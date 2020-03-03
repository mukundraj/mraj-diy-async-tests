//---------------------------------------------------------------------------
//
// diy2-vtk7 parallel particle advection trace classes
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------
#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi/datatypes.hpp>
#include <diy/io/bov.hpp>
#include <diy/pick.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/point.hpp>

using namespace std;

typedef diy::DiscreteBounds            Bounds;
typedef diy::RegularGridLink           RGLink;
typedef diy::RegularDecomposer<Bounds> Decomposer;

typedef     diy::ContinuousBounds       CBounds;
typedef diy::RegularDecomposer<CBounds> CDecomposer;

// incremental stats
// NB, default initialization is 0.0 for all members
struct Stats
{
    double cur_mean_time;
    double prev_mean_time;
    double cur_std_time;
    double prev_std_time;
    double cur_mean_ncalls;
    double prev_mean_ncalls;
    double cur_std_ncalls;
    double prev_std_ncalls;
    double cur_callback_time;
    double cur_mean_callback_time;
    double prev_mean_callback_time;
};

// one point
struct Pt
{
    diy::Point<float, 3>    coords;
};

// whether a point is inside given bounds
// on the boundary is considered inside
bool inside(const Pt& pt, const Bounds bounds)
{
    for (int i = 0; i < 3; i++)
        if (pt.coords[i] < (float)(bounds.min[i]) || pt.coords[i] > (float)(bounds.max[i]))
            return false;
    return true;
}

// continuous domain version of above function
bool cinside(const Pt& pt, const CBounds bounds)
{
    for (int i = 0; i < 3; i++)
        if (pt.coords[i] < (float)(bounds.min[i]) || pt.coords[i] > (float)(bounds.max[i]))
            return false;
    return true;
}



// one end point of a particle trace segment
struct EndPt
{
    int  pid;                                // particle ID, unique within a block
    Pt   pt;                                 // end pointof the trace
    int  gid;                                // block gid of seed particle (start) of this trace
    int  nsteps;                             // number of steps this particle went so far
    bool predonly;                        // whether point is advected in prediction round only
    double st_time;                         // start time of a particle
    int esteps;                             // expected number of steps (decided priority)
    Pt pt_home;                             // coordinates of home point (only used during prediction run)

    const float& operator [](int i) const { return pt.coords[i]; }
    float& operator [](int i)             { return pt.coords[i]; }

    EndPt()
        {
            pid      = 0;
            gid      = 0;
            nsteps   = 0;
            predonly = 0;
            esteps   = 0;
        }
    EndPt(struct Segment& s);                // extract the end point of a segment
};

struct CompareEndPt
{
    bool operator()(EndPt &e1, EndPt &e2)
    {
        return e1.esteps < e2.esteps;
    }
};

// one segment of a particle trace (trajectory)
struct Segment
{
    int pid;        // particle ID, unique within a block
    vector<Pt> pts; // points along trace
    int gid;        // block gid of seed particle (start) of this trace

    Segment()
    {
        pid = 0;
        gid = 0;
    }
    Segment(EndPt& p)                        // construct a segment from one point
        {
            pid      = p.pid;
            gid      = p.gid;
            pid      = p.pid;
            Pt pt    { { p[0], p[1], p[2] } };
            pts.push_back(pt);
        }
     
    Segment(EndPt *p) // construct a segment from one point
        {
            pid = p->pid;
            gid = p->gid;
            pid = p->pid;
            Pt pt{{p->pt.coords[0], p->pt.coords[1], p->pt.coords[2]}};
            pts.push_back(pt);
        }

        // whether end point is inside given bounds
        // on the boundary is considered inside
        bool inside(const int lb[3], const int ub[3]) const
        {
            for (int i = 0; i < 3; i++)
                if (pts.back().coords[i] < (float)(lb[i]) || pts.back().coords[i] > (float)(ub[i]))
                    return false;
            return true;
        }
};

// following constructor defined out of line because references Segment, which needed
// to be defined first
EndPt::
EndPt(Segment& s)                       // extract the end point of a segment
{
    pid = s.pid;
    gid = s.gid;
    pt.coords[0] = s.pts.back().coords[0];
    pt.coords[1] = s.pts.back().coords[1];
    pt.coords[2] = s.pts.back().coords[2];
}

// specialize the serialization of a segment
namespace diy
{
    template<>
    struct Serialization<Segment>
    {
        static
        void save(diy::BinaryBuffer& bb, const Segment& x)
            {
                diy::Serialization<int>::save(bb, x.pid);
                diy::Serialization< vector <Pt> >::
                    save(bb, static_cast< const vector<Pt>& >(x.pts));
                diy::Serialization<int>::save(bb, x.gid);
            }
        static
        void load(diy::BinaryBuffer& bb, Segment& x)
            {
                diy::Serialization<int>::load(bb, x.pid);
                diy::Serialization< vector<Pt> >::
                    load(bb, static_cast< vector<Pt>& >(x.pts));
                diy::Serialization<int>::load(bb, x.gid);
            }
    };
}
