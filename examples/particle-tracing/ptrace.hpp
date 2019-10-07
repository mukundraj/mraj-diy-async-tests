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
        if (pt.coords[i] < bounds.min[i] || pt.coords[i] > bounds.max[i])
            return false;
    return true;
}

// one end point of a particle trace segment
struct EndPt
{
    int  pid;                                // particle ID, for now unique only within a block, not globally
    Pt   pt;                                 // end pointof the trace
    int  sid;                                // segment ID of this part of the trace, for now same as pid
    int  nsteps;                             // number of steps this particle went so far

    const float& operator [](int i) const { return pt.coords[i]; }
    float& operator [](int i)             { return pt.coords[i]; }

    EndPt()
        {
            pid      = 0;
            sid      = 0;
            nsteps   = 0;
        }
    EndPt(struct Segment& s);                // extract the end point of a segment
};

// one segment of a particle trace (trajectory)
struct Segment
{
    int        pid;                          // particle ID, for now unique only within a block, not globally
    vector<Pt> pts;                          // points along trace
    int        sid;                          // segment ID of this part of the trace, for now same as pid

    Segment()
        {
            pid      = 0;
            sid      = 0;
        }
    Segment(EndPt& p)                        // construct a segment from one point
        {
            pid      = p.pid;
            sid      = p.sid;
            Pt pt    { { p[0], p[1], p[2] } };
            pts.push_back(pt);
        }

    // whether end point is inside given bounds
    // on the boundary is considered inside
    bool inside(const int lb[3], const int ub[3]) const
        {
            for (int i = 0; i < 3; i++)
                if (pts.back().coords[i] < lb[i] || pts.back().coords[i] > ub[i])
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
    sid = s.sid;
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
                diy::Serialization<int>::save(bb, x.sid);
            }
        static
        void load(diy::BinaryBuffer& bb, Segment& x)
            {
                diy::Serialization<int>::load(bb, x.pid);
                diy::Serialization< vector<Pt> >::
                    load(bb, static_cast< vector<Pt>& >(x.pts));
                diy::Serialization<int>::load(bb, x.sid);
            }
    };
}

// TODO: uncomment and convert from segment to vtk
// specialize the serialization of a segment
// namespace diy
// {
//     template<>
//     struct Serialization<Segment>
//     {
//         static
//         void save(diy::BinaryBuffer& bb, const Segment& x)
//             {
//                 diy::Serialization<int>::save(bb, x.pid);
//                 diy::Serialization< vector <Pt> >::
//                     save(bb, static_cast< const vector<Pt>& >(x.pts));
//                 diy::Serialization<int>::save(bb, x.sid);
//             }
//         static
//         void load(diy::BinaryBuffer& bb, Segment& x)
//             {
//                 diy::Serialization<int>::load(bb, x.pid);
//                 diy::Serialization< vector<Pt> >::
//                     load(bb, static_cast< vector<Pt>& >(x.pts));
//                 diy::Serialization<int>::load(bb, x.sid);
//             }
//     };
// }
