//---------------------------------------------------------------------------
//
// diy2-vtk7 parallel particle advection
//
// original advection kernel courtesy Hanqi Guo
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
//
// copied and modified with permission by Tom Peterka
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

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkStreamTracer.h>
#include <vtkSOADataArrayTemplate.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkSynchronizedRenderWindows.h>
#include <vtkCompositedSynchronizedRenderers.h>
#include <vtkMPIController.h>
#include <vtkProcessIdScalars.h>
#include <vtkVersion.h>

#include <cassert>
#include <cstring>

#include "../opts.h"
#include "ptrace.hpp"
#include "block.hpp"

#include "advect.h"

#include <pnetcdf.h>

using namespace std;

#define IEXCHANGE 1


static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

// add a block to the master
struct AddBlock
{
    AddBlock(diy::Master &master_) :
        master(master_)
    {}

    Block* operator()(int gid,
                      const Bounds& core,
                      const Bounds& bounds,
                      const Bounds& domain,
                      const RGLink& link) const
    {
        Block *b       = new Block;
        RGLink *l      = new RGLink(link);
        diy::Master &m = const_cast<diy::Master&>(master);
        m.add(gid, b, l);
        return b;
    }

    diy::Master& master;
};

// add a block to the master and read input data
struct AddAndRead : public AddBlock
{
    AddAndRead(diy::Master& m,
               const char*  infile_,
               diy::mpi::communicator& world_,
               const float vec_scale_,
               const int hdr_bytes_) :
        AddBlock(m),
        infile(infile_),
        world(world_),
        vec_scale(vec_scale_),
        hdr_bytes(hdr_bytes_) {}

    void operator()(int gid,
                    const Bounds& core,
                    const Bounds& bounds,
                    const Bounds& domain,
                    const RGLink& link) const
    {
        Block* b = AddBlock::operator()(gid, core, bounds, domain, link);
        MPI_Offset *start, *count;
        float *data_u=NULL, *data_v=NULL, *data_w=NULL;

        int ncfile, ndims, nvars, ngatts, unlimited;
        int ret;
        ret = ncmpi_open(world, infile, NC_NOWRITE, MPI_INFO_NULL,&ncfile);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);


        // reversed order of shape and bounds needed because the sample data file
        // is linearized in row-major (C) order
        vector<int> shape(3);
        for (size_t i = 0; i < 3; i++)
            shape[2 - i] = domain.max[i] - domain.min[i] + 1;
//        diy::io::BOV reader(in, shape, hdr_bytes);

        Bounds r_bounds;
        r_bounds.min[0] = bounds.min[2];
        r_bounds.max[0] = bounds.max[2];
        r_bounds.min[1] = bounds.min[1];
        r_bounds.max[1] = bounds.max[1];
        r_bounds.min[2] = bounds.min[0];
        r_bounds.max[2] = bounds.max[0];
        printf("%d %d %d %d %d %d \n", r_bounds.min[0], r_bounds.max[0], r_bounds.min[1],r_bounds.max[1],
                r_bounds.min[2], r_bounds.max[2]);


        start = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
        count = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));

        if (ndims==4){
            count[0] = 1;
            count[1] = r_bounds.max[0] - r_bounds.min[0]+1;
            count[2] = r_bounds.max[1] - r_bounds.min[1]+1;
            count[3] = r_bounds.max[2] - r_bounds.min[2]+1;

            start[0] =  0; start[1] = r_bounds.min[0]; start[2] = r_bounds.min[1]; start[3] = r_bounds.min[2];
        }else if(ndims==3){
            printf("in ndims 3\n");
            count[0] = r_bounds.max[0] - r_bounds.min[0]+1;
            count[1] = r_bounds.max[1] - r_bounds.min[1]+1;
            count[2] = r_bounds.max[2] - r_bounds.min[2]+1;

            start[0] = r_bounds.min[0]; start[1] = r_bounds.min[1]; start[2] = r_bounds.min[2];
        }

        std::cout<<"counts"<<count[0]<<" "<<count[1]<<" "<<count[2]<<"\n";
        std::cout<<"starts"<<start[0]<<" "<<start[1]<<" "<<start[2]<<"\n";

        size_t nvecs =
                (bounds.max[0] - bounds.min[0] + 1) *
                (bounds.max[1] - bounds.min[1] + 1) *
                (bounds.max[2] - bounds.min[2] + 1);
        vector<float> values(nvecs * 3); // temporary contiguous buffer of input vector values

        data_u = (float*) calloc(nvecs, sizeof(float));
        data_v = (float*) calloc(nvecs, sizeof(float));
        data_w = (float*) calloc(nvecs, sizeof(float));
        ret = ncmpi_get_vara_float_all(ncfile, 0, start, count, data_u);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        ret = ncmpi_get_vara_float_all(ncfile, 1, start, count, data_v);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        ret = ncmpi_get_vara_float_all(ncfile, 2, start, count, data_w);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);



        // copy from temp values into block
        b->vel[0] = new float[nvecs];
        b->vel[1] = new float[nvecs];
        b->vel[2] = new float[nvecs];
        b->nvecs = nvecs;
        for (size_t i = 0; i < nvecs; i++)
        {

            b->vel[0][i] = data_u[i] * vec_scale;
            b->vel[1][i] = data_v[i] * vec_scale;
            b->vel[2][i] = data_w[i] * vec_scale;

        }

        ret = ncmpi_close(ncfile);
        free(start);
        free(count);
        free(data_u);
        free(data_v);
        free(data_w);
    }


    const char*	infile;
    diy::mpi::communicator world;
    float vec_scale;
    int hdr_bytes;
};

#if IEXCHANGE==0
void TraceBlock(Block *b,
                const diy::Master::ProxyWithLink &cp,
                const Decomposer&            decomposer,
                const diy::Assigner&         assigner,
                const int                    max_steps,
                const int                    seed_rate,
                const Decomposer::BoolVector share_face)
{
    const int rank              = cp.master()->communicator().rank();
    const int gid               = cp.gid();
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());
    map<diy::BlockID, vector<Pt> > outgoing_pts;

    vector<EndPt> particles;
    map<diy::BlockID, vector<EndPt> > outgoing_endpts;

    const float *vec[3] = {b->vel[0],
                           b->vel[1],
                           b->vel[2]};
    const int   st[3]   = {l->core().min[0],
                           l->core().min[1],
                           l->core().min[2]};
    const int   sz[3]   = {l->core().max[0] - l->core().min[0] + 1,
                           l->core().max[1] - l->core().min[1] + 1,
                           l->core().max[2] - l->core().min[2] + 1};
    const int   gst[3]  = {l->bounds().min[0],
                           l->bounds().min[1],
                           l->bounds().min[2]};
    const int   gsz[3]  = {l->bounds().max[0] - l->bounds().min[0] + 1,
                           l->bounds().max[1] - l->bounds().min[1] + 1,
                           l->bounds().max[2] - l->bounds().min[2] + 1};

    // debug
    // fprintf(stderr, "bounds min[%d %d %d] max[%d %d %d]\n",
    //         l->bounds().min[0], l->bounds().min[1], l->bounds().min[2],
    //         l->bounds().max[0], l->bounds().max[1], l->bounds().max[2]);

    // initialize seed particles first time
    bool first_time = false;
    if (b->init == 0)
    {
        first_time = true;
        // seed particles at every so many grid points
        int sr = (seed_rate < 1 ? 1 : seed_rate);


        for (int i = st[0]; i < st[0] + sz[0]; i += sr)
        {
            // don't duplicate points on block boundaries
            if (share_face[0] && i < decomposer.domain.max[0] && i == l->core().max[0])
                continue;
            for (int j = st[1]; j < st[1] + sz[1]; j += sr)
            {
                // don't duplicate points on block boundaries
                if (share_face[1] && i < decomposer.domain.max[1] && j == l->core().max[1])
                    continue;
                for (int k = st[2]; k < st[2] + sz[2]; k += sr)
                {
                    // don't duplicate points on block boundaries
                    if (share_face[2] && i < decomposer.domain.max[2] && k == l->core().max[2])
                        continue;
                    //                    Pt p;
                    //                    p.coords[0] = i;  p.coords[1] = j;  p.coords[2] = k;
                    //                    b->points->InsertNextPoint(p.coords);

                    EndPt p;
                    p.pid = b->init;
                    p.sid = b->init;
                    p[0] = i;  p[1] = j;  p[2] = k;
                    particles.push_back(p);

                    b->init++; // needed for both

                }
            }
        }
    }

    // dequeue vectors of endpoints, add to seed particles
    vector<int> in;
    cp.incoming(in);
    for (int i = 0; i < in.size(); i++)
    {
        if (cp.incoming(in[i]).buffer.size() > 0)
        {
            vector<EndPt> incoming_endpts;
            cp.dequeue(in[i], incoming_endpts);
            for (size_t j = 0; j < incoming_endpts.size(); j++) {
                incoming_endpts[j].sid++;
                particles.push_back(incoming_endpts[j]);
            }

        }
    }

    // trace particles
    for (int i = 0; i < particles.size(); i++)
    {
        Pt&     cur_p = particles[i].pt; // current end point
        Segment s(particles[i]);         // segment with one point p
        Pt      next_p;                  // coordinates of next end point
        bool    finished = false;

        // trace this segment as far as it will go in the local vector field
        while (trace_3D_rk1(gst, gsz, st, sz, vec, cur_p.coords, 0.5, next_p.coords))
        {
            s.pts.push_back(next_p);
            cur_p = next_p;
            if (s.pts.size() >= max_steps)
            {
                finished = true;
                break;
            }
        }
        b->segments.push_back(s);

        //        Pt        end_pt;
        //        end_pt.coords[0] = next_p[0];
        //        end_pt.coords[1] = next_p[1];
        //        end_pt.coords[2] = next_p[2];


        //        if (!s.inside(decomposer.domain.min, decomposer.domain.max)) // out of global domain
        // if (!inside(end_pt, decomposer.domain))
        if (!inside(next_p, decomposer.domain))
            finished = true;


        if (finished)                    // this segment is done
            b->done++;
        else                             // package segment endpoint in vector for enqueueing
        {
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int> > insert_it(dests, it);
            diy::in(*l, next_p.coords, insert_it, decomposer.domain);
            EndPt out_pt(s);

            // TODO: deal with multiple dests, also match with the ASYNC case
            // either send to first block that is not me or perturb the point along velocity
//          for (size_t j = 0; j < dests.size(); j++)
            for (size_t j = 0; j < 1; j++)
            {
                diy::BlockID bid = l->target(dests[j]);
                outgoing_endpts[bid].push_back(out_pt);
                // fprintf(stderr, "gid %d enqueue [%.3f %.3f %.3f] to gid %d\n",
                //         gid, out_pt[0], out_pt[1], out_pt[2], bid.gid);
            }
        }
    }

    // enqueue the vectors of endpoints
    for (map<diy::BlockID, vector<EndPt> >::const_iterator it =
         outgoing_endpts.begin(); it != outgoing_endpts.end(); it++)
        cp.enqueue(it->first, it->second);



    // stage all_reduce of total initialized and total finished particle traces
    cp.all_reduce(b->init, plus<int>());
    cp.all_reduce(b->done, plus<int>());


}

#endif

#if IEXCHANGE==1
bool trace_segment(Block *b,
                   const diy::Master::ProxyWithLink &cp,
                   const Decomposer&            decomposer,
                   const diy::Assigner&         assigner,
                   const int                    max_steps,
                   const int                    seed_rate,
                   const Decomposer::BoolVector share_face)

{

    const int rank              = cp.master()->communicator().rank();
    const int gid               = cp.gid();
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());
    map<diy::BlockID, vector<Pt> > outgoing_pts;

    vector<EndPt> particles;
    map<diy::BlockID, vector<EndPt> > outgoing_endpts;

    const float *vec[3] = {b->vel[0],
                           b->vel[1],
                           b->vel[2]};
    const int   st[3]   = {l->core().min[0],
                           l->core().min[1],
                           l->core().min[2]};
    const int   sz[3]   = {l->core().max[0] - l->core().min[0] + 1,
                           l->core().max[1] - l->core().min[1] + 1,
                           l->core().max[2] - l->core().min[2] + 1};
    const int   gst[3]  = {l->bounds().min[0],
                           l->bounds().min[1],
                           l->bounds().min[2]};
    const int   gsz[3]  = {l->bounds().max[0] - l->bounds().min[0] + 1,
                           l->bounds().max[1] - l->bounds().min[1] + 1,
                           l->bounds().max[2] - l->bounds().min[2] + 1};

    // debug
    // fprintf(stderr, "bounds min[%d %d %d] max[%d %d %d]\n",
    //         l->bounds().min[0], l->bounds().min[1], l->bounds().min[2],
    //         l->bounds().max[0], l->bounds().max[1], l->bounds().max[2]);

    // initialize seed particles first time
    bool first_time = false;
    if (b->init == 0)
    {
        first_time = true;
        // seed particles at every so many grid points
        int sr = (seed_rate < 1 ? 1 : seed_rate);


        for (int i = st[0]; i < st[0] + sz[0]; i += sr)
        {
            // don't duplicate points on block boundaries
            if (share_face[0] && i < decomposer.domain.max[0] && i == l->core().max[0])
                continue;
            for (int j = st[1]; j < st[1] + sz[1]; j += sr)
            {
                // don't duplicate points on block boundaries
                if (share_face[1] && i < decomposer.domain.max[1] && j == l->core().max[1])
                    continue;
                for (int k = st[2]; k < st[2] + sz[2]; k += sr)
                {
                    // don't duplicate points on block boundaries
                    if (share_face[2] && i < decomposer.domain.max[2] && k == l->core().max[2])
                        continue;
                    //                    Pt p;
                    //                    p.coords[0] = i;  p.coords[1] = j;  p.coords[2] = k;
                    //                    b->points->InsertNextPoint(p.coords);

                    EndPt p;
                    p.pid = b->init;
                    p.sid = b->init;
                    p[0] = i;  p[1] = j;  p[2] = k;
                    particles.push_back(p);

                    b->init++; // needed for both

                }
            }
        }
    }

    //    diy::Link* l = cp.link();

    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        if (cp.incoming(nbr_gid))      // FIXME: make this while
        {
            EndPt incoming_endpt;
            cp.dequeue(nbr_gid, incoming_endpt);
            particles.push_back(incoming_endpt);
        }


    }



    // trace particles
    for (int i = 0; i < particles.size(); i++)
    {
        Pt&     cur_p = particles[i].pt; // current end point
        Segment s(particles[i]);         // segment with one point p
        Pt      next_p;                  // coordinates of next end point
        bool    finished = false;

        // trace this segment as far as it will go in the local vector field
        while (trace_3D_rk1(gst, gsz, st, sz, vec, cur_p.coords, 0.5, next_p.coords))
        {
            s.pts.push_back(next_p);
            cur_p = next_p;
            if (s.pts.size() >= max_steps)
            {
                finished = true;
                break;
            }
        }
        b->segments.push_back(s);

        if (!inside(next_p, decomposer.domain))
            finished = true;


        if (finished){                    // this segment is done
            b->done++;
        }
        else{                               // asyncronously send out segment

            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int> > insert_it(dests, it);
            diy::in(*l, next_p.coords, insert_it, decomposer.domain);
            EndPt out_pt(s);
            //for (size_t j = 0; j < dests.size(); j++)
            for (size_t j = 0; j < 1; j++)
            {
                diy::BlockID bid = l->target(dests[j]);
                cp.enqueue(bid, out_pt);
                //                outgoing_endpts[bid].push_back(out_pt);
                // fprintf(stderr, "gid %d enqueue [%.3f %.3f %.3f] to gid %d\n",
                //         gid, out_pt[0], out_pt[1], out_pt[2], bid.gid);
            }

        }

    }

    return true;
}

#endif

// merge traces at the root block
void merge_traces(void* b_, const diy::ReduceProxy& rp, const diy::RegularMergePartners&)
{
    Block* b = static_cast<Block*>(b_);


    // dequeue and merge
    for (unsigned i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid())               // skip self
            continue;



        vector<Segment> in_traces;
        rp.dequeue(nbr_gid, in_traces);


        // append in_traces to segments
        // TODO: the right way is to sort segments with the same pid into increasing sid order
        // and renumber them into a single trace (streamline, pathline, etc.)
        // for now, we're leaving the traces segmented and disorganized
        b->segments.insert(b->segments.end(), in_traces.begin(), in_traces.end());
    }

    // enqueue
    // NB, for a merge, the out_link size is 1; ie, there is only one target

    //    printf("size %d\n", rp.out_link().size());
    if (rp.out_link().size()){
        int nbr_gid = rp.out_link().target(0).gid;
        if (rp.out_link().size() && nbr_gid != rp.gid()) // skip self
            rp.enqueue(rp.out_link().target(0), b->segments);
    }

}



int main(int argc, char **argv)
{
    string infile;                           // input file name
    Bounds domain;                           // global domain bounds
    int max_steps;                           // max number of steps a particle is allowed to take
    int seed_rate;                           // seed particle every this many grid pts in each dim

    diy::mpi::environment  env(argc, argv);
    diy::mpi::communicator world;

    using namespace opts;

    // defaults
    int nblocks     = world.size();           // total number of global blocks
    int nthreads    = 1;                      // number of threads diy can use
    int mblocks     = -1;                     // number of blocks in memory (-1 = all)
    string prefix   = "./DIY.XXXXXX";         // storage of temp files
    int ndims       = 3;                      // domain dimensions
    float vec_scale = 1.0;                    // vector field scaling factor
    int hdr_bytes   = 0;                      // num bytes header before start of data in infile
    int max_rounds  = 0;                      // max number of rounds to trace (0 = no limit)

    // print vtk version
    if (world.rank() == 0)
        std::cerr << vtkVersion::GetVTKSourceVersion() << std::endl;

    Options ops(argc, argv);
    ops
            >> Option('b', "blocks",     nblocks,    "Total number of blocks to use")
            >> Option('t', "threads",    nthreads,   "Number of threads to use")
            >> Option('m', "in-memory",  mblocks,    "Number of blocks to keep in memory")
            >> Option('s', "storage",    prefix,     "Path for out-of-core storage")
            >> Option('v', "vec-scale",  vec_scale,  "Vector field scaling factor")
            >> Option('h', "hdr-bytes",  hdr_bytes,  "Skip this number bytes header in infile")
            >> Option('r', "max-rounds", max_rounds, "Max number of rounds to trace")
               ;

    if (ops >> Present('h', "help", "show help") ||
            !(ops >> PosOption(infile) >> PosOption(max_steps) >> PosOption(seed_rate)
              >> PosOption(domain.min[0])  >> PosOption(domain.min[1])  >> PosOption(domain.min[2])
              >> PosOption(domain.max[0])  >> PosOption(domain.max[1])  >> PosOption(domain.max[2])))
    {
        if (world.rank() == 0)
        {
            fprintf(stderr, "Usage: %s [OPTIONS] infile mins maxs\n", argv[0]);
            cout << ops;
        }
        return 1;
    }
    diy::FileStorage             storage(prefix);
    diy::Master                  master(world,
                                        nthreads,
                                        mblocks,
                                        &Block::create,
                                        &Block::destroy,
                                        &storage,
                                        &Block::save,
                                        &Block::load);
    diy::RoundRobinAssigner      assigner(world.size(), nblocks);
    Decomposer::BoolVector       share_face;
    Decomposer::BoolVector       wrap;       // defaults to false
    Decomposer::CoordinateVector ghosts;
    ghosts.push_back(2); ghosts.push_back(2); ghosts.push_back(2);
    share_face.push_back(true); share_face.push_back(true); share_face.push_back(true);
    Decomposer decomposer(ndims,
                          domain,
                          assigner.nblocks(),
                          share_face,
                          wrap,
                          ghosts);

    AddAndRead addblock(master, infile.c_str(), world, vec_scale, hdr_bytes);
    decomposer.decompose(world.rank(), assigner, addblock); // does decomposer pass the link to addblock?

    if (world.rank() == 0)
    {
        fprintf(stderr, "input vectors read from file %s\n", infile.c_str());
        fprintf(stderr, "starting particle tracing\n");
    }

    double time_start = MPI_Wtime();

#if IEXCHANGE==0
    // particle tracing for either a maximum number of rounds or, if max_rounds == 0,
    // then for inifinitely many rounds until breaking out when done is true
     int stop = (max_rounds ? max_rounds : 1);
    int incr = (max_rounds ? 1 : 0);
    for (int round = 0; round < stop; round += incr)
    {
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
        { TraceBlock(b,
                     cp,
                     decomposer,
                     assigner,
                     max_steps,
                     seed_rate,
                     share_face); });
         master.exchange();


        int init, done;
        for (int i = 0; i < master.size(); i++)
        {
            init = master.proxy(i).get<int>();
            done = master.proxy(i).get<int>();
        }

        if (world.rank() == 0)
            fprintf(stderr, "round=%d, INIT=%d, DONE=%d\n", round, init, done);

        if (init == done && done != 0)
            break;
    }
#endif

#if IEXCHANGE==1
    master.iexchange([&](Block* b, const diy::Master::ProxyWithLink& cp) -> bool
    { bool val = trace_segment(b,
                               cp,
                               decomposer,
                               assigner,
                               max_steps,
                               seed_rate,
                               share_face);
        return val;
    });

#endif

    double time_end = MPI_Wtime();
    if (world.rank() == 0)
        fprintf(stderr, "wtime=%lf\n", time_end - time_start);


    // merge-reduce traces to one block
    int k = 2;                               // the radix of the k-ary reduction tree
    diy::RegularMergePartners  partners(decomposer, k);
    diy::reduce(master, assigner, partners, &merge_traces);

    if (world.rank() == 0)
    {
        fprintf(stderr, "converting particle traces to vtk polylines and rendering\n");
        ((Block*)master.block(0))->render();
    }




    return 0;
}
