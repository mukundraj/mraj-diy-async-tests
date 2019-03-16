//---------------------------------------------------------------------------
//
// diy parallel particle advection
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
#include <diy/io/shared.hpp>

#ifdef WITH_VTK

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
#include <vtkRenderWindow.h>
#include <vtkProcessIdScalars.h>
#include <vtkVersion.h>

#endif

#include <cassert>
#include <cstring>

#include "../opts.h"
#include "ptrace.hpp"
#include "block.hpp"

#include "advect.h"

#include <fstream>
#include <string.h>
#include <thread>

#ifdef MPE

#include    "mpe.h"

#endif

using namespace std;

#if !defined(IEXCHANGE)
#define IEXCHANGE 0
#endif

int counter = 0;

void InitSeeds(Block*                       b,
               int                          gid,
               const Decomposer&            decomposer,
               const Decomposer::BoolVector share_face,
               diy::RegularLink<Bounds>*    l,
               const int                    sr,
               const int*                   st,
               const int*                   sz,
               int                          synth,
               vector<EndPt>&               particles)
{
    // for synthetic data, see only blocks at -x side of domain, and skip others
    std::vector<int> coords;
    decomposer.gid_to_coords(gid, coords);
    if (synth && coords[0])
        return;

        // debug
//         fmt::print(stderr, "gid={} st=[{} {} {}] sz=[{} {} {}] coords=[{} {} {}]\n",
//                 gid, st[0], st[1], st[2], sz[0], sz[1], sz[2], coords[0], coords[1], coords[2]);

    // for synthetic data, seed only -x side of the block
    int end = synth ? st[0] + 1: st[0] + sz[0];
    for (int i = st[0]; i < end; i += sr)
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

                // debug
//                 fmt::print(stderr, "p=[{} {} {}]\n", p[0], p[1], p[2]);

                b->init++; // needed for both

            }
        }
    }
}

#if IEXCHANGE == 0                        // callback for synchronous exchange version

void TraceBlock(Block*                              b,
                const diy::Master::ProxyWithLink&   cp,
                const Decomposer&                   decomposer,
                const diy::Assigner&                assigner,
                const int                           max_steps,
                const int                           seed_rate,
                const Decomposer::BoolVector        share_face,
                bool                                synth)
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

    // initialize seed particles first time
    if (b->init == 0)
    {
        int sr = (seed_rate < 1 ? 1 : seed_rate);
        InitSeeds(b, gid, decomposer, share_face, l, sr, st, sz, synth, particles);
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
        while (trace_3D_rk1(gst, gsz, st, sz, vec, cur_p.coords.data(), 0.5, next_p.coords.data()))
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

        if (finished)                    // this segment is done
            b->done++;
        else                             // package segment endpoint in vector for enqueueing
        {
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int> > insert_it(dests, it);
            diy::in(*l, next_p.coords, insert_it, decomposer.domain);
            EndPt out_pt(s);
            if (dests.size())
            {
                diy::BlockID bid = l->target(dests[0]); // in case of multiple dests, send to first dest only
                outgoing_endpts[bid].push_back(out_pt);
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

#if IEXCHANGE == 1                                // callback for asynchronous iexchange version

bool trace_segment(Block*                               b,
                   const diy::Master::ProxyWithLink&    icp,
                   const Decomposer&                    decomposer,
                   const diy::Assigner&                 assigner,
                   const int                            max_steps,
                   const int                            seed_rate,
                   const Decomposer::BoolVector         share_face,
                   int                                  synth)

{
    counter++;
    const int rank              = icp.master()->communicator().rank();
    const int gid               = icp.gid();
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(icp.link());
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

    // initialize seed particles first time
    if (b->init == 0)
    {
        int sr = (seed_rate < 1 ? 1 : seed_rate);
        InitSeeds(b, gid, decomposer, share_face, l, sr, st, sz, synth, particles);
    }

    // get incoming points
    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        while (icp.incoming(nbr_gid))
        {
            EndPt incoming_endpt;
            icp.dequeue(nbr_gid, incoming_endpt);
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
        while (trace_3D_rk1(gst, gsz, st, sz, vec, cur_p.coords.data(), 0.5, next_p.coords.data()))
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

        if (finished)                    // this segment is done
            b->done++;
        else                             // enqueue endpoint
        {
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int> > insert_it(dests, it);
            diy::in(*l, next_p.coords, insert_it, decomposer.domain);
            EndPt out_pt(s);
            if (dests.size())
            {
                diy::BlockID bid = l->target(dests[0]); // in case of multiple dests, send to first dest only
                icp.enqueue(bid, out_pt);
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
    Bounds domain {3};                       // global domain bounds
    int max_steps;                           // max number of steps a particle is allowed to take
    int seed_rate;                           // seed particle every this many grid pts in each dim

    diy::mpi::environment  env(argc, argv);
    diy::mpi::communicator world;

    using namespace opts;

    // defaults
    int nblocks             = world.size();     // total number of global blocks
    int nthreads            = 1;                // number of threads diy can use
    int mblocks             = -1;               // number of blocks in memory (-1 = all)
    string prefix           = "./DIY.XXXXXX";   // storage of temp files
    int ndims               = 3;                // domain dimensions
    float vec_scale         = 1.0;              // vector field scaling factor
    int hdr_bytes           = 0;                // num bytes header before start of data in infile
    int max_rounds          = 0;                // max number of rounds to trace (0 = no limit)
    size_t min_queue_size   = 0;                // min queue size (bytes) for iexchange
    size_t max_hold_time    = 0;                // max hold time (microsec) for iexchange
    int synth               = 0;                // generate various synthetic input datasets
    float slow_vel          = 1.0;              // slow velocity for synthetic data
    float fast_vel          = 10.0;             // fast velocity for synthetic data
    int check               = 0;                // write out traces for checking
    std::string log_level   = "info";           // logging level
    int ntrials             = 1;                // number of trials

    Options ops(argc, argv);
    ops
        >> Option('b', "blocks",        nblocks,        "Total number of blocks to use")
        >> Option('t', "threads",       nthreads,       "Number of threads to use")
        >> Option('m', "in-memory",     mblocks,        "Number of blocks to keep in memory")
        >> Option('s', "storage",       prefix,         "Path for out-of-core storage")
        >> Option('v', "vec-scale",     vec_scale,      "Vector field scaling factor")
        >> Option('h', "hdr-bytes",     hdr_bytes,      "Skip this number bytes header in infile")
        >> Option('r', "max-rounds",    max_rounds,     "Max number of rounds to trace")
        >> Option('q', "min-q-size",    min_queue_size, "Minimum queue size (bytes) for iexchange")
        >> Option('o', "max-hold-time", max_hold_time,  "Maximum queue hold time (ms) for iexchange")
        >> Option('x', "synthetic",     synth,          "Generate various synthetic flows")
        >> Option('w', "slow-vel",      slow_vel,       "Slow velocity for synthetic data")
        >> Option('f', "fast-vel",      fast_vel,       "Fast velocity for synthetic data")
        >> Option('c', "check",         check,          "Write out traces for checking")
        >> Option('l', "log",           log_level,      "log level")
        >> Option('n', "trials",        ntrials,        "number of trials")
        ;
    bool fine = ops >> Present("fine", "Use fine-grain icommunicate");

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
//     diy::create_logger(log_level);
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

    // decompose domain
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

    if (synth == 1)
    {
        AddSynthetic1 addsynth(master, slow_vel, fast_vel, decomposer);
        decomposer.decompose(world.rank(), assigner, addsynth);
    }
    else if (synth == 2)
    {
        AddSynthetic2 addsynth(master, slow_vel, fast_vel, decomposer);
        decomposer.decompose(world.rank(), assigner, addsynth);
    }
    else
    {
        AddAndRead addblock(master, infile.c_str(), world, vec_scale, hdr_bytes);
        decomposer.decompose(world.rank(), assigner, addblock);
    }

    if (world.rank() == 0)
    {
        if (synth)
            fprintf(stderr, "input vectors created synthetically\n");
        else
            fprintf(stderr, "input vectors read from file %s\n", infile.c_str());
        fprintf(stderr, "starting particle tracing\n");
    }

    // incremental stats
    double cur_mean_time            = 0.0;
    double prev_mean_time           = 0.0;
    double cur_std_time             = 0.0;
    double prev_std_time            = 0.0;
    double cur_mean_ncalls          = 0.0;
    double prev_mean_ncalls         = 0.0;
    double cur_std_ncalls           = 0.0;
    double prev_std_ncalls          = 0.0;

    size_t nrounds = 0;

    for (int trial = 0; trial < ntrials; trial++)           // number of trials
    {
        // reset the block particle traces, but leave the vector field intact
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
                {
                    b->init     = 0;
                    b->done     = 0;
                    b->nvecs    = 0;
                    b->segments.clear();
                });

        MPI_Barrier(world);
        double time_start = MPI_Wtime();

#if IEXCHANGE == 1

        master.iexchange([&](Block* b, const diy::Master::ProxyWithLink& icp) -> bool
                {
                bool val = trace_segment(b,
                        icp,
                        decomposer,
                        assigner,
                        max_steps,
                        seed_rate,
                        share_face,
                        synth);
                return val;
                }, min_queue_size, max_hold_time, fine);

#else

        // particle tracing for either a maximum number of rounds or, if max_rounds == 0,
        // then for inifinitely many rounds until breaking out when done is true
        int stop = (max_rounds ? max_rounds : 1);
        int incr = (max_rounds ? 1 : 0);

        for (nrounds = 0; nrounds < stop; nrounds += incr)
        {
            master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
                    {
                    TraceBlock(b,
                            cp,
                            decomposer,
                            assigner,
                            max_steps,
                            seed_rate,
                            share_face,
                            synth);
                    });
            master.exchange();

            int init, done;
            for (int i = 0; i < master.size(); i++)
            {
                init = master.proxy(i).get<int>();
                done = master.proxy(i).get<int>();
            }

            if (init == done && done != 0)
            {
                nrounds++;
                break;
            }
        }

#endif

        MPI_Barrier(world);

        double cur_time = MPI_Wtime() - time_start;

        int cur_ncalls = 0;
        MPI_Reduce(&counter, &cur_ncalls, 1, MPI_INT, MPI_SUM, 0, world);

        // update incremental stats
        // Knuth "The Art of Computer Programming, Volume 2: Seminumerical Algorithms", section 4.2.2
        // originally B.P. Welford, Technometrics, 4,(1962), 419-420
        if (trial == 0)
        {
            cur_mean_time   = cur_time;
            prev_mean_time  = cur_time;
            cur_mean_ncalls = cur_ncalls;
            prev_mean_ncalls= cur_ncalls;
            cur_std_time    = 0.0;
            cur_std_ncalls  = 0.0;
        }
        else
        {
            cur_mean_time   = prev_mean_time   + (cur_time   - prev_mean_time)   / (trial + 1);
            cur_mean_ncalls = prev_mean_ncalls + (cur_ncalls - prev_mean_ncalls) / (trial + 1);
            cur_std_time    = prev_std_time    + (cur_time   - prev_mean_time)   * (cur_time   - cur_mean_time);
            cur_std_ncalls  = prev_std_ncalls  + (cur_ncalls - prev_mean_ncalls) * (cur_ncalls - cur_mean_ncalls);
        }
        prev_mean_time      = cur_mean_time;
        prev_mean_ncalls    = cur_mean_ncalls;
        prev_std_time       = cur_std_time;
        prev_std_ncalls     = cur_std_ncalls;

        // debug
        if (world.rank() == 0)
        {
            fmt::print(stderr, "trial {} time {} nrounds {} ncalls {}\n",
                    trial, cur_time, nrounds, cur_ncalls);
        }

        // merge-reduce traces to one block
        int k = 2;                               // the radix of the k-ary reduction tree
        diy::RegularMergePartners  partners(decomposer, k);
        diy::reduce(master, assigner, partners, &merge_traces);

        // rendering
#ifdef WITH_VTK
#if 0                       // render one block with all traces
        if (world.rank() == 0)
        {
            fprintf(stderr, "converting particle traces to vtk polylines and rendering 1 block only\n");
            ((Block*)master.block(0))->render();
        }
#else                       // debugging: render individual blocks
        if (world.rank() == 0)
            fprintf(stderr, "converting particle traces to vtk polylines and rendering all blocks\n");
        master.foreach(&Block::render_block);
#endif
#endif

        // output profile
#ifdef DIY_PROFILE
#if IEXCHANGE == 1
            diy::io::SharedOutFile prof_out(fmt::format("profile-iexchange-p{}-b{}.txt", world.size(), nblocks), world);
#else
            diy::io::SharedOutFile prof_out(fmt::format("profile-exchange-p{}-b{}.txt", world.size(), nblocks), world);
#endif
            master.prof.output(prof_out, std::to_string(world.rank()));
            prof_out.close();
#endif

    }           // number of trials

    // print stats
    if (world.rank() == 0)
    {
        // ncalls: number of calls in case of IEXCHANGE
        // nrounds: number of rounds in case of EXCHANGE
        fmt::print(stderr, "---------- stats ----------\n");
#if IEXCHANGE == 1
        fmt::print(stderr, "using iexchange\n");
        fmt::print(stderr, "min queue size (bytes):\t{}\n",     min_queue_size);
        fmt::print(stderr, "max hold time (micro s):\t{}\n",    max_hold_time);
#else
        fmt::print(stderr, "using exchange\n");
#endif
        fmt::print(stderr, "seed rate:\t\t\t{}\n",              seed_rate);
        fmt::print(stderr, "nprocs:\t\t\t{}\n",                 world.size());
        fmt::print(stderr, "nblocks:\t\t\t{}\n",                nblocks);
        fmt::print(stderr, "ntrials:\t\t\t{}\n",                ntrials);
        fmt::print(stderr, "mean time (s):\t\t{}\n",            cur_mean_time);
        fmt::print(stderr, "std dev time (s):\t\t{}\n",         ntrials > 1 ? sqrt(cur_std_time / (ntrials - 1)) : 0.0);
#if IEXCHANGE == 1
        fprintf(stderr, "mean # callbacks:\t\t%.0lf\n",         cur_mean_ncalls);
        fprintf(stderr, "std dev # callbacks:\t%.0lf\n",        ntrials > 1 ? sqrt(cur_std_ncalls / (ntrials - 1)) : 0.0);
#else
        fmt::print(stderr, "# rounds:\t\t\t{}\n",               nrounds);
#endif

        char infile[256];           // profile file name
#if IEXCHANGE == 1
        sprintf(infile, "profile-iexchange-p%d-b%d.txt",        world.size(), nblocks);
#else
        sprintf(infile, "profile-exchange-p%d-b%d.txt",         world.size(), nblocks);
#endif

        // count number of occurences of send-same-rank plus send-different rank in profile
        // NB: only for the last trial; profile is overwritten for each trial
#ifdef DIY_PROFILE
            char cmd[256];
            sprintf(cmd, "grep -c '>send-.*-rank' %s", infile);
            char buf[256];
            size_t ncalls;
            FILE* pipe(popen(cmd, "r"));
            if (!pipe)
            {
                fprintf(stderr, "Error: popen failed\n");
                abort();
            }
            while (!feof(pipe))
                if (fgets(buf, 128, pipe) != NULL)
                    ncalls = atol(buf);
            pclose(pipe);
            fmt::print(stderr, "# send calls:\t\t{}\n",             ncalls);
#endif

        fmt::print(stderr, "---------------------------\n");
    }

    // write trajectory segments out in order to validate that they are identical
    if (check)
    {
        std::string filename;
        if(IEXCHANGE==0)
            filename = "exchange.txt";
        else
            filename = "iexchange.txt";
        ((Block*)master.block(0))->write_segments(filename);
    }

    return 0;
}
