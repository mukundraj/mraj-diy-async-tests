//---------------------------------------------------------------------------
//
// diy parallel particle advection
//
// original advection kernel courtesy Hanqi Guo
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
//
// copied and modified with permission by Tom Peterka and Mukund Raj
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
#include "lerp.hpp"
#include "utils.hpp"

#include <fstream>
#include <string.h>
#include <thread>

using namespace std;

void InitSeeds(Block*                       b,
               int                          gid,
               const Decomposer&            decomposer,
               diy::RegularLink<Bounds>*    l,
               float                        sr,
               int                          synth)
{
    // for synthetic data, seed only blocks at -x side of domain, and skip others
    std::vector<int> coords;
    decomposer.gid_to_coords(gid, coords);
    if (synth && coords[0])
        return;

    // seed interior of block (1 step inside of all block core boundaries)
    // particles at block boundaries are ambiguous as to which block they belong
    // best to avoid this case when seeding them

    // for synthetic data, seed only -x side of the block
    float epsilon = 1.0e-10;            // a small increment
    float end = synth ? l->core().min[0] + 1.0 + epsilon: l->core().max[0] - 1.0;

    // seed the block
    for (float i = float(l->core().min[0]) + 1.0; i <= end; i += sr)
    {
        for (float j = l->core().min[1] + 1.0; j <= l->core().max[1] - 1.0; j += sr)
        {
            for (float k = l->core().min[2] + 1.0; k <= l->core().max[2] - 1.0; k += sr)
            {
                EndPt p;
                p.pid = b->init;
                p.gid = gid;
                p[0] = i;  p[1] = j;  p[2] = k;
                b->particles.push_back(p);
                b->init++;
            }
        }
    }
}

// common to both exchange and iexchange
void trace_particles(Block*                             b,
                     const diy::Master::ProxyWithLink&  cp,
                     const Decomposer&                  decomposer,
                     const int                          max_steps,
                     map<diy::BlockID, vector<EndPt> >& outgoing_endpts)
{
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());

    const float *vec[3] = {b->vel[0],           // shallow pointer copy
                           b->vel[1],
                           b->vel[2]};
    const int   st[3]   = {l->core().min[0],
                           l->core().min[1],
                           l->core().min[2]};
    const int   sz[3]   = {l->core().max[0] - l->core().min[0] + 1,
                           l->core().max[1] - l->core().min[1] + 1,
                           l->core().max[2] - l->core().min[2] + 1};

    for (auto i = 0; i < b->particles.size(); i++)
    {
        Pt&     cur_p = b->particles[i].pt;     // current end point
        Segment s(b->particles[i]);             // segment with one point p
        Pt      next_p;                         // coordinates of next end point
        bool    finished = false;

        // trace this segment until it leaves the block core (no ghost)
        while (advect_rk4(st, sz, vec, cur_p.coords.data(), 0.5, next_p.coords.data()))
        {
            b->particles[i].nsteps++;
            s.pts.push_back(next_p);
            cur_p = next_p;
            if (b->particles[i].nsteps >= max_steps)
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
        else                             // find destination of segment endpoint
        {
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int> > insert_it(dests, it);

            utl::in(*l, next_p.coords, insert_it, decomposer.domain, 1);

            EndPt out_pt(s);
            out_pt.nsteps = b->particles[i].nsteps;
            if (dests.size())
            {
                diy::BlockID bid = l->target(dests[0]); // in case of multiple dests, send to first dest only

                // debug
//                 fmt::print(stderr, "gid {} enq to gid {}\n", cp.gid(), bid.gid);

                if (IEXCHANGE)                          // enqueuing single endpoint allows fine-grain iexchange if desired
                    cp.enqueue(bid, out_pt);
                else
                    outgoing_endpts[bid].push_back(out_pt); // vector of endpoints
            }
        }
    }
}

void deq_incoming_exchange(Block*                               b,
                           const diy::Master::ProxyWithLink&    cp)
{
    vector<int> in;
    cp.incoming(in);
    for (int i = 0; i < in.size(); i++)
    {
        if (cp.incoming(in[i]).buffer.size() > 0)
        {
            vector<EndPt> incoming_endpts;
            cp.dequeue(in[i], incoming_endpts);
            for (size_t j = 0; j < incoming_endpts.size(); j++)
                b->particles.push_back(incoming_endpts[j]);
        }
    }
}

void deq_incoming_iexchange(Block*                              b,
                            const diy::Master::ProxyWithLink&   cp)
{
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());
    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        while (cp.incoming(nbr_gid))
        {
            EndPt incoming_endpt;
            cp.dequeue(nbr_gid, incoming_endpt);
            b->particles.push_back(incoming_endpt);
        }
    }
}

// common to both exchange and iexchange
void trace_block(Block*                              b,
                 const diy::Master::ProxyWithLink&   cp,
                 const Decomposer&                   decomposer,
                 const diy::Assigner&                assigner,
                 const int                           max_steps,
                 const float                         seed_rate,
                 const Decomposer::BoolVector        share_face,
                 bool                                synth,
                 map<diy::BlockID, vector<EndPt>>&   outgoing_endpts)
{
    const int gid               = cp.gid();
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds>*>(cp.link());
    b->particles.clear();

    const int   st[3]   = {l->core().min[0],
                           l->core().min[1],
                           l->core().min[2]};
    const int   sz[3]   = {l->core().max[0] - l->core().min[0] + 1,
                           l->core().max[1] - l->core().min[1] + 1,
                           l->core().max[2] - l->core().min[2] + 1};

    // initialize seed particles first time
    if (b->init == 0)
        InitSeeds(b, gid, decomposer, l, seed_rate, synth);

    // dequeue incoming points and trace particles
    if (IEXCHANGE)
    {
        do
        {
            deq_incoming_iexchange(b, cp);
            trace_particles(b, cp, decomposer, max_steps, outgoing_endpts);
            b->particles.clear();
        } while (cp.fill_incoming());
    }
    else
    {
        deq_incoming_exchange(b, cp);
        trace_particles(b, cp, decomposer, max_steps, outgoing_endpts);
    }
}

void trace_block_exchange(Block*                              b,
                          const diy::Master::ProxyWithLink&   cp,
                          const Decomposer&                   decomposer,
                          const diy::Assigner&                assigner,
                          const int                           max_steps,
                          const float                         seed_rate,
                          const Decomposer::BoolVector        share_face,
                          bool                                synth)
{
    map<diy::BlockID, vector<EndPt> > outgoing_endpts;

    trace_block(b, cp, decomposer, assigner, max_steps, seed_rate, share_face, synth, outgoing_endpts);

    // enqueue the vectors of endpoints
    for (map<diy::BlockID, vector<EndPt> >::const_iterator it = outgoing_endpts.begin(); it != outgoing_endpts.end(); it++)
        cp.enqueue(it->first, it->second);

    // stage all_reduce of total initialized and total finished particle traces
    cp.all_reduce(b->particles.size(), plus<size_t>());
}

bool trace_block_iexchange(Block*                               b,
                           const diy::Master::ProxyWithLink&    cp,
                           const Decomposer&                    decomposer,
                           const diy::Assigner&                 assigner,
                           const int                            max_steps,
                           const float                          seed_rate,
                           const Decomposer::BoolVector         share_face,
                           int                                  synth)
{
    map<diy::BlockID, vector<EndPt> > outgoing_endpts;  // needed to call trace_particles() but otherwise unused in iexchange
    trace_block(b, cp, decomposer, assigner, max_steps, seed_rate, share_face, synth, outgoing_endpts);
    return true;
}

// merge traces at the root block
void merge_traces(void* b_, const diy::ReduceProxy& rp, const diy::RegularMergePartners&)
{
    Block* b = static_cast<Block*>(b_);

    // dequeue and merge
    for (unsigned i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid())                    // skip self
            continue;

        vector<Segment> in_traces;
        rp.dequeue(nbr_gid, in_traces);

        // append in_traces to segments, leaving trajectories segmented and disorganized
        // eventually could sort into continuous long trajectories, but not necessary at this time
        b->segments.insert(b->segments.end(), in_traces.begin(), in_traces.end());
    }

    // enqueue
    if (rp.out_link().size())
    {
        int nbr_gid = rp.out_link().target(0).gid;  // for a merge, the out_link size is 1; ie, there is only one target
        if (nbr_gid != rp.gid())                    // skip self
            rp.enqueue(rp.out_link().target(0), b->segments);
    }
}

// debugging: catches segfaults and dumps a backtrace
// copied from Dmitriy's henson/henson-chai.cpp

#include <cxxabi.h>
#include <signal.h>
#include <execinfo.h>
bool abort_on_segfault_ = true;

void catch_sig(int signum)
{
    fmt::print(stderr, "caught signal {}\n", signum);

    // print backtrace
    void*   callstack[128];
    int     frames      = backtrace(callstack, 128);
    char**  strs        = backtrace_symbols(callstack, frames);

    size_t funcnamesize = 256;
    char*  funcname     = (char*) malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < frames; i++)
    {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for (char *p = strs[i]; *p; ++p)
        {
            if (*p == '(')
                begin_name = p;
            else if (*p == '+')
                begin_offset = p;
            else if (*p == ')' && begin_offset)
            {
                end_offset = p;
                break;
            }
        }

        if (begin_name && begin_offset && end_offset && begin_name < begin_offset)
        {
            *begin_name++   = '\0';
            *begin_offset++ = '\0';
            *end_offset     = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply __cxa_demangle():

            int status;
            char* ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
            if (status == 0)
            {
                funcname = ret; // use possibly realloc()-ed string
                fmt::print(stderr, "  {} : {}+{}", strs[i], funcname, begin_offset);
            } else
            {
                // demangling failed. Output function name as a C function with no arguments.
                fmt::print(stderr, "  {} : {}()+{}", strs[i], begin_name, begin_offset);
            }
        }
        else
        {
            // couldn't parse the line? print the whole line.
            fmt::print(stderr, "  {}", strs[i]);
        }
    }

    free(funcname);
    free(strs);

    //for (int i = 0; i < frames; ++i)
    //    logger->critical("{}", strs[i]);

    signal(signum, SIG_DFL);    // restore the default signal
    if (abort_on_segfault_)
        MPI_Abort(MPI_COMM_WORLD, 1);
}

// update incremental stats
// ref: Knuth "The Art of Computer Programming, Volume 2: Seminumerical Algorithms", section 4.2.2
// originally B.P. Welford, Technometrics, 4,(1962), 419-420
void update_stats(
        int                             trial,
        double                          time_start,
        int                             ncalls,
        const diy::mpi::communicator&   world,
        Stats&                          stats)
{
    double cur_time = MPI_Wtime() - time_start;
    int cur_ncalls  = 0;
    MPI_Reduce(&ncalls, &cur_ncalls, 1, MPI_INT, MPI_SUM, 0, world);

    if (trial == 0)
    {
        stats.cur_mean_time               = cur_time;
        stats.prev_mean_time              = cur_time;
        stats.cur_mean_ncalls             = cur_ncalls;
        stats.prev_mean_ncalls            = cur_ncalls;
        stats.cur_mean_callback_time      = stats.cur_callback_time;
        stats.prev_mean_callback_time     = stats.cur_callback_time;
        stats.cur_std_time                = 0.0;
        stats.cur_std_ncalls              = 0.0;
    }
    else
    {
        stats.cur_mean_time           = stats.prev_mean_time            +
                                        (cur_time                 - stats.prev_mean_time)           / (trial + 1);
        stats.cur_mean_ncalls         = stats.prev_mean_ncalls          +
                                        (cur_ncalls         - stats.prev_mean_ncalls)         / (trial + 1);
        stats.cur_mean_callback_time  = stats.prev_mean_callback_time   +
                                        (stats.cur_callback_time  - stats.prev_mean_callback_time)  / (trial + 1);
        stats.cur_std_time            = stats.prev_std_time             +
                                        (cur_time                 - stats.prev_mean_time)   * (cur_time         - stats.cur_mean_time);
        stats.cur_std_ncalls          = stats.prev_std_ncalls           +
                                        (cur_ncalls         - stats.prev_mean_ncalls) * (cur_ncalls - stats.cur_mean_ncalls);
    }
    stats.prev_mean_time              = stats.cur_mean_time;
    stats.prev_mean_ncalls            = stats.cur_mean_ncalls;
    stats.prev_mean_callback_time     = stats.cur_mean_callback_time;
    stats.prev_std_time               = stats.cur_std_time;
    stats.prev_std_ncalls             = stats.cur_std_ncalls;

    // debug
//     if (world.rank() == 0)
//     {
//         fmt::print(stderr, "trial {} time {} callback time {} ncalls {}\n",
//                 trial, cur_time, stats.cur_callback_time, cur_ncalls);
//     }
}

void print_results(
        float           seed_rate,
        int             nprocs,
        int             nblocks,
        int             ntrials,
        int             nrounds,
        const Stats&    stats)
{
    fmt::print(stderr, "---------- stats ----------\n");
    if (IEXCHANGE)
        fmt::print(stderr, "using iexchange\n");
    else
        fmt::print(stderr, "using exchange\n");
    fmt::print(stderr, "seed rate:                       {}\n", seed_rate);
    fmt::print(stderr, "nprocs:                          {}\n", nprocs);
    fmt::print(stderr, "nblocks:                         {}\n", nblocks);
    fmt::print(stderr, "ntrials:                         {}\n", ntrials);
    fmt::print(stderr, "mean time (s):                   {}\n", stats.cur_mean_time);
    fmt::print(stderr, "std dev time (s):                {}\n", ntrials > 1 ? sqrt(stats.cur_std_time / (ntrials - 1)) : 0.0);
    if (IEXCHANGE)
    {
        fprintf(stderr,    "mean # callbacks:                %.0lf\n",  stats.cur_mean_ncalls);
        fprintf(stderr,    "std dev # callbacks:             %.0lf\n",  ntrials > 1 ? sqrt(stats.cur_std_ncalls / (ntrials - 1)) : 0.0);
    }
    else
    {
        fmt::print(stderr, "# rounds:                        {}\n",     nrounds);
        fmt::print(stderr, "mean callback (advect) time (s): {}\n",     stats.cur_mean_callback_time);
    }
    fmt::print(stderr, "---------------------------\n");
}

void print_exceeded_max_rounds(diy::Master& master)
{
    if (master.communicator().rank() == 0)
        fmt::print(stderr, "*** Warning: max # rounds for exchange has been reached. ***\n");

    // debug: print unterminated particles
    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        if (b->particles.size() > 0)
        {
            fmt::print(stderr, "gid = {}, particles size = {}\n", cp.gid(), b->particles.size());
            auto* l = static_cast<RGLink*>(cp.link());
            fmt::print(stderr, "  core = {} - {}, bounds = {} - {}\n",
                    l->core().min,   l->core().max,
                    l->bounds().min, l->bounds().max);
            for (auto& p : b->particles)
                fmt::print(stderr, "  {}\n", p.pt.coords);
        }
    });
}

#ifdef WITH_VTK

void render_traces(
        diy::Master&        master,
        diy::Assigner&      assigner,
        Decomposer&         decomposer,
        bool                merge)                      // merge traces to one block
{
    if (merge)
    {
        // merge-reduce traces to one block
        int k = 2;                               // the radix of the k-ary reduction tree
        diy::RegularMergePartners  partners(decomposer, k);
        diy::reduce(master, assigner, partners, &merge_traces);

        if (master.communicator().rank() == 0)
        {
            fprintf(stderr, "converting particle traces to vtk polylines and rendering 1 block only\n");
            ((Block*)master.block(0))->render();
        }
    }
    else
    {
        if (master.communicator().rank() == 0)
            fprintf(stderr, "converting particle traces to vtk polylines and rendering all blocks\n");
        master.foreach(&Block::render_block);
    }
}

#endif

void write_traces(
        diy::Master&        master,
        diy::Assigner&      assigner,
        Decomposer&         decomposer)
{
    // merge-reduce traces to one block
    int k = 2;                               // the radix of the k-ary reduction tree
    diy::RegularMergePartners  partners(decomposer, k);
    diy::reduce(master, assigner, partners, &merge_traces);

    if (master.communicator().rank() == 0)
    {
        fprintf(stderr, "Check is turned on: merging traces to one block and writing them to disk\n");
        std::string filename;
        if (IEXCHANGE)
            filename = "iexchange.txt";
        else
            filename = "exchange.txt";
        ((Block*)master.block(0))->write_segments(filename);
    }
}

void output_profile(
        diy::Master&            master,
        int                     nblocks)
{
    if (IEXCHANGE)
    {
        diy::io::SharedOutFile prof_out(fmt::format("profile-iexchange-p{}-b{}.txt",
                    master.communicator().size(), nblocks), master.communicator());
        master.prof.output(prof_out, std::to_string(master.communicator().rank()));
        prof_out.close();
    }
    else
    {
        diy::io::SharedOutFile prof_out(fmt::format("profile-exchange-p{}-b{}.txt",
                    master.communicator().size(), nblocks), master.communicator());
        master.prof.output(prof_out, std::to_string(master.communicator().rank()));
        prof_out.close();
    }
}

int main(int argc, char **argv)
{
    signal(SIGSEGV, catch_sig);                 // catch segfault

    string infile;                              // input file name
    Bounds domain {3};                          // global domain bounds
    int max_steps;                              // max number of steps a particle is allowed to take
    float seed_rate;                            // seed particle every this many grid pts in each dim

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
    bool merged_traces      = false;            // traces have already been merged to one block

    // command-line ags
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
    }

    Stats stats;                        // incremental stats, default initialized to 0's
    int nrounds;

    // run the trials
    for (int trial = 0; trial < ntrials; trial++)
    {
        int ncalls = 0;

        // debug
//         if (world.rank() == 0)
//             fprintf(stderr, "started particle tracing trial %d\n", trial);

        // reset the block particle traces, but leave the vector field intact
        master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
                {
                    b->init     = 0;
                    b->done     = 0;
                    b->segments.clear();
                    b->particles.clear();
                });

        world.barrier();
        double time_start = MPI_Wtime();

        if (IEXCHANGE)
        {
            // combined advection and exchange
            master.iexchange([&](Block* b, const diy::Master::ProxyWithLink& icp) -> bool
            {
                ncalls++;
                bool val = trace_block_iexchange(b,
                           icp,
                           decomposer,
                           assigner,
                           max_steps,
                           seed_rate,
                           share_face,
                           synth);
                return val;
            });
        }
        else    // exchange
        {
            // particle tracing for either a maximum number of rounds or, if max_rounds == 0,
            // then for inifinitely many rounds until breaking out when done is true
            int stop = (max_rounds ? max_rounds : 1);
            int incr = (max_rounds ? 1 : 0);

            nrounds                 = 0;
            stats.cur_callback_time = 0.0;
            for (int round = 0; round < stop; round += incr)
            {
                nrounds++;

                // advect
                double t0 = MPI_Wtime();
                master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
                {
                    trace_block_exchange(b,
                                         cp,
                                         decomposer,
                                         assigner,
                                         max_steps,
                                         seed_rate,
                                         share_face,
                                         synth);
                });
                stats.cur_callback_time += (MPI_Wtime() - t0);

                // exchange
                master.exchange();

                // determine if all particles are done
                size_t remaining;
                for (int i = 0; i < master.size(); i++)
                    remaining = master.proxy(i).get<size_t>();
                if (remaining == 0)
                    break;
            }   // rounds

            // debug: exceeded number of rounds for exchange
            if (nrounds == max_rounds)
                print_exceeded_max_rounds(master);
        }

        world.barrier();

        // debug
//         if (world.rank() == 0)
//             fprintf(stderr, "finished particle tracing trial %d\n", trial);
//         master.prof.totals().output(std::cerr);

        update_stats(trial, time_start, ncalls, world, stats);

#ifdef WITH_VTK
        render_traces(master, assigner, decomposer, true);
#endif

#ifdef DIY_PROFILE
        output_profile(master, nblocks);
#endif

    }           // number of trials

    if (world.rank() == 0)
        print_results(seed_rate, world.size(), nblocks, ntrials, nrounds, stats);

    // write trajectory segments for validation
    if (check)
        write_traces(master, assigner, decomposer);

    // debug
//     master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
//     {
//         b->show_geometry(cp);
//     });

    return 0;
}
