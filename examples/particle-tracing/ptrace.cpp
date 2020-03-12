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

#include <diy/algorithms.hpp>

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

#include "misc.h"
#include <random>
#include <string>
#include <iterator>
#include <algorithm>

using namespace std;

void InitSeeds(Block *b,
               int gid,
               const Decomposer &decomposer,
               diy::RegularLink<Bounds> *l,
               float sr,
               int synth)
{
    // for synthetic data, seed only blocks at -x side of domain, and skip others
    std::vector<int> coords;
    decomposer.gid_to_coords(gid, coords);
    if (synth && coords[0])
        return;

    // for synthetic data, seed only -x side of the block
    float end = synth ? l->core().min[0] + 1.0 : l->core().max[0];

    // seed the block
    for (float i = float(l->core().min[0]); i < end; i += sr)
    {
        for (float j = l->core().min[1]; j < l->core().max[1]; j += sr)
        {
            for (float k = l->core().min[2]; k < l->core().max[2]; k += sr)
            {
                EndPt p;
                p.pid = (gid + 1) * 100 + b->init;
                p.gid = gid;
                p[0] = i;
                p[1] = j;
                p[2] = k;
                b->particles.push_back(p);
                b->init++;
                // dprint("init %d: (%f %f %f)", p.pid, p[0], p[1], p[2]);
            }
        }
    }
}

void CInitSeeds(Block *b,
                int gid,
                const CBounds &cdomain,
                diy::RegularLink<CBounds> *l,
                float sr,
                int synth)
{

    // for synthetic data, seed only -x side of the block
    // float end = synth ? l->core().min[0] + 1.0: l->core().max[0];

    // std::vector<int> lastseeds ={65885, 606568, 1198076, 98640, 787319, 524793, 491894, 836084, 147814, 508280, 147919, 902008, 901752, 180578, 820088, 98764, 984059, 1201146, 590201, 688745};
    // std::vector<int> lastseeds ={508385, 147814, 131426, 886262, 458970, 754413, 1233786, 131531, 148039, 458971, 328020, 410068, 475355, 820335, 409843, 836853, 393684, 1017849, 164077, 540896};
    // std::vector<int> lastseeds ={721259, 542821, 82754, 98882, 66368, 98780, 148044, 426200, 115038, 131402, 987372, 1066482, 672596, 1331707, 1330686, 1098613, 65872, 17185, 1332091, 1329663};
    // std::vector<int> lastseeds ={180559, 1065980, 1198076, 65880, 836853, 131413, 164302, 82262, 66008, 426353, 393818, 147796, 147800, 115029, 196953, 131414, 262365, 278881, 98644, 213332};
    std::vector<int> lastseeds ={82254, 409970, 1198076, 393684, 1118587, 540896, 1366253, 17185, 1366377, 1167102, 622953, 1000701, 902008, 886262, 1035510, 721259, 1299188, 475355, 164200, 869627};

    // seed the block
    for (float i = float(cdomain.min[0]); i < cdomain.max[0]; i += sr)
    {
        for (float j = cdomain.min[1]; j < cdomain.max[1]; j += sr)
        {
            for (float k = cdomain.min[2]; k < cdomain.max[2]; k += sr)
            {
                EndPt p;
                p.pid = b->init;
                p.gid = gid;
                p[0] = i; p.pt_home.coords[0] = i;
                p[1] = j; p.pt_home.coords[1] = j;
                p[2] = k; p.pt_home.coords[2] = k;
                
                for (size_t ii=0; ii<lastseeds.size(); ii++){
                    if (lastseeds[ii] == p.pid){
                        b->particles.push_back(p); 
                        break;
                    }
                }
                    
                b->init++;
                // dprint("init %d: (%f %f %f)", p.pid, p[0], p[1], p[2]);
            }
        }
    }
}

// common to both exchange and iexchange
void trace_particles(Block *b,
                     const diy::Master::ProxyWithLink &cp,
                     const Decomposer &decomposer,
                     const int max_steps,
                     map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                     size_t &nsteps)
{
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds> *>(cp.link());

    const float *vec[3] = {b->vel[0], // shallow pointer copy
                           b->vel[1],
                           b->vel[2]};
    const int st[3] = {l->bounds().min[0],
                       l->bounds().min[1],
                       l->bounds().min[2]};
    const int sz[3] = {l->bounds().max[0] - l->bounds().min[0] + 1,
                       l->bounds().max[1] - l->bounds().min[1] + 1,
                       l->bounds().max[2] - l->bounds().min[2] + 1};

    for (auto i = 0; i < b->particles.size(); i++)
    {
        Pt &cur_p = b->particles[i].pt; // current end point
        Segment s(b->particles[i]);     // segment with one point p
        Pt next_p;                      // coordinates of next end point
        bool finished = false;
        if (b->particles[i].pid == 800)
            dprint("%f %f %f @ %d", cur_p.coords[0], cur_p.coords[1], cur_p.coords[2], cp.gid());

        // trace this segment until it leaves the block
        while (advect_rk1(st, sz, vec, cur_p.coords.data(), 0.05, next_p.coords.data()))
        {
            nsteps++;
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

        if (finished)
        { // this segment is done
            b->done++;
            //  dprint("pid: %d, step %d, %f %f %f", b->particles[i].pid, b->particles[i].nsteps, cur_p.coords[0], cur_p.coords[1], cur_p.coords[2]);
        }
        else // find destination of segment endpoint
        {
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int>> insert_it(dests, it);

            utl::in(*l, next_p.coords, insert_it, decomposer.domain, 1);

            EndPt out_pt(s);
            // out_pt.pid = b->particles[i].pid;
            out_pt.nsteps = b->particles[i].nsteps;
            if (dests.size())
            {
                diy::BlockID bid = l->target(dests[0]); // in case of multiple dests, send to first dest only

                // debug
                // fmt::print(stderr, "{}: gid {} enq to gid {}, steps {}, ({}, {}, {})\n", out_pt.pid, cp.gid(), bid.gid,  out_pt.nsteps, out_pt.pt.coords[0], out_pt.pt.coords[1], out_pt.pt.coords[2]);

                if (IEXCHANGE) // enqueuing single endpoint allows fine-grain iexchange if desired
                    cp.enqueue(bid, out_pt);
                else
                    outgoing_endpts[bid].push_back(out_pt); // vector of endpoints
            }
        }
    }
}

void deq_incoming_iexchange(Block *b,
                            const diy::Master::ProxyWithLink &cp, 
                            bool prediction)
{
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds> *>(cp.link());
    for (size_t i = 0; i < l->size(); ++i)
    {
        int nbr_gid = l->target(i).gid;
        // int cnt = 0;
        while (cp.incoming(nbr_gid))
        {   //cnt ++;
            EndPt incoming_endpt;
            cp.dequeue(nbr_gid, incoming_endpt);
            // b->particles.push_back(incoming_endpt);
            // if (prediction)
            //     b->particles.push_back(incoming_endpt);
            // else
                b->particles_pqueue.push(incoming_endpt);
        }
        // if (cnt>0)
        //  dprint("psize %d, gid %d", cnt, cp.gid());
    }
   
}

// Only for iexchange
bool trace_particles_iex(Block *b,
                         const diy::Master::ProxyWithLink &cp,
                         const CBounds &cdomain,
                         const int max_steps,
                         map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                         size_t &nsteps,
                         size_t &ntransfers,
                         bool prediction,
                         double &time_trace, 
                         const double &time_start, 
                         std::vector<int> &step_vs_time, 
                         std::vector<int> &csteps_vs_ftime, 
                         std::vector<int> &fhops_vs_ftime, 
                         size_t &np_core, 
                         std::vector<int> &esteps, 
                         std::vector<int> &final_ids, 
                         std::vector<double> &final_ids_time)
{
    diy::RegularLink<CBounds> *l = static_cast<diy::RegularLink<CBounds> *>(cp.link());

    const float *vec[3] = {b->vel[0], // shallow pointer copy
                           b->vel[1],
                           b->vel[2]};
    const float st[3] = {l->bounds().min[0],
                         l->bounds().min[1],
                         l->bounds().min[2]};
    const float sz[3] = {l->bounds().max[0] - l->bounds().min[0],
                         l->bounds().max[1] - l->bounds().min[1],
                         l->bounds().max[2] - l->bounds().min[2]};

    
    //for (auto i = 0; i < b->particles.size(); i++)

    np_core = b->particles.size() + b->particles_pqueue.size(); 
    // else 
    //     np_core = -1*(b->particles.size() + b->particles_pqueue.size()); 
    
     
    std::unique_ptr<EndPt> par (nullptr);

    if (b->particles.size()>0){
        // Pt &cur_p = b->particles[i].pt; // current end point
        // EndPt par = b->particles.back();
        par = std::unique_ptr<EndPt>(new EndPt(b->particles.back()));
        b->particles.pop_back();
    }else if(b->particles_pqueue.size()>0){
        par = std::unique_ptr<EndPt>(new EndPt(b->particles_pqueue.top()));
        b->particles_pqueue.pop();
        
    }
    
    if (par)
    {   

        Pt &cur_p = par->pt;

        // for vis analytics
        double time_pstart = MPI_Wtime();
        int bin = (int) std::floor(time_pstart - time_start);
        step_vs_time[bin*max_steps/32+par->nsteps/32] += 1;
        int old_step_cnt = par->nsteps;


        Segment s(par.get());     // segment with one point p
        Pt next_p;                      // coordinates of next end point
        bool finished = false;

        // trace this segment until it leaves the block
        double time_start_loc = MPI_Wtime();

        while (cadvect_rk1(st, sz, vec, cur_p.coords.data(), 0.05, next_p.coords.data()))
        {
            nsteps++;

            par->nsteps++;
            s.pts.push_back(next_p);
            cur_p = next_p;

            // if (par.pid == 800)
            //     dprint("%f %f %f @ %d", cur_p.coords[0], cur_p.coords[1], cur_p.coords[2], cp.gid());

            if (par->nsteps >= max_steps)
            {
                finished = true;
                break;
            }

            // // if predicting, add copy coordinates to EndPt and add to b->particles_store
            // if (prediction == true && cinside(cur_p, cdomain) && nsteps % 2 == 0)
            // {   
            //     EndPt way_pt;
            //     way_pt[0] = cur_p.coords[0];
            //     way_pt[1] = cur_p.coords[1];
            //     way_pt[2] = cur_p.coords[2];
            //     way_pt.predonly = true;
            //     b->particles_store.push_back(way_pt);
            // }

        }

        time_trace += MPI_Wtime() - time_start_loc;

        b->segments.push_back(s);

        if (!cinside(next_p, cdomain))
            finished = true;


        double time_pend = MPI_Wtime();
        int binf = (int) std::floor(time_pend - time_start);
        int new_step_cnt = par->nsteps - old_step_cnt;
        double time_diff = time_pend - time_pstart;
        // int cidx = binf*max_steps/32+new_step_cnt/32;
        // csteps_vs_ftime[cidx] += 1; 

        if (finished)
        { // this segment is done
            b->done++;
            // dprint("pid: %d, step %d, %f %f %f", par.pid, par.nsteps, cur_p.coords[0], cur_p.coords[1], cur_p.coords[2]);
            //   if (prediction){
                EndPt start_pt;

                start_pt[0] = par->pt_home.coords[0];
                start_pt[1] = par->pt_home.coords[1];
                start_pt[2] = par->pt_home.coords[2];
                start_pt.predonly = true;
                start_pt.esteps = par->nsteps;

                esteps[par->pid] = par->nsteps;

                // dprint("storing %f %f %f , esteps %d", start_pt[0], start_pt[1], start_pt[2], start_pt.esteps);
                // b->particles_store.push_back(start_pt);

            // }

            int cidx = binf*max_steps/32+par->nsteps/32;
            csteps_vs_ftime[cidx] += 1; 

            cidx = binf*max_steps/8+par->nhops/8;
            fhops_vs_ftime[cidx] += 1;


            if (final_ids.size()<20){
                final_ids.push_back(par->pid);
                final_ids_time.push_back(time_pend);
            }else{
                final_ids[b->done%20] = par->pid;
                final_ids_time[b->done%20] = time_pend;
            }


        }
        else // find destination of segment endpoint
        {   
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int>> insert_it(dests, it);

            // utl::in(*l, next_p.coords, insert_it, cdomain, 1);
            utl::in(*l, next_p.coords, insert_it, cdomain, false);

            EndPt out_pt(s);
            if (prediction){
                out_pt.pt_home.coords[0] = par->pt_home.coords[0]; 
                out_pt.pt_home.coords[1] = par->pt_home.coords[1]; 
                out_pt.pt_home.coords[2] = par->pt_home.coords[2]; 
            }
            // out_pt.pid = par.pid;
            out_pt.nsteps = par->nsteps;
            out_pt.nhops  = par->nhops + 1;  
            if (dests.size())
            {   
                ntransfers ++;
                diy::BlockID bid = l->target(dests[0]); // in case of multiple dests, send to first dest only

                // debug
                // fmt::print(stderr, "{}: gid {} enq to gid {}, steps {}, ({}, {}, {})\n", out_pt.pid, cp.gid(), bid.gid,  out_pt.nsteps, out_pt.pt.coords[0], out_pt.pt.coords[1], out_pt.pt.coords[2]);

                if (IEXCHANGE) // enqueuing single endpoint allows fine-grain iexchange if desired
                    cp.enqueue(bid, out_pt);
                else
                    outgoing_endpts[bid].push_back(out_pt); // vector of endpoints
            }
        }
        // dprint("BREAKING HERE");
        // break;

        // deq_incoming_iexchange(b, cp);
    }

    if (prediction==false)
        b->particles_store.clear();

    if (b->particles.size()>0  || b->particles_pqueue.size()>0){
        np_core = b->particles.size() + b->particles_pqueue.size(); 
        return false;
    }
    else{
        np_core = 0;
        return true;
    }
}

void deq_incoming_exchange(Block *b,
                           const diy::Master::ProxyWithLink &cp)
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



// common to both exchange and iexchange
void trace_block(Block *b,
                 const diy::Master::ProxyWithLink &cp,
                 const Decomposer &decomposer,
                 const diy::Assigner &assigner,
                 const int max_steps,
                 const float seed_rate,
                 const Decomposer::BoolVector share_face,
                 bool synth,
                 map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                 size_t &nsteps)
{
    const int gid = cp.gid();
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds> *>(cp.link());
    b->particles.clear();

    // const int   st[3]   = {l->core().min[0],
    //                        l->core().min[1],
    //                        l->core().min[2]};
    // const int   sz[3]   = {l->core().max[0] - l->core().min[0] + 1,
    //                        l->core().max[1] - l->core().min[1] + 1,
    //                        l->core().max[2] - l->core().min[2] + 1};
    // initialize seed particles first time
    if (b->init == 0)
        InitSeeds(b, gid, decomposer, l, seed_rate, synth);

    // dequeue incoming points and trace particles
    if (IEXCHANGE)
    {
        do
        {
            deq_incoming_iexchange(b, cp, 0);
            trace_particles(b, cp, decomposer, max_steps, outgoing_endpts, nsteps);
            b->particles.clear();
        } while (cp.fill_incoming());
    }
    else
    {
        deq_incoming_exchange(b, cp);
        trace_particles(b, cp, decomposer, max_steps, outgoing_endpts, nsteps);
    }
}

// will be called by trace_block_iexchange
bool trace_block_iex(Block *b,
                     const diy::Master::ProxyWithLink &cp,
                     const CBounds &cdomain,
                     const diy::Assigner &assigner,
                     const int max_steps,
                     const float seed_rate,
                     const Decomposer::BoolVector share_face,
                     bool synth,
                     map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                     size_t &nsteps,
                     size_t &ntransfers,
                     bool prediction,
                     double &time_trace, 
                     const double &time_start, 
                     std::vector<int> &step_vs_time, 
                      std::vector<int> &csteps_vs_ftime, 
                      std::vector<int> &fhops_vs_ftime, 
                      size_t &np_core, 
                      std::vector<int> &esteps, 
                      std::vector<int> &final_ids, 
                      std::vector<double> &final_ids_time)
{
    const int gid = cp.gid();
    diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds> *>(cp.link());
    // b->particles.clear();

    // const int   st[3]   = {l->core().min[0],
    //                        l->core().min[1],
    //                        l->core().min[2]};
    // const int   sz[3]   = {l->core().max[0] - l->core().min[0] + 1,
    //                        l->core().max[1] - l->core().min[1] + 1,
    //                        l->core().max[2] - l->core().min[2] + 1};

    // // initialize seed particles first time
    // if (b->init == 0)
    //     InitSeeds(b, gid, decomposer, l, seed_rate, synth);
    // dequeue incoming points and trace particles
    bool val = true;
    if (IEXCHANGE)
    {
        // do
        // {
            deq_incoming_iexchange(b, cp, prediction);
            // trace_particles(b, cp, decomposer, max_steps, outgoing_endpts, nsteps);
            val = trace_particles_iex(b, cp, cdomain, max_steps, outgoing_endpts, nsteps, ntransfers, prediction, time_trace, time_start, step_vs_time, csteps_vs_ftime, fhops_vs_ftime, np_core, esteps, final_ids, final_ids_time);
            // b->particles.clear();
        // } while (cp.fill_incoming());
    }
    return val;
}

void trace_block_exchange(Block *b,
                          const diy::Master::ProxyWithLink &cp,
                          const Decomposer &decomposer,
                          const diy::Assigner &assigner,
                          const int max_steps,
                          const float seed_rate,
                          const Decomposer::BoolVector share_face,
                          bool synth,
                          size_t &nsteps)
{
    map<diy::BlockID, vector<EndPt>> outgoing_endpts;

    trace_block(b, cp, decomposer, assigner, max_steps, seed_rate, share_face, synth, outgoing_endpts, nsteps);

    // enqueue the vectors of endpoints
    for (map<diy::BlockID, vector<EndPt>>::const_iterator it = outgoing_endpts.begin(); it != outgoing_endpts.end(); it++)
        cp.enqueue(it->first, it->second);

    // stage all_reduce of total initialized and total finished particle traces
    cp.all_reduce(b->particles.size(), plus<size_t>());
}

bool trace_block_iexchange(Block *b,
                           const diy::Master::ProxyWithLink &cp,
                           const CBounds &cdomain,
                           const diy::Assigner &assigner,
                           const int max_steps,
                           const float seed_rate,
                           const Decomposer::BoolVector share_face,
                           int synth,
                           size_t &nsteps,
                           size_t &ntransfers,
                           bool prediction, 
                           double &time_trace,
                           const double &time_start, 
                           std::vector<int> &step_vs_time, 
                           std::vector<int> &csteps_vs_ftime, 
                           std::vector<int> &fhops_vs_ftime, 
                           size_t &np_core, 
                           std::vector<int> &esteps, 
                           std::vector<int> &final_ids, 
                           std::vector<double> &final_ids_time)
{
    map<diy::BlockID, vector<EndPt>> outgoing_endpts; // needed to call trace_particles() but otherwise unused in iexchange
    // trace_block(b, cp, decomposer, assigner, max_steps, seed_rate, share_face, synth, outgoing_endpts, nsteps);
    bool val = trace_block_iex(b, cp, cdomain, assigner, max_steps, seed_rate, share_face, synth, outgoing_endpts, nsteps, ntransfers, prediction, time_trace, time_start, step_vs_time, csteps_vs_ftime, fhops_vs_ftime, np_core, esteps, final_ids, final_ids_time);
    return val;
}

// merge traces at the root block
void merge_traces(void *b_, const diy::ReduceProxy &rp, const diy::RegularMergePartners &)
{
    Block *b = static_cast<Block *>(b_);

    // dequeue and merge
    for (unsigned i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid()) // skip self
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
        int nbr_gid = rp.out_link().target(0).gid; // for a merge, the out_link size is 1; ie, there is only one target
        if (nbr_gid != rp.gid())                   // skip self
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
    void *callstack[128];
    int frames = backtrace(callstack, 128);
    char **strs = backtrace_symbols(callstack, frames);

    size_t funcnamesize = 256;
    char *funcname = (char *)malloc(funcnamesize);

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
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply __cxa_demangle():

            int status;
            char *ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
            if (status == 0)
            {
                funcname = ret; // use possibly realloc()-ed string
                fmt::print(stderr, "  {} : {}+{}", strs[i], funcname, begin_offset);
            }
            else
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

    signal(signum, SIG_DFL); // restore the default signal
    if (abort_on_segfault_)
        MPI_Abort(MPI_COMM_WORLD, 1);
}

// update incremental stats
// ref: Knuth "The Art of Computer Programming, Volume 2: Seminumerical Algorithms", section 4.2.2
// originally B.P. Welford, Technometrics, 4,(1962), 419-420
void update_stats(
    int trial,
    double time_start,
    int ncalls,
    const diy::mpi::communicator &world,
    Stats &stats)
{
    double cur_time = MPI_Wtime() - time_start;
    int cur_ncalls = 0;
    MPI_Reduce(&ncalls, &cur_ncalls, 1, MPI_INT, MPI_SUM, 0, world);

    if (trial == 0)
    {
        stats.cur_mean_time = cur_time;
        stats.prev_mean_time = cur_time;
        stats.cur_mean_ncalls = cur_ncalls;
        stats.prev_mean_ncalls = cur_ncalls;
        stats.cur_mean_callback_time = stats.cur_callback_time;
        stats.prev_mean_callback_time = stats.cur_callback_time;
        stats.cur_std_time = 0.0;
        stats.cur_std_ncalls = 0.0;
    }
    else
    {
        stats.cur_mean_time = stats.prev_mean_time +
                              (cur_time - stats.prev_mean_time) / (trial + 1);
        stats.cur_mean_ncalls = stats.prev_mean_ncalls +
                                (cur_ncalls - stats.prev_mean_ncalls) / (trial + 1);
        stats.cur_mean_callback_time = stats.prev_mean_callback_time +
                                       (stats.cur_callback_time - stats.prev_mean_callback_time) / (trial + 1);
        stats.cur_std_time = stats.prev_std_time +
                             (cur_time - stats.prev_mean_time) * (cur_time - stats.cur_mean_time);
        stats.cur_std_ncalls = stats.prev_std_ncalls +
                               (cur_ncalls - stats.prev_mean_ncalls) * (cur_ncalls - stats.cur_mean_ncalls);
    }
    stats.prev_mean_time = stats.cur_mean_time;
    stats.prev_mean_ncalls = stats.cur_mean_ncalls;
    stats.prev_mean_callback_time = stats.cur_mean_callback_time;
    stats.prev_std_time = stats.cur_std_time;
    stats.prev_std_ncalls = stats.cur_std_ncalls;

    // debug
    //     if (world.rank() == 0)
    //     {
    //         fmt::print(stderr, "trial {} time {} callback time {} ncalls {}\n",
    //                 trial, cur_time, stats.cur_callback_time, cur_ncalls);
    //     }
}

void print_results(
    float seed_rate,
    int nprocs,
    int nblocks,
    int tot_nsynth,
    int ntrials,
    int nrounds,
    const Stats &stats)
{
    fmt::print(stderr, "---------- stats ----------\n");
    if (IEXCHANGE)
        fmt::print(stderr, "using iexchange\n");
    else
        fmt::print(stderr, "using exchange\n");
    fmt::print(stderr, "seed rate:                       {}\n", seed_rate);
    fmt::print(stderr, "nprocs:                          {}\n", nprocs);
    fmt::print(stderr, "nblocks:                         {}\n", nblocks);
    fmt::print(stderr, "nsynth:                          {}\n", tot_nsynth);
    fmt::print(stderr, "ntrials:                         {}\n", ntrials);
    fmt::print(stderr, "mean time (s):                   {}\n", stats.cur_mean_time);
    fmt::print(stderr, "std dev time (s):                {}\n", ntrials > 1 ? sqrt(stats.cur_std_time / (ntrials - 1)) : 0.0);
    if (IEXCHANGE)
    {
        fprintf(stderr, "mean # callbacks:                %.0lf\n", stats.cur_mean_ncalls);
        fprintf(stderr, "std dev # callbacks:             %.0lf\n", ntrials > 1 ? sqrt(stats.cur_std_ncalls / (ntrials - 1)) : 0.0);
    }
    else
    {
        fmt::print(stderr, "# rounds:                        {}\n", nrounds);
        fmt::print(stderr, "mean callback (advect) time (s): {}\n", stats.cur_mean_callback_time);
    }
    fmt::print(stderr, "---------------------------\n");
}

void print_exceeded_max_rounds(diy::Master &master)
{
    if (master.communicator().rank() == 0)
        fmt::print(stderr, "*** Warning: max # rounds for exchange has been reached. ***\n");

    // debug: print unterminated particles
    master.foreach ([](Block *b, const diy::Master::ProxyWithLink &cp) {
        if (b->particles.size() > 0)
        {
            fmt::print(stderr, "gid = {}, particles size = {}\n", cp.gid(), b->particles.size());
            auto *l = static_cast<RGLink *>(cp.link());
            fmt::print(stderr, "  core = {} - {}, bounds = {} - {}\n",
                       l->core().min, l->core().max,
                       l->bounds().min, l->bounds().max);
            for (auto &p : b->particles)
                fmt::print(stderr, "  {}\n", p.pt.coords);
        }
    });
}

#ifdef WITH_VTK

void render_traces(
    diy::Master &master,
    diy::Assigner &assigner,
    Decomposer &decomposer,
    bool merge) // merge traces to one block
{
    if (merge)
    {
        // merge-reduce traces to one block
        int k = 2; // the radix of the k-ary reduction tree
        diy::RegularMergePartners partners(decomposer, k);
        diy::reduce(master, assigner, partners, &merge_traces);

        if (master.communicator().rank() == 0)
        {
            fprintf(stderr, "converting particle traces to vtk polylines and rendering 1 block only\n");
            ((Block *)master.block(0))->render();
        }
    }
    else
    {
        if (master.communicator().rank() == 0)
            fprintf(stderr, "converting particle traces to vtk polylines and rendering all blocks\n");
        master.foreach (&Block::render_block);
    }
}

#endif

void write_traces(
    diy::Master &master,
    diy::Assigner &assigner,
    Decomposer &decomposer)
{
    // merge-reduce traces to one block
    int k = 2; // the radix of the k-ary reduction tree
    diy::RegularMergePartners partners(decomposer, k);
    diy::reduce(master, assigner, partners, &merge_traces);

    if (master.communicator().rank() == 0)
    {
        fprintf(stderr, "Check is turned on: merging traces to one block and writing them to disk\n");
        std::string filename;
        if (IEXCHANGE)
            filename = "iexchange.txt";
        else
            filename = "exchange.txt";
        ((Block *)master.block(0))->write_segments(filename);
    }
}

void output_profile(
    diy::Master &master,
    int nblocks)
{
    if (IEXCHANGE)
    {
        diy::io::SharedOutFile prof_out(fmt::format("profile-iexchange-p{}-b{}.txt",
                                                    master.communicator().size(), nblocks),
                                        master.communicator());
        master.prof.output(prof_out, std::to_string(master.communicator().rank()));
        prof_out.close();
    }
    else
    {
        diy::io::SharedOutFile prof_out(fmt::format("profile-exchange-p{}-b{}.txt",
                                                    master.communicator().size(), nblocks),
                                        master.communicator());
        master.prof.output(prof_out, std::to_string(master.communicator().rank()));
        prof_out.close();
    }
}

void pairsort(std::vector<double> &a, std::vector<int> &b, int n) 
{ 
    pair<double, int> pairt[n]; 
  
    // Storing the respective array 
    // elements in pairs. 
    for (int i = 0; i < n; i++)  
    { 
        pairt[i].first = a[i]; 
        pairt[i].second = b[i]; 
    } 
  
    // Sorting the pair array. 
    sort(pairt, pairt + n); 
      
    // Modifying original arrays 
    for (int i = 0; i < n; i++)  
    { 
        a[i] = pairt[i].first; 
        b[i] = pairt[i].second; 
    } 
} 

int main(int argc, char **argv)
{
    signal(SIGSEGV, catch_sig); // catch segfault

    string infile;    // input file name
    Bounds domain{3}; // global domain bounds
    CBounds cdomain{3};
    int max_steps;   // max number of steps a particle is allowed to take
    float seed_rate; // seed particle every this many grid pts in each dim

    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    using namespace opts;

    // defaults
    int nblocks = world.size();     // total number of global blocks
    int nthreads = 2;               // number of threads diy can use
    int mblocks = -1;               // number of blocks in memory (-1 = all)
    string prefix = "./DIY.XXXXXX"; // storage of temp files
    int ndims = 3;                  // domain dimensions
    float vec_scale = 1.0;          // vector field scaling factor
    int hdr_bytes = 0;              // num bytes header before start of data in infile
    int max_rounds = 0;             // max number of rounds to trace (0 = no limit)
    size_t min_queue_size = 0;      // min queue size (bytes) for iexchange
    size_t max_hold_time = 0;       // max hold time (microsec) for iexchange
    int synth = 0;                  // generate various synthetic input datasets
    float slow_vel = 1.0;           // slow velocity for synthetic data
    float fast_vel = 10.0;          // fast velocity for synthetic data
    int check = 0;                  // write out traces for checking
    std::string log_level = "info"; // logging level
    int ntrials = 1;                // number of trials
    bool merged_traces = false;     // traces have already been merged to one block
    int tot_nsynth = nblocks;       // total number of synthetic slow velocity regions
    bool prediction = true;

    double time_total = 0, time_overhead = 0;
    double time_prep=0, time_predrun=0, time_kdtree=0, time_readdata=0, time_filter=0, time_final=0, time_predrun_loc=0, time_final_loc, time_trace = 0;

    // command-line ags
    Options ops(argc, argv);
    ops >> Option('b', "blocks", nblocks, "Total number of blocks to use") >> Option('t', "threads", nthreads, "Number of threads to use") >> Option('m', "in-memory", mblocks, "Number of blocks to keep in memory") >> Option('s', "storage", prefix, "Path for out-of-core storage") >> Option('v', "vec-scale", vec_scale, "Vector field scaling factor") >> Option('h', "hdr-bytes", hdr_bytes, "Skip this number bytes header in infile") >> Option('r', "max-rounds", max_rounds, "Max number of rounds to trace") >> Option('x', "synthetic", synth, "Generate various synthetic flows") >> Option('w', "slow-vel", slow_vel, "Slow velocity for synthetic data") >> Option('f', "fast-vel", fast_vel, "Fast velocity for synthetic data") >> Option('c', "check", check, "Write out traces for checking") >> Option('l', "log", log_level, "log level") >> Option('n', "trials", ntrials, "number of trials") >> Option('o', "nsynth", tot_nsynth, "total number of synthetic velocity regions");
    bool fine = ops >> Present("fine", "Use fine-grain icommunicate");

    if (ops >> Present('h', "help", "show help") ||
        !(ops >> PosOption(infile) >> PosOption(max_steps) >> PosOption(seed_rate) >> PosOption(domain.min[0]) >> PosOption(domain.min[1]) >> PosOption(domain.min[2]) >> PosOption(domain.max[0]) >> PosOption(domain.max[1]) >> PosOption(domain.max[2]) >> PosOption(prediction)))
    {
        if (world.rank() == 0)
        {
            fprintf(stderr, "Usage: %s [OPTIONS] infile mins maxs\n", argv[0]);
            cout << ops;
        }
        return 1;
    }

    cdomain.min[0] = domain.min[0];
    cdomain.min[1] = domain.min[1];
    cdomain.min[2] = domain.min[2];
    cdomain.max[0] = domain.max[0];
    cdomain.max[1] = domain.max[1];
    cdomain.max[2] = domain.max[2];

    //     diy::create_logger(log_level);
    diy::FileStorage storage(prefix);
    diy::Master master(world,
                       nthreads,
                       mblocks,
                       &Block::create,
                       &Block::destroy,
                       &storage,
                       &Block::save,
                       &Block::load);
    diy::RoundRobinAssigner assigner(world.size(), nblocks);

    // decompose domain
    Decomposer::BoolVector share_face;
    Decomposer::BoolVector wrap; // defaults to false
    Decomposer::CoordinateVector ghosts;
    ghosts.push_back(1);
    ghosts.push_back(1);
    ghosts.push_back(1);
    share_face.push_back(true);
    share_face.push_back(true);
    share_face.push_back(true);

    Decomposer decomposer(ndims,
                          domain,
                          assigner.nblocks(),
                          share_face,
                          wrap,
                          ghosts);
    if (synth == 1)
    {
        AddConsistentSynthetic addsynth(master, slow_vel, fast_vel, tot_nsynth);
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

    Stats stats; // incremental stats, default initialized to 0's
    int nrounds;
    size_t nsteps = 0, ntransfers = 0, nsteps_lagged = 0;
    std::vector<size_t> steps_per_interval;
    std::atomic<bool> done{false};
    std::mutex mutex;
    std::vector<int> step_vs_time((max_steps/32) * 400);
    std::vector<int> csteps_vs_ftime(max_steps/32 * 400);
    std::vector<int> fhops_vs_ftime(max_steps/8 * 200);
    std::vector<size_t> particles_in_core;
    size_t np_core = 0;
    std::vector<int> esteps((512/seed_rate)*(512/seed_rate)*(512/seed_rate));
    std::vector<int> esteps_all(esteps.size());
    int ndone=0;
    std::vector<int> done_counts; // for tracking number of particles done each second
    std::vector<int> final_ids; // ids of particles finishing in the end
    std::vector<double> final_ids_time; // finishting times of particles finishing in the end

    // check if clocks are synchronized by printing the value of MPI_WTIME_IS_GLOBAL and timing an initial barrier
    // barrier also has the effect of removing any skew in generating or reading the data
    // all ranks starting synchronized makes the performance profiles easier to understand
    int flag, *get_val;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_WTIME_IS_GLOBAL, &get_val, &flag);
    if (world.rank() == 0)
        fmt::print(stderr, "MPI_WTIME_IS_GLOBAL = {} flag = {}\n", *get_val, flag);
    master.prof << "initial barrier";
    world.barrier();
    master.prof >> "initial barrier";


    diy::Master master_iex(world,
                                   nthreads,
                                   mblocks,
                                   &Block::create,
                                   &Block::destroy,
                                   &storage,
                                   &Block::save,
                                   &Block::load);

    // run the trials
    for (int trial = 0; trial < ntrials; trial++)
    {
        int ncalls = 0;

        // debug
        if (world.rank() == 0)
            fprintf(stderr, "started particle tracing trial %d\n", trial);

        // reset the block particle traces, but leave the vector field intact
        master.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
            b->init = 0;
            b->done = 0;
            b->segments.clear();
            b->particles.clear();
        });

        world.barrier();
        double time_start = MPI_Wtime();

         

        if (IEXCHANGE)
        {

           
            diy::ContiguousAssigner cassigner(world.size(), nblocks);

            CDecomposer cdecomposer(ndims,
                                    cdomain,
                                    cassigner.nblocks(),
                                    share_face,
                                    wrap);

            CAddAndRead caddblock(master_iex, infile.c_str(), world, vec_scale, hdr_bytes);
            cdecomposer.decompose(world.rank(), cassigner, caddblock);

            //   std::vector<int> gids;
            //   cassigner.local_gids(world.rank(), gids);
            //   for (unsigned i = 0; i < gids.size(); ++i)
            //   {
            //      int             gid = gids[i];
            //      CAddAndRead caddblock(master_iex, infile.c_str(), world, vec_scale, hdr_bytes, cdomain, gid);
            //   }

            master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                b->init = 0;
                b->done = 0;
                b->segments.clear();
                b->particles.clear();
            });

            master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                const int gid = cp.gid();

                diy::RegularLink<CBounds> *l = static_cast<diy::RegularLink<CBounds> *>(cp.link());
                if (b->init == 0 && gid == 0)
                    CInitSeeds(b, gid, cdomain, l, seed_rate, synth);
            });

            diy::kdtree(master_iex, cassigner, ndims, cdomain, &Block::particles, 2 * 512, false);

            bool verbose = false;
            master_iex.foreach ([verbose](Block *b, const diy::Master::ProxyWithLink &cp) { print_block(b, cp, verbose); });

            master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                int gid = cp.gid();
                RCLink *l = static_cast<RCLink *>(cp.link());

                caddblock.read_data(b, l->bounds(), gid);
            });

            // for gantt chart
            std::thread steps_ctr([&]{

                while(!done){

                        std::this_thread::sleep_for( std::chrono::milliseconds(1000));
                        {
                             std::lock_guard<std::mutex> guard(mutex);

                             particles_in_core.push_back(np_core);
                             if (steps_per_interval.size()>0){
                                steps_per_interval.push_back(nsteps - nsteps_lagged);
                                nsteps_lagged = nsteps;
                                
                             }
                             else
                                steps_per_interval.push_back(nsteps);

                        }
                }
                
            });


            
            if (prediction)
            {
                world.barrier();
                double time0 = MPI_Wtime(); 

                // sample prediction points
                master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                    // std::random_device rd;
                    // std::mt19937 g(4);
                    // std::shuffle(b->particles.begin(), b->particles.end(), g);

                    size_t pred_size = b->particles.size();
                    for (size_t i = 0; i < pred_size; i++)
                    {
                        b->particles_store.push_back(b->particles[i]);
                    }

                    // // init store unpredicted particles?; particles_store used to tag them along 
                    // b->particles_store.insert(std::end(b->particles_store), std::begin(b->particles) + pred_size, std::end(b->particles));
                    // b->particles.resize(pred_size);
                });

               

                // master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &icp) -> bool {
                //     dprint("psizes %ld %ld", b->particles.size(), b->particles_store.size());
                // });
                //   prediction run
                world.barrier();
                double time1 = MPI_Wtime();
                time_prep = time1 - time0;

                // advect prediction subset
                master_iex.iexchange([&](Block *b, const diy::Master::ProxyWithLink &icp) -> bool {
                    ncalls++;

                    bool val = trace_block_iexchange(b,
                                                     icp,
                                                     cdomain,
                                                     cassigner,
                                                     max_steps,
                                                     seed_rate,
                                                     share_face,
                                                     synth,
                                                     nsteps,
                                                     ntransfers, 
                                                     true, 
                                                     time_trace, 
                                                     time_start, 
                                                     step_vs_time, 
                                                     csteps_vs_ftime, 
                                                     fhops_vs_ftime,
                                                     np_core, 
                                                     esteps, 
                                                     final_ids, 
                                                     final_ids_time);
                    return val;
                });
                np_core = -1;

                time_predrun_loc = MPI_Wtime() - time1;
                world.barrier();
                double time2 = MPI_Wtime(); 
                time_predrun = time2 - time1;

                // master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                //     b->particles.clear();
                // });

                // get esteps from all
                diy::mpi::reduce(world, esteps, esteps_all, 0, std::plus<int>());

                


                // swap b->particles and b->particles_store
                master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &icp) -> bool {
                    b->particles = std::move(b->particles_store);

                    for (size_t i=0; i<b->particles.size(); i++)
                        b->particles[i].esteps = esteps_all[i];
                });

                //rebalance
                diy::kdtree(master_iex, cassigner, ndims, cdomain, &Block::particles, 2 * 512, false);

                 world.barrier();
                 double time3 = MPI_Wtime();
                 time_kdtree = time3 - time2;


                verbose = false;
                master_iex.foreach ([verbose](Block *b, const diy::Master::ProxyWithLink &cp) { print_block(b, cp, verbose); });

                // read data again
                master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                    dprint("sizes %ld", b->particles.size());
                    int gid = cp.gid();
                    RCLink *l = static_cast<RCLink *>(cp.link());

                    caddblock.read_data(b, l->bounds(), gid);
                });

                world.barrier();
                double time4 = MPI_Wtime();
                time_readdata = time4 - time3;

                // // filter out non prediction particles
                // master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &icp) -> bool {
                //     b->particles_store.clear();

                //     int block_esteps=0; 
                //     int num_pred = 0; // number of predicted particles
                //     for (size_t i = 0; i < b->particles.size(); i++)
                //     {
                //         if (b->particles[i].predonly == false) // store the to be advected points
                //             b->particles_store.push_back(b->particles[i]);
                //         // else store the expected #steps of an arbitrary predicted particle
                //         else{
                            
                //                 block_esteps += b->particles[i].esteps;
                //         }
                //     }

                //     // block's expected prediction length is average of esteps of all its predicted points
                //     if (num_pred>0)
                //         block_esteps /= num_pred;

                //     for (size_t i=0; i<b->particles_store.size(); i++){
                //         b->particles_store[i].esteps = block_esteps;
                //     }

                //     b->particles = std::move(b->particles_store);
                // });

                world.barrier();
                double time5 = MPI_Wtime();
                time_filter = time5 - time4;

                // master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &icp) -> bool {
                //     dprint("psizes %ld %ld", b->particles.size(), b->particles_store.size());
                // });
            }
           
            
            world.barrier();
            double time6 = MPI_Wtime();


            // post prediction run
            master_iex.iexchange([&](Block *b, const diy::Master::ProxyWithLink &icp) -> bool {
                ncalls++;

                bool val = trace_block_iexchange(b,
                                                 icp,
                                                 cdomain,
                                                 cassigner,
                                                 max_steps,
                                                 seed_rate,
                                                 share_face,
                                                 synth,
                                                 nsteps,
                                                 ntransfers,
                                                 false, 
                                                 time_trace, 
                                                 time_start, 
                                                 step_vs_time, 
                                                 csteps_vs_ftime, 
                                                 fhops_vs_ftime,
                                                 np_core, esteps, 
                                                 final_ids, 
                                                 final_ids_time);
                return val;
            });

            np_core = -1;

            time_final_loc = MPI_Wtime() - time6;
            world.barrier();
            double time7 = MPI_Wtime();

            time_final = time7 - time6;


            done = true;
            steps_ctr.join();


        }
        else // exchange
        {
            // particle tracing for either a maximum number of rounds or, if max_rounds == 0,
            // then for inifinitely many rounds until breaking out when done is true
            int stop = (max_rounds ? max_rounds : 1);
            int incr = (max_rounds ? 1 : 0);

            nrounds = 0;
            stats.cur_callback_time = 0.0;

            world.barrier();
            double time1 = MPI_Wtime();
            for (int round = 0; round < stop; round += incr)
            {
                nrounds++;

                // advect
                double t0 = MPI_Wtime();
                master.foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
                    trace_block_exchange(b,
                                         cp,
                                         decomposer,
                                         assigner,
                                         max_steps,
                                         seed_rate,
                                         share_face,
                                         synth,
                                         nsteps);
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
            } // rounds

            world.barrier();
            double time2 = MPI_Wtime();

            time_total = time2 - time1;

            // debug: exceeded number of rounds for exchange
            if (nrounds == max_rounds)
                print_exceeded_max_rounds(master);
        }


        

        size_t nsteps_global=0;
        diy::mpi::reduce(world, nsteps, nsteps_global, 0, std::plus<size_t>());
        size_t maxsteps_global=0;
        diy::mpi::reduce(world, nsteps, maxsteps_global, 0, diy::mpi::maximum<size_t>());
        size_t minsteps_global;
        diy::mpi::reduce(world, nsteps, minsteps_global, 0, diy::mpi::minimum<size_t>());
        size_t ntransfers_global;
        diy::mpi::reduce(world, ntransfers, ntransfers_global, 0, std::plus<size_t>());

        double time_fin_loc_max;
        diy::mpi::reduce(world, time_final_loc, time_fin_loc_max, 0, diy::mpi::maximum<double>());
        double time_fin_loc_avg;
        diy::mpi::reduce(world, time_final_loc, time_fin_loc_avg, 0, std::plus<double>());
        time_fin_loc_avg = time_fin_loc_avg/world.size(); 

        double time_predrun_loc_max;
        diy::mpi::reduce(world, time_predrun_loc, time_predrun_loc_max, 0, diy::mpi::maximum<double>());
        double time_predrun_loc_avg;
        diy::mpi::reduce(world, time_predrun_loc, time_predrun_loc_avg, 0, std::plus<double>());
        time_predrun_loc_avg = time_predrun_loc_avg/world.size();


        double time_trace_max;
        diy::mpi::reduce(world, time_trace, time_trace_max, 0, diy::mpi::maximum<double>()); 

        double time_trace_avg;
        diy::mpi::reduce(world, time_trace, time_trace_avg, 0, std::plus<double>());
        time_trace_avg = time_trace_avg/world.size();

        std::vector<std::vector<size_t>> all_steps_per_interval;
        diy::mpi::gather (world, steps_per_interval, all_steps_per_interval, 0 );

        std::vector<int> step_vs_time_all(step_vs_time.size());
        diy::mpi::reduce(world, step_vs_time, step_vs_time_all, 0, std::plus<int>());

        std::vector<int> csteps_vs_ftime_all(csteps_vs_ftime.size());
        diy::mpi::reduce(world, csteps_vs_ftime, csteps_vs_ftime_all, 0, std::plus<int>());
    
        std::vector<size_t> particles_in_core_all(particles_in_core.size());
        diy::mpi::reduce(world, particles_in_core, particles_in_core_all, 0, std::plus<size_t>());

        std::vector<int> fhops_vs_ftime_all (fhops_vs_ftime.size());
        diy::mpi::reduce(world, fhops_vs_ftime, fhops_vs_ftime_all, 0, std::plus<int>());
        
        std::vector<std::vector<int>> final_ids_all;
        diy::mpi::gather(world, final_ids, final_ids_all, 0);


        std::vector<std::vector<double>> final_ids_time_all;
        diy::mpi::gather(world, final_ids_time, final_ids_time_all, 0);


        float avg = float(nsteps_global) / world.size();
        // float balance = float(maxsteps_global) / avg;
        float balance = (float(maxsteps_global))/ float(avg);
        // debug
        if (world.rank() == 0)
        {
            fprintf(stderr, "finished particle tracing trial %d\n", trial);
            fprintf(stderr, "predd , %d, nsteps_global , %ld, maxsteps_global , %ld, bal , %f, time_tot , %f, time_overhead, %f, worldsize, %d,minsteps, %ld,\n", prediction, nsteps_global, maxsteps_global, balance, time_total, time_overhead, world.size(), minsteps_global);
            dprint("times: predrun, %f, kdtree , %f, readdata, %f, filter ,%f, final , %f, prediction, %d, max, %ld, min, %ld, nsteps, %ld, wsize, %d, time_pred, %f, tot_transfers, %ld, prdrun_local(max avg), %f, %f, fin_local (max avg), %f, %f, max_steps, %d, time_trace_max, %f, time_trace_avg, %f, ", time_predrun, time_kdtree, time_readdata, time_filter, time_final, prediction, maxsteps_global, minsteps_global, nsteps_global, world.size(), time_prep, ntransfers_global, time_predrun_loc_max, time_predrun_loc_avg, time_fin_loc_max, time_fin_loc_avg, max_steps, time_trace_max, time_trace_avg);

            for (size_t i=0; i<all_steps_per_interval.size(); i++)
                {   fprintf(stderr, "ganttrank, %d, p, %d, ws, %d, ", i, prediction, world.size());
                    for (size_t j=0; j<all_steps_per_interval[i].size(); j++){
                        fprintf(stderr, "%ld, ", all_steps_per_interval[i][j]);
                    }
                    fprintf(stderr, "\n");
                }  

            fprintf(stderr, "step_vs_time, p, %d, ws, %d, ", prediction, world.size());
            for (size_t i=0; i<step_vs_time_all.size(); i++){
                fprintf(stderr, "%d ", step_vs_time_all[i]);
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "csteps_vs_ftime, p, %d, ws, %d, ", prediction, world.size());
            for (size_t i=0; i<csteps_vs_ftime.size(); i++){
                fprintf(stderr, "%d ", csteps_vs_ftime_all[i]);
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "particles_in_core, p, %d, ws, %d, ", prediction, world.size());
            for (size_t i=0; i<particles_in_core_all.size(); i++){
                fprintf(stderr, "%ld ", particles_in_core_all[i]);
            }
            fprintf(stderr, "\n");


            fprintf(stderr, "fhops_vs_ftime, p, %d, ws, %d, ", prediction, world.size());
            for (size_t i=0; i<fhops_vs_ftime_all.size(); i++){
                fprintf(stderr, "%ld ", fhops_vs_ftime_all[i]);
            }
            fprintf(stderr, "\n");

            for (size_t i=1; i<final_ids_all.size(); i++){
                for (size_t j=0; j<final_ids_all[i].size(); j++){
                    final_ids.push_back(final_ids_all[i][j]);
                    final_ids_time.push_back(final_ids_time_all[i][j]);
                }

            }            
            pairsort(final_ids_time, final_ids, (int)final_ids.size());
            dprint("finalids %ld", final_ids_all.size());
            int numitems = 20;
            if (final_ids.size()>numitems){
                for (int i=0; i<numitems;i++){
                    dprint("lastpart, %d, %d", i, final_ids[final_ids.size()-1-i]);
                }
            }


            

        }

        //         master.prof.totals().output(std::cerr);

        update_stats(trial, time_start, ncalls, world, stats);

#ifdef WITH_VTK
        render_traces(master, assigner, decomposer, true);
#endif

#ifdef DIY_PROFILE
        output_profile(master, nblocks);
#endif

    } // number f trials

    if (world.rank() == 0)
        print_results(seed_rate, world.size(), nblocks, tot_nsynth, ntrials, nrounds, stats);

    // write trajectory segments for validation
    if (check)
        write_traces(master_iex, assigner, decomposer);

    // debug
    //     master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
    //     {
    //         b->show_geometry(cp);
    //     });

    return 0;
}
