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

#include "zoltan.h"

#include "btrace/partitioning.h"
#include "btrace/message.h"
#include "btrace/io.h"

#include "btrace/bblock.hpp"
#include "btrace/advection.h"

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

// dequeue remote data
// there is still a link, but exchange(remote = true) exchanged messages from any block
void remote_deq(BBlock* b, const diy::Master::ProxyWithLink& cp)
{
    std::vector<int> incoming_gids;
    cp.incoming(incoming_gids);
    for (size_t i = 0; i < incoming_gids.size(); i++)
        if (cp.incoming(incoming_gids[i]).size())
        {
            // int recvd_data;
            datablock recvd_data;
            cp.dequeue(incoming_gids[i], recvd_data);
            // fmt::print(stderr, "Remote dequeue: gid {} received value {} from gid {}\n", cp.gid(), recvd_data.data.size(), recvd_data.from_proc);
            for (size_t i=0; i<recvd_data.cgid.size(); i++){
                b->data[recvd_data.cgid[i]] = recvd_data.data[i];
                b->particles[recvd_data.cgid[i]] = std::move(recvd_data.particles[i]);
                b->bounds[recvd_data.cgid[i]] = std::move(recvd_data.bounds[i]);
            }   
        }

    update_mesh_data(b);

}

// merge traces at the root block
void merge_traces(void *b_, const diy::ReduceProxy &rp, const diy::RegularMergePartners &)
{
    BBlock *b = static_cast<BBlock *>(b_);

    // dequeue and merge
    for (unsigned i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid()) // skip self
            continue;

        vector<BSegment> in_traces;
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

void write_segments(std::string filename, std::vector<Segment> &segments)
    {
        ofstream f;
        f.open(filename);

        // debug
        fmt::print("writing {} segments\n", segments.size());

        for (size_t i = 0; i < segments.size(); i++)
        {
            for (size_t j = 0; j < segments[i].pts.size(); j++)
            {
                // debug
//                 fprintf(f, "%ld %f %f %f, ", segments[i].pts.size(), segments[i].pts[0].coords[0], segments[i].pts[0].coords[1], segments[i].pts[0].coords[2]);
                f << std::setprecision(8) << segments[i].pts[j].coords[0] << " " << segments[i].pts[j].coords[1] << " " << segments[i].pts[j].coords[2] << " ";
            }
            f << endl;
        }
        f.close();
    }

void write_traces(
    diy::Master &master,
    diy::Assigner &assigner,
    Decomposer &decomposer, 
    diy::mpi::communicator &world)
{
    dprint("merging traces");
    // merge-reduce traces to one block
    // int k = 2; // the radix of the k-ary reduction tree
    // diy::RegularMergePartners partners(decomposer, k);
    // diy::reduce(master, assigner, partners, &merge_traces);
    std::vector<float> segs;
  
    
    std::vector<std::vector<float>> segs2;
   


    std::vector<int> sizes;
    std::vector<std::vector<int>> all_sizes;


    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {

        // if (world.rank()==0)
            for (size_t i=0; i<b->segments.size(); i++){
                for (size_t j=0; j<b->segments[i].pts.size(); j++){
                    segs.push_back(b->segments[i].pts[j].coords[0]);
                    segs.push_back(b->segments[i].pts[j].coords[1]);
                    segs.push_back(b->segments[i].pts[j].coords[2]);

                    // dprint("segg %f %f %f| %ld", b->segments[i].pts[j].coords[0], b->segments[i].pts[j].coords[1], b->segments[i].pts[j].coords[2], b->segments[i].pts.size());
                }
                sizes.push_back(b->segments[i].pts.size());
                // dprint("%ld (%f %f %f)", i, b->segments[i].pts[0].coords[0], b->segments[i].pts[0].coords[1], b->segments[i].pts[0].coords[2]);
                // pvi(sizes);
            }

    });

    diy::mpi::gather(world, segs, segs2, 0);
    diy::mpi::gather(world, sizes, all_sizes, 0);

   
   


    if (master.communicator().rank() == 0)
    {

         dprint("segs2size %ld", segs2.size());
    
        std::vector<Segment> all_segs;
        
        for(size_t i=0; i<all_sizes.size(); i++){
            size_t idx=0;
            for (size_t j=0; j<all_sizes[i].size(); j++){
                Segment seg;
                
                
                for(int k=0; k<all_sizes[i][j]; k++){
                    // dprint("all_sizes[i][j] %d, idx %ld, %f %f %f", all_sizes[i][j], idx, segs2[i][idx*3], segs2[i][idx*3+1], segs2[i][idx*3+2]);
                    Pt p;
                    p.coords[0] = segs2[i][idx*3];
                    p.coords[1] = segs2[i][idx*3+1];
                    p.coords[2] = segs2[i][idx*3+2];
                    seg.pts.push_back(p);
                    idx++;
                    // dprint("op %f %f %f", p.coords[0], p.coords[1], p.coords[2]);
                }
                // dprint("i %d j %d sizes %ld, csize %ld ", i, j, all_sizes[i][j], segs2[i].size());
                all_segs.push_back(seg);

            }
            
        }

        fprintf(stderr, "Check is turned on: merging traces to one block and writing them to disk\n");
        std::string filename;
       
        filename = "baseline.txt";
        // ((Block *)master.block(0))->write_segments(filename);
        write_segments(filename, all_segs);
    }
}




int main(int argc, char **argv){

    signal(SIGSEGV, catch_sig); // catch segfault

    string infile;    // input file name
    Bounds domain{3}; // global domain bounds
    CBounds cdomain{3};
    int max_steps;   // max number of steps a particle is allowed to take
    float seed_rate; // seed particle every this many grid pts in each dim

    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    using namespace opts;

   


    // read in arguments

    // defaults
    int nblocks = world.size();     // total number of global blocks
    int nthreads = 1;               // number of threads diy can use
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

    diy::FileStorage storage(prefix);
    diy::Master master(world,
                       nthreads,
                       mblocks,
                       &BBlock::create,
                       &BBlock::destroy,
                       &storage,
                       &BBlock::save,
                       &BBlock::load);
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

    std::vector<int> gids;                     // global ids of local blocks
    assigner.local_gids(world.rank(), gids);   // get the gids of local blocks
    for (size_t i = 0; i < gids.size(); ++i)   // for the local blocks in this processor
    {
        int gid = gids[i];

        diy::Link*   link = new diy::Link;   // link is this block's neighborhood
        master.add(gid, new BBlock, link);    // add the current local block to the master
    }

    size_t nsteps = 0, init = 0, done = 0;

    // definitions
    // element: 1x1x1 pixel
    // cell: smallest unit that can be moved around (formerly block)
    // block: data consisting of multiple cells in a process

    int C[3] = {4, 4, 4}; // cells per side of domain
    set_C(nblocks, C);


    bbounds dom;
    dom.max[0] = domain.max[0];
    dom.max[1] = domain.max[1];
    dom.max[2] = domain.max[2];
    dom.min[0] = domain.min[0];
    dom.min[1] = domain.min[1];
    dom.min[2] = domain.min[2];

    int N = 1;// depth of ghost region
    
    
    assert((domain.max[0] - domain.min[0]+1)%(C[0]*C[1]*C[2]) == 0);
    assert((domain.max[1] - domain.min[1]+1)%(C[0]*C[1]*C[2]) == 0);
    assert((domain.max[2] - domain.min[2]+1)%(C[0]*C[1]*C[2]) == 0);

    // use the partitioner to identify the boundary int value for bottom corner of cell
    // std::map<int, std::vector<float>> data;
    // std::vector<int> bid_to_rank(C*C*C,0);
    // std::vector<int> weights(C*C*C);
    // int bside[3];
    // std::vector<int> partn; // map of gid to current partn

    // MESH_DATA mesh_data;

    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
        // dprint("in master. C={%d %d %d}, nblocks %d", C[0], C[1], C[2], nblocks);
        b->cid_to_rank.resize(C[0]*C[1]*C[2]);
        b->weights.resize(C[0]*C[1]*C[2]);
        partition(world, dom, &C[0], b->data,  world.rank(), world.size(), b->cid_to_rank, &b->cside[0], b->partn, b->mesh_data, b);

        read_data(world, infile.c_str(), b->data, b->weights, &C[0], &b->cside[0], b, dom);

        // assign and send
        assign(world, b->data, b->particles, b->weights, b->partn, b->mesh_data, b, cp, assigner, b->bounds );

    });

    // receive and update data
    bool remote = true;
    master.exchange(remote);
    master.foreach(&remote_deq);

   master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) { 

        update_cid_to_rank_vector(world, b->data, b->cid_to_rank);

   });


    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {

        //  if (world.rank()==0){
        //     for (auto it=b->data.begin(); it!=b->data.end() ; ++it){
        //         dprint("ASDF %d", it->first);
        //     }
        //     dprint("HERE!!!! %ld, %ld", b->data.size(), b->particles.size());
        // }

      
        init = seed(b, dom, &C[0], seed_rate, world.rank());
        
        // dprint("initialized %ld", init);

        // std::map<int, std::vector<BEndPt>>::iterator it = b->particles.begin();
        // while (it != b->particles.end()){
        //     // if(cp.gid()==0 ){
        //         // dprint("rank [%d %d], cid %d, particles %ld, seed_rate %f", world.rank(), cp.gid(), it->first, it->second.size(), seed_rate);
            
        //     // }
        //     it++;
        // }

        // // interpolate
        // float pt[3] = {100.3,100.4,100.5};

        // // dprint("inblock %d, rank %d", in_block(&pt[0], b->data, dom, C), world.rank());
        // int cblock = in_block(&pt[0], b->data, b->data_ghost, dom, C);
        // if (cblock>-1){

        //     float v[3];
        //     lerp(&pt[0], b->bounds[cblock], b->data[cblock], &v[0]);

        //     // dprint("interpolated %f %f %f | cblock %d", v[0], v[1], v[2], cblock);
        // }

        // 

    });

    // update cid_to_rank
    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
        update_cid_to_rank_vector(world, b->data, b->cid_to_rank);
    });

    // get ghost cells based on new partition
    get_ghost_cells(master, assigner, N, dom, C, world.rank());

    //  master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
    //         dprint("myrank %d, data_ghost %ld", cp.gid(), b->data_ghost.size());
    //         // b->data_ghost.clear();
    //  });

    

    size_t init_global=0;
    diy::mpi::all_reduce(world, init, init_global, std::plus<size_t>());
    size_t done_global=0;
   

    int nrounds = 20;
    for (int round=0; round<nrounds; round++){

        done_global = 0;
        diy::mpi::all_reduce(world, done, done_global, std::plus<size_t>());
        if (done_global==init_global){
            if (world.rank() == 0)
            dprint("done %ld init %ld, breaking now", done_global, init_global);
            break;
        }

        if (world.rank()==0)
            dprint("!! starting round %d", round);

        master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {

            std::map<int, std::vector<float>>::iterator it = b->data.begin();
            while (it != b->data.end()){

                int cid = it->first;
                
                // dprint("rank %d, particles %d", world.rank(), b->particles.size());
                // continue;

                 //  if (round>0){
                        // dprint("rank %d cid %d numcells %ld, %ld", world.rank(), cid, b->data.size(), b->particles.size());
                        
                    //  }
                

                // iterate over particles in cell i
                if (b->particles.find(cid)!=b->particles.end())
                for (size_t j=0; j<b->particles[cid].size(); j++){
                    
                    BEndPt &cur_p = b->particles[cid][j];
                    bool finished = false;

                    // if (cur_p.pid != 50430080)
                    //     continue;
                    // else{
                        dprint("pid starting %d (%f %f %f ) in rank %d, round %d pcid %d, nsteps %d", cur_p.pid, cur_p[0], cur_p[1], cur_p[2], world.rank(), round, cur_p.cid, cur_p.nsteps);
                    // }
                   
                    BEndPt next_p;

                    BSegment s(cur_p);
                    while (badvect_rk1(cur_p, b, dom, &C[0], 0.05, next_p, world.rank()))
                    { // returns false if post cid is not in block

                        nsteps++;

                        // next_p.nsteps ++;
                        cur_p = next_p;

                        if (check)
                        {
                            BPt p;
                            p.coords[0] = cur_p[0];
                            p.coords[1] = cur_p[1];
                            p.coords[2] = cur_p[2];
                            s.pts.push_back(p);
                        }

                        // check for exit conditions : global bounds and max steps
                        if (cur_p.nsteps > max_steps || in_global_dom(dom, cur_p) == false)
                        {
                            finished = true;
                            // if (round > 2)
                                dprint("ending pid %d, %f %f %f, nsteps %d, pcid %d, cid %d, rank %d", cur_p.pid, cur_p[0], cur_p[1], cur_p[2], cur_p.nsteps, cur_p.cid, cid, world.rank());
                            break;
                        }

                        // dprint("step pid %d, %f %f %f, nsteps %d, pcid %d, cid %d, rank %d", cur_p.pid, cur_p[0], cur_p[1], cur_p[2], cur_p.nsteps, cur_p.cid, cid, world.rank());
                        
                    }

                        // push back into segment
                        b->segments.push_back(s);
                    // }
                    // if finished done++ else put in unfinised of the new cell
                    if (finished == true){
                        done++;
                    }else{

                        // // if in current rank's cell: key is a cid
                        if (b->data.find(next_p.cid) != b->data.end()){
                            b->unfinished_local[next_p.cid].push_back(next_p);
                        }
                        // else for export : key is a rank
                        else{
                            // if (next_p.pid==401)
                            // if (world.rank()==4)
                                dprint("unfin nonloc cid for rank %d, cid %d, pid %d, nsteps %d", world.rank(), next_p.cid, next_p.pid, next_p.nsteps);
                            int to_rank = b->cid_to_rank[next_p.cid];
                            b->unfinished_nonlocal[to_rank].push_back(next_p);
                        }

                    }


                }
                it++;
            }

            

            

            // reenqueue local particles and clears finished particles
            enqueue_local_particles(b->unfinished_local, b->particles);



             // update cell weights
            update_weights(world, b->particles, b->weights);



        });

        // exchange non local particles
        exchange_nonlocal_particles(master, assigner);

       
        // rebalance cells

        master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
            assign(world, b->data, b->particles, b->weights, b->partn, b->mesh_data, b, cp, assigner, b->bounds);
        });

        // receive and update data
        bool remote = true;
        master.exchange(remote);
        master.foreach (&remote_deq);

        master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
            update_cid_to_rank_vector(world, b->data, b->cid_to_rank);

            if (world.rank()==0)
                pvi(b->cid_to_rank);
        });
       

        if (world.rank() == 0)
            dprint("getting ghost for round %d", round);
        // get ghost cells based on new partition

       master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
            // dprint("myrankinloop %d, data_ghost %ld", cp.gid(), b->data_ghost.size());
            // dprint("!myrankinloop %d, data %ld", cp.gid(), b->data.size());
            b->data_ghost.clear();
            b->bounds_ghost.clear();
       });
        
       
       get_ghost_cells(master, assigner, N, dom, &C[0], world.rank());
        
    //    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
    //         // dprint("myrankinloop2 %d, data_ghost %ld", cp.gid(), b->data_ghost.size());
    //         dprint("!myrankinloop2 %d, data %ld", cp.gid(), b->data.size());
    //         b->data_ghost.clear();
    //         b->bounds_ghost.clear();
    //    });

    }

    
    size_t nsteps_global=0;
    diy::mpi::reduce(world, nsteps, nsteps_global, 0, std::plus<size_t>());

    

    if (world.rank()==0){
        dprint("nsteps_global %ld, init_global %ld, done_global %ld", nsteps_global, init_global, done_global);
    }
  
   

    // write trajectory segments for validation
    if (check){
         write_traces(master, assigner, decomposer, world);

    }

    // clean up
    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
            if(b->mesh_data.numMyPoints > 0){
                free(b->mesh_data.myGlobalIDs);
                free(b->mesh_data.x);
                free(b->mesh_data.y);
                free(b->mesh_data.z);
            }
    });

    

    if (world.rank() ==0)
        dprint("done");
    return 0;
}