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
            fmt::print(stderr, "Remote dequeue: gid {} received value {} from gid {}\n", cp.gid(), recvd_data.data.size(), recvd_data.from_proc);
            for (size_t i=0; i<recvd_data.cgid.size(); i++){
                b->data[recvd_data.cgid[i]] = recvd_data.data[i];
                b->particles[recvd_data.cgid[i]] = std::move(recvd_data.particles[i]);
                b->bounds[recvd_data.cgid[i]] = std::move(recvd_data.bounds[i]);
            }   
        }
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
        for (size_t i=0; i<b->segments.size(); i++){
            for (size_t j=0; j<b->segments[i].pts.size(); j++){
                segs.push_back(b->segments[i].pts[j].coords[0]);
                segs.push_back(b->segments[i].pts[j].coords[1]);
                segs.push_back(b->segments[i].pts[j].coords[2]);
            }
            sizes.push_back(b->segments[i].pts.size());
        }

    });

    diy::mpi::gather(world, segs, segs2, 0);
    diy::mpi::gather(world, sizes, all_sizes, 0);

   



    if (master.communicator().rank() == 0)
    {

        
        std::vector<Segment> all_segs;
        
        for(size_t i=0; i<all_sizes.size(); i++){
            for (size_t j=0; j<all_sizes[i].size(); j++){
                Segment seg;
                size_t idx=0;
                for(int k=0; k<all_sizes[i][j]; k++){
                    Pt p;
                    p.coords[0] = segs2[i][idx*3];
                    p.coords[1] = segs2[i][idx*3+1];
                    p.coords[2] = segs2[i][idx*3+2];
                    seg.pts.push_back(p);
                    idx++;
                }

            }
            dprint("sizes %ld ", all_sizes[i].size());
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
    // definitions
    // element: 1x1x1 pixel
    // cell: smallest unit that can be moved around (formerly block)
    // block: data consisting of multiple cells in a process

    int C = 4; // cells per side of domain
    bbounds dom;
    dom.max[0] = domain.max[0];
    dom.max[1] = domain.max[1];
    dom.max[2] = domain.max[2];
    dom.min[0] = domain.min[0];
    dom.min[1] = domain.min[1];
    dom.min[2] = domain.min[2];

    int N = 1;// depth of ghost region
    
    
    assert((domain.max[0] - domain.min[0]+1)%(C*C*C) == 0);
    assert((domain.max[1] - domain.min[1]+1)%(C*C*C) == 0);
    assert((domain.max[2] - domain.min[2]+1)%(C*C*C) == 0);

    // use the partitioner to identify the boundary int value for bottom corner of cell
    // std::map<int, std::vector<float>> data;
    // std::vector<int> bid_to_rank(C*C*C,0);
    // std::vector<int> weights(C*C*C);
    // int bside[3];
    // std::vector<int> partn; // map of gid to current partn

    // MESH_DATA mesh_data;

    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
        dprint("in master");
        b->cid_to_rank.resize(C*C*C);
        b->weights.resize(C*C*C);
        partition(world, dom, C, b->data,  world.rank(), world.size(), b->cid_to_rank, &b->cside[0], b->partn, b->mesh_data, b);

        read_data(world, infile.c_str(), b->data, b->weights, C, &b->cside[0], b, dom);

        // assign and send
        assign(world, b->data, b->particles, b->weights, b->partn, b->mesh_data, b, cp, assigner, b->bounds );


        // int gid = pos2cgid(130,10,10, dom, C);
        // dprint("gid %d", gid);
        
        // bbounds bnd;
        // gid2bounds(gid, &dom.cside[0], C, bnd);
        // dprint("bnd [%d %d] [%d %d] [%d %d], rank %d", bnd.min[0], bnd.max[0], bnd.min[1], bnd.max[1], bnd.min[2], bnd.max[2], world.rank());

        

        // dprint("seeded %ld, seed_rate %f, rank %d", b->particles.size(), seed_rate, world.rank());

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

      
        seed(b, dom, C, seed_rate, world.rank());

        std::map<int, std::vector<BEndPt>>::iterator it = b->particles.begin();
        while (it != b->particles.end()){
            // if(cp.gid()==0 ){
                // dprint("rank [%d %d], cid %d, particles %ld, seed_rate %f", world.rank(), cp.gid(), it->first, it->second.size(), seed_rate);
            
            // }
            it++;
        }

        // interpolate
        float pt[3] = {100.3,100.4,100.5};

        // dprint("inblock %d, rank %d", in_block(&pt[0], b->data, dom, C), world.rank());
        int cblock = in_block(&pt[0], b->data, b->data_ghost, dom, C);
        if (cblock>-1){

            float v[3];
            lerp(&pt[0], b->bounds[cblock], b->data[cblock], &v[0]);

            // dprint("interpolated %f %f %f | cblock %d", v[0], v[1], v[2], cblock);
        }

        // 

    });

    // update cid_to_rank
    master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {
        update_cid_to_rank_vector(world, b->data, b->cid_to_rank);
    });

    // get ghost cells based on new partition
    get_ghost_cells(master, assigner, N, dom, C, world.rank());

    int nrounds = 1;
    for (int i=0; i<nrounds; i++){

        master.foreach ([&](BBlock *b, const diy::Master::ProxyWithLink &cp) {

            std::map<int, std::vector<float>>::iterator it = b->data.begin();
            while (it != b->data.end()){

                int i = it->first;

                // iterate over particles in cell i
                for (size_t j=0; j<b->particles[i].size(); j++){
                    
                    BEndPt &cur_p = b->particles[i][j];
                    bool finished = false;

                    BEndPt next_p;
                   
                    // if (cur_p.pid == 821){
                        // dprint("*************");
                        // print_cellids_in_block(b->data);
                        // print_cellids_in_block(b->data_ghost);
                        // dprint("startingg %d in rank %d, cid %d", cur_p.pid, world.rank(), cur_p.cid);
                        // dprint("stepsinit %f %f %f, nsteps %d, cid %d", cur_p[0], cur_p[1], cur_p[2], cur_p.nsteps, cur_p.cid);
                    // badvect_rk1(cur_p, b, dom, C, 0.05, next_p);
                        BSegment s;
                        while(badvect_rk1(cur_p, b, dom, C, 0.05, next_p)){ // returns false if post cid is not in block

                            // next_p.nsteps ++;
                            cur_p = next_p;

                            // dprint("steps %f %f %f, nsteps %d, cid %d", cur_p[0], cur_p[1], cur_p[2], cur_p.nsteps, cur_p.cid);

                            if (check){
                                BPt p;
                                p.coords[0] = cur_p[0];
                                p.coords[1] = cur_p[1];
                                p.coords[2] = cur_p[2];
                                s.pts.push_back(p);
                            }
                            
                            if (cur_p.nsteps > max_steps){
                                    finished = true;
                                    break;
                            }

                        }
                    // }
                    // push back into segment
                    b->segments.push_back(s);

                    // if finished done++ else put in unfinised of the new cell



                }

                it++;
            }

            // get ghosts
            b->data_ghost.clear();
            b->bounds_ghost.clear();

             // update cell weights
            update_weights(world, b->particles, b->weights);
                
        });

       
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
        });
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
            }
    });

    if (world.rank() ==0)
        dprint("done");
    return 0;
}