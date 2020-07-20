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

    float ver;
	int rc = Zoltan_Initialize(argc, argv, &ver);
	
	if (rc != ZOLTAN_OK)
	{
		printf("sorry...\n");
		MPI_Finalize();
		exit(0);
	}


    // read in arguments

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


   
    
    int C = 4; // blocks per side of domain
    bbounds dom = {domain.max[0], domain.max[1], domain.max[2], domain.min[0], domain.min[1], domain.min[2]};
    
    
    assert((domain.max[0] - domain.min[0]+1)%(C*C*C) == 0);
    assert((domain.max[1] - domain.min[1]+1)%(C*C*C) == 0);
    assert((domain.max[2] - domain.min[2]+1)%(C*C*C) == 0);

    // use the partitioner to identify the boundary int value for bottom corner of cell
    std::map<int, std::vector<float>> data;
    std::vector<int> bid_to_rank(C*C*C,0);
    std::vector<int> weights(C*C*C);

    partition(dom, C, data, world.rank(), world.size(), bid_to_rank);

    // read data blocks for current process, set weights to be uniform
    read_data(data, weights);

    



    // use the assigner to assigne the cells to partitions (uses Zoltan), also move data blocks around accordingly




    // perform advection and repartitioning over rounds

    if (world.rank() ==0)
        dprint("done");
    return 0;
}