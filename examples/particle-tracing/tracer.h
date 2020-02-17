#ifndef TRACER_H
#define TRACER_H

#include "readerwriterqueue.h"
#include "concurrentqueue.h"

#include "atomicops.h"

#include "async_comm.h"
#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include "state_exchanger.h"
#include <diy/master.hpp>
#include <thread>

#include "ptrace.hpp"
#include "block.hpp"
#include "advect.h"
#include "utils.hpp"


#include <mutex>

using namespace moodycamel;

template <typename MsgType>
struct Worker
{

    using AtomicQueue = ConcurrentQueue<MsgType>;
    AtomicQueue incoming_queue = AtomicQueue(10000);
    AtomicQueue outgoing_queue = AtomicQueue(10000);
};

class Tracer : public Worker<std::unique_ptr<tracer_message>>
{

public:
    AsyncComm acomm;
    StateExchanger state;
    diy::mpi::communicator world;
    std::unique_ptr<std::thread> worker;

    std::mutex mutexx;



    Tracer(diy::mpi::communicator &world, diy::Master &master, diy::Assigner &assigner);
    void exec(diy::Master &master_iex, const CBounds &cdomain,
                         const int max_steps,
                        //  map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                         size_t &nsteps,
                         size_t &ntransfers,
                         bool prediction,
                         double &time_trace, 
                         Tracer &tracer, 
                         diy::Assigner &cassigner);      // main function
    void exec_comm(size_t &nsteps); // handles communication in the main thread
    // void tracer_join();
};

#endif