#ifndef ASYNC_COMM_H
#define ASYNC_COMM_H

#include "readerwriterqueue.h"
#include "atomicops.h"

// #include <mpi.h>
#include <diy/mpi.hpp>
#include <diy/master.hpp>
#include <queue>
#include <memory> 
#include "state_exchanger.h"

using namespace moodycamel;

struct tracer_message {
  int pid;             // particle id
  int gid;            // seed gid
  int nsteps = 0;
  int predonly;
  int efs;            // estimated final steps
  float p[3];
  int dest_gid;       // destination of message; used only locally; will not be sent out.
};


class AsyncComm{

    diy::mpi::communicator *world;
    diy::Master *master;
    diy::Assigner *assigner;
    // std::queue<MPI_Status*> stats;
    std::queue<MPI_Request> stats;
    std::queue<std::unique_ptr<tracer_message>> in_transit_msgs;


    /* create MPI type for tracer message */
    const int nitems=7;
    int          blocklengths[7] = {1, 1, 1, 1, 1, 3, 1};
    MPI_Datatype types[7] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_INT};
    MPI_Datatype mpi_tracer_msg;
    MPI_Aint     offsets[7];


  public:
    void init(diy::mpi::communicator& world, diy::Master &master, diy::Assigner &assigner);

    void send_to_a_nbr(std::unique_ptr<tracer_message> &uptr_msg);
    bool check_nbrs_for_incoming(ReaderWriterQueue<std::unique_ptr<tracer_message>> &incoming_queue, StateExchanger& state);
    bool check_sends_complete(StateExchanger& state);


};

#endif