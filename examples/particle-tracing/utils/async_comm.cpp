#include "async_comm.h"
#include "../ptrace.hpp"
#include "../block.hpp"
#include "../misc.h"
#include "../tracer.h"



void AsyncComm::init(diy::mpi::communicator& world, diy::Master &master, diy::Assigner &assigner){
    this->world = &world;
    this->master = &master;
    this->assigner = &assigner;

    offsets[0] = offsetof( tracer_message , pid);
    offsets[1] = offsetof( tracer_message , gid);
    offsets[2] = offsetof( tracer_message , nsteps);
    offsets[3] = offsetof( tracer_message , predonly);
    offsets[4] = offsetof( tracer_message , efs);
    offsets[5] = offsetof( tracer_message , p);
    offsets[6] = offsetof( tracer_message , dest_gid);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_tracer_msg);
    MPI_Type_commit(&mpi_tracer_msg);

}

void AsyncComm::send_to_a_nbr(std::unique_ptr<tracer_message> &uptr_msg){

        // int MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
            //    MPI_Comm comm, MPI_Request *request);
        
        // tracer_message* msg = new tracer_message(); 
        // msg->gid = 77;
        // msg->steps = 99;
        // msg->efs = 101;
        // msg->p[0] = 0.2; msg->p[1] = 0.3; msg->p[2] = 0.4;

        // push request object handle to queue
        MPI_Request request;
        MPI_Issend(uptr_msg.get(), 1, mpi_tracer_msg, uptr_msg->dest_gid, 0, *world, &request);
        // tracer_message msgg;
        // MPI_Issend(&msgg, 1, mpi_tracer_msg, uptr_msg->dest_gid, 0, *world, &request);
        // int msgg = 1;
        // MPI_Issend(&msgg, 1, MPI_INT, uptr_msg->dest_gid, 0, *world, &request);

        stats.push(request);

        // push message buffer to queue
        // std::unique_ptr<tracer_message> uptr_msg (msg);
        in_transit_msgs.push(std::move(uptr_msg));
        
}
bool AsyncComm::check_nbrs_for_incoming(ReaderWriterQueue<std::unique_ptr<tracer_message>> &incoming_queue, StateExchanger& state){


    MPI_Status status;
     
    
    master->foreach ([&](Block *b, const diy::Master::ProxyWithLink &cp) {
        diy::RegularLink<Bounds> *l = static_cast<diy::RegularLink<Bounds> *>(cp.link());
        for (size_t i = 0; i < l->size(); ++i)
        {
            int nbr_gid = l->target(i).gid;
            int nbr_proc = assigner->rank(nbr_gid);

            
                int flag=0;
                MPI_Status status;
                MPI_Iprobe( nbr_proc, 0, *world, &flag, &status );
                // if (world->rank()==2){
                //     // dprint("probe 2 received %d", flag);
                //     // dprint("nbr rank %d, nbr gid %d", nbr_proc, nbr_gid);
                // }
                while(flag==1){
                        state.add_work();
                        dprint("flagged in %d", world->rank());
                        std::unique_ptr<tracer_message> msg(new tracer_message());
                        // tracer_message msgg;
                        
                        MPI_Recv( msg.get(), 1, mpi_tracer_msg, nbr_proc, 0, *world, &status );

                        // int msgg;
                        // MPI_Recv( &msgg, 1, MPI_INT, nbr_proc, 0, *world, &status );
                        
                        dprint("Received from %d, pid: %d, nsteps: %d, (%f %f %f)", nbr_proc, msg->pid, msg->nsteps, msg->p[0], msg->p[1], msg->p[2]);

                        incoming_queue.enqueue(std::move(msg));
                        // state.dec_work();
                        flag = 0;
                        MPI_Iprobe( nbr_proc, 0, *world, &flag, &status );
                }
        }

     });

    return true;
}
bool AsyncComm::check_sends_complete(StateExchanger &state){

    MPI_Request request;
    bool remaining = false;
    MPI_Status status;



        do{
            if (!stats.empty()){

                request = stats.front();

                int flag = false;
                MPI_Test(&request, &flag, &status);
                // if complete, continue with next on queue
                if (flag){
                    
                    stats.pop();
                    in_transit_msgs.pop();
                    state.dec_work();

                    if (stats.empty())
                       remaining = false;
                    else
                       remaining = true;
                }
               
            }

        }while(remaining);

    if (stats.size()>0)
        return false;
    else
        return true;
}