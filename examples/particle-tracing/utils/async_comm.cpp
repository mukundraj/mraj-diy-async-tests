#include "async_comm.h"
#include "../ptrace.hpp"
#include "../block.hpp"
#include "../misc.h"
#include "../tracer.h"
#include <mutex>



void AsyncComm::init(diy::mpi::communicator& world, diy::Master &master_iex, diy::Assigner &assigner){
    this->world = &world;
    this->master_iex = &master_iex;
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


        // push request object handle to queue
        MPI_Request request;
        MPI_Issend(uptr_msg.get(), 1, mpi_tracer_msg, uptr_msg->dest_gid, 0, *world, &request);

        send_dests.push_back(uptr_msg->dest_gid);
        // send_dests.push(77);
        stats.push_back(request);
        // push message buffer to queue
        in_transit_msgs.push_back(std::move(uptr_msg));

        
        
}
bool AsyncComm::check_nbrs_for_incoming(ConcurrentQueue<std::unique_ptr<tracer_message>> &incoming_queue, StateExchanger& state){


    MPI_Status status;

    for (size_t i=0; i<nbr_procs.size(); i++){
        
       
        int nbr_proc = nbr_procs[i];

            
                int flag=0;
                MPI_Status status;
                MPI_Iprobe( nbr_proc, 0, *world, &flag, &status );
                // if (world->rank() == 14)
                //     pvi(nbr_procs);


                
               
                while(flag==1){
                        
                        state.add_work();
                        // dprint("flagged in %d", world->rank());
                        std::unique_ptr<tracer_message> msg(new tracer_message());
                        
                        
                        MPI_Recv( msg.get(), 1, mpi_tracer_msg, nbr_proc, 0, *world, &status );

                        // tracer_message msgg;
                        // int msgg;
                        // MPI_Recv( &msgg, 1, MPI_INT, nbr_proc, 0, *world, &status );
                        
                        // dprint("Received from %d, pid: %d, nsteps: %d, (%f %f %f)", nbr_proc, msg->pid, msg->nsteps, msg->p[0], msg->p[1], msg->p[2]);

                        incoming_queue.enqueue(std::move(msg));
                        flag = 0;
                        MPI_Iprobe( nbr_proc, 0, *world, &flag, &status );
                }

    }


    return true;
}
bool AsyncComm::check_sends_complete(StateExchanger &state){

    MPI_Request request;
    bool remaining = false;
    MPI_Status status;



        // do{
        //     if (!stats.empty()){

        //         request = stats.front();
                

        //         int flag = false;
        //         MPI_Test(&request, &flag, &status);
        //         // if complete, continue with next on queue
        //         if (flag){
                    
        //             // MPI_Request_free(&request);
        //             stats.pop();
        //             in_transit_msgs.pop();
        //             send_dests.pop();
        //             state.dec_work();


        //             if (stats.empty())
        //                remaining = false;
        //             else
        //                remaining = true;
        //         }
               
        //     }

        // }while(remaining);

        for (size_t i=0; i<stats.size(); i++){

            int flag;
            MPI_Test(&stats[i], &flag, &status); 
            if (flag){
                stats.erase(stats.begin()+i);
                in_transit_msgs.erase(in_transit_msgs.begin()+i);
                send_dests.erase(send_dests.begin()+i);
                state.dec_work();
            }

        }



    if (stats.size()>0)
        return false;
    else
        return true;
}


size_t AsyncComm::get_stat_size(){
      return stats.size();
    }