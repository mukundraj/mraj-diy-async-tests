#include "tracer.h"
#include <thread>
#include "misc.h"
#include "state_exchanger.h"





// Only for iexchange
void trace_particles_iex(Block *b,
                         const diy::Master::ProxyWithLink &cp,
                         const CBounds &cdomain,
                         const int max_steps,
                        //  map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                         size_t &nsteps,
                         size_t &ntransfers,
                         bool prediction,
                         double &time_trace, 
                         Tracer &tracer, 
                         diy::Assigner &cassigner)
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

    // check if any particles in tracer.incoming_queue and load to b->particles

    std::unique_ptr<tracer_message> *uptrptr_msg = tracer.incoming_queue.peek();

    while (uptrptr_msg != nullptr)
    {
        dprint("incoming");
        std::unique_ptr<tracer_message> uptr_msg;
        tracer.incoming_queue.try_dequeue(uptr_msg); 
        EndPt in_pt;
        in_pt.pid = uptr_msg->pid;
        in_pt.gid = uptr_msg->gid;
        in_pt.nsteps = uptr_msg->nsteps;
        in_pt.predonly = uptr_msg->predonly;
        // in_pt.efs = uptr_msg->efs;
        in_pt[0] = uptr_msg->p[0]; in_pt[1] = uptr_msg->p[1]; in_pt[2] = uptr_msg->p[2];

        b->particles.push_back(in_pt);
        uptrptr_msg = tracer.incoming_queue.peek();
        dprint("incomeing pid %d, (%f %f %f)", in_pt.pid, in_pt[1], in_pt[2], in_pt[2]);

    }
    
    for (auto i = 0; i < b->particles.size(); i++)
    {
        Pt &cur_p = b->particles[i].pt; // current end point
        Segment s(b->particles[i]);     // segment with one point p
        Pt next_p;                      // coordinates of next end point
        bool finished = false;

        if (b->particles[i].pid == 800)
            dprint("%f %f %f @ %d", cur_p.coords[0], cur_p.coords[1], cur_p.coords[2], cp.gid());
        // trace this segment until it leaves the block
        double time_start = MPI_Wtime();

        while (cadvect_rk1(st, sz, vec, cur_p.coords.data(), 0.5, next_p.coords.data()))
        {
            nsteps++;
            b->particles[i].nsteps++;
            s.pts.push_back(next_p);
            cur_p = next_p;

            // if (b->particles[i].pid == 800)
            //     dprint("%f %f %f @ %d", cur_p.coords[0], cur_p.coords[1], cur_p.coords[2], cp.gid());

            if (b->particles[i].nsteps >= max_steps)
            {
                finished = true;
                break;
            }

            // if predicting, add copy coordinates to EndPt and add to b->particles_store
            if (prediction == true && cinside(cur_p, cdomain) && nsteps % 2 == 0)
            {   
                EndPt way_pt;
                way_pt[0] = cur_p.coords[0];
                way_pt[1] = cur_p.coords[1];
                way_pt[2] = cur_p.coords[2];
                way_pt.predonly = true;
                b->particles_store.push_back(way_pt);
            }

        }

        time_trace += MPI_Wtime() - time_start;

        b->segments.push_back(s);

        if (!cinside(next_p, cdomain))
            finished = true;

        if (finished)
        { // this segment is done
            b->done++;
            dprint("pid: %d, step %d, %f %f %f", b->particles[i].pid, b->particles[i].nsteps, cur_p.coords[0], cur_p.coords[1], cur_p.coords[2]);
            tracer.state.dec_work();
        }
        else // find destination of segment endpoint
        {   
            vector<int> dests;
            vector<int>::iterator it = dests.begin();
            insert_iterator<vector<int>> insert_it(dests, it);

            // utl::in(*l, next_p.coords, insert_it, cdomain, 1);
            utl::in(*l, next_p.coords, insert_it, cdomain, false);

            EndPt out_pt(s);
            // out_pt.pid = b->particles[i].pid;
            out_pt.nsteps = b->particles[i].nsteps;
            if (dests.size())
            {   
                ntransfers ++;
                diy::BlockID bid = l->target(dests[0]); // in case of multiple dests, send to first dest only

                // debug
                fmt::print(stderr, "{}: gid {} enq to gid {}, steps {}, ({}, {}, {})\n", out_pt.pid, cp.gid(), bid.gid,  out_pt.nsteps, out_pt.pt.coords[0], out_pt.pt.coords[1], out_pt.pt.coords[2]);

                int nbr_proc = cassigner.rank(bid.gid);
                std::unique_ptr<tracer_message> uptr_msg(new tracer_message()); 
                uptr_msg->pid = out_pt.pid;
                uptr_msg->gid = out_pt.gid;
                uptr_msg->nsteps = out_pt.nsteps;
                uptr_msg->predonly = out_pt.predonly;
                uptr_msg->p[0] = out_pt[0]; uptr_msg->p[1] = out_pt[1]; uptr_msg->p[2] = out_pt[2];
                uptr_msg->dest_gid = nbr_proc;
                // dprint("pid: %d (%f %f %f)", uptr_msg->pid , uptr_msg->p[0], uptr_msg->p[1], uptr_msg->p[2]);

                tracer.outgoing_queue.enqueue(std::move(uptr_msg));

                // if (IEXCHANGE) // enqueuing single endpoint allows fine-grain iexchange if desired
                //     cp.enqueue(bid, out_pt);
                // else
                //     outgoing_endpts[bid].push_back(out_pt); // vector of endpoints
            }
        }
        // dprint("BREAKING HERE");
        // break;

    }
    b->particles.clear();

    if (prediction==false)
        b->particles_store.clear();
}

Tracer::Tracer(diy::mpi::communicator& world, diy::Master &master, diy::Assigner &assigner){
    
        acomm.init(world, master, assigner);
        this->world = world;
        this->state.rank = world.rank();
}

void Tracer::exec(diy::Master &master_iex, const CBounds &cdomain,
                         const int max_steps,
                        //  map<diy::BlockID, vector<EndPt>> &outgoing_endpts,
                         size_t &nsteps,
                         size_t &ntransfers,
                         bool prediction,
                         double &time_trace, 
                         Tracer &tracer, 
                         diy::Assigner &cassigner){


    
    // create worker thread, keeps checking atomic status to decide whether to end
    worker = std::unique_ptr<std::thread>(new std::thread([&]{
        while (!state.all_done())
            {
                    
                master_iex.foreach ([&](Block *b, const diy::Master::ProxyWithLink &icp) {
                    trace_particles_iex(b, icp, cdomain, max_steps, nsteps, ntransfers, prediction, time_trace, tracer, cassigner);
                });
               
            }
    }));

    dprint("worker working");

    // enter communication loop, local work count has already been updated
    exec_comm();

   


    // wait for worker threads to join
    worker->join();
}

// void Tracer::tracer_join(){

//     worker->join();
// }

void Tracer::exec_comm(){

    // while not finished
    while(!state.all_done())
    {   
        state.control(); // update state exchanger

        // send out any message to be sent out
        bool success = false;
        do{

            std::unique_ptr<tracer_message> uptr_msg;
            success = outgoing_queue.try_dequeue(uptr_msg); 
            
            if (success){
                dprint("sending from %d to %d", world.rank(), uptr_msg->dest_gid);
                // state.dirty = true;
                acomm.send_to_a_nbr(uptr_msg);
            }
        }while(success);
        
        // check for incoming and receive and enqueue if any
        acomm.check_nbrs_for_incoming(incoming_queue, state);

            // success = false;
            // do{
            //     tracer_message msg;
            //     success = incoming_queue.try_dequeue(uptr_msg);

            // }while(success);


        // check if send messages reached destination and update dirty accordingly
        bool val = acomm.check_sends_complete(state);
            // state.dirty = false;
        // if (val==false)
        //     dprint("sends comple %d, rank %d", val, world.rank());
    }

    
}