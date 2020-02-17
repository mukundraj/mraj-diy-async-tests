#ifndef _STATE_EXCHANGER_H
#define _STATE_EXCHANGER_H

#include <diy/mpi.hpp>
#include "../misc.h"

// modified DIY's diy/detail/master/iexchange-collective.hpp
struct StateExchanger
{
  StateExchanger(diy::mpi::communicator comm_ = MPI_COMM_WORLD) : comm(comm_) 
  {
    local_work_ = 0;
    dirty = 0;
    state = 0;
    this->rank = rank;
  }

  inline bool       all_done();                    // get global all done status
  inline void       add_work(int work=1);            // add work to global work counter
  inline void       dec_work(int work=1) {add_work(-work);}
  inline void       control();

  std::atomic<int>  local_work_;
  std::atomic<int>  dirty;
  int               local_dirty, all_dirty;

  std::atomic<int> state;
  diy::mpi::request r;
  diy::mpi::communicator comm;
  int rank;
};


bool
StateExchanger::
all_done()
{
  return state == 3;
}

void
StateExchanger::
add_work(int work)
{
  local_work_ += work;
  if (local_work_ > 0)
    dirty = 1;
}

void
StateExchanger::
control()
{
  if (state == 0 && local_work_ == 0) { // local work done, inform others if haven't yet
    r = ibarrier(comm);
    dirty = 0;
    state = 1;
  } else if (state == 1) { // informed others, wait to hear from those working
    diy::mpi::optional<diy::mpi::status> ostatus = r.test();
    if (ostatus)
    {
      local_dirty = dirty;
      r = diy::mpi::iall_reduce(comm, local_dirty, all_dirty, std::logical_or<int>());
      state = 2;
    }
  } else if (state == 2) { // heard back from those working, check if all work over
    diy::mpi::optional<diy::mpi::status> ostatus = r.test();
    if (ostatus)
    {
      // dprint("dirty %d all_dirty %d rank %d", int(dirty), int(all_dirty), rank);
      if (all_dirty == 0)    // done
        state = 3;
      else
        state = 0;          // reset
    }
  }
}

#endif
