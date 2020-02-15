// #include "state_exchanger.h"

// inline bool
// StateExchanger::
// all_done()
// {
//   return state == 3;
// }

// inline void
// StateExchanger::
// add_work(int work)
// {
//   local_work_ += work;
//   if (local_work_ > 0)
//     dirty = 1;
// }

// void  
// StateExchanger::
// dec_work(int work) {add_work(-work);}

// void
// StateExchanger::
// control()
// {
//   if (state == 0 && local_work_ == 0) { // local work done, inform others if haven't yet
//     r = ibarrier(comm);
//     dirty = 0;
//     state = 1;
//   } else if (state == 1) { // informed others, wait to hear from those working
//     diy::mpi::optional<diy::mpi::status> ostatus = r.test();
//     if (ostatus)
//     {
//       local_dirty = dirty;
//       r = diy::mpi::iall_reduce(comm, local_dirty, all_dirty, std::logical_or<int>());
//       state = 2;
//     }
//   } else if (state == 2) { // heard back from those working, check if all work over
//     diy::mpi::optional<diy::mpi::status> ostatus = r.test();
//     if (ostatus)
//     {
//       if (all_dirty == 0)    // done
//         state = 3;
//       else
//         state = 0;          // reset
//     }
//   }
// }
