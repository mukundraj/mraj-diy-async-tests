#include "tracer.h"
#include <thread>


void Tracer::exec(){

    // create worker thread, keeps checking atomic status to decide whether to end
    std::thread worker([]{});

    // enter communication loop, updates atomic status
    exec_comm();


    // wait for worker threads to join
    worker.join();
}

void Tracer::exec_comm(){



}