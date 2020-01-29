#ifndef TRACER_H
#define TRACER_H

#include "readerwriterqueue.h"

class Tracer{


    public:
        void exec(); // main function
        void exec_comm(); // handles communication in the main thread
};





#endif