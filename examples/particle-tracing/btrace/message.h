#ifndef BTRACE_HPP
#define BTRACE_HPP

#include <vector>

struct bbounds{

    int max[3];
    int min[3];

};

struct datablock{

    std::vector<float> particles;
    int cellid;
    std::vector<float> data;


    


};




#endif