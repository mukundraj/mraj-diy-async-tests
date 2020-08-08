#ifndef BMISC_HPP
#define BMISC_HPP

struct bbounds{

    int max[3];
    int min[3];

    int cside[3]; // elements per cell side along each dimension

};


struct BEndPt{

    int pid;
    int gid;
    int nsteps;
    // std::vector<float> pt(3);
    float pt[3];

};

#endif