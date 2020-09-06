#ifndef ADVECTION_H
#define ADVECTION_H

#include "bmisc.h"
#include "bblock.hpp"

bool in_global_dom(bbounds &dom, BEndPt &p);

bool badvect_rk1(BEndPt &pt, BBlock *b, bbounds &dom, int *C, float h, BEndPt &next_p, int rank);

int pos2cgid(float px, float py, float pz, bbounds &dom, int *C);

void gid2bounds(int gid, int *cside, int *C, bbounds &dom, bbounds &bnd);

bool inside(float *pt, bbounds &bnd);

int seed(BBlock *b, bbounds &dom, int *C, float sr, int rank);

int in_block(float *pt, std::map<int, std::vector<float>> &data, std::map<int, std::vector<float>> &data_ghost, bbounds &dom, int *C);

void lerp(float *pt, bbounds &bnd,  std::vector<float> &vec, float *vars);

void enqueue_local_particles(std::map<int, std::vector<BEndPt>> &unfinished_local, 
std::map<int, std::vector<BEndPt>> &particles);

void exchange_nonlocal_particles(diy::Master &master, const diy::Assigner& assigner);

#endif