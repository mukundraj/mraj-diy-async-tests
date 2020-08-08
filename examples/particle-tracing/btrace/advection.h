#ifndef ADVECTION_H
#define ADVECTION_H

#include "bmisc.h"
#include "bblock.hpp"


int pos2cgid(float px, float py, float pz, bbounds &dom, int C);

void gid2bounds(int gid, int *cside, int C, bbounds &bnd);

bool inside(float *pt, bbounds &bnd);

void seed(BBlock *b, bbounds &dom, int C, float sr, int rank);


#endif