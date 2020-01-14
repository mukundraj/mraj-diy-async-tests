//---------------------------------------------------------------------------
//
// particle advection
//
// courtesy Hanqi Guo
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
//
// copied with permission by Tom Peterka
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------

#ifndef _ADVECT_H
#define _ADVECT_H

#include <stdbool.h>
#include <functional>

bool advect_rk1(
        const int   *st,
        const int   *sz,
        const float **vec,
        float       *X,
        float        h,
        float       *Y);

// rk1 for continuous domain
bool cadvect_rk1(
        const float   *st,
        const float   *sz,
        const float **vec,
        float       *X,
        float        h,
        float       *Y);

bool advect_brown(
        const int   *st,
        const int   *sz,
        const float **vec,
        float       *X,
        float        h,
        float       *Y);


bool advect_rk4(
        const int   *st,
        const int   *sz,
        const float **vec,
        float       *pt,
        float       h,
        float       *Y );

#endif
