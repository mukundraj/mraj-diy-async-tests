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

#include "advect.h"
#include "lerp.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string.h>
#include "misc.h"

bool advect_brown(
        const int *st,      // min. corner of block
        const int *sz,      // size (number of points) in block
        const float **vec,  // vector field
        float *X,           // input point
        float h,            // step size
        float *Y = NULL)    // output point if not NULL, otherwise in X
{
    if (!inside(3, st, sz, X)) return false;

    if (Y)
    {
        Y[0] = X[0] + (float)rand()/RAND_MAX - 0.5;
        Y[1] = X[1] + (float)rand()/RAND_MAX - 0.5;
        Y[2] = X[2] + (float)rand()/RAND_MAX - 0.5;
    }
    else
    {
        X[0] = X[0] + (float)rand()/RAND_MAX - 0.5;
        X[1] = X[1] + (float)rand()/RAND_MAX - 0.5;
        X[2] = X[2] + (float)rand()/RAND_MAX - 0.5;
    }

    return true;
}

bool advect_rk1(
        const int *st,      // min. corner of block
        const int *sz,      // size (number of points) in block
        const float **vec,  // vector field
        float *X,           // input point
        float h,            // step size
        float *Y = NULL)    // output point if not NULL, otherwise in X
{   
    if (!inside(3, st, sz, X)) return false;

    float v[3];
    if (!lerp3D(X, st, sz, 3, vec, v))
        return false;

    

    if (Y)
    {
        Y[0] = X[0] + h * v[0];
        Y[1] = X[1] + h * v[1];
        Y[2] = X[2] + h * v[2];
    }
    else
    {
        X[0] = X[0] + h * v[0];
        X[1] = X[1] + h * v[1];
        X[2] = X[2] + h * v[2];
    }

    return true;
}

// rk1 for continuous domain
bool cadvect_rk1(
        const float *st,      // min. corner of block
        const float *sz,      // size (number of points) in block
        const float **vec,  // vector field
        float *X,           // input point
        float h,            // step size
        float *Y = NULL)    // output point if not NULL, otherwise in X
{
    // X[0] = 100.3;
    // X[1] = 100.4;
    // X[2] = 100.5;

    if (!cinside(3, st, sz, X)) return false;

    float v[3];
    if (!clerp3D(X, st, sz, 3, vec, v))
        return false;
    
    // dprint("vel %f %f %f", v[0], v[1], v[2]);    
    // exit(0);

    if (Y)
    {
        Y[0] = X[0] + h * v[0];
        Y[1] = X[1] + h * v[1];
        Y[2] = X[2] + h * v[2];
    }
    else
    {
        X[0] = X[0] + h * v[0];
        X[1] = X[1] + h * v[1];
        X[2] = X[2] + h * v[2];
    }

    return true;
}

bool advect_rk4(
        const int *st,      // min. corner of block
        const int *sz,      // size (number of points) in block
        const float **vec,  // vector field
        float *pt,          // input point
        float h,            // step size
        float* Y)           // output point
{
    int num_dims = 3;

    float p0[num_dims];
    memcpy(p0, pt, num_dims * sizeof(float));

    float v[num_dims];

    // 1st rk step
    if (!lerp3D(pt, st, sz, 3, vec, v))
        return false;
    float k1[num_dims];
    for (int i = 0; i < num_dims; i++)
        k1[i] = h * v[i];
    for (int i = 0; i < num_dims; i++)
        Y[i] = p0[i] + 0.5 * k1[i];

    // 2nd rk step
    if (!lerp3D(Y, st, sz, 3, vec, v))
        return true;
    float k2[num_dims];
    for (int i = 0; i < num_dims; i++)
        k2[i] = h * v[i];
    for (int i = 0; i < num_dims; i++)
        Y[i] = Y[i] + 0.5 * k2[i];

    // 3rd rk step
    if (!lerp3D(Y, st, sz, 3, vec, v))
        return true;
    float k3[num_dims];
    for (int i = 0; i < num_dims; i++)
        k3[i] = h * v[i];
    for (int i = 0; i < num_dims; i++)
        Y[i] = Y[i] + k3[i];

    // 4th rk step
    if (!lerp3D(Y, st, sz, 3, vec, v))
        return true;
    for (int i = 0; i < num_dims; i++)
        Y[i] = Y[i] + (k1[i] + 2.0 * (k2[i] + k3[i]) + h * v[i]) / 6.0;

    return true;
}




