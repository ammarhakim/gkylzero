#include <cstdio>
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
#include <gkyl_bc_twistshift.h>
}

#include <cstdio>
#include <cassert>

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)


__global__ void
gkyl_bc_twistshift_inc_cu_kernel(double* ftar, int n, double* matdata)
{
   ftar[START_ID] += matdata[START_ID];
}


__global__ void
gkyl_bc_twistshift_clear_cu_kernel(double* ftar, int n)
{
   ftar[START_ID] = 0.0;
}


void
gkyl_bc_twistshift_inc_cu(double* ftar, int n, struct gkyl_mat* mat)
{
  gkyl_bc_twistshift_inc_cu_kernel<<<1,n>>>(ftar, n, mat->data);
}

void
gkyl_bc_twistshift_clear_cu(double* ftar, int n)
{
  gkyl_bc_twistshift_clear_cu_kernel<<<1,n>>>(ftar, n);
}
