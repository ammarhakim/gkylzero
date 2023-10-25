#pragma once

// Private header for array_integrate updater, not for direct use by user.

#include <gkyl_util.h>
#include <math.h>
#include <gkyl_array_integrate.h>
#include <assert.h>

GKYL_CU_DH void gkyl_array_integrate_op_none(double vol, int num_comp, int num_basis, const double *arr, double *out)
{
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += arr[c*num_basis]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs(double vol, int num_comp, int num_basis, const double *arr, double *out)
{
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(arr[c*num_basis])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq(double vol, int num_comp, int num_basis, const double *arr, double *out)
{
  for (unsigned c=0; c<num_comp; ++c) {
    for (unsigned b=0; b<num_basis; ++b)
      out[c] += pow(arr[c*num_basis+b],2)*vol;
  }
}

// Function pointer type for array_integrate kernels.
typedef void (*array_integrate_t)(double vol, int num_comp, int num_basis,
  const double *arr, double *out);

typedef struct { array_integrate_t kernels[3]; } array_integrate_kern_list;  // For use in kernel tables.

GKYL_CU_D
static const array_integrate_kern_list gkyl_array_integrate_ker_list = {
  .kernels = {
    gkyl_array_integrate_op_none,
    gkyl_array_integrate_op_abs,
    gkyl_array_integrate_op_sq,
  }
};

struct gkyl_array_integrate_kernels {
  array_integrate_t intker;  // single cell integration kernel.
};

// Primary struct in this updater.
struct gkyl_array_integrate {
  int num_basis, num_comp;
  bool use_gpu;
  double vol;  // Single-cell volume factor.
  struct gkyl_array_integrate_kernels *kernels;  // Single-cell integration kernel.
  struct gkyl_array_integrate_kernels *kernels_cu;  // device copy.
};

void gkyl_array_integrate_choose_kernel_cu(enum gkyl_array_integrate_op op, struct gkyl_array_integrate_kernels *kers);

void gkyl_array_integrate_choose_kernel(enum gkyl_array_integrate_op op, struct gkyl_array_integrate_kernels *kers)
{
  kers->intker = gkyl_array_integrate_ker_list.kernels[op];
}
