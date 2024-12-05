#pragma once

// Private header for array_average updater, not for direct use by user.

#include <gkyl_util.h>
#include <math.h>
#include <gkyl_array_average.h>
#include <gkyl_array_average_kernels.h>
#include <assert.h>
#include <gkyl_dg_bin_ops.h>

// Function pointer type for array_average kernels.
typedef void (*array_average_t)( const double subvol, const double *win, const double *fin, double *out);

// For use in kernel tables.
typedef struct { array_average_t kernels[2]; } array_average_kern_list;
typedef struct { array_average_kern_list list[7]; } dim_array_average_kern_list;

GKYL_CU_D
static const dim_array_average_kern_list gkyl_array_average_ker_list[] = {
  { // Kernel list for 1x integration
    .list = 
      {
        {gkyl_array_average_1x_ser_p1_ker, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
      }
  },
  { // Kernel list for 2x integration
    .list = 
      {
        {gkyl_array_average_2x_ser_p1_ker, NULL},
        {gkyl_array_average_2x_ser_p1_x_ker, NULL},
        {gkyl_array_average_2x_ser_p1_y_ker, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
      }
  },
  { // Kernel list for 3x integration
    .list = 
      {
        {gkyl_array_average_3x_ser_p1_ker, NULL},
        {gkyl_array_average_3x_ser_p1_x_ker, NULL},
        {gkyl_array_average_3x_ser_p1_y_ker, NULL},
        {gkyl_array_average_3x_ser_p1_z_ker, NULL},
        {gkyl_array_average_3x_ser_p1_xy_ker, NULL},
        {gkyl_array_average_3x_ser_p1_xz_ker, NULL},
        {gkyl_array_average_3x_ser_p1_yz_ker, NULL},
      }
  }
};

// Primary struct in this updater.
struct gkyl_array_average {
  // dimensionality of the full array
  int ndim;
  struct gkyl_basis tot_basis;
  struct gkyl_basis sub_basis;
  // the updater stores the ranges of the input and output
  struct gkyl_range tot_rng;
  struct gkyl_range sub_rng;

  // if we use gpu or not
  bool use_gpu;

  // array that indicates if the dimension is averaged or not
  int isavg_dim[GKYL_MAX_CDIM]; 

  // array that maps the sub dimensions to the full one
  // examples:
  // - for 3x op_yz, then sub_dir = {1,2},
  // - for 2x op_x, then sub_dir = {0}
  int sub_dir[GKYL_MAX_CDIM]; 

  // Single-cell sub-volume element (length for 1D avg, area for 2D)
  double subvol;

  // Single cell average kernel.
  array_average_t kernel;  

  // Pointer to itself on device.
  struct gkyl_array_average *on_dev;

  // weighted integral
  bool isweighted;
  struct gkyl_array *weights;
  struct gkyl_array *integral_weights;

  // memory for the weak division at the end of averaging
  gkyl_dg_bin_op_mem *div_mem;
};

GKYL_CU_D static
void gkyl_array_average_choose_kernel(struct gkyl_array_average *up, const struct gkyl_basis *basis,  enum gkyl_array_average_op op)
{
  int ndim = basis->ndim, poly_order = basis->poly_order;

  up->kernel = gkyl_array_average_ker_list[ndim-1].list[op].kernels[poly_order-1];

}

struct gkyl_array_average*
gkyl_array_average_cu_dev_new(const struct gkyl_basis *basis,
  const struct gkyl_array *GKYL_RESTRICT weights, enum gkyl_array_average_op op);

void gkyl_array_average_advance_cu(gkyl_array_average *up, const struct gkyl_array *arr,
  const struct gkyl_range *range, double *out);
