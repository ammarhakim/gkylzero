#pragma once

// Private header for bc_twistshift updater, not for direct use in user code.

#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_gyrokinetic_kernels.h>
#include <assert.h>
#include <gkyl_mat.h>

// Function pointer type for twistshift kernels.
typedef void (*twistshift_xlimdg_t)(double sFac, const double *xLimLo,
  const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

typedef void (*twistshift_ylimdg_t)(double sFac, double xLimLo,
  double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

typedef void (*twistshift_fullcell_t)(double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *tsmat);

typedef struct { twistshift_xlimdg_t kernels[3]; }   twistshift_xlimdg_kern_list;  // For use in kernel tables.
typedef struct { twistshift_ylimdg_t kernels[3]; }   twistshift_ylimdg_kern_list;  // For use in kernel tables.
typedef struct { twistshift_fullcell_t kernels[3]; } twistshift_fullcell_kern_list;  // For use in kernel tables.

// Serendipity  kernels.
static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_0v[] = {
  {NULL, twistshift_xlimdg_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_xlimdg_3x_ser_p1_yshift_p1, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_0v[] = {
  {NULL, twistshift_ylimdg_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_ylimdg_3x_ser_p1_yshift_p1, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_0v[] = {
  {NULL, twistshift_fullcell_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_fullcell_3x_ser_p1_yshift_p1, NULL,},
};

static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_2v[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_xlimdg_3x2v_ser_p1_yshift_p1, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_2v[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_ylimdg_3x2v_ser_p1_yshift_p1, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_2v[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_fullcell_3x2v_ser_p1_yshift_p1, NULL,},
};

struct gkyl_bc_twistshift_kernels {
  twistshift_xlimdg_t xlimdg;
  twistshift_ylimdg_t ylimdg;
  twistshift_fullcell_t fullcell;
};

// Primary struct in this updater.
struct gkyl_bc_twistshift {
  int dir;
  int do_dir;
  int shift_dir;
  enum gkyl_edge_loc edge;
  const struct gkyl_basis *basis;
  bool use_gpu;
  struct gkyl_bc_twistshift_kernels *kernels;  // kernels.
  struct gkyl_bc_twistshift_kernels *kernels_cu;  // device copy.
  const struct gkyl_rect_grid *grid;
  const int *ndonors;
  int donor_factor; // total number of donors is sum(ndonors)*donor_factor
  int *ndonors_cum_cu;
  const int *cells_do; // y indices of donor cells for each x and y
  int *remDir;
  int *locDir;
  int *remDir_do;
  int *locDir_do;
  long *locs;
  long *locs_cu;
  long *tar_locs;
  long *tar_locs_cu;
  const struct gkyl_range *local_range_ext;
  const struct gkyl_range *local_range_update;
  struct gkyl_range *yrange;
  struct gkyl_range *xrange;
  struct gkyl_nmat *matsdo;
  struct gkyl_nmat *matsdo_ho;
  struct gkyl_nmat *vecsdo;
  struct gkyl_nmat *vecstar;
};

void gkyl_bc_twistshift_choose_kernels_cu(const struct gkyl_basis *basis, int cdim,
  struct gkyl_bc_twistshift_kernels *kers)
{
  int dim = basis->ndim;
  int vdim = dim - cdim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->xlimdg   = vdim==0?   ser_twistshift_xlimdg_list_0v[cdim-2].kernels[poly_order] :   ser_twistshift_xlimdg_list_2v[cdim-2].kernels[poly_order];
      kers->ylimdg   = vdim==0?   ser_twistshift_ylimdg_list_0v[cdim-2].kernels[poly_order] :   ser_twistshift_ylimdg_list_2v[cdim-2].kernels[poly_order];
      kers->fullcell = vdim==0? ser_twistshift_fullcell_list_0v[cdim-2].kernels[poly_order] : ser_twistshift_fullcell_list_2v[cdim-2].kernels[poly_order];
      return;
    default:
      assert(false);
      break;
  }
}

void gkyl_bc_twistshift_choose_kernels(const struct gkyl_basis *basis, int cdim,
  struct gkyl_bc_twistshift_kernels *kers)
{
  int dim = basis->ndim;
  int vdim = dim - cdim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->xlimdg   = vdim==0?   ser_twistshift_xlimdg_list_0v[cdim-2].kernels[poly_order] :   ser_twistshift_xlimdg_list_2v[cdim-2].kernels[poly_order];
      kers->ylimdg   = vdim==0?   ser_twistshift_ylimdg_list_0v[cdim-2].kernels[poly_order] :   ser_twistshift_ylimdg_list_2v[cdim-2].kernels[poly_order];
      kers->fullcell = vdim==0? ser_twistshift_fullcell_list_0v[cdim-2].kernels[poly_order] : ser_twistshift_fullcell_list_2v[cdim-2].kernels[poly_order];
      return;
    default:
      assert(false);
      break;
  }
}
