#pragma once

#include <gkyl_maxwell_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in Maxwell DG equation object creation
// functions

// context for use in Wall BCs
struct maxwell_wall_bc_ctx {
  int dir; // direction for BCs
  const struct gkyl_basis *basis; // basis function
};

enum { M_EX, M_EY, M_EZ, M_BX, M_BY, M_BZ }; // components of EM field
static const int m_flip_even[3][3] = { // zero tangent E and zero normal B
  {M_BX, M_EY, M_EZ},
  {M_BY, M_EX, M_EZ},
  {M_BZ, M_EX, M_EY},
};
static const int m_flip_odd[3][3] = { // zero gradient
  { M_EX, M_BY, M_BZ },
  { M_EY, M_BX, M_BZ },
  { M_EZ, M_BX, M_BY },
};

// Types for various kernels
typedef double (*maxwell_vol_t)(const gkyl_maxwell_inp *meq,
  const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out);

typedef void (*maxwell_surf_t)(const gkyl_maxwell_inp *meq, const double *w, const double *dx,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

typedef void (*bc_funcf_t)(size_t nc, double *out, const double *inp, void *ctx);

// for use in kernel tables
typedef struct { maxwell_vol_t kernels[3]; } gkyl_dg_maxwell_vol_kern_list;
typedef struct { maxwell_surf_t kernels[3]; } gkyl_dg_maxwell_surf_kern_list;

//
// Serendipity basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ser_vol_kernels[] = {
  { NULL, maxwell_vol_1x_ser_p1, maxwell_vol_1x_ser_p2 }, // 0
  { NULL, maxwell_vol_2x_ser_p1, maxwell_vol_2x_ser_p2 }, // 1
  { NULL, maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_ser_p2 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_ser_p2 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, NULL }, // 2
};

//
// Tensor-product basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ten_vol_kernels[] = {
  { NULL, maxwell_vol_1x_ser_p1, maxwell_vol_1x_tensor_p2 }, // 0
  { NULL, maxwell_vol_2x_ser_p1, maxwell_vol_2x_tensor_p2 }, // 1
  { NULL, maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_tensor_p2 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_tensor_p2 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_tensor_p2 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, NULL }, // 2
};

struct dg_maxwell {
  struct gkyl_dg_eqn eqn; // Base object    
  gkyl_maxwell_inp maxwell_data; // Parameters needed by kernels
  maxwell_vol_t vol; // pointer to volume kernel
  maxwell_surf_t surf[3]; // pointers to surface kernels
  bc_funcf_t wall_bc; // wall BCs function
};


/**
 * Free Maxwell equation object
 *
 * @param ref Reference counter for Maxwell equation
 */
void gkyl_maxwell_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx, 
  const int*  idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell->vol(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  maxwell->surf[dir](&maxwell->maxwell_data, xcC, dxC,
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  
}

static void
maxwell_wall_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct maxwell_wall_bc_ctx *mc = ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  const int *feven = m_flip_even[dir];
  const int *fodd = m_flip_odd[dir];

  for (int i=0; i<3; ++i) {
    int eloc = nbasis*feven[i], oloc = nbasis*fodd[i];
    mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
    mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
  }
  // correction potentials
  int eloc = nbasis*6, oloc = nbasis*7;
  mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
  mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
}