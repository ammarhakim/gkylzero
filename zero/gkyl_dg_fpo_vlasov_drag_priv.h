#pragma once

// Private header, not for direct use in user code
#include <gkyl_fpo_vlasov_kernels.h>
#include <string.h>

// Types for various kernels
typedef void (*fpo_vlasov_drag_surf_t)(const double *dxv, 
  const double* drag_coeff_stencil[9], const double* f_stencil[9], double* out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_vol_kern_list;
typedef struct { fpo_vlasov_drag_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_drag_surf_kern_list;

struct dg_fpo_vlasov_drag {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  fpo_vlasov_drag_surf_t surf[3][3]; // Generic surface term domain stencil
  struct gkyl_range phase_range; // Configuration space range.
  struct gkyl_dg_fpo_vlasov_drag_auxfields auxfields; // Auxiliary fields.
};

// Serendipity volume kernels
// Only 1x3v for now

GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  
  return fpo_vlasov_drag_vol_1x3v_ser_p1( dx,
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff, pidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_fpo_vlasov_drag_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, const int* idx, const double *qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  
  long pidx = gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx);
  
  return fpo_vlasov_drag_vol_1x3v_ser_p2( dx,
    (const double*) gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff, pidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_vol_kern_list ser_vol_kernels[] = {
  { kernel_fpo_vlasov_drag_vol_1x3v_ser_p1, kernel_fpo_vlasov_drag_vol_1x3v_ser_p2, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernels are handled generally using the domain decomp stencil.
// These lists contain both surface and boundary surface
// kernels, and the proper kernel is selected with the index returned from idx_to_inloup_ker 
// in hyper_dg. 

// vx-direction surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_1x3v_vx_kernels[] = {
  { fpo_vlasov_drag_surfvx_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1_lovx, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p1_upvx },
  {  fpo_vlasov_drag_surfvx_1x3v_ser_p2, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p2_lovx, fpo_vlasov_drag_boundary_surfvx_1x3v_ser_p2_upvx},
  { NULL, NULL, NULL },
};

// vy-direction surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_1x3v_vy_kernels[] = {
  { fpo_vlasov_drag_surfvy_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p1_lovy, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p1_upvy },
  {  fpo_vlasov_drag_surfvy_1x3v_ser_p2, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p2_lovy, fpo_vlasov_drag_boundary_surfvy_1x3v_ser_p2_upvy},
  { NULL, NULL, NULL },
};

// vz-direction surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list ser_surf_1x3v_vz_kernels[] = {
  { fpo_vlasov_drag_surfvz_1x3v_ser_p1, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p1_lovz, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p1_upvz },
  {  fpo_vlasov_drag_surfvz_1x3v_ser_p2, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p2_lovz, fpo_vlasov_drag_boundary_surfvz_1x3v_ser_p2_upvz},
  { NULL, NULL, NULL },
};

// Lists of domain stencil kernel lists, ordered by cdim.
GKYL_CU_DH
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list *ser_surf_vx_kernels[3] = {
  ser_surf_1x3v_vx_kernels,
  NULL,
  NULL
};

GKYL_CU_DH
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list *ser_surf_vy_kernels[3] = {
  ser_surf_1x3v_vy_kernels,
  NULL,
  NULL
};

GKYL_CU_DH
static const gkyl_dg_fpo_vlasov_drag_surf_kern_list *ser_surf_vz_kernels[3] = {
  ser_surf_1x3v_vz_kernels,
  NULL,
  NULL
};

GKYL_CU_D
static void
choose_fpo_vlasov_drag_surf_kern(fpo_vlasov_drag_surf_t surf[3],
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list** surf_kernel_list, int cdim, int poly_order)
{
  memcpy(surf, surf_kernel_list[cdim]->kernels[poly_order], 3 * sizeof(fpo_vlasov_drag_surf_t));
}

/* Free fpo_vlasov_diff equation object
 *
 * @param ref Reference counter for constant fpo_vlasov_diff equation
*/
void gkyl_fpo_vlasov_drag_free(const struct gkyl_ref_count* ref);

// Gen surface term called by hyper_dg_gen_stencil_advance
GKYL_CU_D
static void
fpo_drag_gen_surf_term(const struct gkyl_dg_eqn* eqn, int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  int keri, const int idx[9][GKYL_MAX_DIM], const double* qIn[9],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  long sz_dim = 9;
  int cdim = fpo_vlasov_drag->cdim;
  const double* drag_coeff_d[9];
  for (int i=0; i<sz_dim; ++i) {
    if (idx[i]) {
      drag_coeff_d[i] = gkyl_array_cfetch(fpo_vlasov_drag->auxfields.drag_coeff,
          gkyl_range_idx(&fpo_vlasov_drag->phase_range, idx[i]));
    }
  }

  if (dir1 >= cdim && dir1 == dir2) {
    fpo_vlasov_drag->surf[dir1-cdim][keri](dxc, drag_coeff_d, qIn, qRhsOut);
  }
}
