#include <gkyl_util.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_maxwellian_pots_on_basis.h>

GKYL_CU_DH
static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

GKYL_CU_DH
static inline void
edge_idx_to_phase_idx(int ndim, int dir, const int *surf_idx, int edge_idx, int *phase_idx) {
  phase_idx[dir] = edge_idx;
  for (int i=0; i<ndim; ++i) {
    if (i < dir) phase_idx[i] = surf_idx[i];
    else if (i == dir) phase_idx[i] = edge_idx;
    else phase_idx[i] = surf_idx[i-1];
  }
}
static inline void
copy_idx_arrays(int cdim, int pdim, const int *cidx, const int *vidx, int *out)
{
  for (int i=0; i<cdim; ++i)
    out[i] = cidx[i];
  for (int i=cdim; i<pdim; ++i)
    out[i] = vidx[i-cdim];
}


struct gkyl_proj_maxwellian_pots_on_basis {
  struct gkyl_rect_grid grid;
  int cdim; // Configuration space dimension
  int pdim; // Phase space dimension
  int num_quad; // Number of 1D quadrature points

  int num_conf_basis; // number of configuration space basis functions
  int num_phase_basis; // number of phase space basis functions
  int num_surf_basis; // Number of surface basis functions

  struct gkyl_range conf_qrange; // Range of configuration space ordinates
  struct gkyl_range phase_qrange; // Range of phase space ordinates
  struct gkyl_range surf_qrange;

  // quadrature in phase space
  int tot_quad;
  struct gkyl_array *ordinates;
  struct gkyl_array *weights;
  struct gkyl_array *basis_at_ords;

  // quadrature in configuration space
  int tot_conf_quad;
  struct gkyl_array *conf_ordinates; // configuration space ordinates
  struct gkyl_array *conf_weights; // weights for configuration space quadrature
  struct gkyl_array *conf_basis_at_ords; // configuration space basis functions at ordinates

  // Surface expansion quadrature
  int tot_surf_quad;
  struct gkyl_array *surf_ordinates; // Surface ordinates
  struct gkyl_array *surf_weights; // Weights for surface quadrature
  struct gkyl_array *surf_basis_at_ords; // Surface basis functions at ordinates

  struct gkyl_array *fpo_h_at_ords; // Potential H at ordinates
  struct gkyl_array *fpo_g_at_ords; // Potential G at ordinates

  // Arrays for quantities defined at boundary surface quadrature points
  struct gkyl_array *fpo_h_at_surf_ords;
  struct gkyl_array *fpo_g_at_surf_ords;
  struct gkyl_array *fpo_dhdv_at_surf_ords;
  struct gkyl_array *fpo_dgdv_at_surf_ords;
  struct gkyl_array *fpo_d2gdv2_at_surf_ords;
};
