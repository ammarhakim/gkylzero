#include "gkyl_rect_grid.h"
#include <gkyl_proj_maxwellian_pots_on_basis.h>

GKYL_CU_DH
static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
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

  int num_conf_basis; // number of configuration space basis functions
  int num_phase_basis; // number of phase space basis functions

  struct gkyl_range conf_qrange; // Range of configuration space ordinates
  struct gkyl_range phase_qrange; // Range of phase space ordinates

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

  struct gkyl_array *fpo_h_at_ords; // Potential H at ordinates
  struct gkyl_array *fpo_g_at_ords;
};
