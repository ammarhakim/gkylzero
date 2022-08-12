#include <gkyl_proj_maxwellian_on_basis.h>

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

struct gkyl_proj_maxwellian_on_basis {
  struct gkyl_rect_grid grid;
  int num_quad; // number of quadrature points to use in each direction
  int cdim; // Configuration-space dimension
  int pdim; // Phase-space dimension

  int num_conf_basis; // number of conf-space basis functions
  int num_phase_basis; // number of phase-space basis functions

  bool use_gpu;

  // for quadrature in phase-space
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // ordinates for quadrature
  struct gkyl_array *weights; // weights for quadrature
  struct gkyl_array *basis_at_ords; // basis functions at ordinates

  // for quadrature in conf space
  int tot_conf_quad; // total number of quadrature points
  struct gkyl_array *conf_ordinates; // conf-space ordinates for quadrature
  struct gkyl_array *conf_weights; // weights for conf-space quadrature  
  struct gkyl_array *conf_basis_at_ords; // conf-space basis functions at ordinates

  struct gkyl_array *fun_at_ords; // function (Maxwellian) evaluated at
                                  // ordinates in a cell.
  int *p2c_qidx;  // Mapping between conf-space and phase-space ordinates.
};

void
gkyl_proj_maxwellian_on_basis_lab_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *M0, const struct gkyl_array *M1i, const struct gkyl_array *M2,
  struct gkyl_array *fmax);

void
gkyl_proj_maxwellian_on_basis_prim_mom_cu(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_r, const struct gkyl_range *conf_r,
  const struct gkyl_array *m0, const struct gkyl_array *udrift, const struct gkyl_array *vtsq,
  struct gkyl_array *fmax);
