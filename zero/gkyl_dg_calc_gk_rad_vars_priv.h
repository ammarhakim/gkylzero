#include <math.h>
#include <string.h>

#include <gkyl_const.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <gkyl_rect_grid.h>

struct gkyl_dg_calc_gk_rad_vars {
  struct gkyl_rect_grid phase_grid;
  int cdim; // Configuration space dimension
  int pdim; // Phase space dimension

  const struct gkyl_basis *phase_basis;
  const struct gkyl_basis *conf_basis;

  int num_conf_basis; // number of configuration space basis functions
  int num_phase_basis; // number of phase space basis functions
   
  struct gkyl_array *nodes;

  double charge, mass;
  double a, alpha, beta, gamma, v0; // Fitting parameters for radiation drag coefficient
  const struct gk_geometry *gk_geom; // Pointer to geometry struct
};

GKYL_CU_DH
static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static void
nod2mod(int num_ret_vals, const struct gkyl_basis *basis, const struct gkyl_array *fun_at_nodes, double *f) {
  const double *fao = gkyl_array_cfetch(fun_at_nodes, 0);

  int num_basis = basis->num_basis;
  double fnodal[num_basis];
  for (int i=0; i<num_ret_vals; ++i) {
    for (int k=0; k<num_basis; ++k) {
      fnodal[k] = fao[num_ret_vals*k+i];
    }

    basis->nodal_to_modal(fnodal, &f[num_basis*i]);
  }
}

static inline double 
eval_vnu(double charge, double mass, 
  double a, double alpha, double beta, double gamma, double v0, 
  double vpar, double mu, double bmag) 
{
  double scaled_v0 = v0/sqrt(mass/(2.0*fabs(charge)));
  double c_const = 8.0*sqrt(M_PI)*pow(fabs(charge),5.0/2.0)/mass;
  double const_mult = a*(alpha+beta)/c_const;
  double vmag = sqrt(vpar*vpar+2.0*bmag*mu/mass);
  double nu = 0.0;
  if (vmag == 0.0) {
    return 0.0;
  } else {
    nu = a*(alpha+beta)*pow(vmag,gamma)/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta));
    return 2.0/(M_PI*c_const)*nu*vmag;
  }
}

static inline double 
eval_vsqnu(double charge, double mass, 
  double vpar, double mu, double bmag, double vnu) 
{
  double vmag = sqrt(vpar*vpar+2.0*bmag*mu/mass);
  if (vmag == 0.0) {
    return 0.0;
  } else {
    return vmag*sqrt(mu)*M_PI*pow(mass/(2.0*bmag),3.0/2.0)*vnu;
  }
}
