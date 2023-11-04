#include <math.h>

#include <gkyl_const.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
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
surf_comp_to_phys(int dir, int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) {
    if (d < dir)
      xout[d] = 0.5*dx[d]*eta[d]+xc[d];
    else if (d > dir)
      xout[d] = 0.5*dx[d]*eta[d-1]+xc[d];
    else 
      xout[d] = 0.0;
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

static inline double eval_fpo_h(double gamma_q, double den_q, 
  double rel_speed_q, double vtsq) 
{
  double sqrt2temp_over_m_q = sqrt(2.0*vtsq);
  return gamma_q*den_q/rel_speed_q * erf(rel_speed_q/sqrt2temp_over_m_q);
}

static inline double eval_fpo_g(double gamma_q, double den_q, double rel_speed_q, 
  double vtsq) 
{
  double rel_speedsq_q = pow(rel_speed_q, 2);
  double sqrt2temp_over_m_q = sqrt(2.0*vtsq);
  return  gamma_q*den_q*sqrt2temp_over_m_q*
    (1.0/(sqrt(GKYL_PI))*exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q, 2)) + 
    erf(rel_speed_q/sqrt2temp_over_m_q)*(sqrt2temp_over_m_q/(2.0*rel_speed_q) + 
    rel_speed_q/sqrt2temp_over_m_q));
}

static inline double eval_fpo_dhdv(double gamma_q, double den_q, 
  double rel_vel_in_dir_q, double vtsq, double rel_speed_q) 
{
  double rel_speedsq_q = pow(rel_speed_q, 2);
  double sqrt2temp_over_m_q = sqrt(2.0*vtsq);
  return gamma_q*den_q*rel_vel_in_dir_q * (
    2.0*exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))/(sqrt(GKYL_PI)*sqrt2temp_over_m_q*rel_speedsq_q) -
    erf(rel_speed_q/sqrt2temp_over_m_q)/pow(rel_speedsq_q, 1.5)
  );
}

static inline double eval_fpo_dgdv(double gamma_q, double den_q, double rel_vel_in_dir_q,
  double vtsq, double rel_speed_q) {
  double rel_speedsq_q = pow(rel_speed_q,2);
  double sqrt2temp_over_m_q = sqrt(2.0*vtsq);
  return gamma_q*den_q*rel_vel_in_dir_q * (
    exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))*sqrt2temp_over_m_q/(sqrt(GKYL_PI)*rel_speedsq_q) -
    erf(rel_speed_q/sqrt2temp_over_m_q)*(pow(sqrt2temp_over_m_q,2)/(2.0*pow(rel_speed_q,3)) - 1.0/rel_speed_q)
  );
}

static inline double eval_fpo_d2gdv2(double gamma_q, double den_q, 
  double rel_vel_in_dir_q, double vtsq, double rel_speed_q) {
  double rel_speedsq_q = pow(rel_speed_q,2);
  double sqrt2temp_over_m_q = sqrt(2.0*vtsq);
  return gamma_q*den_q*(exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))/sqrt(GKYL_PI)*(sqrt2temp_over_m_q/rel_speedsq_q - 
    pow(rel_vel_in_dir_q,2)*2.0*sqrt2temp_over_m_q*(1.0/pow(rel_speedsq_q, 2) +
    1.0/(pow(sqrt2temp_over_m_q,2)*rel_speedsq_q) + (pow(sqrt2temp_over_m_q,2) - 
    2.0*rel_speedsq_q)/(2.0*pow(sqrt2temp_over_m_q,2)*pow(rel_speedsq_q, 2)))) +
    erf(rel_speed_q/sqrt2temp_over_m_q)*(pow(rel_vel_in_dir_q,2)*
    (3.0*pow(sqrt2temp_over_m_q,2)-2.0*rel_speedsq_q)/(2.0*pow(rel_speed_q, 5)) - 
    pow(sqrt2temp_over_m_q, 2)/(2.0*pow(rel_speed_q, 3)) + 1.0/rel_speed_q));
}

static inline double eval_fpo_d2gdv2_cross(double gamma_q, double den_q,
  double rel_vel_in_dir1_q, double rel_vel_in_dir2_q, double rel_speed_q,
  double vtsq) {
  double rel_speedsq_q = pow(rel_speed_q,2);
  double sqrt2temp_over_m_q = sqrt(2.0*vtsq);
  return gamma_q*den_q*rel_vel_in_dir1_q*rel_vel_in_dir2_q/pow(rel_speedsq_q,2)*(
    erf(rel_speed_q/sqrt2temp_over_m_q)*(3.0*pow(sqrt2temp_over_m_q,2)-2.0*rel_speedsq_q)/(2.0*rel_speed_q) -
    3.0*exp(-rel_speedsq_q/pow(sqrt2temp_over_m_q,2))*sqrt2temp_over_m_q/sqrt(GKYL_PI)
  );
}

struct gkyl_proj_maxwellian_pots_on_basis {
  struct gkyl_rect_grid grid;
  int cdim; // Configuration space dimension
  int pdim; // Phase space dimension
  int num_quad; // Number of 1D quadrature points

  const struct gkyl_basis *phase_basis;
  const struct gkyl_basis *conf_basis;
  const struct gkyl_basis *surf_basis;

  int num_conf_basis; // number of configuration space basis functions
  int num_phase_basis; // number of phase space basis functions
  int num_surf_basis; // Number of surface basis functions
                      
  struct gkyl_array *nodes;
  struct gkyl_array *surf_nodes;
};
