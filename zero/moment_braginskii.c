#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_braginskii_priv.h>

// Makes indexing cleaner
static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

static const unsigned P11 = 4;
static const unsigned P12 = 5;
static const unsigned P13 = 6;
static const unsigned P22 = 7;
static const unsigned P23 = 8;
static const unsigned P33 = 9;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;
static const unsigned PHIE = 6;
static const unsigned PHIM = 7;

// 1D stencil locations (L: lower, C: center, U: upper)
enum loc_1d {
  L_1D, C_1D, U_1D
};

// 2D stencil locations (L: lower, C: center, U: upper)
enum loc_2d {
  LL_2D, LC_2D, LU_2D,
  CL_2D, CC_2D, CU_2D,
  UL_2D, UC_2D, UU_2D
};

// 3D stencil locations (L: lower, C: center, U: upper)
enum loc_3d {
  LLL_3D, LLC_3D, LLU_3D,
  LCL_3D, LCC_3D, LCU_3D,
  LUL_3D, LUC_3D, LUU_3D,
  
  CLL_3D, CLC_3D, CLU_3D,
  CCL_3D, CCC_3D, CCU_3D,
  CUL_3D, CUC_3D, CUU_3D,

  ULL_3D, ULC_3D, ULU_3D,
  UCL_3D, UCC_3D, UCU_3D,
  UUL_3D, UUC_3D, UUU_3D
};

struct gkyl_moment_braginskii {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int nfluids; // number of fluids in multi-fluid system
  struct gkyl_moment_braginskii_data param[GKYL_MAX_SPECIES]; // struct of fluid parameters
  enum gkyl_braginskii_type type_brag; // which Braginskii equations (magnetized versus unmagnetized)
  double epsilon0; // permittivity of free space
};

static void
create_offsets(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

// Calculate the magnitude of the local magnetic field
static inline double
calc_mag_b(double em_tot[8])
{
  return sqrt(em_tot[BX]*em_tot[BX] + em_tot[BY]*em_tot[BY] + em_tot[BZ]*em_tot[BZ]);
}

// Calculate the cyclotron frequency based on the species' parameters
static inline double
calc_omega_c(double charge, double mass, double em_tot[8])
{
  double omega_c = 0.0;
  double Bmag = calc_mag_b(em_tot);
  if (Bmag > 0.0)
    omega_c = charge*Bmag/mass;
  return omega_c;
}

// Calculate magnetic field unit vector
static void
calc_bhat(double em_tot[8], double b[3])
{
  double Bx = em_tot[BX];
  double By = em_tot[BY];
  double Bz = em_tot[BZ];
  double Bmag = calc_mag_b(em_tot);
  // get magnetic field unit vector 
  if (Bmag > 0.0) {
    b[0] = Bx/Bmag;
    b[1] = By/Bmag;
    b[2] = Bz/Bmag;
  }  
}

// Calculate viscous stress tensor Pi (magnetized)
// In 1D, computes quantity at cell edge of two-cell interface
static void
calc_pi_1D(double dx, double fluid_l[5], double fluid_u[5], double em_tot_l[8], double em_tot_u[8], double pi[6])
{
  double rho_l = fluid_l[RHO];
  double rho_u = fluid_u[RHO];
  double u_l[3] = { fluid_l[MX]/rho_l, fluid_l[MY]/rho_l, fluid_l[MZ]/rho_l };
  double u_u[3] = { fluid_u[MX]/rho_u, fluid_u[MY]/rho_u, fluid_u[MZ]/rho_u };
  double p_l = fluid_l[ER] - (fluid_l[MX]*fluid_l[MX] + fluid_l[MY]*fluid_l[MY] + fluid_l[MZ]*fluid_l[MZ])/rho_l;
  double p_u = fluid_u[ER] - (fluid_u[MX]*fluid_u[MX] + fluid_u[MY]*fluid_u[MY] + fluid_u[MZ]*fluid_u[MZ])/rho_u;
  
  double w[6];
  calc_ros_1D(dx, u_l, u_u, w);

  double b_l[3] = { 0.0, 0.0, 0.0 };
  double b_u[3] = { 0.0, 0.0, 0.0 };
  calc_bhat(em_tot_l, b_l);
  calc_bhat(em_tot_u, b_u);
  
  double b_avg[3] = { 0.0, 0.0, 0.0 };

  b_avg[0] = calc_arithm_avg_1D(b_l[0], b_u[0]);
  b_avg[1] = calc_arithm_avg_1D(b_l[1], b_u[1]);
  b_avg[2] = calc_arithm_avg_1D(b_l[2], b_u[2]);

  // pi_parallel = (bb - 1/3 I) (bb - 1/3 I) : W
  
  // parallel rate of strain = (bb - 1/3 I) : W
  double par_ros = (b_avg[0]*b_avg[0] - 1.0/3.0)*w[0] + 2.0*b_avg[0]*b_avg[1]*w[1] + 2.0*b_avg[0]*b_avg[2]*w[2] 
                  + (b_avg[1]*b_avg[1] - 1.0/3.0)*w[3] + 2.0*b_avg[1]*b_avg[2]*w[4] 
                  + (b_avg[2]*b_avg[2] - 1.0/3.0)*w[5];

  double pi_par[6];
  pi_par[0] = (b_avg[0]*b_avg[0] - 1.0/3.0)*par_ros;
  pi_par[1] = b_avg[0]*b_avg[1]*par_ros;
  pi_par[2] = b_avg[0]*b_avg[2]*par_ros;
  pi_par[3] = (b_avg[1]*b_avg[1] - 1.0/3.0)*par_ros;
  pi_par[4] = b_avg[1]*b_avg[2]*par_ros;
  pi_par[5] = (b_avg[2]*b_avg[2] - 1.0/3.0)*par_ros;

  // pi_perp = (I - bb) . W . (I + 3bb) + (I + 3bb) . W . (I - bb)
  
  // (b . W . I)_x = b_x W_xx + b_y W_xy + b_z W_xz
  double bWIx = w[0]*b_avg[0] + w[1]*b_avg[1] + w[2]*b_avg[2];
  // (b . W . I)_y = b_x W_xy + b_y W_yy + b_z W_yz  
  double bWIy = w[1]*b_avg[0] + w[3]*b_avg[1] + w[4]*b_avg[2];
  // (b . W . I)_z = b_x W_xz + b_y W_yz + b_z W_zz
  double bWIz = w[0]*b_avg[2] + w[4]*b_avg[1] + w[5]*b_avg[2];

  // b . W . b
  double bWb = b_avg[0]*b_avg[0]*w[0] + 2.0*b_avg[0]*b_avg[1]*w[1] + 2.0*b_avg[0]*b_avg[2]*w[2] 
              + b_avg[1]*b_avg[1]*w[3] + 2.0*b_avg[1]*b_avg[2]*w[4] 
              + b_avg[2]*b_avg[2]*w[5];
  
  double pi_perp[6];
  pi_perp[0] = 2.0*(w[0] + 2.0*b_avg[0]*bWIx - 3.0*b_avg[0]*b_avg[0]*bWb);
  pi_perp[1] = 2.0*w[1] + 2.0*(b_avg[1]*bWIx + b_avg[0]*bWIy) - 6.0*b_avg[0]*b_avg[1]*bWb;
  pi_perp[2] = 2.0*w[2] + 2.0*(b_avg[2]*bWIx + b_avg[0]*bWIz) - 6.0*b_avg[0]*b_avg[2]*bWb;
  pi_perp[3] = 2.0*(w[3] + 2.0*b_avg[1]*bWIy - 3.0*b_avg[1]*b_avg[1]*bWb);
  pi_perp[4] = 2.0*w[4] + 2.0*(b_avg[2]*bWIy + b_avg[1]*bWIz) - 6.0*b_avg[1]*b_avg[2]*bWb;
  pi_perp[5] = 2.0*(w[5] + 2.0*b_avg[2]*bWIz - 3.0*b_avg[2]*b_avg[2]*bWb);
}

// Compute the RHS for the Braginskii transport terms in 1D
// If isothermal Euler, only update momentum
// If Euler, update pressure due to viscous heating and heat conduction
// If Ten moment, update momentum *only* due to inter-species collisions
//
// Logic is thus as follows: first compute RHS due to inter-species collisions,
// then viscosity in momentum equation if equation system supports viscosity. 
// Then add the corresponding contribution to the pressure variable from heating
// due to the momentum evolution (Joule heating, viscous heating). 
// Finally add thermal conduction and temperature equilibriation. 
static void
mag_braginskii_update(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES], const double *em_tot_d[],
  double *cflrate, double *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids; 
}

static void
unmag_braginskii_update(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double *cflrate, double *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids; 
}

gkyl_moment_braginskii*
gkyl_moment_braginskii_new(struct gkyl_moment_braginskii_inp inp)
{
  gkyl_moment_braginskii *up = gkyl_malloc(sizeof(gkyl_moment_braginskii));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->nfluids = inp.nfluids;
  for (int n=0; n<inp.nfluids; ++n) up->param[n] = inp.param[n];
  up->type_brag = inp.type_brag;
  up->epsilon0 = inp.epsilon0;

  return up;
}

void
gkyl_moment_braginskii_advance(const gkyl_moment_braginskii *bes, struct gkyl_range update_range,
  const struct gkyl_array *fluid[], const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *rhs[])
{
  int nfluids = bes->nfluids;
  int ndim = update_range.ndim;
  long sz[] = { 3, 9, 27 };

  long offsets[sz[ndim-1]];
  create_offsets(&update_range, offsets);

  const double* fluid_d[sz[ndim-1]][GKYL_MAX_SPECIES];
  const double* em_tot_d[sz[ndim-1]];
  double *rhs_d[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_range);
  while (gkyl_range_iter_next(&iter)) {
    
    long linc = gkyl_range_idx(&update_range, iter.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc + offsets[i]); 
      for (int n=0; n<nfluids; ++n)
        fluid_d[i][n] = gkyl_array_cfetch(fluid[n], linc + offsets[i]);
    }
    for (int n=0; n<nfluids; ++n)
      rhs_d[n] = gkyl_array_fetch(rhs[n], linc);

    if (bes->type_brag == GKYL_MAG_BRAG)
      mag_braginskii_update(bes, fluid_d, em_tot_d, gkyl_array_fetch(cflrate, linc), rhs_d);
    else if (bes->type_brag == GKYL_UNMAG_BRAG)
      unmag_braginskii_update(bes, fluid_d, gkyl_array_fetch(cflrate, linc), rhs_d);
  }
}

void
gkyl_moment_braginskii_release(gkyl_moment_braginskii* up)
{
  free(up);
}
