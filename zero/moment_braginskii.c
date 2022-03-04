#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_braginskii_priv.h>
#include <gkyl_prim_euler.h>

// 1D stencil locations (L: lower, U: upper)
enum loc_1d {
  L_1D, U_1D
};

// 2D stencil locations (L: lower, U: upper)
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

// 3D stencil locations (L: lower, U: upper)
enum loc_3d {
  LLL_3D, LLU_3D,
  LUL_3D, LUU_3D,
  ULL_3D, ULU_3D,
  UUL_3D, UUU_3D
};

struct gkyl_moment_braginskii {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int nfluids; // number of fluids in multi-fluid system
  struct gkyl_moment_braginskii_data param[GKYL_MAX_SPECIES]; // struct of fluid parameters
  enum gkyl_braginskii_type type_brag; // which Braginskii equations (magnetized versus unmagnetized)
  double epsilon0; // permittivity of free space
  double coll_fac; // constant multiplicative factor for collision time to increase or decrease collisionality
};

static void
create_offsets_vertices(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

// Fetch input quantities and compute derived quantities for magnetized Braginskii
static void
mag_var_setup(const gkyl_moment_braginskii *bes, int start, int end,
  const double *fluid_d[][GKYL_MAX_SPECIES], const double *em_tot_d[],
  double u[][2][3], double b[][3],
  double T[][2], double tau[][2],
  double eta_par[][2], double eta_perp[][2], double eta_cross[][2],
  double kappa_par[][2], double kappa_perp[][2], double kappa_cross[][2],
  double current_par[][3], double current_cross[][3])
{
  int nfluids = bes->nfluids;

  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };
  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double rho[2] = {0.0}; // Mass density for each species
  double p[2] = {0.0}; // Pressure for each species
  double omega_c[2] = {0.0}; // Cyclotron frequency for each species

  for (int j = start; j <= end; ++j) {
    // Input density and flow in each cell
    rho[ELC] = fluid_d[j][ELC][RHO];
    rho[ION] = fluid_d[j][ION][RHO];
    for (int k = MX; k <= MZ; ++k) {
      u[j][ELC][k - MX] = fluid_d[j][ELC][k] / rho[ELC];
      u[j][ION][k - MX] = fluid_d[j][ION][k] / rho[ION];
    }
    // Input magnetic field in each cell
    calc_bhat(em_tot_d[j], b[j]);
    // Derived cyclotron frequency and collision time
    omega_c[ELC] = calc_omega_c(q[ELC], m[ELC], em_tot_d[j]);
    omega_c[ION] = calc_omega_c(q[ION], m[ION], em_tot_d[j]);
    // Pressure information is different for each equation type
    for (int n = 0; n < nfluids; ++n) {
      if (bes->param[n].type_eqn == GKYL_EQN_EULER)
        p[n] = gkyl_euler_pressure(bes->param[n].p_fac, fluid_d[j][n]); // Euler needs to divide out gas_gamma factor to obtain pressure
      else if (bes->param[n].type_eqn == GKYL_EQN_ISO_EULER)
        p[n] = rho[n]*bes->param[n].p_fac*bes->param[n].p_fac; // isothermal Euler input is vth, pressure = rho*vth^2
    }
    T[j][ELC] = m[ELC] * p[ELC] / rho[ELC];
    T[j][ION] = m[ION] * p[ION] / rho[ION];

    tau[j][ELC] = calc_tau(1.0, bes->coll_fac, bes->epsilon0, q[ELC], q[ION], m[ELC], m[ION], rho[ION], T[j][ELC]);
    tau[j][ION] = calc_tau(1.0, sqrt(2.0)*bes->coll_fac, bes->epsilon0, q[ION], q[ION], m[ION], m[ION], rho[ION], T[j][ION]);

    eta_par[j][ELC] = 1.5*0.73*p[ELC]*tau[j][ELC];
    eta_perp[j][ELC] = 0.51*p[ELC]/(tau[j][ELC]*omega_c[ELC]*omega_c[ELC]);
    eta_cross[j][ELC] = 0.25*p[ELC]/omega_c[ELC];
    kappa_par[j][ELC] = 3.16*p[ELC]*tau[j][ELC]/m[ELC];
    kappa_perp[j][ELC] = 4.66*p[ELC]/(m[ELC]*tau[j][ELC]*omega_c[ELC]*omega_c[ELC]);
    kappa_cross[j][ELC] = 2.5*p[ELC]/(m[ELC]*omega_c[ELC]);

    eta_par[j][ION] = 1.5*0.96*p[ION]*tau[j][ION];
    eta_perp[j][ION] = 0.3*p[ION]/(tau[j][ION]*omega_c[ION]*omega_c[ION]);
    eta_cross[j][ION] = 0.25*p[ION]/omega_c[ION];
    kappa_par[j][ION] = 3.91*p[ION]*tau[j][ION]/m[ION];
    kappa_perp[j][ION] = 2*p[ION]/(m[ION]*tau[j][ION]*omega_c[ION]*omega_c[ION]);
    kappa_cross[j][ION] = 2.5*p[ION]/(m[ION]*omega_c[ION]);

    double thermal_par = -0.71*p[ELC]; // Parallel thermal force coefficient (same for each species)
    double thermal_perp = 1.5*p[ELC]/(omega_c[ELC]*tau[j][ELC]); // Perpendicular thermal force coefficient (same for each species)
    double b_dot_j = b[j][0]*(u[j][ION][0] - u[j][ELC][0]) + b[j][1]*(u[j][ION][1] - u[j][ELC][1]) + b[j][2]*(u[j][ION][2] - u[j][ELC][2]);
    current_par[j][0] = thermal_par*b[j][0]*b_dot_j;
    current_par[j][1] = thermal_par*b[j][1]*b_dot_j;
    current_par[j][2] = thermal_par*b[j][2]*b_dot_j;
    current_cross[j][0] = thermal_perp*(b[j][1]*(u[j][ION][2] - u[j][ELC][2]) - b[j][2]*(u[j][ION][1] - u[j][ELC][1]));
    current_cross[j][1] = thermal_perp*(b[j][2]*(u[j][ION][0] - u[j][ELC][0]) - b[j][0]*(u[j][ION][2] - u[j][ELC][2]));
    current_cross[j][2] = thermal_perp*(b[j][0]*(u[j][ION][1] - u[j][ELC][1]) - b[j][1]*(u[j][ION][0] - u[j][ELC][0]));
  }
}

static void
mag_brag_calc_vars(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES], const double *em_tot_d[],
  double *cflrate, double *brag_d[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  const int ndim = bes->ndim;
  // Allocate some memory on stack regardless of dimensionality or equation type (Euler vs. isothermal Euler)
  // These allocations allow us to make final construction of Braginskii variables dimensionally independent
  // Note: Braginskii implementation currently assumes only two fluids (electrons and an ion species)
  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };
  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double b_avg[3] = {0.0};
  double eta_par_avg[2], eta_perp_avg[2], eta_cross_avg[2] = {0.0};
  double w[2][6] = {0.0};

  double pi_par[2][6] = {0.0};
  double pi_perp[2][6] = {0.0};
  double pi_cross[2][6] = {0.0};

  double kappa_par_avg[2], kappa_perp_avg[2], kappa_cross_avg[2] = {0.0};
  // Parallel current multiplied by parallel thermal force coefficient at cell edges (using arithmetic average)
  double current_par_avg[3] = {0.0};
  // Cross current multiplied by perpendicular thermal force coefficient at cell edges (using arithmetic average)
  double current_cross_avg[3] = {0.0};
  // Temperature gradient, parallel, perp, and cross at cell edges
  double gradxT[2] = {0.0};
  double gradyT[2] = {0.0};
  double gradzT[2] = {0.0};
  // Parallel temperature gradient (b_hat b_hat dot grad T)
  double bbgradT[2][3] = {0.0};
  // Perpendicular temperature gradient (grad T - b_hat b_hat dot grad T)
  double perp_gradT[2][3] = {0.0};
  // Cross temperature gradient (b x grad T)
  double cross_gradT[2][3] = {0.0};
  // Average velocity at cell edges
  double u_avg[2][3] = {0.0};

  if (ndim == 1) {
    double dx = bes->grid.dx[0]; 
    // Derived quantities
    double u[2][2][3] = {0.0}; // Flow for each species, ux, uy, & uz
    double b[2][3] = {0.0}; // Magnetic field unit vector, Bx/|B|, By/|B|, & Bz/|B|
    double T[2][2] = {0.0}; // Temperature for each species (m*p/rho)
    double tau[2][2] = {0.0}; // Collision times for each species
    double eta_par[2][2] = {0.0}; // Parallel viscosity for each species
    double eta_perp[2][2] = {0.0}; // Perpendicular viscosity for each species
    double eta_cross[2][2] = {0.0}; // Gyro-viscosity for each species
    double kappa_par[2][2] = {0.0}; // Parallel conductivity for each species
    double kappa_perp[2][2] = {0.0}; // Perpendicular conductivity for each species
    double kappa_cross[2][2] = {0.0}; // Gyro-conductivity for each species
    double current_par[2][3] = {0.0}; // b_hat b_hat dot (u_i - u_e) 
    double current_cross[2][3] = {0.0}; // b x (u_i - u_e)
    // Compute derived quantities
    mag_var_setup(bes, L_1D, U_1D,
      fluid_d, em_tot_d,
      u, b, T, tau,
      eta_par, eta_perp, eta_cross,
      kappa_par, kappa_perp, kappa_cross,
      current_par, current_cross);
    // Magnetic field at cell edges (using arithmetic average)
    b_avg[0] = calc_arithm_avg_1D(b[L_1D][0], b[U_1D][0]);
    b_avg[1] = calc_arithm_avg_1D(b[L_1D][1], b[U_1D][1]);
    b_avg[2] = calc_arithm_avg_1D(b[L_1D][2], b[U_1D][2]);
    
    // Parallel viscosity coefficients at cell edges (using harmonic average)
    eta_par_avg[ELC] = calc_harmonic_avg_1D(eta_par[L_1D][ELC], eta_par[U_1D][ELC]);
    eta_par_avg[ION] = calc_harmonic_avg_1D(eta_par[L_1D][ION], eta_par[U_1D][ION]);
    // Perpendicular viscosity coefficients at cell edges (using harmonic average)
    eta_perp_avg[ELC] = calc_harmonic_avg_1D(eta_perp[L_1D][ELC], eta_perp[U_1D][ELC]);
    eta_perp_avg[ION] = calc_harmonic_avg_1D(eta_perp[L_1D][ION], eta_perp[U_1D][ION]);   
    // Gyro-viscosity coefficients at cell edges (using harmonic average)
    eta_cross_avg[ELC] = calc_harmonic_avg_1D(eta_cross[L_1D][ELC], eta_cross[U_1D][ELC]);
    eta_cross_avg[ION] = calc_harmonic_avg_1D(eta_cross[L_1D][ION], eta_cross[U_1D][ION]);
      
    // Rate of strain tensor at cell edges for electrons and ions
    calc_ros_1D(dx, u[L_1D][ELC], u[U_1D][ELC], w[ELC]);
    calc_ros_1D(dx, u[L_1D][ION], u[U_1D][ION], w[ION]);
    
    // Calculate heat flux and viscous heating if energy variable exists
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER) {
      // Parallel conductivity coefficients at cell edges (using harmonic average)
      kappa_par_avg[ELC] = calc_harmonic_avg_1D(kappa_par[L_1D][ELC], kappa_par[U_1D][ELC]);
      // Perpendicular conductivity coefficients at cell edges (using harmonic average)
      kappa_perp_avg[ELC] = calc_harmonic_avg_1D(kappa_perp[L_1D][ELC], kappa_perp[U_1D][ELC]);
      // Gyro-conductivity coefficients at cell edges (using harmonic average)
      kappa_cross_avg[ELC] = calc_harmonic_avg_1D(kappa_cross[L_1D][ELC], kappa_cross[U_1D][ELC]);
        
      // Parallel current multiplied by parallel thermal force coefficient at cell edges (using arithmetic average)
      current_par_avg[0] = calc_arithm_avg_1D(current_par[L_1D][0], current_par[U_1D][0]);
      current_par_avg[1] = calc_arithm_avg_1D(current_par[L_1D][1], current_par[U_1D][1]);
      current_par_avg[2] = calc_arithm_avg_1D(current_par[L_1D][2], current_par[U_1D][2]);
      // Cross current multiplied by perpendicular thermal force coefficient at cell edges (using arithmetic average)
      current_cross_avg[0] = calc_arithm_avg_1D(current_cross[L_1D][0], current_cross[U_1D][0]);
      current_cross_avg[1] = calc_arithm_avg_1D(current_cross[L_1D][1], current_cross[U_1D][1]);
      current_cross_avg[2] = calc_arithm_avg_1D(current_cross[L_1D][2], current_cross[U_1D][2]);
    
      gradxT[ELC] = calc_sym_grad_1D(dx, T[L_1D][ELC], T[U_1D][ELC]);
    
      bbgradT[ELC][0] = b_avg[0]*b_avg[0]*gradxT[ELC];
      bbgradT[ELC][1] = b_avg[1]*b_avg[0]*gradxT[ELC];
      bbgradT[ELC][2] = b_avg[2]*b_avg[0]*gradxT[ELC];
    
      perp_gradT[ELC][0] = gradxT[ELC] - bbgradT[ELC][0];
    
      cross_gradT[ELC][1] = b_avg[2]*gradxT[ELC];
      cross_gradT[ELC][2] = -b_avg[1]*gradxT[ELC];

      // Average velocity at cell edges (using arithmetic average)
      u_avg[ELC][0] = calc_arithm_avg_1D(u[L_1D][ELC][0], u[U_1D][ELC][0]);
      u_avg[ELC][1] = calc_arithm_avg_1D(u[L_1D][ELC][1], u[U_1D][ELC][1]);
      u_avg[ELC][2] = calc_arithm_avg_1D(u[L_1D][ELC][2], u[U_1D][ELC][2]);
    }
    if (bes->param[ION].type_eqn == GKYL_EQN_EULER) {
      // Parallel conductivity coefficients at cell edges (using harmonic average)
      kappa_par_avg[ION] = calc_harmonic_avg_1D(kappa_par[L_1D][ION], kappa_par[U_1D][ION]);
      // Perpendicular conductivity coefficients at cell edges (using harmonic average)
      kappa_perp_avg[ION] = calc_harmonic_avg_1D(kappa_perp[L_1D][ION], kappa_perp[U_1D][ION]);
      // Gyro-conductivity coefficients at cell edges (using harmonic average)
      kappa_cross_avg[ION] = calc_harmonic_avg_1D(kappa_cross[L_1D][ION], kappa_cross[U_1D][ION]);
    
      gradxT[ION] = calc_sym_grad_1D(dx, T[L_1D][ION], T[U_1D][ION]);
    
      bbgradT[ION][0] = b_avg[0]*b_avg[0]*gradxT[ION];
      bbgradT[ION][1] = b_avg[1]*b_avg[0]*gradxT[ION];
      bbgradT[ION][2] = b_avg[2]*b_avg[0]*gradxT[ION];
    
      perp_gradT[ION][0] = gradxT[ION] - bbgradT[ION][0];
    
      cross_gradT[ION][1] = b_avg[2]*gradxT[ION];
      cross_gradT[ION][2] = -b_avg[1]*gradxT[ION];

      // Average velocity at cell edges (using arithmetic average)
      u_avg[ION][0] = calc_arithm_avg_1D(u[L_1D][ION][0], u[U_1D][ION][0]);
      u_avg[ION][1] = calc_arithm_avg_1D(u[L_1D][ION][1], u[U_1D][ION][1]);
      u_avg[ION][2] = calc_arithm_avg_1D(u[L_1D][ION][2], u[U_1D][ION][2]);
    }
  }
  else if (ndim == 2) {
    double dx = bes->grid.dx[0];
    double dy = bes->grid.dx[1];
    // Derived quantities
    double u[4][2][3] = {0.0}; // Flow for each species, ux, uy, & uz
    double b[4][3] = {0.0}; // Magnetic field unit vector, Bx/|B|, By/|B|, & Bz/|B|
    double T[4][2] = {0.0}; // Temperature for each species (m*p/rho)
    double tau[4][2] = {0.0}; // Collision times for each species
    double eta_par[4][2] = {0.0}; // Parallel viscosity for each species
    double eta_perp[4][2] = {0.0}; // Perpendicular viscosity for each species
    double eta_cross[4][2] = {0.0}; // Gyro-viscosity for each species
    double kappa_par[4][2] = {0.0}; // Parallel conductivity for each species
    double kappa_perp[4][2] = {0.0}; // Perpendicular conductivity for each species
    double kappa_cross[4][2] = {0.0}; // Gyro-conductivity for each species
    double current_par[4][3] = {0.0}; // b_hat b_hat dot (u_i - u_e) 
    double current_cross[4][3] = {0.0}; // b x (u_i - u_e)
    // Compute derived quantities
    mag_var_setup(bes, LL_2D, UU_2D,
      fluid_d, em_tot_d,
      u, b, T, tau,
      eta_par, eta_perp, eta_cross,
      kappa_par, kappa_perp, kappa_cross,
      current_par, current_cross);
    // Magnetic field at cell vertices (using arithmetic average)
    b_avg[0] = calc_arithm_avg_2D(b[LL_2D][0], b[LU_2D][0], b[UL_2D][0], b[UU_2D][0]);
    b_avg[1] = calc_arithm_avg_2D(b[LL_2D][1], b[LU_2D][1], b[UL_2D][1], b[UU_2D][1]);
    b_avg[2] = calc_arithm_avg_2D(b[LL_2D][2], b[LU_2D][2], b[UL_2D][2], b[UU_2D][2]);

    // Parallel viscosity coefficients at cell edges (using harmonic average)
    eta_par_avg[ELC] = calc_harmonic_avg_2D(eta_par[LL_2D][ELC], eta_par[LU_2D][ELC], eta_par[UL_2D][ELC], eta_par[UU_2D][ELC]);
    eta_par_avg[ION] = calc_harmonic_avg_2D(eta_par[LL_2D][ION], eta_par[LU_2D][ION], eta_par[UL_2D][ION], eta_par[UU_2D][ION]);

    // Perpendicular viscosity coefficients at cell edges (using harmonic average)
    eta_perp_avg[ELC] = calc_harmonic_avg_2D(eta_perp[LL_2D][ELC], eta_perp[LU_2D][ELC], eta_perp[UL_2D][ELC], eta_perp[UU_2D][ELC]);
    eta_perp_avg[ION] = calc_harmonic_avg_2D(eta_perp[LL_2D][ION], eta_perp[LU_2D][ION], eta_perp[UL_2D][ION], eta_perp[UU_2D][ION]);

    // Gyro-viscosity coefficients at cell edges (using harmonic average)
    eta_cross_avg[ELC] = calc_harmonic_avg_2D(eta_cross[LL_2D][ELC], eta_cross[LU_2D][ELC], eta_cross[UL_2D][ELC], eta_cross[UU_2D][ELC]);
    eta_cross_avg[ION] = calc_harmonic_avg_2D(eta_cross[LL_2D][ION], eta_cross[LU_2D][ION], eta_cross[UL_2D][ION], eta_cross[UU_2D][ION]);

    // Rate of strain tensor at cell vertices for electrons and ions
    calc_ros_2D(dx, dy, u[LL_2D][ELC], u[LU_2D][ELC], u[UL_2D][ELC], u[UU_2D][ELC], w[ELC]);
    calc_ros_2D(dx, dy, u[LL_2D][ION], u[LU_2D][ION], u[UL_2D][ION], u[UU_2D][ION], w[ION]);
    
    // Calculate heat flux and viscous heating if energy variable exists
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER) {
      // Parallel conductivity coefficients at cell edges (using harmonic average)
      kappa_par_avg[ELC] = calc_harmonic_avg_2D(kappa_par[LL_2D][ELC], kappa_par[LU_2D][ELC], kappa_par[UL_2D][ELC], kappa_par[UU_2D][ELC]);
      // Perpendicular conductivity coefficients at cell edges (using harmonic average)
      kappa_perp_avg[ELC] = calc_harmonic_avg_2D(kappa_perp[LL_2D][ELC], kappa_perp[LU_2D][ELC], kappa_perp[UL_2D][ELC], kappa_perp[UU_2D][ELC]);
      // Gyro-conductivity coefficients at cell edges (using harmonic average)
      kappa_cross_avg[ELC] = calc_harmonic_avg_2D(kappa_cross[LL_2D][ELC], kappa_cross[LU_2D][ELC], kappa_cross[UL_2D][ELC], kappa_cross[UU_2D][ELC]);

      // Parallel current multiplied by parallel thermal force coefficient at cell edges (using arithmetic average)
      current_par_avg[0] = calc_arithm_avg_2D(current_par[LL_2D][0], current_par[LU_2D][0], current_par[UL_2D][0], current_par[UU_2D][0]);
      current_par_avg[1] = calc_arithm_avg_2D(current_par[LL_2D][1], current_par[LU_2D][1], current_par[UL_2D][1], current_par[UU_2D][1]);
      current_par_avg[2] = calc_arithm_avg_2D(current_par[LL_2D][2], current_par[LU_2D][2], current_par[UL_2D][2], current_par[UU_2D][2]);
      
      // Cross current multiplied by perpendicular thermal force coefficient at cell edges (using arithmetic average)
      current_cross_avg[0] = calc_arithm_avg_2D(current_cross[LL_2D][0], current_cross[LU_2D][0], current_cross[UL_2D][0], current_cross[UU_2D][0]);
      current_cross_avg[1] = calc_arithm_avg_2D(current_cross[LL_2D][1], current_cross[LU_2D][1], current_cross[UL_2D][1], current_cross[UU_2D][1]);
      current_cross_avg[2] = calc_arithm_avg_2D(current_cross[LL_2D][2], current_cross[LU_2D][2], current_cross[UL_2D][2], current_cross[UU_2D][2]);
    
      gradxT[ELC] = calc_sym_gradx_2D(dx, T[LL_2D][ELC], T[LU_2D][ELC], T[UL_2D][ELC], T[UU_2D][ELC]);
      gradyT[ELC] = calc_sym_grady_2D(dy, T[LL_2D][ELC], T[LU_2D][ELC], T[UL_2D][ELC], T[UU_2D][ELC]);
    
      bbgradT[ELC][0] = b_avg[0]*(b_avg[0]*gradxT[ELC] + b_avg[1]*gradyT[ELC]);
      bbgradT[ELC][1] = b_avg[1]*(b_avg[0]*gradxT[ELC] + b_avg[1]*gradyT[ELC]);
      bbgradT[ELC][2] = b_avg[2]*(b_avg[0]*gradxT[ELC] + b_avg[1]*gradyT[ELC]);
    
      perp_gradT[ELC][0] = gradxT[ELC] - bbgradT[ELC][0];
      perp_gradT[ELC][1] = gradyT[ELC] - bbgradT[ELC][1];
    
      cross_gradT[ELC][0] = -b_avg[2]*gradyT[ELC];
      cross_gradT[ELC][1] = b_avg[2]*gradxT[ELC];
      cross_gradT[ELC][2] = b_avg[0]*gradyT[ELC] - b_avg[1]*gradxT[ELC];

      // Average velocity at cell edges (using arithmetic average)
      u_avg[ELC][0] = calc_arithm_avg_2D(u[LL_2D][ELC][0], u[LU_2D][ELC][0], u[UL_2D][ELC][0], u[UU_2D][ELC][0]);
      u_avg[ELC][1] = calc_arithm_avg_2D(u[LL_2D][ELC][1], u[LU_2D][ELC][1], u[UL_2D][ELC][1], u[UU_2D][ELC][1]);
      u_avg[ELC][2] = calc_arithm_avg_2D(u[LL_2D][ELC][2], u[LU_2D][ELC][2], u[UL_2D][ELC][2], u[UU_2D][ELC][2]);
    }
    if (bes->param[ION].type_eqn == GKYL_EQN_EULER) {
      // Parallel conductivity coefficients at cell edges (using harmonic average)
      kappa_par_avg[ION] = calc_harmonic_avg_2D(kappa_par[LL_2D][ION], kappa_par[LU_2D][ION], kappa_par[UL_2D][ION], kappa_par[UU_2D][ION]);
      // Perpendicular conductivity coefficients at cell edges (using harmonic average)
      kappa_perp_avg[ION] = calc_harmonic_avg_2D(kappa_perp[LL_2D][ION], kappa_perp[LU_2D][ION], kappa_perp[UL_2D][ION], kappa_perp[UU_2D][ION]);
      // Gyro-conductivity coefficients at cell edges (using harmonic average)
      kappa_cross_avg[ION] = calc_harmonic_avg_2D(kappa_cross[LL_2D][ION], kappa_cross[LU_2D][ION], kappa_cross[UL_2D][ION], kappa_cross[UU_2D][ION]);
    
      gradxT[ION] = calc_sym_gradx_2D(dx, T[LL_2D][ION], T[LU_2D][ION], T[UL_2D][ION], T[UU_2D][ION]);
      gradyT[ION] = calc_sym_grady_2D(dy, T[LL_2D][ION], T[LU_2D][ION], T[UL_2D][ION], T[UU_2D][ION]);
    
      bbgradT[ION][0] = b_avg[0]*(b_avg[0]*gradxT[ION] + b_avg[1]*gradyT[ION]);
      bbgradT[ION][1] = b_avg[1]*(b_avg[0]*gradxT[ION] + b_avg[1]*gradyT[ION]);
      bbgradT[ION][2] = b_avg[2]*(b_avg[0]*gradxT[ION] + b_avg[1]*gradyT[ION]);
    
      perp_gradT[ION][0] = gradxT[ION] - bbgradT[ION][0];
      perp_gradT[ION][1] = gradyT[ION] - bbgradT[ION][1];
    
      cross_gradT[ION][0] = -b_avg[2]*gradyT[ION];
      cross_gradT[ION][1] = b_avg[2]*gradxT[ION];
      cross_gradT[ION][2] = b_avg[0]*gradyT[ION] - b_avg[1]*gradxT[ION];

      // Velocity at cell edges (using arithmetic average)
      u_avg[ION][0] = calc_arithm_avg_2D(u[LL_2D][ION][0], u[LU_2D][ION][0], u[UL_2D][ION][0], u[UU_2D][ION][0]);
      u_avg[ION][1] = calc_arithm_avg_2D(u[LL_2D][ION][1], u[LU_2D][ION][1], u[UL_2D][ION][1], u[UU_2D][ION][1]);
      u_avg[ION][2] = calc_arithm_avg_2D(u[LL_2D][ION][2], u[LU_2D][ION][2], u[UL_2D][ION][2], u[UU_2D][ION][2]);
    }
  }
  // Parallel viscous stress tensor at cell edges for electrons and ions
  calc_pi_par(eta_par_avg[ELC], b_avg, w[ELC], pi_par[ELC]);
  calc_pi_par(eta_par_avg[ION], b_avg, w[ION], pi_par[ION]);

  // Perpendicular viscous stress tensor at cell edges for electrons and ions
  calc_pi_perp(eta_perp_avg[ELC], b_avg, w[ELC], pi_perp[ELC]);
  calc_pi_perp(eta_perp_avg[ION], b_avg, w[ION], pi_perp[ION]);

  // Gyro-viscous stress tensor at cell edges for electrons and ions
  calc_pi_cross(eta_cross_avg[ELC], b_avg, w[ELC], pi_cross[ELC]);
  calc_pi_cross(eta_cross_avg[ION], b_avg, w[ION], pi_cross[ION]);
  
  // Total viscous stress tensor
  brag_d[ELC][PIXX] = pi_par[ELC][0]+pi_perp[ELC][0]+pi_cross[ELC][0];
  brag_d[ELC][PIXY] = pi_par[ELC][1]+pi_perp[ELC][1]+pi_cross[ELC][1];
  brag_d[ELC][PIXZ] = pi_par[ELC][2]+pi_perp[ELC][2]+pi_cross[ELC][2];
  brag_d[ELC][PIYY] = pi_par[ELC][3]+pi_perp[ELC][3]+pi_cross[ELC][3];
  brag_d[ELC][PIYZ] = pi_par[ELC][4]+pi_perp[ELC][4]+pi_cross[ELC][4];
  brag_d[ELC][PIZZ] = pi_par[ELC][5]+pi_perp[ELC][5]+pi_cross[ELC][5];
  // Total heat flux + viscous heating 
  brag_d[ELC][QX] = -kappa_par_avg[ELC]*bbgradT[ELC][0] + current_par_avg[0]
                    -kappa_perp_avg[ELC]*perp_gradT[ELC][0] + current_cross_avg[0]
                    +kappa_cross_avg[ELC]*cross_gradT[ELC][0]
                    +brag_d[ELC][PIXX]*u_avg[ELC][0] + brag_d[ELC][PIXY]*u_avg[ELC][1] + brag_d[ELC][PIXZ]*u_avg[ELC][2];
  
  brag_d[ELC][QY] = -kappa_par_avg[ELC]*bbgradT[ELC][1] + current_par_avg[1]
                    -kappa_perp_avg[ELC]*perp_gradT[ELC][1] + current_cross_avg[1]
                    +kappa_cross_avg[ELC]*cross_gradT[ELC][1]
                    +brag_d[ELC][PIXY]*u_avg[ELC][0] + brag_d[ELC][PIYY]*u_avg[ELC][1] + brag_d[ELC][PIYZ]*u_avg[ELC][2];
  
  brag_d[ELC][QZ] = -kappa_par_avg[ELC]*bbgradT[ELC][2] + current_par_avg[2]
                    -kappa_perp_avg[ELC]*perp_gradT[ELC][2] + current_cross_avg[2]
                    +kappa_cross_avg[ELC]*cross_gradT[ELC][2]
                    +brag_d[ELC][PIXZ]*u_avg[ELC][0] + brag_d[ELC][PIYZ]*u_avg[ELC][1] + brag_d[ELC][PIZZ]*u_avg[ELC][2];
  
  brag_d[ION][PIXX] = pi_par[ION][0]+pi_perp[ION][0]+pi_cross[ION][0];
  brag_d[ION][PIXY] = pi_par[ION][1]+pi_perp[ION][1]+pi_cross[ION][1];
  brag_d[ION][PIXZ] = pi_par[ION][2]+pi_perp[ION][2]+pi_cross[ION][2];
  brag_d[ION][PIYY] = pi_par[ION][3]+pi_perp[ION][3]+pi_cross[ION][3];
  brag_d[ION][PIYZ] = pi_par[ION][4]+pi_perp[ION][4]+pi_cross[ION][4];
  brag_d[ION][PIZZ] = pi_par[ION][5]+pi_perp[ION][5]+pi_cross[ION][5];
  // Total heat flux + viscous heating 
  brag_d[ION][QX] = -kappa_par_avg[ION]*bbgradT[ION][0] - kappa_perp_avg[ION]*perp_gradT[ION][0] + kappa_cross_avg[ION]*cross_gradT[ION][0] +
                    brag_d[ION][PIXX]*u_avg[ION][0] + brag_d[ION][PIXY]*u_avg[ION][1] + brag_d[ION][PIXZ]*u_avg[ION][2];
  
  brag_d[ION][QY] = -kappa_par_avg[ION]*bbgradT[ION][1] - kappa_perp_avg[ION]*perp_gradT[ION][1] + kappa_cross_avg[ION]*cross_gradT[ION][1] +
                    brag_d[ION][PIXY]*u_avg[ION][0] + brag_d[ION][PIYY]*u_avg[ION][1] + brag_d[ION][PIYZ]*u_avg[ION][2];
  
  brag_d[ION][QZ] = -kappa_par_avg[ION]*bbgradT[ION][2] - kappa_perp_avg[ION]*perp_gradT[ION][2] + kappa_cross_avg[ION]*cross_gradT[ION][2] +
                    brag_d[ION][PIXZ]*u_avg[ION][0] + brag_d[ION][PIYZ]*u_avg[ION][1] + brag_d[ION][PIZZ]*u_avg[ION][2];
}

// Fetch input quantities and compute derived quantities for UNmagnetized Braginskii
static void
unmag_var_setup(const gkyl_moment_braginskii *bes, int start, int end,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double u[][2][3],
  double T[][2], double tau[][2],
  double eta[][2], double kappa[][2],
  double current[][3])
{
  int nfluids = bes->nfluids;

  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };
  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double rho[2] = {0.0}; // Mass density for each species
  double p[2] = {0.0}; // Pressure for each species
  for (int j = start; j <= end; ++j) {
    // Input density and flow in each cell
    rho[ELC] = fluid_d[j][ELC][RHO];
    rho[ION] = fluid_d[j][ION][RHO];
    for (int k = MX; k <= MZ; ++k) {
      u[j][ELC][k - MX] = fluid_d[j][ELC][k] / rho[ELC];
      u[j][ION][k - MX] = fluid_d[j][ION][k] / rho[ION];
    }
    // Pressure information is different for each equation type
    for (int n = 0; n < nfluids; ++n) {
      if (bes->param[n].type_eqn == GKYL_EQN_EULER)
        p[n] = gkyl_euler_pressure(bes->param[n].p_fac, fluid_d[j][n]); // Euler needs to divide out gas_gamma factor to obtain pressure
      else if (bes->param[n].type_eqn == GKYL_EQN_ISO_EULER)
        p[n] = rho[n]*bes->param[n].p_fac*bes->param[n].p_fac; // isothermal Euler input is vth, pressure = rho*vth^2
    }
    T[j][ELC] = m[ELC] * p[ELC] / rho[ELC];
    T[j][ION] = m[ION] * p[ION] / rho[ION];

    tau[j][ELC] = calc_tau(1.0, bes->coll_fac, bes->epsilon0, q[ELC], q[ION], m[ELC], m[ION], rho[ION], T[j][ELC]);
    tau[j][ION] = calc_tau(1.0, sqrt(2.0)*bes->coll_fac, bes->epsilon0, q[ION], q[ION], m[ION], m[ION], rho[ION], T[j][ION]);

    eta[j][ELC] = 0.73*p[ELC]*tau[j][ELC];
    kappa[j][ELC] = 3.16*p[ELC]*tau[j][ELC]/m[ELC];

    eta[j][ION] = 0.96*p[ION]*tau[j][ION];
    kappa[j][ION] = 3.91*p[ION]*tau[j][ION]/m[ION];

    double thermal = -0.71*p[ELC]; // thermal force coefficient (same for each species)      
    current[j][0] = thermal*(u[j][ION][0] - u[j][ELC][0]);
    current[j][1] = thermal*(u[j][ION][1] - u[j][ELC][1]);
    current[j][2] = thermal*(u[j][ION][2] - u[j][ELC][2]);
  }
}
    
static void
unmag_brag_calc_vars(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double *cflrate, double *brag_d[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  const int ndim = bes->ndim;
  // Allocate some memory on stack regardless of dimensionality or equation type (Euler vs. isothermal Euler)
  // These allocations allow us to make final construction of Braginskii variables dimensionally independent
  // Note: Braginskii implementation currently assumes only two fluids (electrons and an ion species)
  double m[2] = { bes->param[0].mass, bes->param[1].mass };
  double q[2] = { bes->param[0].charge, bes->param[1].charge };
  // Grab indices of electron and ion fluid arrays
  int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
  int ION = (ELC + 1) % 2;

  double eta_avg[2] = {0.0};
  double w[2][6] = {0.0};

  double kappa_avg[2] = {0.0};
  // Current multiplied by thermal force coefficient at cell edges (using arithmetic average)
  double current_avg[3] = {0.0};
  // Temperature gradient at cell edges
  double gradxT[2] = {0.0};
  double gradyT[2] = {0.0};
  double gradzT[2] = {0.0};
  // Average velocity at cell edges
  double u_avg[2][3] = {0.0};

  if (ndim == 1) {
    double dx = bes->grid.dx[0]; 
    // Derived quantities
    double u[2][2][3] = {0.0}; // Flow for each species, ux, uy, & uz
    double T[2][2] = {0.0}; // Temperature for each species (m*p/rho)
    double tau[2][2] = {0.0}; // Collision times for each species
    double eta[2][2] = {0.0}; // Viscosity for each species
    double kappa[2][2] = {0.0}; // Conductivity for each species
    double current[2][3] = {0.0}; // (u_i - u_e) 
    // Compute derived quantities
    unmag_var_setup(bes, L_1D, U_1D, fluid_d,
      u, T, tau,
      eta, kappa, current);
    // Viscosity coefficients at cell edges (using harmonic average)
    eta_avg[ELC] = calc_harmonic_avg_1D(eta[L_1D][ELC], eta[U_1D][ELC]);
    eta_avg[ION] = calc_harmonic_avg_1D(eta[L_1D][ION], eta[U_1D][ION]);

    // Rate of strain tensor at cell edges for electrons and ions
    double w[2][6] = {0.0};
    calc_ros_1D(dx, u[L_1D][ELC], u[U_1D][ELC], w[ELC]);
    calc_ros_1D(dx, u[L_1D][ION], u[U_1D][ION], w[ION]);

    // Calculate heat flux and viscous heating if energy variable exists
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER) {
      // Conductivity coefficients at cell edges (using harmonic average)
      kappa_avg[ELC] = calc_harmonic_avg_1D(kappa[L_1D][ELC], kappa[U_1D][ELC]);
        
      // Current multiplied by thermal force coefficient at cell edges (using arithmetic average)
      current_avg[0] = calc_arithm_avg_1D(current[L_1D][0], current[U_1D][0]);
      current_avg[1] = calc_arithm_avg_1D(current[L_1D][1], current[U_1D][1]);
      current_avg[2] = calc_arithm_avg_1D(current[L_1D][2], current[U_1D][2]);
    
      gradxT[ELC] = calc_sym_grad_1D(dx, T[L_1D][ELC], T[U_1D][ELC]);

      // Average velocity at cell edges (using arithmetic average)
      u_avg[ELC][0] = calc_arithm_avg_1D(u[L_1D][ELC][0], u[U_1D][ELC][0]);
      u_avg[ELC][1] = calc_arithm_avg_1D(u[L_1D][ELC][1], u[U_1D][ELC][1]);
      u_avg[ELC][2] = calc_arithm_avg_1D(u[L_1D][ELC][2], u[U_1D][ELC][2]);
    }
    if (bes->param[ION].type_eqn == GKYL_EQN_EULER) {
      // Conductivity coefficients at cell edges (using harmonic average)
      kappa_avg[ION] = calc_harmonic_avg_1D(kappa[L_1D][ION], kappa[U_1D][ION]);
    
      gradxT[ION] = calc_sym_grad_1D(dx, T[L_1D][ION], T[U_1D][ION]);
      
      // Average velocity at cell edges (using arithmetic average)
      u_avg[ION][0] = calc_arithm_avg_1D(u[L_1D][ION][0], u[U_1D][ION][0]);
      u_avg[ION][1] = calc_arithm_avg_1D(u[L_1D][ION][1], u[U_1D][ION][1]);
      u_avg[ION][2] = calc_arithm_avg_1D(u[L_1D][ION][2], u[U_1D][ION][2]);
    }
  }
  else if (ndim == 2) {
    double dx = bes->grid.dx[0];
    double dy = bes->grid.dx[1]; 
    // Derived quantities
    double u[4][2][3] = {0.0}; // Flow for each species, ux, uy, & uz
    double T[4][2] = {0.0}; // Temperature for each species (m*p/rho)
    double tau[4][2] = {0.0}; // Collision times for each species
    double eta[4][2] = {0.0}; // Viscosity for each species
    double kappa[4][2] = {0.0}; // Conductivity for each species
    double current[4][3] = {0.0}; // (u_i - u_e) 
    // Compute derived quantities
    unmag_var_setup(bes, LL_2D, UU_2D, fluid_d,
      u, T, tau,
      eta, kappa, current);

    // Viscosity coefficients at cell edges (using harmonic average)
    eta_avg[ELC] = calc_harmonic_avg_2D(eta[LL_2D][ELC], eta[LU_2D][ELC], eta[UL_2D][ELC], eta[UU_2D][ELC]);
    eta_avg[ION] = calc_harmonic_avg_2D(eta[LL_2D][ION], eta[LU_2D][ION], eta[UL_2D][ION], eta[UU_2D][ION]);

    // Rate of strain tensor at cell vertices for electrons and ions
    calc_ros_2D(dx, dy, u[LL_2D][ELC], u[LU_2D][ELC], u[UL_2D][ELC], u[UU_2D][ELC], w[ELC]);
    calc_ros_2D(dx, dy, u[LL_2D][ION], u[LU_2D][ION], u[UL_2D][ION], u[UU_2D][ION], w[ION]);
    
    // Calculate heat flux and viscous heating if energy variable exists
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER) {
      // Conductivity coefficients at cell edges (using harmonic average)
      kappa_avg[ELC] = calc_harmonic_avg_2D(kappa[LL_2D][ELC], kappa[LU_2D][ELC], kappa[UL_2D][ELC], kappa[UU_2D][ELC]);

      // Current multiplied by thermal force coefficient at cell edges (using arithmetic average)
      current_avg[0] = calc_arithm_avg_2D(current[LL_2D][0], current[LU_2D][0], current[UL_2D][0], current[UU_2D][0]);
      current_avg[1] = calc_arithm_avg_2D(current[LL_2D][1], current[LU_2D][1], current[UL_2D][1], current[UU_2D][1]);
      current_avg[2] = calc_arithm_avg_2D(current[LL_2D][2], current[LU_2D][2], current[UL_2D][2], current[UU_2D][2]);
    
      gradxT[ELC] = calc_sym_gradx_2D(dx, T[LL_2D][ELC], T[LU_2D][ELC], T[UL_2D][ELC], T[UU_2D][ELC]);
      gradyT[ELC] = calc_sym_grady_2D(dy, T[LL_2D][ELC], T[LU_2D][ELC], T[UL_2D][ELC], T[UU_2D][ELC]);

      // Average velocity at cell edges (using arithmetic average)
      u_avg[ELC][0] = calc_arithm_avg_2D(u[LL_2D][ELC][0], u[LU_2D][ELC][0], u[UL_2D][ELC][0], u[UU_2D][ELC][0]);
      u_avg[ELC][1] = calc_arithm_avg_2D(u[LL_2D][ELC][1], u[LU_2D][ELC][1], u[UL_2D][ELC][1], u[UU_2D][ELC][1]);
      u_avg[ELC][2] = calc_arithm_avg_2D(u[LL_2D][ELC][2], u[LU_2D][ELC][2], u[UL_2D][ELC][2], u[UU_2D][ELC][2]);
    }
    if (bes->param[ION].type_eqn == GKYL_EQN_EULER) {
      // Conductivity coefficients at cell edges (using harmonic average)
      kappa_avg[ION] = calc_harmonic_avg_2D(kappa[LL_2D][ION], kappa[LU_2D][ION], kappa[UL_2D][ION], kappa[UU_2D][ION]);
    
      gradxT[ION] = calc_sym_gradx_2D(dx, T[LL_2D][ION], T[LU_2D][ION], T[UL_2D][ION], T[UU_2D][ION]);
      gradyT[ION] = calc_sym_grady_2D(dy, T[LL_2D][ION], T[LU_2D][ION], T[UL_2D][ION], T[UU_2D][ION]);

      // Velocity at cell edges (using arithmetic average)
      u_avg[ION][0] = calc_arithm_avg_2D(u[LL_2D][ION][0], u[LU_2D][ION][0], u[UL_2D][ION][0], u[UU_2D][ION][0]);
      u_avg[ION][1] = calc_arithm_avg_2D(u[LL_2D][ION][1], u[LU_2D][ION][1], u[UL_2D][ION][1], u[UU_2D][ION][1]);
      u_avg[ION][2] = calc_arithm_avg_2D(u[LL_2D][ION][2], u[LU_2D][ION][2], u[UL_2D][ION][2], u[UU_2D][ION][2]);
    }
  }
  // Total viscous stress tensor for electrons
  brag_d[ELC][PIXX] = -eta_avg[ELC]*w[ELC][0];
  brag_d[ELC][PIXY] = -eta_avg[ELC]*w[ELC][1];
  brag_d[ELC][PIXZ] = -eta_avg[ELC]*w[ELC][2];
  brag_d[ELC][PIYY] = -eta_avg[ELC]*w[ELC][3];
  brag_d[ELC][PIYZ] = -eta_avg[ELC]*w[ELC][4];
  brag_d[ELC][PIZZ] = -eta_avg[ELC]*w[ELC][5];

  // Total heat flux + viscous heating for electrons
  brag_d[ELC][QX] = -kappa_avg[ELC]*gradxT[ELC] + current_avg[0] +
                    brag_d[ELC][PIXX]*u_avg[ELC][0] + brag_d[ELC][PIXY]*u_avg[ELC][1] + brag_d[ELC][PIXZ]*u_avg[ELC][2];
  brag_d[ELC][QY] = -kappa_avg[ELC]*gradyT[ELC] + current_avg[1] +
                    brag_d[ELC][PIXY]*u_avg[ELC][0] + brag_d[ELC][PIYY]*u_avg[ELC][1] + brag_d[ELC][PIYZ]*u_avg[ELC][2];
  brag_d[ELC][QZ] = -kappa_avg[ELC]*gradzT[ELC] + current_avg[2] +
                    brag_d[ELC][PIXZ]*u_avg[ELC][0] + brag_d[ELC][PIYZ]*u_avg[ELC][1] + brag_d[ELC][PIZZ]*u_avg[ELC][2];

  // Total viscous stress tensor for ions
  brag_d[ION][PIXX] = -eta_avg[ION]*w[ION][0];
  brag_d[ION][PIXY] = -eta_avg[ION]*w[ION][1];
  brag_d[ION][PIXZ] = -eta_avg[ION]*w[ION][2];
  brag_d[ION][PIYY] = -eta_avg[ION]*w[ION][3];
  brag_d[ION][PIYZ] = -eta_avg[ION]*w[ION][4];
  brag_d[ION][PIZZ] = -eta_avg[ION]*w[ION][5];

  // Total heat flux + viscous heating for ions
  brag_d[ION][QX] = -kappa_avg[ION]*gradxT[ION] +
                    brag_d[ION][PIXX]*u_avg[ION][0] + brag_d[ION][PIXY]*u_avg[ION][1] + brag_d[ION][PIXZ]*u_avg[ION][2];
  brag_d[ION][QY] = -kappa_avg[ION]*gradxT[ION] +
                    brag_d[ION][PIXY]*u_avg[ION][0] + brag_d[ION][PIYY]*u_avg[ION][1] + brag_d[ION][PIYZ]*u_avg[ION][2];
  brag_d[ION][QZ] = -kappa_avg[ION]*gradxT[ION] +
                    brag_d[ION][PIXZ]*u_avg[ION][0] + brag_d[ION][PIYZ]*u_avg[ION][1] + brag_d[ION][PIZZ]*u_avg[ION][2];
}

static void
brag_calc_update(const gkyl_moment_braginskii *bes,
  const double *brag_d[][GKYL_MAX_SPECIES], double *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  const int ndim = bes->ndim;
  double div_pi[GKYL_MAX_SPECIES][3] = {0.0};
  double div_q[GKYL_MAX_SPECIES] = {0.0};
  if (ndim == 1) {
    const double dx = bes->grid.dx[0];
    double pi[2][GKYL_MAX_SPECIES][6] = {0.0};
    double q[2][GKYL_MAX_SPECIES][3] = {0.0};
    for (int n=0; n < nfluids; ++n) {
      for (int j = L_1D; j <= U_1D; ++j) {
        pi[j][n][0] = brag_d[j][n][PIXX];
        pi[j][n][1] = brag_d[j][n][PIXY];
        pi[j][n][2] = brag_d[j][n][PIXZ];
        pi[j][n][3] = brag_d[j][n][PIYY];
        pi[j][n][4] = brag_d[j][n][PIYZ];
        pi[j][n][5] = brag_d[j][n][PIZZ];
      }
      div_pi[n][0] = calc_sym_grad_1D(dx, pi[L_1D][n][0], pi[U_1D][n][0]);
      div_pi[n][1] = calc_sym_grad_1D(dx, pi[L_1D][n][1], pi[U_1D][n][1]);
      div_pi[n][2] = calc_sym_grad_1D(dx, pi[L_1D][n][2], pi[U_1D][n][2]);
      rhs[n][RHO] = 0.0;
      rhs[n][MX] = -div_pi[n][0];
      rhs[n][MY] = -div_pi[n][1];
      rhs[n][MZ] = -div_pi[n][2];
      // If energy variable exists, increment heat flux and viscous heating
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        for (int j = L_1D; j <= U_1D; ++j) {
          q[j][n][0] = brag_d[j][n][QX];
          q[j][n][1] = brag_d[j][n][QY];
          q[j][n][2] = brag_d[j][n][QZ];
        }
        div_q[n] = calc_sym_grad_1D(dx, q[L_1D][n][0], q[U_1D][n][0]);
        rhs[n][ER] = -div_q[n];
      }
    }
  }
  else if (ndim == 2) {
    const double dx = bes->grid.dx[0];
    const double dy = bes->grid.dx[1];
    double pi[4][GKYL_MAX_SPECIES][6] = {0.0};
    double q[4][GKYL_MAX_SPECIES][3] = {0.0};
    for (int n=0; n < nfluids; ++n) {
      for (int j = LL_2D; j <= UU_2D; ++j) {
        pi[j][n][0] = brag_d[j][n][PIXX];
        pi[j][n][1] = brag_d[j][n][PIXY];
        pi[j][n][2] = brag_d[j][n][PIXZ];
        pi[j][n][3] = brag_d[j][n][PIYY];
        pi[j][n][4] = brag_d[j][n][PIYZ];
        pi[j][n][5] = brag_d[j][n][PIZZ];
      }
      div_pi[n][0] = calc_sym_gradx_2D(dx, pi[LL_2D][n][0], pi[LU_2D][n][0], pi[UL_2D][n][0], pi[UU_2D][n][0])
                       + calc_sym_grady_2D(dy, pi[LL_2D][n][1], pi[LU_2D][n][1], pi[UL_2D][n][1], pi[UU_2D][n][1]);
      
      div_pi[n][1] = calc_sym_gradx_2D(dx, pi[LL_2D][n][1], pi[LU_2D][n][1], pi[UL_2D][n][1], pi[UU_2D][n][1])
                       + calc_sym_grady_2D(dy, pi[LL_2D][n][3], pi[LU_2D][n][3], pi[UL_2D][n][3], pi[UU_2D][n][3]);
      
      div_pi[n][2] = calc_sym_gradx_2D(dx, pi[LL_2D][n][2], pi[LU_2D][n][2], pi[UL_2D][n][2], pi[UU_2D][n][2])
                       + calc_sym_grady_2D(dy, pi[LL_2D][n][4], pi[LU_2D][n][4], pi[UL_2D][n][4], pi[UU_2D][n][4]);
      rhs[n][RHO] = 0.0;
      rhs[n][MX] = -div_pi[n][0];
      rhs[n][MY] = -div_pi[n][1];
      rhs[n][MZ] = -div_pi[n][2];
      // If energy variable exists, increment heat flux and viscous heating
      if (bes->param[n].type_eqn == GKYL_EQN_EULER) {
        for (int j = L_1D; j <= U_1D; ++j) {
          q[j][n][0] = brag_d[j][n][QX];
          q[j][n][1] = brag_d[j][n][QY];
          q[j][n][2] = brag_d[j][n][QZ];
        }
        div_q[n] = calc_sym_gradx_2D(dx, q[LL_2D][n][0], q[LU_2D][n][0], q[UL_2D][n][0], q[UU_2D][n][0])
                 + calc_sym_grady_2D(dy, q[LL_2D][n][1], q[LU_2D][n][1], q[UL_2D][n][1], q[UU_2D][n][1]);
        rhs[n][ER] = -div_q[n];
      }
    }
  }
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
  up->coll_fac = inp.coll_fac;

  return up;
}

void
gkyl_moment_braginskii_advance(const gkyl_moment_braginskii *bes,
  struct gkyl_range brag_vars_range, struct gkyl_range update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *brag_vars[GKYL_MAX_SPECIES], struct gkyl_array *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  int ndim = update_range.ndim;
  long sz[] = { 2, 4, 8 };

  long offsets_vertices[sz[ndim-1]];
  create_offsets_vertices(&brag_vars_range, offsets_vertices);

  long offsets_centers[sz[ndim-1]];
  create_offsets_centers(&update_range, offsets_centers);
  
  const double* fluid_d[sz[ndim-1]][GKYL_MAX_SPECIES];
  const double* em_tot_d[sz[ndim-1]];
  double *brag_vars_d[GKYL_MAX_SPECIES];
  const double* brag_vars_up[sz[ndim-1]][GKYL_MAX_SPECIES];
  double *rhs_d[GKYL_MAX_SPECIES];

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, &brag_vars_range);
  while (gkyl_range_iter_next(&iter_vertex)) {
    
    long linc_vertex = gkyl_range_idx(&brag_vars_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(&update_range, iter_vertex.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]); 
      for (int n=0; n<nfluids; ++n)
        fluid_d[i][n] = gkyl_array_cfetch(fluid[n], linc_center + offsets_vertices[i]);
    }
    for (int n=0; n<nfluids; ++n)
      brag_vars_d[n] = gkyl_array_fetch(brag_vars[n], linc_vertex);

    if (bes->type_brag == GKYL_MAG_BRAG)
      mag_brag_calc_vars(bes, fluid_d, em_tot_d, gkyl_array_fetch(cflrate, linc_center), brag_vars_d);
    else if (bes->type_brag == GKYL_UNMAG_BRAG)
      unmag_brag_calc_vars(bes, fluid_d, gkyl_array_fetch(cflrate, linc_center), brag_vars_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, &update_range);
  while (gkyl_range_iter_next(&iter_center)) {
    
    long linc_vertex = gkyl_range_idx(&brag_vars_range, iter_center.idx);
    long linc_center = gkyl_range_idx(&update_range, iter_center.idx);
    
    for (int i=0; i<sz[ndim-1]; ++i)
      for (int n=0; n<nfluids; ++n)
        brag_vars_up[i][n] = gkyl_array_fetch(brag_vars[n], linc_vertex + offsets_centers[i]);
    
    for (int n=0; n<nfluids; ++n)
      rhs_d[n] = gkyl_array_fetch(rhs[n], linc_center);

    brag_calc_update(bes, brag_vars_up, rhs_d);
  }
}

void
gkyl_moment_braginskii_release(gkyl_moment_braginskii* up)
{
  free(up);
}
