#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_braginskii_priv.h>
#include <gkyl_prim_euler.h>

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
  double coll_fac; // constant multiplicative factor for collision time to increase or decrease collisionality
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

// Calculate magnetized parallel viscous stress tensor
static void
calc_pi_par(double eta_par, double b_avg[3], double w[6], double pi_par[6])
{ 
  // parallel rate of strain = (bb - 1/3 I) : W
  double par_ros = (b_avg[0]*b_avg[0] - 1.0/3.0)*w[0] + 2.0*b_avg[0]*b_avg[1]*w[1] + 2.0*b_avg[0]*b_avg[2]*w[2] 
                  + (b_avg[1]*b_avg[1] - 1.0/3.0)*w[3] + 2.0*b_avg[1]*b_avg[2]*w[4] 
                  + (b_avg[2]*b_avg[2] - 1.0/3.0)*w[5];

  // pi_par = -eta_par * (bb - 1/3 I) (bb - 1/3 I) : W
  pi_par[0] = -eta_par*(b_avg[0]*b_avg[0] - 1.0/3.0)*par_ros;
  pi_par[1] = -eta_par*b_avg[0]*b_avg[1]*par_ros;
  pi_par[2] = -eta_par*b_avg[0]*b_avg[2]*par_ros;
  pi_par[3] = -eta_par*(b_avg[1]*b_avg[1] - 1.0/3.0)*par_ros;
  pi_par[4] = -eta_par*b_avg[1]*b_avg[2]*par_ros;
  pi_par[5] = -eta_par*(b_avg[2]*b_avg[2] - 1.0/3.0)*par_ros;
}

// Calculate magnetized perpendicular viscous stress tensor
static void
calc_pi_perp(double eta_perp, double b_avg[3], double w[6], double pi_perp[6])
{   
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
  
  // pi_perp = -eta_perp * ((I - bb) . W . (I + 3bb) + (I + 3bb) . W . (I - bb))
  pi_perp[0] = -eta_perp*(2.0*(w[0] + 2.0*b_avg[0]*bWIx - 3.0*b_avg[0]*b_avg[0]*bWb));
  pi_perp[1] = -eta_perp*(2.0*w[1] + 2.0*(b_avg[1]*bWIx + b_avg[0]*bWIy) - 6.0*b_avg[0]*b_avg[1]*bWb);
  pi_perp[2] = -eta_perp*(2.0*w[2] + 2.0*(b_avg[2]*bWIx + b_avg[0]*bWIz) - 6.0*b_avg[0]*b_avg[2]*bWb);
  pi_perp[3] = -eta_perp*(2.0*(w[3] + 2.0*b_avg[1]*bWIy - 3.0*b_avg[1]*b_avg[1]*bWb));
  pi_perp[4] = -eta_perp*(2.0*w[4] + 2.0*(b_avg[2]*bWIy + b_avg[1]*bWIz) - 6.0*b_avg[1]*b_avg[2]*bWb);
  pi_perp[5] = -eta_perp*(2.0*(w[5] + 2.0*b_avg[2]*bWIz - 3.0*b_avg[2]*b_avg[2]*bWb));
}
// Calculate magnetized gyroviscous viscous stress tensor
static void
calc_pi_cross(double eta_cross, double b_avg[3], double w[6], double pi_cross[6])
{
}

// Compute the RHS for the Braginskii transport terms in 1D
// If isothermal Euler, only update momentum
// If Euler, update pressure due to viscous heating and heat conduction
//
// Logic is thus as follows: first compute RHS due to inter-species collisions,
// and viscosity in momentum equation. Then add the corresponding contributions 
// to the pressure variable (if it exists) from heating due to the momentum evolution,
// i.e., Joule heating and viscous heating, and finally add thermal conduction and temperature equilibriation. 
static void
mag_braginskii_update(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES], const double *em_tot_d[],
  double *cflrate, double *rhs[GKYL_MAX_SPECIES])
{
  int nfluids = bes->nfluids;
  const int ndim = bes->ndim;
  if (ndim == 1 && nfluids == 2) {
    const double dx = bes->grid.dx[0];
    const double m[2] = { bes->param[0].mass, bes->param[1].mass };
    const double q[2] = { bes->param[0].charge, bes->param[1].charge };
    // Grab indices of electron and ion fluid arrays
    const int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
    const int ION = (ELC + 1) % 2;

    // Input quantities
    double rho[3][2] = {0.0}; // Mass density for each species
    double u[3][2][3] = {0.0}; // Flow for each species, ux, uy, & uz
    double p[3][2] = {0.0}; // Pressure for each species
    double T[3][2] = {0.0}; // Temperature for each species (m*p/rho)
    double b[3][3] = {0.0}; // Magnetic field unit vector, Bx/|B|, By/|B|, & Bz/|B|

    // Derived quantities
    double omega_c[3][2] = {0.0}; // Cyclotron frequency for each species
    double tau[3][2] = {0.0}; // Collision times for each species
    double eta_par[3][2] = {0.0}; // Parallel viscosity for each species
    double eta_perp[3][2] = {0.0}; // Perpendicular viscosity for each species
    double eta_cross[3][2] = {0.0}; // Gyro-viscosity for each species
    double kappa_par[3][2] = {0.0}; // Parallel conductivity for each species
    double kappa_perp[3][2] = {0.0}; // Perpendicular conductivity for each species
    double kappa_cross[3][2] = {0.0}; // Gyro-conductivity for each species
    double thermal_par[3] = {0.0}; // Parallel thermal force coefficient (same for each species)
    double thermal_perp[3] = {0.0}; // Perpendicular thermal force coefficient (same for each species)
    double current_par[3][3] = {0.0}; // b_hat b_hat dot (u_i - u_e) 
    double current_cross[3][3] = {0.0}; // b x (u_i - u_e)
    for (int j = L_1D; j <= U_1D; ++j)
    {
      // Input density and flow in each cell
      rho[j][ELC] = fluid_d[j][ELC][RHO];
      rho[j][ION] = fluid_d[j][ION][RHO];
      for (int k = MX; k <= MZ; ++k)
      {
        u[j][ELC][k - MX] = fluid_d[j][ELC][k] / rho[j][ELC];
        u[j][ION][k - MX] = fluid_d[j][ION][k] / rho[j][ION];
      }
      // Input magnetic field in each cell
      calc_bhat(em_tot_d[j], b[j]);
      // Derived cyclotron frequency and collision time
      omega_c[j][ELC] = calc_omega_c(q[ELC], m[ELC], em_tot_d[j]);
      omega_c[j][ION] = calc_omega_c(q[ION], m[ION], em_tot_d[j]);
      tau[j][ELC] = calc_tau(1.0, bes->coll_fac, bes->epsilon0, q[ELC], q[ION], m[ELC], m[ION], rho[j][ION], T[j][ELC]);
      tau[j][ION] = calc_tau(1.0, sqrt(2.0)*bes->coll_fac, bes->epsilon0, q[ION], q[ION], m[ION], m[ION], rho[j][ION], T[j][ION]);
      // Pressure information is different for each equation type
      for (int n = 0; n < nfluids; ++n)
      {
        if (bes->param[n].type_eqn == GKYL_EQN_EULER)
          p[j][n] = gkyl_euler_pressure(bes->param[n].p_fac, fluid_d[j][n]); // Euler needs to divide out gas_gamma factor to obtain pressure
        else if (bes->param[n].type_eqn == GKYL_EQN_ISO_EULER)
          p[j][n] = rho[j][n]*bes->param[n].p_fac; // isothermal Euler input is vth, pressure = rho*vth
      }
      T[j][ELC] = m[ELC] * p[j][ELC] / rho[j][ELC];
      eta_par[j][ELC] = 1.5*0.73*p[j][ELC]*tau[j][ELC];
      eta_perp[j][ELC] = 0.51*p[j][ELC]/(tau[j][ELC]*omega_c[j][ELC]*omega_c[j][ELC]);
      eta_cross[j][ELC] = 0.25*p[j][ELC]/omega_c[j][ELC];
      kappa_par[j][ELC] = 3.16*p[j][ELC]*tau[j][ELC]/m[ELC];
      kappa_perp[j][ELC] = 4.66*p[j][ELC]/(m[ELC]*tau[j][ELC]*omega_c[j][ELC]*omega_c[j][ELC]);

      T[j][ION] = m[ION] * p[j][ION] / rho[j][ION];
      eta_par[j][ION] = 1.5*0.96*p[j][ION]*tau[j][ION];
      eta_perp[j][ION] = 0.3*p[j][ION]/(tau[j][ION]*omega_c[j][ION]*omega_c[j][ION]);
      eta_cross[j][ION] = 0.25*p[j][ION]/omega_c[j][ION];
      kappa_par[j][ION] = 3.91*p[j][ION]*tau[j][ION]/m[ION];
      kappa_perp[j][ION] = 2*p[j][ION]/(m[ION]*tau[j][ION]*omega_c[j][ION]*omega_c[j][ION]);
      kappa_cross[j][ION] = 2.5*p[j][ION]/(m[ION]*omega_c[j][ION]);

      thermal_par[j] = -0.71*rho[j][ELC]/m[ELC];
      thermal_perp[j] = 1.5*rho[j][ELC]/(m[ELC]*omega_c[j][ELC]*tau[j][ELC]);
      double b_dot_j = b[j][0]*(u[j][ION][0] - u[j][ELC][0]) + b[j][1]*(u[j][ION][1] - u[j][ELC][1]) + b[j][2]*(u[j][ION][2] - u[j][ELC][2]);
      current_par[j][0] = b[j][0]*b_dot_j;
      current_par[j][1] = b[j][1]*b_dot_j;
      current_par[j][2] = b[j][2]*b_dot_j;
      current_cross[j][0] = b[j][1]*(u[j][ION][2] - u[j][ELC][2]) - b[j][2]*(u[j][ION][1] - u[j][ELC][1]);
      current_cross[j][1] = b[j][2]*(u[j][ION][0] - u[j][ELC][0]) - b[j][0]*(u[j][ION][2] - u[j][ELC][2]);
      current_cross[j][2] = b[j][0]*(u[j][ION][1] - u[j][ELC][1]) - b[j][1]*(u[j][ION][0] - u[j][ELC][0]);
    }

    // Magnetic field at cell edges (using arithmetic average)
    double b_avg_l[3] = {0.0};
    double b_avg_u[3] = {0.0};
    b_avg_l[0] = calc_arithm_avg_1D(b[L_1D][0], b[C_1D][0]);
    b_avg_l[1] = calc_arithm_avg_1D(b[L_1D][1], b[C_1D][1]);
    b_avg_l[2] = calc_arithm_avg_1D(b[L_1D][2], b[C_1D][2]);
    b_avg_u[0] = calc_arithm_avg_1D(b[C_1D][0], b[U_1D][0]);
    b_avg_u[1] = calc_arithm_avg_1D(b[C_1D][1], b[U_1D][1]);
    b_avg_u[2] = calc_arithm_avg_1D(b[C_1D][2], b[U_1D][2]);

    // Parallel thermal force coefficient at cell edges (using arithmetic average)
    double thermal_par_l = calc_arithm_avg_1D(thermal_par[L_1D], thermal_par[C_1D]);
    double thermal_par_u = calc_arithm_avg_1D(thermal_par[C_1D], thermal_par[U_1D]);
    // Perpendicular thermal force coefficient at cell edges (using arithmetic average)
    double thermal_perp_l = calc_arithm_avg_1D(thermal_perp[L_1D], thermal_perp[C_1D]);
    double thermal_perp_u = calc_arithm_avg_1D(thermal_perp[C_1D], thermal_perp[U_1D]);
    
    // Parallel viscosity coefficients at cell edges (using harmonic average)
    double eta_par_l[2] = {0.0};
    double eta_par_u[2] = {0.0};
    eta_par_l[ELC] = calc_harmonic_avg_1D(eta_par[L_1D][ELC], eta_par[C_1D][ELC]);
    eta_par_l[ION] = calc_harmonic_avg_1D(eta_par[L_1D][ION], eta_par[C_1D][ION]);
    eta_par_u[ELC] = calc_harmonic_avg_1D(eta_par[C_1D][ELC], eta_par[U_1D][ELC]);
    eta_par_u[ION] = calc_harmonic_avg_1D(eta_par[C_1D][ION], eta_par[U_1D][ION]);
    // Perpendicular viscosity coefficients at cell edges (using harmonic average)
    double eta_perp_l[2] = {0.0};
    double eta_perp_u[2] = {0.0};
    eta_perp_l[ELC] = calc_harmonic_avg_1D(eta_perp[L_1D][ELC], eta_perp[C_1D][ELC]);
    eta_perp_l[ION] = calc_harmonic_avg_1D(eta_perp[L_1D][ION], eta_perp[C_1D][ION]);
    eta_perp_u[ELC] = calc_harmonic_avg_1D(eta_perp[C_1D][ELC], eta_perp[U_1D][ELC]);
    eta_perp_u[ION] = calc_harmonic_avg_1D(eta_perp[C_1D][ION], eta_perp[U_1D][ION]);    

    // Rate of strain tensor at cell edges for electrons and ions
    double w_l[2][6] = {0.0};
    double w_u[2][6] = {0.0};
    calc_ros_1D(dx, u[L_1D][ELC], u[C_1D][ELC], w_l[ELC]);
    calc_ros_1D(dx, u[L_1D][ION], u[C_1D][ION], w_l[ION]);
    calc_ros_1D(dx, u[C_1D][ELC], u[U_1D][ELC], w_u[ELC]);
    calc_ros_1D(dx, u[C_1D][ION], u[U_1D][ION], w_u[ION]);

    // Parallel viscous stress tensor at cell edges for electrons and ions
    double pi_par_l[2][6] = {0.0};
    double pi_par_u[2][6] = {0.0};
    calc_pi_par(eta_par_l[ELC], b_avg_l, w_l[ELC], pi_par_l[ELC]);
    calc_pi_par(eta_par_l[ION], b_avg_l, w_l[ION], pi_par_l[ION]);
    calc_pi_par(eta_par_u[ELC], b_avg_u, w_u[ELC], pi_par_u[ELC]);
    calc_pi_par(eta_par_u[ION], b_avg_u, w_u[ION], pi_par_u[ION]); 

    // Perpendicular viscous stress tensor at cell edges for electrons and ions
    double pi_perp_l[2][6] = {0.0};
    double pi_perp_u[2][6] = {0.0};
    calc_pi_perp(eta_perp_l[ELC], b_avg_l, w_l[ELC], pi_perp_l[ELC]);
    calc_pi_perp(eta_perp_l[ION], b_avg_l, w_l[ION], pi_perp_l[ION]);
    calc_pi_perp(eta_perp_u[ELC], b_avg_u, w_u[ELC], pi_perp_u[ELC]);
    calc_pi_perp(eta_perp_u[ION], b_avg_u, w_u[ION], pi_perp_u[ION]); 

    // Temperature gradient, parallel, perp, and cross at cell edges
    double gradxT_l[2] = {0.0};
    double gradxT_u[2] = {0.0};
    // Parallel temperature gradient (b_hat b_hat dot grad T)
    double bbgradT_l[2][3] = {0.0};
    double bbgradT_u[2][3] = {0.0};
    // Perpendicular temperature gradient (grad T - b_hat b_hat dot grad T)
    double perp_gradT_l[2][3] = {0.0};
    double perp_gradT_u[2][3] = {0.0};
    // Cross temperature gradient (b x grad T)
    double cross_gradT_l[2][3] = {0.0};
    double cross_gradT_u[2][3] = {0.0};
    
    // Only electrons needed for friction force; ions may be different equation system
    gradxT_l[ELC] = calc_sym_grad_1D(dx, T[L_1D][ELC], T[C_1D][ELC]);
    gradxT_u[ELC] = calc_sym_grad_1D(dx, T[C_1D][ELC], T[U_1D][ELC]);
    
    bbgradT_l[ELC][0] = b_avg_l[0]*b_avg_l[0]*gradxT_l[ELC];
    bbgradT_l[ELC][1] = b_avg_l[1]*b_avg_l[0]*gradxT_l[ELC];
    bbgradT_l[ELC][2] = b_avg_l[2]*b_avg_l[0]*gradxT_l[ELC];
    
    bbgradT_u[ELC][0] = b_avg_u[0]*b_avg_u[0]*gradxT_u[ELC];
    bbgradT_u[ELC][1] = b_avg_u[1]*b_avg_u[0]*gradxT_u[ELC];
    bbgradT_u[ELC][2] = b_avg_u[2]*b_avg_u[0]*gradxT_u[ELC];
    
    perp_gradT_l[ELC][0] = gradxT_l[ELC] - b_avg_l[0]*b_avg_l[0]*gradxT_l[ELC];
    perp_gradT_u[ELC][0] = gradxT_u[ELC] - b_avg_u[0]*b_avg_u[0]*gradxT_u[ELC];
    
    cross_gradT_l[ELC][1] = b_avg_l[2]*gradxT_l[ELC];
    cross_gradT_l[ELC][2] = -b_avg_l[1]*gradxT_l[ELC];
    
    cross_gradT_u[ELC][1] = b_avg_l[2]*gradxT_u[ELC];
    cross_gradT_u[ELC][2] = -b_avg_l[1]*gradxT_u[ELC];
      
    double div_pi[2][3] = {0.0};
    double perp_current[3] = {0.0};
    double par_force_thermal[3] = {0.0};
    double perp_force_thermal[3] = {0.0};
    double friction_force[3] = {0.0};

    // Compute momentum update -div(Pi) +/- F_friction (+ for electron, - for ion)
    div_pi[ELC][0] = calc_sym_grad_1D(dx, pi_par_l[ELC][0]+pi_perp_l[ELC][0], pi_par_u[ELC][0]+pi_perp_u[ELC][0]);
    div_pi[ELC][1] = calc_sym_grad_1D(dx, pi_par_l[ELC][1]+pi_perp_l[ELC][1], pi_par_u[ELC][1]+pi_perp_u[ELC][1]);
    div_pi[ELC][2] = calc_sym_grad_1D(dx, pi_par_l[ELC][2]+pi_perp_l[ELC][2], pi_par_u[ELC][2]+pi_perp_u[ELC][2]);

    div_pi[ION][0] = calc_sym_grad_1D(dx, pi_par_l[ION][0]+pi_perp_l[ION][0], pi_par_u[ION][0]+pi_perp_u[ION][0]);
    div_pi[ION][1] = calc_sym_grad_1D(dx, pi_par_l[ION][1]+pi_perp_l[ION][1], pi_par_u[ION][1]+pi_perp_u[ION][1]);
    div_pi[ION][2] = calc_sym_grad_1D(dx, pi_par_l[ION][2]+pi_perp_l[ION][2], pi_par_u[ION][2]+pi_perp_u[ION][2]);

    // Perpendicular current j_perp/n = (u_i - u_e) - b_hat (bhat dot (u_i - u_e))
    // Note: density factors are included in friction force coefficient    
    perp_current[0] = u[C_1D][ION][0] - u[C_1D][ELC][0] - current_par[C_1D][0];
    perp_current[1] = u[C_1D][ION][1] - u[C_1D][ELC][1] - current_par[C_1D][1];
    perp_current[2] = u[C_1D][ION][2] - u[C_1D][ELC][2] - current_par[C_1D][2];

    par_force_thermal[0] = 0.5*(thermal_par_l*bbgradT_l[ELC][0] + thermal_par_u*bbgradT_u[ELC][0]);
    par_force_thermal[1] = 0.5*(thermal_par_l*bbgradT_l[ELC][1] + thermal_par_u*bbgradT_u[ELC][1]);
    par_force_thermal[2] = 0.5*(thermal_par_l*bbgradT_l[ELC][2] + thermal_par_u*bbgradT_u[ELC][2]);

    perp_force_thermal[0] = 0.0;
    perp_force_thermal[1] = 0.5*(thermal_perp_l*cross_gradT_l[ELC][1] + thermal_perp_u*cross_gradT_u[ELC][1]);
    perp_force_thermal[2] = 0.5*(thermal_perp_l*cross_gradT_l[ELC][2] + thermal_perp_u*cross_gradT_u[ELC][2]);

    // F_friction = j_parallel/sigma_parallel + j_perp/sigma_perp - b_hat grad_par T - b_hat x grad T
    // Thermal force signs included in respective coefficients thermal_par/perp
    friction_force[0] = 0.51*rho[C_1D][ELC]/tau[C_1D][ELC]*current_par[C_1D][0] + rho[C_1D][ELC]/tau[C_1D][ELC]*perp_current[0] + par_force_thermal[0] + perp_force_thermal[0];
    friction_force[1] = 0.51*rho[C_1D][ELC]/tau[C_1D][ELC]*current_par[C_1D][1] + rho[C_1D][ELC]/tau[C_1D][ELC]*perp_current[1] + par_force_thermal[1] + perp_force_thermal[1];
    friction_force[2] = 0.51*rho[C_1D][ELC]/tau[C_1D][ELC]*current_par[C_1D][2] + rho[C_1D][ELC]/tau[C_1D][ELC]*perp_current[2] + par_force_thermal[2] + perp_force_thermal[2];
        
    rhs[ELC][RHO] = 0;
    rhs[ELC][MX] = -div_pi[ELC][0] + friction_force[0];
    rhs[ELC][MY] = -div_pi[ELC][1] + friction_force[1];
    rhs[ELC][MZ] = -div_pi[ELC][2] + friction_force[2];

    rhs[ION][RHO] = 0;
    rhs[ION][MX] = -div_pi[ION][0] - friction_force[0];
    rhs[ION][MY] = -div_pi[ION][1] - friction_force[1];
    rhs[ION][MZ] = -div_pi[ION][2] - friction_force[2];
    
    // Add contributions to the pressure variable if it exists
    if (bes->param[ELC].type_eqn == GKYL_EQN_EULER)
    {
      // Parallel conductivity coefficients at cell edges (using harmonic average)
      double kappa_par_l = calc_harmonic_avg_1D(kappa_par[L_1D][ELC], kappa_par[C_1D][ELC]);
      double kappa_par_u = calc_harmonic_avg_1D(kappa_par[C_1D][ELC], kappa_par[U_1D][ELC]);
      // Perpendicular conductivity coefficients at cell edges (using harmonic average)
      double kappa_perp_l = calc_harmonic_avg_1D(kappa_perp[L_1D][ELC], kappa_perp[C_1D][ELC]);
      double kappa_perp_u = calc_harmonic_avg_1D(kappa_perp[C_1D][ELC], kappa_perp[U_1D][ELC]); 

      // Parallel current multiplied by electron temperature at cell edges (using arithmetic average)
      double current_par_l[3] = {calc_arithm_avg_1D(T[L_1D][ELC]*current_par[L_1D][0], T[C_1D][ELC]*current_par[C_1D][0]),
                                 calc_arithm_avg_1D(T[L_1D][ELC]*current_par[L_1D][1], T[C_1D][ELC]*current_par[C_1D][1]),
                                 calc_arithm_avg_1D(T[L_1D][ELC]*current_par[L_1D][2], T[C_1D][ELC]*current_par[C_1D][2]) };
      double current_par_u[3] = {calc_arithm_avg_1D(T[C_1D][ELC]*current_par[C_1D][0], T[U_1D][ELC]*current_par[U_1D][0]),
                                 calc_arithm_avg_1D(T[C_1D][ELC]*current_par[C_1D][1], T[U_1D][ELC]*current_par[U_1D][1]),
                                 calc_arithm_avg_1D(T[C_1D][ELC]*current_par[C_1D][2], T[U_1D][ELC]*current_par[U_1D][2]) };
      // Cross current multiplied by electron temperature at cell edges (using arithmetic average)
      double current_cross_l[3] = {calc_arithm_avg_1D(T[L_1D][ELC]*current_cross[L_1D][0], T[C_1D][ELC]*current_cross[C_1D][0]),
                                   calc_arithm_avg_1D(T[L_1D][ELC]*current_cross[L_1D][1], T[C_1D][ELC]*current_cross[C_1D][1]),
                                   calc_arithm_avg_1D(T[L_1D][ELC]*current_cross[L_1D][2], T[C_1D][ELC]*current_cross[C_1D][2]) };
      double current_cross_u[3] = {calc_arithm_avg_1D(T[C_1D][ELC]*current_cross[C_1D][0], T[U_1D][ELC]*current_cross[U_1D][0]),
                                   calc_arithm_avg_1D(T[C_1D][ELC]*current_cross[C_1D][1], T[U_1D][ELC]*current_cross[U_1D][1]),
                                   calc_arithm_avg_1D(T[C_1D][ELC]*current_cross[C_1D][2], T[U_1D][ELC]*current_cross[U_1D][2]) };
      
      double q_par_l[3] = {kappa_par_l*bbgradT_l[ELC][0] + thermal_par_l*current_par_l[0],
                           kappa_par_l*bbgradT_l[ELC][1] + thermal_par_l*current_par_l[1],
                           kappa_par_l*bbgradT_l[ELC][2] + thermal_par_l*current_par_l[2] };
      double q_par_u[3] = {kappa_par_u*bbgradT_u[ELC][0] + thermal_par_u*current_par_u[0],
                           kappa_par_u*bbgradT_u[ELC][1] + thermal_par_u*current_par_u[1],
                           kappa_par_u*bbgradT_u[ELC][2] + thermal_par_u*current_par_u[2] };
      double q_perp_l[3] = {kappa_perp_l*perp_gradT_l[ELC][0] + thermal_perp_l*current_cross_l[0],
                            kappa_perp_l*perp_gradT_l[ELC][1] + thermal_perp_l*current_cross_l[1],
                            kappa_perp_l*perp_gradT_l[ELC][2] + thermal_perp_l*current_cross_l[2] };
      double q_perp_u[3] = {kappa_perp_u*perp_gradT_u[ELC][0] + thermal_perp_u*current_cross_u[0],
                            kappa_perp_u*perp_gradT_u[ELC][1] + thermal_perp_u*current_cross_u[1],
                            kappa_perp_u*perp_gradT_u[ELC][2] + thermal_perp_u*current_cross_u[2] };
      
      double div_q[3] = {calc_sym_grad_1D(dx, q_par_l[0] + q_perp_l[0], q_par_u[0] + q_perp_u[0]), 0.0, 0.0};
      double friction_heating = (u[C_1D][ION][0] - u[C_1D][ELC][0])*friction_force[0]
        + (u[C_1D][ION][1] - u[C_1D][ELC][1])*friction_force[1]
        + (u[C_1D][ION][2] - u[C_1D][ELC][2])*friction_force[2]
        - 3.0*m[ELC]/m[ION]*rho[C_1D][ELC]*(T[C_1D][ELC] - T[C_1D][ION])/tau[C_1D][ELC];
      double viscous_heating[6] = {calc_arithm_avg_1D(pi_par_l[ELC][0]+pi_perp_l[ELC][0], pi_par_u[ELC][0]+pi_perp_u[ELC][0]),
                                   calc_arithm_avg_1D(pi_par_l[ELC][1]+pi_perp_l[ELC][1], pi_par_u[ELC][1]+pi_perp_u[ELC][1]),
                                   calc_arithm_avg_1D(pi_par_l[ELC][2]+pi_perp_l[ELC][2], pi_par_u[ELC][0]+pi_perp_u[ELC][2]),
                                   0.0, 0.0, 0.0};
      rhs[ELC][ER] = -div_q[0] - div_q[1] - div_q[2]
        - viscous_heating[0] - viscous_heating[1] - viscous_heating[2] - viscous_heating[3] - viscous_heating[4] - viscous_heating[5]
        + friction_heating;
    }
    if (bes->param[ION].type_eqn == GKYL_EQN_EULER)
    {
      // Parallel conductivity coefficients at cell edges (using harmonic average)
      double kappa_par_l = calc_harmonic_avg_1D(kappa_par[L_1D][ION], kappa_par[C_1D][ION]);
      double kappa_par_u = calc_harmonic_avg_1D(kappa_par[C_1D][ION], kappa_par[U_1D][ION]);
      // Perpendicular conductivity coefficients at cell edges (using harmonic average)
      double kappa_perp_l = calc_harmonic_avg_1D(kappa_perp[L_1D][ION], kappa_perp[C_1D][ION]);
      double kappa_perp_u = calc_harmonic_avg_1D(kappa_perp[C_1D][ION], kappa_perp[U_1D][ION]);

      // Ion Temperature gradient, parallel, perp and cross at cell edges
      gradxT_l[ION] = calc_sym_grad_1D(dx, T[L_1D][ION], T[C_1D][ION]);
      gradxT_u[ION] = calc_sym_grad_1D(dx, T[C_1D][ION], T[U_1D][ION]);
      
      bbgradT_l[ION][0] = b_avg_l[0]*b_avg_l[0]*gradxT_l[ION];
      bbgradT_l[ION][1] = b_avg_l[1]*b_avg_l[0]*gradxT_l[ION];
      bbgradT_l[ION][2] = b_avg_l[2]*b_avg_l[0]*gradxT_l[ION];
    
      bbgradT_u[ION][0] = b_avg_u[0]*b_avg_u[0]*gradxT_u[ION];
      bbgradT_u[ION][1] = b_avg_u[1]*b_avg_u[0]*gradxT_u[ION];
      bbgradT_u[ION][2] = b_avg_u[2]*b_avg_u[0]*gradxT_u[ION];
    
      perp_gradT_l[ION][0] = gradxT_l[ION] - b_avg_l[0]*b_avg_l[0]*gradxT_l[ION];
      perp_gradT_u[ION][0] = gradxT_u[ION] - b_avg_u[0]*b_avg_u[0]*gradxT_u[ION];
    
      cross_gradT_l[ION][1] = b_avg_l[2]*gradxT_l[ION];
      cross_gradT_l[ION][2] = -b_avg_l[1]*gradxT_l[ION];
    
      cross_gradT_u[ION][1] = b_avg_l[2]*gradxT_u[ION];
      cross_gradT_u[ION][2] = -b_avg_l[1]*gradxT_u[ION];

      double q_par_l[3] = {kappa_par_l*bbgradT_l[ION][0], kappa_par_l*bbgradT_l[ION][1], kappa_par_l*bbgradT_l[ION][2] };
      double q_par_u[3] = {kappa_par_u*bbgradT_u[ION][0], kappa_par_u*bbgradT_u[ION][1], kappa_par_u*bbgradT_u[ION][2] };
      double q_perp_l[3] = {kappa_perp_l*perp_gradT_l[ION][0], kappa_perp_l*perp_gradT_l[ION][1], kappa_perp_l*perp_gradT_l[ION][2] };
      double q_perp_u[3] = {kappa_perp_u*perp_gradT_u[ION][0], kappa_perp_u*perp_gradT_u[ION][1], kappa_perp_u*perp_gradT_u[ION][2] };
      
      double div_q[3] = {calc_sym_grad_1D(dx, q_par_l[0] + q_perp_l[0], q_par_u[0] + q_perp_u[0]), 0.0, 0.0};
      // Temperature equilibriation has opposite sign for ions
      double friction_heating = 3.0*m[ELC]/m[ION]*rho[C_1D][ELC]*(T[C_1D][ELC] - T[C_1D][ION])/tau[C_1D][ELC];
      double viscous_heating[6] = {calc_arithm_avg_1D(pi_par_l[ION][0]+pi_perp_l[ION][0], pi_par_u[ION][0]+pi_perp_u[ION][0]),
                                   calc_arithm_avg_1D(pi_par_l[ION][1]+pi_perp_l[ION][1], pi_par_u[ION][1]+pi_perp_u[ION][1]),
                                   calc_arithm_avg_1D(pi_par_l[ION][2]+pi_perp_l[ION][2], pi_par_u[ION][0]+pi_perp_u[ION][2]),
                                   0.0, 0.0, 0.0};
      rhs[ION][ER] = -div_q[0] - div_q[1] - div_q[2]
        - viscous_heating[0] - viscous_heating[1] - viscous_heating[2] - viscous_heating[3] - viscous_heating[4] - viscous_heating[5]
        + friction_heating;
    }
  }
  else if (ndim == 2 && nfluids == 2) {
    const double dx = bes->grid.dx[0];
    const double dy = bes->grid.dx[1];
    const double m[2] = { bes->param[0].mass, bes->param[1].mass };
  }
}

static void
unmag_braginskii_update(const gkyl_moment_braginskii *bes,
  const double *fluid_d[][GKYL_MAX_SPECIES],
  double *cflrate, double *rhs[GKYL_MAX_SPECIES])
{
  const int ndim = bes->ndim;
  const int nfluids = bes->nfluids;
  const double dx = bes->grid.dx[0];

  // For now we will asume one ion species
  if (ndim == 1 && nfluids == 2)
  {
    const double m[2] = { bes->param[0].mass, bes->param[1].mass };

    double rho[3][2];
    for (int j = L_1D; j <= U_1D; ++j)
      for (int n = 0; n < nfluids; ++n)
        rho[j][n] = fluid_d[j][n][RHO];

    double v[3][2][3];
    for (int j = L_1D; j <= U_1D; ++j)
      for (int n = 0; n < nfluids; ++n)
        for (int k = MX; k <= MZ; ++k)
          v[j][n][k - MX] = fluid_d[j][n][k] / fluid_d[j][n][RHO];

    double p[3][2];
    const double gas_gamma = 5.0 / 3.0;  // Note: Hard-coded value
    for (int j = L_1D; j <= U_1D; ++j)
      for (int n = 0; n < nfluids; ++n)
      {
        double KE = 0.5 * (fluid_d[j][n][MX] * fluid_d[j][n][MX] + fluid_d[j][n][MY] * fluid_d[j][n][MY] + fluid_d[j][n][MZ] * fluid_d[j][n][MZ]) / fluid_d[j][n][RHO];
        p[j][n] = (gas_gamma - 1.0) * fluid_d[j][n][ER] - KE;
      }

    double T[3][2];
    for (int j = L_1D; j <= U_1D; ++j)
      for (int n = 0; n < nfluids; ++n)
        T[j][n] = m[n] * p[j][n] / rho[j][n];

    // Grab indices of electron and ion fluid arrays
    const int ELC = bes->param[0].charge < 0.0 ? 0 : 1;
    const int ION = (ELC + 1) % 2;

    const double e = bes->param[ELC].charge;
    const int Z = bes->param[ION].charge / fabs(bes->param[ELC].charge);  // Note: May always round down on conversion

    double lambda = 10.0; // Note: Change this

    // Collision times for each species
    double tau[3][2];
    const double e4Z2 = e * e * e * e * Z * Z;
    const double e4Z4 = e4Z2 * Z * Z;
    for (int j = L_1D; j <= U_1D; ++j)
    {
      tau[j][ELC] = 3.0 * sqrt(m[ELC]) * m[ION] * pow(T[j][ELC], 3.0 / 2.0) / (4.0 * sqrt(2.0 * M_PI) * lambda * e4Z2 * rho[j][ION]);
      tau[j][ION] = 3.0 * sqrt(m[ION]) * m[ION] * pow(T[j][ION], 3.0 / 2.0) / (4.0 * sqrt(M_PI) * lambda * e4Z4 * rho[j][ION]);
    }

    // Viscosity coefficients for each species
    double eta[3][2];
    for (int j = L_1D; j <= U_1D; ++j)
    {
      eta[j][ELC] = 0.73 * p[j][ELC] * tau[j][ELC];
      eta[j][ION] = 0.96 * p[j][ION] * tau[j][ION];
    }

    // Calculate rhs
    for (int n = 0; n < nfluids; ++n)
    {
      // Compute velocity derivatives at left and right edges
      double dvdxL[3] = { calc_sym_grad_1D(dx, v[L_1D][n][0], v[C_1D][n][0]),
                          calc_sym_grad_1D(dx, v[L_1D][n][1], v[C_1D][n][1]),
                          calc_sym_grad_1D(dx, v[L_1D][n][2], v[C_1D][n][2]) };
      double dvdxR[3] = { calc_sym_grad_1D(dx, v[C_1D][n][0], v[U_1D][n][0]),
                          calc_sym_grad_1D(dx, v[C_1D][n][1], v[U_1D][n][1]),
                          calc_sym_grad_1D(dx, v[C_1D][n][2], v[U_1D][n][2]) };

      // Compute eta at left and right edges using harmonic average
      double etaL = calc_harmonic_avg_1D(eta[L_1D][n], eta[C_1D][n]);
      double etaR = calc_harmonic_avg_1D(eta[C_1D][n], eta[U_1D][n]);

      // Viscous stress tensor
      double piL[3] = { (4.0 / 3.0) * dvdxL[0], dvdxL[1], dvdxL[2] };
      double piR[3] = { (4.0 / 3.0) * dvdxR[0], dvdxR[1], dvdxR[2] };

      // Update momentum
      rhs[n][MX] = calc_sym_grad_1D(dx, etaL * piL[0], etaR * piR[0]);
      rhs[n][MY] = calc_sym_grad_1D(dx, etaL * piL[1], etaR * piR[1]);
      rhs[n][MZ] = calc_sym_grad_1D(dx, etaL * piL[2], etaR * piR[2]);

      // Compute velocity derivatives at center
      double dvdxC[3] = { calc_harmonic_avg_1D(dvdxL[0], dvdxR[0]),
                          calc_harmonic_avg_1D(dvdxL[1], dvdxR[1]),
                          calc_harmonic_avg_1D(dvdxL[2], dvdxR[2]) };

      // Update energy
      rhs[n][ER] = -1.0 * eta[C_1D][n] * (4.0 / 3.0 * dvdxC[0] * dvdxC[0] + dvdxC[1] * dvdxC[1] + dvdxC[2] * dvdxC[2]);
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
gkyl_moment_braginskii_advance(const gkyl_moment_braginskii *bes, struct gkyl_range update_range,
  struct gkyl_array *fluid[GKYL_MAX_SPECIES], const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *rhs[GKYL_MAX_SPECIES])
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

    if (bes->type_brag == GKYL_MAG_BRAG || bes->type_brag == GKYL_MAG_BRAG_NO_CROSS)
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
