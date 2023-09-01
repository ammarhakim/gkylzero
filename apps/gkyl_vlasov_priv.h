// Private header for use in Vlasov app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_prim_vars.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_dist_vars.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_pkpm.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_updater_fluid.h>
#include <gkyl_dg_updater_diffusion.h>
#include <gkyl_dg_updater_lbo_vlasov.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_dg_vlasov_sr.h>
#include <gkyl_dynvec.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_bcorr_lbo_vlasov_pkpm.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_pkpm.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_pkpm.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_util.h>
#include <gkyl_vlasov.h>

// Definitions of private structs and APIs attached to these objects
// for use in Vlasov app.

// data for moments
struct vm_species_moment {
  bool use_gpu; // should we use GPU (if present)
  struct gkyl_dg_updater_moment *mcalc; // moment update

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

// forward declare species struct
struct vm_species;

struct vm_lbo_collisions {
  enum gkyl_model_id model_id; // type of Vlasov equation (e.g., Vlasov vs. SR vs. PKPM model)
  
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *prim_moms, *nu_prim_moms; // LBO primitive moments
  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  double betaGreenep1; // value of Greene's factor beta + 1
  double other_m[GKYL_MAX_SPECIES]; // masses of species being collided with
  struct gkyl_array *other_prim_moms[GKYL_MAX_SPECIES]; // self-primitive moments of species being collided with
  struct gkyl_array *cross_prim_moms[GKYL_MAX_SPECIES]; // LBO cross-primitive moments
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // LBO cross-species collision frequencies
  struct gkyl_array *other_nu[GKYL_MAX_SPECIES];
  struct gkyl_array *cross_nu_prim_moms; // weak multiplication of collision frequency and primitive moments
  
  struct gkyl_array *self_nu, *self_nu_prim_moms; // LBO self-primitive moments

  struct vm_species_moment moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)
  struct gkyl_array *m0;
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  gkyl_dg_bin_op_mem *greene_factor_mem; // memory needed in computing Greene factor
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct vm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

struct vm_boundary_fluxes {
  struct vm_species_moment integ_moms[2*GKYL_MAX_CDIM]; // integrated moments
  gkyl_ghost_surf_calc *flux_slvr; // boundary flux solver
};

struct vm_source {
  struct vm_species_moment moms; // source moments

  bool calc_bflux; // flag for calculating boundary fluxes
  struct vm_boundary_fluxes bflux; // boundary flux object

  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *source_proj; // projector for source

  struct vm_species *source_species; // species to use for the source
  int source_species_idx; // index of source species
  
  double scale_factor; // factor to scale source function
  double source_length; // length used to scale the source function
  double *scale_ptr;
};


// context for use in computing applied acceleration
struct vm_eval_accel_ctx { evalf_t accel_func; void *accel_ctx; };

// species data
struct vm_species {
  struct gkyl_vlasov_species info; // data for species
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    
  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // fixed buffers for time independent BCs 

  struct gkyl_array *f_host; // host copy for use IO and initialization

  struct vm_species_moment m1i; // for computing currents
  struct vm_species_moment m0; // for computing charge density
  struct vm_species_moment pkpm_moms; // for computing pkpm moments needed in update
  struct vm_species_moment pkpm_moms_diag; // for computing pkpm moments diagnostics
  struct vm_species_moment *moms; // diagnostic moments

  struct gkyl_array *L2_f; // L2 norm f^2
  struct vm_species_moment integ_moms; // integrated moments
  double *red_L2_f; // for reduction of integrated L^2 norm on GPU
  double *red_integ_diag; // for reduction of integrated moments on GPU
  gkyl_dynvec integ_L2_f; // integrated L^2 norm reduced across grid
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_L2_write_call; // flag for integrated L^2 norm dynvec written first time
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  enum gkyl_field_id field_id; // type of field equation 
  enum gkyl_model_id model_id; // type of Vlasov equation (e.g., Vlasov vs. SR)
  struct gkyl_array *qmem; // array for q/m*(E,B)

  // Special relativistic Vlasov arrays
  struct gkyl_array *p_over_gamma; // array for p/gamma (velocity) in special relativistic equation
  struct gkyl_array *p_over_gamma_host; // host copy for use in projecting before copying over to GPU
  struct gkyl_array *gamma; // array for gamma = sqrt(1 + p^2) 
  struct gkyl_array *gamma_host; // host copy for use in projecting before copying over to GPU
  struct gkyl_array *gamma_inv; // array for gamma = 1.0/sqrt(1 + p^2) 
  struct gkyl_array *gamma_inv_host; // host copy for use in projecting before copying over to GPU
  // Special relativistic derived quantities
  struct gkyl_array *V_drift; // bulk fluid velocity (computed from M0*V_drift = M1i with weak division)
  struct gkyl_array *GammaV2; // Gamma^2 = 1/(1 - V_drift^2/c^2), Lorentz boost factor squared from bulk fluid velocity
  struct gkyl_array *GammaV_inv; // Gamma_inv = sqrt(1 - V_drift^2/c^2), inverse Lorentz boost factor from bulk fluid velocity

  struct gkyl_dg_bin_op_mem *V_drift_mem; // memory used in the div-op for V_drift from M1i and M0

  // Data for PKPM model
  struct vm_fluid_species *pkpm_fluid_species; // pointers to cross-species we collide with
  int pkpm_fluid_index; // index of the fluid species being collided with as part of PKPM model
                        // index corresponds to location in fluid_species array (size num_fluid_species)
  // PKPM distribution function variables
  struct gkyl_array *g_dist_source; // g_dist_source = [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)), 
                                    //                 (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ]
  struct gkyl_array *F_k_p_1; // k+1 distribution function (first NP components are F_2) 
  struct gkyl_array *F_k_m_1; // k-1 distribution function (first NP components are F_1)

  // PKPM variables
  struct gkyl_array *m1i_pkpm; // "M1i" in the PKPM model for use in current coupling
                               // Used to copy over fluid variables from pkpm fluid_species, which solves for [rho ux, rho uy, rho uz]
  struct gkyl_array *pkpm_div_ppar; // div(p_parallel b_hat) used for computing self-consistent total pressure force 
  struct gkyl_array *pkpm_prim; // [ux, uy, uz, 1/rho*div(p_par b), T_perp/m, m/T_perp]
  struct gkyl_array *pkpm_prim_surf; // Surface primitive variables. Ordered as:
                                     // [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 3.0*Txx_xl/m, 3.0*Txx_xr/m, 
                                     //  ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 3.0*Tyy_yl/m, 3.0*Tyy_yr/m, 
                                     //  ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr, 3.0*Tzz_zl/m, 3.0*Tzz_zr/m] 
  struct gkyl_array *pkpm_p_ij; // (p_par - p_perp) b_i b_j + p_perp g_ij
  struct gkyl_array *pkpm_p_ij_surf; // (p_par - p_perp) b_i b_j + p_perp g_ij at needed surfaces
                                     // [Pxx_xl, Pxx_xr, Pxy_xl, Pxy_xr, Pxz_xl, Pxz_xr,
                                     //  Pxy_yl, Pxy_yr, Pyy_yl, Pyy_yr, Pyz_yl, Pyz_yr,
                                     //  Pxz_zl, Pxz_zr, Pyz_zl, Pyz_zr, Pzz_zl, Pzz_zr]
  struct gkyl_array *pkpm_lax; // Surface expansion of Lax penalization lambda_i = |u_i| + sqrt(3.0*T_ii/m)
  struct gkyl_array *cell_avg_prim; // Integer array for whether e.g., rho *only* uses cell averages for weak division
                                    // Determined when constructing the matrix if rho or p_perp < 0.0 at control points
  struct gkyl_array *pkpm_accel; // Acceleration variables for PKPM, pkpm_accel:
                                 // 0: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))
                                 // 1: bb_grad_u (bb : grad(u))
                                 // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
                                 // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - 2*nu)
  struct gkyl_array *integ_pkpm_mom; // integrated PKPM variables [rho, rho ux, rho uy, rho uz, rho ux^2, rho uy^2, rho uz^2, p_par, p_perp]
  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars; // Updater to compute PKPM variables (primitive and acceleration variables)
  struct gkyl_dg_calc_pkpm_vars *calc_pkpm_vars_ext; // Updater to compute PKPM variables (primitive and acceleration variables)
                                                     // over extended range (used when BCs are not absorbing to minimize apply BCs calls)
  struct gkyl_dg_calc_pkpm_dist_vars *calc_pkpm_dist_vars; // Updater to compute PKPM distribution function variables 
                                                           // div(p_parallel b_hat) and distribution function sources

  // Pointers for io for PKPM fluid variables, handled by kinetic species because of fluid-kinetic coupling.
  // For PKPM we construct the 10 moment conserved variables for ease of analysis 
  // along with an array of the various update variables, primitive and acceleration
  struct gkyl_array *fluid_io;
  struct gkyl_array *fluid_io_host;
  struct gkyl_array *pkpm_vars_io;
  struct gkyl_array *pkpm_vars_io_host;

  gkyl_dg_updater_vlasov *slvr; // Vlasov solver 
  struct gkyl_dg_eqn *eqn_vlasov; // Vlasov equation object

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  bool bc_is_absorb; // boolean for absorbing BCs since 1/rho is undefined in absorbing BCs
                     // If BCs are *not* absorbing, primitive variables can be calculated on *extended* range 

  bool has_accel; // flag to indicate there is applied acceleration
  struct gkyl_array *accel; // applied acceleration
  struct gkyl_array *accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *accel_proj; // projector for acceleration
  struct vm_eval_accel_ctx accel_ctx; // context for applied acceleration

  enum gkyl_source_id source_id; // type of source
  struct vm_source src; // applied source

  enum gkyl_collision_id collision_id; // type of collisions
  struct vm_lbo_collisions lbo; // collisions object

  double *omegaCfl_ptr;
};

// context for use in computing external electromagnetic fields
struct vm_eval_ext_em_ctx { evalf_t ext_em_func; void *ext_em_ctx; };

// context for use in computing applied current
struct vm_eval_app_current_ctx { evalf_t app_current_func; void *app_current_ctx; };

// field data
struct vm_field {
  struct gkyl_vlasov_field info; // data for field

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *em, *em1, *emnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *em_host;  // host copy for use IO and initialization

  bool has_ext_em; // flag to indicate there is external electromagnetic field
  bool ext_em_evolve; // flag to indicate external electromagnetic field is time dependent
  struct gkyl_array *ext_em; // external electromagnetic field
  struct gkyl_array *ext_em_host; // host copy for use in IO and projecting
  struct gkyl_array *tot_em; // total electromagnetic field
  gkyl_proj_on_basis *ext_em_proj; // projector for external electromagnetic field 
  struct vm_eval_ext_em_ctx ext_em_ctx; // context for external electromagnetic field 

  bool has_app_current; // flag to indicate there is an applied current 
  bool app_current_evolve; // flag to indicate applied current is time dependent
  struct gkyl_array *app_current; // applied current
  struct gkyl_array *app_current_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_current_proj; // projector for applied current 
  struct vm_eval_app_current_ctx app_current_ctx; // context for applied current

  struct gkyl_array *cell_avg_magB2; // Integer array for whether |B|^2 *only* uses cell averages for weak division
                                     // Determined when constructing the matrix if |B|^2 < 0.0 at control points
  struct gkyl_array *bvar; // magnetic field unit vector and tensor (diagnostic and for use in pkpm model)
  struct gkyl_array *ExB; // E x B velocity = E x B/|B|^2 (diagnostic and for use in relativistic pkpm model)
  struct gkyl_array *bvar_surf; // Surface expansion magnetic field unit vector and tensor (for use in pkpm model)
  struct gkyl_array *div_b; // Volume expansion of div(b) (for use in pkpm model)
  struct gkyl_array *max_b; // max(|b_i|) penalization (for use in pkpm model)
  struct gkyl_dg_calc_em_vars *calc_bvar; // Updater to compute magnetic field unit vector and tensor
  struct gkyl_dg_calc_em_vars *calc_ExB; // Updater to compute ExB velocity

  gkyl_hyper_dg *slvr; // Maxwell solver

  struct gkyl_array *em_energy; // EM energy components in each cell
  double *em_energy_red; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components

  bool is_first_energy_write_call; // flag for energy dynvec written first time

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_field_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  double* omegaCfl_ptr;
};

struct vm_fluid_source {
  struct vm_species_moment moms; // source moments

  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *source_proj; // projector for source
};

// fluid species data
struct vm_fluid_species {
  struct gkyl_vlasov_fluid_species info; // data for fluid

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *fluid, *fluid1, *fluidnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)

  struct gkyl_array *fluid_host;  // host copy for use IO and initialization

  enum gkyl_eqn_type eqn_id; // type of fluid system (e.g., scalar advection vs. Euler vs. isothermal Euler)
  double param; // Input parameter for fluid species (vt for isothermal Euler, gas_gamma for Euler)

  // applied advection
  struct gkyl_array *app_advect; // applied advection
  struct gkyl_array *app_advect_host; // host copy for use in IO and projecting

  // Pointers to primitive variables, pressure, and boolean array for if we are only using the cell average for primitive variables
  // For isothermal Euler, prim : (ux, uy, uz), p : (vth*rho)
  // For Euler, prim : (ux, uy, uz, T/m), p : (gamma - 1)*(E - 1/2 rho u^2)
  struct gkyl_array *prim; 
  struct gkyl_array *p; 
  struct gkyl_array *cell_avg_prim; // Integer array for whether e.g., rho *only* uses cell averages for weak division
                                    // Determined when constructing the matrix if rho < 0.0 at control points

  struct vm_species *pkpm_species; // pointer to coupling species in pkpm model
  int species_index; // index of the kinetic species being coupled to in pkpm model
                     // index corresponds to location in vm_species array (size num_species)

  gkyl_dg_updater_fluid *advect_slvr; // Fluid equation solver
  gkyl_dg_updater_diffusion *diff_slvr; // Fluid equation solver

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  bool bc_is_absorb; // boolean for absorbing BCs since 1/rho is undefined in absorbing BCs
                     // If BCs are *not* absorbing, primitive variables can be calculated on *extended* range 

  // fluid diffusion
  bool has_diffusion; // flag to indicate there is applied diffusion
  struct gkyl_array *Dij; // array for diffusion tensor
  struct gkyl_array *Dij_host; // host copy of diffusion tensor
  enum gkyl_diffusion_id diffusion_id; // type of diffusion (e.g., isotropic vs. anisotropic)

  struct gkyl_array *integ_mom; // Integrated moments
  double *red_integ_diag; // for reduction on GPU
  gkyl_dynvec integ_diag; // Integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for int-moments dynvec written first time

  // fluid source
  enum gkyl_source_id source_id; // type of source
  struct vm_fluid_source src; // applied source

  double* omegaCfl_ptr;
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_app {
  char name[128]; // name of app
  struct gkyl_job_pool *job_pool; // Job pool
  
  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  double tcurr; // current time
  double cfl; // CFL number

  bool use_gpu; // should we use GPU (if present)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
    
  struct gkyl_rect_grid grid; // config-space grid
  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges  
  struct gkyl_basis basis, confBasis, velBasis; // phase-space, conf-space basis, vel-space basis

  struct gkyl_comm *comm;   // communicator object for conf-space arrays
  int nghost[GKYL_MAX_CDIM]; // number of ghost-cells in each direction  
  
  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *confBasis;
  } basis_on_dev;

  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  bool has_field; // has field
  bool calc_bvar; // boolean for if simulation needs magnetic field unit vector/tensor
  struct vm_field *field; // pointer to field object

  // species data
  int num_species;
  struct vm_species *species; // data for each species
  
  // fluid data
  int num_fluid_species;
  struct vm_fluid_species *fluid_species; // data for each fluid species
  
  struct gkyl_vlasov_stat stat; // statistics
};

/** gkyl_vlasov_app private API */

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct vm_species* vm_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of species, -1 if not found
 */
int vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm);

/**
 * Find fluid species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of fluid species
 * @return Pointer to fluid species with given name. NULL if not found.o
 */
struct vm_fluid_species *vm_find_fluid_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Return index fluid species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of fluid species
 * @return Index of species, -1 if not found
 */
int vm_find_fluid_species_idx(const gkyl_vlasov_app *app, const char *nm);


/** vm_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param sm Species moment object
 * @param nm Name string indicating moment type
 */
void vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void vm_species_moment_calc(const struct vm_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release species moment object.
 *
 * @param app Vlasov app object
 * @param sm Species moment object to release
 */
void vm_species_moment_release(const struct gkyl_vlasov_app *app,
  const struct vm_species_moment *sm);

/** vm_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 * @param collides_with_fluid Boolean for if kinetic species collides with a fluid species
 */
void vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void vm_species_lbo_cross_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vm_species_lbo_moms(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param collides_with_fluid Boolean for if kinetic species collides with a fluid species
 * @param fluidin Input fluid array (size: num_fluid_species)
 */
void vm_species_lbo_cross_moms(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute RHS from LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 * @return Maximum stable time-step
 */
double vm_species_lbo_rhs(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species LBO object.
 *
 * @param app Vlasov app object
 * @param sm Species LBO object to release
 */
void vm_species_lbo_release(const struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo);

/** vm_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param bflux Species boundary flux object
 */
void vm_species_bflux_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_boundary_fluxes *bflux);

/**
 * Compute boundary flux from rhs
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bflux Species boundary flux object
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void vm_species_bflux_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species boundary flux object.
 *
 * @param app Vlasov app object
 * @param bflux Species boundary flux object to release
 */
void vm_species_bflux_release(const struct gkyl_vlasov_app *app, const struct vm_boundary_fluxes *bflux);

/** vm_species_source API */

/**
 * Initialize species source object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param src Species source object
 */
void vm_species_source_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_source *src);

/**
 * Compute species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param tm Time for use in source
 */
void vm_species_source_calc(gkyl_vlasov_app *app, struct vm_species *species, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function
 */
void vm_species_source_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[]);

/**
 * Release species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vm_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_source *src);

/** vm_species API */

/**
 * Initialize species.
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param s On output, initialized species object
 */
void vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void vm_species_apply_ic(gkyl_vlasov_app *app, struct vm_species *species, double t0);

/**
 * Compute species applied acceleration term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param tm Time for use in acceleration
 */
void vm_species_calc_accel(gkyl_vlasov_app *app, struct vm_species *species, double tm);

/**
 * Compute parallel-kinetic-perpendicular-moment (pkpm) model variables
 * These are the coupling moments [rho, p_par, p_perp], the self-consistent
 * pressure force (div(p_par b_hat)), and the primitive variables
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param fin Input distribution function
 * @param fluidin Input fluid species array (size: num_fluid_species)
 */
void vm_species_calc_pkpm_vars(gkyl_vlasov_app *app, struct vm_species *species, 
  const struct gkyl_array *fin, const struct gkyl_array *fluidin[]);

/**
 * Compute parallel-kinetic-perpendicular-moment (pkpm) model update variables
 * These are the acceleration variables in the kinetic equation and 
 * the source distribution functions for Laguerre couplings.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param fin Input distribution function
 */
void vm_species_calc_pkpm_update_vars(gkyl_vlasov_app *app, struct vm_species *species, const struct gkyl_array *fin);

/**
 * Compute RHS from species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param em EM field
 * @param rhs On output, the RHS from the species object
 * @param fluidin Input fluid array for potential fluid force (size: num_fluid_species)
 * @return Maximum stable time-step
 */
double vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *em, 
  struct gkyl_array *rhs);

/**
 * Apply periodic BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param dir Direction to apply BCs
 * @param f Field to apply BCs
 */
void vm_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f);

/**
 * Apply BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param f Field to apply BCs
 */
void vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f);

/**
 * Compute L2 norm (f^2) of the distribution function diagnostic
 *
 * @param app Vlasov app object
 * @param tm Time at which diagnostic is computed
 * @param species Pointer to species
 */
void vm_species_calc_L2(gkyl_vlasov_app *app, double tm, const struct vm_species *species);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers
 */
void vm_species_coll_tm(gkyl_vlasov_app *app);

/**
 * Fill stat object in app with collisionless timers.
 *
 * @param app App object to update stat timers
 */
void vm_species_tm(gkyl_vlasov_app *app);

/**
 * Delete resources used in species.
 *
 * @param app Vlasov app object
 * @param species Species object to delete
 */
void vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s);

/** vm_field API */

/**
 * Create new field object
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @return Newly created field
 */
struct vm_field* vm_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app);

/**
 * Compute field initial conditions.
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param t0 Time for use in ICs
 */
void vm_field_apply_ic(gkyl_vlasov_app *app, struct vm_field *field, double t0);

/**
 * Compute external electromagnetic fields
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param tm Time for use in external electromagnetic fields computation
 */
void vm_field_calc_ext_em(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Compute applied currents
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param tm Time for use in applied current computation
 */
void vm_field_calc_app_current(gkyl_vlasov_app *app, struct vm_field *field, double tm);

/**
 * Compute magnetic field unit vector and unit tensor
 *
 * @param app Vlasov app object
 * @param field Field object (output bvar is stored in field object)
 * @param em Input electromagnetic fields
 */
void vm_field_calc_bvar(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em);

/**
 * Compute E x B velocity
 *
 * @param app Vlasov app object
 * @param field Field object (output ExB is stored in field object)
 * @param em Input electromagnetic fields
 */
void vm_field_calc_ExB(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em);

/**
 * Accumulate current density onto RHS from field equations
 *
 * @param app Vlasov app object
 * @param fin[] Input distribution function (num_species size)
 * @param fluidin[] Input fluid array (num_fluid_species size)
 * @param emout On output, the RHS from the field solver *with* accumulated current density
 */
void vm_field_accumulate_current(gkyl_vlasov_app *app, 
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], struct gkyl_array *emout);

/**
 * Compute RHS from field equations
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param em Input field
 * @param rhs On output, the RHS from the field solver
 * @return Maximum stable time-step
 */
double vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em, struct gkyl_array *rhs);

/**
 * Apply periodic BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param dir Direction to apply BCs
 * @param f Field to apply BCs
 */
void vm_field_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, struct gkyl_array *f);

/**
 * Apply BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param f Field to apply BCs
 */
void vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  struct gkyl_array *f);

/**
 * Compute field energy diagnostic
 *
 * @param app Vlasov app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 */
void vm_field_calc_energy(gkyl_vlasov_app *app, double tm, const struct vm_field *field);

/**
 * Release resources allocated by field
 *
 * @param app Vlasov app object
 * @param f Field object to release
 */
void vm_field_release(const gkyl_vlasov_app* app, struct vm_field *f);

/** vm_fluid_species_source API */

/**
 * Initialize fluid species source object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param src Species source object
 */
void vm_fluid_species_source_init(struct gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct vm_fluid_source *src);

/**
 * Compute fluid species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param tm Time for use in source
 */
void vm_fluid_species_source_calc(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function RHS
 */
void vm_fluid_species_source_rhs(gkyl_vlasov_app *app, const struct vm_fluid_species *species,
  struct vm_fluid_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[]);

/**
 * Release fluid species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vm_fluid_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_fluid_source *src);

/** vm_fluid_species API */

/**
 * Create new fluid species object
 *
 * @param vm Input VM data
 * @param app Vlasov app object
 * @param f On output, initialized fluid species object
 */
void vm_fluid_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_fluid_species *f);

/**
 * Compute fluid species initial conditions.
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param t0 Time for use in ICs
 */
void vm_fluid_species_apply_ic(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double t0);

/**
 * Compute primitive variables (bulk velocity, u, and pressure, p, if pressure present)
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species (where primitive variables are stored)
 * @param fluid Input array fluid species
 */
void vm_fluid_species_prim_vars(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species,
  const struct gkyl_array *fluid);

/**
 * Compute RHS from fluid species equations
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param fluid Input fluid species
 * @param em EM field
 * @param rhs On output, the RHS from the fluid species solver
 * @return Maximum stable time-step
 */
double vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, 
  const struct gkyl_array *fluid, const struct gkyl_array *em, 
  struct gkyl_array *rhs);

/**
 * Apply periodic BCs to fluid species
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param dir Direction to apply BCs
 * @param f Fluid Species to apply BCs
 */
void vm_fluid_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, 
  int dir, struct gkyl_array *f);

/**
 * Apply BCs to fluid species
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param f Fluid Species to apply BCs
 */
void vm_fluid_species_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species, struct gkyl_array *f);

/**
 * Release resources allocated by fluid species
 *
 * @param app Vlasov app object
 * @param f Fluid_Species object to release
 */
void vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f);
