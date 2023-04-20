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
#include <gkyl_dg_calc_sr_vars.h>
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
#include <gkyl_util.h>
#include <gkyl_vlasov.h>

// Definitions of private structs and APIs attached to these objects
// for use in Vlasov app.
// context for use in special relativistic simulations
struct gamma_ctx {
  double mass; // species mass
};

// Labels for lower, upper edge of domain
enum vm_domain_edge { VM_EDGE_LOWER, VM_EDGE_UPPER };

// data for moments
struct vm_species_moment {
  bool use_gpu; // should we use GPU (if present)
  struct gkyl_dg_updater_moment *mcalc; // moment update

  // Special relativistic Vlasov arrays
  struct gkyl_array *p_over_gamma; // array for p/gamma (velocity) 
  struct gkyl_array *gamma; // array for gamma = sqrt(1 + p^2) 
  struct gkyl_array *gamma_inv; // array for gamma = 1.0/sqrt(1 + p^2) 
  struct gkyl_array *V_drift; // bulk fluid velocity (computed from M0*V_drift = M1i with weak division)
  struct gkyl_array *GammaV2; // Gamma^2 = 1/(1 - V_drift^2/c^2), Lorentz boost factor squared from bulk fluid velocity
  struct gkyl_array *GammaV_inv; // Gamma_inv = sqrt(1 - V_drift^2/c^2), inverse Lorentz boost factor from bulk fluid velocity

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
  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

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
  struct vm_species_moment integ_moms; // integrated moments

  double *red_integ_diag; // for reduction on GPU
  gkyl_dynvec integ_diag; // integrated moments reduced across grid

  bool is_first_integ_write_call; // flag for int-moments dynvec written first time

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
  int pkpm_fluid_index; // index of the fluid species being collided with as part of pkpm model
                        // index corresponds to location in fluid_species array (size num_fluid_species)
  struct gkyl_array *m1i_pkpm; // "M1i" in the pkpm model for use in current coupling
                               // Used to copy over fluid variable from pkpm_fluid_species (first three components are momentum)
  struct gkyl_array *g_dist_source; // 2*T_perp/m*G - T_perp/m*(F_0 - F_2)
  struct gkyl_array *F_k_p_1; // k+1 distribution function (first NP components are F_2) 
  struct gkyl_array *F_k_m_1; // k-1 distribution function (first NP components are F_1)

  gkyl_dg_updater_vlasov *slvr; // Vlasov solver 
  struct gkyl_dg_eqn *eqn_vlasov; // Vlasov equation object

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];

  bool has_accel; // flag to indicate there is applied acceleration
  struct gkyl_array *accel; // applied acceleration
  struct gkyl_array *accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *accel_proj; // projector for acceleration
  struct vm_eval_accel_ctx accel_ctx; // context for applied acceleration

  enum gkyl_source_id source_id; // type of source
  struct vm_source src; // applied source
  
  bool has_magB; // flag to indicate Vlasov equation solved along field line
  struct gkyl_array *magB; // magnitude of magnetic field (J = 1/B)
  // host copy for use in IO and projecting
  struct gkyl_array *magB_host;

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
  gkyl_proj_on_basis *ext_em_proj; // projector for external electromagnetic field 
  struct vm_eval_ext_em_ctx ext_em_ctx; // context for external electromagnetic field 

  bool has_app_current; // flag to indicate there is an applied current 
  bool app_current_evolve; // flag to indicate applied current is time dependent
  struct gkyl_array *app_current; // applied current
  struct gkyl_array *app_current_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *app_current_proj; // projector for applied current 
  struct vm_eval_app_current_ctx app_current_ctx; // context for applied current

  struct gkyl_array *bvar; // magnetic field unit vector and tensor (diagnostic and for use in pkpm model)
  struct gkyl_array *ExB; // E x B velocity = E x B/|B|^2 (diagnostic and for use in relativistic pkpm model)
  struct gkyl_array *kappa_inv_b; // b_i/kappa; magnetic field unit vector divided by Lorentz boost factor
                                  // for E x B velocity, b_i/kappa = sqrt((B_i)^2/|B|^2*(1 - |E x B|^2/(c^2 |B|^4)))

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

  double* omegaCfl_ptr;
};

// context for use in computing applied advection
struct vm_eval_advect_ctx { evalf_t advect_func; void *advect_ctx; };

// context for use in computing applied diffusion
struct vm_eval_diffusion_ctx { evalf_t diff_func; void* diff_ctx; };

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

  struct gkyl_array *u; // array for fluid/advection velocity
  struct gkyl_array *p; // array for pressure (used by Euler (1 component) and pkpm Euler (6 components))
  struct gkyl_array *rho_inv; // array for 1/rho
  struct gkyl_array *T_perp_over_m; // array for p_perp/rho = T_perp/m
  struct gkyl_array *T_perp_over_m_inv; // array for (T_perp/m)^-1 
  struct gkyl_array *T_ij; // Temperature tensor for penalization T_ij = 3.0*P_ij/rho

  struct gkyl_array *u_bc_buffer; // buffer for applying BCs to flow
  struct gkyl_array *p_bc_buffer; // buffer for applying BCs to pressure
  
  struct gkyl_array *u_host; // array for host-side fluid/advection velocity (for I/O)
  struct gkyl_array *p_host; // array for host-side pressure (for I/O)
  // pkpm variables
  struct gkyl_array *div_p; // array for divergence of the pressure tensor
  struct gkyl_array *pkpm_accel_vars;  // Acceleration variables for pkpm, pkpm_accel_vars:
                                       // 0: div_b (divergence of magnetic field unit vector)
                                       // 1: bb_grad_u (bb : grad(u))
                                       // 2: p_force (total pressure forces in kinetic equation 1/rho div(p_parallel b_hat) - T_perp/m*div(b)
                                       // 3: p_perp_source (pressure source for higher Laguerre moments -> bb : grad(u) - div(u) - nu + nu rho vth^2/p_perp)
                                       // 4: p_perp_div_b (p_perp/rho*div(b) = T_perp/m*div(b))

  double nuHyp; // Hyper-diffusion coefficient

  struct gkyl_array *D; // array for diffusion tensor
  struct gkyl_array *D_host; // host copy of diffusion tensor

  gkyl_dg_updater_fluid *advect_slvr; // Fluid equation solver
  gkyl_dg_updater_diffusion *diff_slvr; // Fluid equation solver

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];

  // Pointers to updaters that apply BCs to velocity and pressure (and Tij in pkpm)
  bool bc_is_absorb; // boolean for absorbing BCs since 1/rho is undefined in absorbing BCs
  struct gkyl_bc_basic *bc_u_lo[3];
  struct gkyl_bc_basic *bc_u_up[3];
  struct gkyl_bc_basic *bc_p_lo[3];
  struct gkyl_bc_basic *bc_p_up[3];

  // fluid advection
  bool has_advect; // flag to indicate there is advection of fluid equation
  enum gkyl_eqn_type eqn_id; // type of advection (e.g., scalar advection vs. Euler vs. isothermal Euler)
  double param; // Input parameter for fluid species (vt for isothermal Euler, gas_gamma for Euler)

  struct gkyl_dg_bin_op_mem *u_mem; // memory needed in computing flow velocity 
                                    // needed for weak division rho*u = rhou

  // pkpm model
  struct vm_species *pkpm_species; // pointer to coupling species in pkpm model
  int species_index; // index of the kinetic species being coupled to in pkpm model
                     // index corresponds to location in vm_species array (size num_species)

  // applied advection
  struct gkyl_array *advect; // applied advection
  struct gkyl_array *advect_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *advect_proj; // projector for advection
  struct vm_eval_advect_ctx advect_ctx; // context for applied advection

  // advection with another species
  bool advects_with_species; // flag to indicate we are advecting with another species
  struct vm_species *advection_species; // pointer to species we advect with
  struct gkyl_array *other_advect; // pointer to that species drift velocity

  // fluid diffusion
  bool has_diffusion; // flag to indicate there is applied diffusion
  enum gkyl_diffusion_id diffusion_id; // type of diffusion (e.g., isotropic vs. anisotropic)
  gkyl_proj_on_basis* diff_proj; // projector for diffusion
  struct vm_eval_diffusion_ctx diff_ctx; // context for applied diffusion

  // fluid source
  enum gkyl_source_id source_id; // type of source
  struct vm_fluid_source src; // applied source
  
  // collisions with another species present
  enum gkyl_collision_id collision_id; // type of collisions
  struct gkyl_array *other_nu; // pointer to that species collision frequency
  struct gkyl_array *other_m0; // pointer to that species density
  struct gkyl_array *other_nu_vthsq; // pointer to that species nu*vth_sq

  struct gkyl_array *nu_fluid; // collision frequency multiplying fluid_species (nu*nT_perp or nu*nT_z)
  struct gkyl_array *nu_n_vthsq; // nu*n*vthsq (what collisions relax auxiliary temperature to)

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
  
  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *confBasis;
  } basis_on_dev;

  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  bool has_field; // has field
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
 * Apply BCs to primitive moments
 *
 * @param app Vlasov app object
 * @param lbo Species LBO object to apply boundary conditions to primitive moments
 */
void vm_species_lbo_apply_bc(struct gkyl_vlasov_app *app, const struct vm_lbo_collisions *lbo);

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
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param fin Input distribution function
 * @param em Input EM field (needed for bvar, and ExB and kappa_inv_b in relativistic pkpm)
 */
void vm_species_calc_pkpm_vars(gkyl_vlasov_app *app, struct vm_species *species, 
  const struct gkyl_array *fin, const struct gkyl_array *em);

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
 * @param dir Direction to apply BCs
 * @param f Field to apply BCs
 */
void vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f);

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
 * Compute special relativistic electromagnetic variables
 * bvar = magnetic field unit vector (first 3 components) and unit tensor (last 6 components)
 * ExB = E x B velocity, E x B/|B|^2
 * kappa_inv_b = b_i/kappa; magnetic field unit vector divided by Lorentz boost factor
 *               for E x B velocity, b_i/kappa = sqrt((B_i)^2/|B|^2*(1 - |E x B|^2/(c^2 |B|^4)))
 *
 * @param app Vlasov app object
 * @param field Field object (output bvar is stored in field object)
 * @param em Input electromagnetic fields
 */
void vm_field_calc_sr_pkpm_vars(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em);

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
 * @param f Field array
 */
void vm_field_calc_energy(gkyl_vlasov_app *app, double tm, const struct vm_field *field,
  struct gkyl_array *f);

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
 * Compute species applied advection term
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param tm Time for use in advection
 */
void vm_fluid_species_calc_advect(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, double tm);

/**
 * Compute species applied diffusion term
 *
 * @param app Vlasov app object
 * @param fluid_species Fluid Species object
 * @param tm Time for use in advection
 */
void vm_fluid_species_calc_diff(gkyl_vlasov_app* app, struct vm_fluid_species* fluid_species, double tm);

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
 * Apply BCs to primitive variables (bulk velocity, u, and pressure, p, if pressure present)
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 */
void vm_fluid_species_prim_vars_apply_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species);

/**
 * Release resources allocated by fluid species
 *
 * @param app Vlasov app object
 * @param f Fluid_Species object to release
 */
void vm_fluid_species_release(const gkyl_vlasov_app* app, struct vm_fluid_species *f);
