// Private header for use in Vlasov app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_basic.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_prim_vars.h>
#include <gkyl_dg_updater_lbo_vlasov.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_dg_updater_vlasov_poisson.h>
#include <gkyl_dg_vlasov_poisson.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_poisson.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>

// Definitions of private structs and APIs attached to these objects
// for use in Vlasov app.

// Meta-data for IO
struct vlasov_poisson_output_meta {
  int frame; // frame number
  double stime; // output time
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
  char basis_type_nm[64]; // used during read
};

// list of valid moment names
static const char *const valid_moment_names[] = {
  "M0",
  "M1i",
  "M2ij",
  "M2",
  "M3i",
  "M3ijk",
  "FiveMoments", // M0, M1i, M2.
  "LTEMoments", // this is an internal flag for computing moments (n, V_drift, T/m)
                // of the LTE (local thermodynamic equilibrium) distribution
  "Integrated", // this is an internal flag, not for passing to moment type
};

// check if name of moment is valid or not
static bool
is_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_moment_names)/sizeof(valid_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// data for moments
struct vp_species_moment {
  bool is_integrated; // boolean for if computing integrated moments
  int num_mom; // number of moments

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
  // Options for moment calculation: 
  // 1. Compute the moment directly with dg_updater_moment
  // 2. Compute the moments of the equivalent LTE (local thermodynamic equilibrium)
  //    distribution (n, V_drift, T/m) with specialized updater
  union {
    struct {
      struct gkyl_vlasov_lte_moments *vlasov_lte_moms; // updater for computing LTE moments
    };
    struct {
      struct gkyl_dg_updater_moment *mcalc; // moment update
    };
  };

  bool is_vlasov_lte_moms;
};

// forward declare species struct
struct vp_species;

struct vp_lbo_collisions {  
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *prim_moms, *nu_prim_moms; // LBO primitive moments
  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  double self_nu_fac; // Self collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double cross_nu_fac[GKYL_MAX_SPECIES]; // Cross collision freqs without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double vtsq_min; // minimum vtsq
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

  struct vp_species_moment moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)

  struct gkyl_array *m0;
  struct gkyl_array *vtsq;
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  gkyl_dg_bin_op_mem *greene_factor_mem; // memory needed in computing Greene factor
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct vp_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

struct vp_bgk_collisions {  
  struct gkyl_array *nu_sum; // BGK collision frequency 
  struct gkyl_array *nu_sum_host; // BGK collision frequency host-side for I/O
  struct gkyl_array *self_nu; // BGK self-collision frequency

  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  struct gkyl_array *f_lte;
  struct gkyl_array *nu_f_lte;

  enum gkyl_model_id model_id;
  struct vp_species_moment moms; // moments needed in BGK (n, V_drift, T/m) for LTE distribution

  // LTE distribution function projection object
  // also corrects the density of projected distribution function
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte; 

  // Correction updater for insuring LTE distribution has desired LTE (n, V_drift, T/m) moments
  bool correct_all_moms; // boolean if we are correcting all the moments
  struct gkyl_vlasov_lte_correct *corr_lte; 
  gkyl_dynvec corr_stat;
  bool is_first_corr_status_write_call;

  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep)
};

struct vp_boundary_fluxes {
  struct vp_species_moment integ_moms[2*GKYL_MAX_CDIM]; // integrated moments
  gkyl_ghost_surf_calc *flux_slvr; // boundary flux solver
};

struct vp_proj {
  enum gkyl_projection_id proj_id; // type of projection

  // organization of the different projection objects and the required data and solvers
  union {
    // function projection
    struct {
      struct gkyl_proj_on_basis *proj_func; // projection operator for specified function
      struct gkyl_array *proj_host; // array for projection on host-side if running on GPUs
    };
    // LTE (Local thermodynamic equilibrium) distribution function project with moment correction
    // (Maxwellian for non-relativistic, Maxwell-Juttner for relativistic)
    struct {
      struct gkyl_array *dens; // host-side density
      struct gkyl_array *V_drift; // host-side V_drift
      struct gkyl_array *T_over_m; // host-side T/m (temperature/mass)

      struct gkyl_array *vlasov_lte_moms_host; // host-side LTE moms (n, V_drift, T/m)
      struct gkyl_array *vlasov_lte_moms; // LTE moms (n, V_drift, T/m) for passing to updaters

      struct gkyl_proj_on_basis *proj_dens; // projection operator for density
      struct gkyl_proj_on_basis *proj_V_drift; // projection operator for V_drift
      struct gkyl_proj_on_basis *proj_temp; // projection operator for temperature
      
      // LTE distribution function projection object
      // also corrects the density of projected distribution function
      struct gkyl_vlasov_lte_proj_on_basis *proj_lte; 

      // Correction updater for insuring LTE distribution has desired LTE (n, V_drift, T/m) moments
      bool correct_all_moms; // boolean if we are correcting all the moments
      struct gkyl_vlasov_lte_correct *corr_lte;    
    };
  };
};

struct vp_source {
  struct vp_species_moment moms; // source moments

  bool calc_bflux; // flag for calculating boundary fluxes
  struct vp_boundary_fluxes bflux; // boundary flux object

  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  struct vp_proj proj_source; // projector for source

  struct vp_species *source_species; // species to use for the source
  int source_species_idx; // index of source species
  
  double scale_factor; // factor to scale source function
  double source_length; // length used to scale the source function
  double *scale_ptr;
};

// species data
struct vp_species {
  struct gkyl_vlasov_poisson_species info; // data for species

  enum gkyl_vpfield_id field_id; // type of field equation 
  enum gkyl_vpmodel_id model_id; // type of Vlasov equation (e.g., Vlasov vs. SR)
  
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

  double qbym; // Charge (q) divided by mass (m).
  struct gkyl_array *qmem; // array for (q/m)*(phi_tot,A_ext), phi_tot=phi+phi_ext

  struct vp_species_moment m0; // for computing charge density
  struct vp_species_moment integ_moms; // integrated moments
  struct vp_species_moment *moms; // diagnostic moments
  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  gkyl_dg_updater_vlasov_poisson *slvr; // Vlasov solver 
  struct gkyl_dg_eqn *eqn_vlasov; // Vlasov equation object

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
  bool bc_is_np[3]; // whether BC is nonperiodic.
  
  // boundary conditions on lower/upper edges in each direction  
  struct gkyl_vlasov_poisson_bc lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  struct vp_proj proj_init; // projector for initial conditions

  enum gkyl_source_id source_id; // type of source
  struct vp_source src; // applied source

  // Boundary fluxes.
  struct vp_boundary_fluxes bflux;

  // collisions
  struct {
    enum gkyl_collision_id collision_id; // type of collisions
    union {
      struct vp_lbo_collisions lbo; // LBO collisions object
      struct vp_bgk_collisions bgk; // BGK collisions object
    };
  }; 

  double *omega_cfl;
};

// context for use in computing external electromagnetic fields
struct vp_eval_ext_em_ctx { evalf_t ext_em_func; void *ext_em_ctx; };

// context for use in computing applied current
struct vp_eval_app_current_ctx { evalf_t app_current_func; void *app_current_ctx; };

// field data
struct vp_field {
  struct gkyl_vlasov_poisson_field info; // data for field

  enum gkyl_vpfield_id field_id; // Type of field.

  struct gkyl_array *epsilon;  // Permittivity in Poisson equation.

  struct gkyl_job_pool *job_pool; // Job pool  
  // arrays for local and global charge density and potential.
  struct gkyl_array *rho_c, *rho_c_global;
  struct gkyl_array *phi, *phi_global;

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
                                      // for solving subset of Poisson solves with parallelization in z

  struct gkyl_fem_poisson *fem_poisson; // Poisson solver for - nabla . (epsilon * nabla phi) - kSq * phi = rho.

  bool has_ext_em; // flag to indicate there is external electromagnetic field
  bool ext_em_evolve; // flag to indicate external electromagnetic field is time dependent
  struct gkyl_array *ext_em; // external electromagnetic field
  struct gkyl_array *ext_em_host; // host copy for use in IO and projecting
  struct gkyl_array *tot_em; // total electromagnetic field
  gkyl_eval_on_nodes *ext_em_proj; // projector for external electromagnetic field 
  struct vp_eval_ext_em_ctx ext_em_ctx; // context for external electromagnetic field 

  struct gkyl_array *es_energy_fac; // Factor in calculation of ES energy diagnostic.
  struct gkyl_array_integrate *calc_es_energy;
  double *es_energy_red, *es_energy_red_global; // Memory for use in GPU reduction of ES energy.
  gkyl_dynvec integ_energy; // integrated energy components

  bool is_first_energy_write_call; // flag for energy dynvec written first time
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_poisson_app {
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
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis, vel-space basis

  struct gkyl_comm *comm;   // communicator object for conf-space arrays

  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *confBasis;
  } basis_on_dev;

  bool has_field; // has field
  struct vp_field *field; // pointer to field object

  // species data
  int num_species;
  struct vp_species *species; // data for each species
  
  struct gkyl_vlasov_poisson_stat stat; // statistics
};

/** gkyl_vlasov_poisson_app private API */

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct vp_species* vp_find_species(const gkyl_vlasov_poisson_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of species, -1 if not found
 */
int vp_find_species_idx(const gkyl_vlasov_poisson_app *app, const char *nm);

/** vp_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param sm Species moment object
 * @param nm Name string indicating moment type
 */
void vp_species_moment_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s,
  struct vp_species_moment *sm, const char *nm);

/**
 * Calculate moment, given distribution function @a fin.
 *
 * @param sm vp species moment
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void vp_species_moment_calc(const struct vp_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release species moment object.
 *
 * @param app Vlasov app object
 * @param sm Species moment object to release
 */
void vp_species_moment_release(const struct gkyl_vlasov_poisson_app *app,
  const struct vp_species_moment *sm);

/** vp_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void vp_species_lbo_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s,
  struct vp_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void vp_species_lbo_cross_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s,
  struct vp_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vp_species_lbo_moms(gkyl_vlasov_poisson_app *app,
  const struct vp_species *species,
  struct vp_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void vp_species_lbo_cross_moms(gkyl_vlasov_poisson_app *app,
  const struct vp_species *species,
  struct vp_lbo_collisions *lbo,
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
void vp_species_lbo_rhs(gkyl_vlasov_poisson_app *app,
  const struct vp_species *species,
  struct vp_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species LBO object.
 *
 * @param app Vlasov app object
 * @param sm Species LBO object to release
 */
void vp_species_lbo_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_lbo_collisions *lbo);

/** vp_species_bgk API */

/**
 * Initialize species BGK collisions object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param bgk Species BGK object
 */
void vp_species_bgk_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s,
  struct vp_bgk_collisions *bgk);

/**
 * Compute necessary moments for BGK collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void vp_species_bgk_moms(gkyl_vlasov_poisson_app *app,
  const struct vp_species *species,
  struct vp_bgk_collisions *bgk,
  const struct gkyl_array *fin);

/**
 * Compute RHS from BGK collisions
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 * @param rhs On output, the RHS from bgk
 */
void vp_species_bgk_rhs(gkyl_vlasov_poisson_app *app,
  const struct vp_species *species,
  struct vp_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species BGK object.
 *
 * @param app Vlasov app object
 * @param bgk Species BGK object to release
 */
void vp_species_bgk_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_bgk_collisions *bgk);

/** vp_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param bflux Species boundary flux object
 */
void vp_species_bflux_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s,
  struct vp_boundary_fluxes *bflux);

/**
 * Compute boundary flux from rhs
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param bflux Species boundary flux object
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void vp_species_bflux_rhs(gkyl_vlasov_poisson_app *app, const struct vp_species *species,
  struct vp_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species boundary flux object.
 *
 * @param app Vlasov app object
 * @param bflux Species boundary flux object to release
 */
void vp_species_bflux_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_boundary_fluxes *bflux);

/** vp_species_projection API */

/**
 * Initialize species projection object.
 *
 * @param app vlasov app object
 * @param s Species object 
 * @param inp Input struct for projection (contains functions pointers for type of projection)
 * @param proj Species projection object
 */
void vp_species_projection_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s, 
  struct gkyl_vlasov_poisson_projection inp, struct vp_proj *proj);

/**
 * Compute species projection
 *
 * @param app vlasov app object
 * @param species Species object
 * @param proj Species projection object
 * @param f Output distribution function from projection
 * @param tm Time for use in projection
 */
void vp_species_projection_calc(gkyl_vlasov_poisson_app *app, const struct vp_species *species, 
  struct vp_proj *proj, struct gkyl_array *f, double tm);

/**
 * Release species projection object.
 *
 * @param app vlasov app object
 * @param proj Species projection object to release
 */
void vp_species_projection_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_proj *proj);

/** vp_species_source API */

/**
 * Initialize species source object.
 *
 * @param app Vlasov app object
 * @param s Species object 
 * @param src Species source object
 */
void vp_species_source_init(struct gkyl_vlasov_poisson_app *app, struct vp_species *s, struct vp_source *src);

/**
 * Compute species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param src Pointer to source
 * @param tm Time for use in source
 */
void vp_species_source_calc(gkyl_vlasov_poisson_app *app, struct vp_species *species, 
  struct vp_source *src, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function
 */
void vp_species_source_rhs(gkyl_vlasov_poisson_app *app, const struct vp_species *species,
  struct vp_source *src, const struct gkyl_array *fin[], struct gkyl_array *rhs[]);

/**
 * Release species source object.
 *
 * @param app Vlasov app object
 * @param src Species source object to release
 */
void vp_species_source_release(const struct gkyl_vlasov_poisson_app *app, const struct vp_source *src);

/** vp_species API */

/**
 * Initialize species.
 *
 * @param vp Input VP data
 * @param app Vlasov app object
 * @param s On output, initialized species object
 */
void vp_species_init(struct gkyl_vp *vp, struct gkyl_vlasov_poisson_app *app, struct vp_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void vp_species_apply_ic(gkyl_vlasov_poisson_app *app, struct vp_species *species, double t0);

/**
 * Compute RHS from species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param em EM field
 * @param rhs On output, the RHS from the species object
 * @return Maximum stable time-step
 */
double vp_species_rhs(gkyl_vlasov_poisson_app *app, struct vp_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Apply BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param f Field to apply BCs
 */
void vp_species_apply_bc(gkyl_vlasov_poisson_app *app, const struct vp_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers
 */
void vp_species_coll_tm(gkyl_vlasov_poisson_app *app);

/**
 * Fill stat object in app with collisionless timers.
 *
 * @param app App object to update stat timers
 */
void vp_species_tm(gkyl_vlasov_poisson_app *app);

/**
 * Delete resources used in species.
 *
 * @param app Vlasov app object
 * @param species Species object to delete
 */
void vp_species_release(const gkyl_vlasov_poisson_app* app, const struct vp_species *s);

/** vp_field API */

/**
 * Create new field object
 *
 * @param vp Input VP data
 * @param app Vlasov app object
 * @return Newly created field
 */
struct vp_field* vp_field_new(struct gkyl_vp *vp, struct gkyl_vlasov_poisson_app *app);

/**
 * Compute field initial conditions.
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param t0 Time for use in ICs
 */
void vp_field_apply_ic(gkyl_vlasov_poisson_app *app, struct vp_field *field, double t0);

/**
 * Compute external electromagnetic fields
 *
 * @param app Vlasov app object
 * @param field Field object
 * @param tm Time for use in external electromagnetic fields computation
 */
void vp_field_calc_ext_em(gkyl_vlasov_poisson_app *app, struct vp_field *field, double tm);

/**
 * Accumulate charge density for Poisson solve.
 *
 * @param app Vlasov-Poisson app object
 * @param field Pointer to field
 * @param fin[] Input distribution function (num_species size)
 */
void vp_field_accumulate_rho_c(gkyl_vlasov_poisson_app *app, struct vp_field *field,
  const struct gkyl_array *fin[]);

/**
 * Compute RHS from field equations
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param em Input field
 * @param rhs On output, the RHS from the field solver
 * @return Maximum stable time-step
 */
void vp_field_rhs(gkyl_vlasov_poisson_app *app, struct vp_field *field);

/**
 * Compute field energy diagnostic
 *
 * @param app Vlasov app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 */
void vp_field_calc_energy(gkyl_vlasov_poisson_app *app, double tm, const struct vp_field *field);

/**
 * Release resources allocated by field
 *
 * @param app Vlasov app object
 * @param f Field object to release
 */
void vp_field_release(const gkyl_vlasov_poisson_app* app, struct vp_field *f);
