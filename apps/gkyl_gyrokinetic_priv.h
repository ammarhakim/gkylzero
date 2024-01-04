// Private header for use in Gyrokinetic app: do not include in user-facing
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
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_gyrokinetic_vars.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_dg_updater_diffusion_gyrokinetic.h>
#include <gkyl_dg_updater_lbo_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_gk_geometry_fromfile.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_tok_geo.h>
#include <gkyl_util.h>

#include <gkyl_gyrokinetic.h>

// Definitions of private structs and APIs attached to these objects
// for use in Gyrokinetic app.

// list of valid moment names
static const char *const valid_moment_names[] = {
  "M0",
  "M1",
  "M2",
  "M2par",
  "M2perp",
  "M3par",
  "M3perp",
  "ThreeMoments",
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
struct gk_species_moment {
  struct gk_geometry *gk_geom; // geometry struct for dividing moments by Jacobian
  struct gkyl_dg_bin_op_mem *mem_geo; // memory needed in dividing moments by Jacobian

  struct gkyl_dg_updater_moment *mcalc; // moment update

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

struct gk_rad_drag {  
  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  int collide_with_idx[GKYL_MAX_SPECIES]; // index of species we collide with

  // drag coefficients in vparallel and mu for each species being collided with
  struct gkyl_array *vnu[GKYL_MAX_SPECIES]; // vnu = 2/pi*|v|*nu(v)
  struct gkyl_array *vsqnu[GKYL_MAX_SPECIES]; // vsqnu = 1/2*(m/B)^(3/2)*sqrt(mu)*|v|^2*nu(v)
  struct gkyl_array *vnu_host[GKYL_MAX_SPECIES]; // host-side copy of vnu
  struct gkyl_array *vsqnu_host[GKYL_MAX_SPECIES]; // host-side copy of vsqnu
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars[GKYL_MAX_SPECIES]; 

  struct gk_species_moment moms[GKYL_MAX_SPECIES]; // moments needed in radiation update (need number density)

  gkyl_dg_updater_collisions *drag_slvr; // radiation solver
};

// forward declare species struct
struct gk_species;

struct gk_lbo_collisions {  
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *prim_moms, *nu_prim_moms; // LBO primitive moments
  struct gkyl_array *nu_sum_host, *prim_moms_host, *nu_prim_moms_host; // LBO primitive moments host-side for I/O
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

  struct gk_species_moment moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)

  struct gkyl_array *m0;
  struct gkyl_array *m2self; // m2self used for robustness of LBO
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  gkyl_dg_bin_op_mem *greene_factor_mem; // memory needed in computing Greene factor
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

struct gk_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // optional parameter to write out source distribution
  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *source_proj; // projector for source
};

// species data
struct gk_species {
  struct gkyl_gyrokinetic_species info; // data for species

  enum gkyl_gkmodel_id gkmodel_id;
  enum gkyl_gkfield_id gkfield_id;
  
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

  struct gkyl_array *alpha_surf; // array for surface phase space flux
  struct gkyl_array *sgn_alpha_surf; // array for the sign of the surface phase space flux at quadrature points
                                     // utilized for numerical flux function
                                     // F = alpha_surf/2 ( (f^+ + f^-) - sign_alpha_surf*(f^+ - f^-) )
  struct gkyl_array *const_sgn_alpha; // boolean array for if the surface phase space flux is single signed
                                      // if true, numerical flux function inside kernels simplifies to
                                      // F = alpha_surf*f^- (if sign_alpha_surf = 1), 
                                      // F = alpha_surf*f^+ (if sign_alpha_surf = -1)
  
  struct gkyl_array *phi; // array for electrostatic potential
  // organization of the different equation objects and the required data and solvers
  union {
    // EM GK model
    struct {
      struct gkyl_array *apar; // array for A_parallel
      struct gkyl_array *apardot; // array for d/dt A_parallel
    };
  };

  struct gkyl_dg_calc_gyrokinetic_vars *calc_gk_vars;

  struct gk_species_moment m0; // for computing charge density
  struct gk_species_moment integ_moms; // integrated moments
  struct gk_species_moment *moms; // diagnostic moments
  struct gkyl_array *L2_f; // L2 norm f^2
  double *red_L2_f; // for reduction of integrated L^2 norm on GPU
  double *red_integ_diag; // for reduction of integrated moments on GPU
  gkyl_dynvec integ_L2_f; // integrated L^2 norm reduced across grid
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_L2_write_call; // flag for integrated L^2 norm dynvec written first time
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  gkyl_dg_updater_gyrokinetic *slvr; // Gyrokinetic solver 
  struct gkyl_dg_eqn *eqn_gyrokinetic; // Gyrokinetic equation object
  
  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // gyrokinetic sheath boundary conditions
  struct gkyl_bc_sheath_gyrokinetic *bc_sheath_lo;
  struct gkyl_bc_sheath_gyrokinetic *bc_sheath_up;
  // Pointers to updaters that apply (non-sheath) BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  enum gkyl_source_id source_id; // type of source
  struct gk_source src; // applied source

  enum gkyl_collision_id collision_id; // type of collisions
  struct gk_lbo_collisions lbo; // collisions object

  enum gkyl_radiation_id radiation_id; // type of radiation
  struct gk_rad_drag rad; // radiation object

  double *omegaCfl_ptr;
};

// field data
struct gk_field {
  struct gkyl_gyrokinetic_field info; // data for field

  enum gkyl_gkfield_id gkfield_id;

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *rho_c, *rho_c_smooth; // arrays for charge density and smoothed charge density
  struct gkyl_array *phi_fem, *phi_smooth; // arrays for updates

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  // organization of the different equation objects and the required data and solvers
  union {
    // EM GK model
    struct {
      struct gkyl_array *apar_fem; // array for A_parallel
      struct gkyl_array *apardot_fem; // array for d/dt A_parallel
    };
  };

  struct gkyl_array *weight; 
  struct gkyl_array *es_energy_fac; 
  struct gkyl_array *epsilon; 
  struct gkyl_array *kSq; 

  struct gkyl_fem_parproj *fem_parproj; // FEM smoother for projecting DG functions onto continuous FEM basis
                                        // weight*phi_{fem} = phi_{dg} 
  struct gkyl_fem_poisson_perp *fem_poisson_perp; // perpendicular Poisson solve
                                                  // - nabla . (epsilon * nabla phi) - kSq * phi = rho

  struct gkyl_array_integrate *calc_em_energy;
  double *em_energy_red; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components

  bool is_first_energy_write_call; // flag for energy dynvec written first time

  bool has_phi_wall_lo; // flag to indicate there is biased wall potential on lower wall
  bool phi_wall_lo_evolve; // flag to indicate biased wall potential on lower wall is time dependent
  struct gkyl_array *phi_wall_lo; // biased wall potential on lower wall
  struct gkyl_array *phi_wall_lo_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *phi_wall_lo_proj; // projector for biased wall potential on lower wall 

  bool has_phi_wall_up; // flag to indicate there is biased wall potential on upper wall
  bool phi_wall_up_evolve; // flag to indicate biased wall potential on upper wall is time dependent
  struct gkyl_array *phi_wall_up; // biased wall potential on upper wall
  struct gkyl_array *phi_wall_up_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *phi_wall_up_proj; // projector for biased wall potential on upper wall 
};

// gyrokinetic object: used as opaque pointer in user code
struct gkyl_gyrokinetic_app {
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

  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis
  
  struct gkyl_comm *comm;   // communicator object for conf-space arrays

  bool has_mapc2p; // flag to indicate if we have mapc2p
  bool tokamak; // flag to indicate if it is a tokamak geometry
  bool geo_fromfile; // flag to indicate if we should just read the geometry
  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);
  void *bmag_ctx; // context for bmag function
  // pointer to bmag function
  void (*bmag_func)(double t, const double *xc, double *xp, void *ctx);

  void *tok_rz_ctx; // context with RZ data such as efit file for a tokamak
  void *tok_comp_ctx; // context for tokamak geometry with computational domain info

  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *confBasis;
  } basis_on_dev;

  struct gk_geometry *gk_geom;

  struct gk_field *field; // pointer to field object

  // species data
  int num_species;
  struct gk_species *species; // data for each species
  
  struct gkyl_gyrokinetic_stat stat; // statistics
};

/** gkyl_gyrokinetic_app private API */

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct gk_species* gk_find_species(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of species, -1 if not found
 */
int gk_find_species_idx(const gkyl_gyrokinetic_app *app, const char *nm);


/** gk_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param sm Species moment object
 * @param nm Name string indicating moment type
 */
void gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, const char *nm);

/**
 * Calculate moment, given distribution function @a fin.
 * 
 * @param sm Species moment object
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void gk_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release species moment object.
 *
 * @param app gyrokinetic app object
 * @param sm Species moment object to release
 */
void gk_species_moment_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species_moment *sm);

/** gk_species_radiation API */

/**
 * Initialize species radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param rad Species radiation drag object
 */
void gk_species_radiation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_rad_drag *rad);

/**
 * Compute RHS from radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void gk_species_radiation_rhs(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_rad_drag *rad,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param rad Species radiation drag object to release
 */
void gk_species_radiation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_rad_drag *rad);

/** gk_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void gk_species_lbo_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void gk_species_lbo_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void gk_species_lbo_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void gk_species_lbo_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute RHS from LBO collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void gk_species_lbo_rhs(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species LBO object.
 *
 * @param app gyrokinetic app object
 * @param lbo Species LBO object to release
 */
void gk_species_lbo_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lbo_collisions *lbo);

/** gk_species_source API */

/**
 * Initialize species source object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param src Species source object
 */
void gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_source *src);

/**
 * Compute species applied source term
 *
 * @param app gyrokinetic app object
 * @param species Species object
 * @param tm Time for use in source
 */
void gk_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_species *species, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function
 */
void gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species source object.
 *
 * @param app gyrokinetic app object
 * @param src Species source object to release
 */
void gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_species API */

/**
 * Initialize species.
 *
 * @param gk Input gk data
 * @param app gyrokinetic app object
 * @param s On output, initialized species object
 */
void gk_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app gyrokinetic app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute RHS from species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param rhs On output, the RHS from the species object
 * @param fluidin Input fluid array for potential fluid force (size: num_fluid_species)
 * @return Maximum stable time-step
 */
double gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Apply BCs to species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param f Field to apply BCs
 */
void gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f);

/**
 * Compute L2 norm (f^2) of the distribution function diagnostic
 *
 * @param app gyrokinetic app object
 * @param tm Time at which diagnostic is computed
 * @param species Pointer to species
 */
void gk_species_calc_L2(gkyl_gyrokinetic_app *app, double tm, const struct gk_species *species);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers
 */
void gk_species_coll_tm(gkyl_gyrokinetic_app *app);

/**
 * Fill stat object in app with collisionless timers.
 *
 * @param app App object to update stat timers
 */
void gk_species_tm(gkyl_gyrokinetic_app *app);

/**
 * Delete resources used in species.
 *
 * @param app gyrokinetic app object
 * @param species Species object to delete
 */
void gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *s);

/** gk_field API */

/**
 * Create new field object
 *
 * @param gk Input gk data
 * @param app gyrokinetic app object
 * @return Newly created field
 */
struct gk_field* gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app);

/**
 * Compute biased wall potentials 
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 * @param tm Time to compute biased wall potentials at
 */
void gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm);

/**
 * Accumulate charge density for Poisson solve
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 * @param fin[] Input distribution function (num_species size)
 */
void gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *fin[]);

/**
 * Compute EM field 
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 * @param em Output field
 */
void gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Compute field energy diagnostic
 *
 * @param app gyrokinetic app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 */
void gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field);

/**
 * Release resources allocated by field
 *
 * @param app gyrokinetic app object
 * @param f Field object to release
 */
void gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f);
