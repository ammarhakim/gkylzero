// Private header for use in Vlasov app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_updater_lbo_vlasov.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_poisson.h>
#include <gkyl_dg_vlasov_sr.h>
#include <gkyl_dynvec.h>
#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_with_fluid.h>
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

// ranges for use in BCs
struct vm_skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// data for moments
struct vm_species_moment {
  bool use_gpu; // should we use GPU (if present)
  gkyl_mom_calc *mcalc; // moment update
  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

// forward declare species struct
struct vm_species;

struct vm_lbo_collisions {
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_type *bcorr_type; // LBO boundary corrections moment type
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *u_drift, *vth_sq, *nu_u, *nu_vthsq; // LBO primitive moments

  double betaGreenep1; // value of Greene's factor beta + 1
  double other_m[GKYL_MAX_SPECIES]; // masses of species being collided with
  struct gkyl_array *other_u_drift[GKYL_MAX_SPECIES], *other_vth_sq[GKYL_MAX_SPECIES]; // self-primitive moments of species being collided with
  struct gkyl_array *cross_u_drift[GKYL_MAX_SPECIES], *cross_vth_sq[GKYL_MAX_SPECIES]; // LBO cross-primitive moments
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // LBO cross-species collision frequencies
  struct gkyl_array *other_nu[GKYL_MAX_SPECIES];
  struct gkyl_array *cross_nu_u, *cross_nu_vthsq; // weak multiplication of collision frequency and primitive moments
  
  struct gkyl_array *self_nu, *self_nu_u, *self_nu_vthsq; // LBO self-primitive moments
  struct gkyl_prim_lbo_type *coll_prim; // LBO primitive moments type

  struct vm_species_moment moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)
  struct gkyl_array *m0;
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct vm_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_lbo_vlasov *coll_slvr; // collision solver
};

struct vm_bgk_collisions {
  struct gkyl_array *u, *vthsq; // BGK primitive moments
  // BGK Collisions should probably own a project on Maxwellian updater
  // so it can compute its contribution to the RHS
  // struct proj_maxwellian;
};


// context for use in computing applied acceleration
struct vm_eval_accel_ctx { evalf_t accel_func; void *accel_ctx; };

// species data
struct vm_species {
  struct gkyl_vlasov_species info; // data for species
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct vm_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *f_host; // host copy for use IO and initialization

  struct vm_species_moment m1i; // for computing currents
  struct vm_species_moment *moms; // diagnostic moments
  struct vm_species_moment integ_moms; // integrated moments

  double *red_integ_diag; // for reduction on GPU
  gkyl_dynvec integ_diag; // integrated moments reduced across grid

  bool is_first_integ_write_call; // flag for int-moments dynvec written first time

  enum gkyl_field_id field_id; // type of Vlasov equation (based on type of field solve)
  struct gkyl_array *qmem; // array for q/m*(E,B)
  struct gkyl_array *fac_phi; // array for potential (electrostatic or gravitational)
  struct gkyl_array *vecA; // array for vector potential
  struct gkyl_array *p_over_gamma; // array for p/gamma (velocity) in special relativistic equation
  struct gkyl_array *p_over_gamma_host; // host copy for use in projecting before copying over to GPU

  gkyl_dg_updater_vlasov *slvr; // Vlasov solver 
  struct gkyl_dg_eqn *eqn_vlasov; // Vlasov equation object

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Note: we need to store pointers to the struct as these may
  // actually be on the GPUs. Seems ugly, but I am not sure how else
  // to ensure the function and context lives on the GPU
  struct gkyl_array_copy_func *wall_bc_func[3]; // for wall BCs
  struct gkyl_array_copy_func *absorb_bc_func[3]; // for absorbing BCs

  bool has_accel; // flag to indicate there is applied acceleration
  struct gkyl_array *accel; // applied acceleration
  struct gkyl_array *accel_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *accel_proj; // projector for acceleration
  struct vm_eval_accel_ctx accel_ctx; // context for applied acceleration

  bool has_source; // flag to indicate there is applied source
  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *source_proj; // projector for source

  bool has_mirror_force; // flag to indicate Vlasov includes mirror force from external magnetic field
  struct gkyl_array *gradB; // gradient of magnetic field
  struct gkyl_array *magB; // magnitude of magnetic field (J = 1/B)
  struct gkyl_array *n; // array storing density (no Jacobian)
  struct gkyl_array *Tperp; // array storing J*Tperp (J*p_perp/n)
  struct gkyl_array *mirror_force; // array storing full mirror force (J*T_perp*gradB)
  struct gkyl_array *m1i_no_J; // current density without Jacobian (for coupling to EM fields)

  // host copy for use in IO and projecting
  struct gkyl_array *gradB_host;
  struct gkyl_array *magB_host;
  struct gkyl_array *n_host;
  struct gkyl_array *Tperp_host;
  struct gkyl_array *mirror_force_host;
  struct gkyl_array *m1i_no_J_host;

  enum gkyl_collision_id collision_id; // type of collisions
  bool collides_with_fluid; // boolean for if kinetic species collides with a fluid speceis
  int fluid_index; // index of the fluid species being collided with
                   // index corresponds to location in fluid_species array (size num_fluid_species)
  struct vm_lbo_collisions lbo; // collisions object

  double *omegaCfl_ptr;
};

// field data
struct vm_field {
  struct gkyl_vlasov_field info; // data for field

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *em, *em1, *emnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *em_host;  // host copy for use IO and initialization

  gkyl_hyper_dg *slvr; // Maxwell solver

  struct gkyl_array *em_energy; // EM energy components in each cell
  double *em_energy_red; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components

  bool is_first_energy_write_call; // flag for energy dynvec written first time

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_field_bc_type lower_bc[3], upper_bc[3];
  // Note: we need to store pointers to the struct as these may
  // actually be on the GPUs. Seems ugly, but I am not sure how else
  // to ensure the function and context lives on the GPU
  struct gkyl_array_copy_func *wall_bc_func[3]; // for wall BCs

  double* omegaCfl_ptr;
};

// context for use in computing applied advection
struct vm_eval_advect_ctx { evalf_t advect_func; void *advect_ctx; };

// fluid species data
struct vm_fluid_species {
  struct gkyl_vlasov_fluid_species info; // data for fluid

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *fluid, *fluid1, *fluidnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *fluid_host;  // host copy for use IO and initialization

  struct gkyl_array *u; // array for advection flow

  struct gkyl_dg_eqn *eqn; // Fluid equation  
  gkyl_hyper_dg *slvr; // Fluid equation solver

  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_fluid_species_bc_type lower_bc[3], upper_bc[3];
  // Note: we need to store pointers to the struct as these may
  // actually be on the GPUs. Seems ugly, but I am not sure how else
  // to ensure the function and context lives on the GPU
  struct gkyl_array_copy_func *absorb_bc_func[3]; // for absorbing BCs

  // specified advection
  bool has_advect; // flag to indicate there is applied advection
  struct gkyl_array *advect; // applied advection
  struct gkyl_array *advect_host; // host copy for use in IO and projecting
  gkyl_proj_on_basis *advect_proj; // projector for advection
  struct vm_eval_advect_ctx advect_ctx; // context for applied advection

  // advection with another species present
  bool advects_with_species; // flag to indicate we are advecting with another species
  struct vm_species *advection_species; // pointer to species we advect with
  struct gkyl_array *other_advect; // pointer to that species drift velocity

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
  struct gkyl_basis basis, confBasis, velBasis; // phase-space, conf-space basis, vel-space basis

  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *confBasis;
  } basis_on_dev;

  struct vm_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

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

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Compute out = c1*arr1 + c2*arr2
static inline struct gkyl_array*
array_combine(struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range rng)
{
  return gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct vm_skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

/** gkyl_vlasov_app private API */

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.o
 */
struct vm_species* vm_find_species(const gkyl_vlasov_app *app, const char *nm);

/**
 * Find fluid species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of fluid species
 * @return Pointer to fluid species with given name. NULL if not found.o
 */
struct vm_fluid_species* vm_find_fluid_species(const gkyl_vlasov_app *app, const char *nm);

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
  struct vm_lbo_collisions *lbo, bool collides_with_fluid);

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
 * @param collides_with_fluid Boolean for if kinetic species collides with a fluid species
 * @param fluidin Input fluid array (size: num_fluid_species)
 */
void vm_species_lbo_moms(gkyl_vlasov_app *app,
  const struct vm_species *species,
  struct vm_lbo_collisions *lbo,
  const struct gkyl_array *fin,
  bool collides_with_fluid, const struct gkyl_array *fluidin[]);

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
  const struct gkyl_array *fin,
  bool collides_with_fluid, const struct gkyl_array *fluidin[]);

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
 * Compute species applied source term
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param tm Time for use in source
 */
void vm_species_calc_source(gkyl_vlasov_app *app, struct vm_species *species, double tm);

/**
 * Compute gradient and magnitude of magnetic field
 *
 * @param app Vlasov app object
 * @param species Species object
 * @param tm Time for use in source
 */
void vm_species_calc_magB_gradB(gkyl_vlasov_app *app, struct vm_species *species, double tm);

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
  struct gkyl_array *rhs,
  const struct gkyl_array *fluidin[]);

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
 * Apply copy BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param dir Direction to apply BCs
 * @param edge Edge to apply BCs
 * @param f Field to apply BCs
 */
void
vm_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

/**
 * Apply wall BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param dir Direction to apply BCs
 * @param edge Edge to apply BCs
 * @param f Field to apply BCs
 */
void
vm_species_apply_wall_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

/**
 * Apply absorbing BCs to species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param dir Direction to apply BCs
 * @param edge Edge to apply BCs
 * @param f Field to apply BCs
 */
void
vm_species_apply_absorb_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

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
 * Apply copy BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param dir Direction to apply BCs
 * @param f Field to apply BCs
 */
void vm_field_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

/**
 * Apply PEC wall BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param dir Direction to apply BCs
 * @param f Field to apply BCs
 */
void vm_field_apply_pec_bc(gkyl_vlasov_app *app, const struct vm_field *field,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

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
 * Compute RHS from fluid species equations
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param fluid Input fluid species
 * @param rhs On output, the RHS from the fluid species solver
 * @return Maximum stable time-step
 */
double vm_fluid_species_rhs(gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, const struct gkyl_array *fluid, struct gkyl_array *rhs);

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
 * Apply copy BCs to fluid species
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param dir Direction to apply BCs
 * @param f Fluid Species to apply BCs
 */
void vm_fluid_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

/**
 * Apply absorbing BCs to fluid species
 *
 * @param app Vlasov app object
 * @param fluid_species Pointer to fluid species
 * @param dir Direction to apply BCs
 * @param f Fluid Species to apply BCs
 */
void vm_fluid_species_apply_absorb_bc(gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  int dir, enum vm_domain_edge edge, struct gkyl_array *f);

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
