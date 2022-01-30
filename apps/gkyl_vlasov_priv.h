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
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_lbo_updater.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_lbo_mom_bcorr.h>
#include <gkyl_mom_bcorr.h>
#include <gkyl_mom_calc.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_mom.h>


// Definitions of private structs and APIs attached to these objects
// for use in Vlasov app.

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

struct vm_lbo_collisions {
  struct gkyl_array *cM, *cE; // LBO boundary corrections
  struct gkyl_mom_type *cM_mom, *cE_mom; // moments needed in update
  struct gkyl_mom_bcorr *cM_bcorr, *cE_bcorr; // moments needed in update
  struct gkyl_array *nu_sum, *u_drift, *vth_sq, *nu_u, *nu_vthsq; // LBO primitive moments
  struct gkyl_prim_lbo *coll_prim; // primitive moments

  struct vm_species_moment m0, m1i, m2; // moments needed in LBO
  
  gkyl_prim_lbo_calc *coll_pcalc; // primitive moment calculator
  gkyl_dg_lbo_updater *coll_slvr; // collision solver
};

struct vm_bgk_collisions {
  struct gkyl_array *u, *vthsq; // BGK primitive moments
  // BGK Collisions should probably own a project on Maxwellian updater
  // so it can compute its contribution to the RHS
  // struct proj_maxwellian;
};

// species data
struct vm_species {
  struct gkyl_vlasov_species info; // data for species
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct vm_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *f_host; // host copy for use IO and initialization

  struct vm_species_moment m1i; // for computing currents
  struct vm_species_moment *moms; // diagnostic moments

  struct gkyl_dg_eqn *eqn; // Vlasov equation
  gkyl_hyper_dg *slvr; // solver 

  enum gkyl_collision_id collision_id; // type of collisions
  struct vm_lbo_collisions lbo;
  double *omegaCfl_ptr;
};

// field data
struct vm_field {
  struct gkyl_vlasov_field info; // data for field

  struct gkyl_job_pool *job_pool; // Job pool  
  struct gkyl_array *em, *em1, *emnew; // arrays for updates
  struct gkyl_array *qmem; // array for q/m*(E,B)
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *em_host;  // host copy for use IO and initialization

  gkyl_hyper_dg *slvr; // Maxwell solver

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
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  struct vm_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  bool has_field; // has field
  struct vm_field *field; // pointer to field object

  // species data
  int num_species;
  struct vm_species *species; // data for each species

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
 */
void vm_species_lbo_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_lbo_collisions *lbo);

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
double vm_species_lbo_rhs(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lbo_collisions *lbo, const struct gkyl_array *fin, struct gkyl_array *rhs);

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
 * Compute RHS from species distribution function
 *
 * @param app Vlasov app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param qmem EM field scaled by q/m
 * @param rhs On output, the RHS from the species object
 * @return Maximum stable time-step
 */
double vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *qmem, struct gkyl_array *rhs);

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
 * @param f Field to apply BCs
 */
void vm_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_species *species,
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
void vm_field_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_field *field, int dir, struct gkyl_array *f);

/**
 * Apply copy BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param dir Direction to apply BCs
 * @param f Field to apply BCs
 */
void vm_field_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_field *field, int dir, struct gkyl_array *f);

/**
 * Apply BCs to field
 *
 * @param app Vlasov app object
 * @param field Pointer to field
 * @param f Field to apply BCs
 */
void vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f);

/**
 * Release resources allocated by field
 *
 * @param app Vlasov app object
 * @param f Field object to release
 */
void vm_field_release(const gkyl_vlasov_app* app, struct vm_field *f);
