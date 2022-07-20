#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fv_proj.h>
#include <gkyl_kep_scheme.h>
#include <gkyl_moment.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_ten_moment.h>

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Species data
struct moment_species {
  int ndim;
  char name[128]; // species name
  double charge, mass;
  double k0; // closure parameter (default is 0.0, used by 10 moment)

  int evolve; // evolve species? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
    
  struct gkyl_array *fdup, *f[4]; // arrays for updates
  struct gkyl_array *cflrate; // cflrate in each cell used by KEP scheme
  struct gkyl_array *app_accel; // array for applied acceleration/forces
  // pointer to projection operator for applied acceleration/forces function
  gkyl_fv_proj *proj_app_accel;
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  enum gkyl_eqn_type eqn_type; // type ID of equation
  enum gkyl_braginskii_type type_brag; // which Braginskii equations
  int num_equations; // number of equations in species

  const struct gkyl_wv_eqn *equation; // equation object for initializing solvers

  enum gkyl_moment_fluid_scheme mom_fluid_scheme; // particular scheme to be employed (KEP vs. wave propagation)
  gkyl_kep_scheme *kep_slvr; // kep slvr
  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary condition type
  enum gkyl_species_bc_type lower_bct[3], upper_bct[3];
  // boundary condition solvers on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];
};

// Field data
struct moment_field {
  int ndim;
  double epsilon0, mu0;

  int evolve; // evolve species? 1-yes, 0-no
    
  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);    
    
  struct gkyl_array *fdup, *f[4]; // arrays for updates
  struct gkyl_array *app_current, *ext_em; // arrays for applied currents/external fields
  // pointer to projection operator for applied current function
  gkyl_fv_proj *proj_app_current;
  // pointer to projection operator for external fields
  gkyl_fv_proj *proj_ext_em;
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary condition type
  enum gkyl_field_bc_type lower_bct[3], upper_bct[3];
  // boundary conditions on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];    
};

// Source data
struct moment_coupling {
  gkyl_moment_em_coupling *slvr; // source solver function
};

struct moment_non_ideal {
  struct gkyl_rect_grid non_ideal_grid; // grid for braginskii variables (braginskii variables located at cell nodes)
  struct gkyl_range non_ideal_local, non_ideal_local_ext; // local, local-ext ranges for braginskii variables (loop over nodes)

  gkyl_moment_braginskii *brag_slvr; // Braginskii solver (if present)
  struct gkyl_array *non_ideal_cflrate[GKYL_MAX_SPECIES]; // array for stable time-step from non-ideal terms
  struct gkyl_array *non_ideal_vars[GKYL_MAX_SPECIES]; // array for braginskii variables
};

// Moment app object: used as opaque pointer in user code
struct gkyl_moment_app {
  char name[128]; // name of app
  int ndim; // space dimensions
  double tcurr; // current time
  double cfl; // CFL number

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int is_dir_skipped[3]; // flags to tell if update in direction are skipped

  enum gkyl_moment_fluid_scheme fluid_scheme; // scheme to update fluid equations
    
  struct gkyl_rect_grid grid; // grid
  struct gkyl_range local, local_ext; // local, local-ext ranges

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom *geom; // geometry needed for species and field solvers

  struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  int has_field; // flag to indicate if we have a field
  struct moment_field field; // field data
    
  // species data
  int num_species;
  struct moment_species *species; // species data

  int update_sources; // flag to indicate if sources are to be updated
  struct moment_coupling sources; // sources

  bool has_brag;
  int has_non_ideal; // flag to indicate if non-ideal terms are present
  double coll_fac; // multiplicative collisionality factor for Braginskii
  struct moment_non_ideal non_ideal; // non-ideal terms
    
  struct gkyl_moment_stat stat; // statistics
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
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


// function for copy BC
static void
bc_copy(double t, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  for (int c=0; c<nc; ++c) ghost[c] = skin[c];
}

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
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

// apply periodic BCs
static void
moment_apply_periodic_bc(const gkyl_moment_app *app, struct gkyl_array *bc_buffer,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(bc_buffer->data, f, app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, bc_buffer->data, app->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(bc_buffer->data, f, app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
}

static void
moment_apply_periodic_corner_sync(const gkyl_moment_app *app, struct gkyl_array *f)
{
  long idx_src, idx_dest;
  double *out;
  const double *inp;
  if (app->ndim == 2) {
    // LL skin cell -> UU ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.lower[0],app->local.lower[1]});
    idx_dest = gkyl_range_idx(&app->local, (int[]) {app->local.upper[0]+1,app->local.upper[1]+1});
    
    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);

    // LU skin cell -> UL ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.lower[0],app->local.upper[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.upper[0]+1,app->local.lower[1]-1});
    
    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);

    // UL skin cell -> LU ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.upper[0],app->local.lower[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.lower[0]-1,app->local.upper[1]+1});
    
    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);

    // UU skin cell -> LL ghost cell
    idx_src = gkyl_range_idx(&app->local, (int[]) {app->local.upper[0],app->local.upper[1]});
    idx_dest = gkyl_range_idx(&app->local_ext, (int[]) {app->local.lower[0]-1,app->local.lower[1]-1});
    
    out = (double*) gkyl_array_fetch(f, idx_dest);
    inp = (const double*) gkyl_array_cfetch(f, idx_src);
    gkyl_copy_double_arr(f->ncomp, inp, out);
  }
}

// apply wedge BCs
static void
moment_apply_wedge_bc(const gkyl_moment_app *app, double tcurr,
  const struct gkyl_range *update_rng, struct gkyl_array *bc_buffer,
  int dir, const struct gkyl_wv_apply_bc *lo, const struct gkyl_wv_apply_bc *up,
  struct gkyl_array *f)
{
  gkyl_wv_apply_bc_to_buff(lo, tcurr, update_rng, f, bc_buffer->data);
  gkyl_array_copy_from_buffer(f, bc_buffer->data, app->skin_ghost.upper_ghost[dir]);

  gkyl_wv_apply_bc_to_buff(up, tcurr, update_rng, f, bc_buffer->data);
  gkyl_array_copy_from_buffer(f, bc_buffer->data, app->skin_ghost.lower_ghost[dir]);
}

/** moment_species functions */

// initialize species
static void
moment_species_init(const struct gkyl_moment *mom, const struct gkyl_moment_species *mom_sp,
  struct gkyl_moment_app *app, struct moment_species *sp)
{
  sp->ndim = mom->ndim;
  strcpy(sp->name, mom_sp->name);
  sp->charge = mom_sp->charge;
  sp->mass = mom_sp->mass;
  sp->ctx = mom_sp->ctx;
  sp->init = mom_sp->init;

  sp->eqn_type = mom_sp->equation->type;
  sp->type_brag = mom_sp->type_brag;
  sp->num_equations = mom_sp->equation->num_equations;
  sp->equation = gkyl_wv_eqn_acquire(mom_sp->equation);
  // closure parameter, used by 10 moment
  sp->k0 = mom_sp->equation->type == GKYL_EQN_TEN_MOMENT ? gkyl_wv_ten_moment_k0(mom_sp->equation) : 0.0;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  int ndim = mom->ndim;

  // create kep updater
  sp->kep_slvr = gkyl_kep_scheme_new ( (struct gkyl_kep_scheme_inp) {
      .grid = &app->grid,
      .equation = sp->equation,
      .num_up_dirs = ndim,
      .update_dirs = {0, 1, 2},
      .cfl = app->cfl,
    }
  );
  // create wave propagation updaters for each directional update
  for (int d=0; d<ndim; ++d)
    sp->slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = sp->equation,
        .limiter = limiter,
        .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
        .update_dirs = { d },
        .cfl = app->cfl,
        .geom = app->geom,
      }
    );

  // which scheme are we using (KEP vs. wave propagation)
  sp->mom_fluid_scheme = app->fluid_scheme;

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    sp->lower_bc[i] = 0;
    sp->upper_bc[i] = 0;
  }

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_species_bc_type *bc;
      if (dir == 0)
        bc = mom_sp->bcx;
      else if (dir == 1)
        bc = mom_sp->bcy;
      else
        bc = mom_sp->bcz;

      sp->lower_bct[dir] = bc[0];
      sp->upper_bct[dir] = bc[1];

      // lower BCs
      switch (bc[0]) {
        case GKYL_SPECIES_WALL:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE: // wedge also uses bc_copy
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
          break;

        default:
          // can't happen
          break;
      }
      
      // upper BCs
      switch (bc[1]) {
        case GKYL_SPECIES_WALL:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;
            
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
          break;
      }
    }
  }

  int meqn = sp->num_equations;
  sp->fdup = mkarr(meqn, app->local_ext.volume);
  // allocate arrays
  for (int d=0; d<ndim+1; ++d)
    sp->f[d] = mkarr(meqn, app->local_ext.volume);

  sp->cflrate = mkarr(1, app->local_ext.volume);

  // allocate array for applied acceleration/forces for each species
  sp->app_accel = mkarr(3, app->local_ext.volume);
  sp->proj_app_accel = 0;
  if (mom_sp->app_accel_func)
    sp->proj_app_accel = gkyl_fv_proj_new(&app->grid, 2, 3, mom_sp->app_accel_func, sp->ctx);
  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  sp->bc_buffer = mkarr(meqn, buff_sz);
}

// apply BCs to species
static void
moment_species_apply_bc(const gkyl_moment_app *app, double tcurr,
  const struct moment_species *sp, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    moment_apply_periodic_bc(app, sp->bc_buffer, app->periodic_dirs[d], f);
    is_non_periodic[app->periodic_dirs[d]] = 0;
  }
  moment_apply_periodic_corner_sync(app, f);  
  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (sp->lower_bct[d] != GKYL_SPECIES_WEDGE)
        gkyl_wv_apply_bc_advance(sp->lower_bc[d], tcurr, &app->local, f);
      if (sp->upper_bct[d] != GKYL_SPECIES_WEDGE)
        gkyl_wv_apply_bc_advance(sp->upper_bc[d], tcurr, &app->local, f);

      // wedge BCs for upper/lower must be handled in one shot
      if (sp->lower_bct[d] == GKYL_SPECIES_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          sp->bc_buffer, d, sp->lower_bc[d], sp->upper_bc[d], f);
    }
}

// maximum stable time-step
static double
moment_species_max_dt(const gkyl_moment_app *app, const struct moment_species *sp)
{
  double max_dt = DBL_MAX;
  for (int d=0; d<app->ndim; ++d)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(sp->slvr[d], &app->local, sp->f[0]));
  return max_dt;
}

// update solution: initial solution is in sp->f[0] and updated
// solution in sp->f[ndim]
static struct gkyl_update_status
moment_species_update(const gkyl_moment_app *app,
  const struct moment_species *sp, double tcurr, double dt)
{
  int ndim = sp->ndim;
  double dt_suggested = DBL_MAX;
  struct gkyl_wave_prop_status stat;

  for (int d=0; d<ndim; ++d) {
    stat = gkyl_wave_prop_advance(sp->slvr[d], tcurr, dt, &app->local, sp->f[d], sp->f[d+1]);

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };
    
    dt_suggested = fmin(dt_suggested, stat.dt_suggested);
    moment_species_apply_bc(app, tcurr, sp, sp->f[d+1]);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

// free species
static void
moment_species_release(const struct moment_species *sp)
{
  gkyl_wv_eqn_release(sp->equation);
  gkyl_kep_scheme_release(sp->kep_slvr);
  for (int d=0; d<sp->ndim; ++d)
    gkyl_wave_prop_release(sp->slvr[d]);

  for (int d=0; d<sp->ndim; ++d) {
    if (sp->lower_bc[d])
      gkyl_wv_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d])    
      gkyl_wv_apply_bc_release(sp->upper_bc[d]);
  }

  gkyl_array_release(sp->fdup);
  for (int d=0; d<sp->ndim+1; ++d)
    gkyl_array_release(sp->f[d]);

  gkyl_array_release(sp->cflrate);

  gkyl_array_release(sp->app_accel);
  if (sp->proj_app_accel)
    gkyl_fv_proj_release(sp->proj_app_accel);

  gkyl_array_release(sp->bc_buffer);
}

/** moment_field functions */

// initialize field
static void
moment_field_init(const struct gkyl_moment *mom, const struct gkyl_moment_field *mom_fld,
  struct gkyl_moment_app *app, struct moment_field *fld)
{
  fld->ndim = mom->ndim;
  double epsilon0 = fld->epsilon0 = mom_fld->epsilon0;
  double mu0 = fld->mu0 = mom_fld->mu0;

  fld->ctx = mom_fld->ctx;
  fld->init = mom_fld->init;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_fld->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_fld->limiter;

  double c = 1/sqrt(epsilon0*mu0);
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c,
    mom_fld->elc_error_speed_fact, mom_fld->mag_error_speed_fact);

  int ndim = mom->ndim;
  // create updaters for each directional update
  for (int d=0; d<ndim; ++d)
    fld->slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = maxwell,
        .limiter = limiter,
        .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
        .update_dirs = { d },
        .cfl = app->cfl,
        .geom = app->geom,
      }
    );

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  for (int i=0; i<3; ++i) {
    fld->lower_bc[i] = 0;
    fld->upper_bc[i] = 0;
  }  

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0)
        bc = mom_fld->bcx;
      else if (dir == 1)
        bc = mom_fld->bcy;
      else
        bc = mom_fld->bcz;

      fld->lower_bct[dir] = bc[0];
      fld->upper_bct[dir] = bc[1];

      switch (bc[0]) {
        case GKYL_FIELD_PEC_WALL:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
          &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost, maxwell->wall_bc_func, 0);
          break;

        case GKYL_FIELD_COPY:
        case GKYL_FIELD_WEDGE:
          fld->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
          break;
      }

      switch (bc[1]) {
        case GKYL_FIELD_PEC_WALL:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost, maxwell->wall_bc_func, 0);
          break;

        case GKYL_FIELD_COPY:
        case GKYL_FIELD_WEDGE:
          fld->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, maxwell, app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
      }
    }
  }

  // allocate arrays
  fld->fdup = mkarr(8, app->local_ext.volume);
  for (int d=0; d<ndim+1; ++d)
    fld->f[d] = mkarr(8, app->local_ext.volume);

  // allocate arrays for applied current/external fields
  fld->app_current = mkarr(3, app->local_ext.volume);
  fld->proj_app_current = 0;
  if (mom_fld->app_current_func)
    fld->proj_app_current = gkyl_fv_proj_new(&app->grid, 2, 3, mom_fld->app_current_func, fld->ctx);
  fld->ext_em = mkarr(6, app->local_ext.volume);
  fld->proj_ext_em = 0;
  if (mom_fld->ext_em_func)
    fld->proj_ext_em = gkyl_fv_proj_new(&app->grid, 2, 6, mom_fld->ext_em_func, fld->ctx);

  // allocate buffer for applying BCs (used for periodic BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->ndim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  fld->bc_buffer = mkarr(8, buff_sz);

  gkyl_wv_eqn_release(maxwell);
}

// apply BCs to EM field
static void
moment_field_apply_bc(const gkyl_moment_app *app, double tcurr,
  const struct moment_field *field, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    moment_apply_periodic_bc(app, field->bc_buffer, app->periodic_dirs[d], f);
    is_non_periodic[app->periodic_dirs[d]] = 0;
  }

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {
      // handle non-wedge BCs
      if (field->lower_bct[d] != GKYL_FIELD_WEDGE)
        gkyl_wv_apply_bc_advance(field->lower_bc[d], tcurr, &app->local, f);
      if (field->upper_bct[d] != GKYL_FIELD_WEDGE)      
        gkyl_wv_apply_bc_advance(field->upper_bc[d], tcurr, &app->local, f);

      // wedge BCs for upper/lower must be handled in one shot
      if (field->lower_bct[d] == GKYL_FIELD_WEDGE)
        moment_apply_wedge_bc(app, tcurr, &app->local,
          field->bc_buffer, d, field->lower_bc[d], field->upper_bc[d], f);
    }  
}

static double
moment_field_max_dt(const gkyl_moment_app *app, const struct moment_field *fld)
{
  double max_dt = DBL_MAX;
  for (int d=0; d<app->ndim; ++d)
    max_dt = fmin(max_dt, gkyl_wave_prop_max_dt(fld->slvr[d], &app->local, fld->f[0]));
  return max_dt;
}

// update solution: initial solution is in fld->f[0] and updated
// solution in fld->f[ndim]
static struct gkyl_update_status
moment_field_update(const gkyl_moment_app *app,
  const struct moment_field *fld, double tcurr, double dt)
{
  int ndim = fld->ndim;
  struct gkyl_wave_prop_status stat = { 1, DBL_MAX };

  for (int d=0; d<ndim; ++d) {
    // update solution
    stat = gkyl_wave_prop_advance(fld->slvr[d], tcurr, dt, &app->local, fld->f[d], fld->f[d+1]);

    if (!stat.success)
      return (struct gkyl_update_status) {
        .success = false,
        .dt_suggested = stat.dt_suggested
      };
    // apply BC
    moment_field_apply_bc(app, tcurr, fld, fld->f[d+1]);
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = stat.dt_suggested
  };
}

// free field
static void
moment_field_release(const struct moment_field *fld)
{
  for (int d=0; d<fld->ndim; ++d)
    gkyl_wave_prop_release(fld->slvr[d]);

  for (int d=0; d<fld->ndim; ++d) {
    if (fld->lower_bc[d])
      gkyl_wv_apply_bc_release(fld->lower_bc[d]);
    if (fld->upper_bc[d])    
      gkyl_wv_apply_bc_release(fld->upper_bc[d]);
  }

  gkyl_array_release(fld->fdup);
  for (int d=0; d<fld->ndim+1; ++d)
    gkyl_array_release(fld->f[d]);
  
  gkyl_array_release(fld->app_current);
  if (fld->proj_app_current)
    gkyl_fv_proj_release(fld->proj_app_current);
  gkyl_array_release(fld->ext_em);
  if (fld->proj_ext_em)
    gkyl_fv_proj_release(fld->proj_ext_em);

  gkyl_array_release(fld->bc_buffer);
}

/** moment_coupling functions */

// initialize source solver: this should be called after all species
// and fields are initialized
static 
void
moment_coupling_init(const struct gkyl_moment_app *app, struct moment_coupling *src)
{
  struct gkyl_moment_em_coupling_inp src_inp = {
    .grid = &app->grid,
    .nfluids = app->num_species,
    .epsilon0 = app->field.epsilon0,
  };

  for (int i=0; i<app->num_species; ++i)
    src_inp.param[i] = (struct gkyl_moment_em_coupling_data) {
      .type = app->species[i].eqn_type,
      .charge = app->species[i].charge,
      .mass = app->species[i].mass,
      .k0 = app->species[i].k0
    };

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);
}

// update sources: 'nstrang' is 0 for the first Strang step and 1 for
// the second step
static 
void
moment_coupling_update(const gkyl_moment_app *app, struct moment_coupling *src,
  int nstrang, double tcurr, double dt)
{
  int sidx[] = { 0, app->ndim };
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];
  
  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = app->species[i].f[sidx[nstrang]];
    
    if (app->species[i].proj_app_accel)
      gkyl_fv_proj_advance(app->species[i].proj_app_accel, tcurr, &app->local, app->species[i].app_accel);
    app_accels[i] = app->species[i].app_accel;
  }
  
  if (app->field.proj_app_current)
    gkyl_fv_proj_advance(app->field.proj_app_current, tcurr, &app->local, app->field.app_current);
  
  if (app->field.proj_ext_em)
    gkyl_fv_proj_advance(app->field.proj_ext_em, tcurr, &app->local, app->field.ext_em);

  gkyl_moment_em_coupling_advance(src->slvr, dt, &app->local,
    fluids, app_accels,
    app->field.f[sidx[nstrang]], app->field.app_current, app->field.ext_em);

  for (int i=0; i<app->num_species; ++i)
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);

  moment_field_apply_bc(app, tcurr, &app->field, app->field.f[sidx[nstrang]]);
}

// free sources
static 
void
moment_coupling_release(const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
}

// initialize non-ideal solver: this should be called after all species
// and fields are initialized
static 
void
moment_non_ideal_init(const struct gkyl_moment_app *app, struct moment_non_ideal *non_ideal)
{
  for (int n=0; n<app->num_species; ++n)
    non_ideal->non_ideal_cflrate[n] = mkarr(1, app->local_ext.volume); 

  // create grid and ranges for non-ideal variables (grid is in computational space)
  // this grid is the grid of node values
  int ghost[3] = { 2, 2, 2 };
  double non_ideal_lower[3] = {0.0};
  double non_ideal_upper[3] = {0.0};
  int non_ideal_cells[3] = {0};
  // braginskii grid has one "extra" cell and is half a grid cell larger past the lower and upper domain
  for (int d=0; d<app->ndim; ++d) {
    non_ideal_lower[d] = app->grid.lower[d] - (app->grid.upper[d]-app->grid.lower[d])/(2.0* (double) app->grid.cells[d]);
    non_ideal_upper[d] = app->grid.upper[d] + (app->grid.upper[d]-app->grid.lower[d])/(2.0* (double) app->grid.cells[d]);
    non_ideal_cells[d] = app->grid.cells[d] + 1;
  }
  gkyl_rect_grid_init(&non_ideal->non_ideal_grid, app->ndim, non_ideal_lower, non_ideal_upper, non_ideal_cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &non_ideal->non_ideal_local_ext, &non_ideal->non_ideal_local);
  // In Braginskii case, non-ideal variables are 6 pressure tensor + 3 heat flux
  for (int n=0;  n<app->num_species; ++n) 
    non_ideal->non_ideal_vars[n] = mkarr(9, non_ideal->non_ideal_local_ext.volume);

  // check if Braginskii terms present
  struct gkyl_moment_braginskii_inp brag_inp = {
    .grid = &app->grid,
    .nfluids = app->num_species,
    .epsilon0 = app->field.epsilon0,
    // Check for multiplicative collisionality factor, default is 1.0
    .coll_fac = app->coll_fac ? 1.0 : app->coll_fac,
  };
  for (int i=0; i<app->num_species; ++i) {
    // Braginskii coefficients depend on pressure and coefficient to obtain
    // pressure is different for different equation systems (gasGamma, vt, Tr(P))
    double p_fac = 1.0;
    if (app->species[i].eqn_type == GKYL_EQN_EULER)
      p_fac =  gkyl_wv_euler_gas_gamma(app->species[i].equation);
    else if (app->species[i].eqn_type == GKYL_EQN_ISO_EULER)
      p_fac =  gkyl_wv_iso_euler_vt(app->species[i].equation);
    brag_inp.param[i] = (struct gkyl_moment_braginskii_data) {
      .type_eqn = app->species[i].eqn_type,
      .type_brag = app->species[i].type_brag,
      .charge = app->species[i].charge,
      .mass = app->species[i].mass,
      .p_fac = p_fac,
    };
  }
  non_ideal->brag_slvr = gkyl_moment_braginskii_new(brag_inp);
}

// compute braginskii variables
static 
void
moment_non_ideal_update(const gkyl_moment_app *app, struct moment_non_ideal *non_ideal, int rk_idx)
{
  int sidx[] = { 0, app->ndim };
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];

  for (int i=0; i<app->num_species; ++i)
    fluids[i] = app->species[i].f[rk_idx];

  // app->field.f[0] are the EM fields at the known time-step
  gkyl_moment_braginskii_advance(non_ideal->brag_slvr,
    non_ideal->non_ideal_local, app->local,
    fluids, app->field.f[0],
    non_ideal->non_ideal_cflrate, non_ideal->non_ideal_vars);
}

// free non-ideal terms
static 
void
moment_non_ideal_release(const gkyl_moment_app *app, const struct moment_non_ideal *non_ideal)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_release(non_ideal->non_ideal_vars[i]);
    gkyl_array_release(non_ideal->non_ideal_cflrate[i]);
  }

  gkyl_moment_braginskii_release(non_ideal->brag_slvr);
}

// rk3 method for KEP + non-ideal terms scheme
static struct gkyl_update_status
kep_rk3(gkyl_moment_app *app, double tcurr, double dt)
{
  int ndim = app->ndim;
  double dt_suggested = DBL_MAX;
  struct gkyl_wave_prop_status stat;

  const struct gkyl_array *fin;
  struct gkyl_array *fout;

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:

        moment_non_ideal_update(app, &app->non_ideal, 0);
        for (int i=0; i<app->num_species; ++i) { 
          gkyl_array_clear(app->species[i].cflrate, 0.0);
          fin = app->species[i].f[0];
          fout = app->species[i].f[1];
          stat = gkyl_kep_scheme_advance(app->species[i].kep_slvr, dt, &app->local, &app->non_ideal.non_ideal_local, 
            fin, app->non_ideal.non_ideal_vars[i], app->species[i].cflrate, fout);
          dt_suggested = fmin(dt_suggested, stat.dt_suggested);
          gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dt_suggested, app->local), 1.0, fin, app->local);
          moment_species_apply_bc(app, tcurr, &app->species[i], fout);
        }

        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:

        moment_non_ideal_update(app, &app->non_ideal, 1);
        for (int i=0; i<app->num_species; ++i) { 
          gkyl_array_clear(app->species[i].cflrate, 0.0);
          fin = app->species[i].f[1];
          fout = app->species[i].f[2];
          stat = gkyl_kep_scheme_advance(app->species[i].kep_slvr, dt, &app->local, &app->non_ideal.non_ideal_local, 
            fin, app->non_ideal.non_ideal_vars[i], app->species[i].cflrate, fout);
          dt_suggested = fmin(dt_suggested, stat.dt_suggested);
          gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dt_suggested, app->local), 1.0, fin, app->local);
          moment_species_apply_bc(app, tcurr+dt, &app->species[i], fout);

          array_combine(app->species[i].f[1], 3.0/4.0, app->species[i].f[0], 1.0/4.0, app->species[i].f[2], app->local_ext);
        }

        state = RK_STAGE_3;

        break;

      case RK_STAGE_3:

        moment_non_ideal_update(app, &app->non_ideal, 1);
        for (int i=0; i<app->num_species; ++i) { 
          gkyl_array_clear(app->species[i].cflrate, 0.0);
          fin = app->species[i].f[1];
          fout = app->species[i].f[2];
          stat = gkyl_kep_scheme_advance(app->species[i].kep_slvr, dt, &app->local, &app->non_ideal.non_ideal_local, 
            fin, app->non_ideal.non_ideal_vars[i], app->species[i].cflrate, fout);
          dt_suggested = fmin(dt_suggested, stat.dt_suggested);
          gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dt_suggested, app->local), 1.0, fin, app->local);
          moment_species_apply_bc(app, tcurr+dt/2, &app->species[i], fout);

          array_combine(app->species[i].f[1], 1.0/3.0, app->species[i].f[0], 2.0/3.0, app->species[i].f[2], app->local_ext);
          gkyl_array_copy_range(app->species[i].f[ndim], app->species[i].f[1], app->local_ext);
        }

        state = RK_COMPLETE;

        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }
  
  return (struct gkyl_update_status) {
    .success = true,
    .dt_suggested = dt_suggested
  };
}

/** app methods */

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment *mom)
{
  struct gkyl_moment_app *app = gkyl_malloc(sizeof(gkyl_moment_app));

  int ndim = app->ndim = mom->ndim;
  strcpy(app->name, mom->name);
  app->tcurr = 0.0; // reset on init

  // create grid and ranges (grid is in computational space)
  int ghost[3] = { 2, 2, 2 };
  gkyl_rect_grid_init(&app->grid, ndim, mom->lower, mom->upper, mom->cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);

  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = mom->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = mom->c2p_ctx;
    app->mapc2p = mom->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, ndim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(ndim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, ndim, mom->mapc2p, mom->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    const char *fmt = "%s-mapc2p.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name);
    gkyl_grid_sub_array_write(&app->grid, &app->local, c2p, fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object
  app->geom = gkyl_wave_geom_new(&app->grid, &app->local_ext,
    app->mapc2p, app->c2p_ctx);

  double cfl_frac = mom->cfl_frac == 0 ? 0.95 : mom->cfl_frac;
  app->cfl = 1.0*cfl_frac;

  app->fluid_scheme = mom->fluid_scheme;

  app->num_periodic_dir = mom->num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = mom->periodic_dirs[d];

  // construct list of directions to skip
  for (int d=0; d<3; ++d)
    app->is_dir_skipped[d] = 0;
  for (int i=0; i<mom->num_skip_dirs; ++i)
    app->is_dir_skipped[mom->skip_dirs[i]] = 1;

  app->has_field = 0;
  // initialize field if we have one
  if (mom->field.init) {
    app->has_field = 1;
    moment_field_init(mom, &mom->field, app, &app->field);
  }

  int ns = app->num_species = mom->num_species;
  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct moment_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i)
    moment_species_init(mom, &mom->species[i], app, &app->species[i]);

  // check if we should update sources
  app->update_sources = 0;
  if (app->has_field && ns>0) {
    app->update_sources = 1; // only update if field and species are present
    moment_coupling_init(app, &app->sources);
  }

  // Initialize Braginskii terms (only used if KEP scheme)
  app->coll_fac = mom->coll_fac;
  moment_non_ideal_init(app, &app->non_ideal);
  // initialize stat object to all zeros
  app->stat = (struct gkyl_moment_stat) {
  };

  return app;
}

double
gkyl_moment_app_max_dt(gkyl_moment_app* app)
{
  double max_dt = DBL_MAX;
  for (int i=0;  i<app->num_species; ++i) 
    max_dt = fmin(max_dt, moment_species_max_dt(app, &app->species[i]));

  if (app->has_field)
    max_dt = fmin(max_dt, moment_field_max_dt(app, &app->field));

  return max_dt;
}

void
gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_moment_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_moment_app_apply_ic_species(app, i, t0);
}

void
gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0)
{
  if (app->has_field != 1) return;

  app->tcurr = t0;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, 2, 8, app->field.init, app->field.ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->field.f[0]);
  gkyl_fv_proj_release(proj);

  moment_field_apply_bc(app, t0, &app->field, app->field.f[0]);
}

void
gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, 2, app->species[sidx].num_equations,
    app->species[sidx].init, app->species[sidx].ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->species[sidx].f[0]);
  gkyl_fv_proj_release(proj);

  moment_species_apply_bc(app, t0, &app->species[sidx], app->species[sidx].f[0]);
}

void
gkyl_moment_app_write(const gkyl_moment_app* app, double tm, int frame)
{
  gkyl_moment_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_app_write_species(app, i, tm, frame);
}

void
gkyl_moment_app_write_field(const gkyl_moment_app* app, double tm, int frame)
{
  if (app->has_field != 1) return;

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, "field", frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "field", frame);
  
  gkyl_grid_sub_array_write(&app->grid, &app->local, app->field.f[0], fileNm);
}

void
gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].name, frame);
  
  gkyl_grid_sub_array_write(&app->grid, &app->local, app->species[sidx].f[0], fileNm);
}

// internal function that takes a single time-step
static
struct gkyl_update_status
moment_update(gkyl_moment_app* app, double dt0)
{
  int ns = app->num_species, ndim = app->ndim;

  double dt_suggested = DBL_MAX;
  
  // time-stepper states
  enum {
    UPDATE_DONE = 0,
    PRE_UPDATE,
    POST_UPDATE,
    FIRST_COUPLING_UPDATE,
    FIELD_UPDATE,
    SPECIES_UPDATE,
    SECOND_COUPLING_UPDATE,
    UPDATE_REDO,
  } state = PRE_UPDATE;

  double tcurr = app->tcurr, dt = dt0;
  while (state != UPDATE_DONE) {
    switch (state) {
      case PRE_UPDATE:
        state = FIRST_COUPLING_UPDATE; // next state
          
        // copy old solution in case we need to redo this step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].fdup, app->species[i].f[0]);
        if (app->has_field)
          gkyl_array_copy(app->field.fdup, app->field.f[0]);

        break;
          
      
      case FIRST_COUPLING_UPDATE:
        state = FIELD_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src1_tm = gkyl_wall_clock();
          moment_coupling_update(app, &app->sources, 0, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src1_tm);
        }
            
        break;

      case FIELD_UPDATE:
        state = SPECIES_UPDATE; // next state

        if (app->has_field) {
          struct timespec fl_tm = gkyl_wall_clock();
          struct gkyl_update_status s = moment_field_update(app, &app->field, tcurr, dt);
          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
            
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
          app->stat.field_tm += gkyl_time_diff_now_sec(fl_tm);
        }
          
        break;

      case SPECIES_UPDATE:
        state = SECOND_COUPLING_UPDATE; // next state

        struct timespec sp_tm = gkyl_wall_clock();
        if (app->fluid_scheme == GKYL_MOMENT_FLUID_KEP) {
          struct gkyl_update_status s = kep_rk3(app, tcurr, dt);

          if (!s.success) {
            app->stat.nfail += 1;
            dt = s.dt_suggested;
            state = UPDATE_REDO;
            break;
          }
          dt_suggested = fmin(dt_suggested, s.dt_suggested);
        }
        else {
          for (int i=0; i<ns; ++i) { 
            struct gkyl_update_status s =
              moment_species_update(app, &app->species[i], tcurr, dt);

            if (!s.success) {
              app->stat.nfail += 1;
              dt = s.dt_suggested;
              state = UPDATE_REDO;
              break;
            }
            dt_suggested = fmin(dt_suggested, s.dt_suggested);
          }
        }

        app->stat.species_tm += gkyl_time_diff_now_sec(sp_tm);
         
        break;

      case SECOND_COUPLING_UPDATE:
        state = POST_UPDATE; // next state

        if (app->update_sources) {
          struct timespec src2_tm = gkyl_wall_clock();
          moment_coupling_update(app, &app->sources, 1, tcurr, dt/2);
          app->stat.sources_tm += gkyl_time_diff_now_sec(src2_tm);
        }

        break;

      case POST_UPDATE:
        state = UPDATE_DONE;

        // copy solution in prep for next time-step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].f[0], app->species[i].f[ndim]);
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.f[ndim]);
          
        break;

      case UPDATE_REDO:
        state = PRE_UPDATE; // start all-over again
          
        // restore solution and retake step
        for (int i=0; i<ns; ++i)
          gkyl_array_copy(app->species[i].f[0], app->species[i].fdup);
        if (app->has_field)
          gkyl_array_copy(app->field.f[0], app->field.fdup);
          
        break;

      case UPDATE_DONE: // unreachable code! (suppresses warning)
        break;
    }
  }

  return (struct gkyl_update_status) {
    .success = true,
    .dt_actual = dt,
    .dt_suggested = dt_suggested,
  };
}

struct gkyl_update_status
gkyl_moment_update(gkyl_moment_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = moment_update(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  
  return status;
}

struct gkyl_moment_stat
gkyl_moment_app_stat(gkyl_moment_app* app)
{
  return app->stat;
}

void
gkyl_moment_app_stat_write(const gkyl_moment_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  with_file (fp, fileNm, "a") {
    fprintf(fp, "{\n");

    if (strftime(buff, sizeof buff, "%c", &curr_tm))
      fprintf(fp, " \"date\" : \"%s\",\n", buff);

    fprintf(fp, " \"nup\" : \"%ld\",\n", app->stat.nup);
    fprintf(fp, " \"nfail\" : \"%ld\",\n", app->stat.nfail);
    fprintf(fp, " \"total_tm\" : \"%lg\",\n", app->stat.total_tm);
    fprintf(fp, " \"species_tm\" : \"%lg\",\n", app->stat.species_tm);
    fprintf(fp, " \"field_tm\" : \"%lg\",\n", app->stat.field_tm);
    fprintf(fp, " \"sources_tm\" : \"%lg\"\n", app->stat.sources_tm);
  
    fprintf(fp, "}\n");
  }
}

void
gkyl_moment_app_release(gkyl_moment_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    moment_species_release(&app->species[i]);
  gkyl_free(app->species);

  if (app->has_field)
    moment_field_release(&app->field);

  if (app->update_sources)
    moment_coupling_release(&app->sources);

  gkyl_wave_geom_release(app->geom);

  gkyl_free(app);
}
