#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_moment_braginskii.h>
#include <gkyl_moment_braginskii_priv.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_range.h>
#include <gkyl_rect_apply_bc.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_maxwell.h>

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
  bool has_grad_closure; // has gradient-based closure (only for 10 moment)

  int evolve; // evolve species? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
    
  struct gkyl_array *fdup, *f[4]; // arrays for updates
  struct gkyl_array *app_accel; // array for applied acceleration/forces
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  enum gkyl_eqn_type eqn_type; // type ID of equation
  int num_equations; // number of equations in species
  const struct gkyl_wv_eqn *equation; // equation object for initializing solvers
  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary conditions on lower/upper edges in each direction
  gkyl_rect_apply_bc *lower_bc[3], *upper_bc[3];
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
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary conditions on lower/upper edges in each direction
  gkyl_rect_apply_bc *lower_bc[3], *upper_bc[3];    
};

// Source data
struct moment_coupling {
  gkyl_moment_em_coupling *slvr; // source solver function
  gkyl_moment_braginskii *brag_slvr; // Braginskii solver (if present)
  gkyl_ten_moment_grad_closure *grad_closure_slvr[GKYL_MAX_SPECIES]; // Gradient-based closure solver (if present)
  struct gkyl_array *cflrate; // array for stable time-step from non-ideal terms
  struct gkyl_array *rhs[GKYL_MAX_SPECIES]; // array for storing RHS of each species from non-ideal term updates (Braginskii/Gradient-based closure)
};

// Moment app object: used as opaque pointer in user code
struct gkyl_moment_app {
  char name[128]; // name of app
  int ndim; // space dimensions
  double tcurr; // current time
  double cfl; // CFL number

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  enum gkyl_moment_fluid_scheme fluid_scheme; // scheme to update fluid equations
    
  struct gkyl_rect_grid grid; // grid
  struct gkyl_range local, local_ext; // local, local-ext ranges

  struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  int has_field; // flag to indicate if we have a field
  struct moment_field field; // field data
    
  // species data
  int num_species;
  struct moment_species *species; // species data

  int update_sources; // flag to indicate if sources are to be updated
  enum gkyl_braginskii_type type_brag; // enum for Braginskii type (if present)
  struct moment_coupling sources; // sources
    
  struct gkyl_moment_stat stat; // statistics
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// for use in direction suffle in BCs
static const int maxwell_dir_shuffle[][6] = {
  {0, 1, 2, 3, 4, 5},
  {1, 2, 0, 4, 5, 3},
  {2, 0, 1, 5, 3, 4}
};
static const int euler_dir_shuffle[][3] = {
  {1, 2, 3},
  {2, 3, 1},
  {3, 1, 2}
};
// direction shuffle for off-diagonal components of pressure tensor
static const int ten_moment_dir_shuffle[][3] = {
  {8, 6, 5},
  {6, 5, 8},
  {5, 8, 6}
};

// function for copy BC
static void
bc_copy(double t, int dir, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  for (int c=0; c<nc; ++c) ghost[c] = skin[c];
}

// SHOULD THESE REALLY BE HERE? PERHAPS PUT THEM IN EQN THEMSELVES.

// Maxwell perfectly conducting BC
static void
bc_maxwell_wall(double t, int dir, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  const int *d = maxwell_dir_shuffle[dir];
  
  // zero-tangent for E field
  ghost[d[0]] = skin[d[0]];
  ghost[d[1]] = -skin[d[1]];
  ghost[d[2]] = -skin[d[2]];

  // zero-normal for B field
  ghost[d[3]] = -skin[d[3]];
  ghost[d[4]] = skin[d[4]];
  ghost[d[5]] = skin[d[5]];

  // correction potential
  ghost[6] = -skin[6];
  ghost[7] = skin[7];
}

// Euler perfectly reflecting wall
static void
bc_euler_wall(double t, int dir, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  const int *d = euler_dir_shuffle[dir];

  // copy density and pressure
  ghost[0] = skin[0];
  ghost[4] = skin[4];

  // zero-normal for momentum
  ghost[d[0]] = -skin[d[0]];
  ghost[d[1]] = skin[d[1]];
  ghost[d[2]] = skin[d[2]];
}

// Ten moment perfectly reflecting wall
static void
bc_ten_moment_wall(double t, int dir, int nc, const double *skin, double * GKYL_RESTRICT ghost, void *ctx)
{
  const int *d = euler_dir_shuffle[dir];
  const int *pd = ten_moment_dir_shuffle[dir];

  // copy density and Pxx, Pyy, and Pzz
  ghost[0] = skin[0];
  ghost[4] = skin[4];
  ghost[7] = skin[7];
  ghost[9] = skin[9];

  // zero-normal for momentum
  ghost[d[0]] = -skin[d[0]];
  ghost[d[1]] = skin[d[1]];
  ghost[d[2]] = skin[d[2]];

  // zero-tangent for off-diagonal components of pressure tensor
  ghost[pd[0]] = skin[pd[0]];
  ghost[pd[1]] = -skin[pd[1]];
  ghost[pd[2]] = -skin[pd[2]];
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
  // set closure parameter (default is 0.0, used by 10 moment)
  sp->k0 = mom_sp->k0 == 0 ? 0.0 : mom_sp->k0;
  // check if we are running with gradient-based closure
  sp->has_grad_closure = mom_sp->has_grad_closure == 0 ? 0 : mom_sp->has_grad_closure;
  sp->ctx = mom_sp->ctx;
  sp->init = mom_sp->init;

  sp->eqn_type = mom_sp->equation->type;
  sp->num_equations = mom_sp->equation->num_equations;
  sp->equation = gkyl_wv_eqn_aquire(mom_sp->equation);

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  int ndim = mom->ndim;
  // create updaters for each directional update
  for (int d=0; d<ndim; ++d)
    sp->slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = sp->equation,
        .limiter = limiter,
        .num_up_dirs = 1,
        .update_dirs = { d },
        .cfl = app->cfl
      }
    );

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_moment_bc_type *bc;
      if (dir == 0)
        bc = mom_sp->bcx;
      else if (dir == 1)
        bc = mom_sp->bcy;
      else
        bc = mom_sp->bcz;

      // lower BCs in X
      if (bc[0] == GKYL_MOMENT_SPECIES_WALL) {
        if (sp->eqn_type == GKYL_EQN_EULER)
          sp->lower_bc[dir] = gkyl_rect_apply_bc_new(
            &app->grid, dir, GKYL_LOWER_EDGE, nghost, bc_euler_wall, 0);
        else if (sp->eqn_type == GKYL_EQN_TEN_MOMENT)
          sp->lower_bc[dir] = gkyl_rect_apply_bc_new(
            &app->grid, dir, GKYL_LOWER_EDGE, nghost, bc_ten_moment_wall, 0);
      }
      else {
        sp->lower_bc[dir] = gkyl_rect_apply_bc_new(
          &app->grid, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
      }
      
      // upper BCs in X
      if (bc[1] == GKYL_MOMENT_SPECIES_WALL) {
        if (sp->eqn_type == GKYL_EQN_EULER)
          sp->upper_bc[dir] = gkyl_rect_apply_bc_new(
            &app->grid, dir, GKYL_UPPER_EDGE, nghost, bc_euler_wall, 0);
        else if (sp->eqn_type == GKYL_EQN_TEN_MOMENT)
          sp->upper_bc[dir] = gkyl_rect_apply_bc_new(
            &app->grid, dir, GKYL_UPPER_EDGE, nghost, bc_ten_moment_wall, 0);
      }
      else {
        sp->upper_bc[dir] = gkyl_rect_apply_bc_new(
          &app->grid, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
      }
    }
    else {
      sp->upper_bc[dir] = sp->lower_bc[dir] = 0;
    }
  }
    
  int meqn = sp->num_equations;
  sp->fdup = mkarr(meqn, app->local_ext.volume);
  // allocate arrays
  for (int d=0; d<ndim+1; ++d)
    sp->f[d] = mkarr(meqn, app->local_ext.volume);

  // allocate array for applied acceleration/forces for each species
  sp->app_accel = mkarr(3, app->local_ext.volume);
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
      gkyl_rect_apply_bc_advance(sp->lower_bc[d], tcurr, &app->local, f);
      gkyl_rect_apply_bc_advance(sp->upper_bc[d], tcurr, &app->local, f);
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
  for (int d=0; d<sp->ndim; ++d)
    gkyl_wave_prop_release(sp->slvr[d]);

  for (int d=0; d<sp->ndim; ++d) {
    if (sp->lower_bc[d])
      gkyl_rect_apply_bc_release(sp->lower_bc[d]);
    if (sp->upper_bc[d])    
      gkyl_rect_apply_bc_release(sp->upper_bc[d]);
  }

  gkyl_array_release(sp->fdup);
  for (int d=0; d<sp->ndim+1; ++d)
    gkyl_array_release(sp->f[d]);

  gkyl_array_release(sp->app_accel);

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
        .num_up_dirs = 1,
        .update_dirs = { d },
        .cfl = app->cfl
      }
    );

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d)
    is_np[app->periodic_dirs[d]] = 0;

  int nghost[3] = {2, 2, 2};
  for (int dir=0; dir<app->ndim; ++dir) {
    if (is_np[dir]) {
      const enum gkyl_moment_bc_type *bc;
      if (dir == 0)
        bc = mom_fld->bcx;
      else if (dir == 1)
        bc = mom_fld->bcy;
      else
        bc = mom_fld->bcz;

      // lower BCs in X
      if (bc[0] == GKYL_MOMENT_FIELD_COND)
        fld->lower_bc[dir] = gkyl_rect_apply_bc_new(
          &app->grid, dir, GKYL_LOWER_EDGE, nghost, bc_maxwell_wall, 0);
      else
        fld->lower_bc[dir] = gkyl_rect_apply_bc_new(
          &app->grid, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
      
      // upper BCs in X
      if (bc[1] == GKYL_MOMENT_FIELD_COND)
        fld->upper_bc[dir] = gkyl_rect_apply_bc_new(
          &app->grid, dir, GKYL_UPPER_EDGE, nghost, bc_maxwell_wall, 0);
      else
        fld->upper_bc[dir] = gkyl_rect_apply_bc_new(
          &app->grid, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
    }
    else {
      fld->upper_bc[dir] = fld->lower_bc[dir] = 0;
    }
  }

  // allocate arrays
  fld->fdup = mkarr(8, app->local_ext.volume);
  for (int d=0; d<ndim+1; ++d)
    fld->f[d] = mkarr(8, app->local_ext.volume);

  // allocate arrays for applied current/external fields
  fld->app_current = mkarr(3, app->local_ext.volume);
  fld->ext_em = mkarr(6, app->local_ext.volume);

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
      gkyl_rect_apply_bc_advance(field->lower_bc[d], tcurr, &app->local, f);
      gkyl_rect_apply_bc_advance(field->upper_bc[d], tcurr, &app->local, f);
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
      gkyl_rect_apply_bc_release(fld->lower_bc[d]);
    if (fld->upper_bc[d])    
      gkyl_rect_apply_bc_release(fld->upper_bc[d]);
  }

  gkyl_array_release(fld->fdup);
  for (int d=0; d<fld->ndim+1; ++d)
    gkyl_array_release(fld->f[d]);
  
  gkyl_array_release(fld->app_current);
  gkyl_array_release(fld->ext_em);

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
      .k0 = app->species[i].has_grad_closure ? 0.0 : app->species[i].k0,
    };

  // create updater to solve for sources
  src->slvr = gkyl_moment_em_coupling_new(src_inp);

  for (int n=0; n<app->num_species; ++n) {
    int meqn = app->species[n].num_equations;
    src->rhs[n] = mkarr(meqn, app->local_ext.volume);
  }
  src->cflrate = mkarr(1, app->local_ext.volume); 
  
  // check if Braginskii terms present
  if (app->type_brag) {
    struct gkyl_moment_braginskii_inp brag_inp = {
      .grid = &app->grid,
      .nfluids = app->num_species,
      .epsilon0 = app->field.epsilon0,
      .type_brag = app->type_brag,
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
        .charge = app->species[i].charge,
        .mass = app->species[i].mass,
        .p_fac = p_fac,
      };
    }
    src->brag_slvr = gkyl_moment_braginskii_new(brag_inp);
  }

  // check if gradient-closure is present
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      struct gkyl_ten_moment_grad_closure_inp grad_closure_inp = {
        .grid = &app->grid,
        .k0 = app->species[i].k0,
      };
      src->grad_closure_slvr[i] = gkyl_ten_moment_grad_closure_new(grad_closure_inp);      
    }
  }
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
  const struct gkyl_array *rhs_const[GKYL_MAX_SPECIES];
  
  for (int i=0; i<app->num_species; ++i) {
    fluids[i] = app->species[i].f[sidx[nstrang]];
    app_accels[i] = app->species[i].app_accel;
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure) {
      gkyl_ten_moment_grad_closure_advance(src->grad_closure_slvr[i], app->local, app->species[i].f[sidx[nstrang]], app->field.f[sidx[nstrang]], src->cflrate, src->rhs[i]);
    }
  }

  if (app->type_brag) {
    gkyl_moment_braginskii_advance(src->brag_slvr, app->local, fluids, app->field.f[sidx[nstrang]], src->cflrate, src->rhs);
  }

  for (int i=0; i<app->num_species; ++i)
    rhs_const[i] = src->rhs[i];
  
  gkyl_moment_em_coupling_advance(src->slvr, dt, app->local,
    fluids, app_accels, rhs_const,
    app->field.f[sidx[nstrang]], app->field.app_current, app->field.ext_em);

  for (int i=0; i<app->num_species; ++i)
    moment_species_apply_bc(app, tcurr, &app->species[i], fluids[i]);

  moment_field_apply_bc(app, tcurr, &app->field, app->field.f[sidx[nstrang]]);
}

// free sources
static 
void
moment_coupling_release(const gkyl_moment_app *app, const struct moment_coupling *src)
{
  gkyl_moment_em_coupling_release(src->slvr);
  if (app->type_brag)
    gkyl_moment_braginskii_release(src->brag_slvr);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_release(src->rhs[i]);
    if (app->species[i].eqn_type == GKYL_EQN_TEN_MOMENT && app->species[i].has_grad_closure)
      gkyl_ten_moment_grad_closure_release(src->grad_closure_slvr[i]);
  }
  gkyl_array_release(src->cflrate);
}

/** app methods */

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment mom)
{
  struct gkyl_moment_app *app = gkyl_malloc(sizeof(gkyl_moment_app));

  int ndim = app->ndim = mom.ndim;

  // create grid and ranges
  int ghost[3] = { 2, 2, 2 };
  gkyl_rect_grid_init(&app->grid, ndim, mom.lower, mom.upper, mom.cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);

  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  double cfl_frac = mom.cfl_frac == 0 ? 0.95 : mom.cfl_frac;
  app->cfl = 1.0*cfl_frac;

  app->fluid_scheme = mom.fluid_scheme;

  app->num_periodic_dir = mom.num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = mom.periodic_dirs[d];

  app->has_field = 0;
  // initialize field if we have one
  if (mom.field.init) {
    app->has_field = 1;
    moment_field_init(&mom, &mom.field, app, &app->field);
  }

  int ns = app->num_species = mom.num_species;
  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct moment_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i)
    moment_species_init(&mom, &mom.species[i], app, &app->species[i]);

  // check if we should update sources
  app->update_sources = 0;
  if (app->has_field && ns>0) {
    app->update_sources = 1; // only update if field and species are present
    moment_coupling_init(app, &app->sources);
  }

  strcpy(app->name, mom.name);
  app->tcurr = 0.0; // reset on init

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
  free(app->species);

  if (app->has_field)
    moment_field_release(&app->field);

  if (app->update_sources)
    moment_coupling_release(app, &app->sources);
  
  gkyl_free(app);
}
