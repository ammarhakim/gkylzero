#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_maxwell.h>
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
  struct gkyl_array *app_accel; // array for applied acceleration/forces
  // pointer to projection operator for applied acceleration/forces function
  gkyl_fv_proj *proj_app_accel;
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  enum gkyl_eqn_type eqn_type; // type ID of equation
  int num_equations; // number of equations in species
  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary condition type
  enum gkyl_species_bc_type lower_bct[3], upper_bct[3];
  // boundary condition solvers on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_q; // integrated conserved quantities
  bool is_first_q_write_call; // flag for dynvec written first time  
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
  struct gkyl_array *app_current; // arrays for applied currents
  // pointer to projection operator for applied current function
  gkyl_fv_proj *proj_app_current;


  bool is_ext_em_static; // flag to indicate if external field is time-independent
  struct gkyl_array *ext_em; // array external fields  
  gkyl_fv_proj *proj_ext_em;   // pointer to projection operator for external fields
  bool was_ext_em_computed; // flag to indicate if we already computed external EM field
  
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary condition type
  enum gkyl_field_bc_type lower_bct[3], upper_bct[3];
  // boundary conditions on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for dynvec written first time
};

// Source data
struct moment_coupling {
  gkyl_moment_em_coupling *slvr; // source solver function
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
    
  struct gkyl_moment_stat stat; // statistics
};

// Function pointer to compute integrated quantities from input
typedef void (*integ_func)(int nc, const double *qin, double *integ_out);

static inline void
integ_unit(int nc, const double *qin, double *integ_out)
{
  for (int i=0; i<nc; ++i) integ_out[i] = qin[i];
}
static inline void
integ_sq(int nc, const double *qin, double *integ_out)
{
  for (int i=0; i<nc; ++i) integ_out[i] = qin[i]*qin[i];
}

// Compute the nc intergated values over the update_rgn, storing the
// result in the integ_q
static void
calc_integ_quant(int nc, double vol, const struct gkyl_array *q, const struct gkyl_wave_geom *geom,
  struct gkyl_range update_rng, integ_func i_func, double *integ_q)
{
  double integ_out[nc];
  for (int i=0; i<nc; ++i) integ_q[i] = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_rng);
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(geom, iter.idx);
    const double *qcell = gkyl_array_cfetch(q, gkyl_range_idx(&update_rng, iter.idx));

    i_func(nc, qcell, integ_out);
    for (int i=0; i<nc; ++i) integ_q[i] += vol*cg->kappa*integ_out[i];
  }
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
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
  sp->num_equations = mom_sp->equation->num_equations;
  // closure parameter, used by 10 moment
  sp->k0 = mom_sp->equation->type == GKYL_EQN_TEN_MOMENT ? gkyl_wv_ten_moment_k0(mom_sp->equation) : 0.0;

  // choose default limiter
  enum gkyl_wave_limiter limiter =
    mom_sp->limiter == 0 ? GKYL_MONOTONIZED_CENTERED : mom_sp->limiter;

  int ndim = mom->ndim;
  // create updaters for each directional update
  for (int d=0; d<ndim; ++d)
    sp->slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &app->grid,
        .equation = mom_sp->equation,
        .limiter = limiter,
        .num_up_dirs = app->is_dir_skipped[d] ? 0 : 1,
        .force_low_order_flux = mom_sp->force_low_order_flux,
        .check_inv_domain = true,
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
        case GKYL_SPECIES_REFLECT:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost,
            mom_sp->bc_lower_func, mom_sp->ctx);
          break;
        
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE: // wedge also uses bc_copy
          sp->lower_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
          break;

        default:
          assert(false);
          break;
      }
      
      // upper BCs
      switch (bc[1]) {
        case GKYL_SPECIES_REFLECT:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->wall_bc_func, 0);
          break;

        case GKYL_SPECIES_NO_SLIP:      
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->equation->no_slip_bc_func, 0);
          break;

        case GKYL_SPECIES_FUNC:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost,
            mom_sp->bc_upper_func, mom_sp->ctx);
          break;
          
        case GKYL_SPECIES_COPY:
        case GKYL_SPECIES_WEDGE:
          sp->upper_bc[dir] = gkyl_wv_apply_bc_new(
            &app->grid, mom_sp->equation, app->geom, dir, GKYL_UPPER_EDGE, nghost, bc_copy, 0);
          break;

        default:
          assert(false);
          break;
      }
    }
  }

  int meqn = sp->num_equations;
  sp->fdup = mkarr(meqn, app->local_ext.volume);
  // allocate arrays
  for (int d=0; d<ndim+1; ++d)
    sp->f[d] = mkarr(meqn, app->local_ext.volume);

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

  sp->integ_q = gkyl_dynvec_new(GKYL_DOUBLE, meqn);
  sp->is_first_q_write_call = true;
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

  gkyl_array_release(sp->app_accel);
  if (sp->proj_app_accel)
    gkyl_fv_proj_release(sp->proj_app_accel);

  gkyl_array_release(sp->bc_buffer);

  gkyl_dynvec_release(sp->integ_q);
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
        .check_inv_domain = false,
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
  fld->is_ext_em_static = mom_fld->is_ext_em_static;

  fld->was_ext_em_computed = false;
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

  fld->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  fld->is_first_energy_write_call = true;
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

  gkyl_dynvec_release(fld->integ_energy);
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
    // if there is a field, need to update electric field too, otherwise just updating fluid
    .epsilon0 = app->field.epsilon0 ? app->field.epsilon0 : 0.0, 
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
moment_coupling_update(gkyl_moment_app *app, struct moment_coupling *src,
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
  
  if (app->field.proj_ext_em) {

    if (!app->field.was_ext_em_computed)
      gkyl_fv_proj_advance(app->field.proj_ext_em, tcurr, &app->local, app->field.ext_em);

    if (app->field.is_ext_em_static)
      app->field.was_ext_em_computed = true;
    else
      app->field.was_ext_em_computed = false;
  }

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

/** app methods */

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment *mom)
{
  disable_denorm_float();
  
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
gkyl_moment_app_write_field_energy(gkyl_moment_app *app)
{
  // write out field energy
  cstr fileNm = cstr_from_fmt("%s-field-energy.gkyl", app->name);

  if (app->field.is_first_energy_write_call) {
    // write to a new file (this ensure previous output is removed)
    gkyl_dynvec_write(app->field.integ_energy, fileNm.str);
    app->field.is_first_energy_write_call = false;
  }
  else {
    // append to existing file
    gkyl_dynvec_awrite(app->field.integ_energy, fileNm.str);
  }
  gkyl_dynvec_clear(app->field.integ_energy);

  cstr_drop(&fileNm);
}

void
gkyl_moment_app_write_integrated_mom(gkyl_moment_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    // write out diagnostic moments
    cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, app->species[i].name,
      "imom");
    
    if (app->species[i].is_first_q_write_call) {
      gkyl_dynvec_write(app->species[i].integ_q, fileNm.str);
      app->species[i].is_first_q_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(app->species[i].integ_q, fileNm.str);
    }
    gkyl_dynvec_clear(app->species[i].integ_q);

    cstr_drop(&fileNm);
  }
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

void
gkyl_moment_app_calc_field_energy(gkyl_moment_app* app, double tm)
{
  double energy[6] = { 0.0 };
  calc_integ_quant(6, app->grid.cellVolume, app->field.f[0], app->geom,
    app->local, integ_sq, energy);
  gkyl_dynvec_append(app->field.integ_energy, tm, energy);
}

void
gkyl_moment_app_calc_integrated_mom(gkyl_moment_app *app, double tm)
{
  for (int sidx=0; sidx<app->num_species; ++sidx) {
    int meqn = app->species[sidx].num_equations;
    double q_integ[meqn];
    calc_integ_quant(meqn, app->grid.cellVolume, app->species[sidx].f[0], app->geom,
      app->local, integ_unit, q_integ);
    gkyl_dynvec_append(app->species[sidx].integ_q, tm, q_integ);
  }
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
      fprintf(fp, " date : %s\n", buff);

    fprintf(fp, " nup : %ld,\n", app->stat.nup);
    fprintf(fp, " nfail : %ld,\n", app->stat.nfail);
    fprintf(fp, " total_tm : %lg,\n", app->stat.total_tm);
    fprintf(fp, " species_tm : %lg,\n", app->stat.species_tm);
    fprintf(fp, " field_tm : %lg,\n", app->stat.field_tm);
    fprintf(fp, " sources_tm : %lg\n", app->stat.sources_tm);

    for (int i=0; i<app->num_species; ++i) {
      for (int d=0; d<app->ndim; ++d) {
        struct gkyl_wave_prop_stats wvs = gkyl_wave_prop_stats(app->species[i].slvr[d]);
        fprintf(fp, " %s_n_bad_advance_calls[%d] = %ld\n", app->species[i].name, d, wvs.n_bad_advance_calls);
        fprintf(fp, " %s_n_bad_cells[%d] = %ld\n", app->species[i].name, d, wvs.n_bad_cells);
        fprintf(fp, " %s_n_max_bad_cells[%d] = %ld\n", app->species[i].name, d, wvs.n_max_bad_cells);
      }
    }
  
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
