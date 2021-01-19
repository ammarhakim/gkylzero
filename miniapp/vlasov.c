#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_mom.h>

// ranges for use in BCs
struct skin_ghost_ranges {
    struct gkyl_range lower_skin[GKYL_MAX_DIM];
    struct gkyl_range lower_ghost[GKYL_MAX_DIM];

    struct gkyl_range upper_skin[GKYL_MAX_DIM];
    struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// data for moments
struct vm_species_moment {
    struct gkyl_mom_type *mtype;
    gkyl_mom_calc *mcalc;
    struct gkyl_array *marr;
};

// species data
struct vm_species {
    struct gkyl_vlasov_species info; // data for species
    
    struct gkyl_rect_grid grid;
    struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
    struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost

    struct gkyl_array *f, *f1, *fnew; // arrays for updates
    struct gkyl_array *cflrate; // CFL rate in each cell
    struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

    struct vm_species_moment m1i; // for computing currents
    struct vm_species_moment *moms; // diagnostic moments

    struct gkyl_dg_eqn *eqn; // Vlasov equation
    gkyl_hyper_dg *slvr; // solver
};

// field data
struct vm_field {
    struct gkyl_em_field info; // data for field
    
    struct gkyl_array *em, *em1, *emnew; // arrays for updates
    struct gkyl_array *qmem; // array for q/m*(E,B)
    struct gkyl_array *cflrate; // CFL rate in each cell
    struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)
    
    struct gkyl_dg_eqn *eqn; // Maxwell equation
    gkyl_hyper_dg *slvr; // solver
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_app {
    char name[128]; // name of app
    int cdim, vdim; // conf, velocity space dimensions
    int poly_order; // polynomial order
    double tcurr; // current time
    double cfl; // CFL number

    int num_periodic_dir; // number of periodic directions
    int periodic_dirs[3]; // list of periodic directions
    
    struct gkyl_rect_grid grid; // config-space grid
    struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
    struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

    struct skin_ghost_ranges skin_ghost; // conf-space skin/ghost

    struct vm_field field; // field data

    // species data
    int num_species;
    struct vm_species *species; // species data

    struct gkyl_vlasov_stat stat; // statistics
};

static inline double
diff_now_tm(struct timespec tm)
{
  return gkyl_time_sec(gkyl_time_diff(tm, gkyl_wall_clock()));
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(sizeof(double[nc]), size);
  gkyl_array_clear(a, 0.0); 
  return a;
}

// Compute out = c1*arr1 + c2*arr2
static inline struct gkyl_array*
array_combine(struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2)
{
  return gkyl_array_accumulate(gkyl_array_set(out, c1, arr1), c2, arr2);
}

// list of valid moment names
static const char *const valid_moment_names[] = { "M0", "M1i", "M2ij", "M2", "M3i" };

// check if name of moment is valid or not
static int
is_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_moment_names)/sizeof(valid_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// initialize local and local-ext ranges on conf-space.
static void
init_conf_ranges(int cdim, const int *cells,
  struct gkyl_range *local, struct gkyl_range *local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(local_ext, cdim, lower_ext, upper_ext);
  gkyl_sub_range_init(local, local_ext, lower, upper);
}

// initialize local and local-ext ranges on phase-space: each species
// can have different phase-space grid and hence ranges; note that
// this function assumes there are no ghost-cells in velocity space.
static void
init_phase_ranges(int cdim, int pdim, const int *cells,
  struct gkyl_range *local, struct gkyl_range *local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  for (int i=cdim; i<pdim; ++i) {
    lower_ext[i] = 0;
    upper_ext[i] = cells[i]-1;

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(local_ext, pdim, lower_ext, upper_ext);
  gkyl_sub_range_init(local, local_ext, lower, upper);
}

// reduce
static inline void
incr_int_array(int ndim, int fact, const int *restrict del,
  const int *restrict inp, int *restrict out)
{
  for (int i=0; i<ndim; ++i)
    out[i] = inp[i] + fact*del[i];
}

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  // create skin ranges in each direction
  for (int d=0; d<ndim; ++d) {

    incr_int_array(ndim, 1, ghost, parent->lower, lo);
    incr_int_array(ndim, -1, ghost, parent->upper, up);
    // adjust for lower skin in direction 'd'
    up[d] = lo[d];
    gkyl_sub_range_init(&sgr->lower_skin[d], parent, lo, up);

    incr_int_array(ndim, 1, ghost, parent->lower, lo);
    incr_int_array(ndim, -1, ghost, parent->upper, up);
    // adjust for upper skin in direction 'd'
    lo[d] = up[d];
    gkyl_sub_range_init(&sgr->upper_skin[d], parent, lo, up);
  }

  // create ghost ranges in each direction
  for (int d=0; d<ndim; ++d) {

    incr_int_array(ndim, 1, ghost, parent->lower, lo);
    incr_int_array(ndim, -1, ghost, parent->upper, up);
    // adjust for lower ghost in direction 'd'
    lo[d] = lo[d]-ghost[d];
    up[d] = lo[d];
    gkyl_sub_range_init(&sgr->lower_ghost[d], parent, lo, up);

    incr_int_array(ndim, 1, ghost, parent->lower, lo);
    incr_int_array(ndim, -1, ghost, parent->upper, up);
    // adjust for upper ghost in direction 'd'
    up[d] = up[d]+ghost[d];
    lo[d] = up[d];
    gkyl_sub_range_init(&sgr->upper_ghost[d], parent, lo, up);
  }  
}

// initialize field object
static void
vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm)
{
  assert(is_moment_name_valid(nm));
  
  sm->mtype = gkyl_vlasov_mom_new(&app->confBasis, &app->basis, nm);
  sm->mcalc = gkyl_mom_calc_new(&s->grid, sm->mtype);
  sm->marr = mkarr(sm->mtype->num_mom*app->confBasis.numBasis,
    app->local_ext.volume);
}

// release memory for moment data object
static void
vm_species_moment_release(struct vm_species_moment *sm)
{
  gkyl_mom_type_release(sm->mtype);
  gkyl_mom_calc_release(sm->mcalc);
  gkyl_array_release(sm->marr);
}

// initialize field object
static void
vm_field_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_field *f)
{
  // allocate EM arrays
  f->em = mkarr(8*app->confBasis.numBasis, app->local_ext.volume);
  f->em1 = mkarr(8*app->confBasis.numBasis, app->local_ext.volume);
  f->emnew = mkarr(8*app->confBasis.numBasis, app->local_ext.volume);
  f->qmem = mkarr(8*app->confBasis.numBasis, app->local_ext.volume);

  // allocate buffer for applying BCs (used for both periodic and copy BCs)
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<app->cdim; ++d) {
    long vol = app->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(8*app->confBasis.numBasis, buff_sz);

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(1, app->local_ext.volume);

  // equation object
  double c = 1/sqrt(f->info.epsilon0*f->info.mu0);
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;
  f->eqn = gkyl_dg_maxwell_new(&app->confBasis, c, ef, mf);

  int up_dirs[] = {0, 1, 2}, zero_flux_flags[] = {0, 0, 0};
  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(&app->grid, &app->confBasis, f->eqn,
    app->cdim, up_dirs, zero_flux_flags, 1);
}

// Compute the RHS for field update, returning maximum stable
// time-step.
static double
vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field,
  const struct gkyl_array *em, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_array_clear(field->cflrate, 0.0);

  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);

  double omegaCfl = gkyl_array_reduce(field->cflrate, GKYL_MAX);

  app->stat.field_rhs_tm += diff_now_tm(wst);
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on EM fields
static void
vm_field_apply_periodic_bc(gkyl_vlasov_app *app, struct vm_field *field, int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, &app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, &app->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, &app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, &app->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on EM fields
static void
vm_field_apply_copy_bc(gkyl_vlasov_app *app, struct vm_field *field, int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, &app->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, &app->skin_ghost.lower_ghost[dir]);

  gkyl_array_copy_to_buffer(field->bc_buffer->data, f, &app->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, field->bc_buffer->data, &app->skin_ghost.upper_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for EM fields
static void
vm_field_apply_bc(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim, is_copy[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_field_apply_periodic_bc(app, field, app->periodic_dirs[d], f);
    is_copy[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_copy[d])
      vm_field_apply_copy_bc(app, field, d, f);
}

// release resources for field
static void
vm_field_release(struct vm_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
  gkyl_array_release(f->qmem);
  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  
  gkyl_dg_eqn_release(f->eqn);
  gkyl_hyper_dg_release(f->slvr);
}

// initialize species object
static void
vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s)
{
  int cdim = app->cdim, vdim = app->vdim, poly_order = app->poly_order;
  int pdim = cdim+vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<cdim; ++d) {
    cells[d] = vm->cells[d];
    lower[d] = vm->lower[d];
    upper[d] = vm->upper[d];
    ghost[d] = 1;
  }
  for (int d=0; d<vdim; ++d) {
    cells[cdim+d] = s->info.cells[d];
    lower[cdim+d] = s->info.lower[d];
    upper[cdim+d] = s->info.upper[d];
    ghost[cdim+d] = 0; // no ghost-cells in velocity space
  }
  gkyl_rect_grid_init(&s->grid, pdim, lower, upper, cells);
  init_phase_ranges(cdim, pdim, cells, &s->local, &s->local_ext);

  // create skin/ghost range
  skin_ghost_ranges_init(&s->skin_ghost, &s->local_ext, ghost);
  
  // allocate distribution function arrays
  s->f = mkarr(app->basis.numBasis, s->local_ext.volume);
  s->f1 = mkarr(app->basis.numBasis, s->local_ext.volume);
  s->fnew = mkarr(app->basis.numBasis, s->local_ext.volume);

  // allocate cflrate (scalar array)
  s->cflrate = mkarr(1, s->local_ext.volume);

  // allocate buffer for applying periodic BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int d=0; d<cdim; ++d) {
    long vol = s->skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  s->bc_buffer = mkarr(app->basis.numBasis, buff_sz);

  // allocate data for momentum (for use in current accumulation)
  vm_species_moment_init(app, s, &s->m1i, "M1i");

  int ndm = s->info.num_diag_moments;
  // allocate data for diagnostic moments
  s->moms = gkyl_malloc(sizeof(struct vm_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    vm_species_moment_init(app, s, &s->moms[m], s->info.diag_moments[m]);

  // create equation object
  s->eqn = gkyl_dg_vlasov_new(&app->confBasis, &app->basis, &app->local);

  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
  }
  for (int d=cdim; d<pdim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = 1; // zero-flux BCs in vel-space
  }
  
  // create solver
  s->slvr = gkyl_hyper_dg_new(&s->grid, &app->basis, s->eqn,
    pdim, up_dirs, zero_flux_flags, 1);
}

// Compute the RHS for species update, returning maximum stable
// time-step.
static double
vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *qmem, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_vlasov_set_qmem(species->eqn, qmem); // must set EM fields to use

  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(species->slvr, &species->local, fin, species->cflrate, rhs);

  double omegaCfl = gkyl_array_reduce(species->cflrate, GKYL_MAX);

  app->stat.species_rhs_tm += diff_now_tm(wst);
  
  return app->cfl/omegaCfl;
}

// Apply periodic BCs on distribution function
static void
vm_species_apply_periodic_bc(gkyl_vlasov_app *app, struct vm_species *species, int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, &species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, &species->skin_ghost.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, &species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, &species->skin_ghost.lower_ghost[dir]);
}

// Apply copy BCs on distribution function
static void
vm_species_apply_copy_bc(gkyl_vlasov_app *app, struct vm_species *species, int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, &species->skin_ghost.lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, &species->skin_ghost.lower_ghost[dir]);

  gkyl_array_copy_to_buffer(species->bc_buffer->data, f, &species->skin_ghost.upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, species->bc_buffer->data, &species->skin_ghost.upper_ghost[dir]);
}

// Determine which directions are periodic and which directions are copy,
// and then apply boundary conditions for distribution function
static void
vm_species_apply_bc(gkyl_vlasov_app *app, struct vm_species *species, struct gkyl_array *f)
{
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim, is_copy[3] = {1, 1, 1};
  for (int d=0; d<num_periodic_dir; ++d) {
    vm_species_apply_periodic_bc(app, species, app->periodic_dirs[d], f);
    is_copy[app->periodic_dirs[d]] = 0;
  }
  for (int d=0; d<cdim; ++d)
    if (is_copy[d])
      vm_species_apply_copy_bc(app, species, d, f);
}

// release resources for species
static void
vm_species_release(struct vm_species *s)
{
  // release various arrays
  gkyl_array_release(s->f);
  gkyl_array_release(s->f1);
  gkyl_array_release(s->fnew);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->bc_buffer);

  // release moment data
  vm_species_moment_release(&s->m1i);
  for (int i=0; i<s->info.num_diag_moments; ++i)
    vm_species_moment_release(&s->moms[i]);
  gkyl_free(s->moms);

  gkyl_dg_eqn_release(s->eqn);
  gkyl_hyper_dg_release(s->slvr);
}

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm vm)
{
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));
  
  int cdim = app->cdim = vm.cdim;
  int vdim = app->vdim = vm.vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm.poly_order;
  int ns = app->num_species = vm.num_species;

  double cfl_frac = vm.cfl_frac == 0 ? 1.0 : vm.cfl_frac;
  app->cfl = cfl_frac/(2*poly_order+1);

  app->num_periodic_dir = vm.num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vm.periodic_dirs[d];

  strcpy(app->name, vm.name);
  app->tcurr = 0.0; // reset on init

  // basis functions
  gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);

  // create config grid & ranges
  gkyl_rect_grid_init(&app->grid, cdim, vm.lower, vm.upper, vm.cells);
  init_conf_ranges(cdim, vm.cells, &app->local, &app->local_ext);

  int ghost[] = { 1, 1, 1 };
  // create config skin/ghost ranges
  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  // initialize EM field
  app->field.info = vm.field;
  vm_field_init(&vm, app, &app->field);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct vm_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i) {
    app->species[i].info = vm.species[i];
    vm_species_init(&vm, app, &app->species[i]);
  }

  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat) {
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };
  
  return app;
}

void
gkyl_vlasov_app_init(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_vlasov_app_init_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_vlasov_app_init_species(app, i, t0);
}

void
gkyl_vlasov_app_init_field(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    poly_order+1, 8, app->field.info.init, app->field.info.ctx);
  
  gkyl_proj_on_basis_advance(proj, t0, &app->local, app->field.em);
  gkyl_proj_on_basis_release(proj);
  vm_field_apply_bc(app, &app->field, app->field.em);
}

void
gkyl_vlasov_app_init_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  int poly_order = app->poly_order;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->species[sidx].grid, &app->basis,
    poly_order+1, 1, app->species[sidx].info.init, app->species[sidx].info.ctx);
  
  gkyl_proj_on_basis_advance(proj, t0, &app->species[sidx].local, app->species[sidx].f);
  gkyl_proj_on_basis_release(proj);
  vm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
}

void
gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m)
      gkyl_mom_calc_advance(s->moms[m].mcalc, &s->local, &app->local,
        s->f, s->moms[m].marr);
  }
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  gkyl_vlasov_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_app_write_species(app, i, tm, frame);
}

void
gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame)
{
  char fileNm[256];
  sprintf(fileNm, "%s-field_%d.gkyl", app->name, frame);
  gkyl_grid_array_write(&app->grid, &app->local, app->field.em, fileNm);
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  char fileNm[256];
  sprintf(fileNm, "%s-%s_%d.gkyl", app->name, app->species[sidx].info.name, frame);
  gkyl_grid_array_write(&app->species[sidx].grid, &app->species[sidx].local,
    app->species[sidx].f, fileNm);
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      char fileNm[256];
      sprintf(fileNm, "%s-%s-%s_%d.gkyl", app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      gkyl_grid_array_write(&app->grid, &app->local, app->species[i].moms[m].marr, fileNm);
    }
  }
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *emout, struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  
  double dtmin = dt;

  // compute RHS of Vlasov equations
  for (int i=0; i<app->num_species; ++i) {
    double qbym = app->species[i].info.charge/app->species[i].info.mass;
    gkyl_array_set(app->field.qmem, qbym, emin);
    
    double dt1 = vm_species_rhs(app, &app->species[i], fin[i], app->field.qmem, fout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // compute RHS of Maxwell equations
  double dt1 = vm_field_rhs(app, &app->field, emin, emout);
  dtmin = fmin(dtmin, dt1);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
    vm_species_apply_bc(app, &app->species[i], fout[i]);
  }

  struct timespec wst = gkyl_wall_clock();
  // accumulate current contribution to electric field terms
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];    
    gkyl_mom_calc_advance(s->m1i.mcalc, &s->local, &app->local,
      fin[i], s->m1i.marr);
    double qbyeps = s->info.charge/app->field.info.epsilon0;
    gkyl_array_accumulate_range(emout, -qbyeps, s->m1i.marr, &app->local);
  }
  app->stat.current_tm += diff_now_tm(wst);
  
  // complete update of field
  gkyl_array_accumulate(gkyl_array_scale(emout, dta), 1.0, emin);
  vm_field_apply_bc(app, &app->field, emout);
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_vlasov_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  struct gkyl_update_status st = { .success = 1 };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
          for (int i=0; i<app->num_species; ++i) {
            fin[i] = app->species[i].f;
            fout[i] = app->species[i].f1;
          }
          forward_euler(app, tcurr, dt, fin, app->field.em, fout, app->field.em1, &st);
          dt = st.dt_actual;
          state = RK_STAGE_2;
          break;

      case RK_STAGE_2:
          for (int i=0; i<app->num_species; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
          forward_euler(app, tcurr+dt, dt, fin, app->field.em1, fout, app->field.emnew, &st);
          if (st.dt_actual < dt) {

            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_2_fail += 1;
            
            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1

          }
          else {
            for (int i=0; i<app->num_species; ++i)
              array_combine(app->species[i].f1, 3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew);
            array_combine(app->field.em1, 3.0/4.0, app->field.em, 1.0/4.0, app->field.emnew);

            state = RK_STAGE_3;
          }
          break;

      case RK_STAGE_3:
          for (int i=0; i<app->num_species; ++i) {
            fin[i] = app->species[i].f1;
            fout[i] = app->species[i].fnew;
          }
          forward_euler(app, tcurr+dt/2, dt, fin, app->field.em1, fout, app->field.emnew, &st);
          if (st.dt_actual < dt) {

            // collect stats
            double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
            app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
              dt_rel_diff);
            app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
              dt_rel_diff);
            app->stat.nstage_3_fail += 1;
            
            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1
            
            dt = st.dt_actual;
            state = RK_STAGE_1; // restart from stage 1
            app->stat.nstage_2_fail += 1;
          }
          else {
            for (int i=0; i<app->num_species; ++i) {
              array_combine(app->species[i].f1, 1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew);
              gkyl_array_copy(app->species[i].f, app->species[i].f1);
            }
            array_combine(app->field.em1, 1.0/3.0, app->field.em, 2.0/3.0, app->field.emnew);
            gkyl_array_copy(app->field.em, app->field.em1);

            state = RK_COMPLETE;
          }
          break;

      case RK_COMPLETE: // can't happen: suppresses warning
          break;
    }
  }
  
  return st;
}

struct gkyl_update_status
gkyl_vlasov_update(gkyl_vlasov_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = rk3(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += diff_now_tm(wst);
  
  return status;
}

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app* app)
{
  return app->stat;
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vm_species_release(&app->species[i]);
  gkyl_free(app->species);
  
  vm_field_release(&app->field);

  gkyl_free(app);
}
