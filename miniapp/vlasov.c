#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_mom_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_mom.h>

// species data
struct vm_species {
    struct gkyl_rect_grid grid;
    struct gkyl_range local, local_ext; // local, local-ext phase-space ranges

    // arrays for distribution function updates
    struct gkyl_array *f;
};

struct gkyl_vlasov_app {
    struct gkyl_vm vm;
    struct gkyl_rect_grid grid; // config-space grid
    struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
    struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

    // arrays for field-updates
    struct gkyl_array *em;

    // species data
    struct vm_species *species;
};

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
  gkyl_range_init(local, cdim, lower, upper);    
  gkyl_range_init(local_ext, cdim, lower_ext, upper_ext);
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
  gkyl_range_init(local, pdim, lower, upper);
  gkyl_range_init(local_ext, pdim, lower_ext, upper_ext);
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(sizeof(double)*nc, size);
  gkyl_array_clear(a, 0.0);
  return a;
}

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm vm)
{
  int cdim = vm.cdim, vdim = vm.vdim, pdim = cdim+vdim;
  int poly_order = vm.poly_order, ns = vm.num_species;

  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));  
  app->vm = vm;

  // basis functions
  gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);

  if (ns > 0)
    app->species = gkyl_malloc(ns*sizeof(struct vm_species));
  else
    app->species = 0; // so free won't barf

  // create config grid & ranges
  gkyl_rect_grid_init(&app->grid, cdim, vm.lower, vm.upper, vm.cells);
  init_conf_ranges(cdim, vm.cells, &app->local, &app->local_ext);

  // allocate field arrays (6 field components + 2 error fields)
  app->em = mkarr(8*app->confBasis.numBasis, app->local_ext.volume);
    
  // create species grid & ranges
  for (int i=0; i<ns; ++i) {
    int cells[GKYL_MAX_DIM];
    double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

    for (int d=0; d<cdim; ++d) {
      cells[d] = vm.cells[d];
      lower[d] = vm.lower[d];
      upper[d] = vm.upper[d];
    }
    for (int d=0; d<vdim; ++d) {
      cells[cdim+d] = vm.species[i].cells[d];
      lower[cdim+d] = vm.species[i].lower[d];
      upper[cdim+d] = vm.species[i].upper[d];
    }
    gkyl_rect_grid_init(&app->species[i].grid, pdim, lower, upper, cells);
    init_phase_ranges(cdim, pdim, cells,
      &app->species[i].local, &app->species[i].local_ext);

    // allocate distribution function arrays
    app->species[i].f = mkarr(app->basis.numBasis, app->species[i].local_ext.volume);
  }
  
  return app;
}

void
gkyl_vlasov_app_init_sim(gkyl_vlasov_app* app)
{
  int poly_order = app->vm.poly_order;

  do { // initialize EM field
    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
      poly_order+1, 8, app->vm.field.init, app->vm.field.ctx);
    gkyl_proj_on_basis_advance(proj, 0.0, &app->local, app->em);
    gkyl_proj_on_basis_release(proj);
  } while(0);

  // initialize species
  for (int i=0;  i<app->vm.num_species; ++i) {
    gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->species[i].grid, &app->basis,
      poly_order+1, 1, app->vm.species[i].init, app->vm.species[i].ctx);
    gkyl_proj_on_basis_advance(proj, 0.0, &app->species[i].local, app->species[i].f);
    gkyl_proj_on_basis_release(proj);
  }
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  do { // write EM field
    char fileNm[256];
    sprintf(fileNm, "%s-field_%d.gkyl", app->vm.name, frame);
    gkyl_grid_array_write(&app->grid, &app->local, app->em, fileNm);
  } while(0);

  // write species distribution function
  for (int i=0; i<app->vm.num_species; ++i) {
    char fileNm[256];
    sprintf(fileNm, "%s-%s_%d.gkyl", app->vm.name, app->vm.species[i].name, frame);
    gkyl_grid_array_write(&app->species[i].grid, &app->species[i].local,
      app->species[i].f, fileNm);
  }
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->vm.num_species; ++i) {
    gkyl_array_release(app->species[i].f);
  }
  gkyl_free(app->species);

  gkyl_array_release(app->em);
  gkyl_free(app);
}
