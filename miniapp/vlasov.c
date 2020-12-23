#include <gkyl_alloc.h>
#include <gkyl_array.h>
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
};

struct gkyl_vlasov_app {
    struct gkyl_vm vm;
    struct gkyl_rect_grid confGrid;
    struct gkyl_basis basis, confBasis;

    struct vm_species *species;
};

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm vm)
{
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));

  app->vm = vm;
  int cdim = vm.cdim, vdim = vm.vdim, pdim = cdim+vdim;
  int ns = vm.num_species;

  if (ns > 0)
    app->species = gkyl_malloc(ns*sizeof(struct vm_species));
  else
    app->species = 0; // so free won't barf

  // create config grid
  gkyl_rect_grid_init(&app->confGrid, cdim, vm.lower, vm.upper, vm.cells);

  // create species grids
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
  }

  // basis functions
  int poly_order = vm.poly_order;
  gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
    
  return app;
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  gkyl_free(app->species);
  gkyl_free(app);
}
