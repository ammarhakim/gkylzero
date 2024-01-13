#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>


void
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = z*cos(xn[0]);
}

void
test_1(){
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { M_PI, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  // project initial function
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_func, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field);
  gkyl_proj_on_basis_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, field, "xzfunc.gkyl");

  // create deflated grid, ranges, basis, and field
  // create xz grid
  double deflated_lower[] = { -M_PI}, deflated_upper[] = { M_PI};
  int deflated_cells[] = { 12};
  struct gkyl_rect_grid deflated_grid;
  gkyl_rect_grid_init(&deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);

  //ranges
  struct gkyl_range deflated_local, deflated_local_ext;
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);

  // basis function
  int deflated_poly_order = 1;
  struct gkyl_basis deflated_basis;
  gkyl_cart_modal_serendip(&deflated_basis, 1, deflated_poly_order);
  
  //field
  struct gkyl_array *deflated_field = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_local_ext.volume);

  // now deflate
  int edge = 1; //lower = 0
  gkyl_deflate_zsurf *deflator = gkyl_deflate_zsurf_new(&basis, &deflated_basis, &grid, &deflated_grid, edge, false);

  int zidx = 1;
  gkyl_deflate_zsurf_advance(deflator, zidx, &local, &deflated_local, field, deflated_field, 1);




  gkyl_grid_sub_array_write(&deflated_grid, &deflated_local, deflated_field, "xzfunc_deflated.gkyl");

}


TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
