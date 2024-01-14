#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_nodal_ops.h>


void
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = cos(xn[0]);
}

void
test_p1(){
  // create RZ grid
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 64, 128 };


  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rzpoly_order = 1;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, funcdg);
  gkyl_eval_on_nodes_release(eon);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg, "proj_func.gkyl");

  int poly_order = rzpoly_order;
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<rzgrid.ndim; ++d)
      nodes[d] = rzgrid.cells[d] + 1;
  }
                   
  if (poly_order == 2){
    for (int d=0; d<rzgrid.ndim; ++d)
      nodes[d] = 2*(rzgrid.cells[d]) + 1;
  }

  for(int d=0; d<rzgrid.ndim; d++){
    printf("d[%d] = %d\n", d, nodes[d]);
  }
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, rzgrid.ndim, nodes);
  printf("nrange lower = %d %d\n", nrange.lower[0], nrange.lower[1]);
  printf("nrange upper = %d %d\n", nrange.upper[0], nrange.upper[1]);


  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, rzgrid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(&rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg);

  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_nodal_ops_n2m(&rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg2);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg2, "proj_func2.gkyl");

}

void
test_p2(){
  // create RZ grid
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 64, 128 };


  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rzpoly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, funcdg);
  gkyl_eval_on_nodes_release(eon);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg, "proj_func.gkyl");

  int poly_order = rzpoly_order;
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<rzgrid.ndim; ++d)
      nodes[d] = rzgrid.cells[d] + 1;
  }
                   
  if (poly_order == 2){
    for (int d=0; d<rzgrid.ndim; ++d)
      nodes[d] = 2*(rzgrid.cells[d]) + 1;
  }

  for(int d=0; d<rzgrid.ndim; d++){
    printf("d[%d] = %d\n", d, nodes[d]);
  }
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, rzgrid.ndim, nodes);
  printf("nrange lower = %d %d\n", nrange.lower[0], nrange.lower[1]);
  printf("nrange upper = %d %d\n", nrange.upper[0], nrange.upper[1]);


  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, rzgrid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(&rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg);

  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_nodal_ops_n2m(&rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg2);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg2, "proj_func2.gkyl");

}




TEST_LIST = {
  //{ "test_p1", test_p1},
  { "test_p2", test_p2},
  { NULL, NULL },
};
