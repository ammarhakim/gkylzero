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
#include <gkyl_deflate_zsurf.h>


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

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, rzgrid.ndim, nodes);


  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, rzgrid.ndim, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&rzbasis, &rzgrid, false);
  gkyl_nodal_ops_m2n(n2m, &rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg);

  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_nodal_ops_n2m(n2m, &rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg2);
  gkyl_nodal_ops_release(n2m);
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

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, rzgrid.ndim, nodes);


  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, rzgrid.ndim, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&rzbasis, &rzgrid, false);
  gkyl_nodal_ops_m2n(n2m, &rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg);

  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_nodal_ops_n2m(n2m, &rzbasis, &rzgrid, &nrange, &rzlocal, 1, nodal_fld, funcdg2);
  gkyl_nodal_ops_release(n2m);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg2, "proj_func2.gkyl");

}

void
test_p1_cu(){
  // create RZ grid
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 8, 16 };


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

  struct gkyl_basis *rzbasis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(rzbasis_on_dev, 2, rzpoly_order);

  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid, &rzbasis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, funcdg);
  gkyl_eval_on_nodes_release(eon);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg, "proj_func.gkyl");
  struct gkyl_array *funcdg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_array_copy(funcdg_dev, funcdg);

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

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, rzgrid.ndim, nodes);


  struct gkyl_array* nodal_fld_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, rzgrid.ndim, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&rzbasis, &rzgrid, true);
  gkyl_nodal_ops_m2n(n2m, rzbasis_on_dev, &rzgrid, &nrange, &rzlocal, 1, nodal_fld_dev, funcdg_dev);

  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  struct gkyl_array *funcdg2_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  gkyl_nodal_ops_n2m(n2m, rzbasis_on_dev, &rzgrid, &nrange, &rzlocal, 1, nodal_fld_dev, funcdg2_dev);
  gkyl_nodal_ops_release(n2m);
  gkyl_array_copy(funcdg2, funcdg2_dev);
  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, funcdg2, "proj_func2.gkyl");

}

void
test_p1_deflated(){
  // create  grid
  double lower[] = { 0.0, -1.5 }, upper[] = { 1.5, 1.5 };
  // as ellipitical surfaces are exact, we only need 1 cell in each
  // direction
  int cells[] = { 4, 6 };


  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //  ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  //  basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  struct gkyl_array *funcdg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_func, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, funcdg);
  gkyl_eval_on_nodes_release(eon);
  gkyl_grid_sub_array_write(&grid, &local, funcdg, "proj_func.gkyl");

  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = grid.cells[d] + 1;
  }
                   
  if (poly_order == 2){
    for (int d=0; d<grid.ndim; ++d)
      nodes[d] = 2*(grid.cells[d]) + 1;
  }

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);


  struct gkyl_array* nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);

  // Create deflated 1d grid, ranges, basis, and nodal range
  struct gkyl_rect_grid deflated_grid;
  struct gkyl_basis deflated_basis;
  struct gkyl_range deflated_local, deflated_local_ext;
  double deflated_lower[1] = { grid.lower[0]}, deflated_upper[1] = { grid.upper[0]};
  int deflated_cells[1] = { grid.cells[0] };
  gkyl_rect_grid_init(&deflated_grid, 1, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = { 1 };
  gkyl_create_grid_ranges(&deflated_grid, deflated_nghost, &deflated_local_ext, &deflated_local);
  gkyl_cart_modal_serendip(&deflated_basis, 1, poly_order);

  struct gkyl_basis *deflated_basis_on_dev = &deflated_basis;
  nodal_fld = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);

  int deflated_nodes[1];
  if (poly_order == 1){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = deflated_grid.cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<deflated_grid.ndim; ++d)
      deflated_nodes[d] = 2*(deflated_grid.cells[d]) + 1;
  }
  struct gkyl_range deflated_nrange;
  gkyl_range_init_from_shape(&deflated_nrange, deflated_grid.ndim, deflated_nodes);
  struct gkyl_array *deflated_phi = gkyl_array_new(GKYL_DOUBLE, deflated_basis.num_basis, deflated_nrange.volume);

  struct gkyl_deflate_zsurf *deflator_lo = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 0, false);
  struct gkyl_deflate_zsurf *deflator_up = gkyl_deflate_zsurf_new(&basis, &deflated_basis, 1, false);
  struct gkyl_nodal_ops *n2m_1d = gkyl_nodal_ops_new(&deflated_basis, &deflated_grid, false);

  for(int zidx = local.lower[1]; zidx <= local.upper[1]; zidx++){
    gkyl_deflate_zsurf_advance(deflator_lo, zidx, &local, &deflated_local, funcdg, deflated_phi, 1);
    gkyl_nodal_ops_m2n_deflated(n2m_1d, deflated_basis_on_dev, &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, nodal_fld, deflated_phi, zidx-1);
    if (zidx == local.upper[1]) {
      gkyl_deflate_zsurf_advance(deflator_up, zidx, &local, &deflated_local, funcdg, deflated_phi, 1);
      gkyl_nodal_ops_m2n_deflated(n2m_1d, deflated_basis_on_dev, &deflated_grid, &nrange, &deflated_nrange, &deflated_local, 1, nodal_fld, deflated_phi, zidx);
    }
  }


  struct gkyl_array *funcdg2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_nodal_ops_n2m(n2m, &basis, &grid, &nrange, &local, 1, nodal_fld, funcdg2);
  gkyl_nodal_ops_release(n2m);
  gkyl_grid_sub_array_write(&grid, &local, funcdg2, "proj_func2.gkyl");

}








TEST_LIST = {
  //{ "test_p1", test_p1},
  { "test_p1_deflated", test_p1_deflated},
  //{ "test_p2", test_p2},
#ifdef GKYL_HAVE_CUDA
  { "test_p1_cu", test_p1_cu},
#endif
  { NULL, NULL },
};
