#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_fem_parproj.h>


void
proj_func(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = cos(5*z)*cos(x);
}

void
proj_func2(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 4*z*cos(2*x);
}

void evalFunc1x_neumannx_dirichletx(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double a = 5.0;
  double c0 = 0.;
  double c1 = a/12. - 1./2.;
  fout[0] = -(1.-a*pow(x,2));
}

void check_continuity(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field)
{
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid.ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));
  // Check continuity
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  int nidx;
  long lin_nidx;
  int idx[3];
  const double *node_i  ;

  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[1] != range.upper[1]) {
      long lidx = gkyl_range_idx(&range, iter.idx);
      idx[0] = iter.idx[0];
      idx[1] = iter.idx[1] + 1;
      long lidx_up = gkyl_range_idx(&range, idx);
      double *arr = gkyl_array_fetch(field, lidx);
      double *arr_up = gkyl_array_fetch(field, lidx_up);
      node_i  = gkyl_array_cfetch(nodes, 2);
      double temp1 = basis.eval_expand(node_i, arr);
      node_i  = gkyl_array_cfetch(nodes, 3);
      double temp2 = basis.eval_expand(node_i, arr);
      node_i  = gkyl_array_cfetch(nodes, 0);
      double temp_up1 = basis.eval_expand(node_i, arr_up);
      node_i  = gkyl_array_cfetch(nodes, 1);
      double temp_up2 = basis.eval_expand(node_i, arr_up);
      TEST_CHECK( gkyl_compare(temp1, temp_up1, 1e-12) );
      TEST_CHECK( gkyl_compare(temp2, temp_up2, 1e-12) );
    }
  }

}

void check_bc(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field1, struct gkyl_array *field2)
{
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid.ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));
  // Check continuity
  int nidx;
  long lin_nidx;
  const double *node_i;
  int remDir[2] = {0,1};
  int locDir[2] = {0, range.upper[1]};
  struct gkyl_range defr;
  gkyl_range_deflate(&defr, &range, remDir, locDir);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &defr);
  while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&defr, iter.idx);
      double *arr1 = gkyl_array_fetch(field1, lidx);
      double *arr2 = gkyl_array_fetch(field2, lidx);
      node_i  = gkyl_array_cfetch(nodes, 2);
      double temp_11 = basis.eval_expand(node_i, arr1);
      node_i  = gkyl_array_cfetch(nodes, 3);
      double temp_12 = basis.eval_expand(node_i, arr1);
      node_i  = gkyl_array_cfetch(nodes, 2);
      double temp_21 = basis.eval_expand(node_i, arr2);
      node_i  = gkyl_array_cfetch(nodes, 3);
      double temp_22 = basis.eval_expand(node_i, arr2);
      TEST_CHECK( gkyl_compare(temp_11, temp_21, 1e-12) );
      TEST_CHECK( gkyl_compare(temp_12, temp_22, 1e-12) );
  }

}


void
test_1(){
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { 3*M_PI/4, 1.0 };
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
  
  // project initial function on 2d field
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_func, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field_discont);
  gkyl_proj_on_basis_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, field_discont, "in_field.gkyl");

  // smooth output
  struct gkyl_array *field= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // smooth
  struct gkyl_array *weight=0;
  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &local_ext, &basis, GKYL_FEM_PARPROJ_DIRICHLET, weight, false);
  gkyl_fem_parproj_set_rhs(parproj, field_discont, field_discont);
  gkyl_fem_parproj_solve(parproj, field);
  gkyl_grid_sub_array_write(&grid, &local, field, "smooth_field.gkyl");

  check_continuity(grid, local, basis, field);
  check_bc(grid, local, basis, field, field_discont);
}

void
test_1_cu(){
  double lower[] = { -M_PI, 0.0 }, upper[] = { 3*M_PI/4, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  // Project initial function on 2d field
  struct gkyl_array *field_discont = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_func, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, field_discont);
  gkyl_proj_on_basis_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, field_discont, "in_field.gkyl");
  struct gkyl_array *field_discont_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(field_discont_dev, field_discont);

  // Allocate output fields and smooth
  struct gkyl_array *field = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *weight = 0;
  struct gkyl_fem_parproj *parproj = gkyl_fem_parproj_new(&local, &local_ext, &basis, GKYL_FEM_PARPROJ_DIRICHLET, weight, true);
  gkyl_fem_parproj_set_rhs(parproj, field_discont_dev, field_discont_dev); // Second argument is bc. Skin cells will be used
  gkyl_fem_parproj_solve(parproj, field_dev);
  gkyl_array_copy(field, field_dev);
  gkyl_grid_sub_array_write(&grid, &local, field, "smooth_field.gkyl");

  check_continuity(grid, local, basis, field);
  check_bc(grid, local, basis, field, field_discont);

  gkyl_array_release(field);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_discont);
  gkyl_array_release(field_discont_dev);

}


TEST_LIST = {
  { "test_1", test_1},
#ifdef GKYL_HAVE_CUDA
  {"test_1_cu", test_1_cu},
#endif
  { NULL, NULL },
};
