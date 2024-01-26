
#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>

#include <gkyl_bc_block_tensor.h>

void
proj_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1;
  fout[1] = 0;
  fout[2] = 0;
  fout[3] = 0;
  fout[4] = 1;
  fout[5] = 0;
  fout[6] = 0;
  fout[7] = 0;
  fout[8] = 1;
}

void
proj_cyl(double t, const double *xn, double *fout, void *ctx)
{
  double theta = atan(xn[1]/xn[0]);
  fout[0] = cos(theta);
  fout[1] = sin(theta);
  fout[2] = 0;
  fout[3] = -sin(theta);
  fout[4] = cos(theta);
  fout[5] = 0;
  fout[6] = 0;
  fout[7] = 0;
  fout[8] = 1;
}



void test_cartesian_2x()
{
  double lower[] = { -1, -1 }, upper[] = { 1, 1 };
  int cells[] = { 8, 16 };
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

  // Create cartesian tangent vectors and duals
  struct gkyl_array *dxdz= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  struct gkyl_array *dzdx= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dxdz);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dzdx);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, dzdx, "dzdx.gkyl");

  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid, &local, &local_ext, &basis, false);

  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 1; // second direction
  int idx1[2] = {4,16};
  int idx2[2] = {4,1};
  long loc1 = gkyl_range_idx(&local, idx1);
  long loc2 = gkyl_range_idx(&local, idx2);
  const double *dxdz_i = gkyl_array_fetch(dxdz, loc1);
  const double *dzdx_j = gkyl_array_fetch(dzdx, loc2);

  int num_nodes = 2;
  double tj_i[2*2*num_nodes];

  calc_tensor(bt, dir, edge1, edge2, dzdx_j, dxdz_i, tj_i);
  for(int n = 0; n<num_nodes; n++)
    for (int i = 0; i < 4; i++){
      printf("tj_i[%d][%d] = %g\n", i,n, tj_i[i*num_nodes + n]);
    }

  //printf("\n OR \n");
  //for (int i = 0; i < 8; i++)
  //  printf("tj_i[%d] = %g\n", i, tj_i[i]);
}

void test_cartesian_3x()
{
  double lower[] = { -1, -1 ,-1}, upper[] = { 1, 1 ,1};
  int cells[] = { 8, 16,4 };
  int dim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 ,1};
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // Create cartesian tangent vectors and duals
  struct gkyl_array *dxdz= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  struct gkyl_array *dzdx= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dxdz);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local, dzdx);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, dzdx, "dzdx.gkyl");

  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid, &local, &local_ext, &basis, false);

  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 2; // third direction
  int idx1[3] = {4,8,4};
  int idx2[3] = {4,8,4};
  long loc1 = gkyl_range_idx(&local, idx1);
  long loc2 = gkyl_range_idx(&local, idx2);
  const double *dxdz_i = gkyl_array_fetch(dxdz, loc1);
  const double *dzdx_j = gkyl_array_fetch(dzdx, loc2);

  int num_nodes = 4;
  double tj_i[dim*dim*num_nodes];

  calc_tensor(bt, dir, edge1, edge2, dzdx_j, dxdz_i, tj_i);
  for(int n = 0; n<num_nodes; n++)
    for (int i = 0; i < 9; i++){
      printf("tj_i[%d][%d] = %g\n", i,n, tj_i[i*num_nodes + n]);
    }

  //printf("\n OR \n");
  //for (int i = 0; i < 36; i++)
  //  printf("tj_i[%d] = %g\n", i, tj_i[i]);


}

void test_cyl_cart_2x()
{
  double lower[] = { 0.5, 1 }, upper[] = { 1.5, 1 };
  int cells[] = { 128, 64 };
  struct gkyl_rect_grid grid;
  int dim = sizeof(lower)/sizeof(lower[0]);
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  double lower1[] = { 0.5, 0.5 }, upper1[] = { 1, 1 };
  int cells1[] = { 128, 64 };
  struct gkyl_rect_grid grid1;
  gkyl_rect_grid_init(&grid1, 2, lower1, upper1, cells1);

  //ranges
  struct gkyl_range local1, local_ext1;
  int nghost1[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid1, nghost1, &local_ext1, &local1);

  // Create cartesian tangent vectors and duals
  struct gkyl_array *dxdz= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  struct gkyl_array *dzdx= gkyl_array_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *proj = gkyl_eval_on_nodes_new(&grid1, &basis, 9, &proj_one, 0);
  gkyl_eval_on_nodes *projcyl = gkyl_eval_on_nodes_new(&grid, &basis, 9, &proj_cyl, 0);
  gkyl_eval_on_nodes_advance(projcyl, 0.0, &local, dxdz);
  gkyl_eval_on_nodes_advance(proj, 0.0, &local1, dzdx);
  gkyl_eval_on_nodes_release(proj);
  gkyl_grid_sub_array_write(&grid, &local, dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_write(&grid1, &local1, dzdx, "dzdx.gkyl");

  struct bc_block_tensor *bt = gkyl_bc_block_tensor_new(&grid, &local, &local_ext, &basis, false);

  int edge1 = 1; //upper edge
  int edge2 = 0; //lower edge
  int dir = 1; // second direction
  int idx1[2] = {1,64};
  int idx2[2] = {1,1};
  long loc1 = gkyl_range_idx(&local1, idx1);
  long loc2 = gkyl_range_idx(&local, idx2);
  const double *dxdz_i = gkyl_array_fetch(dxdz, loc1);
  const double *dzdx_j = gkyl_array_fetch(dzdx, loc2);

  int num_nodes = 2;
  double tj_i[2*2*num_nodes];

  calc_tensor(bt, dir, edge1, edge2, dzdx_j, dxdz_i, tj_i);
  for(int n = 0; n<num_nodes; n++)
    for (int i = 0; i < 4; i++){
      printf("tj_i[%d][%d] = %g\n", i,n, tj_i[i*num_nodes + n]);
    }

  //printf("\n OR \n");
  //for (int i = 0; i < 8; i++)
  //  printf("tj_i[%d] = %g\n", i, tj_i[i]);
}


TEST_LIST = {
  { "test_cartesian_2x", test_cartesian_2x},
  { "test_cartesian_3x", test_cartesian_3x},
  { "test_cyl_cart_2x", test_cyl_cart_2x},
  { NULL, NULL },
};
