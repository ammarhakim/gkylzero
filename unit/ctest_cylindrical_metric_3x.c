#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_kernels.h>


void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  fout[0] = r*cos(theta); fout[1] = r*sin(theta); fout[2] = z;
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = r*r;
  fout[4] = 0.0;
  fout[5] = 1.0;
}

void
test_3x_p1()
{
struct gkyl_basis basis;
int poly_order = 1;
gkyl_cart_modal_serendip(&basis, 3, poly_order);

double lower[3] = { 0.5, -M_PI/3, -1 }, upper[3] = { 1.0, M_PI/3, 1 };
int cells[3] = { 8, 80, 8 };
struct gkyl_rect_grid grid;
gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

struct gkyl_range ext_range, range;
int nghost[3] = { 1,1,1};
gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
struct gkyl_range_iter iter;

//calculate metrics directly
gkyl_eval_on_nodes *eval_metric = gkyl_eval_on_nodes_new(&grid, &basis, 6, exact_gij, 0);
struct gkyl_array *gFldExact = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_metric, 0.0, &ext_range, gFldExact);

//calculate the metrics with recovery method
gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 3, mapc2p, 0);
struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume);
struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &ext_range, XYZ);

int bcs[3] = {0,0,0};
gkyl_calc_metric *calculator = gkyl_calc_metric_new(&basis, &grid, bcs, false);
gkyl_calc_metric_advance( calculator, &range, XYZ, gFld);

gkyl_grid_sub_array_write(&grid, &range, gFld, "cylindrical_gFld.gkyl");

gkyl_range_iter_init(&iter, &range);
//while (gkyl_range_iter_next(&iter)) {
//  long loc = gkyl_range_idx(&range, iter.idx);
//  double *gij = gkyl_array_fetch(gFld, loc);
//  double *gijExact = gkyl_array_fetch(gFldExact, loc);
//  double *g22  = &gij[3 * 8];
//  double *g22Exact  = &gijExact[3 * 8];
//  double z[3] = {0.};
//  double val = basis.eval_expand(z, g22);
//  double valExact = basis.eval_expand(z, g22Exact);
//  printf("iter.idx = %d, %d, %d\n", iter.idx[0], iter.idx[1], iter.idx[2]);
//  printf("val, valExact = %g, %g\n", val, valExact);
//  TEST_CHECK( gkyl_compare(val, valExact, 1e-10) );
//  //for(int i = 0; i < basis.num_basis; i++){
//  // //double val = *gij;
//  // //double valExact = *gijExact;
//  // if(iter.idx[1] == 1){
//  //   printf("iter.idx = %d, %d, %d\n", iter.idx[0], iter.idx[1], iter.idx[2]);
//  //   printf("val, valExact = %g, %g\n", g22[i], g22Exact[i]);
//  //   TEST_CHECK( gkyl_compare(g22[i], g22Exact[i], 1e-10) );
//  // }
//  // TEST_CHECK( gkyl_compare(g22[i], g22Exact[i], 1e-12) );
//  //}
//}

gkyl_eval_on_nodes_release(eval_metric);
gkyl_eval_on_nodes_release(eval_mapc2p);
gkyl_array_release(XYZ);
gkyl_array_release(gFld);
gkyl_array_release(gFldExact);
gkyl_calc_metric_release(calculator);
}

//void
//test_3x_p2()
//{
//struct gkyl_basis basis;
//int poly_order = 2;
//gkyl_cart_modal_serendip(&basis, 3, poly_order);
//
//double lower[3] = { 0.5, -M_PI/2, 0 }, upper[3] = { 1.0, M_PI/2, M_PI/2 };
//int cells[3] = { 8, 8, 8 };
//struct gkyl_rect_grid grid;
//gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
//
//struct gkyl_range ext_range, range;
//int nghost[3] = { 1,1,1};
//gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
//struct gkyl_range_iter iter;
//
////calculate metrics directly
//gkyl_eval_on_nodes *eval_metric = gkyl_eval_on_nodes_new(&grid, &basis, 6, exact_gij, 0);
//struct gkyl_array *gFldExact = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, ext_range.volume);
//gkyl_eval_on_nodes_advance(eval_metric, 0.0, &ext_range, gFldExact);
//
////calculate the metrics with recovery method
//gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 3, mapc2p, 0);
//struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume); //on ghosts
//struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, ext_range.volume);
//gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &ext_range, XYZ); //on ghosts with ext_range
//
//gkyl_calc_metric *calculator = gkyl_calc_metric_new(&basis, &grid, false);
//gkyl_calc_metric_advance( calculator, &range, XYZ, gFld);
//
//gkyl_range_iter_init(&iter, &range);
//while (gkyl_range_iter_next(&iter)) {
//  long loc = gkyl_range_idx(&range, iter.idx);
//  double *gij = gkyl_array_fetch(gFld, loc);
//  double *gijExact = gkyl_array_fetch(gFldExact, loc);
//  double val = *gij;
//  double valExact = *gijExact;
//  TEST_CHECK( gkyl_compare(val, valExact, 1e-12) );
//}
//
//gkyl_eval_on_nodes_release(eval_metric);
//gkyl_eval_on_nodes_release(eval_mapc2p);
//gkyl_array_release(XYZ);
//gkyl_array_release(gFld);
//gkyl_array_release(gFldExact);
//gkyl_calc_metric_release(calculator);
//}

TEST_LIST = {
  { "test_3x_p1", test_3x_p1},
  //{ "test_3x_p2", test_3x_p2},
  { NULL, NULL },
};
