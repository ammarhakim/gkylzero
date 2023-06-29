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
  double r = xn[0], theta = xn[1];
  fout[0] = r*cos(theta); fout[1] = r*sin(theta);
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1];
  fout[0] = 1.0; //g11
  fout[1] = 0.0; //g12
  fout[2] = r*r; //g22
}

void
test_2x_p1()
{
struct gkyl_basis basis;
int poly_order = 1;
gkyl_cart_modal_serendip(&basis, 2, poly_order);

double lower[2] = { 0.5, -M_PI/2}, upper[2] = { 1.0, M_PI/2};
int cells[2] = { 8, 8 };
struct gkyl_rect_grid grid;
gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

struct gkyl_range ext_range, range;
int nghost[2] = { 1, 1};
gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
struct gkyl_range_iter iter;

//calculate metrics directly
gkyl_eval_on_nodes *eval_metric = gkyl_eval_on_nodes_new(&grid, &basis, 3, exact_gij, 0);
struct gkyl_array *gFldExact = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_metric, 0.0, &ext_range, gFldExact);

//calculate the metrics with recovery method
gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 2, mapc2p, 0);
struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, ext_range.volume);
struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &ext_range, XYZ);

gkyl_calc_metric *calculator = gkyl_calc_metric_new(&basis, &grid, false);
gkyl_calc_metric_advance( calculator, &range, XYZ, gFld);

gkyl_range_iter_init(&iter, &range);
while (gkyl_range_iter_next(&iter)) {
  long loc = gkyl_range_idx(&range, iter.idx);
  double *gij = gkyl_array_fetch(gFld, loc);
  double *gijExact = gkyl_array_fetch(gFldExact, loc);
  double val = *gij;
  double valExact = *gijExact;
  TEST_CHECK( gkyl_compare(val, valExact, 1e-12) );
}
 

gkyl_eval_on_nodes_release(eval_metric);
gkyl_eval_on_nodes_release(eval_mapc2p);
gkyl_array_release(XYZ);
gkyl_array_release(gFld);
gkyl_array_release(gFldExact);
gkyl_calc_metric_release(calculator);
}

void
test_2x_p2()
{
struct gkyl_basis basis;
int poly_order = 2;
gkyl_cart_modal_serendip(&basis, 2, poly_order);

double lower[2] = { 0.5, -M_PI/2 }, upper[2] = { 1.0, M_PI/2 };
int cells[2] = { 8, 8 };
struct gkyl_rect_grid grid;
gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

struct gkyl_range ext_range, range;
int nghost[2] = { 1, 1};
gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
struct gkyl_range_iter iter;

//calculate metrics directly
gkyl_eval_on_nodes *eval_metric = gkyl_eval_on_nodes_new(&grid, &basis, 3, exact_gij, 0);
struct gkyl_array *gFldExact = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_metric, 0.0, &ext_range, gFldExact);

//calculate the metrics with recovery method
gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 2, mapc2p, 0);
struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, ext_range.volume); //on ghosts
struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &ext_range, XYZ); //on ghosts with ext_range
gkyl_calc_metric *calculator = gkyl_calc_metric_new(&basis, &grid, false);
gkyl_calc_metric_advance( calculator, &range, XYZ, gFld);

gkyl_range_iter_init(&iter, &range);
while (gkyl_range_iter_next(&iter)) {
  long loc = gkyl_range_idx(&range, iter.idx);
  double *gij = gkyl_array_fetch(gFld, loc);
  double *gijExact = gkyl_array_fetch(gFldExact, loc);
  double val = *gij;
  double valExact = *gijExact;
  TEST_CHECK( gkyl_compare(val, valExact, 1e-12) );
}

gkyl_eval_on_nodes_release(eval_metric);
gkyl_eval_on_nodes_release(eval_mapc2p);
gkyl_array_release(XYZ);
gkyl_array_release(gFld);
gkyl_array_release(gFldExact);
gkyl_calc_metric_release(calculator);
}

TEST_LIST = {
  { "test_2x_p1", test_2x_p1},
  { "test_2x_p2", test_2x_p2},
  { NULL, NULL },
};
