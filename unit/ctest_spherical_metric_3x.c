#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_calc_metric.h>
#include <gkyl_calc_metric_kernels.h>

// Helper Functions

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

// Apply periodic BCs along direction
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
}




void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  fout[0] = r*sin(theta)*cos(phi); fout[1] = r*sin(theta)*sin(phi); fout[2] = r*cos(theta);
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = r*r;
  fout[4] = 0.0;
  fout[5] = r*r*sin(theta)*sin(theta);
}

void
test_3x_p1()
{
struct gkyl_basis basis;
int poly_order = 1;
gkyl_cart_modal_serendip(&basis, 3, poly_order);

double lower[3] = { 0.5, -M_PI/2, 0 }, upper[3] = { 1.0, M_PI/2, M_PI/2 };
int cells[3] = { 8, 8, 8 };
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

int bcs[3] = {1,1,1};
// Need to pass appropriate range to eval mapc2p depending on bcsbcs
struct gkyl_range conversion_range;
if(bcs[1]==1){
  int sublower[3] = {range.lower[0], ext_range.lower[1], range.lower[2]};
  int subupper[3] = {range.upper[0], ext_range.upper[1], range.upper[2]};
  gkyl_sub_range_init(&conversion_range, &ext_range, sublower, subupper);
}
else{
  conversion_range = range;
  //int sublower[3] = {range.lower[0], range.lower[1], range.lower[2]};
  //int subupper[3] = {range.upper[0], range.upper[1], range.upper[2]};
  //gkyl_sub_range_init(&conversion_range, &ext_range, sublower, subupper);
}


//calculate the metrics with recovery method
gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 3, mapc2p, 0);
struct gkyl_array *XYZ = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, ext_range.volume);
struct gkyl_array *gFld = gkyl_array_new(GKYL_DOUBLE, 6*basis.num_basis, ext_range.volume);
gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &ext_range, XYZ);

// sync in periodic dirs (y)
struct skin_ghost_ranges skin_ghost; // skin/ghost.
skin_ghost_ranges_init(&skin_ghost, &ext_range, nghost);
for(int i = 0; i<3; i++){
  if(bcs[i]==0){
    struct gkyl_array *perbuff = mkarr(3*basis.num_basis, skin_ghost.lower_skin[i].volume);
    apply_periodic_bc(perbuff, XYZ, i, skin_ghost);
  }
}

gkyl_calc_metric *calculator = gkyl_calc_metric_new(&basis, &grid, bcs, false);
gkyl_calc_metric_advance( calculator, &range, XYZ, gFld);



gkyl_grid_sub_array_write(&grid, &range, gFld, "spherical_gFld.gkyl");


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
