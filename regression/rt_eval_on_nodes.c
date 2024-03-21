#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_array_rio.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

void evalFunc(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + sin(3*y)*cos(z) + z*z;
}

void
test_eval_1(enum gkyl_basis_type type)
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {8, 8, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  if (type == GKYL_BASIS_MODAL_SERENDIPITY)
    gkyl_cart_modal_serendip(&basis, 3, poly_order);
  else
    gkyl_cart_modal_tensor(&basis, 3, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evalDistf = gkyl_eval_on_nodes_new(&grid, &basis,1, evalFunc, 0);

  // create array range: no ghost-cells 
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // run updater
  gkyl_eval_on_nodes_advance(evalDistf, 0.0, &arr_range, distf);

  // construct file name and write data out
  const char *fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "%s-ser.gkyl",
    [GKYL_BASIS_MODAL_TENSOR] = "%s-ten.gkyl",
  };
  int sz = snprintf(0, 0, fmt[type], "rt_eval_on_nodes");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt[type], "rt_eval_on_nodes");
  
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, fileNm);
  
  gkyl_eval_on_nodes_release(evalDistf);
  gkyl_array_release(distf);
  
}

void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], th = xn[1];
  fout[0] = r*cos(th); fout[1] = r*sin(th);
}

void gaussian(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double xp[2], vth = 0.5, xc = 0.0, yc = 2.5;
  mapc2p(t, xn, xp, 0);
  double r2 = (xp[0]-xc)*(xp[0]-xc) + (xp[1]-yc)*(xp[1]-yc);
  fout[0] = exp(-r2/(2*vth*vth));
}

void
test_eval_2(enum gkyl_basis_type type)
{
  int poly_order = 2;
  double lower[] = {1.0, 0.0}, upper[] = {4.0, 2.0*M_PI};
  int cells[] = {8, 8*2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  if (type == GKYL_BASIS_MODAL_SERENDIPITY)
    gkyl_cart_modal_serendip(&basis, 2, poly_order);
  else
    gkyl_cart_modal_tensor(&basis, 2, poly_order);

  // create array range: no ghost-cells 
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);  

  // create DG expansion of mapping
  gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 2, mapc2p, 0);
  struct gkyl_array *rtheta = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
  gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &arr_range, rtheta);

  // ICs on mapped grid
  gkyl_proj_on_basis *proj_ic = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, gaussian, 0);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(proj_ic, 0.0, &arr_range, f);

  // construct file name and write data out
  const char *fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "%s-ser.gkyl",
    [GKYL_BASIS_MODAL_TENSOR] = "%s-ten.gkyl",
  };

  do {
    int sz = snprintf(0, 0, fmt[type], "rt_eval_on_nodes_rtheta");
    char fileNm[sz+1]; // ensures no buffer overflow  
    snprintf(fileNm, sizeof fileNm, fmt[type], "rt_eval_on_nodes_rtheta");
    gkyl_grid_sub_array_write(&grid, &arr_range, 0, rtheta, fileNm);
  } while (0);

  do {
    int sz = snprintf(0, 0, fmt[type], "rt_eval_on_nodes_f");
    char fileNm[sz+1]; // ensures no buffer overflow  
    snprintf(fileNm, sizeof fileNm, fmt[type], "rt_eval_on_nodes_f");
    gkyl_grid_sub_array_write(&grid, &arr_range, 0, f, fileNm);
  } while (0);
  
  gkyl_eval_on_nodes_release(eval_mapc2p);
  gkyl_proj_on_basis_release(proj_ic);
  gkyl_array_release(rtheta);
}

void shock(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0] > 2.5 ? 1.0 : 0.0;
}

void
test_eval_3(enum gkyl_basis_type type)
{
  int poly_order = 1;
  double lower[] = {1.0, 0.0}, upper[] = {4.0, 2.0*M_PI};
  int cells[] = {8, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  if (type == GKYL_BASIS_MODAL_SERENDIPITY)
    gkyl_cart_modal_serendip(&basis, 2, poly_order);
  else
    gkyl_cart_modal_tensor(&basis, 2, poly_order);

  // create array range: no ghost-cells 
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);  

  // create DG expansion of mapping
  gkyl_eval_on_nodes *eval_mapc2p = gkyl_eval_on_nodes_new(&grid, &basis, 2, mapc2p, 0);
  struct gkyl_array *rtheta = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);
  gkyl_eval_on_nodes_advance(eval_mapc2p, 0.0, &arr_range, rtheta);

  // ICs on mapped grid
  gkyl_proj_on_basis *proj_ic = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, shock, 0);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(proj_ic, 0.0, &arr_range, f);

  // construct file name and write data out
  const char *fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "%s-ser.gkyl",
    [GKYL_BASIS_MODAL_TENSOR] = "%s-ten.gkyl",
  };

  do {
    int sz = snprintf(0, 0, fmt[type], "rt_eval_on_nodes_rtheta_shock");
    char fileNm[sz+1]; // ensures no buffer overflow  
    snprintf(fileNm, sizeof fileNm, fmt[type], "rt_eval_on_nodes_rtheta_shock");
    gkyl_grid_sub_array_write(&grid, &arr_range, 0, rtheta, fileNm);
  } while (0);

  do {
    int sz = snprintf(0, 0, fmt[type], "rt_eval_on_nodes_f_shock");
    char fileNm[sz+1]; // ensures no buffer overflow  
    snprintf(fileNm, sizeof fileNm, fmt[type], "rt_eval_on_nodes_f_shock");
    gkyl_grid_sub_array_write(&grid, &arr_range, 0, f, fileNm);
  } while (0);
  
  gkyl_eval_on_nodes_release(eval_mapc2p);
  gkyl_proj_on_basis_release(proj_ic);
  gkyl_array_release(rtheta);
}

int
main(int argc, char **argv)
{
  test_eval_1(GKYL_BASIS_MODAL_SERENDIPITY);
  test_eval_1(GKYL_BASIS_MODAL_TENSOR);

  test_eval_2(GKYL_BASIS_MODAL_SERENDIPITY);
  test_eval_2(GKYL_BASIS_MODAL_TENSOR);

  test_eval_3(GKYL_BASIS_MODAL_SERENDIPITY);
  test_eval_3(GKYL_BASIS_MODAL_TENSOR);  

  return 0;
}
