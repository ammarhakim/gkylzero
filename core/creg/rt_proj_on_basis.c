#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <thpool.h>

void evalFunc(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + sin(3*y)*cos(z) + z*z;
}

void
test_proj_1(enum gkyl_basis_type type)
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
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, 0);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // run updater
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // construct file name and write data out
  const char *fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "%s-ser.gkyl",
    [GKYL_BASIS_MODAL_TENSOR] = "%s-ten.gkyl",
  };
  int sz = snprintf(0, 0, fmt[type], "rt_proj_on_basis");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt[type], "rt_proj_on_basis");
  
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, fileNm);
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  
}

void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], th = xn[1];
  fout[0] = r*cos(th); fout[1] = r*sin(th);
}

void
test_proj_2(enum gkyl_basis_type type)
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

  // projection updater for dist-function
  gkyl_proj_on_basis *proj_mapc2p = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, mapc2p, 0);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create DG expansion of mapping
  struct gkyl_array *rtheta = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, arr_range.volume);

  // run updater
  gkyl_proj_on_basis_advance(proj_mapc2p, 0.0, &arr_range, rtheta);

  // construct file name and write data out
  const char *fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "%s-ser.gkyl",
    [GKYL_BASIS_MODAL_TENSOR] = "%s-ten.gkyl",
  };
  int sz = snprintf(0, 0, fmt[type], "rt_proj_on_basis_rtheta");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt[type], "rt_proj_on_basis_rtheta");
  
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, rtheta, fileNm);
  
  gkyl_proj_on_basis_release(proj_mapc2p);
  gkyl_array_release(rtheta);
}

int
main(int argc, char **argv)
{
  test_proj_1(GKYL_BASIS_MODAL_SERENDIPITY);
  test_proj_1(GKYL_BASIS_MODAL_TENSOR);

  test_proj_2(GKYL_BASIS_MODAL_SERENDIPITY);
  test_proj_2(GKYL_BASIS_MODAL_TENSOR);

  return 0;
}
