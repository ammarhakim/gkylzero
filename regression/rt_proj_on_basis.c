#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_util.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <thpool.h>

void evalFunc_2d(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  //fout[0] = 0.1*(x-0.5)*(x-0.5)*(x-1.5)*sin(x) + 0.1*(y-1)*(y-1)*(y-1.5)*cos(y);
  fout[0] = 0.1*(x-0.5)*(x-0.5) + 0.1*(y-1)*(y-1)*(y-1.5);
}

void
test_proj_2d(enum gkyl_basis_type bt, enum gkyl_quad_type qt)
{
  int ndim = 2;
  int poly_order = 2;
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {4, 4};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  if (bt == GKYL_BASIS_MODAL_SERENDIPITY)
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  else
    gkyl_cart_modal_tensor(&basis, ndim, poly_order);

  const int num_quad[] = {
    [GKYL_QUAD_GAUSS_LEGENDRE] = poly_order+1,
    [GKYL_QUAD_GAUSS_LOBATTO] = poly_order+2
  };

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, qt,
    num_quad[qt], 1, evalFunc_2d, 0);

  // create array range: no ghost-cells 
  struct gkyl_range arr_range;
  gkyl_range_init_from_shape(&arr_range, ndim, cells);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // run updater
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // construct file name and write data out
  const char *b_fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "ser",
    [GKYL_BASIS_MODAL_TENSOR] = "ten",
  };
  const char *q_fmt[] = {
    [GKYL_QUAD_GAUSS_LEGENDRE] = "leg",
    [GKYL_QUAD_GAUSS_LOBATTO] = "lob",
  };  
  int sz = snprintf(0, 0, "rt_proj_on_basis-%s-%s_2d.gkyl", b_fmt[bt], q_fmt[qt]);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, "rt_proj_on_basis-%s-%s_2d.gkyl", b_fmt[bt], q_fmt[qt]);
  
  gkyl_grid_sub_array_write(&grid, &arr_range, distf, fileNm);
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  
}

void evalFunc_3d(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + sin(3*y)*cos(z) + z*z;
}

void
test_proj_3d(enum gkyl_basis_type bt, enum gkyl_quad_type qt)
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {8, 8, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  if (bt == GKYL_BASIS_MODAL_SERENDIPITY)
    gkyl_cart_modal_serendip(&basis, 3, poly_order);
  else
    gkyl_cart_modal_tensor(&basis, 3, poly_order);

  const int num_quad[] = {
    [GKYL_QUAD_GAUSS_LEGENDRE] = poly_order+1,
    [GKYL_QUAD_GAUSS_LOBATTO] = poly_order+2
  };

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, qt,
    num_quad[qt], 1, evalFunc_3d, 0);

  // create array range: no ghost-cells 
  struct gkyl_range arr_range;
  gkyl_range_init_from_shape(&arr_range, 3, cells);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // run updater
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // construct file name and write data out
  const char *b_fmt[] = {
    [GKYL_BASIS_MODAL_SERENDIPITY] = "ser",
    [GKYL_BASIS_MODAL_TENSOR] = "ten",
  };
  const char *q_fmt[] = {
    [GKYL_QUAD_GAUSS_LEGENDRE] = "leg",
    [GKYL_QUAD_GAUSS_LOBATTO] = "lob",
  };  
  int sz = snprintf(0, 0, "rt_proj_on_basis-%s-%s_3d.gkyl", b_fmt[bt], q_fmt[qt]);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, "rt_proj_on_basis-%s-%s_3d.gkyl", b_fmt[bt], q_fmt[qt]);
  
  gkyl_grid_sub_array_write(&grid, &arr_range, distf, fileNm);
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  
}

int
main(int argc, char **argv)
{
  // 2D tests
  test_proj_2d(GKYL_BASIS_MODAL_SERENDIPITY, GKYL_QUAD_GAUSS_LEGENDRE);
  test_proj_2d(GKYL_BASIS_MODAL_TENSOR, GKYL_QUAD_GAUSS_LEGENDRE);

  test_proj_2d(GKYL_BASIS_MODAL_SERENDIPITY, GKYL_QUAD_GAUSS_LOBATTO);
  test_proj_2d(GKYL_BASIS_MODAL_TENSOR, GKYL_QUAD_GAUSS_LOBATTO);
  
  // 3D tests
  test_proj_3d(GKYL_BASIS_MODAL_SERENDIPITY, GKYL_QUAD_GAUSS_LEGENDRE);
  test_proj_3d(GKYL_BASIS_MODAL_TENSOR, GKYL_QUAD_GAUSS_LEGENDRE);

  test_proj_3d(GKYL_BASIS_MODAL_SERENDIPITY, GKYL_QUAD_GAUSS_LOBATTO);
  test_proj_3d(GKYL_BASIS_MODAL_TENSOR, GKYL_QUAD_GAUSS_LOBATTO);

  return 0;
}
