#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_vlasov_alpha_gen_geo.h>

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
// 
void eval_conf_one(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}
void eval_conf_zero(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

//#define mkarr1(use_gpu, nc, size) (fprintf(stderr, "mkarr1: %d\n", __LINE__), mkarr1_(use_gpu, nc, size))

int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr);

void
test_alpha_gen_geo_(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 3, vdim = 3;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8, 8};
  int ghost[] = {0, 0, 0, 0, 0, 0};
  double lower[] = {0., 0., 0., -1., -1., -1.};
  double upper[] = {1., 1., 1., 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // projection updater for geo factors
  gkyl_proj_on_basis *projConfOne = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_conf_one, NULL);
  gkyl_proj_on_basis *projConfZero = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_conf_zero, NULL);

  // initialize arrays
  struct gkyl_array *a0, *a1;
  struct gkyl_array *tv_comp, *gij, *alpha_geo;
  
  a0 = mkarr1(use_gpu, confBasis.num_basis, confRange.volume);
  a1 = mkarr1(use_gpu, confBasis.num_basis, confRange.volume);
  tv_comp = mkarr1(use_gpu, 9*confBasis.num_basis, confRange.volume);
  gij = mkarr1(use_gpu, 6*confBasis.num_basis, confRange.volume);
  alpha_geo = mkarr1(use_gpu, 3*basis.num_basis, phaseRange.volume);

  // initialize arrays of 1 and 0.
  gkyl_proj_on_basis_advance(projConfOne, 0.0, &confRange, a1);
  gkyl_grid_sub_array_write(&confGrid, &confRange, a1, "ctest_alpha_geo_a1.gkyl");
  gkyl_array_clear(a0, 0.0);

  // Set other arrays to 0.
  gkyl_array_clear(tv_comp, 0.0);
  gkyl_array_clear(gij, 0.0);

  // Loop to initialize tv_comp:
  // dX/dx, dY/dx, dZ/dx,
  // dX/dy, dY/dy, dZ/dy,
  // dX/dz, dY/dz, dZ/dz,
  for(int i=0; i< 9; i++) {
    if (i%4 == 0.) gkyl_array_accumulate_offset(tv_comp, 1.0, a1, i*a1->ncomp);
  }
  gkyl_grid_sub_array_write(&confGrid, &confRange, tv_comp, "ctest_alpha_geo_tv_comp.gkyl");
  
  // gij is single array containing gxx, gxy, gxz, gyy, gyz, gzz
  gkyl_array_accumulate_offset(gij, 1.0, a1, 0*a1->ncomp);
  gkyl_array_accumulate_offset(gij, 1.0, a1, 3*a1->ncomp);
  gkyl_array_accumulate_offset(gij, 1.0, a1, 5*a1->ncomp);
  gkyl_grid_sub_array_write(&confGrid, &confRange, gij, "ctest_alpha_geo_gij.gkyl");
  
  // call alpha gen geo method
  gkyl_dg_alpha_gen_geo(&confBasis, &basis, &confRange, &phaseRange, &phaseGrid, tv_comp, gij, alpha_geo);
  // write out alpha_gen_geo
  gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, alpha_geo, "ctest_alpha_geo_3x3v.gkyl");
  
  /* // clean up */
  /* gkyl_array_release(fin); */
  /* gkyl_array_release(rhs); */
  /* gkyl_array_release(rhs_h); */
  /* gkyl_array_release(cflrate); */
  /* gkyl_array_release(qmem); */

  /* gkyl_hyper_dg_release(slvr); */
  /* gkyl_dg_eqn_release(eqn); */
}

void
test_alpha_gen_geo()
{
  test_alpha_gen_geo_(false);
}

TEST_LIST = {
  { "test_alpha_gen_geo", test_alpha_gen_geo },
  { NULL, NULL },
};

