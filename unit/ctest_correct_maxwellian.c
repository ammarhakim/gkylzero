#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

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

void eval_M0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.25;
}

void eval_M1i_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_M2_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double n = 1.0, vth2 = 1.0, ux = 0.5;
  double x = xn[0];
  fout[0] = n*vth2 + n*ux*ux;
}

void eval_udrift_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_vtsq_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void eval_M1i_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void eval_M2_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double n = 1.0, vth2 = 1.0, ux = 0.5, uy = 0.25;
  double x = xn[0];
  fout[0] = 2*n*vth2 + n*(ux*ux+uy*uy);
}

void eval_udrift_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void eval_vtsq_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void
test_1x1v(int poly_order, bool use_gpu)
{
  double lower[] = {0.1, -6.0}, upper[] = {1.0, 6.0};
  int cells[] = {2, 32};
  int vdim = 1, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0_cu, *m1i_cu, *m2_cu;

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_1v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute Maxwellian
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, use_gpu);
  // correction updater
  gkyl_correct_maxwellian *corr_max = gkyl_correct_maxwellian_new(&grid, &confBasis, &basis,
    confLocal.volume, confLocal_ext.volume);
  
  // project the Maxwellian
  gkyl_proj_maxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, m0, m1i, m2, distf);

 // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_correct_maxwellian_test_1x1v_p%d_uc.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  struct gkyl_array *m0_r;
  m0_r = mkarr(confBasis.num_basis, confLocal_ext.volume);  
  gkyl_array_scale(gkyl_array_copy(m0_r, m0), 2.5);
  
  // correct the Maxwellian
  gkyl_correct_maxwellian_fix(corr_max, distf, m0_r, &local, &confLocal);

 // write distribution function to file
  sprintf(fname, "ctest_correct_maxwellian_test_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // compute the number density
  struct gkyl_mom_type *m0_t = gkyl_mom_vlasov_new(&confBasis, &basis, "M0");
  struct gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, m0_t);
  gkyl_mom_type_release(m0_t);

  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0); // m0 = 2.5*orginal m0

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *n0 = gkyl_array_cfetch(m0, gkyl_range_idx(&confLocal, iter.idx));
    const double *nr = gkyl_array_cfetch(m0_r, gkyl_range_idx(&confLocal, iter.idx));
    
    for (int k=0; k<confBasis.num_basis; ++k)
      TEST_CHECK( gkyl_compare_double(n0[k], nr[k], 1e-14) );
  }
  
  gkyl_array_release(m0); gkyl_array_release(m0_r);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);

  gkyl_mom_calc_release(m0calc);
  gkyl_array_release(distf);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_correct_maxwellian_release(corr_max);
}

void test_1x1v_p1() { test_1x1v(1, false); }
void test_1x1v_p2() { test_1x1v(2, false); }

TEST_LIST = {
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x1v_p2", test_1x1v_p2 },
  { NULL, NULL },
};
