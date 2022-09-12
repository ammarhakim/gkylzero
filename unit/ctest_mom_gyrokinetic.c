// Test calculation of gyrokinetic moments of a distribution function.
//
#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_array_rio.h>
#include <math.h>

void
test_mom_gyrokinetic()
{
  double mass = 1.0;
  int poly_order = 2;
  double lower[] = {-M_PI, -2.0, 0.0}, upper[] = {M_PI, 2.0, 2.0};
  int cells[] = {4, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 1;

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

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  struct gkyl_mom_type *m2 = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", false);

  TEST_CHECK( m2->cdim == 1 );
  TEST_CHECK( m2->pdim == 3 );
  TEST_CHECK( m2->poly_order == 2 );
  TEST_CHECK( m2->num_config == confBasis.num_basis );
  TEST_CHECK( m2->num_phase == basis.num_basis );
  TEST_CHECK( m2->num_mom == 1 );

  struct gkyl_mom_type *m3par = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M3par", false);
  TEST_CHECK( m3par->num_mom == 1 );

  gkyl_mom_type_release(m2);
  gkyl_mom_type_release(m3par);
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
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

void bmag_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = cos((2.*M_PI/(2.*2.*M_PI))*x);
}
void bmag_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = cos((2.*M_PI/(2.*2.*M_PI))*x)*exp(-(y*y)/(2.*pow(M_PI/3,2)));
}

void distf_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1];
  double bmag[1];
  bmag_1x(t, xn, &bmag[0], ctx); 
  fout[0] = bmag[0]*(x*x)*(vpar-0.5)*(vpar-0.5);
}
void distf_1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1], mu = xn[2];
  double bmag[1];
  bmag_1x(t, xn, &bmag[0], ctx); 
  fout[0] = bmag[0]*(x*x)*(vpar-0.5)*(vpar-0.5);
}
void distf_2x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], vpar = xn[2], mu = xn[3];
  double bmag[1];
  bmag_2x(t, xn, &bmag[0], ctx); 
  fout[0] = bmag[0]*(x*x+y*y)*(vpar-0.5)*(vpar-0.5);
}

void
test_1x1v(int polyOrder)
{
  double mass = 1.;
  int poly_order = 1;
  double lower[] = {-M_PI, -2.0}, upper[] = {M_PI, 2.0};
  int cells[] = {4, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
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

  // create distribution function array and project distribution function on basis
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, distf_1x1v, NULL);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
//  gkyl_grid_sub_array_write(&grid, &local, distf, "ctest_mom_gyrokinetic_1x1v_p1_distf.gkyl");

  // create bmag array and project magnetic field amplitude function on basis
  struct gkyl_array *bmag;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *projbmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, bmag_1x, NULL);
  gkyl_proj_on_basis_advance(projbmag, 0.0, &confLocal, bmag);
//  gkyl_grid_sub_array_write(&confGrid, &confLocal, bmag, "ctest_mom_gyrokinetic_1x1v_p1_bmag.gkyl");

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", false);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", false);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", false);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, false);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, false);

  // create moment arrays
  struct gkyl_array *m0, *m1, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, distf, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);

  double *m00 = gkyl_array_fetch(m0, 0+confGhost[0]); double *m01 = gkyl_array_fetch(m0, 1+confGhost[0]);
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  double *m10 = gkyl_array_fetch(m1, 0+confGhost[0]); double *m11 = gkyl_array_fetch(m1, 1+confGhost[0]);
  double *m12 = gkyl_array_fetch(m1, 2+confGhost[0]); double *m13 = gkyl_array_fetch(m1, 3+confGhost[0]);
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);

  if (poly_order==1) {
    // Check M0.
    TEST_CHECK( gkyl_compare(  1.52537436025689e+01, m00[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  3.57234787161818e+00, m00[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  6.08286277145113e+00, m01[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -5.10949091314880e+00, m01[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  6.08286277145113e+00, m02[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  5.10949091314880e+00, m02[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.52537436025689e+01, m03[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -3.57234787161818e+00, m03[1], 1e-12) );

    // Check M1.
    TEST_CHECK( gkyl_compare( -1.28452577705844e+01, m10[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -3.00829294452057e+00, m10[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -5.12241075490621e+00, m11[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  4.30272919002004e+00, m11[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -5.12241075490621e+00, m12[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -4.30272919002004e+00, m12[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -1.28452577705844e+01, m13[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  3.00829294452057e+00, m13[1], 1e-12) );

    // Check M2.
    TEST_CHECK( gkyl_compare(  3.31835825740096e+01, m20[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  7.77142344001148e+00, m20[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.32328944501744e+01, m21[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -1.11153837408851e+01, m21[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.32328944501744e+01, m22[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.11153837408851e+01, m22[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  3.31835825740096e+01, m23[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -7.77142344001148e+00, m23[1], 1e-12) );
  } else if (poly_order==2) {
    double m0Correct[] = {
       1.526837339934706e+01,   3.951518219554417e+00,  -3.363344534446567e+00,
       6.052510010088350e+00,  -4.868229034295940e+00,   8.480389975048731e-01,
       6.052510010088350e+00,   4.868229034295939e+00,   8.480389975048728e-01,
       1.526837339934706e+01,  -3.951518219554418e+00,  -3.363344534446568e+00
    };
    double m1Correct[] = {
       -1.285757759945016e+01,  -3.327594290151089e+00,   2.832290134270792e+00,
       -5.096850534811242e+00,   4.099561292038686e+00,  -7.141381031619991e-01,
       -5.096850534811242e+00,  -4.099561292038685e+00,  -7.141381031619988e-01,
       -1.285757759945016e+01,   3.327594290151089e+00,   2.832290134270792e+00
    };
    double m2Correct[] = {
       3.407258063854292e+01,   8.818124868900387e+00,  -7.316749513532877e+00,
       1.350665391724979e+01,  -1.086383742390252e+01,   1.844856766501832e+00,
       1.350665391724979e+01,   1.086383742390252e+01,   1.844856766501832e+00,
       3.407258063854292e+01,  -8.818124868900387e+00,  -7.316749513532877e+00
    };
    // Check M0.
    TEST_CHECK( gkyl_compare( m0Correct[0 ],  m00[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[1 ],  m00[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[2 ],  m00[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[3 ],  m01[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[4 ],  m01[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[5 ],  m01[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[6 ],  m02[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[7 ],  m02[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[8 ],  m02[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[9 ],  m03[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[10], m03[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m0Correct[11], m03[2], 1e-12) );

    // Check M1.
    TEST_CHECK( gkyl_compare( m1Correct[0 ], m10[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[1 ], m10[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[2 ], m10[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[3 ], m11[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[4 ], m11[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[5 ], m11[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[6 ], m12[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[7 ], m12[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[8 ], m12[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[9 ], m13[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[10], m13[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m1Correct[11], m13[2], 1e-12) );

    // Check M2.
    TEST_CHECK( gkyl_compare( m2Correct[0 ], m20[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[1 ], m20[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[2 ], m20[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[3 ], m21[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[4 ], m21[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[5 ], m21[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[6 ], m22[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[7 ], m22[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[8 ], m22[2], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[9 ], m23[0], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[10], m23[1], 1e-12) );
    TEST_CHECK( gkyl_compare( m2Correct[11], m23[2], 1e-12) );
  }

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1calc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(M0_t); gkyl_mom_type_release(M1_t); gkyl_mom_type_release(M2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_1x2v(int poly_order)
{
  double mass = 1.;
  double lower[] = {-M_PI, -2.0, 0.0}, upper[] = {M_PI, 2.0, 2.0};
  int cells[] = {4, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create distribution function array and project distribution function on basis
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, distf_1x2v, NULL);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // create bmag array and project magnetic field amplitude function on basis
  struct gkyl_array *bmag;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *projbmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, bmag_1x, NULL);
  gkyl_proj_on_basis_advance(projbmag, 0.0, &confLocal, bmag);

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", false);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", false);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", false);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, false);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, false);

  // create moment arrays
  struct gkyl_array *m0, *m1, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, distf, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);

  double *m00 = gkyl_array_fetch(m0, gkyl_range_idx(&confLocal, &(int) {0+confGhost[0]})); 
  double *m01 = gkyl_array_fetch(m0, gkyl_range_idx(&confLocal, &(int) {1+confGhost[0]}));
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  double *m10 = gkyl_array_fetch(m1, 0+confGhost[0]); double *m11 = gkyl_array_fetch(m1, 1+confGhost[0]);
  double *m12 = gkyl_array_fetch(m1, 2+confGhost[0]); double *m13 = gkyl_array_fetch(m1, 3+confGhost[0]);
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);
  if (poly_order==1) {
    // Check M0.
    TEST_CHECK( gkyl_compare(  191.6841953662915, m00[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  44.89144731817122, m00[1], 1e-12) );
    TEST_CHECK( gkyl_compare(   76.4395079823428, m01[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  -64.2077564653283, m01[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  76.43950798234282, m02[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  64.20775646532829, m02[1], 1e-12) );
    TEST_CHECK( gkyl_compare( 191.68419536629153, m03[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -44.89144731817124, m03[1], 1e-12) );
  
    // Check M1.
    TEST_CHECK( gkyl_compare( -161.41826978214021, m10[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  -37.80332405740736, m10[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  -64.37011198513079, m11[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  54.069689655013306, m11[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  -64.37011198513079, m12[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  -54.06968965501329, m12[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -161.41826978214021, m13[0], 1e-12) );
    TEST_CHECK( gkyl_compare(   37.80332405740738, m13[1], 1e-12) );
  
    // Check M2.
    TEST_CHECK( gkyl_compare(   578.5971330556217, m20[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  210.75432886828457, m20[1], 1e-12) );
    TEST_CHECK( gkyl_compare(   292.8699479183419, m21[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  -242.1332009483718, m21[1], 1e-12) );
    TEST_CHECK( gkyl_compare(   292.8699479183419, m22[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  242.13320094837178, m22[1], 1e-12) );
    TEST_CHECK( gkyl_compare(   578.5971330556217, m23[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -210.75432886828463, m23[1], 1e-12) );
  } else if (poly_order==2) {
    // Check M0.
    TEST_CHECK( gkyl_compare(  1.918680388146181e+02, m00[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  4.965624243631351e+01, m00[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -4.226503392363477e+01, m00[2], 1e-12) );
    TEST_CHECK( gkyl_compare(  7.605808393388897e+01, m01[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -6.117597028054660e+01, m01[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.065677233807588e+01, m01[2], 1e-12) );
    TEST_CHECK( gkyl_compare(  7.605808393388898e+01, m02[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  6.117597028054661e+01, m02[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.065677233807589e+01, m02[2], 1e-12) );
    TEST_CHECK( gkyl_compare(  1.918680388146182e+02, m03[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -4.965624243631353e+01, m03[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -4.226503392363478e+01, m03[2], 1e-12) );
  
    // Check M1.
    TEST_CHECK( gkyl_compare( -1.615730853175732e+02, m10[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -4.181578310426401e+01, m10[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  3.559160751463980e+01, m10[2], 1e-12) );
    TEST_CHECK( gkyl_compare( -6.404891278643279e+01, m11[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  5.151660655203925e+01, m11[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -8.974124074169159e+00, m11[2], 1e-12) );
    TEST_CHECK( gkyl_compare( -6.404891278643282e+01, m12[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -5.151660655203925e+01, m12[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -8.974124074169165e+00, m12[2], 1e-12) );
    TEST_CHECK( gkyl_compare( -1.615730853175732e+02, m13[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  4.181578310426401e+01, m13[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  3.559160751463981e+01, m13[2], 1e-12) );
  
    // Check M2.
    TEST_CHECK( gkyl_compare(  5.924940346248225e+02, m20[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  2.106245140015062e+02, m20[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -1.080258558862571e+02, m20[2], 1e-09) );
    TEST_CHECK( gkyl_compare(  2.957803339447422e+02, m21[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -2.297438543767285e+02, m21[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  2.952994007812927e+01, m21[2], 1e-12) );
    TEST_CHECK( gkyl_compare(  2.957803339447422e+02, m22[0], 1e-12) );
    TEST_CHECK( gkyl_compare(  2.297438543767285e+02, m22[1], 1e-12) );
    TEST_CHECK( gkyl_compare(  2.952994007812929e+01, m22[2], 1e-12) );
    TEST_CHECK( gkyl_compare(  5.924940346248227e+02, m23[0], 1e-12) );
    TEST_CHECK( gkyl_compare( -2.106245140015063e+02, m23[1], 1e-12) );
    TEST_CHECK( gkyl_compare( -1.080258558862570e+02, m23[2], 1e-12) );
  }

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1calc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(M0_t); gkyl_mom_type_release(M1_t); gkyl_mom_type_release(M2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_2x2v(int poly_order)
{
  double mass = 1.;
  double lower[] = {-M_PI, -M_PI, -2.0, 0.0}, upper[] = {M_PI, M_PI, 2.0, 2.0};
  int cells[] = {4, 4, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 2, cdim = 2;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // create distribution function array and project distribution function on basis
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, distf_2x2v, NULL);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // create bmag array and project magnetic field amplitude function on basis
  struct gkyl_array *bmag;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *projbmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, bmag_2x, NULL);
  gkyl_proj_on_basis_advance(projbmag, 0.0, &confLocal, bmag);

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", false);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", false);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", false);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, false);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, false);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, false);

  // create moment arrays
  struct gkyl_array *m0, *m1, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, distf, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);

  if (poly_order==1) {
    double m0Correct[] = {
       5.674373976691270e+01,  2.201144875962325e+01,  3.651966323006221e+01,  1.314016650990054e+01,
       2.219569636062632e+02,  6.028640616624954e+01,  4.206235100555175e+01,  3.552942178215428e+00,
       2.219569636062633e+02,  6.028640616624952e+01, -4.206235100555176e+01, -3.552942178215411e+00,
       5.674373976691270e+01,  2.201144875962326e+01, -3.651966323006221e+01, -1.314016650990054e+01,
       7.709667278852976e+01, -3.719961743450785e+00,  4.321049456829391e+01, -4.192642488321794e+00,
       1.403753088009537e+02, -5.979190987388544e+01, -2.255457479022396e+01, -2.512741619288884e+01,
       1.403753088009537e+02, -5.979190987388544e+01,  2.255457479022397e+01,  2.512741619288884e+01,
       7.709667278852976e+01, -3.719961743450784e+00, -4.321049456829390e+01,  4.192642488321796e+00,
       7.709667278852976e+01,  3.719961743450789e+00,  4.321049456829391e+01,  4.192642488321798e+00,
       1.403753088009537e+02,  5.979190987388544e+01, -2.255457479022396e+01,  2.512741619288883e+01,
       1.403753088009538e+02,  5.979190987388544e+01,  2.255457479022397e+01, -2.512741619288886e+01,
       7.709667278852974e+01,  3.719961743450781e+00, -4.321049456829391e+01, -4.192642488321795e+00,
       5.674373976691270e+01, -2.201144875962325e+01,  3.651966323006221e+01, -1.314016650990054e+01,
       2.219569636062632e+02, -6.028640616624953e+01,  4.206235100555175e+01, -3.552942178215412e+00,
       2.219569636062632e+02, -6.028640616624953e+01, -4.206235100555175e+01,  3.552942178215411e+00,
       5.674373976691268e+01, -2.201144875962326e+01, -3.651966323006221e+01,  1.314016650990054e+01
    };
    double m1Correct[] = {
      -4.778420190897913e+01, -1.853595685020905e+01, -3.075340061478923e+01, -1.106540337675835e+01,
      -1.869111272473795e+02, -5.076749992947330e+01, -3.542092716256990e+01, -2.991951307970891e+00,
      -1.869111272473796e+02, -5.076749992947327e+01,  3.542092716256991e+01,  2.991951307970870e+00,
      -4.778420190897913e+01, -1.853595685020905e+01,  3.075340061478923e+01,  1.106540337675835e+01,
      -6.492351392718295e+01,  3.132599362905920e+00, -3.638778489961594e+01,  3.530646305955197e+00,
      -1.182107863586978e+02,  5.035108199906144e+01,  1.899332613913597e+01,  2.115992942559060e+01,
      -1.182107863586978e+02,  5.035108199906144e+01, -1.899332613913598e+01, -2.115992942559060e+01,
      -6.492351392718297e+01,  3.132599362905925e+00,  3.638778489961592e+01, -3.530646305955196e+00,
      -6.492351392718297e+01, -3.132599362905928e+00, -3.638778489961594e+01, -3.530646305955196e+00,
      -1.182107863586979e+02, -5.035108199906144e+01,  1.899332613913597e+01, -2.115992942559060e+01,
      -1.182107863586979e+02, -5.035108199906144e+01, -1.899332613913597e+01,  2.115992942559062e+01,
      -6.492351392718295e+01, -3.132599362905919e+00,  3.638778489961592e+01,  3.530646305955196e+00,
      -4.778420190897911e+01,  1.853595685020907e+01, -3.075340061478922e+01,  1.106540337675835e+01,
      -1.869111272473795e+02,  5.076749992947325e+01, -3.542092716256990e+01,  2.991951307970878e+00,
      -1.869111272473796e+02,  5.076749992947332e+01,  3.542092716256988e+01, -2.991951307970872e+00,
      -4.778420190897911e+01,  1.853595685020906e+01,  3.075340061478923e+01, -1.106540337675835e+01
    };
    double m2Correct[] = {
       1.317741859760154e+02,  5.432242387904527e+01,  8.726478357864879e+01,  3.461290078036008e+01,
       6.282598132762612e+02,  2.349967248260282e+02,  1.585668164431439e+02,  5.344698656702761e+01,
       6.282598132762612e+02,  2.349967248260282e+02, -1.585668164431439e+02, -5.344698656702760e+01,
       1.317741859760154e+02,  5.432242387904527e+01, -8.726478357864879e+01, -3.461290078036008e+01,
       1.892086113026537e+02, -7.382118864615809e+00,  1.138208348715406e+02, -8.593088099826790e+00,
       4.706732142611326e+02, -2.016751811606539e+02, -2.651705819234549e+01, -1.091874194410921e+02,
       4.706732142611326e+02, -2.016751811606539e+02,  2.651705819234548e+01,  1.091874194410921e+02,
       1.892086113026537e+02, -7.382118864615817e+00, -1.138208348715406e+02,  8.593088099826790e+00,
       1.892086113026537e+02,  7.382118864615816e+00,  1.138208348715406e+02,  8.593088099826790e+00,
       4.706732142611326e+02,  2.016751811606539e+02, -2.651705819234548e+01,  1.091874194410921e+02,
       4.706732142611327e+02,  2.016751811606539e+02,  2.651705819234547e+01, -1.091874194410922e+02,
       1.892086113026537e+02,  7.382118864615806e+00, -1.138208348715406e+02, -8.593088099826787e+00,
       1.317741859760154e+02, -5.432242387904528e+01,  8.726478357864879e+01, -3.461290078036008e+01,
       6.282598132762611e+02, -2.349967248260282e+02,  1.585668164431439e+02, -5.344698656702759e+01,
       6.282598132762612e+02, -2.349967248260283e+02, -1.585668164431438e+02,  5.344698656702758e+01,
       1.317741859760154e+02, -5.432242387904529e+01, -8.726478357864879e+01,  3.461290078036008e+01
    };
    for (int i=0; i<cells[0]; i++) {
      for (int j=0; j<cells[1]; j++) {
        int idx[] = {i+confGhost[0], j+confGhost[1]};
        // Check M0.
        double *m0ptr = gkyl_array_fetch(m0, gkyl_range_idx(&confLocal, idx)); 
        for (int k=0; k<confBasis.num_basis; k++)
          TEST_CHECK( gkyl_compare( m0Correct[(i*cells[1]+j)*confBasis.num_basis+k], m0ptr[k], 1e-12) );
        // Check M1.
        double *m1ptr = gkyl_array_fetch(m1, gkyl_range_idx(&confLocal, idx)); 
        for (int k=0; k<confBasis.num_basis; k++)
          TEST_CHECK( gkyl_compare( m1Correct[(i*cells[1]+j)*confBasis.num_basis+k], m1ptr[k], 1e-12) );
        // Check M2.
        double *m2ptr = gkyl_array_fetch(m2, gkyl_range_idx(&confLocal, idx)); 
        for (int k=0; k<confBasis.num_basis; k++)
          TEST_CHECK( gkyl_compare( m2Correct[(i*cells[1]+j)*confBasis.num_basis+k], m2ptr[k], 1e-12) );
      }
    }
  
  } else if (poly_order==2) {
  
    double m0Correct[] = {
      5.648341168299724e+01, 2.250387513900885e+01, 3.652502205171847e+01, 1.335144639473130e+01, -7.166399139585262e+00, 7.761081311319963e+00,
      -5.437466693389656e+00, 2.161243931947423e+00, 2.224236975980227e+02, 6.521607861574371e+01, 4.094337722353598e+01, 5.582134934124797e+00,
      -4.387646236656694e+01, -8.976042695293010e+00, -1.237376484382488e+01, -2.753498001080819e+00, 2.224236975980228e+02, 6.521607861574371e+01,
      -4.094337722353598e+01, -5.582134934124794e+00, -4.387646236656691e+01, -8.976042695293000e+00, 1.237376484382486e+01, -2.753498001080819e+00,
      5.648341168299725e+01, 2.250387513900884e+01, -3.652502205171847e+01, -1.335144639473129e+01, -7.166399139585268e+00, 7.761081311319959e+00,
      5.437466693389650e+00, 2.161243931947424e+00, 7.646435675172015e+01, -3.325734694746438e+00, 4.121233955674857e+01, -4.386344478763552e+00,
      1.342837549799004e-01, 4.123263549958735e+00, 5.440652066518922e-01, -2.190342820971765e+00, 1.406410562141511e+02, -5.667019945920631e+01,
      -1.815281995189654e+01, -2.239122325134061e+01, 9.440015133641253e+00, -6.509924435714448e+00, 4.183511326759515e+00, 2.060415698733552e+00,
      1.406410562141512e+02, -5.667019945920630e+01, 1.815281995189653e+01, 2.239122325134061e+01, 9.440015133641275e+00, -6.509924435714441e+00,
      -4.183511326759525e+00, 2.060415698733553e+00, 7.646435675172015e+01, -3.325734694746449e+00, -4.121233955674858e+01, 4.386344478763548e+00,
      1.342837549798988e-01, 4.123263549958738e+00, -5.440652066518882e-01, -2.190342820971765e+00, 7.646435675172016e+01, 3.325734694746450e+00,
      4.121233955674858e+01, 4.386344478763550e+00, 1.342837549798992e-01, 4.123263549958738e+00, 5.440652066518876e-01, 2.190342820971762e+00,
      1.406410562141513e+02, 5.667019945920629e+01, -1.815281995189653e+01, 2.239122325134060e+01, 9.440015133641266e+00, -6.509924435714446e+00,
      4.183511326759531e+00, -2.060415698733550e+00, 1.406410562141513e+02, 5.667019945920630e+01, 1.815281995189654e+01, -2.239122325134062e+01,
      9.440015133641259e+00, -6.509924435714444e+00, -4.183511326759524e+00, -2.060415698733558e+00, 7.646435675172017e+01, 3.325734694746441e+00,
      -4.121233955674857e+01, -4.386344478763546e+00, 1.342837549798972e-01, 4.123263549958737e+00, -5.440652066518941e-01, 2.190342820971768e+00,
      5.648341168299724e+01, -2.250387513900884e+01, 3.652502205171847e+01, -1.335144639473129e+01, -7.166399139585268e+00, 7.761081311319965e+00,
      -5.437466693389656e+00, -2.161243931947421e+00, 2.224236975980227e+02, -6.521607861574371e+01, 4.094337722353598e+01, -5.582134934124787e+00,
      -4.387646236656693e+01, -8.976042695293009e+00, -1.237376484382487e+01, 2.753498001080815e+00, 2.224236975980228e+02, -6.521607861574373e+01,
      -4.094337722353598e+01, 5.582134934124794e+00, -4.387646236656688e+01, -8.976042695292994e+00, 1.237376484382487e+01, 2.753498001080822e+00,
      5.648341168299726e+01, -2.250387513900885e+01, -3.652502205171848e+01, 1.335144639473129e+01, -7.166399139585263e+00, 7.761081311319963e+00,
      5.437466693389656e+00, -2.161243931947425e+00,
    };
    double m1Correct[] = {
      -4.756497825936609e+01, -1.895063169600745e+01, -3.075791330671029e+01, -1.124332327977373e+01, 6.034862433334959e+00, -6.535647420058911e+00,
      4.578919320749178e+00, -1.819994890060987e+00, -1.873041663983350e+02, -5.491880304483679e+01, -3.447863345139874e+01, -4.700745207684037e+00,
      3.694859988763532e+01, 7.558772796036234e+00, 1.042001250006305e+01, 2.318735158804893e+00, -1.873041663983350e+02, -5.491880304483681e+01,
      3.447863345139872e+01, 4.700745207684039e+00, 3.694859988763526e+01, 7.558772796036222e+00, -1.042001250006304e+01, 2.318735158804894e+00,
      -4.756497825936611e+01, -1.895063169600745e+01, 3.075791330671029e+01, 1.124332327977372e+01, 6.034862433334967e+00, -6.535647420058911e+00,
      -4.578919320749177e+00, -1.819994890060993e+00, -6.439103726460642e+01, 2.800618690312792e+00, -3.470512804778827e+01, 3.693763771590361e+00,
      -1.130810568251719e-01, -3.472221936807352e+00, -4.581601740226504e-01, 1.844499217660437e+00, -1.184345736540220e+02, 4.772227322880535e+01,
      1.528658522264971e+01, 1.885576694849735e+01, -7.949486428329483e+00, 5.482041630075337e+00, -3.522956906744856e+00, -1.735086904196671e+00,
      -1.184345736540220e+02, 4.772227322880529e+01, -1.528658522264972e+01, -1.885576694849735e+01, -7.949486428329486e+00, 5.482041630075328e+00,
      3.522956906744861e+00, -1.735086904196669e+00, -6.439103726460642e+01, 2.800618690312799e+00, 3.470512804778828e+01, -3.693763771590359e+00,
      -1.130810568251768e-01, -3.472221936807356e+00, 4.581601740226479e-01, 1.844499217660435e+00, -6.439103726460642e+01, -2.800618690312798e+00,
      -3.470512804778826e+01, -3.693763771590358e+00, -1.130810568251744e-01, -3.472221936807353e+00, -4.581601740226480e-01, -1.844499217660433e+00,
      -1.184345736540221e+02, -4.772227322880531e+01, 1.528658522264971e+01, -1.885576694849735e+01, -7.949486428329481e+00, 5.482041630075334e+00,
      -3.522956906744867e+00, 1.735086904196668e+00, -1.184345736540222e+02, -4.772227322880531e+01, -1.528658522264971e+01, 1.885576694849737e+01,
      -7.949486428329483e+00, 5.482041630075329e+00, 3.522956906744854e+00, 1.735086904196673e+00, -6.439103726460645e+01, -2.800618690312791e+00,
      3.470512804778826e+01, 3.693763771590357e+00, -1.130810568251665e-01, -3.472221936807356e+00, 4.581601740226530e-01, -1.844499217660438e+00,
      -4.756497825936611e+01, 1.895063169600745e+01, -3.075791330671029e+01, 1.124332327977372e+01, 6.034862433334961e+00, -6.535647420058913e+00,
      4.578919320749178e+00, 1.819994890060988e+00, -1.873041663983349e+02, 5.491880304483681e+01, -3.447863345139872e+01, 4.700745207684028e+00,
      3.694859988763531e+01, 7.558772796036230e+00, 1.042001250006304e+01, -2.318735158804888e+00, -1.873041663983350e+02, 5.491880304483683e+01,
      3.447863345139870e+01, -4.700745207684037e+00, 3.694859988763528e+01, 7.558772796036221e+00, -1.042001250006305e+01, -2.318735158804896e+00,
      -4.756497825936611e+01, 1.895063169600746e+01, 3.075791330671029e+01, -1.124332327977372e+01, 6.034862433334961e+00, -6.535647420058913e+00,
      -4.578919320749179e+00, 1.819994890060991e+00,
    };
    double m2Correct[] = {
      1.346810209205576e+02, 5.636184644070686e+01, 9.081689031496188e+01, 3.638800942322894e+01, -1.531915442293607e+01, 2.162068038060000e+01,
      -1.154920978322388e+01, 8.056050757541808e+00, 6.437526165222573e+02, 2.383699952856996e+02, 1.543494567891024e+02, 4.915181136302768e+01,
      -1.058058544975506e+02, -2.407700215608105e+01, -3.508892496659295e+01, -1.032850505731126e+01, 6.437526165222575e+02, 2.383699952856996e+02,
      -1.543494567891024e+02, -4.915181136302770e+01, -1.058058544975505e+02, -2.407700215608102e+01, 3.508892496659290e+01, -1.032850505731125e+01,
      1.346810209205576e+02, 5.636184644070686e+01, -9.081689031496188e+01, -3.638800942322894e+01, -1.531915442293607e+01, 2.162068038059999e+01,
      1.154920978322387e+01, 8.056050757541811e+00, 1.921238511758901e+02, -6.770111068912961e+00, 1.140587007866537e+02, -9.496106164226617e+00,
      -2.226967471797892e-01, 1.915556929571275e+01, 6.944928133740057e-01, -5.022783835435162e+00, 4.823471699659661e+02, -1.893777713359310e+02,
      -1.637529707819997e+01, -9.468924055391847e+01, 2.310532509390621e+01, -3.594756731257781e+01, 1.302051527839764e+01, 3.653556204232846e-01,
      4.823471699659663e+02, -1.893777713359309e+02, 1.637529707819996e+01, 9.468924055391845e+01, 2.310532509390625e+01, -3.594756731257780e+01,
      -1.302051527839767e+01, 3.653556204232817e-01, 1.921238511758900e+02, -6.770111068912972e+00, -1.140587007866538e+02, 9.496106164226617e+00,
      -2.226967471797914e-01, 1.915556929571276e+01, -6.944928133739985e-01, -5.022783835435162e+00, 1.921238511758901e+02, 6.770111068912968e+00,
      1.140587007866538e+02, 9.496106164226621e+00, -2.226967471797906e-01, 1.915556929571276e+01, 6.944928133739963e-01, 5.022783835435160e+00,
      4.823471699659665e+02, 1.893777713359309e+02, -1.637529707819991e+01, 9.468924055391845e+01, 2.310532509390623e+01, -3.594756731257780e+01,
      1.302051527839768e+01, -3.653556204232722e-01, 4.823471699659665e+02, 1.893777713359309e+02, 1.637529707819994e+01, -9.468924055391849e+01,
      2.310532509390620e+01, -3.594756731257780e+01, -1.302051527839765e+01, -3.653556204232893e-01, 1.921238511758901e+02, 6.770111068912962e+00,
      -1.140587007866538e+02, -9.496106164226616e+00, -2.226967471797981e-01, 1.915556929571276e+01, -6.944928133740064e-01, 5.022783835435166e+00,
      1.346810209205576e+02, -5.636184644070686e+01, 9.081689031496187e+01, -3.638800942322894e+01, -1.531915442293607e+01, 2.162068038060000e+01,
      -1.154920978322388e+01, -8.056050757541808e+00, 6.437526165222574e+02, -2.383699952856996e+02, 1.543494567891024e+02, -4.915181136302766e+01,
      -1.058058544975506e+02, -2.407700215608102e+01, -3.508892496659293e+01, 1.032850505731125e+01, 6.437526165222575e+02, -2.383699952856996e+02,
      -1.543494567891024e+02, 4.915181136302769e+01, -1.058058544975505e+02, -2.407700215608101e+01, 3.508892496659291e+01, 1.032850505731126e+01,
      1.346810209205576e+02, -5.636184644070686e+01, -9.081689031496188e+01, 3.638800942322894e+01, -1.531915442293607e+01, 2.162068038060000e+01,
      1.154920978322388e+01, -8.056050757541813e+00,
    };
    for (int i=0; i<cells[0]; i++) {
      for (int j=0; j<cells[1]; j++) {
        int idx[] = {i+confGhost[0], j+confGhost[1]};
        // Check M0.
        double *m0ptr = gkyl_array_fetch(m0, gkyl_range_idx(&confLocal, idx)); 
        for (int k=0; k<confBasis.num_basis; k++)
          TEST_CHECK( gkyl_compare( m0Correct[(i*cells[1]+j)*confBasis.num_basis+k], m0ptr[k], 1e-12) );
        // Check M1.
        double *m1ptr = gkyl_array_fetch(m1, gkyl_range_idx(&confLocal, idx)); 
        for (int k=0; k<confBasis.num_basis; k++)
          TEST_CHECK( gkyl_compare( m1Correct[(i*cells[1]+j)*confBasis.num_basis+k], m1ptr[k], 1e-12) );
        // Check M2.
        double *m2ptr = gkyl_array_fetch(m2, gkyl_range_idx(&confLocal, idx)); 
        for (int k=0; k<confBasis.num_basis; k++)
          TEST_CHECK( gkyl_compare( m2Correct[(i*cells[1]+j)*confBasis.num_basis+k], m2ptr[k], 1e-12) );
      }
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1calc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(M0_t); gkyl_mom_type_release(M1_t); gkyl_mom_type_release(M2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void test_1x1v_p1() { test_1x1v(1); } 
void test_1x1v_p2() { test_1x1v(2); } 
void test_1x2v_p1() { test_1x2v(1); } 
void test_1x2v_p2() { test_1x2v(2); } 
void test_2x2v_p1() { test_2x2v(1); } 
void test_2x2v_p2() { test_2x2v(2); } 

TEST_LIST = {
  { "mom_gyrokinetic", test_mom_gyrokinetic },
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x2v_p1", test_1x2v_p1 },
  { "test_1x1v_p2", test_1x1v_p2 },
  { "test_1x2v_p2", test_1x2v_p2 },
  { "test_2x2v_p1", test_2x2v_p1 },
  { "test_2x2v_p2", test_2x2v_p2 },
//  { "test_3x2v_p1", test_3x2v_p1 },
#ifdef GKYL_HAVE_CUDA
//  { "cu_mom_gyrokinetic", test_cu_mom_gyrokinetic },
//  { "test_1x1v_p1_cu", test_1x1v_p1_cu },
//  { "test_1x2v_p1_cu", test_1x2v_p1_cu },
//  { "test_2x2v_p1_cu", test_2x2v_p1_cu },
//  { "test_3x2v_p1_cu", test_3x2v_p1_cu },
#endif
  { NULL, NULL },
};
