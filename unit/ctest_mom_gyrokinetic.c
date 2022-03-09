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

  struct gkyl_mom_type *m2 = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2");

  TEST_CHECK( m2->cdim == 1 );
  TEST_CHECK( m2->pdim == 3 );
  TEST_CHECK( m2->poly_order == 2 );
  TEST_CHECK( m2->num_config == confBasis.num_basis );
  TEST_CHECK( m2->num_phase == basis.num_basis );
  TEST_CHECK( m2->num_mom == 1 );

  struct gkyl_mom_type *m3par = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M3par");
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

void distf_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1];
  double bmag[1];
  bmag_1x(t, xn, &bmag[0], ctx); 
  fout[0] = bmag[0]*(x*x)*(vpar-0.5)*(vpar-0.5);
}

void
test_1x1v_p1()
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

  // create distribution function array and project distribution function on basis
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, distf_1x1v, NULL);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // create bmag array and project magnetic field amplitude function on basis
  struct gkyl_array *bmag;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *projbmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, bmag_1x, NULL);
  gkyl_proj_on_basis_advance(projbmag, 0.0, &confLocal, bmag);

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0");
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1");
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2");
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t);

  // create moment arrays
  struct gkyl_array *m0, *m1, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, distf, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);

  // Check M0.
  double *m00 = gkyl_array_fetch(m0, 0+confGhost[0]); double *m01 = gkyl_array_fetch(m0, 1+confGhost[0]);
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  1.52537436025689e+01, m00[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  3.57234787161818e+00, m00[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  6.08286277145113e+00, m01[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -5.10949091314880e+00, m01[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  6.08286277145113e+00, m02[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  5.10949091314880e+00, m02[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  1.52537436025689e+01, m03[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -3.57234787161818e+00, m03[1], 1e-12) );

  // Check M1.
  double *m10 = gkyl_array_fetch(m1, 0+confGhost[0]); double *m11 = gkyl_array_fetch(m1, 1+confGhost[0]);
  double *m12 = gkyl_array_fetch(m1, 2+confGhost[0]); double *m13 = gkyl_array_fetch(m1, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare( -1.28452577705844e+01, m10[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -3.00829294452057e+00, m10[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -5.12241075490621e+00, m11[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  4.30272919002004e+00, m11[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -5.12241075490621e+00, m12[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -4.30272919002004e+00, m12[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.28452577705844e+01, m13[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  3.00829294452057e+00, m13[1], 1e-12) );

  // Check M2.
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  3.31835825740096e+01, m20[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  7.77142344001148e+00, m20[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  1.32328944501744e+01, m21[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -1.11153837408851e+01, m21[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  1.32328944501744e+01, m22[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  1.11153837408851e+01, m22[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  3.31835825740096e+01, m23[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -7.77142344001148e+00, m23[1], 1e-12) );

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1calc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(M0_t); gkyl_mom_type_release(M1_t); gkyl_mom_type_release(M2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_1x2v_p1()
{
  double mass = 1.;
  int poly_order = 1;
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
    poly_order+1, 1, distf_1x1v, NULL);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  // create bmag array and project magnetic field amplitude function on basis
  struct gkyl_array *bmag;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *projbmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, bmag_1x, NULL);
  gkyl_proj_on_basis_advance(projbmag, 0.0, &confLocal, bmag);

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0");
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1");
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2");
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t);

  // create moment arrays
  struct gkyl_array *m0, *m1, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // compute the moments
  gkyl_mom_calc_advance(m0calc, &local, &confLocal, distf, m0);
  gkyl_mom_calc_advance(m1calc, &local, &confLocal, distf, m1);
  gkyl_mom_calc_advance(m2calc, &local, &confLocal, distf, m2);

  // Check M0.
  double *m00 = gkyl_array_fetch(m0, 0+confGhost[0]); double *m01 = gkyl_array_fetch(m0, 1+confGhost[0]);
  double *m02 = gkyl_array_fetch(m0, 2+confGhost[0]); double *m03 = gkyl_array_fetch(m0, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(  191.6841953662915, m00[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  44.89144731817122, m00[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   76.4395079823428, m01[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -64.2077564653283, m01[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  76.43950798234282, m02[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  64.20775646532829, m02[1], 1e-12) );
  TEST_CHECK( gkyl_compare( 191.68419536629153, m03[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -44.89144731817124, m03[1], 1e-12) );

  // Check M1.
  double *m10 = gkyl_array_fetch(m1, 0+confGhost[0]); double *m11 = gkyl_array_fetch(m1, 1+confGhost[0]);
  double *m12 = gkyl_array_fetch(m1, 2+confGhost[0]); double *m13 = gkyl_array_fetch(m1, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare( -161.41826978214021, m10[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -37.80332405740736, m10[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -64.37011198513079, m11[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  54.069689655013306, m11[1], 1e-12) );
  TEST_CHECK( gkyl_compare(  -64.37011198513079, m12[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -54.06968965501329, m12[1], 1e-12) );
  TEST_CHECK( gkyl_compare( -161.41826978214021, m13[0], 1e-12) );
  TEST_CHECK( gkyl_compare(   37.80332405740738, m13[1], 1e-12) );

  // Check M2.
  double *m20 = gkyl_array_fetch(m2, 0+confGhost[0]); double *m21 = gkyl_array_fetch(m2, 1+confGhost[0]);
  double *m22 = gkyl_array_fetch(m2, 2+confGhost[0]); double *m23 = gkyl_array_fetch(m2, 3+confGhost[0]);
  TEST_CHECK( gkyl_compare(   578.5971330556217, m20[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  210.75432886828457, m20[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   292.8699479183419, m21[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  -242.1332009483718, m21[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   292.8699479183419, m22[0], 1e-12) );
  TEST_CHECK( gkyl_compare(  242.13320094837178, m22[1], 1e-12) );
  TEST_CHECK( gkyl_compare(   578.5971330556217, m23[0], 1e-12) );
  TEST_CHECK( gkyl_compare( -210.75432886828463, m23[1], 1e-12) );

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1); gkyl_array_release(m2);
  gkyl_mom_calc_release(m0calc); gkyl_mom_calc_release(m1calc); gkyl_mom_calc_release(m2calc);
  gkyl_mom_type_release(M0_t); gkyl_mom_type_release(M1_t); gkyl_mom_type_release(M2_t);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}


TEST_LIST = {
  { "mom_gyrokinetic", test_mom_gyrokinetic },
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x2v_p1", test_1x2v_p1 },
//  { "test_2x2v_p1", test_2x2v_p1 },
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
