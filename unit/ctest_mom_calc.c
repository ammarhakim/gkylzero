// A test of calculation of moments of a distribution function.
//
#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_mom_calc.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_array_rio.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  fout[0] = x*x;
}

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

void
test_1()
{
  int polyOrder = 1;
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {20, 20};
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
  gkyl_cart_modal_serendip(&basis, ndim, polyOrder);
  gkyl_cart_modal_serendip(&confBasis, cdim, polyOrder);

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

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    polyOrder+1, 1, evalFunc, NULL);

  // create distribution function
  struct gkyl_array *distf = mkarr(basis.numBasis, local_ext.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);

  //// bottom left cell
  //double *df0 = gkyl_array_fetch(distf, 0);
  //TEST_CHECK( gkyl_compare( 2.666666666666667, df0[0], 1e-12) );
  //TEST_CHECK( gkyl_compare(-2.309401076758503, df0[1], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df0[2], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df0[3], 1e-12) );

  //// top left cell
  //double *df1 = gkyl_array_fetch(distf, 1);
  //TEST_CHECK( gkyl_compare( 2.666666666666667, df1[0], 1e-12) );
  //TEST_CHECK( gkyl_compare(-2.309401076758503, df1[1], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df1[2], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df1[3], 1e-12) );

  //// bottom right cell
  //double *df2 = gkyl_array_fetch(distf, 2);
  //TEST_CHECK( gkyl_compare( 2.666666666666667, df2[0], 1e-12) );
  //TEST_CHECK( gkyl_compare( 2.309401076758503, df2[1], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df2[2], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df2[3], 1e-12) );

  //// bottom left cell
  //double *df3 = gkyl_array_fetch(distf, 3);
  //TEST_CHECK( gkyl_compare( 2.666666666666667, df3[0], 1e-12) );
  //TEST_CHECK( gkyl_compare( 2.309401076758503, df3[1], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df3[2], 1e-12) );
  //TEST_CHECK( gkyl_compare( 0., df3[3], 1e-12) );

  struct gkyl_mom_type *mtype;
  mtype = gkyl_vlasov_mom_new(&confBasis, &basis, "M0");
  gkyl_mom_calc *mcalc;
  mcalc = gkyl_mom_calc_new(&grid, mtype);
  struct gkyl_array *marr;
  marr = mkarr(mtype->num_mom*confBasis.numBasis, confLocal_ext.volume);

  // compute the moment
  gkyl_mom_calc_advance(mcalc, &local, &confLocal, distf, marr);

  const char *fmt = "%s-%s-%s_%d.gkyl";
  char name[] = "mom_calc";
  char speciesName[] = "distf";
  char momName[] = "M0";
  int frame = 0;
  int sz = snprintf(0, 0, fmt, name, speciesName, momName, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, speciesName, momName, frame);

  gkyl_grid_array_write(&confGrid, &confLocal, marr, fileNm);

  // release memory for moment data object
  gkyl_array_release(marr);
  gkyl_mom_calc_release(mcalc);
  gkyl_mom_type_release(mtype);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
