#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_gyrokinetic_priv.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

void
test_dg_gyrokinetic()
{
  // initialize grid and ranges
  int cdim = 3, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0};
  double lower[] = {0., 0., 0., -1., 0.};
  double upper[] = {1., 1., 1., 1., 1.};

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

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Initialize geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&confGrid, &confRange, &confRange_ext, &confBasis, 
    mapc2p, 0, bmag_func, 0, false);

  double charge = 1.;
  double mass = 1.;

  struct gkyl_dg_eqn* eqn = gkyl_dg_gyrokinetic_new(&confBasis, &basis, &confRange, &phaseRange, charge, mass, 0, gk_geom, false);

  TEST_CHECK( eqn->num_equations == 1 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  TEST_CHECK( gyrokinetic->cdim == 3 );
  TEST_CHECK( gyrokinetic->pdim == 5 );
  TEST_CHECK( gyrokinetic->conf_range.volume == 32768 );

  gkyl_gk_geometry_release(gk_geom);  
  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

/* int cu_gyrokinetic_test(const struct gkyl_dg_eqn *eqn); */

/* void */
/* test_cu_dg_gyrokinetic() */
/* { */
/*   struct gkyl_basis cbasis, pbasis; */
/*   gkyl_cart_modal_serendip(&cbasis, 1, 1); */
/*   gkyl_cart_modal_serendip(&pbasis, 2, 1); */

/*   struct gkyl_range crange; */
/*   gkyl_range_init_from_shape(&crange, 1, (int[]) { 100 } ); */

/*   struct gkyl_dg_eqn* eqn = gkyl_dg_gyrokinetic_cu_dev_new(&cbasis, &pbasis, &crange); */

/*   // this is not possible from user code and should NOT be done. This */
/*   // is for testing only */
/*   struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn); */

/*   TEST_CHECK( gyrokinetic->cdim == 1 ); */
/*   TEST_CHECK( gyrokinetic->pdim == 2 ); */
/*   TEST_CHECK( gyrokinetic->conf_range.volume == 100 ); */

/*   int nfail = cu_gyrokinetic_test(eqn->on_dev); */

/*   TEST_CHECK( nfail == 0 ); */

/*   gkyl_dg_eqn_release(eqn); */
/* } */

#endif

TEST_LIST = {
  { "dg_gyrokinetic", test_dg_gyrokinetic },
#ifdef GKYL_HAVE_CUDA
/*  { "cu_dg_gyrokinetic", test_cu_dg_gyrokinetic }, */
#endif  
  { NULL, NULL },
};
