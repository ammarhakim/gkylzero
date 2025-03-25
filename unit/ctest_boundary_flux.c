#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_boundary_flux.h>
#include <gkyl_boundary_flux_priv.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_gyrokinetic_priv.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_velocity_map.h>
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
test_gk_boundary_flux_init()
{
  // initialize grid and ranges
  bool use_gpu = false;
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

  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }
  struct gkyl_rect_grid velGrid;
  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  // Local skin and ghost ranges for configuration space fields.
  struct gkyl_range lower_skin[cdim], lower_ghost[cdim], upper_skin[cdim], upper_ghost[cdim];
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&lower_skin[dir], &lower_ghost[dir], dir, GKYL_LOWER_EDGE, &phaseRange_ext, ghost); 
    gkyl_skin_ghost_ranges(&upper_skin[dir], &upper_ghost[dir], dir, GKYL_UPPER_EDGE, &phaseRange_ext, ghost);
  }

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
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func, // magnetic field magnitude
      .bmag_ctx =0 ,
      .grid = confGrid,
      .local = confRange,
      .local_ext = confRange_ext,
      .global = confRange,
      .global_ext = confRange_ext,
      .basis = confBasis,
      .geo_grid = confGrid,
      .geo_local = confRange,
      .geo_local_ext = confRange_ext,
      .geo_global = confRange,
      .geo_global_ext = confRange_ext,
      .geo_basis = confBasis,
  };


  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  double charge = 1.;
  double mass = 1.;

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, phaseGrid, velGrid,
    phaseRange, phaseRange_ext, velLocal, velLocal_ext, false);

  struct gkyl_dg_eqn* eqn = gkyl_dg_gyrokinetic_new(&confBasis, &basis, &confRange, &phaseRange, 
    charge, mass, 0, gk_geom, gvm, false);

  struct gkyl_boundary_flux *flux_slvr[2*cdim]; // boundary flux solver
  
  for (int d=0; d<cdim; ++d) {
    for (int e=0; e<2; ++e) {
      // Allocate solver.
      struct gkyl_range *skin_r = e==0? &lower_skin[d] : &upper_skin[d];
      struct gkyl_range *ghost_r = e==0? &lower_ghost[d] : &upper_ghost[d];
      flux_slvr[2*d+e] = gkyl_boundary_flux_new(d, e, &phaseGrid,
        skin_r, ghost_r, eqn, false, use_gpu);
    }
  }
  TEST_CHECK(flux_slvr[0] != 0);
  TEST_CHECK(flux_slvr[1] != 0);

  TEST_CHECK(flux_slvr[0]->dir == 0);
  TEST_CHECK(flux_slvr[0]->edge == 0);
  TEST_CHECK(flux_slvr[0]->grid.ndim == 5);
  TEST_CHECK(flux_slvr[0]->skin_r.ndim == 5);
  TEST_CHECK(flux_slvr[0]->ghost_r.ndim == 5);
  TEST_CHECK(flux_slvr[0]->equation == eqn);
  TEST_CHECK(flux_slvr[0]->use_boundary_surf == false);
  TEST_CHECK(flux_slvr[0]->use_gpu == use_gpu);

  TEST_CHECK(flux_slvr[1]->dir == 0);
  TEST_CHECK(flux_slvr[1]->edge == 1);
  TEST_CHECK(flux_slvr[1]->grid.ndim == 5);
  TEST_CHECK(flux_slvr[1]->skin_r.ndim == 5);
  TEST_CHECK(flux_slvr[1]->ghost_r.ndim == 5);
  TEST_CHECK(flux_slvr[1]->equation == eqn);
  TEST_CHECK(flux_slvr[1]->use_boundary_surf == false);
  TEST_CHECK(flux_slvr[1]->use_gpu == use_gpu);

  // release resources
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_velocity_map_release(gvm);
  gkyl_dg_eqn_release(eqn);
  for (int d=0; d<cdim; ++d) {
    for (int e=0; e<2; ++e) {
      gkyl_boundary_flux_release(flux_slvr[2*d+e]);
    }
  }
}


TEST_LIST = {
  { "test_gk_boundary_flux_init", test_gk_boundary_flux_init },
  { NULL, NULL },
};
