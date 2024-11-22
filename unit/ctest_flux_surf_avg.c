// Test creation and deallocation of updater that applies the
// twist shift BCs.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <mpack.h>
#include <gkyl_array_rio.h>
#include <gkyl_velocity_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

struct test_flux_surf_avg_ctx {
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  double B0;
  double vt;
  double mass;
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

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void eval_bmag_3x(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_flux_surf_avg_ctx *pars = ctx;
  double B0 = pars->B0;

  fout[0] = B0;
}


void
init_field(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double y = xn[1];

  double mu = 0.0;
  double sigma = 0.3;

  fout[0] = ( 1.0/sqrt(2.0*M_PI*pow(sigma,2)) ) * exp( -pow(y-mu,2)/(2.0*pow(sigma,2)) );
}

void
test_flux_surf_avg_3x2v_wcells(const int *cells, enum gkyl_edge_loc edge,
  bool check_distf, bool use_gpu, bool write_f)
{
  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0, -5.0*vt, 0.};
  const double upper[] = { 2.0,  1.50,  3.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  const int vdim = 2;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = cells[cdim+d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis basis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[vdim];
  for (int d=0; d<vdim; d++) ghost_vel[d] = 0;
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  for (int d=cdim; d<ndim; d++) ghost[d] = ghost_vel[d-cdim];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_flux_surf_avg_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]},
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]},
    .cells = {cells[0], cells[1], cells[2], cells[3], cells[4]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_field,
      .ctx = &proj_ctx,
    }
  );

  // Compute the integrated moments of the skin cell and the ghost cell.
  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);


  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id = GKYL_MAPC2P,
    .c2p_ctx = 0,
    .mapc2p = mapc2p,
    .bmag_ctx = &proj_ctx,
    .bmag_func = eval_bmag_3x,
    .grid = grid_conf,
    .local = local_conf,
    .local_ext = local_ext_conf,
    .global = local_conf,
    .global_ext = local_ext_conf,
    .basis = basis_conf,
    .geo_grid = grid_conf,
    .geo_local = local_conf,
    .geo_local_ext = local_ext_conf,
    .geo_global = local_conf,
    .geo_global_ext = local_ext_conf,
    .geo_basis = basis_conf,
  };

  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_inp);
  struct gk_geometry* gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
  if (use_gpu) {  // If we are on the gpu, copy from host
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_inp, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }
  // Need the magnetic field to be initialized in the ghost cell in order to
  // compute integrated moments in the ghost cell (for checking moment
  // conservation).
  gkyl_array_clear(gk_geom->bmag, 0.0);
  gkyl_array_shiftc(gk_geom->bmag, B0*pow(sqrt(2.0),cdim), 0);


  /*
  A:  We compute the total integral of the energy moment
      we will compare it to the sum of each flux surface integral in B.
  */
  //1.--Create an updater for the energy moment M2, (T/m)
  struct gkyl_dg_updater_moment *m2calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &basis_conf,
    &basis, &local_conf, mass, gvm, gk_geom, "M2", true, use_gpu);
  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(m2calc);
  
  //2.--Allocate memory for the flux surface average of the moment over extended local conf space
  struct gkyl_array *m2arr = mkarr(use_gpu, num_mom, local_ext_conf.volume);

  //3.--Allocate memory for the result of the integration
  double *red_total_integ_mom;
  if (use_gpu) {
    red_total_integ_mom = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  else {
    red_total_integ_mom = gkyl_malloc(sizeof(double[vdim+2]));
  }
  //also on the host
  double *red_total_integ_mom_ho = gkyl_malloc(sizeof(double[vdim+2]));

  //4.--Compute the M2 moment over the config space
  gkyl_dg_updater_moment_gyrokinetic_advance(m2calc,
      &local, &local_conf, distf, m2arr);

  //5.--Build up the integrant (TODO)
  /*
  This should compute the integrant to compute the integral and not only the sum over each cell
  */
  struct gkyl_array *integrant = mkarr(use_gpu, num_mom, local_ext_conf.volume);
  gkyl_array_copy(m2arr,integrant);

  //6.--Integration by reduce sum operation of the integrant
  gkyl_array_reduce_range(red_total_integ_mom, integrant, GKYL_SUM, &local_conf);

  if (use_gpu) {
    gkyl_cu_memcpy(red_total_integ_mom_ho, red_total_integ_mom, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(red_total_integ_mom_ho, red_total_integ_mom, sizeof(double[2+vdim]));
  }

  /*
  B:  We compute the flux surface integral of the energy moment for each cell in x
  */
  //0.--Create phase space and configuration subranges to define the flux surface
  struct gkyl_range fs_range_ps, fs_range_conf;
  int *sublower[ndim];
  int *subupper[ndim];
  //set upper = lower = ix where the fs avg must computed
  sublower[0] = 0;
  subupper[0] = 0;
  //fill the other sub limits (from d=1 to avoid x)
  for (int d=1; d<ndim; d++) {
    sublower[d] = 0;
    subupper[d] = cells[d];
  }
  //initialize the fs_range that is a slice of the conf domain at ix=0
  gkyl_sub_range_init(&fs_range_ps, &local, sublower, subupper);
  gkyl_sub_range_init(&fs_range_conf, &local_conf, sublower, subupper);

  //1.--Create an updater for the energy moment M2, (T/m) (on the flux surface only)
  struct gkyl_dg_updater_moment *m2calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &basis_conf,
    &basis, &fs_range_conf, mass, gvm, gk_geom, "M2", true, use_gpu);
  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(m2calc);
  
  //2.--Allocate memory for the flux surface average of the moment over flux surface range
  struct gkyl_array *m2arr = mkarr(use_gpu, num_mom, fs_range_conf.volume);

  //3.--Allocate memory for the result of the integration
  double *red_fs_integ_mom;
  if (use_gpu) {
    red_fs_integ_mom = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  else {
    red_fs_integ_mom = gkyl_malloc(sizeof(double[vdim+2]));
  }
  //also on the host
  double *red_fs_integ_mom_ho = gkyl_malloc(sizeof(double[vdim+2]));

  //4.--Compute the M2 moment over the config space
  gkyl_dg_updater_moment_gyrokinetic_advance(m2calc,
      &fs_range_ps, &fs_range_conf, distf, m2arr);

  //5.--Build up the integrant (TODO)
  /*
  This should compute the integrant to compute the integral and not only the sum over each cell
  */
  struct gkyl_array *integrant = mkarr(use_gpu, num_mom, local_ext_conf.volume);
  gkyl_array_copy(m2arr,integrant);

  //6.--Integration by reduce sum operation of the integrant
  gkyl_array_reduce_range(red_total_integ_mom, integrant, GKYL_SUM, &local_conf);

  if (use_gpu) {
    gkyl_cu_memcpy(red_fs_integ_mom_ho, red_fs_integ_mom, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(red_fs_integ_mom_ho, red_fs_integ_mom, sizeof(double[2+vdim]));
  }


  for (int k=0; k<vdim+2; k++) {
    // TEST_CHECK( gkyl_compare(red_integ_mom_skin_ho[k], red_integ_mom_ghost_ho[k], 1e-12));
    // TEST_MSG( "integ_mom %d | Expected: %.14e | Got: %.14e\n",k,red_integ_mom_skin_ho[k],red_integ_mom_ghost_ho[k]);
  }

  if (check_distf) {
    // Check 0th and 2nd DG coeffs.
    const double f0[] =
    {
      2.5656054446469541e+00, 4.0323377043325592e-01, 2.4549376970832818e-02,
      5.6930416832388616e-04, 5.6930416832388561e-04, 2.4549376970832790e-02,
      4.0323377043325587e-01, 2.5656054446469554e+00, 6.4341275514494081e+00,
      6.4341275514494081e+00,
    };
    const double f2[] =
    {
      -1.0463462350335195e+00, -2.4917980000416121e-01, -1.8802664978950782e-02,
      -4.9044149652606374e-04,  4.9044149652607144e-04,  1.8802664978951094e-02,
       2.4917980000416726e-01,  1.0463462350335588e+00,  9.2229040193425293e-01,
      -9.2229040193415335e-01,
    };

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ghost_rng);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[3]==1 && iter.idx[4]==1) {
        long linidx = gkyl_range_idx(&ghost_rng, iter.idx);
        double *f_c = gkyl_array_fetch(distf_ho, linidx);
        int refidx = (iter.idx[1]-1)*cells[0] + iter.idx[0]-1;
        TEST_CHECK( gkyl_compare(f0[refidx], f_c[0], 1e-13) );
        TEST_CHECK( gkyl_compare(f2[refidx], f_c[2], 1e-13) );
      }
    }
  }

  gkyl_free(red_fs_integ_mom_ho);
  gkyl_free(red_total_integ_mom_ho);
  if (use_gpu) {
    gkyl_cu_free(red_fs_integ_mom);
    gkyl_cu_free(red_total_integ_mom);
  }
  else {
    gkyl_free(red_total_integ_mom);
    gkyl_free(red_total_integ_mom);
  }
  gkyl_dg_updater_moment_gyrokinetic_release(m2calc);
  gkyl_array_release(m2arr);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_velocity_map_release(gvm);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);
}

void
test_flux_surf_avg_3x2v(bool use_gpu)
{
  const int cells0[] = {1, 10, 4, 2, 1};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_flux_surf_avg_3x2v_wcells(cells0, edgelo, true, use_gpu, false);
}


void test_flux_surf_avg_3x2v_ho(){ test_flux_surf_avg_3x2v(false); }

#ifdef GKYL_HAVE_CUDA

void test_flux_surf_avg_3x2v_dev(){ test_flux_surf_avg_3x2v(true); }
#endif

TEST_LIST = {
  { "test_flux_surf_avg_3x2v_ho", test_flux_surf_avg_3x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_flux_surf_avg_3x2v_dev", test_flux_surf_avg_3x2v_dev },
#endif
  { NULL, NULL },
};
