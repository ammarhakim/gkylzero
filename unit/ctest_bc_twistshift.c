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
#include <gkyl_bc_twistshift.h>
#include <gkyl_velocity_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>

// Meta-data for IO
struct test_bc_twistshift_output_meta {
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
};

// returned gkyl_array_meta must be freed using gyrokinetic_array_meta_release
static struct gkyl_array_meta*
test_bc_twistshift_array_meta_new(struct test_bc_twistshift_output_meta meta)
{
  struct gkyl_array_meta *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);

  mpack_write_cstr(&writer, "polyOrder");
  mpack_write_i64(&writer, meta.poly_order);

  mpack_write_cstr(&writer, "basisType");
  mpack_write_cstr(&writer, meta.basis_type);

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

static void
test_bc_twistshift_array_meta_release(struct gkyl_array_meta *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
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
// Apply periodic BCs along parallel direction
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
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

struct test_bc_twistshift_ctx {
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  double B0;
};

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void eval_bmag_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_bc_twistshift_ctx *pars = ctx;
  double B0 = pars->B0;

  fout[0] = B0;
}

void
shift_3x2v_const(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[2] = {pars->upper[0]-pars->lower[0], pars->upper[1]-pars->lower[1]};
  double dx[2] = {Lx[0]/pars->cells[0], Lx[1]/pars->cells[1]};

//  fout[0] = 2.5*dx[1];
//  fout[0] = 0.4*x+2.5*dx[1];
  fout[0] = 0.6*x+1.8;
}

void
init_donor_3x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vpar = xn[3], mu = xn[4];

  struct test_bc_twistshift_ctx *pars = ctx;
  double Lx[2] = {pars->upper[0]-pars->lower[0], pars->upper[1]-pars->lower[1]};

  double beta[2] = {0.0, 0.0};
//  double sigma[2] = {Lx[0]/7.0, Lx[1]/10.0};
  double sigma[2] = {0.6, 0.2};

  fout[0] = exp( -pow(x-beta[0],2)/(2.0*pow(sigma[0],2))
                 -pow(y-beta[1],2)/(2.0*pow(sigma[1],2)) );  
//  fout[0] = 0.;
//  if (-1.6 < x && x < -1.2 && -0.3 < y && y < 0.)
//    fout[0] = 1.;
}

void
test_bc_twistshift_3x2v(bool use_gpu)
{
  assert(!use_gpu); // 2x test only available on CPUs.

  double vt = 1.0; // Thermal speed.
  double mass = 1.0;
  double B0 = 1.0; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.
  enum gkyl_edge_loc edge = GKYL_LOWER_EDGE; // Lower or upper edge.

  const int poly_order = 1;
  const double lower[] = {-2.0, -1.50, -3.0, -5.0*vt, 0.};
  const double upper[] = { 2.0,  1.50,  3.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  const int cells[] = {80, 40, 4, 2, 1};
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

  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[cdim];
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

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]},
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]},
    .cells = {cells[0], cells[1], cells[2], cells[3], cells[4]},
    .B0 = B0,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = init_donor_3x2v,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf);
  struct gkyl_array_meta *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  gkyl_grid_sub_array_write(&grid, &local, mt, distf, "ctest_bc_twistshift_3x2v_do.gkyl");

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .local_ext_r = local_ext,
    .num_ghost = ghost,
    .basis = basis,
    .grid = grid,
    .shift_func = shift_3x2v_const,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_new(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_ghost.lower_skin[bc_dir].volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, distf, distf);

  // Write out the target in the extended range.
  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = cells[d]+2*ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
  gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf, "ctest_bc_twistshift_3x2v_tar.gkyl");

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

  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &basis_conf,
    &basis, &local, mass, gvm, gk_geom, "ThreeMoments", true, use_gpu);
  int num_mom = gkyl_dg_updater_moment_gyrokinetic_num_mom(mcalc);

  struct gkyl_array *marr = mkarr(use_gpu, num_mom, local_ext.volume);
  struct gkyl_array *marr_ho = marr;
  if (use_gpu)
    marr_ho = mkarr(false, num_mom, local_ext.volume);
  double *red_integ_mom_skin, *red_integ_mom_ghost;
  if (use_gpu) {
    red_integ_mom_skin = gkyl_cu_malloc(sizeof(double[vdim+2]));
    red_integ_mom_ghost = gkyl_cu_malloc(sizeof(double[vdim+2]));
  }
  else {
    red_integ_mom_skin = gkyl_malloc(sizeof(double[vdim+2]));
    red_integ_mom_ghost = gkyl_malloc(sizeof(double[vdim+2]));
  }

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &skin_ghost.upper_skin[bc_dir], &skin_ghost_conf.upper_skin[bc_dir], distf, marr);
  gkyl_array_reduce_range(red_integ_mom_skin, marr, GKYL_SUM, &skin_ghost_conf.upper_skin[bc_dir]);

  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc,
      &skin_ghost.lower_ghost[bc_dir], &skin_ghost_conf.lower_ghost[bc_dir], distf, marr);
  gkyl_array_reduce_range(red_integ_mom_ghost, marr, GKYL_SUM, &skin_ghost_conf.lower_ghost[bc_dir]);

  for (int k=0; k<vdim+2; k++) {
    TEST_CHECK( gkyl_compare(red_integ_mom_skin[k], red_integ_mom_ghost[k], 1e-12));
    TEST_MSG( "integ_mom %d | Expected: %.14e | Got: %.14e\n",k,red_integ_mom_skin[k],red_integ_mom_ghost[k]);
  }

  if (use_gpu) {
    gkyl_cu_free(red_integ_mom_skin);
    gkyl_cu_free(red_integ_mom_ghost);
  }
  else {
    gkyl_free(red_integ_mom_skin);
    gkyl_free(red_integ_mom_ghost);
  }
  if (use_gpu)
    gkyl_array_release(marr_ho);
  gkyl_dg_updater_moment_gyrokinetic_release(mcalc);
  gkyl_array_release(marr);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_bc_twistshift_release(tsup);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);

}

void test_bc_twistshift_3x2v_ho(){ test_bc_twistshift_3x2v(false); }

TEST_LIST = {
  { "test_bc_twistshift_3x2v_ho", test_bc_twistshift_3x2v_ho },
  { NULL, NULL },
};
