#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_mirror_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mirror.h>

#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_app_priv.h>

#include <mpack.h>


void
test_lores()
{

  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 1e-3, -0.01, -M_PI+1e-14 };
  double cupper[] = { 3e-3,  0.01,  M_PI-1e-14 };

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_mirror_geo_grid_inp ginp = {
    .rclose = 0.2,
    .zmin = -2.0,
    .zmax =  2.0,
  };

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_MIRROR,
    .efit_info = inp,
    .mirror_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_mirror_new(&geometry_inp); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_hires()
{
  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/wham_hires.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 1e-3, -0.01, -M_PI+1e-14 };
  double cupper[] = { 3e-3,  0.01,  M_PI-1e-14 };

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


struct gkyl_mirror_geo_grid_inp ginp = {
  .rclose = 0.2,
  .zmin = -2.0,
  .zmax =  2.0,
};

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_MIRROR,
    .efit_info = inp,
    .mirror_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_mirror_new(&geometry_inp); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
gyrokinetic_geom_new(struct gkyl_gk *gk)
{
  disable_denorm_float();
  gkyl_gyrokinetic_app *app = gkyl_malloc(sizeof(gkyl_gyrokinetic_app));
  int cdim = app->cdim = gk->cdim;
  int vdim = app->vdim = gk->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = gk->poly_order;
  app->use_gpu = false; // can't use GPUs if we don't have them!
  app->num_periodic_dir = gk->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = gk->periodic_dirs[d];
  strcpy(app->name, gk->name);
  app->tcurr = 0.0; // reset on init
  app->basis_on_dev.basis = &app->basis;
  app->basis_on_dev.neut_basis = &app->neut_basis;
  app->basis_on_dev.confBasis = &app->confBasis;
  gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
  gkyl_cart_modal_gkhybrid(&app->basis, cdim, vdim); // p=2 in vparallel
  gkyl_cart_modal_hybrid(&app->neut_basis, cdim, vdim+1); // p=2 in v for neutral species
  gkyl_rect_grid_init(&app->grid, cdim, gk->lower, gk->upper, gk->cells);
  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);
  memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
  memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));
  int cuts[3] = { 1, 1, 1 };
  struct gkyl_rect_decomp *rect_decomp =
    gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
  app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = rect_decomp,
      .use_gpu = app->use_gpu
    }
  );
  gkyl_rect_decomp_release(rect_decomp);
  // local skin and ghost ranges for configuration space fields
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }
  // Configuration space geometry initialization
  // Initialize the input struct from user side input struct
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = gk->geometry.geometry_id,
    .c2p_ctx = gk->geometry.c2p_ctx,
    .mapc2p = gk->geometry.mapc2p,
    .bmag_ctx = gk->geometry.bmag_ctx,
    .bmag_func = gk->geometry.bmag_func,
    .efit_info = gk->geometry.efit_info,
    .tok_grid_info = gk->geometry.tok_grid_info,
    .mirror_grid_info = gk->geometry.mirror_grid_info,
    .grid = app->grid,
    .local = app->local,
    .local_ext = app->local_ext,
    .global = app->global,
    .global_ext = app->global_ext,
    .basis = app->confBasis,
  };
  for(int i = 0; i<3; i++)
    geometry_inp.world[i] = gk->geometry.world[i];
  if(app->cdim < 3){
    geometry_inp.geo_grid = gkyl_gk_geometry_augment_grid(app->grid, geometry_inp);
    gkyl_cart_modal_serendip(&geometry_inp.geo_basis, 3, poly_order);
    int ghost[] = { 1, 1, 1 };
    gkyl_create_grid_ranges(&geometry_inp.geo_grid, ghost, &geometry_inp.geo_global_ext, &geometry_inp.geo_global);
    if (gk->has_low_inp) {
      // create local and local_ext from user-supplied local range
      gkyl_gk_geometry_augment_local(&gk->low_inp.local_range, ghost, &geometry_inp.geo_local_ext, &geometry_inp.geo_local);
    }
    else {
      // global and local ranges are same, and so just copy
      memcpy(&geometry_inp.geo_local, &geometry_inp.geo_global, sizeof(struct gkyl_range));
      memcpy(&geometry_inp.geo_local_ext, &geometry_inp.geo_global_ext, sizeof(struct gkyl_range));
    }
  }
  else{
    geometry_inp.geo_grid = app->grid;
    geometry_inp.geo_local = app->local;
    geometry_inp.geo_local_ext = app->local_ext;
    geometry_inp.geo_global = app->global;
    geometry_inp.geo_global_ext = app->global_ext;
    geometry_inp.geo_basis = app->confBasis;
  }
  struct gk_geometry* gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mirror_new(&geometry_inp);
  // deflate geometry if necessary
  if(app->cdim < 3)
    app->gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_inp);
  else
    app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
  gkyl_gk_geometry_bmag_mid(app->gk_geom); // set bmag mid
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  int bcast_rank = comm_sz/2;
  gkyl_comm_array_bcast_host(app->comm, app->gk_geom->bmag_mid, app->gk_geom->bmag_mid, bcast_rank);
  gkyl_gyrokinetic_app_write_geometry(app);
}

void
test_geom_2x() {
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  struct gkyl_efit_inp inp = {
    .filepath = "./data/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  struct gkyl_mirror_geo_grid_inp ginp = {
    .rclose = 0.2,
    .zmin = -2.0,
    .zmax =  2.0,
  };
  double psi_min = 2e-4, psi_max = 3e-3;
  double z_min = -M_PI+1e-1, z_max = M_PI-1e-1;
  int cells_psi = 4, cells_z = 32;
  struct gkyl_gk app_inp = {
    .name = "gk_wham_2x2v",
    .cdim = 2,  .vdim = 2,
    .lower = {psi_min, z_min},
    .upper = {psi_max, z_max},
    .cells = { cells_psi, cells_z },
    .poly_order = 1,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {0.0},
      .efit_info = inp,
      .mirror_grid_info = ginp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
  };
  gyrokinetic_geom_new(&app_inp);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}


TEST_LIST = {
  { "test_lores", test_lores },
  { "test_hires", test_hires },
  { "test_geom_2x", test_geom_2x },
  { NULL, NULL },
};
