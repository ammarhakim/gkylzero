#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_comm.h>
#include <gkyl_const.h>
#include <gkyl_dflt.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>
#include <gkyl_mirror_geo.h>
#include <gkyl_mirror_geo_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <rt_arg_parse.h>

struct gkyl_mirror_geo_efit_inp inp = {
  // psiRZ and related inputs
  .filepath = "./efit_data/wham.geqdsk",
  .rzpoly_order = 2,
  .fluxpoly_order = 1,
  .plate_spec = false,
  .quad_param = {  .eps = 1e-10 }
};


struct gkyl_mirror_geo_grid_inp ginp = {
  .rclose = 0.2,
  .zmin = -2.48,
  .zmax =  2.48,
  .write_node_coord_array = true,
  .node_file_nm = "bmag.gkyl",
  .nonuniform_mapping_fraction = 0.7,
};

void mirror_geometry_new(struct gkyl_gk *gk)
{
    // Copy directly from zero/gyrokinetic.c
  disable_denorm_float();

  assert(gk->num_species <= GKYL_MAX_SPECIES);

  gkyl_gyrokinetic_app *app = gkyl_malloc(sizeof(gkyl_gyrokinetic_app));

  int cdim = app->cdim = gk->cdim;
  int vdim = app->vdim = gk->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = gk->poly_order;
  int ns = app->num_species = gk->num_species;
  int neuts = app->num_neut_species = gk->num_neut_species;

  double cfl_frac = gk->cfl_frac == 0 ? 1.0 : gk->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = gk->use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = gk->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = gk->periodic_dirs[d];

  strcpy(app->name, gk->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev.basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.neut_basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.confBasis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev.basis = &app->basis;
    app->basis_on_dev.neut_basis = &app->neut_basis;
    app->basis_on_dev.confBasis = &app->confBasis;
  }

  // basis functions
  switch (gk->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
        gkyl_cart_modal_serendip(&app->neut_basis, pdim+1, poly_order); // neutral species are 3v
      }
      else if (poly_order == 1) {
        gkyl_cart_modal_gkhybrid(&app->basis, cdim, vdim); // p=2 in vparallel
        gkyl_cart_modal_hybrid(&app->neut_basis, cdim, vdim+1); // p=2 in v for neutral species
      }

      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (poly_order > 1) {
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.neut_basis, pdim+1, poly_order); // neutral species are 3v
        }
        else if (poly_order == 1) {
          gkyl_cart_modal_gkhybrid_cu_dev(app->basis_on_dev.basis, cdim, vdim); // p=2 in vparallel
          gkyl_cart_modal_hybrid_cu_dev(app->basis_on_dev.neut_basis, cdim, vdim+1); // p=2 in v for neutral species
        }
      }
      break;
    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, gk->lower, gk->upper, gk->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (gk->has_low_inp) {
    // create local and local_ext from user-supplied local range
    gkyl_create_ranges(&gk->low_inp.local_range, ghost, &app->local_ext, &app->local);
    
    if (gk->low_inp.comm)
      app->comm = gkyl_comm_acquire(gk->low_inp.comm);
    else {
      int cuts[3] = { 1, 1, 1 };
      struct gkyl_rect_decomp *rect_decomp =
        gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
      
      app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
          .decomp = rect_decomp,
          .use_gpu = app->use_gpu
        }
      );

      gkyl_rect_decomp_release(rect_decomp);
    }
  }
  else {
    // global and local ranges are same, and so just copy
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
  }
  // local skin and ghost ranges for configuration space fields
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }


  // Configuration space geometry initialization
  struct gkyl_rect_grid geo_grid;
  struct gkyl_range geo_local;
  struct gkyl_range geo_local_ext;
  struct gkyl_basis geo_basis;
  bool geo_3d_use_gpu = app->use_gpu;

  if(app->cdim < 3){
    geo_grid = agument_grid(app->grid, gk->geometry);
    gkyl_create_grid_ranges(&geo_grid, ghost, &geo_local_ext, &geo_local);
    geo_3d_use_gpu = false;
    switch (gk->basis_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        gkyl_cart_modal_serendip(&geo_basis, 3, poly_order);
        break;
      default:
        assert(false);
        break;
    }
  }
  else{
    geo_grid = app->grid;
    geo_local = app->local;
    geo_local_ext = app->local_ext;
    geo_basis = app->confBasis;
  }

  struct gk_geometry* gk_geom_3d;
  switch (gk->geometry.geometry_id) {
    case GKYL_GEOMETRY_FROMFILE:
      gk_geom_3d = gkyl_gk_geometry_fromfile_new(&app->grid, &app->local, &app->local_ext, &app->confBasis, app->use_gpu);
      break;
    case GKYL_TOKAMAK:
      gk_geom_3d = gkyl_gk_geometry_tok_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
          gk->geometry.tok_efit_info, gk->geometry.tok_grid_info, geo_3d_use_gpu);
      break;
    case GKYL_MIRROR:
      gk_geom_3d = gkyl_gk_geometry_mirror_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
          gk->geometry.mirror_efit_info, gk->geometry.mirror_grid_info, geo_3d_use_gpu);
      break;
    case GKYL_MAPC2P:
      gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
          gk->geometry.mapc2p, gk->geometry.c2p_ctx, gk->geometry.bmag_func,  gk->geometry.bmag_ctx, geo_3d_use_gpu);
      break;
  }
}

void
test_mirror_geometry_new_1x2v_nonuniform()
{
  struct gkyl_mirror_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./efit_data/wham.geqdsk",
    .rzpoly_order = 2,
    .fluxpoly_order = 1,
    .plate_spec = false,
    .quad_param = {  .eps = 1e-10 }
  };


  struct gkyl_mirror_geo_grid_inp ginp = {
    .rclose = 0.2,
    .zmin = -2.48,
    .zmax =  2.48,
    .write_node_coord_array = true,
    .node_file_nm = "bmag.gkyl",
    .nonuniform_mapping_fraction = 0.7, // To do a non-uniform mapping, this parameter should be set to a
    // value between 0 and 1. Zero is for no mapping and 1 is for full mapping. One may not desire full 
    // mapping to maintain cell density in the center of the domain.
  };

  struct gkyl_gk gk = {  // GK app
    .name = "test_mirror_geometry_new_1x2v",
    .cdim = 1,
    .vdim = 2,
    .lower = {-2.48},
    .upper = {2.48},
    .cells = {10},
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {0.02, 0.0},
      .mirror_efit_info = &inp,
      .mirror_grid_info = &ginp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 1,
    .use_gpu = false,
  };
  mirror_geometry_new(&gk);
}

void
test_mirror_geometry_new_1x2v_uniform()
{
  struct gkyl_mirror_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./efit_data/wham.geqdsk",
    .rzpoly_order = 2,
    .fluxpoly_order = 1,
    .plate_spec = false,
    .quad_param = {  .eps = 1e-10 }
  };


  struct gkyl_mirror_geo_grid_inp ginp = {
    .rclose = 0.2,
    .zmin = -2.48,
    .zmax =  2.48,
    .write_node_coord_array = true,
    .node_file_nm = "bmag.gkyl",
    // To do uniform mapping, one may simply neglect the .nonuniform_mapping_fraction parameter, or set it to 0, or 0.0.
  };

  struct gkyl_gk gk = {  // GK app
    .name = "test_mirror_geometry_new_1x2v",
    .cdim = 1,
    .vdim = 2,
    .lower = {-2.48},
    .upper = {2.48},
    .cells = {10},
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {0.02, 0.0},
      .mirror_efit_info = &inp,
      .mirror_grid_info = &ginp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 1,
    .use_gpu = false,
  };
  mirror_geometry_new(&gk);
}

void
test_uniform_grid()
{
  // First read the mapc2p file by defining the appropriate grid and range
  double lower[] = {-2.48};
  double upper[] = {2.48};
  int cells[] = {10};
  int ndim = 1;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range local_ext, local;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Now read the mapc2p file
  struct gkyl_array* mapc2p = gkyl_grid_array_new_from_file(&grid, "mapc2p.gkyl");
}

void
test_mirror_geometry_new_2x2v()
{
  struct gkyl_mirror_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./efit_data/wham.geqdsk",
    .rzpoly_order = 2,
    .fluxpoly_order = 1,
    .plate_spec = false,
    .quad_param = {  .eps = 1e-10 }
  };


  struct gkyl_mirror_geo_grid_inp ginp = {
    .rclose = 0.2,
    .zmin = -2.48,
    .zmax =  2.48,
    .write_node_coord_array = true,
    .node_file_nm = "bmag.gkyl",
    .nonuniform_mapping_fraction = 0.7,
  };

  struct gkyl_gk gk = {  // GK app
    .name = "test_mirror_geometry_new_2x2v",
    .cdim = 2,
    .vdim = 2,
    .lower = {0.0, -2.48},
    .upper = {0.02, 2.48},
    .cells = {10, 10},
    .poly_order = 1,
    .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {0.0},
      .mirror_efit_info = &inp,
      .mirror_grid_info = &ginp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 1,
    .use_gpu = false,
  };
  mirror_geometry_new(&gk);
}

TEST_LIST = {
  // { "test_mirror_geometry_new_1x2v_uniform", test_mirror_geometry_new_1x2v_uniform},
  { "test_uniform_grid", test_uniform_grid},
  // { "test_mirror_geometry_new_1x2v_nonuniform", test_mirror_geometry_new_1x2v_nonuniform},
  // { "test_mirror_geometry_new_2x2v", test_mirror_geometry_new_2x2v},
  { NULL, NULL },
};
