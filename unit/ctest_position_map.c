#include <stdio.h>
#include <stdlib.h>
#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_position_map.h>
#include <gkyl_position_map_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_calc_bmag.h>

void
test_nonuniform_position_map(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double poly_order = 2;
  double z = xn[0];
  double left = 0.25;
  double right = 0.75;
  if (z < -left)
    fout[0] = z;
  else if (z < right)
    fout[0] = - pow(z - right, poly_order)/fabs(pow(left-right, poly_order-1)) + right;
  else
    fout[0] = z;
}

void
test_identity_position_map(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = xn[0];
}


void
test_nonuniform_position_map_3x(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double poly_order = 2;
  double left = 0.25;
  double right = 0.75;
  for (int i = 0; i<3; i++)
  {
    double z = xn[i];
    if (z < -left)
      fout[i] = z;
    else if (z < right)
      fout[i] = - pow(z - right, poly_order)/fabs(pow(left-right, poly_order-1)) + right;
    else
      fout[i] = z;
  }
}

void 
bmag_func(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0];
  double s = 0.6;
  double c = 0.;
  // double B = (4*pow(s*(x-c),2) - 0.3*pow(s*(x-c),4) + 1)*exp(-pow(s*(x-c),2));
  double B = 1/(1+100*pow(x-M_PI/2,2)) + 1/(1+100*pow(x+M_PI/2,2));
  fout[0] = B;
}

void
test_position_map_init_1x()
{
  int cells[] = {32};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .map_x = test_nonuniform_position_map,
    .ctx_x = NULL,
    .map_y = test_nonuniform_position_map,
    .ctx_y = NULL,
    .map_z = test_nonuniform_position_map,
    .ctx_z = NULL,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp,\
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->bmag_ctx->bmag == 0);
  TEST_ASSERT(pos_map->grid.ndim == 1);
  TEST_ASSERT(pos_map->local.ndim == 1);
  TEST_ASSERT(pos_map->local_ext.ndim == 1);
  TEST_ASSERT(pos_map->basis.ndim == 1);
  TEST_ASSERT(pos_map->basis.poly_order == 1);
  TEST_ASSERT(pos_map->flags == 0);

  gkyl_position_map_release(pos_map);
}


void
test_position_map_init_1x_null()
{
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {};
  
  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, \
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->id == GKYL_PMAP_FUNC);
  for (double i = 0; i < 1; i = i+0.1){
    double x[1] = {i};
    double y[1];
    pos_map->map_x(0, x, y, pos_map->map_x_ctx);
    TEST_ASSERT(y[0] == x[0]);
    pos_map->map_y(0, x, y, pos_map->map_y_ctx);
    TEST_ASSERT(y[0] == x[0]);
    pos_map->map_z(0, x, y, pos_map->map_z_ctx);
    TEST_ASSERT(y[0] == x[0]);
  }
  TEST_ASSERT(pos_map->bmag_ctx->bmag == 0);
  TEST_ASSERT(pos_map->grid.ndim == 1);
  TEST_ASSERT(pos_map->local.ndim == 1);
  TEST_ASSERT(pos_map->local_ext.ndim == 1);
  TEST_ASSERT(pos_map->basis.ndim == 1);
  TEST_ASSERT(pos_map->basis.poly_order == 1);
  TEST_ASSERT(pos_map->flags == 0);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_init_2x()
{
  int cells[] = {8,8};
  int poly_order = 1;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .map_x = test_nonuniform_position_map,
    .ctx_x = NULL,
    .map_y = test_nonuniform_position_map,
    .ctx_y = NULL,
    .map_z = test_nonuniform_position_map,
    .ctx_z = NULL,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, \
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->bmag_ctx->bmag == 0);
  TEST_ASSERT(pos_map->grid.ndim == 2);
  TEST_ASSERT(pos_map->local.ndim == 2);
  TEST_ASSERT(pos_map->local_ext.ndim == 2);
  TEST_ASSERT(pos_map->basis.ndim == 2);
  TEST_ASSERT(pos_map->basis.poly_order == 1);
  TEST_ASSERT(pos_map->flags == 0);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_init_3x()
{
  int cells[] = {8, 8, 8};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .map_x = test_nonuniform_position_map,
    .ctx_x = NULL,
    .map_y = test_nonuniform_position_map,
    .ctx_y = NULL,
    .map_z = test_nonuniform_position_map,
    .ctx_z = NULL,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, \
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->bmag_ctx->bmag == 0);
  TEST_ASSERT(pos_map->grid.ndim == 3);
  TEST_ASSERT(pos_map->local.ndim == 3);
  TEST_ASSERT(pos_map->local_ext.ndim == 3);
  TEST_ASSERT(pos_map->basis.ndim == 3);
  TEST_ASSERT(pos_map->basis.poly_order == 1);
  TEST_ASSERT(pos_map->flags == 0);

  gkyl_position_map_release(pos_map);
}

void
test_position_map_set()
{
  int cells[] = {8, 8, 8};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .map_x = test_nonuniform_position_map,
    .ctx_x = NULL,
    .map_y = test_nonuniform_position_map,
    .ctx_y = NULL,
    .map_z = test_nonuniform_position_map,
    .ctx_z = NULL,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, \
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  struct gkyl_array *pmap_arr_set = gkyl_array_new(GKYL_DOUBLE, \
    3*pos_map->basis.num_basis, pos_map->local_ext.volume);
  gkyl_array_clear(pmap_arr_set, 1.0);

  gkyl_position_map_set(pos_map, pmap_arr_set);

  double *pos_map_i  = pos_map->mc2nu->data; 
  for (unsigned i=0; i<pos_map->mc2nu->size; ++i)
    TEST_CHECK( gkyl_compare(pos_map_i[i], 1.0, 1e-14) );

  gkyl_array_release(pmap_arr_set);
  gkyl_position_map_release(pos_map);
}


void
test_gkyl_position_map_eval_mc2nu()
{
  int cells[] = {8, 8, 8};
  int poly_order = 2;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1, 1};
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .map_x = test_nonuniform_position_map,
    .ctx_x = NULL,
    .map_y = test_nonuniform_position_map,
    .ctx_y = NULL,
    .map_z = test_nonuniform_position_map,
    .ctx_z = NULL,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, \
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  struct gkyl_array *pmap_arr_set = gkyl_array_new(GKYL_DOUBLE, \
    3*pos_map->basis.num_basis, pos_map->local_ext.volume);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 3, test_nonuniform_position_map_3x, NULL);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &localRange, pmap_arr_set);
  gkyl_proj_on_basis_release(projDistf);

  gkyl_position_map_set(pos_map, pmap_arr_set);

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<5; k++) {
        double x[3] = {i/10.0, j/10.0, k/10.0};
        double x_fa[3];
        gkyl_position_map_eval_mc2nu(pos_map, x, x_fa);
        double x_analytic[3];
        test_nonuniform_position_map(0.0, &x[0], &x_analytic[0], NULL);
        test_nonuniform_position_map(0.0, &x[1], &x_analytic[1], NULL);
        test_nonuniform_position_map(0.0, &x[2], &x_analytic[2], NULL);
        for (int d=0; d<3; ++d)
          TEST_CHECK( gkyl_compare(x_fa[d], x_analytic[d], 1e-12) );
      }
    }
  }

  gkyl_array_release(pmap_arr_set);
  gkyl_position_map_release(pos_map);
}

void
test_position_polynomial_map_optimize_1x()
{
  int cells[] = {64};
  int poly_order = 1;
  double lower[] = {-M_PI+1e-2}, upper[] = {M_PI-1e-2};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .id = GKYL_PMAP_UNIFORM_B_POLYNOMIAL,
    .map_strength = 1.0,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp,\
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  // Project bmag_func onto bmag_global
  struct gkyl_array *bmag_global = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  gkyl_proj_on_basis *projB = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, bmag_func, NULL);
  gkyl_proj_on_basis_advance(projB, 0.0, &localRange, bmag_global);
  gkyl_proj_on_basis_release(projB);

  // Write bmag_global
  gkyl_grid_sub_array_write(&grid, &localRange,NULL, bmag_global, "bmag_global_1x.gkyl");

  pos_map->bmag_ctx->crange_global = &localRange;
  pos_map->bmag_ctx->cbasis = &basis;
  pos_map->bmag_ctx->cgrid = &grid;
  pos_map->bmag_ctx->bmag = gkyl_array_acquire(bmag_global);

  pos_map->constB_ctx->psi = 0.5;
  pos_map->constB_ctx->psi_min = 0.4;
  pos_map->constB_ctx->psi_max = 0.6;
  pos_map->constB_ctx->alpha = 0.0;
  pos_map->constB_ctx->alpha_min = -0.1;
  pos_map->constB_ctx->alpha_max = 0.1;
  pos_map->constB_ctx->theta_max = upper[0];
  pos_map->constB_ctx->theta_min = lower[0];
  pos_map->constB_ctx->map_x_backup = test_identity_position_map;
  pos_map->constB_ctx->map_x_ctx_backup = NULL;
  pos_map->constB_ctx->map_y_backup = test_identity_position_map;
  pos_map->constB_ctx->map_y_ctx_backup = NULL;

  pos_map->constB_ctx->N_theta_boundaries = cells[0]+1;

  gkyl_position_map_optimize(pos_map);

  TEST_ASSERT(pos_map->bmag_ctx->bmag != 0);
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->theta_throat, 1.565796, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->Bmag_throat, 1.093613, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->psi, 0.5, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->alpha, 0.0, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->map_strength, 1.0, 1e-6) );
  TEST_ASSERT( pos_map->constB_ctx->map_order_center == 2 );
  TEST_ASSERT( pos_map->constB_ctx->map_order_expander == 3 );
  TEST_ASSERT( pos_map->constB_ctx->N_theta_boundaries == 65 );

  gkyl_position_map_release(pos_map);
  gkyl_array_release(bmag_global);
}

void
test_position_map_numeric_optimize_1x()
{
  int cells[] = {64};
  int poly_order = 1;
  double lower[] = {-M_PI+1e-2}, upper[] = {M_PI-1e-2};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);
  // Ranges
  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  
  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);
  struct gkyl_position_map_inp pos_map_inp = {
    .id = GKYL_PMAP_UNIFORM_B_NUMERIC,
    .map_strength = 1.0,
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp,\
    grid, localRange, localRange_ext, localRange, localRange_ext, basis);

  // Project bmag_func onto bmag_global
  struct gkyl_array *bmag_global = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  gkyl_proj_on_basis *projB = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, bmag_func, NULL);
  gkyl_proj_on_basis_advance(projB, 0.0, &localRange, bmag_global);
  gkyl_proj_on_basis_release(projB);

  // Write bmag_global
  gkyl_grid_sub_array_write(&grid, &localRange,NULL, bmag_global, "bmag_global_1x.gkyl");

  pos_map->bmag_ctx->crange_global = &localRange;
  pos_map->bmag_ctx->cbasis = &basis;
  pos_map->bmag_ctx->cgrid = &grid;
  pos_map->bmag_ctx->bmag = gkyl_array_acquire(bmag_global);

  pos_map->constB_ctx->psi = 0.5;
  pos_map->constB_ctx->psi_min = 0.4;
  pos_map->constB_ctx->psi_max = 0.6;
  pos_map->constB_ctx->alpha = 0.0;
  pos_map->constB_ctx->alpha_min = -0.1;
  pos_map->constB_ctx->alpha_max = 0.1;
  pos_map->constB_ctx->theta_max = upper[0];
  pos_map->constB_ctx->theta_min = lower[0];
  pos_map->constB_ctx->map_x_backup = test_identity_position_map;
  pos_map->constB_ctx->map_x_ctx_backup = NULL;
  pos_map->constB_ctx->map_y_backup = test_identity_position_map;
  pos_map->constB_ctx->map_y_ctx_backup = NULL;
  pos_map->constB_ctx->N_theta_boundaries = cells[0]+1;

  gkyl_position_map_optimize(pos_map);

  TEST_ASSERT( pos_map->constB_ctx->num_extrema == 5 );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->theta_extrema[0], -3.131593, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->theta_extrema[1], -1.589886, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->theta_extrema[2], 0.0, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->theta_extrema[3], 1.589886, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->theta_extrema[4], 3.131593, 1e-6) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->bmag_extrema[0], 0.004531427925, 1e-10) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->bmag_extrema[1], 0.957785094052, 1e-10) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->bmag_extrema[2], 0.008057373351, 1e-10) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->bmag_extrema[3], 0.957785094052, 1e-10) );
  TEST_ASSERT( gkyl_compare(pos_map->constB_ctx->bmag_extrema[4], 0.004531427925, 1e-10) );

  gkyl_position_map_release(pos_map);
  gkyl_array_release(bmag_global);
}

TEST_LIST = {
  { "test_position_map_init_1x", test_position_map_init_1x },
  { "test_position_map_init_1x_null", test_position_map_init_1x_null },
  { "test_position_map_init_2x", test_position_map_init_2x },
  { "test_position_map_init_3x", test_position_map_init_3x },
  { "test_position_map_set", test_position_map_set },
  { "test_gkyl_position_map_eval_mc2nu", test_gkyl_position_map_eval_mc2nu }, 
  { "test_position_polynomial_map_optimize_1x", test_position_polynomial_map_optimize_1x },
  { "test_position_map_numeric_optimize_1x", test_position_map_numeric_optimize_1x },
  { NULL, NULL },
};