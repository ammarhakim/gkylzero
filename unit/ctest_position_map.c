#include "gkyl_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_util.h>
#include <acutest.h>
#include <gkyl_position_map.h>

void
mapc2fa(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double poly_order = 2;
  fout[0] = xn[0];
  fout[1] = xn[1];
  double z = xn[2];
  double left = 0.25;
  double right = 0.75;
  if (z < -left)
    fout[2] = z;
  else if (z < right)
    fout[2] = - pow(z - right, poly_order)/pow(left-right, poly_order-1) + right;
  else
    fout[2] = z;
}

void
test_position_map_init_1x()
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
  struct gkyl_position_map_inp pos_map_inp = {
    .ctx = NULL,
    .mapping = mapc2fa,
    .numerical_mapping_fraction = 0.5
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->is_identity == false);
  TEST_ASSERT(pos_map->grid.ndim == 1);
  TEST_ASSERT(pos_map->local.ndim == 1);
  TEST_ASSERT(pos_map->local_ext.ndim == 1);
  TEST_ASSERT(pos_map->pmap_basis->ndim == 1);
  TEST_ASSERT(pos_map->pmap_basis->poly_order == 1);
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
  struct gkyl_position_map_inp pos_map_inp;
  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->is_identity == true);
  TEST_ASSERT(pos_map->grid.ndim == 1);
  TEST_ASSERT(pos_map->local.ndim == 1);
  TEST_ASSERT(pos_map->local_ext.ndim == 1);
  TEST_ASSERT(pos_map->pmap_basis->ndim == 1);
  TEST_ASSERT(pos_map->pmap_basis->poly_order == 1);
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
    .ctx = NULL,
    .mapping = mapc2fa,
    .numerical_mapping_fraction = 0.5
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->is_identity == false);
  TEST_ASSERT(pos_map->grid.ndim == 2);
  TEST_ASSERT(pos_map->local.ndim == 2);
  TEST_ASSERT(pos_map->local_ext.ndim == 2);
  TEST_ASSERT(pos_map->pmap_basis->ndim == 2);
  TEST_ASSERT(pos_map->pmap_basis->poly_order == 1);
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
    .ctx = NULL,
    .mapping = mapc2fa,
    .numerical_mapping_fraction = 0.5
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, localRange, localRange_ext, basis);

  TEST_ASSERT(pos_map->is_identity == false);
  TEST_ASSERT(pos_map->grid.ndim == 3);
  TEST_ASSERT(pos_map->local.ndim == 3);
  TEST_ASSERT(pos_map->local_ext.ndim == 3);
  TEST_ASSERT(pos_map->pmap_basis->ndim == 3);
  TEST_ASSERT(pos_map->pmap_basis->poly_order == 1);
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
    .ctx = NULL,
    .mapping = mapc2fa,
    .numerical_mapping_fraction = 0.5
  };

  struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, localRange, localRange_ext, basis);

  struct gkyl_array *pmap_arr_set = gkyl_array_new(GKYL_DOUBLE, 3*pos_map->pmap_basis->num_basis, pos_map->local_ext.volume);
  gkyl_array_clear(pmap_arr_set, 1.0);

  gkyl_position_map_set(pos_map, pmap_arr_set);

  double *pos_map_i  = pos_map->pmap->data; 
  for (unsigned i=0; i<pos_map->pmap->size; ++i)
    TEST_CHECK( gkyl_compare(pos_map_i[i], 1.0, 1e-14) );

  gkyl_array_release(pmap_arr_set);
  gkyl_position_map_release(pos_map);
}

// Need a test for gkyl_position_map_eval_c2p

TEST_LIST = {
  { "test_position_map_init_1x", test_position_map_init_1x },
  { "test_position_map_init_1x_null", test_position_map_init_1x_null },
  { "test_position_map_init_2x", test_position_map_init_2x },
  { "test_position_map_init_3x", test_position_map_init_3x },
  { "test_position_map_set", test_position_map_set },
  { NULL, NULL },
};