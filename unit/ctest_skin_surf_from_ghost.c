#include "gkyl_array.h"

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_translate_dim_gyrokinetic.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
// Include our new updater
#include <gkyl_skin_surf_from_ghost.h>

// Allocate array (filled with zeros).
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

void eval_field_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
test_3x(int poly_order, bool use_gpu, enum gkyl_edge_loc edge)
{
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 4, 6};
  int dir = 2;

  const int ndim = sizeof(cells)/sizeof(cells[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Ranges
  int ghost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Create the field array.
  struct gkyl_array *field_ho, *field;
  field    = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  field_ho = use_gpu ? mkarr(false, field->ncomp, field->size)
                     : gkyl_array_acquire(field);

  // Fill the inside values
  /* Project the analytical field defined in eval_field_3x to a DG representation 
    and then copy it to the device */
  gkyl_proj_on_basis *proj_field = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_field_3x, NULL);
  gkyl_proj_on_basis_advance(proj_field, 0.0, &local, field_ho);

  // Fill the ghost values
  struct gkyl_range_iter iter_ghost;
  gkyl_range_iter_init(&iter_ghost, &ghost_r);
  while(gkyl_range_iter_next(&iter_ghost)) {

    long linidx = gkyl_range_idx(&ghost_r, iter_ghost.idx);

    double *f_i = gkyl_array_fetch(field_ho, linidx);

    f_i[0] = 1.0;
  }

  // Copy the array from the host to the device
  gkyl_array_copy(field, field_ho);

  // Create a new skin_surf_ghost updater, call it, and release it
  gkyl_skin_surf_from_ghost* up = gkyl_skin_surf_from_ghost_new(dir, edge, &basis, &skin_r, &ghost_r, use_gpu);
  gkyl_skin_surf_from_ghost_advance(up, field);
  gkyl_skin_surf_from_ghost_release(up);

  // Copy the array from the device to the host
  gkyl_array_copy(field_ho, field);

  // Check the result and print the values
  struct gkyl_range_iter iter_skin;
  gkyl_range_iter_init(&iter_skin, &skin_r);
  while (gkyl_range_iter_next(&iter_skin)) {

    long linidx = gkyl_range_idx(&skin_r, iter_skin.idx);

    double *f_i = gkyl_array_fetch(field_ho, linidx);

    // Print the value at each skin location
    // printf("Field value at skin index [%ld]: %f\n", linidx, f_i[0]);
    TEST_CHECK( gkyl_compare( f_i[0], 0.5, 1e-14));
  }

  gkyl_array_release(field);
  gkyl_array_release(field_ho);
}

void test_3x_ho()
{
  test_3x(1, false, GKYL_LOWER_EDGE);
  test_3x(1, false, GKYL_UPPER_EDGE);
}

void test_3x_dev()
{
  test_3x(1, true, GKYL_LOWER_EDGE);
  test_3x(1, true, GKYL_UPPER_EDGE);
}

TEST_LIST = {
  { "test_3x_ho", test_3x_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_3x_dev", test_3x_dev },
#endif
  { NULL, NULL },
};
