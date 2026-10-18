// Additional tests for the PKPM variables updater in 2 configuration-space
// dimensions (gkyl_dg_calc_pkpm_vars, 2x p1).  These complement the 1x tests
// by exercising the 2x kernels and a different magnetic-field orientation.
//
//   u : weak divide rhou_i / rho with uniform fields.
//
// For a 2x p1 orthonormal basis the cell-average basis is 1/sqrt(2^2) = 1/2,
// so a field with physical value V has component-0 equal to V*sqrt(4).
#include <acutest.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// set cell-average physical value `val` of block `blk` for cdim=2 (factor 2).
static void
set_const_block_2x(struct gkyl_array *arr, int blk, int nb, double val)
{
  gkyl_array_shiftc(arr, val * 2.0, blk * nb);
}

static void
mk_setup(
  struct gkyl_rect_grid *grid, struct gkyl_basis *cbasis, struct gkyl_range *local,
  struct gkyl_range *local_ext
)
{
  int poly_order = 1, cdim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {4, 4};
  gkyl_rect_grid_init(grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(cbasis, cdim, poly_order);
  int ghost[] = {1, 1};
  gkyl_create_grid_ranges(grid, ghost, local_ext, local);
}

void
test_u_2x_p1()
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis cbasis;
  struct gkyl_range local, local_ext;
  mk_setup(&grid, &cbasis, &local, &local_ext);
  int nb = cbasis.num_basis;

  struct gkyl_dg_calc_pkpm_vars *up =
    gkyl_dg_calc_pkpm_vars_new(&grid, &cbasis, &local, false);

  struct gkyl_array *moms = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *euler = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *pkpm_u = mkarr(3 * nb, local_ext.volume);
  struct gkyl_array *cell_avg_prim = gkyl_array_new(GKYL_INT, 1, local_ext.volume);
  memset(cell_avg_prim->data, 0, cell_avg_prim->size * cell_avg_prim->esznc);

  double rho = 5.0, rhoux = 10.0, rhouy = -5.0, rhouz = 15.0;
  gkyl_array_clear(moms, 0.0);
  set_const_block_2x(moms, 0, nb, rho);
  gkyl_array_clear(euler, 0.0);
  set_const_block_2x(euler, 0, nb, rhoux);
  set_const_block_2x(euler, 1, nb, rhouy);
  set_const_block_2x(euler, 2, nb, rhouz);

  gkyl_array_clear(pkpm_u, 0.0);
  gkyl_dg_calc_pkpm_vars_u(up, moms, euler, cell_avg_prim, pkpm_u);

  double s = 2.0;
  double ux_exp = rhoux / rho, uy_exp = rhouy / rho, uz_exp = rhouz / rho;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local, iter.idx);
    const double *u = gkyl_array_cfetch(pkpm_u, loc);
    TEST_CHECK(gkyl_compare(u[0 * nb] / s, ux_exp, 1e-12));
    TEST_CHECK(gkyl_compare(u[1 * nb] / s, uy_exp, 1e-12));
    TEST_CHECK(gkyl_compare(u[2 * nb] / s, uz_exp, 1e-12));
  }

  gkyl_array_release(cell_avg_prim);
  gkyl_array_release(pkpm_u);
  gkyl_array_release(euler);
  gkyl_array_release(moms);
  gkyl_dg_calc_pkpm_vars_release(up);
}

TEST_LIST = {
  {"u_2x_p1", test_u_2x_p1},
  {NULL, NULL}
};
