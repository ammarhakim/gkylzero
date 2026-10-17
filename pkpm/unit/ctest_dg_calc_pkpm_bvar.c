// Tests of the direct nodal magnetic field unit vector (dg_calc_pkpm_bvar) and the
// nodal theta = T_perp/m quotient (dg_calc_pkpm_vars_theta) of the refactored PKPM solver.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_ten_moment.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Smooth, non-trivial field with all three components varying in every direction.
static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.5 * cos(x);
  fout[1] = -0.3 * sin(y);
  fout[2] = 0.1 * cos(z);
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.2 * y;
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.1 * x;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

static void
eval_moms(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * z; // rho
  fout[1] = 1.0 + 0.3 * cos(x); // p_par
  fout[2] = 0.7 + 0.25 * cos(y) * sin(x) + 0.1 * z * z; // p_perp
  fout[3] = 0.0; // M1
}

// Fill the 2^cdim Gauss-Legendre nodes of the reference cell.
static int
gl2_nodes(int cdim, double nodes[8][3])
{
  double g = 1.0 / sqrt(3.0);
  int nn = 1 << cdim;
  for (int n = 0; n < nn; ++n) {
    for (int d = 0; d < cdim; ++d) {
      nodes[n][d] = ((n >> d) & 1) ? g : -g;
    }
  }
  return nn;
}

static void
test_bvar(int cdim)
{
  int poly_order = 1;
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {16, 6, 4};
  int ghost[3] = {1, 1, 1};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  int nb = basis.num_basis;

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *em = mkarr(8 * nb, local_ext.volume);
  struct gkyl_array *bvar = mkarr(3 * nb, local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, em);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_dg_calc_pkpm_bvar *up = gkyl_dg_calc_pkpm_bvar_new(&grid, &basis, false);
  gkyl_dg_calc_pkpm_bvar_advance(up, &local_ext, em, bvar);
  TEST_CHECK(gkyl_dg_calc_pkpm_bvar_num_null_b(up) == 0);

  // Legacy unit vector for a closeness check where its fallback is inactive.
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(1.0, 0.0, 0.0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);
  struct gkyl_dg_calc_em_vars *calc_legacy =
    gkyl_dg_calc_em_vars_new(&grid, &basis, &local_ext, maxwell, geom, 0.0, 0, false);
  int nsurf = nb / (poly_order + 1);
  struct gkyl_array *bvar_legacy = mkarr(9 * nb, local_ext.volume);
  struct gkyl_array *bvar_surf_legacy = mkarr(2 * cdim * 4 * nsurf, local_ext.volume);
  struct gkyl_array *cell_avg_magB2 = gkyl_array_new(GKYL_INT, 1, local_ext.volume);
  gkyl_dg_calc_em_vars_advance(calc_legacy, em, cell_avg_magB2, bvar_legacy, bvar_surf_legacy);

  double nodes[8][3];
  int nn = gl2_nodes(cdim, nodes);
  double max_legacy_diff = 0.0;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *em_d = gkyl_array_cfetch(em, loc);
    const double *b_d = gkyl_array_cfetch(bvar, loc);
    const int *fallback = gkyl_array_cfetch(cell_avg_magB2, loc);
    const double *bl_d = gkyl_array_cfetch(bvar_legacy, loc);

    for (int n = 0; n < nn; ++n) {
      double B[3], b[3];
      for (int i = 0; i < 3; ++i) {
        B[i] = basis.eval_expand(nodes[n], &em_d[(3 + i) * nb]);
        b[i] = basis.eval_expand(nodes[n], &b_d[i * nb]);
      }
      double magB = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
      // b is the interpolant of B/|B| at the construction nodes: exact there.
      for (int i = 0; i < 3; ++i) {
        TEST_CHECK(fabs(b[i] - B[i] / magB) < 1e-13);
        TEST_MSG("b[%d] = %.17g, B/|B| = %.17g", i, b[i], B[i] / magB);
      }
      TEST_CHECK(fabs(b[0] * b[0] + b[1] * b[1] + b[2] * b[2] - 1.0) < 1e-13);
    }

    if (fallback[0] == 0) {
      for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < nb; ++k) {
          double diff = fabs(b_d[i * nb + k] - bl_d[i * nb + k]);
          max_legacy_diff = diff > max_legacy_diff ? diff : max_legacy_diff;
        }
      }
    }
  }
  // The legacy vector is sqrt-with-sign of a projected bb; both are second-order
  // accurate representations of the same smooth field, so they agree to the
  // projection error of the grid but not to round-off.
  TEST_CHECK(max_legacy_diff < 2.0e-2);
  TEST_MSG("max |b_new - b_legacy| = %g", max_legacy_diff);

  // Null-field policy: B = 0 in one cell gives b = 0 there and counts every node.
  int idx0[3] = {local.lower[0], cdim > 1 ? local.lower[1] : 0, cdim > 2 ? local.lower[2] : 0};
  long loc0 = gkyl_range_idx(&local_ext, idx0);
  double *em0 = gkyl_array_fetch(em, loc0);
  for (int k = 3 * nb; k < 6 * nb; ++k) {
    em0[k] = 0.0;
  }
  gkyl_dg_calc_pkpm_bvar_advance(up, &local, em, bvar);
  TEST_CHECK(gkyl_dg_calc_pkpm_bvar_num_null_b(up) == nn);
  const double *b0 = gkyl_array_cfetch(bvar, loc0);
  for (int k = 0; k < 3 * nb; ++k) {
    TEST_CHECK(b0[k] == 0.0);
  }

  gkyl_array_release(cell_avg_magB2);
  gkyl_array_release(bvar_surf_legacy);
  gkyl_array_release(bvar_legacy);
  gkyl_dg_calc_em_vars_release(calc_legacy);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(maxwell);
  gkyl_dg_calc_pkpm_bvar_release(up);
  gkyl_array_release(bvar);
  gkyl_array_release(em);
}

static void
test_theta(int cdim)
{
  int poly_order = 1;
  double lower[3] = {-1.0, -1.0, -1.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {12, 6, 4};
  int ghost[3] = {1, 1, 1};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  int nb = basis.num_basis;

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_array *moms = mkarr(4 * nb, local_ext.volume);
  struct gkyl_array *theta = mkarr(nb, local_ext.volume);
  struct gkyl_array *theta_weak = mkarr(nb, local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 8, 4, eval_moms, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &local_ext, moms);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_wv_eqn *eqn = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&grid, &local_ext, 0, 0, false);
  struct gkyl_dg_calc_pkpm_vars *up =
    gkyl_dg_calc_pkpm_vars_new(&grid, &basis, &local_ext, eqn, geom, 0.0, false);

  int num_null = gkyl_dg_calc_pkpm_vars_theta(up, &local_ext, moms, theta);
  TEST_CHECK(num_null == 0);

  // Weak division p_perp/rho (same answer for a positive p=1 density).
  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_new(local_ext.volume, nb);
  gkyl_dg_div_op_range(mem, &basis, 0, theta_weak, 2, moms, 0, moms, &local_ext);
  gkyl_dg_bin_op_mem_release(mem);

  double nodes[8][3];
  int nn = gl2_nodes(cdim, nodes);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local_ext);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&local_ext, iter.idx);
    const double *m_d = gkyl_array_cfetch(moms, loc);
    const double *th_d = gkyl_array_cfetch(theta, loc);
    const double *thw_d = gkyl_array_cfetch(theta_weak, loc);
    for (int n = 0; n < nn; ++n) {
      double rho = basis.eval_expand(nodes[n], &m_d[0]);
      double p_perp = basis.eval_expand(nodes[n], &m_d[2 * nb]);
      double th = basis.eval_expand(nodes[n], th_d);
      TEST_CHECK(fabs(th - p_perp / rho) < 1e-13);
    }
    for (int k = 0; k < nb; ++k) {
      TEST_CHECK(fabs(th_d[k] - thw_d[k]) < 1e-12);
      TEST_MSG("theta[%d] = %.17g, weak = %.17g", k, th_d[k], thw_d[k]);
    }
  }

  // Non-positive density policy: rho = 0 in one cell gives theta = 0 and counts every node.
  int idx0[3] = {local.lower[0], cdim > 1 ? local.lower[1] : 0, cdim > 2 ? local.lower[2] : 0};
  long loc0 = gkyl_range_idx(&local_ext, idx0);
  double *m0 = gkyl_array_fetch(moms, loc0);
  for (int k = 0; k < nb; ++k) {
    m0[k] = 0.0;
  }
  num_null = gkyl_dg_calc_pkpm_vars_theta(up, &local, moms, theta);
  TEST_CHECK(num_null == nn);
  const double *th0 = gkyl_array_cfetch(theta, loc0);
  for (int k = 0; k < nb; ++k) {
    TEST_CHECK(th0[k] == 0.0);
  }

  gkyl_dg_calc_pkpm_vars_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(eqn);
  gkyl_array_release(theta_weak);
  gkyl_array_release(theta);
  gkyl_array_release(moms);
}

void test_bvar_1x() { test_bvar(1); }
void test_bvar_2x() { test_bvar(2); }
void test_bvar_3x() { test_bvar(3); }
void test_theta_1x() { test_theta(1); }
void test_theta_2x() { test_theta(2); }
void test_theta_3x() { test_theta(3); }

TEST_LIST = {
  {"test_bvar_1x", test_bvar_1x},
  {"test_bvar_2x", test_bvar_2x},
  {"test_bvar_3x", test_bvar_3x},
  {"test_theta_1x", test_theta_1x},
  {"test_theta_2x", test_theta_2x},
  {"test_theta_3x", test_theta_3x},
  {NULL, NULL},
};
