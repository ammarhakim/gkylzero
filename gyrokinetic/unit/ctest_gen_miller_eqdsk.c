#include <acutest.h>
#include <gkyl_gen_miller_eqdsk.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_efit.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_position_map.h>
#include <gkyl_tok_geo.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> // mkstemp(): create a uniquely named temporary file.
#include <unistd.h> // unlink(): remove a file's directory entry.

static double
constant_q(double r_omp, void *ctx)
{
  return 2;
}

static double
ltx_q(double r_omp, void *ctx)
{
  double radius = r_omp - 0.406051632, a = 0.61183 - 0.406051632;
  return 4.3153848 * exp(0.5 * 2.6899871 * (radius * radius / (a * a) - 1));
}

static double
d3d_q(double r_omp, void *ctx)
{
  return ((49.46395467479657 * r_omp - 260.79513158768754) * r_omp + 458.42618139184754) * r_omp -
         267.63441353752336;
}

static double
tcv_q(double r_omp, void *ctx)
{
  return ((497.3420166252413 * r_omp - 1408.736172826569) * r_omp + 1331.4134861681464) * r_omp -
         419.00692601227627;
}

static double
invalid_q(double r_omp, void *ctx)
{
  // Parse at runtime so -ffast-math cannot fold a constant NaN return.
  return strtod("nan", NULL);
}

static void
test_circle(void)
{
  struct gkyl_gen_miller_eqdsk_inp inp = {
    .r_axis = 2,
    .b_axis = 1.5,
    .r_lcfs = 2.5,
    .qprofile = constant_q,
  };

  gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inp);
  TEST_ASSERT(up != NULL);

  const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(up);

  for (int j = 0; j <= 100; ++j) {
    double radius = info->r_guard * j / 100;
    double exact = 1.5 * radius * radius / (2 + sqrt(4 - radius * radius));
    double flux = gkyl_gen_miller_eqdsk_psi_at_R_omp(up, 2 + radius);
    TEST_CHECK(fabs(flux - exact) < 1e-8 * info->psi_boundary);
    TEST_CHECK(fabs(gkyl_gen_miller_eqdsk_R_omp_at_psi(up, flux) - (2 + radius)) < 1e-8);
  }

  TEST_CHECK(gkyl_gen_miller_eqdsk_advance(up));

  gkyl_gen_miller_eqdsk_release(up);

  inp.qprofile = invalid_q;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);

  inp.qprofile = constant_q;
  inp.a_shift = 10;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);

  inp.a_shift = 0;
  inp.nr = 64;
  TEST_CHECK(gkyl_gen_miller_eqdsk_inew(&inp) == NULL);
}

static void
check_reader(
  const char *filename, const struct gkyl_gen_miller_eqdsk_inp *inp,
  const gkyl_gen_miller_eqdsk *generator
)
{
  struct gkyl_efit_inp efit_inp = {.rz_poly_order = 2, .flux_poly_order = 1, .reflect = false};
  snprintf(efit_inp.filepath, sizeof efit_inp.filepath, "%s", filename);

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);

  const struct gkyl_gen_miller_eqdsk_info *info = gkyl_gen_miller_eqdsk_get_info(generator);

  TEST_CHECK(fabs(efit->sibry / info->psi_boundary - 1) < 1e-8);
  TEST_CHECK(fabs(efit->zmaxis - inp->z_axis) < 1e-9);

  // The axis is a degenerate critical point. Check any reported X-points lie
  // outside the physical simulation annulus, rather than rejecting the axis.
  for (int j = 0; j < efit->num_xpts_cubic; ++j) {
    double point[2] = {efit->Rxpt_cubic[j], efit->Zxpt_cubic[j]}, value[3];
    efit->evf->eval_cubic_wgrad(0, point, value, efit->evf->ctx);
    double lower = gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, inp->omp_bounds[0]);
    double upper = gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, inp->omp_bounds[1]);
    TEST_CHECK(value[0] < lower || value[0] > upper);
  }

  double max_bp_error = 0, max_flux_error = 0, max_q_error = 0;
  for (int surface = 0; surface < 5; ++surface) {
    double r_omp = inp->omp_bounds[0] + (inp->omp_bounds[1] - inp->omp_bounds[0]) * surface / 4;
    double flux = gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, r_omp);
    double offset = r_omp - inp->r_axis;
    double radius = 2 * offset / (1 + sqrt(1 - 2 * inp->a_shift * offset / inp->r_axis));
    double h = radius * 1e-5;
    double omp_lo =
      inp->r_axis - inp->a_shift * (radius - h) * (radius - h) / (2 * inp->r_axis) + radius - h;
    double omp_hi =
      inp->r_axis - inp->a_shift * (radius + h) * (radius + h) / (2 * inp->r_axis) + radius + h;
    double dpsi = (gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, omp_hi) -
                   gkyl_gen_miller_eqdsk_psi_at_R_omp(generator, omp_lo)) /
                  (2 * h);
    double q_integral = 0;
    for (int j = 0; j < 1024; ++j) {
      double theta = 2 * M_PI * j / 1024, rz[2], value[3];
      gkyl_gen_miller_eqdsk_RZ_at_rtheta(generator, radius, theta, rz);
      efit->evf->eval_cubic_wgrad(0, rz, value, efit->evf->ctx);
      double eta = asin(inp->delta), angle = theta + eta * sin(theta);
      double Rr = cos(angle) - inp->a_shift * radius / inp->r_axis;
      double Rt = -radius * (1 + eta * cos(theta)) * sin(angle);
      double Zr = inp->kappa * sin(theta), Zt = inp->kappa * radius * cos(theta);
      double dl = hypot(Rt, Zt), bp = hypot(value[1], value[2]) / rz[0];
      double expected_bp = dpsi * dl / (rz[0] * (Rr * Zt - Rt * Zr));
      max_bp_error = fmax(max_bp_error, fabs(bp / expected_bp - 1));
      max_flux_error = fmax(max_flux_error, fabs(value[0] - flux) / info->psi_boundary);
      q_integral += inp->r_axis * inp->b_axis * dl / (rz[0] * rz[0] * bp) / 1024;
    }
    max_q_error = fmax(max_q_error, fabs(q_integral / inp->qprofile(r_omp, inp->qctx) - 1));
  }

  TEST_CHECK(max_bp_error < 1e-3);
  TEST_MSG("Bp error = %.9g", max_bp_error);
  TEST_CHECK(max_flux_error < 1e-5);
  TEST_MSG("flux error = %.9g", max_flux_error);
  TEST_CHECK(max_q_error < 1e-3);
  TEST_MSG("q error = %.9g", max_q_error);

  gkyl_efit_release(efit);

  // Check roots of the actual quadratic and cubic EFIT reconstructions.
  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_GEOMETRY_TOKAMAK_IWL,
    .rleft = info->rz_bounds[0],
    .rright = info->rz_bounds[1],
    .rmin = info->rz_bounds[0],
    .rmax = info->rz_bounds[1],
    .zmin = info->rz_bounds[2],
    .zmax = info->rz_bounds[3],
  };

  for (int cubic = 0; cubic < 2; ++cubic) {
    ginp.use_cubics = cubic;
    ginp.use_hyperbolic_numbers = true;
    struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&efit_inp, &ginp);
    geo->rmin = ginp.rmin;
    geo->rmax = ginp.rmax;
    double R[4], slope[4], dR[4], dZ[4];
    int roots = gkyl_tok_geo_R_psiZ(geo, info->psi_boundary, inp->z_axis, 4, R, slope, dR, dZ);
    TEST_CHECK(roots == 2);
    double error = 1e9;
    for (int j = 0; j < roots; ++j) {
      error = fmin(error, fabs(R[j] - inp->r_lcfs));
    }
    TEST_CHECK(error < 1e-4 * info->r_lcfs);
    gkyl_tok_geo_release(geo);
  }

  ginp.use_cubics = false;
  struct gkyl_gk_geometry_inp geometry = {
    .geometry_id = GKYL_GEOMETRY_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .has_LCFS = true,
    .x_LCFS = info->psi_boundary,
  };
  double lower[] = {0.95 * info->psi_boundary, -0.01, -M_PI + 1e-10};
  double upper[] = {1.05 * info->psi_boundary, 0.01, M_PI - 1e-10};
  int cells[] = {2, 1, 2}, ghosts[] = {1, 1, 1};

  gkyl_rect_grid_init(&geometry.grid, 3, lower, upper, cells);
  gkyl_cart_modal_serendip(&geometry.basis, 3, 1);
  gkyl_create_grid_ranges(&geometry.grid, ghosts, &geometry.local_ext, &geometry.local);
  geometry.global = geometry.local;
  geometry.global_ext = geometry.local_ext;
  geometry.geo_grid = geometry.grid;
  geometry.geo_basis = geometry.basis;
  geometry.geo_local = geometry.local;
  geometry.geo_local_ext = geometry.local_ext;
  geometry.geo_global = geometry.global;
  geometry.geo_global_ext = geometry.global_ext;
  geometry.position_map = gkyl_position_map_null_new();

  struct gk_geometry *gk_geo = gkyl_gk_geometry_tok_new(&geometry);
  TEST_CHECK(gkyl_gk_geometry_reduce_bmag(gk_geo, GKYL_MIN) > 0);

  gkyl_gk_geometry_release(gk_geo);
  gkyl_position_map_release(geometry.position_map);
}

static void
test_devices(void)
{
  const char *names[] = {"ltx", "d3d", "tcv"};
  struct gkyl_gen_miller_eqdsk_inp inputs[] = {
    {
      .r_axis = 0.406051632,
      .b_axis = 0.240606108,
      .r_lcfs = 0.61183,
      .kappa = 1.3,
      .delta = 0.4,
      .qprofile = ltx_q,
      .has_omp_bounds = true,
      .omp_bounds = {0.58, 0.7},
    },
    {
      .r_axis = 1.6,
      .z_axis = 0.013055028,
      .b_axis = 2 * 1.6486461 / 1.6,
      .r_lcfs = 2.17,
      .kappa = 1.35,
      .delta = 0.4,
      .qprofile = d3d_q,
      .has_omp_bounds = true,
      .omp_bounds = {2.17 - 5 * 0.15 / 8, 2.17 + 3 * 0.15 / 8},
    },
    {
      .r_axis = 0.8867856264,
      .z_axis = 0.1414361745,
      .b_axis = 1.4,
      .r_lcfs = 1.0870056099999,
      .kappa = 1.4,
      .delta = -0.38,
      .a_shift = 0.5,
      .qprofile = tcv_q,
      .has_omp_bounds = true,
      .omp_bounds = {1.0513586729299085, 1.1555930803275622},
    }
  };

  // Independent Python values at the same tolerance and physical parameters.
  double expected_boundary[] = {0.003518612110703164, 0.28091394499783906, 0.03483576483233173};
  for (int device = 0; device < 3; ++device) {
    char filename[] = "/tmp/gkyl_miller_XXXXXX";
    int descriptor = mkstemp(filename);
    TEST_ASSERT(descriptor >= 0);
    close(descriptor);
    gkyl_gen_miller_eqdsk *up = gkyl_gen_miller_eqdsk_inew(&inputs[device]);

    TEST_ASSERT(up != NULL);
    TEST_CHECK(
      fabs(gkyl_gen_miller_eqdsk_get_info(up)->psi_boundary / expected_boundary[device] - 1) < 1e-8
    );
    TEST_ASSERT(gkyl_gen_miller_eqdsk_advance(up));

    TEST_ASSERT(gkyl_gen_miller_eqdsk_write(up, filename));
    check_reader(filename, &inputs[device], up);
    unlink(filename);
    // Optional Python-generated files exercise the same native reader checks.
    const char *references = getenv("GKYL_MILLER_REFERENCE_DIR");
    if (references) {
      char path[1024];
      snprintf(path, sizeof path, "%s/miller_%s.geqdsk", references, names[device]);
      check_reader(path, &inputs[device], up);
    }
    gkyl_gen_miller_eqdsk_release(up);
  }
}

TEST_LIST = {
  {"miller_circle_and_invalid_inputs", test_circle},
  {"miller_devices_and_reader", test_devices},
  {NULL, NULL}
};
