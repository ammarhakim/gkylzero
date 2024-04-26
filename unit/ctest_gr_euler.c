#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_mild_primitive()
{
  double *q = gkyl_malloc(sizeof(double[29]));
  double *prims = gkyl_malloc(sizeof(double[29]));
  double gas_gamma = 4.0 / 3.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, spacetime, false);

  double rho = 2.0;
  double u = 0.1;
  double v = 0.2;
  double w = 0.3;
  double p = 0.01;

  double v_sq = (u * u) + (v * v) + (w * w);
  double W = 1.0 / sqrt(1.0 - v_sq);
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(pow(10.0, -8.0));
  }

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

  q[0] = rho * W;
  q[1] = rho * h * (W * W) * u;
  q[2] = rho * h * (W * W) * v;
  q[3] = rho * h * (W * W) * w;
  q[4] = rho * h * (W * W) - p - (rho * W);
  
  q[5] = 1.0;
  q[6] = 1.0;
  q[7] = 0.0; q[8] = 0.0; q[9] = 0.0;

  q[10] = 1.0; q[11] = 0.0, q[12] = 0.0;
  q[13] = 0.0; q[14] = 1.0; q[15] = 0.0;
  q[16] = 0.0; q[17] = 0.0; q[18] = 1.0;

  q[19] = 1.0; q[20] = 0.0, q[21] = 0.0;
  q[22] = 0.0; q[23] = 1.0; q[24] = 0.0;
  q[25] = 0.0; q[26] = 0.0; q[27] = 1.0;

  q[28] = 1.0;
  
  gkyl_gr_euler_prim_vars(gas_gamma, q, prims);
  
  TEST_CHECK( gkyl_compare(prims[0], rho, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[1], u, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[2], v, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[3], w, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[4], p, 1e-6) );

  gkyl_free(q);
  gkyl_free(prims);

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_strong_primitive()
{
  double *q = gkyl_malloc(sizeof(double[29]));
  double *prims = gkyl_malloc(sizeof(double[29]));
  double gas_gamma = 4.0 / 3.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, spacetime, false);

  double rho = 3.0;
  double u = 0.4;
  double v = 0.5;
  double w = 0.6;
  double p = 0.001;

  double v_sq = (u * u) + (v * v) + (w * w);
  double W = 1.0 / sqrt(1.0 - v_sq);
  if (v_sq > 1.0 - pow(10.0, -8.0)) {
    W = 1.0 / sqrt(pow(10.0, -8.0));
  }

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

  q[0] = rho * W;
  q[1] = rho * h * (W * W) * u;
  q[2] = rho * h * (W * W) * v;
  q[3] = rho * h * (W * W) * w;
  q[4] = rho * h * (W * W) - p - (rho * W);
  
  q[5] = 1.0;
  q[6] = 1.0;
  q[7] = 0.0; q[8] = 0.0; q[9] = 0.0;

  q[10] = 1.0; q[11] = 0.0, q[12] = 0.0;
  q[13] = 0.0; q[14] = 1.0; q[15] = 0.0;
  q[16] = 0.0; q[17] = 0.0; q[18] = 1.0;

  q[19] = 1.0; q[20] = 0.0, q[21] = 0.0;
  q[22] = 0.0; q[23] = 1.0; q[24] = 0.0;
  q[25] = 0.0; q[26] = 0.0; q[27] = 1.0;

  q[28] = 1.0;
  
  gkyl_gr_euler_prim_vars(gas_gamma, q, prims);
  
  TEST_CHECK( gkyl_compare(prims[0], rho, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[1], u, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[2], v, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[3], w, 1e-6) );
  TEST_CHECK( gkyl_compare(prims[4], p, 1e-6) );

  gkyl_free(q);
  gkyl_free(prims);

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void test_gr_blackhole_mild_primitive()
{
  double gas_gamma = 5.0 / 3.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, spacetime, false);

  double rho = 2.0;
  double u = 0.1;
  double v = 0.2;
  double w = 0.3;
  double p = 0.01;

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 1; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double *q = gkyl_malloc(sizeof(double[29]));
      double *prims = gkyl_malloc(sizeof(double[29]));

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);

      double *vel = gkyl_malloc(sizeof(double[3]));
      double v_sq = 0.0;
      vel[0] = u; vel[1] = v; vel[2] = w;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      if (v_sq > 1.0 - pow(10.0, -8.0)) {
        W = 1.0 / sqrt(pow(10.0, -8.0));
      }

      q[0] = sqrt(spatial_det) * rho * W;
      q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
      q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
      q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
      q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

      q[5] = spatial_det;
      q[6] = lapse;
      q[7] = shift[0]; q[8] = shift[1]; q[9] = shift[2];

      q[10] = spatial_metric[0][0]; q[11] = spatial_metric[0][1]; q[12] = spatial_metric[0][2];
      q[13] = spatial_metric[1][0]; q[14] = spatial_metric[1][1]; q[15] = spatial_metric[1][2];
      q[16] = spatial_metric[2][0]; q[17] = spatial_metric[2][1]; q[18] = spatial_metric[2][2];

      q[19] = inv_spatial_metric[0][0]; q[20] = inv_spatial_metric[0][1]; q[21] = inv_spatial_metric[0][2];
      q[22] = inv_spatial_metric[1][0]; q[23] = inv_spatial_metric[1][1]; q[24] = inv_spatial_metric[1][2];
      q[25] = inv_spatial_metric[2][0]; q[26] = inv_spatial_metric[2][1]; q[27] = inv_spatial_metric[2][2];

      if (in_excision_region) {
        for (int i = 0; i < 28; i++) {
          q[i] = 0.0;
        }

        q[28] = -1.0;
      }
      else {
        q[28] = 1.0;
      }

      gkyl_gr_euler_prim_vars(gas_gamma, q, prims);

      if (in_excision_region) {
        for (int i = 0; i < 28; i++) {
          TEST_CHECK( gkyl_compare(prims[i], 0.0, 1e-10) );
        }
        TEST_CHECK( gkyl_compare(prims[28], -1.0, 1e-10) );
      }
      else {
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-1) );

        TEST_CHECK( gkyl_compare(prims[5], spatial_det, 1e-10) );
        TEST_CHECK( gkyl_compare(prims[6], lapse, 1e-10) );
        TEST_CHECK( gkyl_compare(prims[7], shift[0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[8], shift[1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[9], shift[2], 1e-10) );

        TEST_CHECK( gkyl_compare(prims[10], spatial_metric[0][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[11], spatial_metric[0][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[12], spatial_metric[0][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[13], spatial_metric[1][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[14], spatial_metric[1][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[15], spatial_metric[1][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[16], spatial_metric[2][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[17], spatial_metric[2][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[18], spatial_metric[2][2], 1e-10) );

        TEST_CHECK( gkyl_compare(prims[19], inv_spatial_metric[0][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[20], inv_spatial_metric[0][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[21], inv_spatial_metric[0][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[22], inv_spatial_metric[1][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[23], inv_spatial_metric[1][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[24], inv_spatial_metric[1][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[25], inv_spatial_metric[2][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[26], inv_spatial_metric[2][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[27], inv_spatial_metric[2][2], 1e-10) );

        TEST_CHECK( gkyl_compare(prims[28], 1.0, 1e-10) );
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);
      gkyl_free(shift);
      gkyl_free(vel);

      gkyl_free(q);
      gkyl_free(prims);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void test_gr_blackhole_strong_primitive()
{
  double gas_gamma = 5.0 / 3.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, spacetime, false);

  double rho = 3.0;
  double u = 0.2;
  double v = 0.3;
  double w = 0.4;
  double p = 0.001;

  double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 1; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double *q = gkyl_malloc(sizeof(double[29]));
      double *prims = gkyl_malloc(sizeof(double[29]));

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);

      double *vel = gkyl_malloc(sizeof(double[3]));
      double v_sq = 0.0;
      vel[0] = u; vel[1] = v; vel[2] = w;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      if (v_sq > 1.0 - pow(10.0, -8.0)) {
        W = 1.0 / sqrt(pow(10.0, -8.0));
      }

      q[0] = sqrt(spatial_det) * rho * W;
      q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
      q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
      q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
      q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

      q[5] = spatial_det;
      q[6] = lapse;
      q[7] = shift[0]; q[8] = shift[1]; q[9] = shift[2];

      q[10] = spatial_metric[0][0]; q[11] = spatial_metric[0][1]; q[12] = spatial_metric[0][2];
      q[13] = spatial_metric[1][0]; q[14] = spatial_metric[1][1]; q[15] = spatial_metric[1][2];
      q[16] = spatial_metric[2][0]; q[17] = spatial_metric[2][1]; q[18] = spatial_metric[2][2];

      q[19] = inv_spatial_metric[0][0]; q[20] = inv_spatial_metric[0][1]; q[21] = inv_spatial_metric[0][2];
      q[22] = inv_spatial_metric[1][0]; q[23] = inv_spatial_metric[1][1]; q[24] = inv_spatial_metric[1][2];
      q[25] = inv_spatial_metric[2][0]; q[26] = inv_spatial_metric[2][1]; q[27] = inv_spatial_metric[2][2];

      if (in_excision_region) {
        for (int i = 0; i < 28; i++) {
          q[i] = 0.0;
        }

        q[28] = -1.0;
      }
      else {
        q[28] = 1.0;
      }

      gkyl_gr_euler_prim_vars(gas_gamma, q, prims);

      if (in_excision_region) {
        for (int i = 0; i < 28; i++) {
          TEST_CHECK( gkyl_compare(prims[i], 0.0, 1e-10) );
        }
        TEST_CHECK( gkyl_compare(prims[28], -1.0, 1e-10) );
      }
      else {
        TEST_CHECK( gkyl_compare(prims[0], rho, 3e-1) );
        TEST_CHECK( gkyl_compare(prims[1], u, 3e-1) );
        TEST_CHECK( gkyl_compare(prims[2], v, 3e-1) );
        TEST_CHECK( gkyl_compare(prims[3], w, 3e-1) );
        TEST_CHECK( gkyl_compare(prims[4], p, 3e-1) );

        TEST_CHECK( gkyl_compare(prims[5], spatial_det, 1e-10) );
        TEST_CHECK( gkyl_compare(prims[6], lapse, 1e-10) );
        TEST_CHECK( gkyl_compare(prims[7], shift[0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[8], shift[1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[9], shift[2], 1e-10) );

        TEST_CHECK( gkyl_compare(prims[10], spatial_metric[0][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[11], spatial_metric[0][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[12], spatial_metric[0][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[13], spatial_metric[1][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[14], spatial_metric[1][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[15], spatial_metric[1][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[16], spatial_metric[2][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[17], spatial_metric[2][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[18], spatial_metric[2][2], 1e-10) );

        TEST_CHECK( gkyl_compare(prims[19], inv_spatial_metric[0][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[20], inv_spatial_metric[0][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[21], inv_spatial_metric[0][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[22], inv_spatial_metric[1][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[23], inv_spatial_metric[1][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[24], inv_spatial_metric[1][2], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[25], inv_spatial_metric[2][0], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[26], inv_spatial_metric[2][1], 1e-10) );
        TEST_CHECK( gkyl_compare(prims[27], inv_spatial_metric[2][2], 1e-10) );

        TEST_CHECK( gkyl_compare(prims[28], 1.0, 1e-10) );
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);
      gkyl_free(shift);
      gkyl_free(vel);

      gkyl_free(q);
      gkyl_free(prims);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_mild_primitive", test_gr_mild_primitive },
  { "gr_strong_primitive", test_gr_strong_primitive },
  { "gr_blackhole_mild_primitive", test_gr_blackhole_mild_primitive },
  { "gr_blackhole_strong_primitive", test_gr_blackhole_strong_primitive },
  { NULL, NULL },
};