#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_twofluid.h>
#include <gkyl_wv_gr_twofluid_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_twofluid_basic_minkowski()
{
  double gas_gamma_elc = 5.0 / 3.0;
  double gas_gamma_ion = 5.0 / 3.0;
  double mass_ion = 1.0;
  double charge_ion = 1.0;
  double mass_elc = 1.0 / 1836.2;
  double charge_elc = -1.0;

  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_twofluid = gkyl_wv_gr_twofluid_new(mass_elc, mass_ion, charge_elc, charge_ion, gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_twofluid->num_equations == 84 );
  TEST_CHECK( gr_twofluid->num_waves == 6 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_ion = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;
      double rho_elc = rho_ion * mass_elc / mass_ion;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double he = 1.0 + ((p / rho_elc) * (gas_gamma_elc / (gas_gamma_elc - 1.0)));
      double hi = 1.0 + ((p / rho_ion) * (gas_gamma_ion / (gas_gamma_ion - 1.0)));

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double q[84];
      q[0] = sqrt(spatial_det) * rho_elc * W;
      q[1] = sqrt(spatial_det) * rho_elc * he * (W * W) * u;
      q[2] = sqrt(spatial_det) * rho_elc * he * (W * W) * v;
      q[3] = sqrt(spatial_det) * rho_elc * he * (W * W) * w;
      q[4] = sqrt(spatial_det) * ((rho_elc * he * (W * W)) - p - (rho_elc * W));

      q[5] = sqrt(spatial_det) * rho_ion * W;
      q[6] = sqrt(spatial_det) * rho_ion * hi * (W * W) * u;
      q[7] = sqrt(spatial_det) * rho_ion * hi * (W * W) * v;
      q[8] = sqrt(spatial_det) * rho_ion * hi * (W * W) * w;
      q[9] = sqrt(spatial_det) * ((rho_ion * hi * (W * W)) - p - (rho_ion * W));

      q[10] = Dx; q[11] = Dy; q[12] = Dz;
      q[13] = Bx; q[14] = By; q[15] = Bz;
      q[16] = phi; q[17] = psi;

      q[18] = lapse;
      q[19] = shift[0]; q[20] = shift[1]; q[21] = shift[2];

      q[22] = spatial_metric[0][0]; q[23] = spatial_metric[0][1]; q[24] = spatial_metric[0][2];
      q[25] = spatial_metric[1][0]; q[26] = spatial_metric[1][1]; q[27] = spatial_metric[1][2];
      q[28] = spatial_metric[2][0]; q[29] = spatial_metric[2][1]; q[30] = spatial_metric[2][2];

      q[31] = extrinsic_curvature[0][0]; q[32] = extrinsic_curvature[0][1]; q[33] = extrinsic_curvature[0][2];
      q[34] = extrinsic_curvature[1][0]; q[35] = extrinsic_curvature[1][1]; q[36] = extrinsic_curvature[1][2];
      q[37] = extrinsic_curvature[2][0]; q[38] = extrinsic_curvature[2][1]; q[39] = extrinsic_curvature[2][2];

      q[40] = 1.0;

      q[41] = lapse_der[0]; q[42] = lapse_der[1]; q[43] = lapse_der[2];
      q[44] = shift_der[0][0]; q[45] = shift_der[0][1]; q[46] = shift_der[0][2];
      q[47] = shift_der[1][0]; q[48] = shift_der[1][1]; q[49] = shift_der[1][2];
      q[50] = shift_der[2][0]; q[51] = shift_der[2][1]; q[52] = shift_der[2][2];

      q[53] = spatial_metric_der[0][0][0]; q[54] = spatial_metric_der[0][0][1]; q[55] = spatial_metric_der[0][0][2];
      q[56] = spatial_metric_der[0][1][0]; q[57] = spatial_metric_der[0][1][1]; q[58] = spatial_metric_der[0][1][2];
      q[59] = spatial_metric_der[0][2][0]; q[60] = spatial_metric_der[0][2][1]; q[61] = spatial_metric_der[0][2][2];

      q[62] = spatial_metric_der[1][0][0]; q[63] = spatial_metric_der[1][0][1]; q[64] = spatial_metric_der[1][0][2];
      q[65] = spatial_metric_der[1][1][0]; q[66] = spatial_metric_der[1][1][1]; q[67] = spatial_metric_der[1][1][2];
      q[68] = spatial_metric_der[1][2][0]; q[69] = spatial_metric_der[1][2][1]; q[70] = spatial_metric_der[1][2][2];

      q[71] = spatial_metric_der[2][0][0]; q[72] = spatial_metric_der[2][0][1]; q[73] = spatial_metric_der[2][0][2];
      q[74] = spatial_metric_der[2][1][0]; q[75] = spatial_metric_der[2][1][1]; q[76] = spatial_metric_der[2][1][2];
      q[77] = spatial_metric_der[2][2][0]; q[78] = spatial_metric_der[2][2][1]; q[79] = spatial_metric_der[2][2][2];

      q[80] = 0.0;
      q[81] = x; q[82] = y; q[83] = 0.0;

      double prims[84];
      gkyl_gr_twofluid_prim_vars(gas_gamma_elc, gas_gamma_ion, q, prims);
      
      TEST_CHECK( gkyl_compare(prims[0], rho_elc, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[1], u, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[2], v, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[3], w, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[4], p, 1e-12) );

      TEST_CHECK( gkyl_compare(prims[5], rho_ion, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[6], u, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[7], v, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[8], w, 1e-12) );
      TEST_CHECK( gkyl_compare(prims[9], p, 1e-12) );

      double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
      double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
      double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

      double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
      double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
      double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

      double fluxes[3][18] = {
        { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
          (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])),
          e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
          -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
        { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
          (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])),
          -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
          -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
        { (lapse * sqrt(spatial_det)) * (rho_elc * W * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_elc * he * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (((rho_elc * he * (W * W)) - p - (rho_elc * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
          (lapse * sqrt(spatial_det)) * (rho_ion * W * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho_ion * hi * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (((rho_ion * hi * (W * W)) - p - (rho_ion * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])),
          (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
          e_fact * Dz, b_fact * (light_speed * light_speed) * Bz },
      };

      double norm[3][3] = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
      };

      double tau1[3][3] = {
        { 0.0, 1.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
      };

      double tau2[3][3] = {
        { 0.0, 0.0, 1.0 },
        { 0.0, 0.0, -1.0 },
        { 0.0, 1.0, 0.0 },
      };

      double q_local[84], flux_local[84], flux[84];
      for (int d = 0; d < 3; d++) {
        gr_twofluid->rotate_to_local_func(gr_twofluid, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_twofluid_flux(gas_gamma_elc, gas_gamma_ion, light_speed, e_fact, b_fact, q_local, flux_local);
        gr_twofluid->rotate_to_global_func(gr_twofluid, tau1[d], tau2[d], norm[d], flux_local, flux);

        for (int i = 0; i < 18; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-8) );
        }
      }

      double q_l[84], q_g[84];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_twofluid, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_twofluid, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 18; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
        }

        double w1[84], q1[84];
        gr_twofluid->cons_to_riem(gr_twofluid, q_local, q_local, w1);
        gr_twofluid->riem_to_cons(gr_twofluid, q_local, w1, q1);

        for (int i = 0; i < 18; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_twofluid);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_twofluid_basic_minkowski", test_gr_twofluid_basic_minkowski },
  { NULL, NULL },
};