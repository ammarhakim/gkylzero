#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_maxwell_tetrad.h>
#include <gkyl_wv_gr_maxwell_tetrad_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_maxwell_tetrad_basic_minkowski()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, spacetime, false);

  TEST_CHECK( gr_maxwell_tetrad->num_equations == 22 );
  TEST_CHECK( gr_maxwell_tetrad->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);

      double q[22];
      q[0] = Dx; q[1] = Dy; q[2] = Dz;
      q[3] = Bx; q[4] = By; q[5] = Bz;
      q[6] = phi; q[7] = psi;

      q[8] = lapse;
      q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

      q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
      q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
      q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

      q[21] = 1.0;

      double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
      double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
      double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

      double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
      double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
      double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

      double fluxes[3][8] = {
        { e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
          -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
        { -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
          -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
        { (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
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

      double q_local[22], flux_local_sr[22], flux_local_gr[22], flux[22];
      for (int d = 0; d < 3; d++) {
        gr_maxwell_tetrad->rotate_to_local_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, q_local, flux_local_sr);
        gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
        gr_maxwell_tetrad->rotate_to_global_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

        for (int i = 0; i < 8; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-16) );
        }
      }

      double q_l[6], q_g[6];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 6; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
        }

        double w1[6], q1[6];
        gr_maxwell_tetrad->cons_to_riem(gr_maxwell_tetrad, q_local, q_local, w1);
        gr_maxwell_tetrad->riem_to_cons(gr_maxwell_tetrad, q_local, w1, q1);

        for (int i = 0; i < 6; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(shift);
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_basic_schwarzschild()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, spacetime, false);

  TEST_CHECK( gr_maxwell_tetrad->num_equations == 22 );
  TEST_CHECK( gr_maxwell_tetrad->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);

      if (!in_excision_region) {
        double q[22];
        q[0] = Dx; q[1] = Dy; q[2] = Dz;
        q[3] = Bx; q[4] = By; q[5] = Bz;
        q[6] = phi; q[7] = psi;

        q[8] = lapse;
        q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

        q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
        q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
        q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

        q[21] = 1.0;

        double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
        double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
        double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

        double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
        double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
        double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

        double fluxes[3][8] = {
          { e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
            -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
          { -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
            -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
          { (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
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

        double q_local[22], flux_local_sr[22], flux_local_gr[22], flux[22];
        for (int d = 0; d < 3; d++) {
          gr_maxwell_tetrad->rotate_to_local_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, q_local, flux_local_sr);
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
          gr_maxwell_tetrad->rotate_to_global_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-13) );
          }
        }

        double q_l[6], q_g[6];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 6; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-15) );
          }

          double w1[6], q1[6];
          gr_maxwell_tetrad->cons_to_riem(gr_maxwell_tetrad, q_local, q_local, w1);
          gr_maxwell_tetrad->riem_to_cons(gr_maxwell_tetrad, q_local, w1, q1);

          for (int i = 0; i < 6; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(shift);
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_maxwell_tetrad_basic_kerr()
{
  double light_speed = 1.0;
  double e_fact = 0.0;
  double b_fact = 0.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_maxwell_tetrad = gkyl_wv_gr_maxwell_tetrad_new(light_speed, e_fact, b_fact, spacetime, false);

  TEST_CHECK( gr_maxwell_tetrad->num_equations == 22 );
  TEST_CHECK( gr_maxwell_tetrad->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double Dx = 0.1, Dy = 0.2, Dz = 0.3;
      double Bx = 0.4, By = 0.5, Bz = 0.6;
      double phi = 0.0, psi = 0.0;

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);

      if (!in_excision_region) {
        double q[22];
        q[0] = Dx; q[1] = Dy; q[2] = Dz;
        q[3] = Bx; q[4] = By; q[5] = Bz;
        q[6] = phi; q[7] = psi;

        q[8] = lapse;
        q[9] = shift[0]; q[10] = shift[1]; q[11] = shift[2];

        q[12] = spatial_metric[0][0]; q[13] = spatial_metric[0][1]; q[14] = spatial_metric[0][2];
        q[15] = spatial_metric[1][0]; q[16] = spatial_metric[1][1]; q[17] = spatial_metric[1][2];
        q[18] = spatial_metric[2][0]; q[19] = spatial_metric[2][1]; q[20] = spatial_metric[2][2];

        q[21] = 1.0;

        double Ex = (lapse * Dx) + ((shift[1] * Bz) - (shift[2] * By));
        double Ey = (lapse * Dy) - ((shift[0] * Bz) - (shift[2] * Bx));
        double Ez = (lapse * Dz) + ((shift[0] * By) - (shift[1] * Bx));

        double Hx = (lapse * Bx) - ((shift[1] * Dz) - (shift[2] * Dy));
        double Hy = (lapse * By) + ((shift[0] * Dz) - (shift[2] * Dx));
        double Hz = (lapse * Bz) - ((shift[0] * Dy) - (shift[1] * Dx));

        double fluxes[3][8] = {
          { e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hz, -(light_speed * light_speed) * Hy, b_fact * psi,
            -Ez, Ey, e_fact * Dx, b_fact * (light_speed * light_speed) * Bx },
          { -(light_speed * light_speed) * Hz, e_fact * (light_speed * light_speed) * phi, (light_speed * light_speed) * Hx, Ez, b_fact * psi,
            -Ex, e_fact * Dy, b_fact * (light_speed * light_speed) * By },
          { (light_speed * light_speed) * Hy, -(light_speed * light_speed) * Hx, e_fact * (light_speed * light_speed) * phi, -Ey, Ex, b_fact * psi,
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

        double q_local[22], flux_local_sr[22], flux_local_gr[22], flux[22];
        for (int d = 0; d < 3; d++) {
          gr_maxwell_tetrad->rotate_to_local_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_maxwell_tetrad_flux(light_speed, e_fact, b_fact, q_local, flux_local_sr);
          gkyl_gr_maxwell_tetrad_flux_correction(light_speed, e_fact, b_fact, q_local, flux_local_sr, flux_local_gr);
          gr_maxwell_tetrad->rotate_to_global_func(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], flux_local_gr, flux);

          for (int i = 0; i < 8; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-13) );
          }
        }

        double q_l[6], q_g[6];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_maxwell_tetrad, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 6; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-15) );
          }

          double w1[6], q1[6];
          gr_maxwell_tetrad->cons_to_riem(gr_maxwell_tetrad, q_local, q_local, w1);
          gr_maxwell_tetrad->riem_to_cons(gr_maxwell_tetrad, q_local, w1, q1);

          for (int i = 0; i < 6; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(shift);
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_maxwell_tetrad_basic_minkowski", test_gr_maxwell_tetrad_basic_minkowski },
  { "gr_maxwell_tetrad_basic_schwarzschild", test_gr_maxwell_tetrad_basic_schwarzschild },
  { "gr_maxwell_tetrad_basic_kerr", test_gr_maxwell_tetrad_basic_kerr },
  { NULL, NULL },
};