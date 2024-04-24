#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_gr_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double **spatial_metric = malloc(sizeof(double*) * 3);
      double **inv_spatial_metric = malloc(sizeof(double*) * 3);
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = malloc(sizeof(double) * 3);
        inv_spatial_metric[i] = malloc(sizeof(double) * 3);
      }

      double **spacetime_metric = malloc(sizeof(double*) * 4);
      double **inv_spacetime_metric = malloc(sizeof(double*) * 4);
      for (int i = 0; i < 4; i++) {
        spacetime_metric[i] = malloc(sizeof(double) * 4);
        inv_spacetime_metric[i] = malloc(sizeof(double) * 4);
      }

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
      spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          if (i == j) {
            TEST_CHECK( gkyl_compare(spatial_metric[i][j], 1.0, 1e-10) );
            TEST_CHECK( gkyl_compare(inv_spatial_metric[i][j], 1.0, 1e-10) );
          }
          else {
            TEST_CHECK( gkyl_compare(spatial_metric[i][j], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(inv_spatial_metric[i][j], 0.0, 1e-10) );
          }
        }
      }

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          if (i == j) {
            if (i == 0) {
              TEST_CHECK( gkyl_compare(spacetime_metric[i][j], -1.0, 1e-10) );
              TEST_CHECK( gkyl_compare(inv_spacetime_metric[i][j], -1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric[i][j], 1.0, 1e-10) );
              TEST_CHECK( gkyl_compare(inv_spacetime_metric[i][j], 1.0, 1e-10) );
            }
          }
          else {
            TEST_CHECK( gkyl_compare(spacetime_metric[i][j], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(inv_spacetime_metric[i][j], 0.0, 1e-10) );
          }
        }
      }

      double spatial_metric_det;
      double spacetime_metric_det;

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
      spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);

      TEST_CHECK( gkyl_compare(spatial_metric_det, 1.0, 1e-10) );
      TEST_CHECK( gkyl_compare(spacetime_metric_det, -1.0, 1e-10) );

      double ***spatial_metric_der = malloc(sizeof(double**) * 3);
      double ***spatial_christoffel = malloc(sizeof(double**) * 3);
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = malloc(sizeof(double*) * 3);
        spatial_christoffel[i] = malloc(sizeof(double*) * 3);

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = malloc(sizeof(double) * 3);
          spatial_christoffel[i][j] = malloc(sizeof(double) * 3);
        }
      }

      double ***spacetime_metric_der = malloc(sizeof(double**) * 4);
      double ***spacetime_christoffel = malloc(sizeof(double**) * 4);
      for (int i = 0; i < 4; i++) {
        spacetime_metric_der[i] = malloc(sizeof(double*) * 4);
        spacetime_christoffel[i] = malloc(sizeof(double*) * 4);

        for (int j = 0; j < 4; j++) {
          spacetime_metric_der[i][j] = malloc(sizeof(double) * 4);
          spacetime_christoffel[i][j] = malloc(sizeof(double) * 4);
        }
      }

      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
      spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

      spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
      spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            TEST_CHECK( gkyl_compare(spatial_metric_der[i][j][k], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(spatial_christoffel[i][j][k], 0.0, 1e-10) );
          }
        }
      }

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
            TEST_CHECK( gkyl_compare(spacetime_metric_der[i][j][k], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(spacetime_christoffel[i][j][k], 0.0, 1e-10) );
          }
        }
      }

      double lapse_function;
      double *shift_vector = malloc(sizeof(double) * 3);

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

      TEST_CHECK( gkyl_compare(lapse_function, 1.0, 1e-10) );
      for (int i = 0; i < 3; i++) {
        TEST_CHECK( gkyl_compare(shift_vector[i], 0.0, 1e-10) );
      }

      double *lapse_function_der = malloc(sizeof(double) * 3);
      double **shift_vector_der = malloc(sizeof(double*) * 3);
      double **extrinsic_curvature = malloc(sizeof(double*) * 3);

      for (int i = 0; i < 3; i++) {
        shift_vector_der[i] = malloc(sizeof(double) * 3);
        extrinsic_curvature[i] = malloc(sizeof(double) * 3);
      }

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &lapse_function_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

      for (int i = 0; i < 3; i++) {
        TEST_CHECK( gkyl_compare(lapse_function_der[i], 0.0, 1e-10) );
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          TEST_CHECK( gkyl_compare(shift_vector_der[i][j], 0.0, 1e-10) );
          TEST_CHECK( gkyl_compare(extrinsic_curvature[i][j], 0.0, 1e-10) );
        }
      }

      bool in_excision_region;
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      TEST_CHECK( (in_excision_region == false) );
    }
  }
}

TEST_LIST = {
  { "gr_minkowski", test_gr_minkowski },
  { NULL, NULL },
};