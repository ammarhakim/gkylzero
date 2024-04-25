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

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));

      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
      double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));

      for (int i = 0; i < 4; i++) {
        spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
        inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
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

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));

      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
      double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));

      for (int i = 0; i < 4; i++) {
        spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
        spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));

        for (int j = 0; j < 4; j++) {
          spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
          spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
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
      double *shift_vector = gkyl_malloc(sizeof(double[3]));

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

      TEST_CHECK( gkyl_compare(lapse_function, 1.0, 1e-10) );
      for (int i = 0; i < 3; i++) {
        TEST_CHECK( gkyl_compare(shift_vector[i], 0.0, 1e-10) );
      }

      double *lapse_function_der = gkyl_malloc(sizeof(double[3]));
      double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));

      for (int i = 0; i < 3; i++) {
        shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
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

void test_gr_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double** spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-10) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-10) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }
}

void test_gr_kerr()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.1 * (1.0 + sqrt(0.19))) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double** spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-10) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-10) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }
}

TEST_LIST = {
  { "gr_minkowski", test_gr_minkowski },
  { "gr_schwarzschild", test_gr_schwarzschild },
  { "gr_kerr", test_gr_kerr },
  { NULL, NULL },
};