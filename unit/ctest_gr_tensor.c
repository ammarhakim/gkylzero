#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_tensor.h>

void
test_gr_tensor_rank1()
{
  struct gkyl_gr_spacetime *schwarzschild_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double *spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double *spacetime_tensor = gkyl_malloc(sizeof(double[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor[i] = 0.1 * i;
        }
        for (int i = 0; i < 4; i++) {
          spacetime_tensor[i] = 0.1 * i;
        }

        double *contravariant_spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double *contravariant_spacetime_tensor = gkyl_malloc(sizeof(double[4]));

        bool *covariant_indices = gkyl_malloc(sizeof(bool[1]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[1]));
        covariant_indices[0] = true;
        contravariant_indices[0] = false;

        contravariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spatial_tensor);
        contravariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spacetime_tensor);

        double *covariant_spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double *covariant_spacetime_tensor = gkyl_malloc(sizeof(double[4]));

        covariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spatial_tensor);
        covariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spacetime_tensor);

        for (int i = 0; i < 3; i++) {
          TEST_CHECK( gkyl_compare(spatial_tensor[i], covariant_spatial_tensor[i], 1e-10) );
        }
        for (int i = 0; i < 4; i++) {
          TEST_CHECK( gkyl_compare(spacetime_tensor[i], covariant_spacetime_tensor[i], 1e-10) );
        }

        gkyl_free(spatial_tensor);
        gkyl_free(spacetime_tensor);
        gkyl_free(contravariant_spatial_tensor);
        gkyl_free(contravariant_spacetime_tensor);
        gkyl_free(covariant_indices);
        gkyl_free(contravariant_indices);
        gkyl_free(covariant_spatial_tensor);
        gkyl_free(covariant_spacetime_tensor);
      }
    }
  }

  struct gkyl_gr_spacetime *kerr_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.1 * (1.0 + sqrt(0.19))) {
        double *spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double *spacetime_tensor = gkyl_malloc(sizeof(double[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor[i] = 0.1 * i;
        }
        for (int i = 0; i < 4; i++) {
          spacetime_tensor[i] = 0.1 * i;
        }

        double *contravariant_spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double *contravariant_spacetime_tensor = gkyl_malloc(sizeof(double[4]));

        bool *covariant_indices = gkyl_malloc(sizeof(bool[1]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[1]));
        covariant_indices[0] = true;
        contravariant_indices[0] = false;

        contravariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank1(kerr_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spatial_tensor);
        contravariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank1(kerr_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spacetime_tensor);

        double *covariant_spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double *covariant_spacetime_tensor = gkyl_malloc(sizeof(double[4]));

        covariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank1(kerr_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spatial_tensor);
        covariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank1(kerr_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spacetime_tensor);

        for (int i = 0; i < 3; i++) {
          TEST_CHECK( gkyl_compare(spatial_tensor[i], covariant_spatial_tensor[i], 1e-10) );
        }
        for (int i = 0; i < 4; i++) {
          TEST_CHECK( gkyl_compare(spacetime_tensor[i], covariant_spacetime_tensor[i], 1e-10) );
        }

        gkyl_free(spatial_tensor);
        gkyl_free(spacetime_tensor);
        gkyl_free(contravariant_spatial_tensor);
        gkyl_free(contravariant_spacetime_tensor);
        gkyl_free(covariant_indices);
        gkyl_free(contravariant_indices);
        gkyl_free(covariant_spatial_tensor);
        gkyl_free(covariant_spacetime_tensor);
      }
    }
  }
}

void
test_gr_tensor_rank2()
{
  struct gkyl_gr_spacetime *schwarzschild_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double **spatial_tensor = gkyl_malloc(sizeof(double*[3]));
        double **spacetime_tensor = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor[i][j] = (0.1 * i) + (0.1 * j);
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor[i][j] = (0.1 * i) + (0.1 * j);
          }
        }

        double **contravariant_spatial_tensor = gkyl_malloc(sizeof(double*[3]));
        double **contravariant_spacetime_tensor = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          contravariant_spatial_tensor[i] = gkyl_malloc(sizeof(double[3]));
        }
        for (int i = 0; i < 4; i++) {
          contravariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double[4]));
        }

        bool *covariant_indices = gkyl_malloc(sizeof(bool[2]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[2]));
        covariant_indices[0] = true; covariant_indices[1] = true;
        contravariant_indices[0] = false; contravariant_indices[1] = false;

        contravariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spatial_tensor);
        contravariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spacetime_tensor);

        double **covariant_spatial_tensor = gkyl_malloc(sizeof(double*[3]));
        double **covariant_spacetime_tensor = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          covariant_spatial_tensor[i] = gkyl_malloc(sizeof(double[3]));
        }
        for (int i = 0; i < 3; i++) {
          covariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double[4]));
        }

        covariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spatial_tensor);
        covariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spacetime_tensor);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(spatial_tensor[i][j], covariant_spatial_tensor[i][j], 1e-10) );
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            TEST_CHECK( gkyl_compare(spacetime_tensor[i][j], covariant_spacetime_tensor[i][j], 1e-10) );
          }
        }


        double **spatial_tensor_du = gkyl_malloc(sizeof(double*[3]));
        double **spacetime_tensor_du = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_du[i] = gkyl_malloc(sizeof(double[3]));
        }
        for (int i = 0; i < 4; i++) {
          spacetime_tensor_du[i] = gkyl_malloc(sizeof(double[4]));
        }

        bool *ud_indices = gkyl_malloc(sizeof(bool[2]));
        bool *du_indices = gkyl_malloc(sizeof(bool[2]));
        ud_indices[0] = false; ud_indices[1] = true;
        du_indices[0] = true; du_indices[1] = false;

        spatial_tensor_du = gkyl_gr_spatial_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          ud_indices, du_indices, spatial_tensor);
        spacetime_tensor_du = gkyl_gr_spacetime_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          ud_indices, du_indices, spacetime_tensor);

        double **spatial_tensor_ud = gkyl_malloc(sizeof(double*[3]));
        double **spacetime_tensor_ud = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_ud[i] = gkyl_malloc(sizeof(double[3]));
        }
        for (int i = 0; i < 3; i++) {
          spacetime_tensor_ud[i] = gkyl_malloc(sizeof(double[4]));
        }

        spatial_tensor_ud = gkyl_gr_spatial_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          du_indices, ud_indices, spatial_tensor_du);
        spacetime_tensor_ud = gkyl_gr_spacetime_tensor_reindex_rank2(schwarzschild_spacetime, 0.0, x, y, 0.0,
          du_indices, ud_indices, spacetime_tensor_du);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(spatial_tensor[i][j], spatial_tensor_ud[i][j], 1e-10) );
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            TEST_CHECK( gkyl_compare(spacetime_tensor[i][j], spacetime_tensor_ud[i][j], 1e-10) );
          }
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_tensor[i]);
          gkyl_free(contravariant_spatial_tensor[i]);
          gkyl_free(covariant_spatial_tensor[i]);
          gkyl_free(spatial_tensor_du[i]);
          gkyl_free(spatial_tensor_ud[i]);
        }
        gkyl_free(spatial_tensor);
        gkyl_free(contravariant_spatial_tensor);
        gkyl_free(covariant_spatial_tensor);
        gkyl_free(spatial_tensor_du);
        gkyl_free(spatial_tensor_ud);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_tensor[i]);
          gkyl_free(contravariant_spacetime_tensor[i]);
          gkyl_free(covariant_spacetime_tensor[i]);
          gkyl_free(spacetime_tensor_du[i]);
          gkyl_free(spacetime_tensor_ud[i]);
        }
        gkyl_free(spacetime_tensor);
        gkyl_free(contravariant_spacetime_tensor);
        gkyl_free(covariant_spacetime_tensor);
        gkyl_free(spacetime_tensor_du);
        gkyl_free(spacetime_tensor_ud);

        gkyl_free(covariant_indices);
        gkyl_free(contravariant_indices);
        gkyl_free(ud_indices);
        gkyl_free(du_indices);
      }
    }
  }
}

void
test_gr_covariant_der_rank1()
{
  struct gkyl_gr_spacetime *schwarzschild_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double *spatial_tensor = gkyl_malloc(sizeof(double[3]));
        double **spatial_tensor_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_tensor_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *spacetime_tensor = gkyl_malloc(sizeof(double[4]));
        double **spacetime_tensor_der = gkyl_malloc(sizeof(double*[4]));
        for (int i = 0; i < 4; i++) {
          spacetime_tensor_der[i] = gkyl_malloc(sizeof(double[4]));
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            spatial_tensor_der[i][j] = (pow(10.0, -6.0) * i) + (pow(10.0, -6.0) * j);
          }

          spatial_tensor[i] = 0.1 * i;
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            spacetime_tensor_der[i][j] = (pow(10.0, -6.0) * i) + (pow(10.0, -6.0) * j);
          }

          spacetime_tensor[i] = 0.1 * i;
        }

        bool *covariant_indices = gkyl_malloc(sizeof(bool[1]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[1]));
        covariant_indices[0] = true;
        contravariant_indices[0] = false;

        double **covariant_spatial_der = gkyl_malloc(sizeof(double*[3]));
        double **covariant_spatial_der_comp = gkyl_malloc(sizeof(double*[3]));
        double **contravariant_spatial_der = gkyl_malloc(sizeof(double*[3]));
        double **contravariant_spatial_der_comp = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          covariant_spatial_der[i] = gkyl_malloc(sizeof(double[3]));
          covariant_spatial_der_comp[i] = gkyl_malloc(sizeof(double[3]));
          contravariant_spatial_der[i] = gkyl_malloc(sizeof(double[3]));
          contravariant_spatial_der_comp[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            covariant_spatial_der_comp[i][j] = 0.0;
            contravariant_spatial_der_comp[i][j] = 0.0;
          }
        }

        double **covariant_spacetime_der = gkyl_malloc(sizeof(double*[4]));
        double **covariant_spacetime_der_comp = gkyl_malloc(sizeof(double*[4]));
        double **contravariant_spacetime_der = gkyl_malloc(sizeof(double*[4]));
        double **contravariant_spacetime_der_comp = gkyl_malloc(sizeof(double*[4]));
        for (int i = 0; i < 4; i++) {
          covariant_spacetime_der[i] = gkyl_malloc(sizeof(double[4]));
          covariant_spacetime_der_comp[i] = gkyl_malloc(sizeof(double[4]));
          contravariant_spacetime_der[i] = gkyl_malloc(sizeof(double[4]));
          contravariant_spacetime_der_comp[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            covariant_spacetime_der_comp[i][j] = 0.0;
            contravariant_spacetime_der_comp[i][j] = 0.0;
          }
        }

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

        schwarzschild_spacetime->spatial_metric_tensor_func(schwarzschild_spacetime, 0.0, x, y, 0.0, &spatial_metric);
        schwarzschild_spacetime->spatial_inv_metric_tensor_func(schwarzschild_spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        
        schwarzschild_spacetime->spacetime_metric_tensor_func(schwarzschild_spacetime, 0.0, x, y, 0.0, &spacetime_metric);
        schwarzschild_spacetime->spacetime_inv_metric_tensor_func(schwarzschild_spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        double *spatial_tensor_contr = gkyl_malloc(sizeof(double[3]));
        double **spatial_tensor_der_contr = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_tensor_der_contr[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_der_contr[i][j] = 0.0;
          }
          spatial_tensor_contr[i] = 0.0;
        }

        double *spacetime_tensor_contr = gkyl_malloc(sizeof(double[4]));
        double **spacetime_tensor_der_contr = gkyl_malloc(sizeof(double*[4]));
        for (int i = 0; i < 4; i++) {
          spacetime_tensor_der_contr[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_der_contr[i][j] = 0.0;
          }
          spacetime_tensor_contr[i] = 0.0;
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            spatial_tensor_contr[i] += inv_spatial_metric[i][j] * spatial_tensor[j];

            for (int k = 0; k < 3; k++) {
              spatial_tensor_der_contr[i][j] += inv_spatial_metric[k][j] * spatial_tensor_der[i][k];
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            spacetime_tensor_contr[i] += inv_spacetime_metric[i][j] * spacetime_tensor[j];

            for (int k = 0; k < 4; k++) {
              spacetime_tensor_der_contr[i][j] += inv_spacetime_metric[k][j] * spacetime_tensor_der[i][k];
            }
          }
        }

        covariant_spatial_der = gkyl_gr_spatial_covariant_der_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), covariant_indices, spatial_tensor, spatial_tensor_der);
        contravariant_spatial_der = gkyl_gr_spatial_covariant_der_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), contravariant_indices, spatial_tensor_contr, spatial_tensor_der_contr);

        covariant_spacetime_der = gkyl_gr_spacetime_covariant_der_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), covariant_indices, spacetime_tensor, spacetime_tensor_der);
        contravariant_spacetime_der = gkyl_gr_spacetime_covariant_der_rank1(schwarzschild_spacetime, 0.0, x, y, 0.0,
          pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), contravariant_indices, spacetime_tensor_contr, spacetime_tensor_der_contr);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              covariant_spatial_der_comp[i][j] += spatial_metric[k][j] * contravariant_spatial_der[i][k];
              contravariant_spatial_der_comp[i][j] += inv_spatial_metric[k][j] * covariant_spatial_der[i][k];
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              covariant_spacetime_der_comp[i][j] += spacetime_metric[k][j] * contravariant_spacetime_der[i][k];
              contravariant_spacetime_der_comp[i][j] += inv_spacetime_metric[k][j] * contravariant_spacetime_der[i][k];
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            //TEST_CHECK( gkyl_compare(covariant_spatial_der[i][j], covariant_spatial_der_comp[i][j], 1e-2) );
            //TEST_CHECK( gkyl_compare(contravariant_spatial_der[i][j], contravariant_spatial_der_comp[i][j], 1e-2) );
          }
        }
      }
    }
  }
}

TEST_LIST = {
  { "gr_tensor_rank1", test_gr_tensor_rank1 },
  { "gr_tensor_rank2", test_gr_tensor_rank2 },
  { "gr_covariant_der_rank1", test_gr_covariant_der_rank1 },
  { NULL, NULL },
};