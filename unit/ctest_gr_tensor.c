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
        for (int i = 0; i < 4; i++) {
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
        for (int i = 0; i < 4; i++) {
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

  struct gkyl_gr_spacetime *kerr_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);

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

        contravariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spatial_tensor);
        contravariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spacetime_tensor);

        double **covariant_spatial_tensor = gkyl_malloc(sizeof(double*[3]));
        double **covariant_spacetime_tensor = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          covariant_spatial_tensor[i] = gkyl_malloc(sizeof(double[3]));
        }
        for (int i = 0; i < 4; i++) {
          covariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double[4]));
        }

        covariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spatial_tensor);
        covariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
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

        spatial_tensor_du = gkyl_gr_spatial_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
          ud_indices, du_indices, spatial_tensor);
        spacetime_tensor_du = gkyl_gr_spacetime_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
          ud_indices, du_indices, spacetime_tensor);

        double **spatial_tensor_ud = gkyl_malloc(sizeof(double*[3]));
        double **spacetime_tensor_ud = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_ud[i] = gkyl_malloc(sizeof(double[3]));
        }
        for (int i = 0; i < 4; i++) {
          spacetime_tensor_ud[i] = gkyl_malloc(sizeof(double[4]));
        }

        spatial_tensor_ud = gkyl_gr_spatial_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
          du_indices, ud_indices, spatial_tensor_du);
        spacetime_tensor_ud = gkyl_gr_spacetime_tensor_reindex_rank2(kerr_spacetime, 0.0, x, y, 0.0,
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
test_gr_tensor_rank3()
{
  struct gkyl_gr_spacetime *schwarzschild_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double ***spatial_tensor = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor[i][j] = gkyl_malloc(sizeof(double[3]));

            for (int k = 0; k < 3; k++) {
              spatial_tensor[i][j][k] = (0.1 * i) + (0.1 * j) + (0.1 * k);
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor[i][j] = gkyl_malloc(sizeof(double[4]));

            for (int k = 0; k < 4; k++) {
              spacetime_tensor[i][j][k] = (0.1 * i) + (0.1 * j) + (0.1 * k);
            }
          }
        }

        double ***contravariant_spatial_tensor = gkyl_malloc(sizeof(double**[3]));
        double ***contravariant_spacetime_tensor = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          contravariant_spatial_tensor[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            contravariant_spatial_tensor[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          contravariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            contravariant_spacetime_tensor[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *covariant_indices = gkyl_malloc(sizeof(bool[3]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[3]));
        covariant_indices[0] = true; covariant_indices[1] = true; covariant_indices[2] = true;
        contravariant_indices[0] = false; contravariant_indices[1] = false; contravariant_indices[2] = false;

        contravariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spatial_tensor);
        contravariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spacetime_tensor);

        double ***covariant_spatial_tensor = gkyl_malloc(sizeof(double**[3]));
        double ***covariant_spacetime_tensor = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          covariant_spatial_tensor[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            covariant_spatial_tensor[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 3; i++) {
          covariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            covariant_spacetime_tensor[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        covariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spatial_tensor);
        covariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spacetime_tensor);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], covariant_spatial_tensor[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], covariant_spacetime_tensor[i][j][k], 1e-10) );
            }
          }
        }

        double ***spatial_tensor_ddu = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_ddu = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_ddu[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_ddu[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_ddu[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_ddu[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *uud_indices = gkyl_malloc(sizeof(bool[3]));
        bool *ddu_indices = gkyl_malloc(sizeof(bool[3]));
        uud_indices[0] = false; uud_indices[1] = false; uud_indices[2] = true;
        ddu_indices[0] = true; ddu_indices[1] = true; ddu_indices[2] = false;

        spatial_tensor_ddu = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          uud_indices, ddu_indices, spatial_tensor);
        spacetime_tensor_ddu = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          uud_indices, ddu_indices, spacetime_tensor);

        double ***spatial_tensor_uud = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_uud = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_uud[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_uud[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_uud[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_uud[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spatial_tensor_uud = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          ddu_indices, uud_indices, spatial_tensor_ddu);
        spacetime_tensor_uud = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          ddu_indices, uud_indices, spacetime_tensor_ddu);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], spatial_tensor_uud[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], spacetime_tensor_uud[i][j][k], 1e-10) );
            }
          }
        }

        double ***spatial_tensor_udd = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_udd = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_udd[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_udd[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_udd[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_udd[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *duu_indices = gkyl_malloc(sizeof(bool[3]));
        bool *udd_indices = gkyl_malloc(sizeof(bool[3]));
        duu_indices[0] = true; duu_indices[1] = false; duu_indices[2] = false;
        udd_indices[0] = false; udd_indices[1] = true; udd_indices[2] = true;

        spatial_tensor_udd = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          duu_indices, udd_indices, spatial_tensor);
        spacetime_tensor_udd = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          duu_indices, udd_indices, spacetime_tensor);

        double ***spatial_tensor_duu = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_duu = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_duu[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_duu[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_duu[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_duu[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spatial_tensor_duu = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          udd_indices, duu_indices, spatial_tensor_udd);
        spacetime_tensor_duu = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          udd_indices, duu_indices, spacetime_tensor_udd);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], spatial_tensor_duu[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], spacetime_tensor_duu[i][j][k], 1e-10) );
            }
          }
        }

        double ***spatial_tensor_dud = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_dud = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_dud[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_dud[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_dud[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_dud[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *udu_indices = gkyl_malloc(sizeof(bool[3]));
        bool *dud_indices = gkyl_malloc(sizeof(bool[3]));
        udu_indices[0] = false; udu_indices[1] = true; udu_indices[2] = false;
        dud_indices[0] = true; dud_indices[1] = false; dud_indices[2] = true;

        spatial_tensor_dud = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          udu_indices, dud_indices, spatial_tensor);
        spacetime_tensor_dud = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          udu_indices, dud_indices, spacetime_tensor);

        double ***spatial_tensor_udu = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_udu = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_udu[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_udu[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_udu[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_udu[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spatial_tensor_udu = gkyl_gr_spatial_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          dud_indices, udu_indices, spatial_tensor_dud);
        spacetime_tensor_udu = gkyl_gr_spacetime_tensor_reindex_rank3(schwarzschild_spacetime, 0.0, x, y, 0.0,
          dud_indices, udu_indices, spacetime_tensor_dud);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], spatial_tensor_udu[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], spacetime_tensor_udu[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_tensor[i][j]);
            gkyl_free(contravariant_spatial_tensor[i][j]);
            gkyl_free(covariant_spatial_tensor[i][j]);
            gkyl_free(spatial_tensor_ddu[i][j]);
            gkyl_free(spatial_tensor_uud[i][j]);
            gkyl_free(spatial_tensor_udd[i][j]);
            gkyl_free(spatial_tensor_duu[i][j]);
            gkyl_free(spatial_tensor_dud[i][j]);
            gkyl_free(spatial_tensor_udu[i][j]);
          }
          gkyl_free(spatial_tensor[i]);
          gkyl_free(contravariant_spatial_tensor[i]);
          gkyl_free(covariant_spatial_tensor[i]);
          gkyl_free(spatial_tensor_ddu[i]);
          gkyl_free(spatial_tensor_uud[i]);
          gkyl_free(spatial_tensor_udd[i]);
          gkyl_free(spatial_tensor_duu[i]);
          gkyl_free(spatial_tensor_dud[i]);
          gkyl_free(spatial_tensor_udu[i]);
        }
        gkyl_free(spatial_tensor);
        gkyl_free(contravariant_spatial_tensor);
        gkyl_free(covariant_spatial_tensor);
        gkyl_free(spatial_tensor_ddu);
        gkyl_free(spatial_tensor_uud);
        gkyl_free(spatial_tensor_udd);
        gkyl_free(spatial_tensor_duu);
        gkyl_free(spatial_tensor_dud);
        gkyl_free(spatial_tensor_udu);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_tensor[i][j]);
            gkyl_free(contravariant_spacetime_tensor[i][j]);
            gkyl_free(covariant_spacetime_tensor[i][j]);
            gkyl_free(spacetime_tensor_ddu[i][j]);
            gkyl_free(spacetime_tensor_uud[i][j]);
            gkyl_free(spacetime_tensor_udd[i][j]);
            gkyl_free(spacetime_tensor_duu[i][j]);
            gkyl_free(spacetime_tensor_dud[i][j]);
            gkyl_free(spacetime_tensor_udu[i][j]);
          }
          gkyl_free(spacetime_tensor[i]);
          gkyl_free(contravariant_spacetime_tensor[i]);
          gkyl_free(covariant_spacetime_tensor[i]);
          gkyl_free(spacetime_tensor_ddu[i]);
          gkyl_free(spacetime_tensor_uud[i]);
          gkyl_free(spacetime_tensor_udd[i]);
          gkyl_free(spacetime_tensor_duu[i]);
          gkyl_free(spacetime_tensor_dud[i]);
          gkyl_free(spacetime_tensor_udu[i]);
        }
        gkyl_free(spacetime_tensor);
        gkyl_free(contravariant_spacetime_tensor);
        gkyl_free(covariant_spacetime_tensor);
        gkyl_free(spacetime_tensor_ddu);
        gkyl_free(spacetime_tensor_uud);
        gkyl_free(spacetime_tensor_udd);
        gkyl_free(spacetime_tensor_duu);
        gkyl_free(spacetime_tensor_dud);
        gkyl_free(spacetime_tensor_udu);

        gkyl_free(covariant_indices);
        gkyl_free(contravariant_indices);
        gkyl_free(uud_indices);
        gkyl_free(ddu_indices);
        gkyl_free(duu_indices);
        gkyl_free(udd_indices);
        gkyl_free(udu_indices);
        gkyl_free(dud_indices);
      }
    }
  }

  struct gkyl_gr_spacetime *kerr_spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double ***spatial_tensor = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor[i][j] = gkyl_malloc(sizeof(double[3]));

            for (int k = 0; k < 3; k++) {
              spatial_tensor[i][j][k] = (0.1 * i) + (0.1 * j) + (0.1 * k);
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor[i][j] = gkyl_malloc(sizeof(double[4]));

            for (int k = 0; k < 4; k++) {
              spacetime_tensor[i][j][k] = (0.1 * i) + (0.1 * j) + (0.1 * k);
            }
          }
        }

        double ***contravariant_spatial_tensor = gkyl_malloc(sizeof(double**[3]));
        double ***contravariant_spacetime_tensor = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          contravariant_spatial_tensor[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            contravariant_spatial_tensor[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          contravariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            contravariant_spacetime_tensor[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *covariant_indices = gkyl_malloc(sizeof(bool[3]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[3]));
        covariant_indices[0] = true; covariant_indices[1] = true; covariant_indices[2] = true;
        contravariant_indices[0] = false; contravariant_indices[1] = false; contravariant_indices[2] = false;

        contravariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spatial_tensor);
        contravariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          covariant_indices, contravariant_indices, spacetime_tensor);

        double ***covariant_spatial_tensor = gkyl_malloc(sizeof(double**[3]));
        double ***covariant_spacetime_tensor = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          covariant_spatial_tensor[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            covariant_spatial_tensor[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 3; i++) {
          covariant_spacetime_tensor[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            covariant_spacetime_tensor[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        covariant_spatial_tensor = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spatial_tensor);
        covariant_spacetime_tensor = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          contravariant_indices, covariant_indices, contravariant_spacetime_tensor);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], covariant_spatial_tensor[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], covariant_spacetime_tensor[i][j][k], 1e-10) );
            }
          }
        }

        double ***spatial_tensor_ddu = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_ddu = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_ddu[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_ddu[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_ddu[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_ddu[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *uud_indices = gkyl_malloc(sizeof(bool[3]));
        bool *ddu_indices = gkyl_malloc(sizeof(bool[3]));
        uud_indices[0] = false; uud_indices[1] = false; uud_indices[2] = true;
        ddu_indices[0] = true; ddu_indices[1] = true; ddu_indices[2] = false;

        spatial_tensor_ddu = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          uud_indices, ddu_indices, spatial_tensor);
        spacetime_tensor_ddu = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          uud_indices, ddu_indices, spacetime_tensor);

        double ***spatial_tensor_uud = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_uud = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_uud[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_uud[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_uud[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_uud[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spatial_tensor_uud = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          ddu_indices, uud_indices, spatial_tensor_ddu);
        spacetime_tensor_uud = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          ddu_indices, uud_indices, spacetime_tensor_ddu);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], spatial_tensor_uud[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], spacetime_tensor_uud[i][j][k], 1e-10) );
            }
          }
        }

        double ***spatial_tensor_udd = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_udd = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_udd[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_udd[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_udd[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_udd[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *duu_indices = gkyl_malloc(sizeof(bool[3]));
        bool *udd_indices = gkyl_malloc(sizeof(bool[3]));
        duu_indices[0] = true; duu_indices[1] = false; duu_indices[2] = false;
        udd_indices[0] = false; udd_indices[1] = true; udd_indices[2] = true;

        spatial_tensor_udd = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          duu_indices, udd_indices, spatial_tensor);
        spacetime_tensor_udd = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          duu_indices, udd_indices, spacetime_tensor);

        double ***spatial_tensor_duu = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_duu = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_duu[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_duu[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_duu[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_duu[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spatial_tensor_duu = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          udd_indices, duu_indices, spatial_tensor_udd);
        spacetime_tensor_duu = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          udd_indices, duu_indices, spacetime_tensor_udd);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], spatial_tensor_duu[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], spacetime_tensor_duu[i][j][k], 1e-10) );
            }
          }
        }

        double ***spatial_tensor_dud = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_dud = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_dud[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_dud[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_dud[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_dud[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        bool *udu_indices = gkyl_malloc(sizeof(bool[3]));
        bool *dud_indices = gkyl_malloc(sizeof(bool[3]));
        udu_indices[0] = false; udu_indices[1] = true; udu_indices[2] = false;
        dud_indices[0] = true; dud_indices[1] = false; dud_indices[2] = true;

        spatial_tensor_dud = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          udu_indices, dud_indices, spatial_tensor);
        spacetime_tensor_dud = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          udu_indices, dud_indices, spacetime_tensor);

        double ***spatial_tensor_udu = gkyl_malloc(sizeof(double**[3]));
        double ***spacetime_tensor_udu = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 3; i++) {
          spatial_tensor_udu[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_tensor_udu[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        for (int i = 0; i < 4; i++) {
          spacetime_tensor_udu[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_tensor_udu[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spatial_tensor_udu = gkyl_gr_spatial_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          dud_indices, udu_indices, spatial_tensor_dud);
        spacetime_tensor_udu = gkyl_gr_spacetime_tensor_reindex_rank3(kerr_spacetime, 0.0, x, y, 0.0,
          dud_indices, udu_indices, spacetime_tensor_dud);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              TEST_CHECK( gkyl_compare(spatial_tensor[i][j][k], spatial_tensor_udu[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              TEST_CHECK( gkyl_compare(spacetime_tensor[i][j][k], spacetime_tensor_udu[i][j][k], 1e-10) );
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_tensor[i][j]);
            gkyl_free(contravariant_spatial_tensor[i][j]);
            gkyl_free(covariant_spatial_tensor[i][j]);
            gkyl_free(spatial_tensor_ddu[i][j]);
            gkyl_free(spatial_tensor_uud[i][j]);
            gkyl_free(spatial_tensor_udd[i][j]);
            gkyl_free(spatial_tensor_duu[i][j]);
            gkyl_free(spatial_tensor_dud[i][j]);
            gkyl_free(spatial_tensor_udu[i][j]);
          }
          gkyl_free(spatial_tensor[i]);
          gkyl_free(contravariant_spatial_tensor[i]);
          gkyl_free(covariant_spatial_tensor[i]);
          gkyl_free(spatial_tensor_ddu[i]);
          gkyl_free(spatial_tensor_uud[i]);
          gkyl_free(spatial_tensor_udd[i]);
          gkyl_free(spatial_tensor_duu[i]);
          gkyl_free(spatial_tensor_dud[i]);
          gkyl_free(spatial_tensor_udu[i]);
        }
        gkyl_free(spatial_tensor);
        gkyl_free(contravariant_spatial_tensor);
        gkyl_free(covariant_spatial_tensor);
        gkyl_free(spatial_tensor_ddu);
        gkyl_free(spatial_tensor_uud);
        gkyl_free(spatial_tensor_udd);
        gkyl_free(spatial_tensor_duu);
        gkyl_free(spatial_tensor_dud);
        gkyl_free(spatial_tensor_udu);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_tensor[i][j]);
            gkyl_free(contravariant_spacetime_tensor[i][j]);
            gkyl_free(covariant_spacetime_tensor[i][j]);
            gkyl_free(spacetime_tensor_ddu[i][j]);
            gkyl_free(spacetime_tensor_uud[i][j]);
            gkyl_free(spacetime_tensor_udd[i][j]);
            gkyl_free(spacetime_tensor_duu[i][j]);
            gkyl_free(spacetime_tensor_dud[i][j]);
            gkyl_free(spacetime_tensor_udu[i][j]);
          }
          gkyl_free(spacetime_tensor[i]);
          gkyl_free(contravariant_spacetime_tensor[i]);
          gkyl_free(covariant_spacetime_tensor[i]);
          gkyl_free(spacetime_tensor_ddu[i]);
          gkyl_free(spacetime_tensor_uud[i]);
          gkyl_free(spacetime_tensor_udd[i]);
          gkyl_free(spacetime_tensor_duu[i]);
          gkyl_free(spacetime_tensor_dud[i]);
          gkyl_free(spacetime_tensor_udu[i]);
        }
        gkyl_free(spacetime_tensor);
        gkyl_free(contravariant_spacetime_tensor);
        gkyl_free(covariant_spacetime_tensor);
        gkyl_free(spacetime_tensor_ddu);
        gkyl_free(spacetime_tensor_uud);
        gkyl_free(spacetime_tensor_udd);
        gkyl_free(spacetime_tensor_duu);
        gkyl_free(spacetime_tensor_dud);
        gkyl_free(spacetime_tensor_udu);

        gkyl_free(covariant_indices);
        gkyl_free(contravariant_indices);
        gkyl_free(uud_indices);
        gkyl_free(ddu_indices);
        gkyl_free(duu_indices);
        gkyl_free(udd_indices);
        gkyl_free(udu_indices);
        gkyl_free(dud_indices);
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
  { "gr_tensor_rank3", test_gr_tensor_rank3 },
  { "gr_covariant_der_rank1", test_gr_covariant_der_rank1 },
  { NULL, NULL },
};