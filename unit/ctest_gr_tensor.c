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

        bool *covariant_indices = gkyl_malloc(sizeof(bool[3]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[4]));
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

        bool *covariant_indices = gkyl_malloc(sizeof(bool[3]));
        bool *contravariant_indices = gkyl_malloc(sizeof(bool[4]));
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
      }
    }
  }
}

TEST_LIST = {
  { "gr_tensor_rank1", test_gr_tensor_rank1 },
  { NULL, NULL },
};