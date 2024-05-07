#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_tensor.h>

double*
gkyl_gr_spatial_tensor_reindex_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double* tensor)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);

  double *covariant_tensor = gkyl_malloc(sizeof(double[3]));
  double *new_tensor = gkyl_malloc(sizeof(double[3]));

  for (int i = 0; i < 3; i++) {
    covariant_tensor[i] = 0.0;
    new_tensor[i] = 0.0;
  }

  if (old_indices[0] == true) {
    covariant_tensor = tensor;
  }
  else if (old_indices[0] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        covariant_tensor[i] += spatial_metric[i][j] * tensor[j];
      }
    }
  }

  if (new_indices[0] == true) {
    new_tensor = covariant_tensor;
  }
  else if (new_indices[0] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        new_tensor[i] += inv_spatial_metric[i][j] * covariant_tensor[j];
      }
    }
  }

  return new_tensor;
}

double*
gkyl_gr_spacetime_tensor_reindex_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double* tensor)
{
  double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
  double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
    inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z, &spacetime_metric);
  spacetime->spacetime_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spacetime_metric);

  double *covariant_tensor = gkyl_malloc(sizeof(double[4]));
  double *new_tensor = gkyl_malloc(sizeof(double[4]));

  for (int i = 0; i < 4; i++) {
    covariant_tensor[i] = 0.0;
    new_tensor[i] = 0.0;
  }

  if (old_indices[0] == true) {
    covariant_tensor = tensor;
  }
  else if (old_indices[0] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        covariant_tensor[i] += spacetime_metric[i][j] * tensor[j];
      }
    }
  }

  if (new_indices[0] == true) {
    new_tensor = covariant_tensor;
  }
  else if (new_indices[0] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        new_tensor[i] += inv_spacetime_metric[i][j] * covariant_tensor[j];
      }
    }
  }

  return new_tensor;
}

double**
gkyl_gr_spatial_covariant_der_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, bool* indices, double* tensor, double** tensor_der)
{
  double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
    
    for (int j = 0; j < 3; j++) {
      spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_christoffel_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_christoffel);

  double **covariant_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    covariant_der[i] = gkyl_malloc(sizeof(double[3]));

    for (int j = 0; j < 3; j++) {
      covariant_der[i][j] = tensor_der[i][j];
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if (indices[0] == true) {
          covariant_der[i][j] -= spatial_christoffel[k][i][j] * tensor[k];
        }
        else if (indices[0] == false) {
          covariant_der[i][j] += spatial_christoffel[j][i][k] * tensor[k];
        }
      }
    }
  }

  return covariant_der;
}

double**
gkyl_gr_spacetime_covariant_der_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, bool* indices, double* tensor, double** tensor_der)
{
  double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
    
    for (int j = 0; j < 4; j++) {
      spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
    }
  }

  spacetime->spacetime_christoffel_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_christoffel);

  double **covariant_der = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    covariant_der[i] = gkyl_malloc(sizeof(double[4]));

    for (int j = 0; j < 4; j++) {
      covariant_der[i][j] = tensor_der[i][j];
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        if (indices[0] == true) {
          covariant_der[i][j] -= spacetime_christoffel[k][i][j] * tensor[k];
        }
        else if (indices[0] == false) {
          covariant_der[i][j] += spacetime_christoffel[j][i][k] * tensor[k];
        }
      }
    }
  }

  return covariant_der;
}