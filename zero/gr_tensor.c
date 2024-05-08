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
gkyl_gr_spatial_tensor_reindex_rank2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double** tensor)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);

  double **covariant_tensor = gkyl_malloc(sizeof(double*[3]));
  double **new_tensor = gkyl_malloc(sizeof(double*[3]));

  for (int i = 0; i < 3; i++) {
    covariant_tensor[i] = gkyl_malloc(sizeof(double[3]));
    new_tensor[i] = gkyl_malloc(sizeof(double[3]));

    for (int j = 0; j < 3; j++) {
      covariant_tensor[i][j] = 0.0;
      new_tensor[i][j] = 0.0;
    }
  }

  if (old_indices[0] == true && old_indices[1] == true) {
    covariant_tensor = tensor;
  }
  else if (old_indices[0] == true && old_indices[1] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          covariant_tensor[i][j] += spatial_metric[k][j] * tensor[i][k];
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          covariant_tensor[i][j] += spatial_metric[i][k] * tensor[k][j];
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            covariant_tensor[i][j] += spatial_metric[i][k] * spatial_metric[l][j] * tensor[k][l];
          }
        }
      }
    }
  }

  if (new_indices[0] == true && new_indices[1] == true) {
    new_tensor = covariant_tensor;
  }
  else if (new_indices[0] == true && new_indices[1] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          new_tensor[i][j] += inv_spatial_metric[k][j] * covariant_tensor[i][k];
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          new_tensor[i][j] += inv_spatial_metric[i][k] * covariant_tensor[k][j];
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            new_tensor[i][j] += inv_spatial_metric[i][k] * inv_spatial_metric[l][j] * covariant_tensor[k][l];
          }
        }
      }
    }
  }

  return new_tensor;
}

double**
gkyl_gr_spacetime_tensor_reindex_rank2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double** tensor)
{
  double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
  double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
    inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z, &spacetime_metric);
  spacetime->spacetime_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spacetime_metric);

  double **covariant_tensor = gkyl_malloc(sizeof(double*[4]));
  double **new_tensor = gkyl_malloc(sizeof(double*[4]));

  for (int i = 0; i < 4; i++) {
    covariant_tensor[i] = gkyl_malloc(sizeof(double[4]));
    new_tensor[i] = gkyl_malloc(sizeof(double[4]));

    for (int j = 0; j < 4; j++) {
      covariant_tensor[i][j] = 0.0;
      new_tensor[i][j] = 0.0;
    }
  }

  if (old_indices[0] == true && old_indices[1] == true) {
    covariant_tensor = tensor;
  }
  else if (old_indices[0] == true && old_indices[1] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          covariant_tensor[i][j] += spacetime_metric[k][j] * tensor[i][k];
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          covariant_tensor[i][j] += spacetime_metric[i][k] * tensor[k][j];
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            covariant_tensor[i][j] += spacetime_metric[i][k] * spacetime_metric[l][j] * tensor[k][l];
          }
        }
      }
    }
  }

  if (new_indices[0] == true && new_indices[1] == true) {
    new_tensor = covariant_tensor;
  }
  else if (new_indices[0] == true && new_indices[1] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          new_tensor[i][j] += inv_spacetime_metric[k][j] * covariant_tensor[i][k];
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          new_tensor[i][j] += inv_spacetime_metric[i][k] * covariant_tensor[k][j];
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            new_tensor[i][j] += inv_spacetime_metric[i][k] * inv_spacetime_metric[l][j] * covariant_tensor[k][l];
          }
        }
      }
    }
  }

  return new_tensor;
}

double***
gkyl_gr_spatial_tensor_reindex_rank3(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double*** tensor)
{
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_tensor_func(spacetime, t, x, y, z, &spatial_metric);
  spacetime->spatial_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spatial_metric);

  double ***covariant_tensor = gkyl_malloc(sizeof(double**[3]));
  double ***new_tensor = gkyl_malloc(sizeof(double**[3]));

  for (int i = 0; i < 3; i++) {
    covariant_tensor[i] = gkyl_malloc(sizeof(double*[3]));
    new_tensor[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      covariant_tensor[i][j] = gkyl_malloc(sizeof(double[3]));
      new_tensor[i][j] = gkyl_malloc(sizeof(double[3]));

      for (int k = 0; k < 3; k++) {
        covariant_tensor[i][j][k] = 0.0;
        new_tensor[i][j][k] = 0.0;
      }
    }
  }

  if (old_indices[0] == true && old_indices[1] == true && old_indices[2] == true) {
    covariant_tensor = tensor;
  }
  else if (old_indices[0] == true && old_indices[1] == true && old_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            covariant_tensor[i][j][k] += spatial_metric[l][k] * tensor[i][j][l];
          }
        }
      }
    }
  }
  else if (old_indices[0] == true && old_indices[1] == false && old_indices[2] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            covariant_tensor[i][j][k] += spatial_metric[l][j] * tensor[i][l][k];
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == true && old_indices[2] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            covariant_tensor[i][j][k] += spatial_metric[i][l] * tensor[l][j][k];
          }
        }
      }
    }
  }
  else if (old_indices[0] == true && old_indices[1] == false && old_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              covariant_tensor[i][j][k] += spatial_metric[j][l] * spatial_metric[m][k] * tensor[i][l][m];
            }
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == true && old_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              covariant_tensor[i][j][k] += spatial_metric[i][l] * spatial_metric[m][k] * tensor[l][j][m];
            }
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == false && old_indices[2] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              covariant_tensor[i][j][k] += spatial_metric[i][l] * spatial_metric[m][j] * tensor[l][m][k];
            }
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == false && old_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              for (int n = 0; n < 3; n++) {
                covariant_tensor[i][j][k] += spatial_metric[i][l] * spatial_metric[j][m] * spatial_metric[n][k] * tensor[l][m][n];
              }
            }
          }
        }
      }
    }
  }

  if (new_indices[0] == true && new_indices[1] == true && new_indices[2] == true) {
    new_tensor = covariant_tensor;
  }
  else if (new_indices[0] == true && new_indices[1] == true && new_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            new_tensor[i][j][k] += inv_spatial_metric[l][k] * covariant_tensor[i][j][l];
          }
        }
      }
    }
  }
  else if (new_indices[0] == true && new_indices[1] == false && new_indices[2] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            new_tensor[i][j][k] += inv_spatial_metric[l][j] * covariant_tensor[i][l][k];
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == true && new_indices[2] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            new_tensor[i][j][k] += inv_spatial_metric[i][l] * covariant_tensor[l][j][k];
          }
        }
      }
    }
  }
  else if (new_indices[0] == true && new_indices[1] == false && new_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              new_tensor[i][j][k] += inv_spatial_metric[j][l] * inv_spatial_metric[m][k] * covariant_tensor[i][l][m];
            }
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == true && new_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              new_tensor[i][j][k] += inv_spatial_metric[i][l] * inv_spatial_metric[m][k] * covariant_tensor[l][j][m];
            }
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == false && new_indices[2] == true) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              new_tensor[i][j][k] += inv_spatial_metric[i][l] * inv_spatial_metric[m][j] * covariant_tensor[l][m][k];
            }
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == false && new_indices[2] == false) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              for (int n = 0; n < 3; n++) {
                new_tensor[i][j][k] += inv_spatial_metric[i][l] * inv_spatial_metric[j][m] * inv_spatial_metric[n][k] * covariant_tensor[l][m][n];
              }
            }
          }
        }
      }
    }
  }

  return new_tensor;
}

double***
gkyl_gr_spacetime_tensor_reindex_rank3(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double*** tensor)
{
  double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
  double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));

  for (int i = 0; i < 4; i++) {
    spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
    inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime->spacetime_metric_tensor_func(spacetime, t, x, y, z, &spacetime_metric);
  spacetime->spacetime_inv_metric_tensor_func(spacetime, t, x, y, z, &inv_spacetime_metric);

  double ***covariant_tensor = gkyl_malloc(sizeof(double**[4]));
  double ***new_tensor = gkyl_malloc(sizeof(double**[4]));

  for (int i = 0; i < 4; i++) {
    covariant_tensor[i] = gkyl_malloc(sizeof(double*[4]));
    new_tensor[i] = gkyl_malloc(sizeof(double*[4]));

    for (int j = 0; j < 4; j++) {
      covariant_tensor[i][j] = gkyl_malloc(sizeof(double[4]));
      new_tensor[i][j] = gkyl_malloc(sizeof(double[4]));

      for (int k = 0; k < 4; k++) {
        covariant_tensor[i][j][k] = 0.0;
        new_tensor[i][j][k] = 0.0;
      }
    }
  }

  if (old_indices[0] == true && old_indices[1] == true && old_indices[2] == true) {
    covariant_tensor = tensor;
  }
  else if (old_indices[0] == true && old_indices[1] == true && old_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            covariant_tensor[i][j][k] += spacetime_metric[l][k] * tensor[i][j][l];
          }
        }
      }
    }
  }
  else if (old_indices[0] == true && old_indices[1] == false && old_indices[2] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            covariant_tensor[i][j][k] += spacetime_metric[l][j] * tensor[i][l][k];
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == true && old_indices[2] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            covariant_tensor[i][j][k] += spacetime_metric[i][l] * tensor[l][j][k];
          }
        }
      }
    }
  }
  else if (old_indices[0] == true && old_indices[1] == false && old_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              covariant_tensor[i][j][k] += spacetime_metric[j][l] * spacetime_metric[m][k] * tensor[i][l][m];
            }
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == true && old_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              covariant_tensor[i][j][k] += spacetime_metric[i][l] * spacetime_metric[m][k] * tensor[l][j][m];
            }
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == false && old_indices[2] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              covariant_tensor[i][j][k] += spacetime_metric[i][l] * spacetime_metric[m][j] * tensor[l][m][k];
            }
          }
        }
      }
    }
  }
  else if (old_indices[0] == false && old_indices[1] == false && old_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              for (int n = 0; n < 4; n++) {
                covariant_tensor[i][j][k] += spacetime_metric[i][l] * spacetime_metric[j][m] * spacetime_metric[n][k] * tensor[l][m][n];
              }
            }
          }
        }
      }
    }
  }

  if (new_indices[0] == true && new_indices[1] == true && new_indices[2] == true) {
    new_tensor = covariant_tensor;
  }
  else if (new_indices[0] == true && new_indices[1] == true && new_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            new_tensor[i][j][k] += inv_spacetime_metric[l][k] * covariant_tensor[i][j][l];
          }
        }
      }
    }
  }
  else if (new_indices[0] == true && new_indices[1] == false && new_indices[2] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            new_tensor[i][j][k] += inv_spacetime_metric[l][j] * covariant_tensor[i][l][k];
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == true && new_indices[2] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            new_tensor[i][j][k] += inv_spacetime_metric[i][l] * covariant_tensor[l][j][k];
          }
        }
      }
    }
  }
  else if (new_indices[0] == true && new_indices[1] == false && new_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              new_tensor[i][j][k] += inv_spacetime_metric[j][l] * inv_spacetime_metric[m][k] * covariant_tensor[i][l][m];
            }
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == true && new_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              new_tensor[i][j][k] += inv_spacetime_metric[i][l] * inv_spacetime_metric[m][k] * covariant_tensor[l][j][m];
            }
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == false && new_indices[2] == true) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              new_tensor[i][j][k] += inv_spacetime_metric[i][l] * inv_spacetime_metric[m][j] * covariant_tensor[l][m][k];
            }
          }
        }
      }
    }
  }
  else if (new_indices[0] == false && new_indices[1] == false && new_indices[2] == false) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            for (int m = 0; m < 4; m++) {
              for (int n = 0; n < 4; n++) {
                new_tensor[i][j][k] += inv_spacetime_metric[i][l] * inv_spacetime_metric[j][m] * inv_spacetime_metric[n][k] * covariant_tensor[l][m][n];
              }
            }
          }
        }
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

double***
gkyl_gr_spatial_covariant_der_rank2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, bool* indices, double** tensor, double*** tensor_der)
{
  double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
    
    for (int j = 0; j < 3; j++) {
      spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
    }
  }

  spacetime->spatial_christoffel_func(spacetime, t, x, y, z, dx, dy, dz, &spatial_christoffel);

  double ***covariant_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    covariant_der[i] = gkyl_malloc(sizeof(double*[3]));

    for (int j = 0; j < 3; j++) {
      covariant_der[i][j] = gkyl_malloc(sizeof(double[3]));

      for (int k = 0; k < 3; k++) {
        covariant_der[i][j][k] = tensor_der[i][j][k];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          if (indices[0] == true && indices[1] == true) {
            covariant_der[i][j][k] -= spatial_christoffel[l][i][j] * tensor[l][k];
            covariant_der[i][j][k] -= spatial_christoffel[l][i][k] * tensor[j][l];
          }
          else if (indices[0] == true && indices[1] == false) {
            covariant_der[i][j][k] += spatial_christoffel[k][i][l] * tensor[j][l];
            covariant_der[i][j][k] -= spatial_christoffel[l][i][j] * tensor[l][k];
          }
          else if (indices[0] == false && indices[1] == true) {
            covariant_der[i][j][k] += spatial_christoffel[j][i][l] * tensor[l][k];
            covariant_der[i][j][k] -= spatial_christoffel[l][i][k] * tensor[j][l];
          }
          else if (indices[0] == false && indices[1] == false) {
            covariant_der[i][j][k] += spatial_christoffel[j][i][l] * tensor[l][k];
            covariant_der[i][j][k] += spatial_christoffel[k][i][l] * tensor[j][l];
          }
        }
      }
    }
  }

  return covariant_der;
}

double***
gkyl_gr_spacetime_covariant_der_rank2(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, bool* indices, double** tensor, double*** tensor_der)
{
  double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
    
    for (int j = 0; j < 4; j++) {
      spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
    }
  }

  spacetime->spacetime_christoffel_func(spacetime, t, x, y, z, dt, dx, dy, dz, &spacetime_christoffel);

  double ***covariant_der = gkyl_malloc(sizeof(double**[4]));
  for (int i = 0; i < 4; i++) {
    covariant_der[i] = gkyl_malloc(sizeof(double*[4]));

    for (int j = 0; j < 4; j++) {
      covariant_der[i][j] = gkyl_malloc(sizeof(double[4]));

      for (int k = 0; k < 4; k++) {
        covariant_der[i][j][k] = tensor_der[i][j][k];
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          if (indices[0] == true && indices[1] == true) {
            covariant_der[i][j][k] -= spacetime_christoffel[l][i][j] * tensor[l][k];
            covariant_der[i][j][k] -= spacetime_christoffel[l][i][k] * tensor[j][l];
          }
          else if (indices[0] == true && indices[1] == false) {
            covariant_der[i][j][k] += spacetime_christoffel[k][i][l] * tensor[j][l];
            covariant_der[i][j][k] -= spacetime_christoffel[l][i][j] * tensor[l][k];
          }
          else if (indices[0] == false && indices[1] == true) {
            covariant_der[i][j][k] += spacetime_christoffel[j][i][l] * tensor[l][k];
            covariant_der[i][j][k] -= spacetime_christoffel[l][i][k] * tensor[j][l];
          }
          else if (indices[0] == false && indices[1] == false) {
            covariant_der[i][j][k] += spacetime_christoffel[j][i][l] * tensor[l][k];
            covariant_der[i][j][k] += spacetime_christoffel[k][i][l] * tensor[j][l];
          }
        }
      }
    }
  }

  return covariant_der;
}