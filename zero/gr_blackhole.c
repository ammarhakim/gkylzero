#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_blackhole.h>

double
blackhole_kerrschildscalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

  double mass = blackhole->mass;
  double spin = blackhole->spin;

  double pos_x = blackhole->pos_x;
  double pos_y = blackhole->pos_y;
  double pos_z = blackhole->pos_z;

  double norm_sq = ((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z));

  double rad_dist = sqrt(0.5 * (norm_sq - ((mass * mass) * (spin * spin)) + sqrt(((norm_sq - ((mass * mass) * (spin * spin))) *
    (norm_sq - ((mass * mass) * (spin * spin)))) + (4.0 * ((mass * mass) * (spin * spin) * ((z - pos_z) * (z - pos_z)))))));
  
  return - (mass * (rad_dist * rad_dist * rad_dist)) / ((rad_dist * rad_dist * rad_dist * rad_dist) + ((mass * mass) * (spin * spin) *
    ((z - pos_z) * (z - pos_z))));
}

double*
blackhole_kerrschildvector(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

  double mass = blackhole->mass;
  double spin = blackhole->spin;

  double pos_x = blackhole->pos_x;
  double pos_y = blackhole->pos_y;
  double pos_z = blackhole->pos_z;

  double norm_sq = ((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z));

  double rad_dist = sqrt(0.5 * (norm_sq - ((mass * mass) * (spin * spin)) + sqrt(((norm_sq - ((mass * mass) * (spin * spin))) *
    (norm_sq - ((mass * mass) * (spin * spin)))) + (4.0 * ((mass * mass) * (spin * spin) * ((z - pos_z) * (z - pos_z)))))));
  
  double *kerrschild_vector = malloc(sizeof(double) * 3);
  kerrschild_vector[0] = - ((rad_dist * (x - pos_x)) + (spin * mass * (y - pos_y))) / ((rad_dist * rad_dist) + ((mass * mass) * (spin * spin)));
  kerrschild_vector[1] = - ((rad_dist * (y - pos_y)) - (spin * mass * (x - pos_x))) / ((rad_dist * rad_dist) + ((mass * mass) * (spin * spin)));
  kerrschild_vector[2] = - (z - pos_z) / rad_dist;

  return kerrschild_vector;
}

static void
blackhole_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        (*spatial_metric_tensor)[i][j] = 1.0;
      }
      else {
        (*spatial_metric_tensor)[i][j] = 0.0;
      }
    }
  }

  double kerrschild_scalar = blackhole_kerrschildscalar(spacetime, x, y, z);
  double *kerrschild_vector = blackhole_kerrschildvector(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_metric_tensor)[i][j] = (*spatial_metric_tensor)[i][j] - (2.0 * kerrschild_scalar * kerrschild_vector[i] * kerrschild_vector[j]);
    }
  }
}

static void
blackhole_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        if (i == 0) {
          (*spacetime_metric_tensor)[i][j] = -1.0;
        }
        else {
          (*spacetime_metric_tensor)[i][j] = 1.0;
        }
      }
      else {
        (*spacetime_metric_tensor)[i][j] = 0.0;
      }
    }
  }
}

static void
blackhole_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  double **spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  blackhole_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
  double spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double **spatial_metric_sq = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    spatial_metric_sq[i] = malloc(sizeof(double) * 3);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
  }

  double **euclidean_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    euclidean_metric[i] = malloc(sizeof(double) * 3);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_inv_metric_tensor)[i][j] = (1.0 / spatial_metric_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }
}

static void
blackhole_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  blackhole_spacetime_metric_tensor(spacetime, t, x, y, z, spacetime_inv_metric_tensor);
}

static void
blackhole_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = malloc(sizeof(double*) * 3);
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = malloc(sizeof(double) * 3);
  }

  blackhole_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
}

static void
blackhole_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  *spacetime_metric_det = -1.0;
}

static void
blackhole_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        (*spatial_metric_tensor_der)[i][j][k] = 0.0;
      }
    }
  }
}

static void
blackhole_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        (*spacetime_metric_tensor_der)[i][j][k] = 0.0;
      }
    }
  }
}

static void
blackhole_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  double kerrschild_scalar = blackhole_kerrschildscalar(spacetime, x, y, z);

  *lapse_function = 1.0 / sqrt(1.0 - (2.0 * kerrschild_scalar));
}

static void
blackhole_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  double kerrschild_scalar = blackhole_kerrschildscalar(spacetime, x, y, z);
  double *kerrschild_vector = blackhole_kerrschildvector(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    (*shift_vector)[i] = ((2.0 * kerrschild_scalar) / (1.0 - (2.0 * kerrschild_scalar))) * kerrschild_vector[i];
  }
}

static void
blackhole_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  for (int i = 0; i < 3; i++) {
    (*lapse_function_der)[i] = 0.0;
  }
}

static void
blackhole_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*shift_vector_der)[i][j] = 0.0;
    }
  }
}

static void
blackhole_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*extrinsic_curvature_tensor)[i][j] = 0.0;
    }
  }
}

static void
blackhole_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

  double mass = blackhole->mass;
  double spin = blackhole->spin;

  double pos_x = blackhole->pos_x;
  double pos_y = blackhole->pos_y;
  double pos_z = blackhole->pos_z;

  double r = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));

  if (r <= (mass * (1.0 + sqrt(1.0 - (spin * spin))))) {
    *in_excision_region = true;
  } else {
    *in_excision_region = false;
  }
}

void
gkyl_gr_blackhole_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_blackhole *gr_blackhole = container_of(base->on_dev, struct gr_blackhole, spacetime);
    gkyl_cu_free(gr_blackhole);
  }

  struct gr_blackhole *gr_blackhole = container_of(base, struct gr_blackhole, spacetime);
  gkyl_free(gr_blackhole);
}

struct gkyl_gr_spacetime*
gkyl_gr_blackhole_new(bool use_gpu, double mass, double spin, double pos_x, double pos_y, double pos_z)
{
  return gkyl_gr_blackhole_inew(&(struct gkyl_gr_blackhole_inp) {
      .use_gpu = use_gpu,
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_blackhole_inew(const struct gkyl_gr_blackhole_inp* inp)
{
  struct gr_blackhole *gr_blackhole = gkyl_malloc(sizeof(struct gr_blackhole));

  gr_blackhole->mass = inp->mass;
  gr_blackhole->spin = inp->spin;

  gr_blackhole->pos_x = inp->pos_x;
  gr_blackhole->pos_y = inp->pos_y;
  gr_blackhole->pos_z = inp->pos_z;

  gr_blackhole->spacetime.spatial_metric_tensor_func = blackhole_spatial_metric_tensor;
  gr_blackhole->spacetime.spacetime_metric_tensor_func = blackhole_spacetime_metric_tensor;

  gr_blackhole->spacetime.spatial_inv_metric_tensor_func = blackhole_spatial_inv_metric_tensor;
  gr_blackhole->spacetime.spacetime_inv_metric_tensor_func = blackhole_spacetime_inv_metric_tensor;

  gr_blackhole->spacetime.spatial_metric_det_func = blackhole_spatial_metric_det;
  gr_blackhole->spacetime.spacetime_metric_det_func = blackhole_spacetime_metric_det;

  gr_blackhole->spacetime.spatial_metric_tensor_der_func = blackhole_spatial_metric_tensor_der;
  gr_blackhole->spacetime.spacetime_metric_tensor_der_func = blackhole_spacetime_metric_tensor_der;

  gr_blackhole->spacetime.lapse_function_func = blackhole_lapse_function;
  gr_blackhole->spacetime.shift_vector_func = blackhole_shift_vector;

  gr_blackhole->spacetime.lapse_function_der_func = blackhole_lapse_function_der;
  gr_blackhole->spacetime.shift_vector_der_func = blackhole_shift_vector_der;

  gr_blackhole->spacetime.extrinsic_curvature_tensor_func = blackhole_extrinsic_curvature_tensor;

  gr_blackhole->spacetime.excision_region_func = blackhole_excision_region;

  gr_blackhole->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_blackhole->spacetime.flags);
  gr_blackhole->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_blackhole_free);
  gr_blackhole->spacetime.on_dev = &gr_blackhole->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_blackhole->spacetime;
}