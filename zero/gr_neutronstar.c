#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_neutronstar.h>
#include <gkyl_gr_spacetime_diff.h>

double
neutronstar_A_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double mass_hexadecapole = neutronstar->mass_hexadecapole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (8.0 * (rho * rho) * (z_cylindrical * z_cylindrical) * ((24.0 * (spin * spin) * mass) + (17.0 * (mass * mass) * mass_quadrupole) + (21.0 * mass_hexadecapole))) +
    ((rho * rho * rho * rho) * ((-10.0 * (spin * spin) * mass) + (7.0 * (mass * mass * mass * mass * mass)) + (32.0 * mass_quadrupole * (mass * mass)) - (21.0 * mass_hexadecapole))) +
    (8.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical) * ((20.0 * (spin * spin) * mass) - (7.0 * (mass * mass * mass * mass * mass)) -
      (22.0 * mass_quadrupole * (mass * mass)) - (7.0 * mass_hexadecapole)));
}

double
neutronstar_B_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double spin_octupole = neutronstar->spin_octupole;
  double mass_hexadecapole = neutronstar->mass_hexadecapole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return ((rho * rho * rho * rho) * ((10.0 * (spin * spin) * (mass * mass)) + (10.0 * mass_quadrupole * (mass * mass * mass)) + (21.0 * mass_hexadecapole * mass) +
      (7.0 * (mass_quadrupole * mass_quadrupole)))) +
    (4.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical) * ((-40.0 * (spin * spin) * (mass * mass)) - (14.0 * spin * (spin_octupole * spin_octupole)) +
      (7.0 * (mass * mass * mass * mass * mass * mass)) + (30.0 * mass_quadrupole * (mass * mass * mass)) + (14.0 * mass_hexadecapole * mass) + (7.0 * (mass_quadrupole * mass_quadrupole)))) -
    (4.0 * (rho * rho) * (z_cylindrical * z_cylindrical) * ((27.0 * (spin * spin) * (mass * mass)) - (21.0 * spin * spin_octupole) + (7.0 * (mass * mass * mass * mass * mass * mass)) +
      (48.0 * mass_quadrupole * (mass * mass * mass)) + (42.0 * mass_hexadecapole * mass) + (7.0 * (mass_quadrupole * mass_quadrupole))));
}

double
neutronstar_H_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double spin_octupole = neutronstar->spin_octupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (4.0 * (rho * rho) * (z_cylindrical * z_cylindrical) * ((spin * (mass_quadrupole - (2.0 * (mass * mass * mass)))) - (3.0 * mass * spin_octupole))) +
    ((rho * rho * rho * rho) * ((spin * mass_quadrupole) + (3.0 * mass * spin_octupole)));
}

double
neutronstar_G_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;
  double spin_octupole = neutronstar->spin_octupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (rho * rho) * ((spin * spin * spin) * (- ((rho * rho * rho * rho) + (8.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical)) -
      (12.0 * (rho * rho) * (z_cylindrical * z_cylindrical)))) +
    (spin * mass) * ((((mass * mass * mass) + (2.0 * mass_quadrupole)) * (rho * rho * rho * rho)) - (8.0 * ((3.0 * (mass * mass * mass)) + (2.0 * mass_quadrupole)) *
      (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical)) +
    (4.0 * (((mass * mass * mass) + (10.0 * mass_quadrupole)) * ((rho * rho) * (z_cylindrical * z_cylindrical))))) +
    ((mass * mass) * spin_octupole) * ((3.0 * (rho * rho * rho * rho)) - (40.0 * (z_cylindrical * z_cylindrical * z_cylindrical * z_cylindrical)) +
      (12.0 * (rho * rho) * (z_cylindrical * z_cylindrical))));
}

double
neutronstar_F_scalar(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double spin_octupole = neutronstar->spin_octupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return ((rho * rho * rho * rho) * (spin_octupole - (spin * (mass * mass)))) - ((4.0 * (spin * spin) * (z_cylindrical * z_cylindrical)) * ((spin * (mass * mass)) + spin_octupole));
}

double
neutronstar_f_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  double A_scalar = neutronstar_A_scalar(spacetime, x, y, z);
  double B_scalar = neutronstar_B_scalar(spacetime, x, y, z);

  return 1.0 - ((2.0 * mass) / sqrt((rho * rho) + (z_cylindrical * z_cylindrical))) + ((2.0 * (mass * mass)) / ((rho * rho) + (z_cylindrical * z_cylindrical))) +
    (((mass_quadrupole - (mass * mass * mass)) * (rho * rho)) - (2.0 * ((mass * mass * mass) + mass_quadrupole) * (z_cylindrical * z_cylindrical))) /
      pow((rho * rho) + (z_cylindrical * z_cylindrical), 5.0 / 2.0) +
    ((2.0 * (z_cylindrical * z_cylindrical) * (-(spin * spin) + (mass * mass * mass * mass) + (2.0 * mass_quadrupole * mass))) - (2.0 * mass * mass_quadrupole * (rho * rho))) /
      pow((rho * rho) + (z_cylindrical * z_cylindrical), 3.0) +
    (A_scalar / (28.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 9.0 / 2.0))) + (B_scalar / (14.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 5.0)));
}

double
neutronstar_omega_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  double H_scalar = neutronstar_H_scalar(spacetime, x, y, z);
  double G_scalar = neutronstar_G_scalar(spacetime, x, y, z);
  double F_scalar = neutronstar_F_scalar(spacetime, x, y, z);

  return (-(2.0 * spin * (rho * rho)) / pow((rho * rho) + (z_cylindrical * z_cylindrical), 3.0 / 2.0)) - ((2.0 * spin * mass * (rho * rho)) /
      pow((rho * rho) + (z_cylindrical * z_cylindrical), 2.0)) +
    (F_scalar / pow((rho * rho) + (z_cylindrical * z_cylindrical), 7.0 / 2.0)) + (H_scalar / (2.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 4.0))) +
    (G_scalar / (4.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 11.0 / 2.0)));
}

double
neutronstar_gamma_function(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double mass = neutronstar->mass;
  double spin = neutronstar->spin;

  double mass_quadrupole = neutronstar->mass_quadrupole;

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  return (((rho * rho) * (((spin * spin) * ((rho * rho) - (8.0 * (z_cylindrical * z_cylindrical)))) + (mass * ((mass * mass * mass) +
    (3.0 * mass_quadrupole)) * ((rho * rho) - (4.0 * (z_cylindrical * z_cylindrical)))))) / (4.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 4.0))) -
    (((mass * mass) * (rho * rho)) / (2.0 * pow((rho * rho) + (z_cylindrical * z_cylindrical), 2.0)));
}

static void
neutronstar_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);

  (*spatial_metric_tensor)[0][0] = exp(2.0 * gamma_function) / f_function;
  (*spatial_metric_tensor)[1][1] = ((rho * rho) / f_function) - (f_function * (omega_function * omega_function));
  (*spatial_metric_tensor)[2][2] = exp(2.0 * gamma_function) / f_function;
}

static void
neutronstar_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*spacetime_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);

  (*spacetime_metric_tensor)[0][0] = -f_function;
  (*spacetime_metric_tensor)[1][1] = exp(2.0 * gamma_function) / f_function;
  (*spacetime_metric_tensor)[2][2] = ((rho * rho) / f_function) - (f_function * (omega_function * omega_function));
  (*spacetime_metric_tensor)[3][3] = exp(2.0 * gamma_function) / f_function;

  (*spacetime_metric_tensor)[0][2] = f_function * omega_function;
  (*spacetime_metric_tensor)[2][0] = f_function * omega_function;
}

static void
neutronstar_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*spatial_inv_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);
  
  (*spatial_inv_metric_tensor)[0][0] = exp(-2.0 * gamma_function) * f_function;
  (*spatial_inv_metric_tensor)[1][1] = f_function / ((rho * rho) - ((f_function * f_function) * (omega_function * omega_function)));
  (*spatial_inv_metric_tensor)[2][2] = exp(-2.0 * gamma_function) * f_function;
}

static void
neutronstar_spacetime_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_inv_metric_tensor)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;
  double pos_z = neutronstar->pos_z;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));
  double z_cylindrical = z - pos_z;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      (*spacetime_inv_metric_tensor)[i][j] = 0.0;
    }
  }

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);
  double gamma_function = neutronstar_gamma_function(spacetime, x, y, z);

  (*spacetime_inv_metric_tensor)[0][0] = -(1.0 / f_function) + ((f_function * (omega_function * omega_function)) / (rho * rho));
  (*spacetime_inv_metric_tensor)[1][1] = exp(-2.0 * gamma_function) * f_function;
  (*spacetime_inv_metric_tensor)[2][2] = f_function / (rho * rho);
  (*spacetime_inv_metric_tensor)[3][3] = exp(-2.0 * gamma_function) * f_function;

  (*spacetime_inv_metric_tensor)[0][2] = (f_function * omega_function) / (rho * rho);
  (*spacetime_inv_metric_tensor)[2][0] = (f_function * omega_function) / (rho * rho);
}

static void
neutronstar_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  neutronstar_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static void
neutronstar_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  double spatial_metric_det;
  double lapse_function;
  neutronstar_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);
  neutronstar_lapse_function(spacetime, t, x, y, z, &lapse_function);

  *spacetime_metric_det = - (lapse_function * lapse_function) * spatial_metric_det;
}

static void
neutronstar_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);

  *lapse_function = (sqrt(f_function) * rho) / sqrt((rho * rho) - ((f_function * f_function) * (omega_function * omega_function)));
}

static void
neutronstar_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  const struct gr_neutronstar *neutronstar = container_of(spacetime, struct gr_neutronstar, spacetime);

  double pos_x = neutronstar->pos_x;
  double pos_y = neutronstar->pos_y;

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)));

  double f_function = neutronstar_f_function(spacetime, x, y, z);
  double omega_function = neutronstar_omega_function(spacetime, x, y, z);

  (*shift_vector)[0] = 0.0;
  (*shift_vector)[1] = - ((f_function * f_function) * omega_function) / (-(rho * rho) + ((f_function * f_function) * (omega_function * omega_function)));
  (*shift_vector)[2] = 0.0;
}

void
gkyl_gr_neutronstar_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_neutronstar *gr_neutronstar = container_of(base->on_dev, struct gr_neutronstar, spacetime);
    gkyl_cu_free(gr_neutronstar);
  }

  struct gr_neutronstar *gr_neutronstar = container_of(base, struct gr_neutronstar, spacetime);
  gkyl_free(gr_neutronstar);
}

struct gkyl_gr_spacetime*
gkyl_gr_neutronstar_new(bool use_gpu, double mass, double spin, double mass_quadrupole, double spin_octupole, double mass_hexadecapole, double pos_x, double pos_y, double pos_z)
{
  return gkyl_gr_neutronstar_inew(&(struct gkyl_gr_neutronstar_inp) {
      .use_gpu = use_gpu,
      .mass = mass,
      .spin = spin,
      .mass_quadrupole = mass_quadrupole,
      .spin_octupole = spin_octupole,
      .mass_hexadecapole = mass_hexadecapole,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_neutronstar_inew(const struct gkyl_gr_neutronstar_inp* inp)
{
  struct gr_neutronstar *gr_neutronstar = gkyl_malloc(sizeof(struct gr_neutronstar));

  gr_neutronstar->mass = inp->mass;
  gr_neutronstar->spin = inp->spin;

  gr_neutronstar->pos_x = inp->pos_x;
  gr_neutronstar->pos_y = inp->pos_y;
  gr_neutronstar->pos_z = inp->pos_z;

  gr_neutronstar->spacetime.spatial_metric_tensor_func = neutronstar_spatial_metric_tensor;
  gr_neutronstar->spacetime.spacetime_metric_tensor_func = neutronstar_spacetime_metric_tensor;

  gr_neutronstar->spacetime.spatial_inv_metric_tensor_func = neutronstar_spatial_inv_metric_tensor;
  gr_neutronstar->spacetime.spacetime_inv_metric_tensor_func = neutronstar_spacetime_inv_metric_tensor;

  gr_neutronstar->spacetime.spatial_metric_det_func = neutronstar_spatial_metric_det;
  gr_neutronstar->spacetime.spacetime_metric_det_func = neutronstar_spacetime_metric_det;

  gr_neutronstar->spacetime.lapse_function_func = neutronstar_lapse_function;
  gr_neutronstar->spacetime.shift_vector_func = neutronstar_shift_vector;

  gr_neutronstar->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_neutronstar->spacetime.flags);
  gr_neutronstar->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_neutronstar_free);
  gr_neutronstar->spacetime.on_dev = &gr_neutronstar->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_neutronstar->spacetime;
}