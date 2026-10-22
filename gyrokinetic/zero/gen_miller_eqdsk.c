#include <gkyl_gen_miller_eqdsk.h>
#include <gkyl_alloc.h>

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

enum { MILLER_MAX_NODES = 262145, MILLER_CONTOUR_NODES = 4097 };

struct gkyl_gen_miller_eqdsk {
  struct gkyl_gen_miller_eqdsk_inp inp;
  struct gkyl_gen_miller_eqdsk_info info;
  double eta, fpol, s_guard;
  double *s, *g, *slope, *flux, *mid, *new_s, *new_g, *new_mid;
  double *psi_rz, *qpsi, *guard;
  bool sampled;
};

/**
 * Check IEEE binary64 exponent bits for finiteness under the default -ffast-math build.
 *
 * @param value Floating-point value to classify.
 * @return True for finite values, false for NaN or infinity.
 */
static bool
gen_miller_finite(double value)
{
  uint64_t bits;
  memcpy(&bits, &value, sizeof bits);
  return (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}

/**
 * Construct a quiet NaN from its IEEE binary64 bit representation.
 *
 * @return Quiet NaN used to signal an invalid evaluation.
 */
static double
gen_miller_nan(void)
{
  uint64_t bits = UINT64_C(0x7ff8000000000000);
  double value;
  memcpy(&value, &bits, sizeof value);
  return value;
}

/**
 * Report a Miller generation failure to stderr.
 *
 * @param message Null-terminated diagnostic message.
 * @return False, for use directly in a failure return statement.
 */
static bool
gen_miller_error(const char *message)
{
  fprintf(stderr, "Miller GEQDSK: %s\n", message);
  return false;
}

/**
 * Map Miller minor radius to outboard midplane major radius, including the quadratic shift.
 *
 * @param up Generator.
 * @param radius Miller minor radius in meters.
 * @return Outboard midplane major radius in meters.
 */
static double
gen_miller_R_omp_at_r(const gkyl_gen_miller_eqdsk *up, double radius)
{
  return up->inp.r_axis - up->inp.a_shift * radius * radius / (2 * up->inp.r_axis) + radius;
}

/**
 * Invert the increasing outboard midplane map using a cancellation-resistant quadratic root.
 *
 * @param up Generator.
 * @param r_omp Outboard midplane major radius in meters.
 * @return Miller minor radius in meters, or NaN for an invalid radius or branch.
 */
static double
gen_miller_r_at_R_omp(const gkyl_gen_miller_eqdsk *up, double r_omp)
{
  double offset = r_omp - up->inp.r_axis;
  double disc = 1 - 2 * up->inp.a_shift * offset / up->inp.r_axis;
  if (!gen_miller_finite(offset) || offset < 0 || disc <= 0) {
    return gen_miller_nan();
  }

  return 2 * offset / (1 + sqrt(disc));
}

void
gkyl_gen_miller_eqdsk_RZ_at_rtheta(
  const gkyl_gen_miller_eqdsk *up, double radius, double theta, double rz[2]
)
{
  rz[0] = up->inp.r_axis - up->inp.a_shift * radius * radius / (2 * up->inp.r_axis) +
          radius * cos(theta + up->eta * sin(theta));
  rz[1] = up->inp.z_axis + up->inp.kappa * radius * sin(theta);
}

/**
 * Evaluate (R_r*Z_theta - R_theta*Z_r)/r, including its finite axis limit.
 *
 * @param up Generator.
 * @param radius Miller minor radius in meters.
 * @param theta Poloidal angle in radians.
 * @return Dimensionless R-Z map determinant divided by minor radius.
 */
static double
gen_miller_RZdet_over_r(const gkyl_gen_miller_eqdsk *up, double radius, double theta)
{
  double sint = sin(theta), cost = cos(theta), angle = theta + up->eta * sint;
  return up->inp.kappa * ((cos(angle) - up->inp.a_shift * radius / up->inp.r_axis) * cost +
                          (1 + up->eta * cost) * sin(angle) * sint);
}

/**
 * Evaluate the user safety-factor callback at the outboard major radius of a surface.
 *
 * @param up Generator containing the callback and its context.
 * @param radius Miller minor radius in meters.
 * @return Callback value of q; validity checks are performed by the caller.
 */
static double
gen_miller_q_at_r(const gkyl_gen_miller_eqdsk *up, double radius)
{
  return up->inp.qprofile(gen_miller_R_omp_at_r(up, radius), up->inp.qctx);
}

/**
 * Compute dpsi/ds by refining periodic trapezoidal quadrature in poloidal angle.
 * The smooth angular integrand permits spectral convergence; successive integrals
 * are compared against an explicit tolerance.
 *
 * @param up Generator containing the geometry, q profile, and quadrature tolerance.
 * @param s Squared Miller minor radius in square meters.
 * @return dpsi/ds [Wb/(rad m^2)], or NaN for invalid q, geometry, or failed convergence.
 */
static double
gen_miller_dpsids_at_s(const gkyl_gen_miller_eqdsk *up, double s)
{
  double radius = sqrt(s), q = gen_miller_q_at_r(up, radius), previous = 0;
  if (!gen_miller_finite(q) || q <= 0) {
    return gen_miller_nan();
  }

  for (int count = 32; count <= 16384; count *= 2) {
    double sum = 0;
    for (int j = 0; j < count; ++j) {
      double theta = 2 * M_PI * j / count, rz[2];
      gkyl_gen_miller_eqdsk_RZ_at_rtheta(up, radius, theta, rz);
      double det = gen_miller_RZdet_over_r(up, radius, theta);
      if (rz[0] <= 0 || det <= 0) {
        return gen_miller_nan();
      }
      sum += up->inp.r_axis * det / rz[0];
    }
    double integral = sum * 2 * M_PI / count;
    if (count > 32 && fabs(integral - previous) <= up->inp.rtol * 0.05 * fmax(1, fabs(integral))) {
      return up->fpol * integral / (4 * M_PI * up->inp.r_axis * q);
    }
    previous = integral;
  }

  return gen_miller_nan();
}

/**
 * Limit a one-sided PCHIP endpoint derivative to preserve the local profile shape.
 *
 * @param h0 Width of the interval adjacent to the endpoint.
 * @param h1 Width of the next interval toward the interior.
 * @param m0 Secant slope on the endpoint interval.
 * @param m1 Secant slope on the next interval.
 * @return Limited derivative at the endpoint.
 */
static double
gen_miller_endpoint_slope(double h0, double h1, double m0, double m1)
{
  double result = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
  if (result * m0 <= 0) {
    return 0;
  }

  if (m0 * m1 <= 0 && fabs(result) > fabs(3 * m0)) {
    return 3 * m0;
  }

  return result;
}

/**
 * Compute shape-preserving PCHIP derivatives for the current dpsi/ds table.
 *
 * @param up Generator with populated s and g arrays; its slope array is updated.
 */
static void
gen_miller_slopes(gkyl_gen_miller_eqdsk *up)
{
  int count = up->info.radial_nodes;
  for (int j = 1; j < count - 1; ++j) {
    double h0 = up->s[j] - up->s[j - 1], h1 = up->s[j + 1] - up->s[j];
    double m0 = (up->g[j] - up->g[j - 1]) / h0, m1 = (up->g[j + 1] - up->g[j]) / h1;
    double w0 = 2 * h1 + h0, w1 = h1 + 2 * h0;
    up->slope[j] = m0 * m1 <= 0 ? 0 : (w0 + w1) / (w0 / m0 + w1 / m1);
  }

  double h0 = up->s[1] - up->s[0], h1 = up->s[2] - up->s[1];
  up->slope[0] =
    gen_miller_endpoint_slope(h0, h1, (up->g[1] - up->g[0]) / h0, (up->g[2] - up->g[1]) / h1);
  h0 = up->s[count - 1] - up->s[count - 2];
  h1 = up->s[count - 2] - up->s[count - 3];
  up->slope[count - 1] = gen_miller_endpoint_slope(
    h0, h1, (up->g[count - 1] - up->g[count - 2]) / h0, (up->g[count - 2] - up->g[count - 3]) / h1
  );
}

/**
 * Evaluate the cubic Hermite interpolant of dpsi/ds at an interval midpoint.
 *
 * @param up Generator with a radial table and PCHIP slopes.
 * @param index Index of the interval between s[index] and s[index+1].
 * @return Interpolated dpsi/ds at the midpoint.
 */
static double
gen_miller_midpoint(const gkyl_gen_miller_eqdsk *up, int index)
{
  return 0.5 * (up->g[index] + up->g[index + 1]) +
         (up->s[index + 1] - up->s[index]) * (up->slope[index] - up->slope[index + 1]) / 8;
}

/**
 * Integrate the cubic Hermite interpolant of dpsi/ds over part of one interval.
 *
 * @param up Generator with a radial table and PCHIP slopes.
 * @param index Index of the interval to integrate.
 * @param x Offset in s from the interval start, between zero and the interval width.
 * @return Flux increment from s[index] to s[index]+x [Wb/rad].
 */
static double
gen_miller_integral(const gkyl_gen_miller_eqdsk *up, int index, double x)
{
  double h = up->s[index + 1] - up->s[index], t = x / h;
  double left = up->g[index], right = up->g[index + 1];
  double dl = up->slope[index] * h, dr = up->slope[index + 1] * h;
  double a2 = 3 * (right - left) - 2 * dl - dr, a3 = 2 * (left - right) + dl + dr;
  return h * t * (left + t * (dl / 2 + t * (a2 / 3 + t * a3 / 4)));
}

/**
 * Adaptively refine the radial dpsi/ds table using midpoint quadrature checks,
 * then integrate its PCHIP representation to obtain flux from the axis.
 *
 * @param up Generator with preallocated table storage; tables and radial metadata are updated.
 * @return True on convergence, false for invalid quadrature or exhausted refinement limits.
 */
static bool
gen_miller_build_table(gkyl_gen_miller_eqdsk *up)
{
  up->info.radial_nodes = 65;

  for (int j = 0; j < 65; ++j) {
    up->s[j] = up->s_guard * pow(j / 64.0, 4);
    up->g[j] = gen_miller_dpsids_at_s(up, up->s[j]);
    up->mid[j] = -1;
    if (!gen_miller_finite(up->g[j]) || up->g[j] <= 0) {
      return false;
    }
  }

  for (int pass = 0; pass < 120; ++pass) {
    gen_miller_slopes(up);
    int count = up->info.radial_nodes, next = 0;
    double gmax = 0, error = 0;
    for (int j = 0; j < count; ++j) {
      gmax = fmax(gmax, up->g[j]);
    }

    for (int j = 0; j < count - 1; ++j) {
      double middle = 0.5 * (up->s[j] + up->s[j + 1]);
      if (up->mid[j] <= 0) {
        up->mid[j] = gen_miller_dpsids_at_s(up, middle);
      }
      if (!gen_miller_finite(up->mid[j]) || up->mid[j] <= 0) {
        return false;
      }
      double relative =
        fabs(gen_miller_midpoint(up, j) - up->mid[j]) / fmax(up->mid[j], gmax * 1e-8);
      error = fmax(error, relative);
      if (next + 3 >= MILLER_MAX_NODES) {
        return false;
      }
      up->new_s[next] = up->s[j];
      up->new_g[next] = up->g[j];
      up->new_mid[next++] = relative > up->inp.rtol ? -1 : up->mid[j];
      if (relative > up->inp.rtol) {
        up->new_s[next] = middle;
        up->new_g[next] = up->mid[j];
        up->new_mid[next++] = -1;
      }
    }

    if (next == count - 1) {
      up->info.radial_error = error;
      up->flux[0] = 0;
      for (int j = 0; j < count - 1; ++j) {
        up->flux[j + 1] = up->flux[j] + gen_miller_integral(up, j, up->s[j + 1] - up->s[j]);
      }
      up->info.psi_guard = up->flux[count - 1];
      return true;
    }

    up->new_s[next] = up->s[count - 1];
    up->new_g[next] = up->g[count - 1];
    up->new_mid[next++] = -1;
    memcpy(up->s, up->new_s, next * sizeof(double));
    memcpy(up->g, up->new_g, next * sizeof(double));
    memcpy(up->mid, up->new_mid, next * sizeof(double));
    up->info.radial_nodes = next;
  }
  return false;
}

/**
 * Evaluate integrated radial flux, with continuation linear in s beyond the guard.
 *
 * @param up Generator with a completed radial flux table.
 * @param s Nonnegative squared Miller minor radius in square meters.
 * @return Poloidal flux divided by 2*pi [Wb/rad].
 */
static double
gen_miller_psi_at_s(const gkyl_gen_miller_eqdsk *up, double s)
{
  int lo = 0, hi = up->info.radial_nodes - 1;
  if (s >= up->s_guard) {
    return up->info.psi_guard + up->g[hi] * (s - up->s_guard);
  }

  while (hi - lo > 1) {
    int middle = (lo + hi) / 2;
    if (up->s[middle] <= s) {
      lo = middle;
    } else {
      hi = middle;
    }
  }

  return up->flux[lo] + gen_miller_integral(up, lo, s - up->s[lo]);
}

double
gkyl_gen_miller_eqdsk_psi_at_R_omp(const gkyl_gen_miller_eqdsk *up, double r_omp)
{
  double radius = gen_miller_r_at_R_omp(up, r_omp);
  if (!gen_miller_finite(radius) || radius > up->info.r_guard * (1 + 4 * DBL_EPSILON)) {
    return gen_miller_nan();
  }
  return gen_miller_psi_at_s(up, radius * radius);
}

double
gkyl_gen_miller_eqdsk_R_omp_at_psi(const gkyl_gen_miller_eqdsk *up, double psi)
{
  if (!gen_miller_finite(psi) || psi < 0 || psi > up->info.psi_guard) {
    return gen_miller_nan();
  }

  if (psi == 0) {
    return up->inp.r_axis;
  }

  double lo = 0, hi = up->s_guard;
  for (int iter = 0; iter < 52; ++iter) {
    double middle = (lo + hi) / 2;
    if (gen_miller_psi_at_s(up, middle) < psi) {
      lo = middle;
    } else {
      hi = middle;
    }
  }

  return gen_miller_R_omp_at_r(up, sqrt((lo + hi) / 2));
}

const struct gkyl_gen_miller_eqdsk_info *
gkyl_gen_miller_eqdsk_get_info(const gkyl_gen_miller_eqdsk *up)
{
  return &up->info;
}

void
gkyl_gen_miller_eqdsk_release(gkyl_gen_miller_eqdsk *up)
{
  if (!up) {
    return;
  }
  gkyl_free(up->s);
  gkyl_free(up->g);
  gkyl_free(up->slope);
  gkyl_free(up->flux);
  gkyl_free(up->mid);
  gkyl_free(up->new_s);
  gkyl_free(up->new_g);
  gkyl_free(up->new_mid);
  gkyl_free(up->psi_rz);
  gkyl_free(up->qpsi);
  gkyl_free(up->guard);
  gkyl_free(up);
}

gkyl_gen_miller_eqdsk *
gkyl_gen_miller_eqdsk_inew(const struct gkyl_gen_miller_eqdsk_inp *inp)
{
  if (!inp || !inp->qprofile) {
    gen_miller_error("Miller geometry requires a qprofile callback");
    return 0;
  }

  gkyl_gen_miller_eqdsk *up = gkyl_calloc(1, sizeof *up);

  up->inp = *inp;
  // Set default values for inputs if not given.
  if (up->inp.kappa == 0) {
    up->inp.kappa = 1;
  }
  if (up->inp.rtol == 0) {
    up->inp.rtol = 1e-10;
  }
  if (up->inp.guard_fraction == 0) {
    up->inp.guard_fraction = 0.05;
  }
  if (up->inp.nr == 0) {
    up->inp.nr = 257;
  }
  if (up->inp.nz == 0) {
    up->inp.nz = 257;
  }

  // Check inputs have finite and reasonable values.
  double params[] = {inp->r_axis, inp->z_axis,  inp->b_axis,  inp->r_lcfs,           up->inp.kappa,
                     inp->delta,  inp->a_shift, up->inp.rtol, up->inp.guard_fraction};
  bool valid = true;
  for (int j = 0; j < sizeof params / sizeof params[0]; ++j) {
    valid = valid && gen_miller_finite(params[j]);
  }
  if (!valid || inp->r_axis <= 0 || inp->b_axis <= 0 || inp->r_lcfs <= inp->r_axis ||
      up->inp.kappa <= 0 || fabs(inp->delta) >= 1 || up->inp.rtol < 1e-13 || up->inp.rtol > 1e-3 ||
      up->inp.guard_fraction <= 0 || up->inp.nr < 5 || up->inp.nz < 5 || up->inp.nr > 999 ||
      up->inp.nz > 999 || up->inp.nr % 2 == 0 || up->inp.nz % 2 == 0) {
    gen_miller_error("Invalid Miller parameters or odd R-Z node counts (5..999)");
    goto fail;
  }

  up->eta = asin(inp->delta);
  up->fpol = inp->r_axis * inp->b_axis;
  up->info.nr = up->inp.nr;
  up->info.nz = up->inp.nz;
  up->info.r_lcfs = gen_miller_r_at_R_omp(up, inp->r_lcfs);

  if (!gen_miller_finite(up->info.r_lcfs)) {
    gen_miller_error("LCFS exceeds the increasing midplane branch");
    goto fail;
  }
  up->info.r_required = up->info.r_lcfs;

  if (inp->has_omp_bounds) {
    double lower = gen_miller_r_at_R_omp(up, inp->omp_bounds[0]);
    double upper = gen_miller_r_at_R_omp(up, inp->omp_bounds[1]);
    if (!gen_miller_finite(lower) || !gen_miller_finite(upper) || lower >= upper) {
      gen_miller_error("Invalid outboard midplane coverage bounds");
      goto fail;
    }
    up->info.r_required = fmax(up->info.r_required, upper);
  }

  up->info.r_guard = up->info.r_required * (1 + up->inp.guard_fraction);
  up->s_guard = up->info.r_guard * up->info.r_guard;

  if (!gen_miller_finite(up->s_guard) || up->info.r_guard <= 0) {
    gen_miller_error("LCFS or guard exceeds the increasing midplane branch");
    goto fail;
  }

  for (int ir = 0; ir <= 256; ++ir) {
    double radius = up->info.r_guard * ir / 256, q = gen_miller_q_at_r(up, radius);
    if (!gen_miller_finite(q) || q <= 0 || 1 - inp->a_shift * radius / inp->r_axis <= 0) {
      gen_miller_error("Nonpositive/nonfinite q or folded midplane map before guard");
      goto fail;
    }
    for (int it = 0; it <= 1024; ++it) {
      double theta = 2 * M_PI * it / 1024, rz[2];
      gkyl_gen_miller_eqdsk_RZ_at_rtheta(up, radius, theta, rz);
      if (rz[0] <= 0 || gen_miller_RZdet_over_r(up, radius, theta) <= 0) {
        gen_miller_error("Physical Miller surfaces fold or cross R=0");
        goto fail;
      }
    }
  }

  size_t bytes = MILLER_MAX_NODES * sizeof(double);
  up->s = gkyl_malloc(bytes);
  up->g = gkyl_malloc(bytes);
  up->slope = gkyl_malloc(bytes);
  up->flux = gkyl_malloc(bytes);
  up->mid = gkyl_malloc(bytes);
  up->new_s = gkyl_malloc(bytes);
  up->new_g = gkyl_malloc(bytes);
  up->new_mid = gkyl_malloc(bytes);
  up->psi_rz = gkyl_malloc((size_t)up->info.nr * up->info.nz * sizeof(double));
  up->qpsi = gkyl_malloc(up->info.nr * sizeof(double));
  up->guard = gkyl_malloc(2 * MILLER_CONTOUR_NODES * sizeof(double));

  if (!gen_miller_build_table(up)) {
    gen_miller_error("Miller flux quadrature/interpolation failed; check q and rtol");
    goto fail;
  }

  up->info.psi_boundary = gen_miller_psi_at_s(up, up->info.r_lcfs * up->info.r_lcfs);

  double bounds[] = {DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX};
  for (int j = 0; j < MILLER_CONTOUR_NODES; ++j) {
    double *rz = up->guard + 2 * j;
    gkyl_gen_miller_eqdsk_RZ_at_rtheta(
      up, up->info.r_guard, 2 * M_PI * j / (MILLER_CONTOUR_NODES - 1), rz
    );
    bounds[0] = fmin(bounds[0], rz[0]);
    bounds[1] = fmax(bounds[1], rz[0]);
    bounds[2] = fmin(bounds[2], rz[1]);
    bounds[3] = fmax(bounds[3], rz[1]);
  }
  double rpad = 0.02 * (bounds[1] - bounds[0]), zpad = 0.02 * (bounds[3] - bounds[2]);
  bounds[0] -= rpad;
  bounds[1] += rpad;
  bounds[2] -= zpad;
  bounds[3] += zpad;
  memcpy(up->info.rz_bounds, inp->has_rz_bounds ? inp->rz_bounds : bounds, sizeof bounds);

  for (int j = 0; j < 4; ++j) {
    valid = valid && gen_miller_finite(up->info.rz_bounds[j]);
  }

  const double *box = up->info.rz_bounds;
  if (!valid || !(0 < box[0] && box[0] < inp->r_axis && inp->r_axis < box[1] &&
                  box[2] < inp->z_axis && inp->z_axis < box[3])) {
    gen_miller_error("R-Z bounds must contain axis and have Rmin>0");
    goto fail;
  }

  return up;

fail:
  gkyl_gen_miller_eqdsk_release(up);
  return 0;
}

/**
 * Evaluate the Miller major-radius map in auxiliary coordinates u and v.
 *
 * @param up Generator.
 * @param u Auxiliary coordinate r*cos(theta) in meters.
 * @param v Auxiliary coordinate r*sin(theta) = (Z-Z_axis)/kappa in meters.
 * @return Cylindrical major radius R in meters.
 */
static double
gen_miller_R_uv(const gkyl_gen_miller_eqdsk *up, double u, double v)
{
  double radius = hypot(u, v), beta = radius > 0 ? up->eta * v / radius : 0;
  return up->inp.r_axis - up->inp.a_shift * radius * radius / (2 * up->inp.r_axis) + u * cos(beta) -
         v * sin(beta);
}

/**
 * Invert the R-Z surface map by bracketed bisection in u at fixed v,
 * checking branch monotonicity and the forward-map residual.
 *
 * @param up Generator.
 * @param R Cylindrical major radius in meters.
 * @param Z Vertical coordinate in meters.
 * @return Miller minor radius in meters, or NaN if inversion fails.
 */
static double
gen_miller_r_at_RZ(const gkyl_gen_miller_eqdsk *up, double R, double Z)
{
  double v = (Z - up->inp.z_axis) / up->inp.kappa, bound = up->info.r_guard;

  int expansion = 0;
  while (gen_miller_R_uv(up, -bound, v) > R || gen_miller_R_uv(up, bound, v) < R) {
    if (++expansion > 12) {
      return gen_miller_nan();
    }
    bound *= 1.5;
  }

  for (int j = 0; j <= 64; ++j) {
    double u = bound * (j / 32.0 - 1), radius = hypot(u, v);
    double beta = radius > 0 ? up->eta * v / radius : 0;
    double dbeta = radius > 0 ? -up->eta * v * u / (radius * radius * radius) : 0;
    double derivative =
      cos(beta) - up->inp.a_shift * u / up->inp.r_axis - (u * sin(beta) + v * cos(beta)) * dbeta;
    if (derivative <= 0) {
      return gen_miller_nan();
    }
  }

  double lo = -bound, hi = bound;
  for (int iter = 0; iter < 50; ++iter) {
    double middle = (lo + hi) / 2;
    if (gen_miller_R_uv(up, middle, v) < R) {
      lo = middle;
    } else {
      hi = middle;
    }
  }

  if (fabs(gen_miller_R_uv(up, (lo + hi) / 2, v) - R) > 1e-11 * up->inp.r_axis) {
    return gen_miller_nan();
  }
  return hypot((lo + hi) / 2, v);
}

bool
gkyl_gen_miller_eqdsk_advance(gkyl_gen_miller_eqdsk *up)
{
  up->sampled = false;
  const double *bounds = up->info.rz_bounds;
  double dr = (bounds[1] - bounds[0]) / (up->info.nr - 1);
  double dz = (bounds[3] - bounds[2]) / (up->info.nz - 1);
  up->info.clearance_cells = DBL_MAX;

  for (int contour = 0; contour < 2; ++contour) {
    double radius = contour == 0 ? up->info.r_lcfs : up->info.r_required;
    for (int j = 0; j < MILLER_CONTOUR_NODES; ++j) {
      double rz[2], distance = DBL_MAX;
      gkyl_gen_miller_eqdsk_RZ_at_rtheta(up, radius, 2 * M_PI * j / (MILLER_CONTOUR_NODES - 1), rz);
      distance = fmin((rz[0] - bounds[0]) / dr, (bounds[1] - rz[0]) / dr);
      distance = fmin(distance, fmin((rz[1] - bounds[2]) / dz, (bounds[3] - rz[1]) / dz));
      // Distances to guard samples suffice here: 4097 angular samples resolve
      // the sub-cell distance well on the supported <=999-node R-Z meshes.
      for (int k = 0; k < MILLER_CONTOUR_NODES; ++k) {
        double diffR = (rz[0] - up->guard[2 * k]) / dr, diffZ = (rz[1] - up->guard[2 * k + 1]) / dz;
        distance = fmin(distance, sqrt(diffR * diffR + diffZ * diffZ));
      }
      up->info.clearance_cells = fmin(up->info.clearance_cells, distance);
    }
  }

  if (up->info.clearance_cells < 3) {
    return gen_miller_error("Need 3 cells of clearance: refine R-Z mesh or enlarge guard/bounds");
  }

  for (int iz = 0; iz < up->info.nz; ++iz) {
    for (int ir = 0; ir < up->info.nr; ++ir) {
      double R = bounds[0] + ir * dr, Z = bounds[2] + iz * dz;
      double radius = gen_miller_r_at_RZ(up, R, Z);
      if (!gen_miller_finite(radius)) {
        fprintf(stderr, "Miller GEQDSK: inverse failed at R=%g Z=%g\n", R, Z);
        return false;
      }
      up->psi_rz[(size_t)iz * up->info.nr + ir] = gen_miller_psi_at_s(up, radius * radius);
    }
  }

  for (int j = 0; j < up->info.nr; ++j) {
    double psi = up->info.psi_boundary * j / (up->info.nr - 1);
    double R = gkyl_gen_miller_eqdsk_R_omp_at_psi(up, psi);
    up->qpsi[j] = up->inp.qprofile(R, up->inp.qctx);
    if (!gen_miller_finite(up->qpsi[j]) || up->qpsi[j] <= 0) {
      return gen_miller_error("Invalid q at a GEQDSK flux node");
    }
  }

  up->sampled = true;

  return true;
}

/**
 * Write finite values as GEQDSK records with up to five E16.9 fields per line.
 *
 * @param stream Open output stream; ownership remains with the caller.
 * @param count Number of values to write.
 * @param values Array containing count values.
 * @return True on success, false for nonfinite values, invalid field widths, or write errors.
 */
static bool
gen_miller_write_values(FILE *stream, size_t count, const double *values)
{
  for (size_t j = 0; j < count; ++j) {
    if (!gen_miller_finite(values[j]) || fprintf(stream, "%16.9E", values[j]) != 16) {
      return false;
    }
    if ((j + 1) % 5 == 0 || j + 1 == count) {
      if (fputc('\n', stream) == EOF) {
        return false;
      }
    }
  }
  return true;
}

bool
gkyl_gen_miller_eqdsk_write(const gkyl_gen_miller_eqdsk *up, const char *filename)
{
  if (!up->sampled) {
    return gen_miller_error("Call Miller advance before writing");
  }
  char temporary[2048];
  int length = snprintf(temporary, sizeof temporary, "%s.XXXXXX", filename);
  if (length < 0 || length >= sizeof temporary) {
    return gen_miller_error("Miller output path is too long");
  }
  int descriptor = mkstemp(temporary);
  if (descriptor < 0) {
    return gen_miller_error("Cannot create Miller output file");
  }
  FILE *stream = fdopen(descriptor, "w");
  if (!stream) {
    close(descriptor);
    unlink(temporary);
    return gen_miller_error("Cannot open Miller output stream");
  }
  const double *bounds = up->info.rz_bounds;
  double header[] = {
    bounds[1] - bounds[0],
    bounds[3] - bounds[2],
    up->inp.r_axis,
    bounds[0],
    (bounds[2] + bounds[3]) / 2,
    up->inp.r_axis,
    up->inp.z_axis,
    0,
    up->info.psi_boundary,
    up->inp.b_axis,
    0,
    0,
    0,
    up->inp.r_axis,
    0,
    up->inp.z_axis,
    0,
    up->info.psi_boundary,
    0,
    0
  };

  bool ok =
    fprintf(stream, "%-48s%4d%4d%4d\n", "GKEYLL SYNTHETIC MILLER", 0, up->info.nr, up->info.nz) > 0;
  ok = ok && gen_miller_write_values(stream, 20, header);

  double profile[999];
  for (int array = 0; array < 4; ++array) {
    for (int j = 0; j < up->info.nr; ++j) {
      profile[j] = array == 0 ? up->fpol : 0;
    }
    ok = ok && gen_miller_write_values(stream, up->info.nr, profile);
  }

  ok = ok && gen_miller_write_values(stream, (size_t)up->info.nr * up->info.nz, up->psi_rz);
  ok = ok && gen_miller_write_values(stream, up->info.nr, up->qpsi);
  ok = ok && fprintf(stream, "%5d%5d\n", 257, 0) > 0;

  double boundary[514];
  for (int j = 0; j < 257; ++j) {
    gkyl_gen_miller_eqdsk_RZ_at_rtheta(up, up->info.r_lcfs, 2 * M_PI * j / 256, boundary + 2 * j);
  }
  ok = ok && gen_miller_write_values(stream, 514, boundary);

  if (fclose(stream) != 0) {
    ok = false;
  }
  if (ok && rename(temporary, filename) != 0) {
    ok = false;
  }
  if (!ok) {
    unlink(temporary);
    return gen_miller_error("Failed to write/replace Miller GEQDSK");
  }

  return true;
}
