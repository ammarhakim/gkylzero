#pragma once

#include <math.h>
#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_ref_count.h>

struct gkyl_spectrum_model;

typedef void (*emission_spectrum_dist_func_t)(double t, const double *xn, double *fout, void *ctx);
typedef void (*emission_spectrum_norm_func_t)(double *out, struct gkyl_spectrum_model *spectrum,
  const double *flux, double effective_delta);

struct gkyl_spectrum_model {
  int cdim;
  int vdim;
  double mass;
  double charge;
  emission_spectrum_dist_func_t distribution;
  emission_spectrum_norm_func_t normalization;

  struct gkyl_ref_count ref_count; // reference count
};

struct gkyl_spectrum_chung_everhart {
  struct gkyl_spectrum_model spectrum;
  double phi;
};

struct gkyl_spectrum_model* gkyl_spectrum_chung_everhart_new(double phi);

struct gkyl_spectrum_model* gkyl_spectrum_model_acquire(const struct gkyl_spectrum_model* model);

void gkyl_spectrum_model_release(const struct gkyl_spectrum_model* model);
