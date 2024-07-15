#pragma once

#include <math.h>
#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_ref_count.h>

// Object type
struct gkyl_yield_model;

typedef void (*emission_yield_func_t)(double *out, struct gkyl_yield_model *yield,
  double xc[GKYL_MAX_DIM]);

struct gkyl_yield_model {
  int cdim;
  int vdim;
  double mass;
  double charge;
  emission_yield_func_t function;

  struct gkyl_ref_count ref_count; // reference count
};

struct gkyl_yield_furman_pivi {
  struct gkyl_yield_model yield;
  double deltahat_ts;
  double Ehat_ts;
  double t1;
  double t2;
  double t3;
  double t4;
  double s;
};

struct gkyl_yield_model* gkyl_yield_furman_pivi_new(double deltahat_ts, double Ehat_ts, double t1,
  double t2, double t3, double t4, double s);

struct gkyl_yield_model* gkyl_yield_model_acquire(const struct gkyl_yield_model* model);

void gkyl_yield_model_release(const struct gkyl_yield_model* model);
