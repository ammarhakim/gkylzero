#pragma once

typedef void (*emission_elastic_func_t)(double t, const double *xn, double *fout, void *ctx);

struct gkyl_elastic_model {
  int cdim;
  int vdim;
  double mass;
  double charge;
  emission_elastic_func_t function;
};

struct gkyl_elastic_furman_pivi {
  struct gkyl_elastic_model elastic;
  double P1_inf;
  double P1_hat;
  double E_hat;
  double W;
  double p;
};

struct gkyl_elastic_model* gkyl_elastic_furman_pivi_new(double P1_inf, double P1_hat, double E_hat,
  double W, double p);
