#pragma once

#include <gkyl_basis.h>
#include <gkyl_mat.h>

// kernel pointer type
typedef void (*self_prim_t)(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0,
  const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq);

struct gkyl_prim_mom {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  self_prim_t kernel; // moment calculation kernel
  //struct gkyl_ref_count ref_count; // reference count
};

struct gkyl_prim_mom* gkyl_prim_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);
