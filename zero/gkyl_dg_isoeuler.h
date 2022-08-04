#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_isoeuler_auxfields {
  const struct gkyl_array *uvar; // Pointer to covariant components of B-field unit vector.
};

struct gkyl_dg_eqn* gkyl_dg_isoeuler_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const double vth, bool use_gpu);

void gkyl_isoeuler_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_isoeuler_auxfields auxin);
