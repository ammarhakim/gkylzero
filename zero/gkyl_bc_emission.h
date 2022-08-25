#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_bc_basic.h>
#include <gkyl_vlasov.h>

struct gkyl_bc_emission {
  double gain;
};

void gkyl_bc_emission_init_gain(struct gkyl_bc_emission *bc, void *bc_params);

struct gkyl_bc_emission* gkyl_bc_emission_new(void *bc_params, enum gkyl_bc_basic_type bctype, bool use_gpu);

void gkyl_bc_emission_release(struct gkyl_bc_emission *up);
