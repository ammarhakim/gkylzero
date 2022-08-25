#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_bc_basic.h>
#include <gkyl_vlasov.h>

struct gkyl_bc_external {
  double gain;
};

void gkyl_bc_external_init_gain(struct gkyl_bc_external *bc, void *bc_params);

struct gkyl_bc_external* gkyl_bc_external_new(void *bc_params, enum gkyl_bc_basic_type bctype, bool use_gpu);

void gkyl_bc_external_release(struct gkyl_bc_external *up);
