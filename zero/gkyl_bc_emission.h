#pragma once

#include <gkyl_spectrum_model.h>
#include <gkyl_yield_model.h>
#include <gkyl_elastic_model.h>

// context for use in computing applied acceleration
struct gkyl_bc_emission_ctx {
  int num_species;
  double t_bound;
  bool elastic;
  struct gkyl_spectrum_model *spectrum_model[GKYL_MAX_SPECIES];
  struct gkyl_yield_model *yield_model[GKYL_MAX_SPECIES];
  struct gkyl_elastic_model *elastic_model;
  char in_species[GKYL_MAX_SPECIES][128];
};

struct gkyl_bc_emission_ctx* gkyl_bc_emission_new(int num_species, bool elastic, struct gkyl_spectrum_model *spectrum_model[], struct gkyl_yield_model *yield_model[], struct gkyl_elastic_model *elastic_model, char in_species[][128]);

void gkyl_bc_emission_release(struct gkyl_bc_emission_ctx *ctx);
