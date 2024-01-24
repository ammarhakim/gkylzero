#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_react_init(struct gkyl_gyrokinetic_oneb *app, struct gk_species *s, struct gk_react *react)
{
  
}

void 
gk_species_react_cross_init(struct gkyl_gyrokinetic_oneb *app, struct gk_species *s, struct gk_react *react)
{
  
}

// computes reaction coefficients
void
gk_species_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{

}

// updates the reaction terms in the rhs
void
gk_species_react_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{

}

void 
gk_species_react_release(const struct gkyl_gyrokinetic_oneb *app, const struct gk_react *react)
{

}
