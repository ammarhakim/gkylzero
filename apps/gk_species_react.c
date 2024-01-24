#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_react inp, struct gk_react *react)
{
  react->num_react = inp.num_react; 
  // initialize information about reactions from input struct
  for (int i=0; i<react->num_react; ++i) 
    react->react_type[i] = inp.react_type[i];
}

void 
gk_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_react *react)
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
gk_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react)
{

}
