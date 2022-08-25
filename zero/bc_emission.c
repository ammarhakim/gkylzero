#include <gkyl_bc_emission.h>
#include <gkyl_alloc.h>
#include <gkyl_vlasov_priv.h>
#include <assert.h>

void
gkyl_bc_emission_init_gain(struct gkyl_bc_emission *bc, void *bc_params)
{
  struct bc_ext_ctx *params = (struct bc_ext_ctx*) bc_params;
  assert(params->gain);
  bc->gain = params->gain;
}

struct gkyl_bc_emission*
gkyl_bc_emission_new(void *bc_params, enum gkyl_bc_basic_type bctype, bool use_gpu)
{
  struct gkyl_bc_emission *up = gkyl_malloc(sizeof(struct gkyl_bc_emission));
  
  switch (bctype) {
    case GKYL_BC_GAIN:
      gkyl_bc_emission_init_gain(up, bc_params);
      break;
  }
  
  return up;
}

void gkyl_bc_emission_release(struct gkyl_bc_emission *up)
{
  gkyl_free(up);
}
