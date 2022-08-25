#include <gkyl_bc_external.h>
#include <gkyl_alloc.h>
#include <assert.h>

void
gkyl_bc_external_init_gain(struct gkyl_bc_external *bc, void *bc_params)
{
  struct bc_ext_ctx *params = (struct bc_ext_ctx*) bc_params;
  assert(params->gain);
  bc->gain = params->gain;
}

struct gkyl_bc_external*
gkyl_bc_external_new(void *bc_params, enum gkyl_bc_basic_type bctype, bool use_gpu)
{
  struct gkyl_bc_external *up = gkyl_malloc(sizeof(struct gkyl_bc_external));
  
  switch (bctype) {
    case GKYL_BC_GAIN:
      gkyl_bc_external_init_gain(up, bc_params);
      break;
  }
  
  return up;
}

void gkyl_bc_external_release(struct gkyl_bc_external *up)
{
  gkyl_free(up);
}
