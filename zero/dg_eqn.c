#include <gkyl_dg_eqn.h>

struct gkyl_dg_eqn*
gkyl_dg_eqn_aquire(const struct gkyl_dg_eqn* eqn)
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_dg_eqn*) eqn;
}

void
gkyl_dg_eqn_release(const struct gkyl_dg_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
