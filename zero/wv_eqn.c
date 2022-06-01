#include <gkyl_wv_eqn.h>
#include <gkyl_alloc_flags_priv.h>

bool
gkyl_wvn_eq_is_cu_dev(const struct gkyl_wv_eqn *eqn)
{
  return GKYL_IS_CU_ALLOC(eqn->flags);
}

struct gkyl_wv_eqn*
gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_wv_eqn*) eqn;
}

void
gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
