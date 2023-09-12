#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_alloc_flags_priv.h>

bool
gkyl_dg_prim_vars_type_is_cu_dev(const struct gkyl_dg_prim_vars_type* pvt)
{
  return GKYL_IS_CU_ALLOC(pvt->flags);
}

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_type_acquire(const struct gkyl_dg_prim_vars_type* pvt)
{
  gkyl_ref_count_inc(&pvt->ref_count);
  return (struct gkyl_dg_prim_vars_type*) pvt;
}

void
gkyl_dg_prim_vars_type_release(const struct gkyl_dg_prim_vars_type* pvt)
{
  gkyl_ref_count_dec(&pvt->ref_count);
}
