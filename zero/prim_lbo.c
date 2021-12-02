#include <gkyl_prim_lbo.h>

struct gkyl_prim_lbo*
gkyl_prim_lbo_acquire(const struct gkyl_prim_lbo* prim)
{
  gkyl_ref_count_inc(&prim->ref_count);
  return (struct gkyl_prim_lbo*) prim;
}

void
gkyl_prim_lbo_release(const struct gkyl_prim_lbo* prim)
{
  gkyl_ref_count_dec(&prim->ref_count);
}
