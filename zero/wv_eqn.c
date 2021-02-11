#include <gkyl_wv_eqn.h>

struct gkyl_wv_eqn*
gkyl_wv_eqn_aquire(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_wv_eqn*) eqn;
}

double
gkyl_wv_eqn_waves(const struct gkyl_wv_eqn *eqn,
  int dir, const double *delta, const double *ql, const double *qr, double *waves, double *speeds)
{
  return eqn->waves_func(eqn, dir, delta, ql, qr, waves, speeds);
}

void
gkyl_wv_eqn_qfluct(const struct gkyl_wv_eqn *eqn,
  int dir, const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq)
{
  eqn->qfluct_func(eqn, dir, ql, qr, waves, speeds, amdq, apdq);
}

void
gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
