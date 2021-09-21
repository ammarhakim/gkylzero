#include <gkyl_wv_eqn.h>

// these ensure inline functions are defined only once

extern inline double gkyl_wv_eqn_waves(const struct gkyl_wv_eqn *eqn,
  int dir, const double *delta, const double *ql, const double *qr, double *waves, double *speeds);

extern inline void gkyl_wv_eqn_qfluct(const struct gkyl_wv_eqn *eqn,
  int dir, const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq);

extern inline double gkyl_wv_eqn_max_speed(const struct gkyl_wv_eqn *eqn,
  int dir, const double *q);

extern inline void gkyl_wv_eqn_rotate_to_local(const struct gkyl_wv_eqn *eqn,
  int dir,  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal);

extern inline void gkyl_wv_eqn_rotate_to_global(const struct gkyl_wv_eqn *eqn,
  int dir,  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal);

struct gkyl_wv_eqn*
gkyl_wv_eqn_aquire(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_wv_eqn*) eqn;
}

void
gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
