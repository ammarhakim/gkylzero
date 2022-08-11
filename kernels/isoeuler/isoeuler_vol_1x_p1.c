#include <gkyl_isoeuler_kernels.h>
GKYL_CU_DH double isoeuler_vol_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out)
{
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // statevec: [rho, rho ux, rho uy, rho uz].
  // uvar: [ux, uy, uz].
  // out: Incremented output.

  const double *rho = &statevec[0];
  const double *rhou0 = &statevec[2];
  const double *rhou1 = &statevec[4];
  const double *rhou2 = &statevec[6];
  const double *uvar0 = &uvar[0];
  const double *uvar1 = &uvar[2];
  const double *uvar2 = &uvar[4];
  double *outrho = &out[0];
  double *outrhou0 = &out[2];
  double *outrhou1 = &out[4];
  double *outrhou2 = &out[6];
  double dx10 = 2./dxv[0];

  double vthsq = vth*vth;
  outrho[0] += rho[1]*rhou0[1]*dx10+rho[0]*rhou0[0]*dx10;

  outrhou0[0] += 1.732050807568877*rho[1]*vthsq+rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10;
  outrhou0[1] += rhou0[1]*uvar0[1]*dx10+rhou0[0]*uvar0[0]*dx10;

  return 0.;
}
