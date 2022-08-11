#include <gkyl_isoeuler_kernels.h>
GKYL_CU_DH void isoeuler_surfx_1x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)
{
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells
  const double dxl1 = 2.0/dxv[0];
  const double dxr1 = 2.0/dxv[0];
  const double *rhol = &statevecl[0];
  const double *rhou0l = &statevecl[2];
  const double *rhou1l = &statevecl[4];
  const double *rhou2l = &statevecl[6];
  const double *rhoc = &statevecc[0];
  const double *rhou0c = &statevecc[2];
  const double *rhou1c = &statevecc[4];
  const double *rhou2c = &statevecc[6];
  const double *rhor = &statevecr[0];
  const double *rhou0r = &statevecr[2];
  const double *rhou1r = &statevecr[4];
  const double *rhou2r = &statevecr[6];
  const double *uvar0l = &uvarl[0];
  const double *uvar1l = &uvarl[2];
  const double *uvar2l = &uvarl[4];
  const double *uvar0c = &uvarc[0];
  const double *uvar1c = &uvarc[2];
  const double *uvar2c = &uvarc[4];
  const double *uvar0r = &uvarr[0];
  const double *uvar1r = &uvarr[2];
  const double *uvar2r = &uvarr[4];
  double *outrho = &out[0];
  double *outrhoux = &out[2];
  double *outrhouy = &out[4];
  double *outrhouz = &out[6];
  double incr[2];

  double vthsq = vth*vth;
  //FluxRho;
  outrho[0] += 0.4330127018922193*rhou0r[1]*dxr1-0.4330127018922193*rhou0c[1]*dxr1-0.25*rhou0r[0]*dxr1-0.25*rhou0c[0]*dxr1+0.4330127018922193*rhou0l[1]*dxl1-0.4330127018922193*rhou0c[1]*dxl1+0.25*rhou0l[0]*dxl1+0.25*rhou0c[0]*dxl1;
  outrho[1] += 0.75*rhou0r[1]*dxr1-0.75*rhou0c[1]*dxr1-0.4330127018922193*rhou0r[0]*dxr1-0.4330127018922193*rhou0c[0]*dxr1-0.75*rhou0l[1]*dxl1+0.75*rhou0c[1]*dxl1-0.4330127018922193*rhou0l[0]*dxl1-0.4330127018922193*rhou0c[0]*dxl1;


  printf("in outrhoux[0]: %f\n",outrhoux[0]);

  printf("rhou0l: %f, rhou0c: %f, rhou0r:%f\n",rhou0l[0],rhou0c[0],rhou0r[0]);
  printf("uvar0c: %f, uvar0c: %f, uvar0c:%f\n",uvar0c[0],uvar0c[0],uvar0c[0]);

  //FluxRhoUx;
  outrhoux[0] += (-0.3061862178478972*rhor[1]*dxr1*vthsq)-0.3061862178478972*rhoc[1]*dxr1*vthsq+0.3061862178478972*rhol[1]*dxl1*vthsq+0.3061862178478972*rhoc[1]*dxl1*vthsq-0.5303300858899106*rhor[1]*uvar0r[1]*dxr1+0.3061862178478972*rhor[0]*uvar0r[1]*dxr1-0.5303300858899106*rhoc[1]*uvar0c[1]*dxr1-0.3061862178478972*rhoc[0]*uvar0c[1]*dxr1+0.3061862178478972*uvar0r[0]*rhor[1]*dxr1-0.3061862178478972*uvar0c[0]*rhoc[1]*dxr1-0.1767766952966369*rhor[0]*uvar0r[0]*dxr1-0.1767766952966369*rhoc[0]*uvar0c[0]*dxr1+0.5303300858899106*rhol[1]*uvar0l[1]*dxl1+0.3061862178478972*rhol[0]*uvar0l[1]*dxl1+0.5303300858899106*rhoc[1]*uvar0c[1]*dxl1-0.3061862178478972*rhoc[0]*uvar0c[1]*dxl1+0.3061862178478972*uvar0l[0]*rhol[1]*dxl1-0.3061862178478972*uvar0c[0]*rhoc[1]*dxl1+0.1767766952966369*rhol[0]*uvar0l[0]*dxl1+0.1767766952966369*rhoc[0]*uvar0c[0]*dxl1;
  outrhoux[1] += 0.9185586535436916*rhor[1]*dxr1*vthsq-0.9185586535436916*rhoc[1]*dxr1*vthsq-0.9185586535436916*rhol[1]*dxl1*vthsq+0.9185586535436916*rhoc[1]*dxl1*vthsq-0.9185586535436916*rhor[1]*uvar0r[1]*dxr1+0.5303300858899106*rhor[0]*uvar0r[1]*dxr1-0.9185586535436916*rhoc[1]*uvar0c[1]*dxr1-0.5303300858899106*rhoc[0]*uvar0c[1]*dxr1+0.5303300858899106*uvar0r[0]*rhor[1]*dxr1-0.5303300858899106*uvar0c[0]*rhoc[1]*dxr1-0.3061862178478972*rhor[0]*uvar0r[0]*dxr1-0.3061862178478972*rhoc[0]*uvar0c[0]*dxr1-0.9185586535436916*rhol[1]*uvar0l[1]*dxl1-0.5303300858899106*rhol[0]*uvar0l[1]*dxl1-0.9185586535436916*rhoc[1]*uvar0c[1]*dxl1+0.5303300858899106*rhoc[0]*uvar0c[1]*dxl1-0.5303300858899106*uvar0l[0]*rhol[1]*dxl1+0.5303300858899106*uvar0c[0]*rhoc[1]*dxl1-0.3061862178478972*rhol[0]*uvar0l[0]*dxl1-0.3061862178478972*rhoc[0]*uvar0c[0]*dxl1;

  printf("out outrhoux[0]: %f\n",outrhoux[0]);

}
