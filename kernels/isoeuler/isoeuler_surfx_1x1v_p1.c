#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_1x1v_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[2]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[2]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar0r = &uvarr[0]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[2]; 
  double *outrrho = &out[4]; 
  double *outrrhoux = &out[6]; 
  double incr[2]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.4330127018922193*rhou0r[1])+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = 0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += incr[1]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = 0.3061862178478972*rhor[1]*vthsq+0.3061862178478972*rhol[1]*vthsq+0.5303300858899106*rhor[1]*uvar0r[1]-0.3061862178478972*rhor[0]*uvar0r[1]+0.5303300858899106*rhol[1]*uvar0l[1]+0.3061862178478972*rhol[0]*uvar0l[1]-0.3061862178478972*uvar0r[0]*rhor[1]+0.3061862178478972*uvar0l[0]*rhol[1]+0.1767766952966369*rhor[0]*uvar0r[0]+0.1767766952966369*rhol[0]*uvar0l[0]; 
  incr[1] = 0.9185586535436916*rhor[1]*vthsq-0.9185586535436916*rhol[1]*vthsq-0.9185586535436916*rhor[1]*uvar0r[1]+0.5303300858899106*rhor[0]*uvar0r[1]-0.9185586535436916*rhol[1]*uvar0l[1]-0.5303300858899106*rhol[0]*uvar0l[1]+0.5303300858899106*uvar0r[0]*rhor[1]-0.5303300858899106*uvar0l[0]*rhol[1]-0.3061862178478972*rhor[0]*uvar0r[0]-0.3061862178478972*rhol[0]*uvar0l[0]; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 

  outlrhoux[0] += -1.0*incr[0]*dxl1; 
  outlrhoux[1] += incr[1]*dxl1; 

 
} 
