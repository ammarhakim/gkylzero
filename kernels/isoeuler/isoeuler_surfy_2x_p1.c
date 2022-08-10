#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_2x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[4]; 
  const double *rhou1l = &statevecl[8]; 
  const double *rhou2l = &statevecl[12]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[4]; 
  const double *rhou1r = &statevecr[8]; 
  const double *rhou2r = &statevecr[12]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[4]; 
  const double *uvar2l = &uvarl[8]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[4]; 
  const double *uvar2r = &uvarr[8]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[4]; 
  double *outlrhouy = &out[8]; 
  double *outlrhouz = &out[12]; 
  double *outrrho = &out[16]; 
  double *outrrhoux = &out[20]; 
  double *outrrhouy = &out[24]; 
  double *outrrhouz = &out[28]; 
  double incr[4]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.4330127018922193*rhou0r[1])+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = 0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 
  incr[2] = (-0.4330127018922193*rhou0r[3])+0.4330127018922193*rhou0l[3]+0.25*rhou0r[2]+0.25*rhou0l[2]; 
  incr[3] = 0.75*rhou0r[3]-0.75*rhou0l[3]-0.4330127018922194*rhou0r[2]-0.4330127018922194*rhou0l[2]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = 0.4330127018922193*rhor[1]*vthsq+0.4330127018922193*rhol[1]*vthsq+0.75*rhor[3]*uvar0r[3]-0.4330127018922193*rhor[2]*uvar0r[3]+0.75*rhol[3]*uvar0l[3]+0.4330127018922193*rhol[2]*uvar0l[3]-0.4330127018922193*uvar0r[2]*rhor[3]+0.4330127018922193*uvar0l[2]*rhol[3]+0.25*rhor[2]*uvar0r[2]+0.25*rhol[2]*uvar0l[2]+0.75*rhor[1]*uvar0r[1]-0.4330127018922193*rhor[0]*uvar0r[1]+0.75*rhol[1]*uvar0l[1]+0.4330127018922193*rhol[0]*uvar0l[1]-0.4330127018922193*uvar0r[0]*rhor[1]+0.4330127018922193*uvar0l[0]*rhol[1]+0.25*rhor[0]*uvar0r[0]+0.25*rhol[0]*uvar0l[0]; 
  incr[1] = 1.299038105676658*rhor[1]*vthsq-1.299038105676658*rhol[1]*vthsq-1.299038105676658*rhor[3]*uvar0r[3]+0.75*rhor[2]*uvar0r[3]-1.299038105676658*rhol[3]*uvar0l[3]-0.75*rhol[2]*uvar0l[3]+0.75*uvar0r[2]*rhor[3]-0.75*uvar0l[2]*rhol[3]-0.4330127018922193*rhor[2]*uvar0r[2]-0.4330127018922193*rhol[2]*uvar0l[2]-1.299038105676658*rhor[1]*uvar0r[1]+0.75*rhor[0]*uvar0r[1]-1.299038105676658*rhol[1]*uvar0l[1]-0.75*rhol[0]*uvar0l[1]+0.75*uvar0r[0]*rhor[1]-0.75*uvar0l[0]*rhol[1]-0.4330127018922193*rhor[0]*uvar0r[0]-0.4330127018922193*rhol[0]*uvar0l[0]; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 

  outlrhoux[0] += -1.0*incr[0]*dxl1; 
  outlrhoux[1] += incr[1]*dxl1; 

 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 


 
} 
GKYL_CU_DH void isoeuler_surfy_2x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[1]; 
  const double dxr1 = 2.0/dxv[1]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[4]; 
  const double *rhou1l = &statevecl[8]; 
  const double *rhou2l = &statevecl[12]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[4]; 
  const double *rhou1r = &statevecr[8]; 
  const double *rhou2r = &statevecr[12]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[4]; 
  const double *uvar2l = &uvarl[8]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[4]; 
  const double *uvar2r = &uvarr[8]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[4]; 
  double *outlrhouy = &out[8]; 
  double *outlrhouz = &out[12]; 
  double *outrrho = &out[16]; 
  double *outrrhoux = &out[20]; 
  double *outrrhouy = &out[24]; 
  double *outrrhouz = &out[28]; 
  double incr[4]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.4330127018922193*rhou1r[2])+0.4330127018922193*rhou1l[2]+0.25*rhou1r[0]+0.25*rhou1l[0]; 
  incr[1] = (-0.4330127018922193*rhou1r[3])+0.4330127018922193*rhou1l[3]+0.25*rhou1r[1]+0.25*rhou1l[1]; 
  incr[2] = 0.75*rhou1r[2]-0.75*rhou1l[2]-0.4330127018922193*rhou1r[0]-0.4330127018922193*rhou1l[0]; 
  incr[3] = 0.75*rhou1r[3]-0.75*rhou1l[3]-0.4330127018922194*rhou1r[1]-0.4330127018922194*rhou1l[1]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += -1.0*incr[1]*dxl1; 
  outlrho[2] += incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 

 
  //FluxRhoUx; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 


 
  //FluxRhoUy; 
  incr[0] = 0.4330127018922193*rhor[2]*vthsq+0.4330127018922193*rhol[2]*vthsq+0.75*rhor[3]*uvar1r[3]-0.4330127018922193*rhor[1]*uvar1r[3]+0.75*rhol[3]*uvar1l[3]+0.4330127018922193*rhol[1]*uvar1l[3]-0.4330127018922193*uvar1r[1]*rhor[3]+0.4330127018922193*uvar1l[1]*rhol[3]+0.75*rhor[2]*uvar1r[2]-0.4330127018922193*rhor[0]*uvar1r[2]+0.75*rhol[2]*uvar1l[2]+0.4330127018922193*rhol[0]*uvar1l[2]-0.4330127018922193*uvar1r[0]*rhor[2]+0.4330127018922193*uvar1l[0]*rhol[2]+0.25*rhor[1]*uvar1r[1]+0.25*rhol[1]*uvar1l[1]+0.25*rhor[0]*uvar1r[0]+0.25*rhol[0]*uvar1l[0]; 
  incr[2] = 1.299038105676658*rhor[2]*vthsq-1.299038105676658*rhol[2]*vthsq-1.299038105676658*rhor[3]*uvar1r[3]+0.75*rhor[1]*uvar1r[3]-1.299038105676658*rhol[3]*uvar1l[3]-0.75*rhol[1]*uvar1l[3]+0.75*uvar1r[1]*rhor[3]-0.75*uvar1l[1]*rhol[3]-1.299038105676658*rhor[2]*uvar1r[2]+0.75*rhor[0]*uvar1r[2]-1.299038105676658*rhol[2]*uvar1l[2]-0.75*rhol[0]*uvar1l[2]+0.75*uvar1r[0]*rhor[2]-0.75*uvar1l[0]*rhol[2]-0.4330127018922193*rhor[1]*uvar1r[1]-0.4330127018922193*rhol[1]*uvar1l[1]-0.4330127018922193*rhor[0]*uvar1r[0]-0.4330127018922193*rhol[0]*uvar1l[0]; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 

  outlrhouy[0] += -1.0*incr[0]*dxl1; 
  outlrhouy[2] += incr[2]*dxl1; 

 
} 
