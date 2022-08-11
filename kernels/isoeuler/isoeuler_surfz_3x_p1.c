#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_3x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[8]; 
  const double *rhou1l = &statevecl[16]; 
  const double *rhou2l = &statevecl[24]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[8]; 
  const double *rhou1c = &statevecc[16]; 
  const double *rhou2c = &statevecc[24]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *rhou2r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar2l = &uvarl[16]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[8]; 
  const double *uvar2c = &uvarc[16]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  const double *uvar2r = &uvarr[16]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[8]; 
  double *outrhouy = &out[16]; 
  double *outrhouz = &out[24]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  outrho[0] += 0.4330127018922193*rhou0r[1]*dxr1-0.4330127018922193*rhou0c[1]*dxr1-0.25*rhou0r[0]*dxr1-0.25*rhou0c[0]*dxr1+0.4330127018922193*rhou0l[1]*dxl1-0.4330127018922193*rhou0c[1]*dxl1+0.25*rhou0l[0]*dxl1+0.25*rhou0c[0]*dxl1; 
  outrho[1] += 0.75*rhou0r[1]*dxr1-0.75*rhou0c[1]*dxr1-0.4330127018922193*rhou0r[0]*dxr1-0.4330127018922193*rhou0c[0]*dxr1-0.75*rhou0l[1]*dxl1+0.75*rhou0c[1]*dxl1-0.4330127018922193*rhou0l[0]*dxl1-0.4330127018922193*rhou0c[0]*dxl1; 
  outrho[2] += 0.4330127018922193*rhou0r[4]*dxr1-0.4330127018922193*rhou0c[4]*dxr1-0.25*rhou0r[2]*dxr1-0.25*rhou0c[2]*dxr1+0.4330127018922193*rhou0l[4]*dxl1-0.4330127018922193*rhou0c[4]*dxl1+0.25*rhou0l[2]*dxl1+0.25*rhou0c[2]*dxl1; 
  outrho[3] += 0.4330127018922193*rhou0r[5]*dxr1-0.4330127018922193*rhou0c[5]*dxr1-0.25*rhou0r[3]*dxr1-0.25*rhou0c[3]*dxr1+0.4330127018922193*rhou0l[5]*dxl1-0.4330127018922193*rhou0c[5]*dxl1+0.25*rhou0l[3]*dxl1+0.25*rhou0c[3]*dxl1; 
  outrho[4] += 0.75*rhou0r[4]*dxr1-0.75*rhou0c[4]*dxr1-0.4330127018922194*rhou0r[2]*dxr1-0.4330127018922194*rhou0c[2]*dxr1-0.75*rhou0l[4]*dxl1+0.75*rhou0c[4]*dxl1-0.4330127018922194*rhou0l[2]*dxl1-0.4330127018922194*rhou0c[2]*dxl1; 
  outrho[5] += 0.75*rhou0r[5]*dxr1-0.75*rhou0c[5]*dxr1-0.4330127018922194*rhou0r[3]*dxr1-0.4330127018922194*rhou0c[3]*dxr1-0.75*rhou0l[5]*dxl1+0.75*rhou0c[5]*dxl1-0.4330127018922194*rhou0l[3]*dxl1-0.4330127018922194*rhou0c[3]*dxl1; 
  outrho[6] += 0.4330127018922194*rhou0r[7]*dxr1-0.4330127018922194*rhou0c[7]*dxr1-0.25*rhou0r[6]*dxr1-0.25*rhou0c[6]*dxr1+0.4330127018922194*rhou0l[7]*dxl1-0.4330127018922194*rhou0c[7]*dxl1+0.25*rhou0l[6]*dxl1+0.25*rhou0c[6]*dxl1; 
  outrho[7] += 0.75*rhou0r[7]*dxr1-0.75*rhou0c[7]*dxr1-0.4330127018922193*rhou0r[6]*dxr1-0.4330127018922193*rhou0c[6]*dxr1-0.75*rhou0l[7]*dxl1+0.75*rhou0c[7]*dxl1-0.4330127018922193*rhou0l[6]*dxl1-0.4330127018922193*rhou0c[6]*dxl1; 

 
  //FluxRhoUx; 
  outrhoux[0] += (-0.6123724356957945*rhor[1]*dxr1*vthsq)-0.6123724356957945*rhoc[1]*dxr1*vthsq+0.6123724356957945*rhol[1]*dxl1*vthsq+0.6123724356957945*rhoc[1]*dxl1*vthsq-1.060660171779821*rhor[7]*uvar0r[7]*dxr1+0.6123724356957945*rhor[6]*uvar0r[7]*dxr1-1.060660171779821*rhoc[7]*uvar0c[7]*dxr1-0.6123724356957945*rhoc[6]*uvar0c[7]*dxr1+0.6123724356957945*uvar0r[6]*rhor[7]*dxr1-0.6123724356957945*uvar0c[6]*rhoc[7]*dxr1-0.3535533905932737*rhor[6]*uvar0r[6]*dxr1-0.3535533905932737*rhoc[6]*uvar0c[6]*dxr1-1.060660171779821*rhor[5]*uvar0r[5]*dxr1+0.6123724356957945*rhor[3]*uvar0r[5]*dxr1-1.060660171779821*rhoc[5]*uvar0c[5]*dxr1-0.6123724356957945*rhoc[3]*uvar0c[5]*dxr1+0.6123724356957945*uvar0r[3]*rhor[5]*dxr1-0.6123724356957945*uvar0c[3]*rhoc[5]*dxr1-1.060660171779821*rhor[4]*uvar0r[4]*dxr1+0.6123724356957945*rhor[2]*uvar0r[4]*dxr1-1.060660171779821*rhoc[4]*uvar0c[4]*dxr1-0.6123724356957945*rhoc[2]*uvar0c[4]*dxr1+0.6123724356957945*uvar0r[2]*rhor[4]*dxr1-0.6123724356957945*uvar0c[2]*rhoc[4]*dxr1-0.3535533905932737*rhor[3]*uvar0r[3]*dxr1-0.3535533905932737*rhoc[3]*uvar0c[3]*dxr1-0.3535533905932737*rhor[2]*uvar0r[2]*dxr1-0.3535533905932737*rhoc[2]*uvar0c[2]*dxr1-1.060660171779821*rhor[1]*uvar0r[1]*dxr1+0.6123724356957945*rhor[0]*uvar0r[1]*dxr1-1.060660171779821*rhoc[1]*uvar0c[1]*dxr1-0.6123724356957945*rhoc[0]*uvar0c[1]*dxr1+0.6123724356957945*uvar0r[0]*rhor[1]*dxr1-0.6123724356957945*uvar0c[0]*rhoc[1]*dxr1-0.3535533905932737*rhor[0]*uvar0r[0]*dxr1-0.3535533905932737*rhoc[0]*uvar0c[0]*dxr1+1.060660171779821*rhol[7]*uvar0l[7]*dxl1+0.6123724356957945*rhol[6]*uvar0l[7]*dxl1+1.060660171779821*rhoc[7]*uvar0c[7]*dxl1-0.6123724356957945*rhoc[6]*uvar0c[7]*dxl1+0.6123724356957945*uvar0l[6]*rhol[7]*dxl1-0.6123724356957945*uvar0c[6]*rhoc[7]*dxl1+0.3535533905932737*rhol[6]*uvar0l[6]*dxl1+0.3535533905932737*rhoc[6]*uvar0c[6]*dxl1+1.060660171779821*rhol[5]*uvar0l[5]*dxl1+0.6123724356957945*rhol[3]*uvar0l[5]*dxl1+1.060660171779821*rhoc[5]*uvar0c[5]*dxl1-0.6123724356957945*rhoc[3]*uvar0c[5]*dxl1+0.6123724356957945*uvar0l[3]*rhol[5]*dxl1-0.6123724356957945*uvar0c[3]*rhoc[5]*dxl1+1.060660171779821*rhol[4]*uvar0l[4]*dxl1+0.6123724356957945*rhol[2]*uvar0l[4]*dxl1+1.060660171779821*rhoc[4]*uvar0c[4]*dxl1-0.6123724356957945*rhoc[2]*uvar0c[4]*dxl1+0.6123724356957945*uvar0l[2]*rhol[4]*dxl1-0.6123724356957945*uvar0c[2]*rhoc[4]*dxl1+0.3535533905932737*rhol[3]*uvar0l[3]*dxl1+0.3535533905932737*rhoc[3]*uvar0c[3]*dxl1+0.3535533905932737*rhol[2]*uvar0l[2]*dxl1+0.3535533905932737*rhoc[2]*uvar0c[2]*dxl1+1.060660171779821*rhol[1]*uvar0l[1]*dxl1+0.6123724356957945*rhol[0]*uvar0l[1]*dxl1+1.060660171779821*rhoc[1]*uvar0c[1]*dxl1-0.6123724356957945*rhoc[0]*uvar0c[1]*dxl1+0.6123724356957945*uvar0l[0]*rhol[1]*dxl1-0.6123724356957945*uvar0c[0]*rhoc[1]*dxl1+0.3535533905932737*rhol[0]*uvar0l[0]*dxl1+0.3535533905932737*rhoc[0]*uvar0c[0]*dxl1; 
  outrhoux[1] += 1.837117307087383*rhor[1]*dxr1*vthsq-1.837117307087383*rhoc[1]*dxr1*vthsq-1.837117307087383*rhol[1]*dxl1*vthsq+1.837117307087383*rhoc[1]*dxl1*vthsq-1.837117307087383*rhor[7]*uvar0r[7]*dxr1+1.060660171779821*rhor[6]*uvar0r[7]*dxr1-1.837117307087383*rhoc[7]*uvar0c[7]*dxr1-1.060660171779821*rhoc[6]*uvar0c[7]*dxr1+1.060660171779821*uvar0r[6]*rhor[7]*dxr1-1.060660171779821*uvar0c[6]*rhoc[7]*dxr1-0.6123724356957945*rhor[6]*uvar0r[6]*dxr1-0.6123724356957945*rhoc[6]*uvar0c[6]*dxr1-1.837117307087383*rhor[5]*uvar0r[5]*dxr1+1.060660171779821*rhor[3]*uvar0r[5]*dxr1-1.837117307087383*rhoc[5]*uvar0c[5]*dxr1-1.060660171779821*rhoc[3]*uvar0c[5]*dxr1+1.060660171779821*uvar0r[3]*rhor[5]*dxr1-1.060660171779821*uvar0c[3]*rhoc[5]*dxr1-1.837117307087383*rhor[4]*uvar0r[4]*dxr1+1.060660171779821*rhor[2]*uvar0r[4]*dxr1-1.837117307087383*rhoc[4]*uvar0c[4]*dxr1-1.060660171779821*rhoc[2]*uvar0c[4]*dxr1+1.060660171779821*uvar0r[2]*rhor[4]*dxr1-1.060660171779821*uvar0c[2]*rhoc[4]*dxr1-0.6123724356957945*rhor[3]*uvar0r[3]*dxr1-0.6123724356957945*rhoc[3]*uvar0c[3]*dxr1-0.6123724356957945*rhor[2]*uvar0r[2]*dxr1-0.6123724356957945*rhoc[2]*uvar0c[2]*dxr1-1.837117307087383*rhor[1]*uvar0r[1]*dxr1+1.060660171779821*rhor[0]*uvar0r[1]*dxr1-1.837117307087383*rhoc[1]*uvar0c[1]*dxr1-1.060660171779821*rhoc[0]*uvar0c[1]*dxr1+1.060660171779821*uvar0r[0]*rhor[1]*dxr1-1.060660171779821*uvar0c[0]*rhoc[1]*dxr1-0.6123724356957945*rhor[0]*uvar0r[0]*dxr1-0.6123724356957945*rhoc[0]*uvar0c[0]*dxr1-1.837117307087383*rhol[7]*uvar0l[7]*dxl1-1.060660171779821*rhol[6]*uvar0l[7]*dxl1-1.837117307087383*rhoc[7]*uvar0c[7]*dxl1+1.060660171779821*rhoc[6]*uvar0c[7]*dxl1-1.060660171779821*uvar0l[6]*rhol[7]*dxl1+1.060660171779821*uvar0c[6]*rhoc[7]*dxl1-0.6123724356957945*rhol[6]*uvar0l[6]*dxl1-0.6123724356957945*rhoc[6]*uvar0c[6]*dxl1-1.837117307087383*rhol[5]*uvar0l[5]*dxl1-1.060660171779821*rhol[3]*uvar0l[5]*dxl1-1.837117307087383*rhoc[5]*uvar0c[5]*dxl1+1.060660171779821*rhoc[3]*uvar0c[5]*dxl1-1.060660171779821*uvar0l[3]*rhol[5]*dxl1+1.060660171779821*uvar0c[3]*rhoc[5]*dxl1-1.837117307087383*rhol[4]*uvar0l[4]*dxl1-1.060660171779821*rhol[2]*uvar0l[4]*dxl1-1.837117307087383*rhoc[4]*uvar0c[4]*dxl1+1.060660171779821*rhoc[2]*uvar0c[4]*dxl1-1.060660171779821*uvar0l[2]*rhol[4]*dxl1+1.060660171779821*uvar0c[2]*rhoc[4]*dxl1-0.6123724356957945*rhol[3]*uvar0l[3]*dxl1-0.6123724356957945*rhoc[3]*uvar0c[3]*dxl1-0.6123724356957945*rhol[2]*uvar0l[2]*dxl1-0.6123724356957945*rhoc[2]*uvar0c[2]*dxl1-1.837117307087383*rhol[1]*uvar0l[1]*dxl1-1.060660171779821*rhol[0]*uvar0l[1]*dxl1-1.837117307087383*rhoc[1]*uvar0c[1]*dxl1+1.060660171779821*rhoc[0]*uvar0c[1]*dxl1-1.060660171779821*uvar0l[0]*rhol[1]*dxl1+1.060660171779821*uvar0c[0]*rhoc[1]*dxl1-0.6123724356957945*rhol[0]*uvar0l[0]*dxl1-0.6123724356957945*rhoc[0]*uvar0c[0]*dxl1; 

 
  //FluxRhoUy; 

 
  //FluxRhoUz; 

 
} 
GKYL_CU_DH void isoeuler_surfy_3x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[1]; 
  const double dxr1 = 2.0/dxv[1]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[8]; 
  const double *rhou1l = &statevecl[16]; 
  const double *rhou2l = &statevecl[24]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[8]; 
  const double *rhou1c = &statevecc[16]; 
  const double *rhou2c = &statevecc[24]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *rhou2r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar2l = &uvarl[16]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[8]; 
  const double *uvar2c = &uvarc[16]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  const double *uvar2r = &uvarr[16]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[8]; 
  double *outrhouy = &out[16]; 
  double *outrhouz = &out[24]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  outrho[0] += 0.4330127018922193*rhou1r[2]*dxr1-0.4330127018922193*rhou1c[2]*dxr1-0.25*rhou1r[0]*dxr1-0.25*rhou1c[0]*dxr1+0.4330127018922193*rhou1l[2]*dxl1-0.4330127018922193*rhou1c[2]*dxl1+0.25*rhou1l[0]*dxl1+0.25*rhou1c[0]*dxl1; 
  outrho[1] += 0.4330127018922193*rhou1r[4]*dxr1-0.4330127018922193*rhou1c[4]*dxr1-0.25*rhou1r[1]*dxr1-0.25*rhou1c[1]*dxr1+0.4330127018922193*rhou1l[4]*dxl1-0.4330127018922193*rhou1c[4]*dxl1+0.25*rhou1l[1]*dxl1+0.25*rhou1c[1]*dxl1; 
  outrho[2] += 0.75*rhou1r[2]*dxr1-0.75*rhou1c[2]*dxr1-0.4330127018922193*rhou1r[0]*dxr1-0.4330127018922193*rhou1c[0]*dxr1-0.75*rhou1l[2]*dxl1+0.75*rhou1c[2]*dxl1-0.4330127018922193*rhou1l[0]*dxl1-0.4330127018922193*rhou1c[0]*dxl1; 
  outrho[3] += 0.4330127018922193*rhou1r[6]*dxr1-0.4330127018922193*rhou1c[6]*dxr1-0.25*rhou1r[3]*dxr1-0.25*rhou1c[3]*dxr1+0.4330127018922193*rhou1l[6]*dxl1-0.4330127018922193*rhou1c[6]*dxl1+0.25*rhou1l[3]*dxl1+0.25*rhou1c[3]*dxl1; 
  outrho[4] += 0.75*rhou1r[4]*dxr1-0.75*rhou1c[4]*dxr1-0.4330127018922194*rhou1r[1]*dxr1-0.4330127018922194*rhou1c[1]*dxr1-0.75*rhou1l[4]*dxl1+0.75*rhou1c[4]*dxl1-0.4330127018922194*rhou1l[1]*dxl1-0.4330127018922194*rhou1c[1]*dxl1; 
  outrho[5] += 0.4330127018922194*rhou1r[7]*dxr1-0.4330127018922194*rhou1c[7]*dxr1-0.25*rhou1r[5]*dxr1-0.25*rhou1c[5]*dxr1+0.4330127018922194*rhou1l[7]*dxl1-0.4330127018922194*rhou1c[7]*dxl1+0.25*rhou1l[5]*dxl1+0.25*rhou1c[5]*dxl1; 
  outrho[6] += 0.75*rhou1r[6]*dxr1-0.75*rhou1c[6]*dxr1-0.4330127018922194*rhou1r[3]*dxr1-0.4330127018922194*rhou1c[3]*dxr1-0.75*rhou1l[6]*dxl1+0.75*rhou1c[6]*dxl1-0.4330127018922194*rhou1l[3]*dxl1-0.4330127018922194*rhou1c[3]*dxl1; 
  outrho[7] += 0.75*rhou1r[7]*dxr1-0.75*rhou1c[7]*dxr1-0.4330127018922193*rhou1r[5]*dxr1-0.4330127018922193*rhou1c[5]*dxr1-0.75*rhou1l[7]*dxl1+0.75*rhou1c[7]*dxl1-0.4330127018922193*rhou1l[5]*dxl1-0.4330127018922193*rhou1c[5]*dxl1; 

 
  //FluxRhoUx; 

 
  //FluxRhoUy; 
  outrhouy[0] += (-0.6123724356957945*rhor[2]*dxr1*vthsq)-0.6123724356957945*rhoc[2]*dxr1*vthsq+0.6123724356957945*rhol[2]*dxl1*vthsq+0.6123724356957945*rhoc[2]*dxl1*vthsq-1.060660171779821*rhor[7]*uvar1r[7]*dxr1+0.6123724356957945*rhor[5]*uvar1r[7]*dxr1-1.060660171779821*rhoc[7]*uvar1c[7]*dxr1-0.6123724356957945*rhoc[5]*uvar1c[7]*dxr1+0.6123724356957945*uvar1r[5]*rhor[7]*dxr1-0.6123724356957945*uvar1c[5]*rhoc[7]*dxr1-1.060660171779821*rhor[6]*uvar1r[6]*dxr1+0.6123724356957945*rhor[3]*uvar1r[6]*dxr1-1.060660171779821*rhoc[6]*uvar1c[6]*dxr1-0.6123724356957945*rhoc[3]*uvar1c[6]*dxr1+0.6123724356957945*uvar1r[3]*rhor[6]*dxr1-0.6123724356957945*uvar1c[3]*rhoc[6]*dxr1-0.3535533905932737*rhor[5]*uvar1r[5]*dxr1-0.3535533905932737*rhoc[5]*uvar1c[5]*dxr1-1.060660171779821*rhor[4]*uvar1r[4]*dxr1+0.6123724356957945*rhor[1]*uvar1r[4]*dxr1-1.060660171779821*rhoc[4]*uvar1c[4]*dxr1-0.6123724356957945*rhoc[1]*uvar1c[4]*dxr1+0.6123724356957945*uvar1r[1]*rhor[4]*dxr1-0.6123724356957945*uvar1c[1]*rhoc[4]*dxr1-0.3535533905932737*rhor[3]*uvar1r[3]*dxr1-0.3535533905932737*rhoc[3]*uvar1c[3]*dxr1-1.060660171779821*rhor[2]*uvar1r[2]*dxr1+0.6123724356957945*rhor[0]*uvar1r[2]*dxr1-1.060660171779821*rhoc[2]*uvar1c[2]*dxr1-0.6123724356957945*rhoc[0]*uvar1c[2]*dxr1+0.6123724356957945*uvar1r[0]*rhor[2]*dxr1-0.6123724356957945*uvar1c[0]*rhoc[2]*dxr1-0.3535533905932737*rhor[1]*uvar1r[1]*dxr1-0.3535533905932737*rhoc[1]*uvar1c[1]*dxr1-0.3535533905932737*rhor[0]*uvar1r[0]*dxr1-0.3535533905932737*rhoc[0]*uvar1c[0]*dxr1+1.060660171779821*rhol[7]*uvar1l[7]*dxl1+0.6123724356957945*rhol[5]*uvar1l[7]*dxl1+1.060660171779821*rhoc[7]*uvar1c[7]*dxl1-0.6123724356957945*rhoc[5]*uvar1c[7]*dxl1+0.6123724356957945*uvar1l[5]*rhol[7]*dxl1-0.6123724356957945*uvar1c[5]*rhoc[7]*dxl1+1.060660171779821*rhol[6]*uvar1l[6]*dxl1+0.6123724356957945*rhol[3]*uvar1l[6]*dxl1+1.060660171779821*rhoc[6]*uvar1c[6]*dxl1-0.6123724356957945*rhoc[3]*uvar1c[6]*dxl1+0.6123724356957945*uvar1l[3]*rhol[6]*dxl1-0.6123724356957945*uvar1c[3]*rhoc[6]*dxl1+0.3535533905932737*rhol[5]*uvar1l[5]*dxl1+0.3535533905932737*rhoc[5]*uvar1c[5]*dxl1+1.060660171779821*rhol[4]*uvar1l[4]*dxl1+0.6123724356957945*rhol[1]*uvar1l[4]*dxl1+1.060660171779821*rhoc[4]*uvar1c[4]*dxl1-0.6123724356957945*rhoc[1]*uvar1c[4]*dxl1+0.6123724356957945*uvar1l[1]*rhol[4]*dxl1-0.6123724356957945*uvar1c[1]*rhoc[4]*dxl1+0.3535533905932737*rhol[3]*uvar1l[3]*dxl1+0.3535533905932737*rhoc[3]*uvar1c[3]*dxl1+1.060660171779821*rhol[2]*uvar1l[2]*dxl1+0.6123724356957945*rhol[0]*uvar1l[2]*dxl1+1.060660171779821*rhoc[2]*uvar1c[2]*dxl1-0.6123724356957945*rhoc[0]*uvar1c[2]*dxl1+0.6123724356957945*uvar1l[0]*rhol[2]*dxl1-0.6123724356957945*uvar1c[0]*rhoc[2]*dxl1+0.3535533905932737*rhol[1]*uvar1l[1]*dxl1+0.3535533905932737*rhoc[1]*uvar1c[1]*dxl1+0.3535533905932737*rhol[0]*uvar1l[0]*dxl1+0.3535533905932737*rhoc[0]*uvar1c[0]*dxl1; 
  outrhouy[2] += 1.837117307087383*rhor[2]*dxr1*vthsq-1.837117307087383*rhoc[2]*dxr1*vthsq-1.837117307087383*rhol[2]*dxl1*vthsq+1.837117307087383*rhoc[2]*dxl1*vthsq-1.837117307087383*rhor[7]*uvar1r[7]*dxr1+1.060660171779821*rhor[5]*uvar1r[7]*dxr1-1.837117307087383*rhoc[7]*uvar1c[7]*dxr1-1.060660171779821*rhoc[5]*uvar1c[7]*dxr1+1.060660171779821*uvar1r[5]*rhor[7]*dxr1-1.060660171779821*uvar1c[5]*rhoc[7]*dxr1-1.837117307087383*rhor[6]*uvar1r[6]*dxr1+1.060660171779821*rhor[3]*uvar1r[6]*dxr1-1.837117307087383*rhoc[6]*uvar1c[6]*dxr1-1.060660171779821*rhoc[3]*uvar1c[6]*dxr1+1.060660171779821*uvar1r[3]*rhor[6]*dxr1-1.060660171779821*uvar1c[3]*rhoc[6]*dxr1-0.6123724356957945*rhor[5]*uvar1r[5]*dxr1-0.6123724356957945*rhoc[5]*uvar1c[5]*dxr1-1.837117307087383*rhor[4]*uvar1r[4]*dxr1+1.060660171779821*rhor[1]*uvar1r[4]*dxr1-1.837117307087383*rhoc[4]*uvar1c[4]*dxr1-1.060660171779821*rhoc[1]*uvar1c[4]*dxr1+1.060660171779821*uvar1r[1]*rhor[4]*dxr1-1.060660171779821*uvar1c[1]*rhoc[4]*dxr1-0.6123724356957945*rhor[3]*uvar1r[3]*dxr1-0.6123724356957945*rhoc[3]*uvar1c[3]*dxr1-1.837117307087383*rhor[2]*uvar1r[2]*dxr1+1.060660171779821*rhor[0]*uvar1r[2]*dxr1-1.837117307087383*rhoc[2]*uvar1c[2]*dxr1-1.060660171779821*rhoc[0]*uvar1c[2]*dxr1+1.060660171779821*uvar1r[0]*rhor[2]*dxr1-1.060660171779821*uvar1c[0]*rhoc[2]*dxr1-0.6123724356957945*rhor[1]*uvar1r[1]*dxr1-0.6123724356957945*rhoc[1]*uvar1c[1]*dxr1-0.6123724356957945*rhor[0]*uvar1r[0]*dxr1-0.6123724356957945*rhoc[0]*uvar1c[0]*dxr1-1.837117307087383*rhol[7]*uvar1l[7]*dxl1-1.060660171779821*rhol[5]*uvar1l[7]*dxl1-1.837117307087383*rhoc[7]*uvar1c[7]*dxl1+1.060660171779821*rhoc[5]*uvar1c[7]*dxl1-1.060660171779821*uvar1l[5]*rhol[7]*dxl1+1.060660171779821*uvar1c[5]*rhoc[7]*dxl1-1.837117307087383*rhol[6]*uvar1l[6]*dxl1-1.060660171779821*rhol[3]*uvar1l[6]*dxl1-1.837117307087383*rhoc[6]*uvar1c[6]*dxl1+1.060660171779821*rhoc[3]*uvar1c[6]*dxl1-1.060660171779821*uvar1l[3]*rhol[6]*dxl1+1.060660171779821*uvar1c[3]*rhoc[6]*dxl1-0.6123724356957945*rhol[5]*uvar1l[5]*dxl1-0.6123724356957945*rhoc[5]*uvar1c[5]*dxl1-1.837117307087383*rhol[4]*uvar1l[4]*dxl1-1.060660171779821*rhol[1]*uvar1l[4]*dxl1-1.837117307087383*rhoc[4]*uvar1c[4]*dxl1+1.060660171779821*rhoc[1]*uvar1c[4]*dxl1-1.060660171779821*uvar1l[1]*rhol[4]*dxl1+1.060660171779821*uvar1c[1]*rhoc[4]*dxl1-0.6123724356957945*rhol[3]*uvar1l[3]*dxl1-0.6123724356957945*rhoc[3]*uvar1c[3]*dxl1-1.837117307087383*rhol[2]*uvar1l[2]*dxl1-1.060660171779821*rhol[0]*uvar1l[2]*dxl1-1.837117307087383*rhoc[2]*uvar1c[2]*dxl1+1.060660171779821*rhoc[0]*uvar1c[2]*dxl1-1.060660171779821*uvar1l[0]*rhol[2]*dxl1+1.060660171779821*uvar1c[0]*rhoc[2]*dxl1-0.6123724356957945*rhol[1]*uvar1l[1]*dxl1-0.6123724356957945*rhoc[1]*uvar1c[1]*dxl1-0.6123724356957945*rhol[0]*uvar1l[0]*dxl1-0.6123724356957945*rhoc[0]*uvar1c[0]*dxl1; 

 
  //FluxRhoUz; 

 
} 
GKYL_CU_DH void isoeuler_surfz_3x_ser_p1(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[2]; 
  const double dxr1 = 2.0/dxv[2]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[8]; 
  const double *rhou1l = &statevecl[16]; 
  const double *rhou2l = &statevecl[24]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[8]; 
  const double *rhou1c = &statevecc[16]; 
  const double *rhou2c = &statevecc[24]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *rhou2r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar2l = &uvarl[16]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[8]; 
  const double *uvar2c = &uvarc[16]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  const double *uvar2r = &uvarr[16]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[8]; 
  double *outrhouy = &out[16]; 
  double *outrhouz = &out[24]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  outrho[0] += 0.4330127018922193*rhou2r[3]*dxr1-0.4330127018922193*rhou2c[3]*dxr1-0.25*rhou2r[0]*dxr1-0.25*rhou2c[0]*dxr1+0.4330127018922193*rhou2l[3]*dxl1-0.4330127018922193*rhou2c[3]*dxl1+0.25*rhou2l[0]*dxl1+0.25*rhou2c[0]*dxl1; 
  outrho[1] += 0.4330127018922193*rhou2r[5]*dxr1-0.4330127018922193*rhou2c[5]*dxr1-0.25*rhou2r[1]*dxr1-0.25*rhou2c[1]*dxr1+0.4330127018922193*rhou2l[5]*dxl1-0.4330127018922193*rhou2c[5]*dxl1+0.25*rhou2l[1]*dxl1+0.25*rhou2c[1]*dxl1; 
  outrho[2] += 0.4330127018922193*rhou2r[6]*dxr1-0.4330127018922193*rhou2c[6]*dxr1-0.25*rhou2r[2]*dxr1-0.25*rhou2c[2]*dxr1+0.4330127018922193*rhou2l[6]*dxl1-0.4330127018922193*rhou2c[6]*dxl1+0.25*rhou2l[2]*dxl1+0.25*rhou2c[2]*dxl1; 
  outrho[3] += 0.75*rhou2r[3]*dxr1-0.75*rhou2c[3]*dxr1-0.4330127018922193*rhou2r[0]*dxr1-0.4330127018922193*rhou2c[0]*dxr1-0.75*rhou2l[3]*dxl1+0.75*rhou2c[3]*dxl1-0.4330127018922193*rhou2l[0]*dxl1-0.4330127018922193*rhou2c[0]*dxl1; 
  outrho[4] += 0.4330127018922194*rhou2r[7]*dxr1-0.4330127018922194*rhou2c[7]*dxr1-0.25*rhou2r[4]*dxr1-0.25*rhou2c[4]*dxr1+0.4330127018922194*rhou2l[7]*dxl1-0.4330127018922194*rhou2c[7]*dxl1+0.25*rhou2l[4]*dxl1+0.25*rhou2c[4]*dxl1; 
  outrho[5] += 0.75*rhou2r[5]*dxr1-0.75*rhou2c[5]*dxr1-0.4330127018922194*rhou2r[1]*dxr1-0.4330127018922194*rhou2c[1]*dxr1-0.75*rhou2l[5]*dxl1+0.75*rhou2c[5]*dxl1-0.4330127018922194*rhou2l[1]*dxl1-0.4330127018922194*rhou2c[1]*dxl1; 
  outrho[6] += 0.75*rhou2r[6]*dxr1-0.75*rhou2c[6]*dxr1-0.4330127018922194*rhou2r[2]*dxr1-0.4330127018922194*rhou2c[2]*dxr1-0.75*rhou2l[6]*dxl1+0.75*rhou2c[6]*dxl1-0.4330127018922194*rhou2l[2]*dxl1-0.4330127018922194*rhou2c[2]*dxl1; 
  outrho[7] += 0.75*rhou2r[7]*dxr1-0.75*rhou2c[7]*dxr1-0.4330127018922193*rhou2r[4]*dxr1-0.4330127018922193*rhou2c[4]*dxr1-0.75*rhou2l[7]*dxl1+0.75*rhou2c[7]*dxl1-0.4330127018922193*rhou2l[4]*dxl1-0.4330127018922193*rhou2c[4]*dxl1; 

 
  //FluxRhoUx; 

 
  //FluxRhoUy; 

 
  //FluxRhoUz; 
  outrhouz[0] += (-0.6123724356957945*rhor[3]*dxr1*vthsq)-0.6123724356957945*rhoc[3]*dxr1*vthsq+0.6123724356957945*rhol[3]*dxl1*vthsq+0.6123724356957945*rhoc[3]*dxl1*vthsq-1.060660171779821*rhor[7]*uvar2r[7]*dxr1+0.6123724356957945*rhor[4]*uvar2r[7]*dxr1-1.060660171779821*rhoc[7]*uvar2c[7]*dxr1-0.6123724356957945*rhoc[4]*uvar2c[7]*dxr1+0.6123724356957945*uvar2r[4]*rhor[7]*dxr1-0.6123724356957945*uvar2c[4]*rhoc[7]*dxr1-1.060660171779821*rhor[6]*uvar2r[6]*dxr1+0.6123724356957945*rhor[2]*uvar2r[6]*dxr1-1.060660171779821*rhoc[6]*uvar2c[6]*dxr1-0.6123724356957945*rhoc[2]*uvar2c[6]*dxr1+0.6123724356957945*uvar2r[2]*rhor[6]*dxr1-0.6123724356957945*uvar2c[2]*rhoc[6]*dxr1-1.060660171779821*rhor[5]*uvar2r[5]*dxr1+0.6123724356957945*rhor[1]*uvar2r[5]*dxr1-1.060660171779821*rhoc[5]*uvar2c[5]*dxr1-0.6123724356957945*rhoc[1]*uvar2c[5]*dxr1+0.6123724356957945*uvar2r[1]*rhor[5]*dxr1-0.6123724356957945*uvar2c[1]*rhoc[5]*dxr1-0.3535533905932737*rhor[4]*uvar2r[4]*dxr1-0.3535533905932737*rhoc[4]*uvar2c[4]*dxr1-1.060660171779821*rhor[3]*uvar2r[3]*dxr1+0.6123724356957945*rhor[0]*uvar2r[3]*dxr1-1.060660171779821*rhoc[3]*uvar2c[3]*dxr1-0.6123724356957945*rhoc[0]*uvar2c[3]*dxr1+0.6123724356957945*uvar2r[0]*rhor[3]*dxr1-0.6123724356957945*uvar2c[0]*rhoc[3]*dxr1-0.3535533905932737*rhor[2]*uvar2r[2]*dxr1-0.3535533905932737*rhoc[2]*uvar2c[2]*dxr1-0.3535533905932737*rhor[1]*uvar2r[1]*dxr1-0.3535533905932737*rhoc[1]*uvar2c[1]*dxr1-0.3535533905932737*rhor[0]*uvar2r[0]*dxr1-0.3535533905932737*rhoc[0]*uvar2c[0]*dxr1+1.060660171779821*rhol[7]*uvar2l[7]*dxl1+0.6123724356957945*rhol[4]*uvar2l[7]*dxl1+1.060660171779821*rhoc[7]*uvar2c[7]*dxl1-0.6123724356957945*rhoc[4]*uvar2c[7]*dxl1+0.6123724356957945*uvar2l[4]*rhol[7]*dxl1-0.6123724356957945*uvar2c[4]*rhoc[7]*dxl1+1.060660171779821*rhol[6]*uvar2l[6]*dxl1+0.6123724356957945*rhol[2]*uvar2l[6]*dxl1+1.060660171779821*rhoc[6]*uvar2c[6]*dxl1-0.6123724356957945*rhoc[2]*uvar2c[6]*dxl1+0.6123724356957945*uvar2l[2]*rhol[6]*dxl1-0.6123724356957945*uvar2c[2]*rhoc[6]*dxl1+1.060660171779821*rhol[5]*uvar2l[5]*dxl1+0.6123724356957945*rhol[1]*uvar2l[5]*dxl1+1.060660171779821*rhoc[5]*uvar2c[5]*dxl1-0.6123724356957945*rhoc[1]*uvar2c[5]*dxl1+0.6123724356957945*uvar2l[1]*rhol[5]*dxl1-0.6123724356957945*uvar2c[1]*rhoc[5]*dxl1+0.3535533905932737*rhol[4]*uvar2l[4]*dxl1+0.3535533905932737*rhoc[4]*uvar2c[4]*dxl1+1.060660171779821*rhol[3]*uvar2l[3]*dxl1+0.6123724356957945*rhol[0]*uvar2l[3]*dxl1+1.060660171779821*rhoc[3]*uvar2c[3]*dxl1-0.6123724356957945*rhoc[0]*uvar2c[3]*dxl1+0.6123724356957945*uvar2l[0]*rhol[3]*dxl1-0.6123724356957945*uvar2c[0]*rhoc[3]*dxl1+0.3535533905932737*rhol[2]*uvar2l[2]*dxl1+0.3535533905932737*rhoc[2]*uvar2c[2]*dxl1+0.3535533905932737*rhol[1]*uvar2l[1]*dxl1+0.3535533905932737*rhoc[1]*uvar2c[1]*dxl1+0.3535533905932737*rhol[0]*uvar2l[0]*dxl1+0.3535533905932737*rhoc[0]*uvar2c[0]*dxl1; 
  outrhouz[3] += 1.837117307087383*rhor[3]*dxr1*vthsq-1.837117307087383*rhoc[3]*dxr1*vthsq-1.837117307087383*rhol[3]*dxl1*vthsq+1.837117307087383*rhoc[3]*dxl1*vthsq-1.837117307087383*rhor[7]*uvar2r[7]*dxr1+1.060660171779821*rhor[4]*uvar2r[7]*dxr1-1.837117307087383*rhoc[7]*uvar2c[7]*dxr1-1.060660171779821*rhoc[4]*uvar2c[7]*dxr1+1.060660171779821*uvar2r[4]*rhor[7]*dxr1-1.060660171779821*uvar2c[4]*rhoc[7]*dxr1-1.837117307087383*rhor[6]*uvar2r[6]*dxr1+1.060660171779821*rhor[2]*uvar2r[6]*dxr1-1.837117307087383*rhoc[6]*uvar2c[6]*dxr1-1.060660171779821*rhoc[2]*uvar2c[6]*dxr1+1.060660171779821*uvar2r[2]*rhor[6]*dxr1-1.060660171779821*uvar2c[2]*rhoc[6]*dxr1-1.837117307087383*rhor[5]*uvar2r[5]*dxr1+1.060660171779821*rhor[1]*uvar2r[5]*dxr1-1.837117307087383*rhoc[5]*uvar2c[5]*dxr1-1.060660171779821*rhoc[1]*uvar2c[5]*dxr1+1.060660171779821*uvar2r[1]*rhor[5]*dxr1-1.060660171779821*uvar2c[1]*rhoc[5]*dxr1-0.6123724356957945*rhor[4]*uvar2r[4]*dxr1-0.6123724356957945*rhoc[4]*uvar2c[4]*dxr1-1.837117307087383*rhor[3]*uvar2r[3]*dxr1+1.060660171779821*rhor[0]*uvar2r[3]*dxr1-1.837117307087383*rhoc[3]*uvar2c[3]*dxr1-1.060660171779821*rhoc[0]*uvar2c[3]*dxr1+1.060660171779821*uvar2r[0]*rhor[3]*dxr1-1.060660171779821*uvar2c[0]*rhoc[3]*dxr1-0.6123724356957945*rhor[2]*uvar2r[2]*dxr1-0.6123724356957945*rhoc[2]*uvar2c[2]*dxr1-0.6123724356957945*rhor[1]*uvar2r[1]*dxr1-0.6123724356957945*rhoc[1]*uvar2c[1]*dxr1-0.6123724356957945*rhor[0]*uvar2r[0]*dxr1-0.6123724356957945*rhoc[0]*uvar2c[0]*dxr1-1.837117307087383*rhol[7]*uvar2l[7]*dxl1-1.060660171779821*rhol[4]*uvar2l[7]*dxl1-1.837117307087383*rhoc[7]*uvar2c[7]*dxl1+1.060660171779821*rhoc[4]*uvar2c[7]*dxl1-1.060660171779821*uvar2l[4]*rhol[7]*dxl1+1.060660171779821*uvar2c[4]*rhoc[7]*dxl1-1.837117307087383*rhol[6]*uvar2l[6]*dxl1-1.060660171779821*rhol[2]*uvar2l[6]*dxl1-1.837117307087383*rhoc[6]*uvar2c[6]*dxl1+1.060660171779821*rhoc[2]*uvar2c[6]*dxl1-1.060660171779821*uvar2l[2]*rhol[6]*dxl1+1.060660171779821*uvar2c[2]*rhoc[6]*dxl1-1.837117307087383*rhol[5]*uvar2l[5]*dxl1-1.060660171779821*rhol[1]*uvar2l[5]*dxl1-1.837117307087383*rhoc[5]*uvar2c[5]*dxl1+1.060660171779821*rhoc[1]*uvar2c[5]*dxl1-1.060660171779821*uvar2l[1]*rhol[5]*dxl1+1.060660171779821*uvar2c[1]*rhoc[5]*dxl1-0.6123724356957945*rhol[4]*uvar2l[4]*dxl1-0.6123724356957945*rhoc[4]*uvar2c[4]*dxl1-1.837117307087383*rhol[3]*uvar2l[3]*dxl1-1.060660171779821*rhol[0]*uvar2l[3]*dxl1-1.837117307087383*rhoc[3]*uvar2c[3]*dxl1+1.060660171779821*rhoc[0]*uvar2c[3]*dxl1-1.060660171779821*uvar2l[0]*rhol[3]*dxl1+1.060660171779821*uvar2c[0]*rhoc[3]*dxl1-0.6123724356957945*rhol[2]*uvar2l[2]*dxl1-0.6123724356957945*rhoc[2]*uvar2c[2]*dxl1-0.6123724356957945*rhol[1]*uvar2l[1]*dxl1-0.6123724356957945*rhoc[1]*uvar2c[1]*dxl1-0.6123724356957945*rhol[0]*uvar2l[0]*dxl1-0.6123724356957945*rhoc[0]*uvar2c[0]*dxl1; 

 
} 
