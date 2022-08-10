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
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *rhou2r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar2l = &uvarl[16]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  const double *uvar2r = &uvarr[16]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[8]; 
  double *outlrhouy = &out[16]; 
  double *outlrhouz = &out[24]; 
  double *outrrho = &out[32]; 
  double *outrrhoux = &out[40]; 
  double *outrrhouy = &out[48]; 
  double *outrrhouz = &out[56]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.4330127018922193*rhou0r[1])+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = 0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 
  incr[2] = (-0.4330127018922193*rhou0r[4])+0.4330127018922193*rhou0l[4]+0.25*rhou0r[2]+0.25*rhou0l[2]; 
  incr[3] = (-0.4330127018922193*rhou0r[5])+0.4330127018922193*rhou0l[5]+0.25*rhou0r[3]+0.25*rhou0l[3]; 
  incr[4] = 0.75*rhou0r[4]-0.75*rhou0l[4]-0.4330127018922194*rhou0r[2]-0.4330127018922194*rhou0l[2]; 
  incr[5] = 0.75*rhou0r[5]-0.75*rhou0l[5]-0.4330127018922194*rhou0r[3]-0.4330127018922194*rhou0l[3]; 
  incr[6] = (-0.4330127018922194*rhou0r[7])+0.4330127018922194*rhou0l[7]+0.25*rhou0r[6]+0.25*rhou0l[6]; 
  incr[7] = 0.75*rhou0r[7]-0.75*rhou0l[7]-0.4330127018922193*rhou0r[6]-0.4330127018922193*rhou0l[6]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += -1.0*incr[3]*dxl1; 
  outlrho[4] += incr[4]*dxl1; 
  outlrho[5] += incr[5]*dxl1; 
  outlrho[6] += -1.0*incr[6]*dxl1; 
  outlrho[7] += incr[7]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = 0.6123724356957945*rhor[1]*vthsq+0.6123724356957945*rhol[1]*vthsq+1.060660171779821*rhor[7]*uvar0r[7]-0.6123724356957945*rhor[6]*uvar0r[7]+1.060660171779821*rhol[7]*uvar0l[7]+0.6123724356957945*rhol[6]*uvar0l[7]-0.6123724356957945*uvar0r[6]*rhor[7]+0.6123724356957945*uvar0l[6]*rhol[7]+0.3535533905932737*rhor[6]*uvar0r[6]+0.3535533905932737*rhol[6]*uvar0l[6]+1.060660171779821*rhor[5]*uvar0r[5]-0.6123724356957945*rhor[3]*uvar0r[5]+1.060660171779821*rhol[5]*uvar0l[5]+0.6123724356957945*rhol[3]*uvar0l[5]-0.6123724356957945*uvar0r[3]*rhor[5]+0.6123724356957945*uvar0l[3]*rhol[5]+1.060660171779821*rhor[4]*uvar0r[4]-0.6123724356957945*rhor[2]*uvar0r[4]+1.060660171779821*rhol[4]*uvar0l[4]+0.6123724356957945*rhol[2]*uvar0l[4]-0.6123724356957945*uvar0r[2]*rhor[4]+0.6123724356957945*uvar0l[2]*rhol[4]+0.3535533905932737*rhor[3]*uvar0r[3]+0.3535533905932737*rhol[3]*uvar0l[3]+0.3535533905932737*rhor[2]*uvar0r[2]+0.3535533905932737*rhol[2]*uvar0l[2]+1.060660171779821*rhor[1]*uvar0r[1]-0.6123724356957945*rhor[0]*uvar0r[1]+1.060660171779821*rhol[1]*uvar0l[1]+0.6123724356957945*rhol[0]*uvar0l[1]-0.6123724356957945*uvar0r[0]*rhor[1]+0.6123724356957945*uvar0l[0]*rhol[1]+0.3535533905932737*rhor[0]*uvar0r[0]+0.3535533905932737*rhol[0]*uvar0l[0]; 
  incr[1] = 1.837117307087383*rhor[1]*vthsq-1.837117307087383*rhol[1]*vthsq-1.837117307087383*rhor[7]*uvar0r[7]+1.060660171779821*rhor[6]*uvar0r[7]-1.837117307087383*rhol[7]*uvar0l[7]-1.060660171779821*rhol[6]*uvar0l[7]+1.060660171779821*uvar0r[6]*rhor[7]-1.060660171779821*uvar0l[6]*rhol[7]-0.6123724356957945*rhor[6]*uvar0r[6]-0.6123724356957945*rhol[6]*uvar0l[6]-1.837117307087383*rhor[5]*uvar0r[5]+1.060660171779821*rhor[3]*uvar0r[5]-1.837117307087383*rhol[5]*uvar0l[5]-1.060660171779821*rhol[3]*uvar0l[5]+1.060660171779821*uvar0r[3]*rhor[5]-1.060660171779821*uvar0l[3]*rhol[5]-1.837117307087383*rhor[4]*uvar0r[4]+1.060660171779821*rhor[2]*uvar0r[4]-1.837117307087383*rhol[4]*uvar0l[4]-1.060660171779821*rhol[2]*uvar0l[4]+1.060660171779821*uvar0r[2]*rhor[4]-1.060660171779821*uvar0l[2]*rhol[4]-0.6123724356957945*rhor[3]*uvar0r[3]-0.6123724356957945*rhol[3]*uvar0l[3]-0.6123724356957945*rhor[2]*uvar0r[2]-0.6123724356957945*rhol[2]*uvar0l[2]-1.837117307087383*rhor[1]*uvar0r[1]+1.060660171779821*rhor[0]*uvar0r[1]-1.837117307087383*rhol[1]*uvar0l[1]-1.060660171779821*rhol[0]*uvar0l[1]+1.060660171779821*uvar0r[0]*rhor[1]-1.060660171779821*uvar0l[0]*rhol[1]-0.6123724356957945*rhor[0]*uvar0r[0]-0.6123724356957945*rhol[0]*uvar0l[0]; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 

  outlrhoux[0] += -1.0*incr[0]*dxl1; 
  outlrhoux[1] += incr[1]*dxl1; 

 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 


 
  //FluxRhoUz; 

  outrrhouz[0] += incr[0]*dxr1; 
  outrrhouz[1] += incr[1]*dxr1; 
  outrrhouz[2] += incr[2]*dxr1; 
  outrrhouz[3] += incr[3]*dxr1; 
  outrrhouz[4] += incr[4]*dxr1; 
  outrrhouz[5] += incr[5]*dxr1; 
  outrrhouz[6] += incr[6]*dxr1; 
  outrrhouz[7] += incr[7]*dxr1; 


 
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
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *rhou2r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar2l = &uvarl[16]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  const double *uvar2r = &uvarr[16]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[8]; 
  double *outlrhouy = &out[16]; 
  double *outlrhouz = &out[24]; 
  double *outrrho = &out[32]; 
  double *outrrhoux = &out[40]; 
  double *outrrhouy = &out[48]; 
  double *outrrhouz = &out[56]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.4330127018922193*rhou1r[2])+0.4330127018922193*rhou1l[2]+0.25*rhou1r[0]+0.25*rhou1l[0]; 
  incr[1] = (-0.4330127018922193*rhou1r[4])+0.4330127018922193*rhou1l[4]+0.25*rhou1r[1]+0.25*rhou1l[1]; 
  incr[2] = 0.75*rhou1r[2]-0.75*rhou1l[2]-0.4330127018922193*rhou1r[0]-0.4330127018922193*rhou1l[0]; 
  incr[3] = (-0.4330127018922193*rhou1r[6])+0.4330127018922193*rhou1l[6]+0.25*rhou1r[3]+0.25*rhou1l[3]; 
  incr[4] = 0.75*rhou1r[4]-0.75*rhou1l[4]-0.4330127018922194*rhou1r[1]-0.4330127018922194*rhou1l[1]; 
  incr[5] = (-0.4330127018922194*rhou1r[7])+0.4330127018922194*rhou1l[7]+0.25*rhou1r[5]+0.25*rhou1l[5]; 
  incr[6] = 0.75*rhou1r[6]-0.75*rhou1l[6]-0.4330127018922194*rhou1r[3]-0.4330127018922194*rhou1l[3]; 
  incr[7] = 0.75*rhou1r[7]-0.75*rhou1l[7]-0.4330127018922193*rhou1r[5]-0.4330127018922193*rhou1l[5]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += -1.0*incr[1]*dxl1; 
  outlrho[2] += incr[2]*dxl1; 
  outlrho[3] += -1.0*incr[3]*dxl1; 
  outlrho[4] += incr[4]*dxl1; 
  outlrho[5] += -1.0*incr[5]*dxl1; 
  outlrho[6] += incr[6]*dxl1; 
  outlrho[7] += incr[7]*dxl1; 

 
  //FluxRhoUx; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 


 
  //FluxRhoUy; 
  incr[0] = 0.6123724356957945*rhor[2]*vthsq+0.6123724356957945*rhol[2]*vthsq+1.060660171779821*rhor[7]*uvar1r[7]-0.6123724356957945*rhor[5]*uvar1r[7]+1.060660171779821*rhol[7]*uvar1l[7]+0.6123724356957945*rhol[5]*uvar1l[7]-0.6123724356957945*uvar1r[5]*rhor[7]+0.6123724356957945*uvar1l[5]*rhol[7]+1.060660171779821*rhor[6]*uvar1r[6]-0.6123724356957945*rhor[3]*uvar1r[6]+1.060660171779821*rhol[6]*uvar1l[6]+0.6123724356957945*rhol[3]*uvar1l[6]-0.6123724356957945*uvar1r[3]*rhor[6]+0.6123724356957945*uvar1l[3]*rhol[6]+0.3535533905932737*rhor[5]*uvar1r[5]+0.3535533905932737*rhol[5]*uvar1l[5]+1.060660171779821*rhor[4]*uvar1r[4]-0.6123724356957945*rhor[1]*uvar1r[4]+1.060660171779821*rhol[4]*uvar1l[4]+0.6123724356957945*rhol[1]*uvar1l[4]-0.6123724356957945*uvar1r[1]*rhor[4]+0.6123724356957945*uvar1l[1]*rhol[4]+0.3535533905932737*rhor[3]*uvar1r[3]+0.3535533905932737*rhol[3]*uvar1l[3]+1.060660171779821*rhor[2]*uvar1r[2]-0.6123724356957945*rhor[0]*uvar1r[2]+1.060660171779821*rhol[2]*uvar1l[2]+0.6123724356957945*rhol[0]*uvar1l[2]-0.6123724356957945*uvar1r[0]*rhor[2]+0.6123724356957945*uvar1l[0]*rhol[2]+0.3535533905932737*rhor[1]*uvar1r[1]+0.3535533905932737*rhol[1]*uvar1l[1]+0.3535533905932737*rhor[0]*uvar1r[0]+0.3535533905932737*rhol[0]*uvar1l[0]; 
  incr[2] = 1.837117307087383*rhor[2]*vthsq-1.837117307087383*rhol[2]*vthsq-1.837117307087383*rhor[7]*uvar1r[7]+1.060660171779821*rhor[5]*uvar1r[7]-1.837117307087383*rhol[7]*uvar1l[7]-1.060660171779821*rhol[5]*uvar1l[7]+1.060660171779821*uvar1r[5]*rhor[7]-1.060660171779821*uvar1l[5]*rhol[7]-1.837117307087383*rhor[6]*uvar1r[6]+1.060660171779821*rhor[3]*uvar1r[6]-1.837117307087383*rhol[6]*uvar1l[6]-1.060660171779821*rhol[3]*uvar1l[6]+1.060660171779821*uvar1r[3]*rhor[6]-1.060660171779821*uvar1l[3]*rhol[6]-0.6123724356957945*rhor[5]*uvar1r[5]-0.6123724356957945*rhol[5]*uvar1l[5]-1.837117307087383*rhor[4]*uvar1r[4]+1.060660171779821*rhor[1]*uvar1r[4]-1.837117307087383*rhol[4]*uvar1l[4]-1.060660171779821*rhol[1]*uvar1l[4]+1.060660171779821*uvar1r[1]*rhor[4]-1.060660171779821*uvar1l[1]*rhol[4]-0.6123724356957945*rhor[3]*uvar1r[3]-0.6123724356957945*rhol[3]*uvar1l[3]-1.837117307087383*rhor[2]*uvar1r[2]+1.060660171779821*rhor[0]*uvar1r[2]-1.837117307087383*rhol[2]*uvar1l[2]-1.060660171779821*rhol[0]*uvar1l[2]+1.060660171779821*uvar1r[0]*rhor[2]-1.060660171779821*uvar1l[0]*rhol[2]-0.6123724356957945*rhor[1]*uvar1r[1]-0.6123724356957945*rhol[1]*uvar1l[1]-0.6123724356957945*rhor[0]*uvar1r[0]-0.6123724356957945*rhol[0]*uvar1l[0]; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 

  outlrhouy[0] += -1.0*incr[0]*dxl1; 
  outlrhouy[2] += incr[2]*dxl1; 

 
  //FluxRhoUz; 

  outrrhouz[0] += incr[0]*dxr1; 
  outrrhouz[1] += incr[1]*dxr1; 
  outrrhouz[2] += incr[2]*dxr1; 
  outrrhouz[3] += incr[3]*dxr1; 
  outrrhouz[4] += incr[4]*dxr1; 
  outrrhouz[5] += incr[5]*dxr1; 
  outrrhouz[6] += incr[6]*dxr1; 
  outrrhouz[7] += incr[7]*dxr1; 


 
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
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *rhou2r = &statevecr[24]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar2l = &uvarl[16]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  const double *uvar2r = &uvarr[16]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[8]; 
  double *outlrhouy = &out[16]; 
  double *outlrhouz = &out[24]; 
  double *outrrho = &out[32]; 
  double *outrrhoux = &out[40]; 
  double *outrrhouy = &out[48]; 
  double *outrrhouz = &out[56]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = (-0.4330127018922193*rhou2r[3])+0.4330127018922193*rhou2l[3]+0.25*rhou2r[0]+0.25*rhou2l[0]; 
  incr[1] = (-0.4330127018922193*rhou2r[5])+0.4330127018922193*rhou2l[5]+0.25*rhou2r[1]+0.25*rhou2l[1]; 
  incr[2] = (-0.4330127018922193*rhou2r[6])+0.4330127018922193*rhou2l[6]+0.25*rhou2r[2]+0.25*rhou2l[2]; 
  incr[3] = 0.75*rhou2r[3]-0.75*rhou2l[3]-0.4330127018922193*rhou2r[0]-0.4330127018922193*rhou2l[0]; 
  incr[4] = (-0.4330127018922194*rhou2r[7])+0.4330127018922194*rhou2l[7]+0.25*rhou2r[4]+0.25*rhou2l[4]; 
  incr[5] = 0.75*rhou2r[5]-0.75*rhou2l[5]-0.4330127018922194*rhou2r[1]-0.4330127018922194*rhou2l[1]; 
  incr[6] = 0.75*rhou2r[6]-0.75*rhou2l[6]-0.4330127018922194*rhou2r[2]-0.4330127018922194*rhou2l[2]; 
  incr[7] = 0.75*rhou2r[7]-0.75*rhou2l[7]-0.4330127018922193*rhou2r[4]-0.4330127018922193*rhou2l[4]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += -1.0*incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 
  outlrho[4] += -1.0*incr[4]*dxl1; 
  outlrho[5] += incr[5]*dxl1; 
  outlrho[6] += incr[6]*dxl1; 
  outlrho[7] += incr[7]*dxl1; 

 
  //FluxRhoUx; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 


 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 


 
  //FluxRhoUz; 
  incr[0] = 0.6123724356957945*rhor[3]*vthsq+0.6123724356957945*rhol[3]*vthsq+1.060660171779821*rhor[7]*uvar2r[7]-0.6123724356957945*rhor[4]*uvar2r[7]+1.060660171779821*rhol[7]*uvar2l[7]+0.6123724356957945*rhol[4]*uvar2l[7]-0.6123724356957945*uvar2r[4]*rhor[7]+0.6123724356957945*uvar2l[4]*rhol[7]+1.060660171779821*rhor[6]*uvar2r[6]-0.6123724356957945*rhor[2]*uvar2r[6]+1.060660171779821*rhol[6]*uvar2l[6]+0.6123724356957945*rhol[2]*uvar2l[6]-0.6123724356957945*uvar2r[2]*rhor[6]+0.6123724356957945*uvar2l[2]*rhol[6]+1.060660171779821*rhor[5]*uvar2r[5]-0.6123724356957945*rhor[1]*uvar2r[5]+1.060660171779821*rhol[5]*uvar2l[5]+0.6123724356957945*rhol[1]*uvar2l[5]-0.6123724356957945*uvar2r[1]*rhor[5]+0.6123724356957945*uvar2l[1]*rhol[5]+0.3535533905932737*rhor[4]*uvar2r[4]+0.3535533905932737*rhol[4]*uvar2l[4]+1.060660171779821*rhor[3]*uvar2r[3]-0.6123724356957945*rhor[0]*uvar2r[3]+1.060660171779821*rhol[3]*uvar2l[3]+0.6123724356957945*rhol[0]*uvar2l[3]-0.6123724356957945*uvar2r[0]*rhor[3]+0.6123724356957945*uvar2l[0]*rhol[3]+0.3535533905932737*rhor[2]*uvar2r[2]+0.3535533905932737*rhol[2]*uvar2l[2]+0.3535533905932737*rhor[1]*uvar2r[1]+0.3535533905932737*rhol[1]*uvar2l[1]+0.3535533905932737*rhor[0]*uvar2r[0]+0.3535533905932737*rhol[0]*uvar2l[0]; 
  incr[3] = 1.837117307087383*rhor[3]*vthsq-1.837117307087383*rhol[3]*vthsq-1.837117307087383*rhor[7]*uvar2r[7]+1.060660171779821*rhor[4]*uvar2r[7]-1.837117307087383*rhol[7]*uvar2l[7]-1.060660171779821*rhol[4]*uvar2l[7]+1.060660171779821*uvar2r[4]*rhor[7]-1.060660171779821*uvar2l[4]*rhol[7]-1.837117307087383*rhor[6]*uvar2r[6]+1.060660171779821*rhor[2]*uvar2r[6]-1.837117307087383*rhol[6]*uvar2l[6]-1.060660171779821*rhol[2]*uvar2l[6]+1.060660171779821*uvar2r[2]*rhor[6]-1.060660171779821*uvar2l[2]*rhol[6]-1.837117307087383*rhor[5]*uvar2r[5]+1.060660171779821*rhor[1]*uvar2r[5]-1.837117307087383*rhol[5]*uvar2l[5]-1.060660171779821*rhol[1]*uvar2l[5]+1.060660171779821*uvar2r[1]*rhor[5]-1.060660171779821*uvar2l[1]*rhol[5]-0.6123724356957945*rhor[4]*uvar2r[4]-0.6123724356957945*rhol[4]*uvar2l[4]-1.837117307087383*rhor[3]*uvar2r[3]+1.060660171779821*rhor[0]*uvar2r[3]-1.837117307087383*rhol[3]*uvar2l[3]-1.060660171779821*rhol[0]*uvar2l[3]+1.060660171779821*uvar2r[0]*rhor[3]-1.060660171779821*uvar2l[0]*rhol[3]-0.6123724356957945*rhor[2]*uvar2r[2]-0.6123724356957945*rhol[2]*uvar2l[2]-0.6123724356957945*rhor[1]*uvar2r[1]-0.6123724356957945*rhol[1]*uvar2l[1]-0.6123724356957945*rhor[0]*uvar2r[0]-0.6123724356957945*rhol[0]*uvar2l[0]; 

  outrrhouz[0] += incr[0]*dxr1; 
  outrrhouz[1] += incr[1]*dxr1; 
  outrrhouz[2] += incr[2]*dxr1; 
  outrrhouz[3] += incr[3]*dxr1; 
  outrrhouz[4] += incr[4]*dxr1; 
  outrrhouz[5] += incr[5]*dxr1; 
  outrrhouz[6] += incr[6]*dxr1; 
  outrrhouz[7] += incr[7]*dxr1; 

  outlrhouz[0] += -1.0*incr[0]*dxl1; 
  outlrhouz[3] += incr[3]*dxl1; 

 
} 
