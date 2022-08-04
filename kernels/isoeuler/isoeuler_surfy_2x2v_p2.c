#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_2x2v_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[8]; 
  const double *rhou1l = &statevecl[16]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[8]; 
  double *outlrhouy = &out[16]; 
  double *outrrho = &out[24]; 
  double *outrrhoux = &out[32]; 
  double *outrrhouy = &out[40]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = 0.5590169943749475*rhou0r[4]+0.5590169943749475*rhou0l[4]-0.4330127018922193*rhou0r[1]+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = (-0.9682458365518543*rhou0r[4])-0.9682458365518543*rhou0l[4]+0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 
  incr[2] = 0.5590169943749476*rhou0r[6]+0.5590169943749476*rhou0l[6]-0.4330127018922193*rhou0r[3]+0.4330127018922193*rhou0l[3]+0.25*rhou0r[2]+0.25*rhou0l[2]; 
  incr[3] = (-0.9682458365518543*rhou0r[6])-0.9682458365518543*rhou0l[6]+0.75*rhou0r[3]-0.75*rhou0l[3]-0.4330127018922194*rhou0r[2]-0.4330127018922194*rhou0l[2]; 
  incr[4] = 1.25*rhou0r[4]+1.25*rhou0l[4]-0.9682458365518543*rhou0r[1]+0.9682458365518543*rhou0l[1]+0.5590169943749475*rhou0r[0]+0.5590169943749475*rhou0l[0]; 
  incr[5] = (-0.4330127018922193*rhou0r[7])+0.4330127018922193*rhou0l[7]+0.25*rhou0r[5]+0.25*rhou0l[5]; 
  incr[6] = 1.25*rhou0r[6]+1.25*rhou0l[6]-0.9682458365518543*rhou0r[3]+0.9682458365518543*rhou0l[3]+0.5590169943749476*rhou0r[2]+0.5590169943749476*rhou0l[2]; 
  incr[7] = 0.75*rhou0r[7]-0.75*rhou0l[7]-0.4330127018922194*rhou0r[5]-0.4330127018922194*rhou0l[5]; 

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
  outlrho[3] += incr[3]*dxl1; 
  outlrho[4] += -1.0*incr[4]*dxl1; 
  outlrho[5] += -1.0*incr[5]*dxl1; 
  outlrho[6] += -1.0*incr[6]*dxl1; 
  outlrho[7] += incr[7]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = (-1.677050983124842*rhor[4]*vthsq)+1.677050983124842*rhol[4]*vthsq+0.4330127018922193*rhor[1]*vthsq+0.4330127018922193*rhol[1]*vthsq+0.75*rhor[7]*uvar0r[7]-0.4330127018922194*rhor[5]*uvar0r[7]+0.75*rhol[7]*uvar0l[7]+0.4330127018922194*rhol[5]*uvar0l[7]-0.4330127018922194*uvar0r[5]*rhor[7]+0.4330127018922194*uvar0l[5]*rhol[7]+1.25*rhor[6]*uvar0r[6]-0.9682458365518543*rhor[3]*uvar0r[6]+0.5590169943749472*rhor[2]*uvar0r[6]+1.25*rhol[6]*uvar0l[6]+0.9682458365518543*rhol[3]*uvar0l[6]+0.5590169943749472*rhol[2]*uvar0l[6]-0.9682458365518543*uvar0r[3]*rhor[6]+0.5590169943749472*uvar0r[2]*rhor[6]+0.9682458365518543*uvar0l[3]*rhol[6]+0.5590169943749472*uvar0l[2]*rhol[6]+0.25*rhor[5]*uvar0r[5]+0.25*rhol[5]*uvar0l[5]+1.25*rhor[4]*uvar0r[4]-0.9682458365518541*rhor[1]*uvar0r[4]+0.5590169943749475*rhor[0]*uvar0r[4]+1.25*rhol[4]*uvar0l[4]+0.9682458365518541*rhol[1]*uvar0l[4]+0.5590169943749475*rhol[0]*uvar0l[4]-0.9682458365518541*uvar0r[1]*rhor[4]+0.5590169943749475*uvar0r[0]*rhor[4]+0.9682458365518541*uvar0l[1]*rhol[4]+0.5590169943749475*uvar0l[0]*rhol[4]+0.75*rhor[3]*uvar0r[3]-0.4330127018922193*rhor[2]*uvar0r[3]+0.75*rhol[3]*uvar0l[3]+0.4330127018922193*rhol[2]*uvar0l[3]-0.4330127018922193*uvar0r[2]*rhor[3]+0.4330127018922193*uvar0l[2]*rhol[3]+0.25*rhor[2]*uvar0r[2]+0.25*rhol[2]*uvar0l[2]+0.75*rhor[1]*uvar0r[1]-0.4330127018922193*rhor[0]*uvar0r[1]+0.75*rhol[1]*uvar0l[1]+0.4330127018922193*rhol[0]*uvar0l[1]-0.4330127018922193*uvar0r[0]*rhor[1]+0.4330127018922193*uvar0l[0]*rhol[1]+0.25*rhor[0]*uvar0r[0]+0.25*rhol[0]*uvar0l[0]; 
  incr[1] = (-5.031152949374527*rhor[4]*vthsq)-5.031152949374527*rhol[4]*vthsq+1.299038105676658*rhor[1]*vthsq-1.299038105676658*rhol[1]*vthsq-1.299038105676658*rhor[7]*uvar0r[7]+0.75*rhor[5]*uvar0r[7]-1.299038105676658*rhol[7]*uvar0l[7]-0.75*rhol[5]*uvar0l[7]+0.75*uvar0r[5]*rhor[7]-0.75*uvar0l[5]*rhol[7]-2.165063509461096*rhor[6]*uvar0r[6]+1.677050983124842*rhor[3]*uvar0r[6]-0.9682458365518543*rhor[2]*uvar0r[6]-2.165063509461096*rhol[6]*uvar0l[6]-1.677050983124842*rhol[3]*uvar0l[6]-0.9682458365518543*rhol[2]*uvar0l[6]+1.677050983124842*uvar0r[3]*rhor[6]-0.9682458365518543*uvar0r[2]*rhor[6]-1.677050983124842*uvar0l[3]*rhol[6]-0.9682458365518543*uvar0l[2]*rhol[6]-0.4330127018922193*rhor[5]*uvar0r[5]-0.4330127018922193*rhol[5]*uvar0l[5]-2.165063509461096*rhor[4]*uvar0r[4]+1.677050983124842*rhor[1]*uvar0r[4]-0.9682458365518543*rhor[0]*uvar0r[4]-2.165063509461096*rhol[4]*uvar0l[4]-1.677050983124842*rhol[1]*uvar0l[4]-0.9682458365518543*rhol[0]*uvar0l[4]+1.677050983124842*uvar0r[1]*rhor[4]-0.9682458365518543*uvar0r[0]*rhor[4]-1.677050983124842*uvar0l[1]*rhol[4]-0.9682458365518543*uvar0l[0]*rhol[4]-1.299038105676658*rhor[3]*uvar0r[3]+0.75*rhor[2]*uvar0r[3]-1.299038105676658*rhol[3]*uvar0l[3]-0.75*rhol[2]*uvar0l[3]+0.75*uvar0r[2]*rhor[3]-0.75*uvar0l[2]*rhol[3]-0.4330127018922193*rhor[2]*uvar0r[2]-0.4330127018922193*rhol[2]*uvar0l[2]-1.299038105676658*rhor[1]*uvar0r[1]+0.75*rhor[0]*uvar0r[1]-1.299038105676658*rhol[1]*uvar0l[1]-0.75*rhol[0]*uvar0l[1]+0.75*uvar0r[0]*rhor[1]-0.75*uvar0l[0]*rhol[1]-0.4330127018922193*rhor[0]*uvar0r[0]-0.4330127018922193*rhol[0]*uvar0l[0]; 
  incr[4] = (-8.385254915624213*rhor[4]*vthsq)+8.385254915624213*rhol[4]*vthsq+2.165063509461096*rhor[1]*vthsq+2.165063509461096*rhol[1]*vthsq+1.677050983124842*rhor[7]*uvar0r[7]-0.9682458365518543*rhor[5]*uvar0r[7]+1.677050983124842*rhol[7]*uvar0l[7]+0.9682458365518543*rhol[5]*uvar0l[7]-0.9682458365518543*uvar0r[5]*rhor[7]+0.9682458365518543*uvar0l[5]*rhol[7]+2.795084971874738*rhor[6]*uvar0r[6]-2.165063509461097*rhor[3]*uvar0r[6]+1.25*rhor[2]*uvar0r[6]+2.795084971874738*rhol[6]*uvar0l[6]+2.165063509461097*rhol[3]*uvar0l[6]+1.25*rhol[2]*uvar0l[6]-2.165063509461097*uvar0r[3]*rhor[6]+1.25*uvar0r[2]*rhor[6]+2.165063509461097*uvar0l[3]*rhol[6]+1.25*uvar0l[2]*rhol[6]+0.5590169943749475*rhor[5]*uvar0r[5]+0.5590169943749475*rhol[5]*uvar0l[5]+2.795084971874738*rhor[4]*uvar0r[4]-2.165063509461096*rhor[1]*uvar0r[4]+1.25*rhor[0]*uvar0r[4]+2.795084971874738*rhol[4]*uvar0l[4]+2.165063509461096*rhol[1]*uvar0l[4]+1.25*rhol[0]*uvar0l[4]-2.165063509461096*uvar0r[1]*rhor[4]+1.25*uvar0r[0]*rhor[4]+2.165063509461096*uvar0l[1]*rhol[4]+1.25*uvar0l[0]*rhol[4]+1.677050983124842*rhor[3]*uvar0r[3]-0.9682458365518543*rhor[2]*uvar0r[3]+1.677050983124842*rhol[3]*uvar0l[3]+0.9682458365518543*rhol[2]*uvar0l[3]-0.9682458365518543*uvar0r[2]*rhor[3]+0.9682458365518543*uvar0l[2]*rhol[3]+0.5590169943749475*rhor[2]*uvar0r[2]+0.5590169943749475*rhol[2]*uvar0l[2]+1.677050983124842*rhor[1]*uvar0r[1]-0.9682458365518543*rhor[0]*uvar0r[1]+1.677050983124842*rhol[1]*uvar0l[1]+0.9682458365518543*rhol[0]*uvar0l[1]-0.9682458365518543*uvar0r[0]*rhor[1]+0.9682458365518543*uvar0l[0]*rhol[1]+0.5590169943749475*rhor[0]*uvar0r[0]+0.5590169943749475*rhol[0]*uvar0l[0]; 

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
  outlrhoux[4] += -1.0*incr[4]*dxl1; 

 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 


 
} 
GKYL_CU_DH void isoeuler_surfy_2x2v_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[1]; 
  const double dxr1 = 2.0/dxv[1]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[8]; 
  const double *rhou1l = &statevecl[16]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[8]; 
  const double *rhou1r = &statevecr[16]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[8]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[8]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[8]; 
  double *outlrhouy = &out[16]; 
  double *outrrho = &out[24]; 
  double *outrrhoux = &out[32]; 
  double *outrrhouy = &out[40]; 
  double incr[8]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = 0.5590169943749475*rhou1r[5]+0.5590169943749475*rhou1l[5]-0.4330127018922193*rhou1r[2]+0.4330127018922193*rhou1l[2]+0.25*rhou1r[0]+0.25*rhou1l[0]; 
  incr[1] = 0.5590169943749476*rhou1r[7]+0.5590169943749476*rhou1l[7]-0.4330127018922193*rhou1r[3]+0.4330127018922193*rhou1l[3]+0.25*rhou1r[1]+0.25*rhou1l[1]; 
  incr[2] = (-0.9682458365518543*rhou1r[5])-0.9682458365518543*rhou1l[5]+0.75*rhou1r[2]-0.75*rhou1l[2]-0.4330127018922193*rhou1r[0]-0.4330127018922193*rhou1l[0]; 
  incr[3] = (-0.9682458365518543*rhou1r[7])-0.9682458365518543*rhou1l[7]+0.75*rhou1r[3]-0.75*rhou1l[3]-0.4330127018922194*rhou1r[1]-0.4330127018922194*rhou1l[1]; 
  incr[4] = (-0.4330127018922193*rhou1r[6])+0.4330127018922193*rhou1l[6]+0.25*rhou1r[4]+0.25*rhou1l[4]; 
  incr[5] = 1.25*rhou1r[5]+1.25*rhou1l[5]-0.9682458365518543*rhou1r[2]+0.9682458365518543*rhou1l[2]+0.5590169943749475*rhou1r[0]+0.5590169943749475*rhou1l[0]; 
  incr[6] = 0.75*rhou1r[6]-0.75*rhou1l[6]-0.4330127018922194*rhou1r[4]-0.4330127018922194*rhou1l[4]; 
  incr[7] = 1.25*rhou1r[7]+1.25*rhou1l[7]-0.9682458365518543*rhou1r[3]+0.9682458365518543*rhou1l[3]+0.5590169943749476*rhou1r[1]+0.5590169943749476*rhou1l[1]; 

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
  outlrho[3] += incr[3]*dxl1; 
  outlrho[4] += -1.0*incr[4]*dxl1; 
  outlrho[5] += -1.0*incr[5]*dxl1; 
  outlrho[6] += incr[6]*dxl1; 
  outlrho[7] += -1.0*incr[7]*dxl1; 

 
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
  incr[0] = (-1.677050983124842*rhor[5]*vthsq)+1.677050983124842*rhol[5]*vthsq+0.4330127018922193*rhor[2]*vthsq+0.4330127018922193*rhol[2]*vthsq+1.25*rhor[7]*uvar1r[7]-0.9682458365518543*rhor[3]*uvar1r[7]+0.5590169943749472*rhor[1]*uvar1r[7]+1.25*rhol[7]*uvar1l[7]+0.9682458365518543*rhol[3]*uvar1l[7]+0.5590169943749472*rhol[1]*uvar1l[7]-0.9682458365518543*uvar1r[3]*rhor[7]+0.5590169943749472*uvar1r[1]*rhor[7]+0.9682458365518543*uvar1l[3]*rhol[7]+0.5590169943749472*uvar1l[1]*rhol[7]+0.75*rhor[6]*uvar1r[6]-0.4330127018922194*rhor[4]*uvar1r[6]+0.75*rhol[6]*uvar1l[6]+0.4330127018922194*rhol[4]*uvar1l[6]-0.4330127018922194*uvar1r[4]*rhor[6]+0.4330127018922194*uvar1l[4]*rhol[6]+1.25*rhor[5]*uvar1r[5]-0.9682458365518541*rhor[2]*uvar1r[5]+0.5590169943749475*rhor[0]*uvar1r[5]+1.25*rhol[5]*uvar1l[5]+0.9682458365518541*rhol[2]*uvar1l[5]+0.5590169943749475*rhol[0]*uvar1l[5]-0.9682458365518541*uvar1r[2]*rhor[5]+0.5590169943749475*uvar1r[0]*rhor[5]+0.9682458365518541*uvar1l[2]*rhol[5]+0.5590169943749475*uvar1l[0]*rhol[5]+0.25*rhor[4]*uvar1r[4]+0.25*rhol[4]*uvar1l[4]+0.75*rhor[3]*uvar1r[3]-0.4330127018922193*rhor[1]*uvar1r[3]+0.75*rhol[3]*uvar1l[3]+0.4330127018922193*rhol[1]*uvar1l[3]-0.4330127018922193*uvar1r[1]*rhor[3]+0.4330127018922193*uvar1l[1]*rhol[3]+0.75*rhor[2]*uvar1r[2]-0.4330127018922193*rhor[0]*uvar1r[2]+0.75*rhol[2]*uvar1l[2]+0.4330127018922193*rhol[0]*uvar1l[2]-0.4330127018922193*uvar1r[0]*rhor[2]+0.4330127018922193*uvar1l[0]*rhol[2]+0.25*rhor[1]*uvar1r[1]+0.25*rhol[1]*uvar1l[1]+0.25*rhor[0]*uvar1r[0]+0.25*rhol[0]*uvar1l[0]; 
  incr[2] = (-5.031152949374527*rhor[5]*vthsq)-5.031152949374527*rhol[5]*vthsq+1.299038105676658*rhor[2]*vthsq-1.299038105676658*rhol[2]*vthsq-2.165063509461096*rhor[7]*uvar1r[7]+1.677050983124842*rhor[3]*uvar1r[7]-0.9682458365518543*rhor[1]*uvar1r[7]-2.165063509461096*rhol[7]*uvar1l[7]-1.677050983124842*rhol[3]*uvar1l[7]-0.9682458365518543*rhol[1]*uvar1l[7]+1.677050983124842*uvar1r[3]*rhor[7]-0.9682458365518543*uvar1r[1]*rhor[7]-1.677050983124842*uvar1l[3]*rhol[7]-0.9682458365518543*uvar1l[1]*rhol[7]-1.299038105676658*rhor[6]*uvar1r[6]+0.75*rhor[4]*uvar1r[6]-1.299038105676658*rhol[6]*uvar1l[6]-0.75*rhol[4]*uvar1l[6]+0.75*uvar1r[4]*rhor[6]-0.75*uvar1l[4]*rhol[6]-2.165063509461096*rhor[5]*uvar1r[5]+1.677050983124842*rhor[2]*uvar1r[5]-0.9682458365518543*rhor[0]*uvar1r[5]-2.165063509461096*rhol[5]*uvar1l[5]-1.677050983124842*rhol[2]*uvar1l[5]-0.9682458365518543*rhol[0]*uvar1l[5]+1.677050983124842*uvar1r[2]*rhor[5]-0.9682458365518543*uvar1r[0]*rhor[5]-1.677050983124842*uvar1l[2]*rhol[5]-0.9682458365518543*uvar1l[0]*rhol[5]-0.4330127018922193*rhor[4]*uvar1r[4]-0.4330127018922193*rhol[4]*uvar1l[4]-1.299038105676658*rhor[3]*uvar1r[3]+0.75*rhor[1]*uvar1r[3]-1.299038105676658*rhol[3]*uvar1l[3]-0.75*rhol[1]*uvar1l[3]+0.75*uvar1r[1]*rhor[3]-0.75*uvar1l[1]*rhol[3]-1.299038105676658*rhor[2]*uvar1r[2]+0.75*rhor[0]*uvar1r[2]-1.299038105676658*rhol[2]*uvar1l[2]-0.75*rhol[0]*uvar1l[2]+0.75*uvar1r[0]*rhor[2]-0.75*uvar1l[0]*rhol[2]-0.4330127018922193*rhor[1]*uvar1r[1]-0.4330127018922193*rhol[1]*uvar1l[1]-0.4330127018922193*rhor[0]*uvar1r[0]-0.4330127018922193*rhol[0]*uvar1l[0]; 
  incr[5] = (-8.385254915624213*rhor[5]*vthsq)+8.385254915624213*rhol[5]*vthsq+2.165063509461096*rhor[2]*vthsq+2.165063509461096*rhol[2]*vthsq+2.795084971874738*rhor[7]*uvar1r[7]-2.165063509461097*rhor[3]*uvar1r[7]+1.25*rhor[1]*uvar1r[7]+2.795084971874738*rhol[7]*uvar1l[7]+2.165063509461097*rhol[3]*uvar1l[7]+1.25*rhol[1]*uvar1l[7]-2.165063509461097*uvar1r[3]*rhor[7]+1.25*uvar1r[1]*rhor[7]+2.165063509461097*uvar1l[3]*rhol[7]+1.25*uvar1l[1]*rhol[7]+1.677050983124842*rhor[6]*uvar1r[6]-0.9682458365518543*rhor[4]*uvar1r[6]+1.677050983124842*rhol[6]*uvar1l[6]+0.9682458365518543*rhol[4]*uvar1l[6]-0.9682458365518543*uvar1r[4]*rhor[6]+0.9682458365518543*uvar1l[4]*rhol[6]+2.795084971874738*rhor[5]*uvar1r[5]-2.165063509461096*rhor[2]*uvar1r[5]+1.25*rhor[0]*uvar1r[5]+2.795084971874738*rhol[5]*uvar1l[5]+2.165063509461096*rhol[2]*uvar1l[5]+1.25*rhol[0]*uvar1l[5]-2.165063509461096*uvar1r[2]*rhor[5]+1.25*uvar1r[0]*rhor[5]+2.165063509461096*uvar1l[2]*rhol[5]+1.25*uvar1l[0]*rhol[5]+0.5590169943749475*rhor[4]*uvar1r[4]+0.5590169943749475*rhol[4]*uvar1l[4]+1.677050983124842*rhor[3]*uvar1r[3]-0.9682458365518543*rhor[1]*uvar1r[3]+1.677050983124842*rhol[3]*uvar1l[3]+0.9682458365518543*rhol[1]*uvar1l[3]-0.9682458365518543*uvar1r[1]*rhor[3]+0.9682458365518543*uvar1l[1]*rhol[3]+1.677050983124842*rhor[2]*uvar1r[2]-0.9682458365518543*rhor[0]*uvar1r[2]+1.677050983124842*rhol[2]*uvar1l[2]+0.9682458365518543*rhol[0]*uvar1l[2]-0.9682458365518543*uvar1r[0]*rhor[2]+0.9682458365518543*uvar1l[0]*rhol[2]+0.5590169943749475*rhor[1]*uvar1r[1]+0.5590169943749475*rhol[1]*uvar1l[1]+0.5590169943749475*rhor[0]*uvar1r[0]+0.5590169943749475*rhol[0]*uvar1l[0]; 

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
  outlrhouy[5] += -1.0*incr[5]*dxl1; 

 
} 
