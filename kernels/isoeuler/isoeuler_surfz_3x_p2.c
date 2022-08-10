#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_3x_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[20]; 
  const double *rhou1l = &statevecl[40]; 
  const double *rhou2l = &statevecl[60]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[20]; 
  const double *rhou1r = &statevecr[40]; 
  const double *rhou2r = &statevecr[60]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[20]; 
  const double *uvar2l = &uvarl[40]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[20]; 
  const double *uvar2r = &uvarr[40]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[20]; 
  double *outlrhouy = &out[40]; 
  double *outlrhouz = &out[60]; 
  double *outrrho = &out[80]; 
  double *outrrhoux = &out[100]; 
  double *outrrhouy = &out[120]; 
  double *outrrhouz = &out[140]; 
  double incr[20]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = 0.5590169943749475*rhou0r[7]+0.5590169943749475*rhou0l[7]-0.4330127018922193*rhou0r[1]+0.4330127018922193*rhou0l[1]+0.25*rhou0r[0]+0.25*rhou0l[0]; 
  incr[1] = (-0.9682458365518543*rhou0r[7])-0.9682458365518543*rhou0l[7]+0.75*rhou0r[1]-0.75*rhou0l[1]-0.4330127018922193*rhou0r[0]-0.4330127018922193*rhou0l[0]; 
  incr[2] = 0.5590169943749476*rhou0r[11]+0.5590169943749476*rhou0l[11]-0.4330127018922193*rhou0r[4]+0.4330127018922193*rhou0l[4]+0.25*rhou0r[2]+0.25*rhou0l[2]; 
  incr[3] = 0.5590169943749476*rhou0r[13]+0.5590169943749476*rhou0l[13]-0.4330127018922193*rhou0r[5]+0.4330127018922193*rhou0l[5]+0.25*rhou0r[3]+0.25*rhou0l[3]; 
  incr[4] = (-0.9682458365518543*rhou0r[11])-0.9682458365518543*rhou0l[11]+0.75*rhou0r[4]-0.75*rhou0l[4]-0.4330127018922194*rhou0r[2]-0.4330127018922194*rhou0l[2]; 
  incr[5] = (-0.9682458365518543*rhou0r[13])-0.9682458365518543*rhou0l[13]+0.75*rhou0r[5]-0.75*rhou0l[5]-0.4330127018922194*rhou0r[3]-0.4330127018922194*rhou0l[3]; 
  incr[6] = 0.5590169943749475*rhou0r[17]+0.5590169943749475*rhou0l[17]-0.4330127018922194*rhou0r[10]+0.4330127018922194*rhou0l[10]+0.25*rhou0r[6]+0.25*rhou0l[6]; 
  incr[7] = 1.25*rhou0r[7]+1.25*rhou0l[7]-0.9682458365518543*rhou0r[1]+0.9682458365518543*rhou0l[1]+0.5590169943749475*rhou0r[0]+0.5590169943749475*rhou0l[0]; 
  incr[8] = (-0.4330127018922193*rhou0r[12])+0.4330127018922193*rhou0l[12]+0.25*rhou0r[8]+0.25*rhou0l[8]; 
  incr[9] = (-0.4330127018922193*rhou0r[15])+0.4330127018922193*rhou0l[15]+0.25*rhou0r[9]+0.25*rhou0l[9]; 
  incr[10] = (-0.9682458365518543*rhou0r[17])-0.9682458365518543*rhou0l[17]+0.75*rhou0r[10]-0.75*rhou0l[10]-0.4330127018922193*rhou0r[6]-0.4330127018922193*rhou0l[6]; 
  incr[11] = 1.25*rhou0r[11]+1.25*rhou0l[11]-0.9682458365518543*rhou0r[4]+0.9682458365518543*rhou0l[4]+0.5590169943749476*rhou0r[2]+0.5590169943749476*rhou0l[2]; 
  incr[12] = 0.75*rhou0r[12]-0.75*rhou0l[12]-0.4330127018922194*rhou0r[8]-0.4330127018922194*rhou0l[8]; 
  incr[13] = 1.25*rhou0r[13]+1.25*rhou0l[13]-0.9682458365518543*rhou0r[5]+0.9682458365518543*rhou0l[5]+0.5590169943749476*rhou0r[3]+0.5590169943749476*rhou0l[3]; 
  incr[14] = (-0.4330127018922194*rhou0r[18])+0.4330127018922194*rhou0l[18]+0.25*rhou0r[14]+0.25*rhou0l[14]; 
  incr[15] = 0.75*rhou0r[15]-0.75*rhou0l[15]-0.4330127018922194*rhou0r[9]-0.4330127018922194*rhou0l[9]; 
  incr[16] = (-0.4330127018922194*rhou0r[19])+0.4330127018922194*rhou0l[19]+0.25*rhou0r[16]+0.25*rhou0l[16]; 
  incr[17] = 1.25*rhou0r[17]+1.25*rhou0l[17]-0.9682458365518545*rhou0r[10]+0.9682458365518545*rhou0l[10]+0.5590169943749475*rhou0r[6]+0.5590169943749475*rhou0l[6]; 
  incr[18] = 0.75*rhou0r[18]-0.75*rhou0l[18]-0.4330127018922193*rhou0r[14]-0.4330127018922193*rhou0l[14]; 
  incr[19] = 0.75*rhou0r[19]-0.75*rhou0l[19]-0.4330127018922193*rhou0r[16]-0.4330127018922193*rhou0l[16]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 
  outrrho[8] += incr[8]*dxr1; 
  outrrho[9] += incr[9]*dxr1; 
  outrrho[10] += incr[10]*dxr1; 
  outrrho[11] += incr[11]*dxr1; 
  outrrho[12] += incr[12]*dxr1; 
  outrrho[13] += incr[13]*dxr1; 
  outrrho[14] += incr[14]*dxr1; 
  outrrho[15] += incr[15]*dxr1; 
  outrrho[16] += incr[16]*dxr1; 
  outrrho[17] += incr[17]*dxr1; 
  outrrho[18] += incr[18]*dxr1; 
  outrrho[19] += incr[19]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += -1.0*incr[3]*dxl1; 
  outlrho[4] += incr[4]*dxl1; 
  outlrho[5] += incr[5]*dxl1; 
  outlrho[6] += -1.0*incr[6]*dxl1; 
  outlrho[7] += -1.0*incr[7]*dxl1; 
  outlrho[8] += -1.0*incr[8]*dxl1; 
  outlrho[9] += -1.0*incr[9]*dxl1; 
  outlrho[10] += incr[10]*dxl1; 
  outlrho[11] += -1.0*incr[11]*dxl1; 
  outlrho[12] += incr[12]*dxl1; 
  outlrho[13] += -1.0*incr[13]*dxl1; 
  outlrho[14] += -1.0*incr[14]*dxl1; 
  outlrho[15] += incr[15]*dxl1; 
  outlrho[16] += -1.0*incr[16]*dxl1; 
  outlrho[17] += -1.0*incr[17]*dxl1; 
  outlrho[18] += incr[18]*dxl1; 
  outlrho[19] += incr[19]*dxl1; 

 
  //FluxRhoUx; 
  incr[0] = (-2.371708245126285*rhor[7]*vthsq)+2.371708245126285*rhol[7]*vthsq+0.6123724356957945*rhor[1]*vthsq+0.6123724356957945*rhol[1]*vthsq+1.060660171779821*rhor[19]*uvar0r[19]-0.6123724356957945*rhor[16]*uvar0r[19]+1.060660171779821*rhol[19]*uvar0l[19]+0.6123724356957945*rhol[16]*uvar0l[19]-0.6123724356957945*uvar0r[16]*rhor[19]+0.6123724356957945*uvar0l[16]*rhol[19]+1.060660171779821*rhor[18]*uvar0r[18]-0.6123724356957945*rhor[14]*uvar0r[18]+1.060660171779821*rhol[18]*uvar0l[18]+0.6123724356957945*rhol[14]*uvar0l[18]-0.6123724356957945*uvar0r[14]*rhor[18]+0.6123724356957945*uvar0l[14]*rhol[18]+1.767766952966369*rhor[17]*uvar0r[17]-1.369306393762915*rhor[10]*uvar0r[17]+0.7905694150420948*rhor[6]*uvar0r[17]+1.767766952966369*rhol[17]*uvar0l[17]+1.369306393762915*rhol[10]*uvar0l[17]+0.7905694150420948*rhol[6]*uvar0l[17]-1.369306393762915*uvar0r[10]*rhor[17]+0.7905694150420948*uvar0r[6]*rhor[17]+1.369306393762915*uvar0l[10]*rhol[17]+0.7905694150420948*uvar0l[6]*rhol[17]+0.3535533905932737*rhor[16]*uvar0r[16]+0.3535533905932737*rhol[16]*uvar0l[16]+1.060660171779821*rhor[15]*uvar0r[15]-0.6123724356957945*rhor[9]*uvar0r[15]+1.060660171779821*rhol[15]*uvar0l[15]+0.6123724356957945*rhol[9]*uvar0l[15]-0.6123724356957945*uvar0r[9]*rhor[15]+0.6123724356957945*uvar0l[9]*rhol[15]+0.3535533905932737*rhor[14]*uvar0r[14]+0.3535533905932737*rhol[14]*uvar0l[14]+1.767766952966369*rhor[13]*uvar0r[13]-1.369306393762915*rhor[5]*uvar0r[13]+0.7905694150420947*rhor[3]*uvar0r[13]+1.767766952966369*rhol[13]*uvar0l[13]+1.369306393762915*rhol[5]*uvar0l[13]+0.7905694150420947*rhol[3]*uvar0l[13]-1.369306393762915*uvar0r[5]*rhor[13]+0.7905694150420947*uvar0r[3]*rhor[13]+1.369306393762915*uvar0l[5]*rhol[13]+0.7905694150420947*uvar0l[3]*rhol[13]+1.060660171779821*rhor[12]*uvar0r[12]-0.6123724356957945*rhor[8]*uvar0r[12]+1.060660171779821*rhol[12]*uvar0l[12]+0.6123724356957945*rhol[8]*uvar0l[12]-0.6123724356957945*uvar0r[8]*rhor[12]+0.6123724356957945*uvar0l[8]*rhol[12]+1.767766952966369*rhor[11]*uvar0r[11]-1.369306393762915*rhor[4]*uvar0r[11]+0.7905694150420947*rhor[2]*uvar0r[11]+1.767766952966369*rhol[11]*uvar0l[11]+1.369306393762915*rhol[4]*uvar0l[11]+0.7905694150420947*rhol[2]*uvar0l[11]-1.369306393762915*uvar0r[4]*rhor[11]+0.7905694150420947*uvar0r[2]*rhor[11]+1.369306393762915*uvar0l[4]*rhol[11]+0.7905694150420947*uvar0l[2]*rhol[11]+1.060660171779821*rhor[10]*uvar0r[10]-0.6123724356957945*rhor[6]*uvar0r[10]+1.060660171779821*rhol[10]*uvar0l[10]+0.6123724356957945*rhol[6]*uvar0l[10]-0.6123724356957945*uvar0r[6]*rhor[10]+0.6123724356957945*uvar0l[6]*rhol[10]+0.3535533905932737*rhor[9]*uvar0r[9]+0.3535533905932737*rhol[9]*uvar0l[9]+0.3535533905932737*rhor[8]*uvar0r[8]+0.3535533905932737*rhol[8]*uvar0l[8]+1.767766952966369*rhor[7]*uvar0r[7]-1.369306393762915*rhor[1]*uvar0r[7]+0.7905694150420948*rhor[0]*uvar0r[7]+1.767766952966369*rhol[7]*uvar0l[7]+1.369306393762915*rhol[1]*uvar0l[7]+0.7905694150420948*rhol[0]*uvar0l[7]-1.369306393762915*uvar0r[1]*rhor[7]+0.7905694150420948*uvar0r[0]*rhor[7]+1.369306393762915*uvar0l[1]*rhol[7]+0.7905694150420948*uvar0l[0]*rhol[7]+0.3535533905932737*rhor[6]*uvar0r[6]+0.3535533905932737*rhol[6]*uvar0l[6]+1.060660171779821*rhor[5]*uvar0r[5]-0.6123724356957945*rhor[3]*uvar0r[5]+1.060660171779821*rhol[5]*uvar0l[5]+0.6123724356957945*rhol[3]*uvar0l[5]-0.6123724356957945*uvar0r[3]*rhor[5]+0.6123724356957945*uvar0l[3]*rhol[5]+1.060660171779821*rhor[4]*uvar0r[4]-0.6123724356957945*rhor[2]*uvar0r[4]+1.060660171779821*rhol[4]*uvar0l[4]+0.6123724356957945*rhol[2]*uvar0l[4]-0.6123724356957945*uvar0r[2]*rhor[4]+0.6123724356957945*uvar0l[2]*rhol[4]+0.3535533905932737*rhor[3]*uvar0r[3]+0.3535533905932737*rhol[3]*uvar0l[3]+0.3535533905932737*rhor[2]*uvar0r[2]+0.3535533905932737*rhol[2]*uvar0l[2]+1.060660171779821*rhor[1]*uvar0r[1]-0.6123724356957945*rhor[0]*uvar0r[1]+1.060660171779821*rhol[1]*uvar0l[1]+0.6123724356957945*rhol[0]*uvar0l[1]-0.6123724356957945*uvar0r[0]*rhor[1]+0.6123724356957945*uvar0l[0]*rhol[1]+0.3535533905932737*rhor[0]*uvar0r[0]+0.3535533905932737*rhol[0]*uvar0l[0]; 
  incr[1] = (-7.115124735378854*rhor[7]*vthsq)-7.115124735378854*rhol[7]*vthsq+1.837117307087383*rhor[1]*vthsq-1.837117307087383*rhol[1]*vthsq-1.837117307087383*rhor[19]*uvar0r[19]+1.060660171779821*rhor[16]*uvar0r[19]-1.837117307087383*rhol[19]*uvar0l[19]-1.060660171779821*rhol[16]*uvar0l[19]+1.060660171779821*uvar0r[16]*rhor[19]-1.060660171779821*uvar0l[16]*rhol[19]-1.837117307087383*rhor[18]*uvar0r[18]+1.060660171779821*rhor[14]*uvar0r[18]-1.837117307087383*rhol[18]*uvar0l[18]-1.060660171779821*rhol[14]*uvar0l[18]+1.060660171779821*uvar0r[14]*rhor[18]-1.060660171779821*uvar0l[14]*rhol[18]-3.061862178478972*rhor[17]*uvar0r[17]+2.371708245126284*rhor[10]*uvar0r[17]-1.369306393762915*rhor[6]*uvar0r[17]-3.061862178478972*rhol[17]*uvar0l[17]-2.371708245126284*rhol[10]*uvar0l[17]-1.369306393762915*rhol[6]*uvar0l[17]+2.371708245126284*uvar0r[10]*rhor[17]-1.369306393762915*uvar0r[6]*rhor[17]-2.371708245126284*uvar0l[10]*rhol[17]-1.369306393762915*uvar0l[6]*rhol[17]-0.6123724356957945*rhor[16]*uvar0r[16]-0.6123724356957945*rhol[16]*uvar0l[16]-1.837117307087383*rhor[15]*uvar0r[15]+1.060660171779821*rhor[9]*uvar0r[15]-1.837117307087383*rhol[15]*uvar0l[15]-1.060660171779821*rhol[9]*uvar0l[15]+1.060660171779821*uvar0r[9]*rhor[15]-1.060660171779821*uvar0l[9]*rhol[15]-0.6123724356957945*rhor[14]*uvar0r[14]-0.6123724356957945*rhol[14]*uvar0l[14]-3.061862178478972*rhor[13]*uvar0r[13]+2.371708245126285*rhor[5]*uvar0r[13]-1.369306393762915*rhor[3]*uvar0r[13]-3.061862178478972*rhol[13]*uvar0l[13]-2.371708245126285*rhol[5]*uvar0l[13]-1.369306393762915*rhol[3]*uvar0l[13]+2.371708245126285*uvar0r[5]*rhor[13]-1.369306393762915*uvar0r[3]*rhor[13]-2.371708245126285*uvar0l[5]*rhol[13]-1.369306393762915*uvar0l[3]*rhol[13]-1.837117307087383*rhor[12]*uvar0r[12]+1.060660171779821*rhor[8]*uvar0r[12]-1.837117307087383*rhol[12]*uvar0l[12]-1.060660171779821*rhol[8]*uvar0l[12]+1.060660171779821*uvar0r[8]*rhor[12]-1.060660171779821*uvar0l[8]*rhol[12]-3.061862178478972*rhor[11]*uvar0r[11]+2.371708245126285*rhor[4]*uvar0r[11]-1.369306393762915*rhor[2]*uvar0r[11]-3.061862178478972*rhol[11]*uvar0l[11]-2.371708245126285*rhol[4]*uvar0l[11]-1.369306393762915*rhol[2]*uvar0l[11]+2.371708245126285*uvar0r[4]*rhor[11]-1.369306393762915*uvar0r[2]*rhor[11]-2.371708245126285*uvar0l[4]*rhol[11]-1.369306393762915*uvar0l[2]*rhol[11]-1.837117307087383*rhor[10]*uvar0r[10]+1.060660171779821*rhor[6]*uvar0r[10]-1.837117307087383*rhol[10]*uvar0l[10]-1.060660171779821*rhol[6]*uvar0l[10]+1.060660171779821*uvar0r[6]*rhor[10]-1.060660171779821*uvar0l[6]*rhol[10]-0.6123724356957945*rhor[9]*uvar0r[9]-0.6123724356957945*rhol[9]*uvar0l[9]-0.6123724356957945*rhor[8]*uvar0r[8]-0.6123724356957945*rhol[8]*uvar0l[8]-3.061862178478972*rhor[7]*uvar0r[7]+2.371708245126284*rhor[1]*uvar0r[7]-1.369306393762915*rhor[0]*uvar0r[7]-3.061862178478972*rhol[7]*uvar0l[7]-2.371708245126284*rhol[1]*uvar0l[7]-1.369306393762915*rhol[0]*uvar0l[7]+2.371708245126284*uvar0r[1]*rhor[7]-1.369306393762915*uvar0r[0]*rhor[7]-2.371708245126284*uvar0l[1]*rhol[7]-1.369306393762915*uvar0l[0]*rhol[7]-0.6123724356957945*rhor[6]*uvar0r[6]-0.6123724356957945*rhol[6]*uvar0l[6]-1.837117307087383*rhor[5]*uvar0r[5]+1.060660171779821*rhor[3]*uvar0r[5]-1.837117307087383*rhol[5]*uvar0l[5]-1.060660171779821*rhol[3]*uvar0l[5]+1.060660171779821*uvar0r[3]*rhor[5]-1.060660171779821*uvar0l[3]*rhol[5]-1.837117307087383*rhor[4]*uvar0r[4]+1.060660171779821*rhor[2]*uvar0r[4]-1.837117307087383*rhol[4]*uvar0l[4]-1.060660171779821*rhol[2]*uvar0l[4]+1.060660171779821*uvar0r[2]*rhor[4]-1.060660171779821*uvar0l[2]*rhol[4]-0.6123724356957945*rhor[3]*uvar0r[3]-0.6123724356957945*rhol[3]*uvar0l[3]-0.6123724356957945*rhor[2]*uvar0r[2]-0.6123724356957945*rhol[2]*uvar0l[2]-1.837117307087383*rhor[1]*uvar0r[1]+1.060660171779821*rhor[0]*uvar0r[1]-1.837117307087383*rhol[1]*uvar0l[1]-1.060660171779821*rhol[0]*uvar0l[1]+1.060660171779821*uvar0r[0]*rhor[1]-1.060660171779821*uvar0l[0]*rhol[1]-0.6123724356957945*rhor[0]*uvar0r[0]-0.6123724356957945*rhol[0]*uvar0l[0]; 
  incr[7] = (-11.85854122563142*rhor[7]*vthsq)+11.85854122563142*rhol[7]*vthsq+3.061862178478972*rhor[1]*vthsq+3.061862178478972*rhol[1]*vthsq+2.371708245126285*rhor[19]*uvar0r[19]-1.369306393762915*rhor[16]*uvar0r[19]+2.371708245126285*rhol[19]*uvar0l[19]+1.369306393762915*rhol[16]*uvar0l[19]-1.369306393762915*uvar0r[16]*rhor[19]+1.369306393762915*uvar0l[16]*rhol[19]+2.371708245126285*rhor[18]*uvar0r[18]-1.369306393762915*rhor[14]*uvar0r[18]+2.371708245126285*rhol[18]*uvar0l[18]+1.369306393762915*rhol[14]*uvar0l[18]-1.369306393762915*uvar0r[14]*rhor[18]+1.369306393762915*uvar0l[14]*rhol[18]+3.952847075210474*rhor[17]*uvar0r[17]-3.061862178478972*rhor[10]*uvar0r[17]+1.767766952966369*rhor[6]*uvar0r[17]+3.952847075210474*rhol[17]*uvar0l[17]+3.061862178478972*rhol[10]*uvar0l[17]+1.767766952966369*rhol[6]*uvar0l[17]-3.061862178478972*uvar0r[10]*rhor[17]+1.767766952966369*uvar0r[6]*rhor[17]+3.061862178478972*uvar0l[10]*rhol[17]+1.767766952966369*uvar0l[6]*rhol[17]+0.7905694150420948*rhor[16]*uvar0r[16]+0.7905694150420948*rhol[16]*uvar0l[16]+2.371708245126285*rhor[15]*uvar0r[15]-1.369306393762915*rhor[9]*uvar0r[15]+2.371708245126285*rhol[15]*uvar0l[15]+1.369306393762915*rhol[9]*uvar0l[15]-1.369306393762915*uvar0r[9]*rhor[15]+1.369306393762915*uvar0l[9]*rhol[15]+0.7905694150420948*rhor[14]*uvar0r[14]+0.7905694150420948*rhol[14]*uvar0l[14]+3.952847075210474*rhor[13]*uvar0r[13]-3.061862178478973*rhor[5]*uvar0r[13]+1.767766952966368*rhor[3]*uvar0r[13]+3.952847075210474*rhol[13]*uvar0l[13]+3.061862178478973*rhol[5]*uvar0l[13]+1.767766952966368*rhol[3]*uvar0l[13]-3.061862178478973*uvar0r[5]*rhor[13]+1.767766952966368*uvar0r[3]*rhor[13]+3.061862178478973*uvar0l[5]*rhol[13]+1.767766952966368*uvar0l[3]*rhol[13]+2.371708245126285*rhor[12]*uvar0r[12]-1.369306393762915*rhor[8]*uvar0r[12]+2.371708245126285*rhol[12]*uvar0l[12]+1.369306393762915*rhol[8]*uvar0l[12]-1.369306393762915*uvar0r[8]*rhor[12]+1.369306393762915*uvar0l[8]*rhol[12]+3.952847075210474*rhor[11]*uvar0r[11]-3.061862178478973*rhor[4]*uvar0r[11]+1.767766952966368*rhor[2]*uvar0r[11]+3.952847075210474*rhol[11]*uvar0l[11]+3.061862178478973*rhol[4]*uvar0l[11]+1.767766952966368*rhol[2]*uvar0l[11]-3.061862178478973*uvar0r[4]*rhor[11]+1.767766952966368*uvar0r[2]*rhor[11]+3.061862178478973*uvar0l[4]*rhol[11]+1.767766952966368*uvar0l[2]*rhol[11]+2.371708245126285*rhor[10]*uvar0r[10]-1.369306393762915*rhor[6]*uvar0r[10]+2.371708245126285*rhol[10]*uvar0l[10]+1.369306393762915*rhol[6]*uvar0l[10]-1.369306393762915*uvar0r[6]*rhor[10]+1.369306393762915*uvar0l[6]*rhol[10]+0.7905694150420948*rhor[9]*uvar0r[9]+0.7905694150420948*rhol[9]*uvar0l[9]+0.7905694150420948*rhor[8]*uvar0r[8]+0.7905694150420948*rhol[8]*uvar0l[8]+3.952847075210474*rhor[7]*uvar0r[7]-3.061862178478972*rhor[1]*uvar0r[7]+1.767766952966369*rhor[0]*uvar0r[7]+3.952847075210474*rhol[7]*uvar0l[7]+3.061862178478972*rhol[1]*uvar0l[7]+1.767766952966369*rhol[0]*uvar0l[7]-3.061862178478972*uvar0r[1]*rhor[7]+1.767766952966369*uvar0r[0]*rhor[7]+3.061862178478972*uvar0l[1]*rhol[7]+1.767766952966369*uvar0l[0]*rhol[7]+0.7905694150420948*rhor[6]*uvar0r[6]+0.7905694150420948*rhol[6]*uvar0l[6]+2.371708245126285*rhor[5]*uvar0r[5]-1.369306393762915*rhor[3]*uvar0r[5]+2.371708245126285*rhol[5]*uvar0l[5]+1.369306393762915*rhol[3]*uvar0l[5]-1.369306393762915*uvar0r[3]*rhor[5]+1.369306393762915*uvar0l[3]*rhol[5]+2.371708245126285*rhor[4]*uvar0r[4]-1.369306393762915*rhor[2]*uvar0r[4]+2.371708245126285*rhol[4]*uvar0l[4]+1.369306393762915*rhol[2]*uvar0l[4]-1.369306393762915*uvar0r[2]*rhor[4]+1.369306393762915*uvar0l[2]*rhol[4]+0.7905694150420948*rhor[3]*uvar0r[3]+0.7905694150420948*rhol[3]*uvar0l[3]+0.7905694150420948*rhor[2]*uvar0r[2]+0.7905694150420948*rhol[2]*uvar0l[2]+2.371708245126285*rhor[1]*uvar0r[1]-1.369306393762915*rhor[0]*uvar0r[1]+2.371708245126285*rhol[1]*uvar0l[1]+1.369306393762915*rhol[0]*uvar0l[1]-1.369306393762915*uvar0r[0]*rhor[1]+1.369306393762915*uvar0l[0]*rhol[1]+0.7905694150420948*rhor[0]*uvar0r[0]+0.7905694150420948*rhol[0]*uvar0l[0]; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 
  outrrhoux[8] += incr[8]*dxr1; 
  outrrhoux[9] += incr[9]*dxr1; 
  outrrhoux[10] += incr[10]*dxr1; 
  outrrhoux[11] += incr[11]*dxr1; 
  outrrhoux[12] += incr[12]*dxr1; 
  outrrhoux[13] += incr[13]*dxr1; 
  outrrhoux[14] += incr[14]*dxr1; 
  outrrhoux[15] += incr[15]*dxr1; 
  outrrhoux[16] += incr[16]*dxr1; 
  outrrhoux[17] += incr[17]*dxr1; 
  outrrhoux[18] += incr[18]*dxr1; 
  outrrhoux[19] += incr[19]*dxr1; 

  outlrhoux[0] += -1.0*incr[0]*dxl1; 
  outlrhoux[1] += incr[1]*dxl1; 
  outlrhoux[7] += -1.0*incr[7]*dxl1; 

 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 
  outrrhouy[8] += incr[8]*dxr1; 
  outrrhouy[9] += incr[9]*dxr1; 
  outrrhouy[10] += incr[10]*dxr1; 
  outrrhouy[11] += incr[11]*dxr1; 
  outrrhouy[12] += incr[12]*dxr1; 
  outrrhouy[13] += incr[13]*dxr1; 
  outrrhouy[14] += incr[14]*dxr1; 
  outrrhouy[15] += incr[15]*dxr1; 
  outrrhouy[16] += incr[16]*dxr1; 
  outrrhouy[17] += incr[17]*dxr1; 
  outrrhouy[18] += incr[18]*dxr1; 
  outrrhouy[19] += incr[19]*dxr1; 


 
  //FluxRhoUz; 

  outrrhouz[0] += incr[0]*dxr1; 
  outrrhouz[1] += incr[1]*dxr1; 
  outrrhouz[2] += incr[2]*dxr1; 
  outrrhouz[3] += incr[3]*dxr1; 
  outrrhouz[4] += incr[4]*dxr1; 
  outrrhouz[5] += incr[5]*dxr1; 
  outrrhouz[6] += incr[6]*dxr1; 
  outrrhouz[7] += incr[7]*dxr1; 
  outrrhouz[8] += incr[8]*dxr1; 
  outrrhouz[9] += incr[9]*dxr1; 
  outrrhouz[10] += incr[10]*dxr1; 
  outrrhouz[11] += incr[11]*dxr1; 
  outrrhouz[12] += incr[12]*dxr1; 
  outrrhouz[13] += incr[13]*dxr1; 
  outrrhouz[14] += incr[14]*dxr1; 
  outrrhouz[15] += incr[15]*dxr1; 
  outrrhouz[16] += incr[16]*dxr1; 
  outrrhouz[17] += incr[17]*dxr1; 
  outrrhouz[18] += incr[18]*dxr1; 
  outrrhouz[19] += incr[19]*dxr1; 


 
} 
GKYL_CU_DH void isoeuler_surfy_3x_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[1]; 
  const double dxr1 = 2.0/dxv[1]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[20]; 
  const double *rhou1l = &statevecl[40]; 
  const double *rhou2l = &statevecl[60]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[20]; 
  const double *rhou1r = &statevecr[40]; 
  const double *rhou2r = &statevecr[60]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[20]; 
  const double *uvar2l = &uvarl[40]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[20]; 
  const double *uvar2r = &uvarr[40]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[20]; 
  double *outlrhouy = &out[40]; 
  double *outlrhouz = &out[60]; 
  double *outrrho = &out[80]; 
  double *outrrhoux = &out[100]; 
  double *outrrhouy = &out[120]; 
  double *outrrhouz = &out[140]; 
  double incr[20]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = 0.5590169943749475*rhou1r[8]+0.5590169943749475*rhou1l[8]-0.4330127018922193*rhou1r[2]+0.4330127018922193*rhou1l[2]+0.25*rhou1r[0]+0.25*rhou1l[0]; 
  incr[1] = 0.5590169943749476*rhou1r[12]+0.5590169943749476*rhou1l[12]-0.4330127018922193*rhou1r[4]+0.4330127018922193*rhou1l[4]+0.25*rhou1r[1]+0.25*rhou1l[1]; 
  incr[2] = (-0.9682458365518543*rhou1r[8])-0.9682458365518543*rhou1l[8]+0.75*rhou1r[2]-0.75*rhou1l[2]-0.4330127018922193*rhou1r[0]-0.4330127018922193*rhou1l[0]; 
  incr[3] = 0.5590169943749476*rhou1r[14]+0.5590169943749476*rhou1l[14]-0.4330127018922193*rhou1r[6]+0.4330127018922193*rhou1l[6]+0.25*rhou1r[3]+0.25*rhou1l[3]; 
  incr[4] = (-0.9682458365518543*rhou1r[12])-0.9682458365518543*rhou1l[12]+0.75*rhou1r[4]-0.75*rhou1l[4]-0.4330127018922194*rhou1r[1]-0.4330127018922194*rhou1l[1]; 
  incr[5] = 0.5590169943749475*rhou1r[18]+0.5590169943749475*rhou1l[18]-0.4330127018922194*rhou1r[10]+0.4330127018922194*rhou1l[10]+0.25*rhou1r[5]+0.25*rhou1l[5]; 
  incr[6] = (-0.9682458365518543*rhou1r[14])-0.9682458365518543*rhou1l[14]+0.75*rhou1r[6]-0.75*rhou1l[6]-0.4330127018922194*rhou1r[3]-0.4330127018922194*rhou1l[3]; 
  incr[7] = (-0.4330127018922193*rhou1r[11])+0.4330127018922193*rhou1l[11]+0.25*rhou1r[7]+0.25*rhou1l[7]; 
  incr[8] = 1.25*rhou1r[8]+1.25*rhou1l[8]-0.9682458365518543*rhou1r[2]+0.9682458365518543*rhou1l[2]+0.5590169943749475*rhou1r[0]+0.5590169943749475*rhou1l[0]; 
  incr[9] = (-0.4330127018922193*rhou1r[16])+0.4330127018922193*rhou1l[16]+0.25*rhou1r[9]+0.25*rhou1l[9]; 
  incr[10] = (-0.9682458365518543*rhou1r[18])-0.9682458365518543*rhou1l[18]+0.75*rhou1r[10]-0.75*rhou1l[10]-0.4330127018922193*rhou1r[5]-0.4330127018922193*rhou1l[5]; 
  incr[11] = 0.75*rhou1r[11]-0.75*rhou1l[11]-0.4330127018922194*rhou1r[7]-0.4330127018922194*rhou1l[7]; 
  incr[12] = 1.25*rhou1r[12]+1.25*rhou1l[12]-0.9682458365518543*rhou1r[4]+0.9682458365518543*rhou1l[4]+0.5590169943749476*rhou1r[1]+0.5590169943749476*rhou1l[1]; 
  incr[13] = (-0.4330127018922194*rhou1r[17])+0.4330127018922194*rhou1l[17]+0.25*rhou1r[13]+0.25*rhou1l[13]; 
  incr[14] = 1.25*rhou1r[14]+1.25*rhou1l[14]-0.9682458365518543*rhou1r[6]+0.9682458365518543*rhou1l[6]+0.5590169943749476*rhou1r[3]+0.5590169943749476*rhou1l[3]; 
  incr[15] = (-0.4330127018922194*rhou1r[19])+0.4330127018922194*rhou1l[19]+0.25*rhou1r[15]+0.25*rhou1l[15]; 
  incr[16] = 0.75*rhou1r[16]-0.75*rhou1l[16]-0.4330127018922194*rhou1r[9]-0.4330127018922194*rhou1l[9]; 
  incr[17] = 0.75*rhou1r[17]-0.75*rhou1l[17]-0.4330127018922193*rhou1r[13]-0.4330127018922193*rhou1l[13]; 
  incr[18] = 1.25*rhou1r[18]+1.25*rhou1l[18]-0.9682458365518545*rhou1r[10]+0.9682458365518545*rhou1l[10]+0.5590169943749475*rhou1r[5]+0.5590169943749475*rhou1l[5]; 
  incr[19] = 0.75*rhou1r[19]-0.75*rhou1l[19]-0.4330127018922193*rhou1r[15]-0.4330127018922193*rhou1l[15]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 
  outrrho[8] += incr[8]*dxr1; 
  outrrho[9] += incr[9]*dxr1; 
  outrrho[10] += incr[10]*dxr1; 
  outrrho[11] += incr[11]*dxr1; 
  outrrho[12] += incr[12]*dxr1; 
  outrrho[13] += incr[13]*dxr1; 
  outrrho[14] += incr[14]*dxr1; 
  outrrho[15] += incr[15]*dxr1; 
  outrrho[16] += incr[16]*dxr1; 
  outrrho[17] += incr[17]*dxr1; 
  outrrho[18] += incr[18]*dxr1; 
  outrrho[19] += incr[19]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += -1.0*incr[1]*dxl1; 
  outlrho[2] += incr[2]*dxl1; 
  outlrho[3] += -1.0*incr[3]*dxl1; 
  outlrho[4] += incr[4]*dxl1; 
  outlrho[5] += -1.0*incr[5]*dxl1; 
  outlrho[6] += incr[6]*dxl1; 
  outlrho[7] += -1.0*incr[7]*dxl1; 
  outlrho[8] += -1.0*incr[8]*dxl1; 
  outlrho[9] += -1.0*incr[9]*dxl1; 
  outlrho[10] += incr[10]*dxl1; 
  outlrho[11] += incr[11]*dxl1; 
  outlrho[12] += -1.0*incr[12]*dxl1; 
  outlrho[13] += -1.0*incr[13]*dxl1; 
  outlrho[14] += -1.0*incr[14]*dxl1; 
  outlrho[15] += -1.0*incr[15]*dxl1; 
  outlrho[16] += incr[16]*dxl1; 
  outlrho[17] += incr[17]*dxl1; 
  outlrho[18] += -1.0*incr[18]*dxl1; 
  outlrho[19] += incr[19]*dxl1; 

 
  //FluxRhoUx; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 
  outrrhoux[8] += incr[8]*dxr1; 
  outrrhoux[9] += incr[9]*dxr1; 
  outrrhoux[10] += incr[10]*dxr1; 
  outrrhoux[11] += incr[11]*dxr1; 
  outrrhoux[12] += incr[12]*dxr1; 
  outrrhoux[13] += incr[13]*dxr1; 
  outrrhoux[14] += incr[14]*dxr1; 
  outrrhoux[15] += incr[15]*dxr1; 
  outrrhoux[16] += incr[16]*dxr1; 
  outrrhoux[17] += incr[17]*dxr1; 
  outrrhoux[18] += incr[18]*dxr1; 
  outrrhoux[19] += incr[19]*dxr1; 


 
  //FluxRhoUy; 
  incr[0] = (-2.371708245126285*rhor[8]*vthsq)+2.371708245126285*rhol[8]*vthsq+0.6123724356957945*rhor[2]*vthsq+0.6123724356957945*rhol[2]*vthsq+1.060660171779821*rhor[19]*uvar1r[19]-0.6123724356957945*rhor[15]*uvar1r[19]+1.060660171779821*rhol[19]*uvar1l[19]+0.6123724356957945*rhol[15]*uvar1l[19]-0.6123724356957945*uvar1r[15]*rhor[19]+0.6123724356957945*uvar1l[15]*rhol[19]+1.767766952966369*rhor[18]*uvar1r[18]-1.369306393762915*rhor[10]*uvar1r[18]+0.7905694150420948*rhor[5]*uvar1r[18]+1.767766952966369*rhol[18]*uvar1l[18]+1.369306393762915*rhol[10]*uvar1l[18]+0.7905694150420948*rhol[5]*uvar1l[18]-1.369306393762915*uvar1r[10]*rhor[18]+0.7905694150420948*uvar1r[5]*rhor[18]+1.369306393762915*uvar1l[10]*rhol[18]+0.7905694150420948*uvar1l[5]*rhol[18]+1.060660171779821*rhor[17]*uvar1r[17]-0.6123724356957945*rhor[13]*uvar1r[17]+1.060660171779821*rhol[17]*uvar1l[17]+0.6123724356957945*rhol[13]*uvar1l[17]-0.6123724356957945*uvar1r[13]*rhor[17]+0.6123724356957945*uvar1l[13]*rhol[17]+1.060660171779821*rhor[16]*uvar1r[16]-0.6123724356957945*rhor[9]*uvar1r[16]+1.060660171779821*rhol[16]*uvar1l[16]+0.6123724356957945*rhol[9]*uvar1l[16]-0.6123724356957945*uvar1r[9]*rhor[16]+0.6123724356957945*uvar1l[9]*rhol[16]+0.3535533905932737*rhor[15]*uvar1r[15]+0.3535533905932737*rhol[15]*uvar1l[15]+1.767766952966369*rhor[14]*uvar1r[14]-1.369306393762915*rhor[6]*uvar1r[14]+0.7905694150420947*rhor[3]*uvar1r[14]+1.767766952966369*rhol[14]*uvar1l[14]+1.369306393762915*rhol[6]*uvar1l[14]+0.7905694150420947*rhol[3]*uvar1l[14]-1.369306393762915*uvar1r[6]*rhor[14]+0.7905694150420947*uvar1r[3]*rhor[14]+1.369306393762915*uvar1l[6]*rhol[14]+0.7905694150420947*uvar1l[3]*rhol[14]+0.3535533905932737*rhor[13]*uvar1r[13]+0.3535533905932737*rhol[13]*uvar1l[13]+1.767766952966369*rhor[12]*uvar1r[12]-1.369306393762915*rhor[4]*uvar1r[12]+0.7905694150420947*rhor[1]*uvar1r[12]+1.767766952966369*rhol[12]*uvar1l[12]+1.369306393762915*rhol[4]*uvar1l[12]+0.7905694150420947*rhol[1]*uvar1l[12]-1.369306393762915*uvar1r[4]*rhor[12]+0.7905694150420947*uvar1r[1]*rhor[12]+1.369306393762915*uvar1l[4]*rhol[12]+0.7905694150420947*uvar1l[1]*rhol[12]+1.060660171779821*rhor[11]*uvar1r[11]-0.6123724356957945*rhor[7]*uvar1r[11]+1.060660171779821*rhol[11]*uvar1l[11]+0.6123724356957945*rhol[7]*uvar1l[11]-0.6123724356957945*uvar1r[7]*rhor[11]+0.6123724356957945*uvar1l[7]*rhol[11]+1.060660171779821*rhor[10]*uvar1r[10]-0.6123724356957945*rhor[5]*uvar1r[10]+1.060660171779821*rhol[10]*uvar1l[10]+0.6123724356957945*rhol[5]*uvar1l[10]-0.6123724356957945*uvar1r[5]*rhor[10]+0.6123724356957945*uvar1l[5]*rhol[10]+0.3535533905932737*rhor[9]*uvar1r[9]+0.3535533905932737*rhol[9]*uvar1l[9]+1.767766952966369*rhor[8]*uvar1r[8]-1.369306393762915*rhor[2]*uvar1r[8]+0.7905694150420948*rhor[0]*uvar1r[8]+1.767766952966369*rhol[8]*uvar1l[8]+1.369306393762915*rhol[2]*uvar1l[8]+0.7905694150420948*rhol[0]*uvar1l[8]-1.369306393762915*uvar1r[2]*rhor[8]+0.7905694150420948*uvar1r[0]*rhor[8]+1.369306393762915*uvar1l[2]*rhol[8]+0.7905694150420948*uvar1l[0]*rhol[8]+0.3535533905932737*rhor[7]*uvar1r[7]+0.3535533905932737*rhol[7]*uvar1l[7]+1.060660171779821*rhor[6]*uvar1r[6]-0.6123724356957945*rhor[3]*uvar1r[6]+1.060660171779821*rhol[6]*uvar1l[6]+0.6123724356957945*rhol[3]*uvar1l[6]-0.6123724356957945*uvar1r[3]*rhor[6]+0.6123724356957945*uvar1l[3]*rhol[6]+0.3535533905932737*rhor[5]*uvar1r[5]+0.3535533905932737*rhol[5]*uvar1l[5]+1.060660171779821*rhor[4]*uvar1r[4]-0.6123724356957945*rhor[1]*uvar1r[4]+1.060660171779821*rhol[4]*uvar1l[4]+0.6123724356957945*rhol[1]*uvar1l[4]-0.6123724356957945*uvar1r[1]*rhor[4]+0.6123724356957945*uvar1l[1]*rhol[4]+0.3535533905932737*rhor[3]*uvar1r[3]+0.3535533905932737*rhol[3]*uvar1l[3]+1.060660171779821*rhor[2]*uvar1r[2]-0.6123724356957945*rhor[0]*uvar1r[2]+1.060660171779821*rhol[2]*uvar1l[2]+0.6123724356957945*rhol[0]*uvar1l[2]-0.6123724356957945*uvar1r[0]*rhor[2]+0.6123724356957945*uvar1l[0]*rhol[2]+0.3535533905932737*rhor[1]*uvar1r[1]+0.3535533905932737*rhol[1]*uvar1l[1]+0.3535533905932737*rhor[0]*uvar1r[0]+0.3535533905932737*rhol[0]*uvar1l[0]; 
  incr[2] = (-7.115124735378854*rhor[8]*vthsq)-7.115124735378854*rhol[8]*vthsq+1.837117307087383*rhor[2]*vthsq-1.837117307087383*rhol[2]*vthsq-1.837117307087383*rhor[19]*uvar1r[19]+1.060660171779821*rhor[15]*uvar1r[19]-1.837117307087383*rhol[19]*uvar1l[19]-1.060660171779821*rhol[15]*uvar1l[19]+1.060660171779821*uvar1r[15]*rhor[19]-1.060660171779821*uvar1l[15]*rhol[19]-3.061862178478972*rhor[18]*uvar1r[18]+2.371708245126284*rhor[10]*uvar1r[18]-1.369306393762915*rhor[5]*uvar1r[18]-3.061862178478972*rhol[18]*uvar1l[18]-2.371708245126284*rhol[10]*uvar1l[18]-1.369306393762915*rhol[5]*uvar1l[18]+2.371708245126284*uvar1r[10]*rhor[18]-1.369306393762915*uvar1r[5]*rhor[18]-2.371708245126284*uvar1l[10]*rhol[18]-1.369306393762915*uvar1l[5]*rhol[18]-1.837117307087383*rhor[17]*uvar1r[17]+1.060660171779821*rhor[13]*uvar1r[17]-1.837117307087383*rhol[17]*uvar1l[17]-1.060660171779821*rhol[13]*uvar1l[17]+1.060660171779821*uvar1r[13]*rhor[17]-1.060660171779821*uvar1l[13]*rhol[17]-1.837117307087383*rhor[16]*uvar1r[16]+1.060660171779821*rhor[9]*uvar1r[16]-1.837117307087383*rhol[16]*uvar1l[16]-1.060660171779821*rhol[9]*uvar1l[16]+1.060660171779821*uvar1r[9]*rhor[16]-1.060660171779821*uvar1l[9]*rhol[16]-0.6123724356957945*rhor[15]*uvar1r[15]-0.6123724356957945*rhol[15]*uvar1l[15]-3.061862178478972*rhor[14]*uvar1r[14]+2.371708245126285*rhor[6]*uvar1r[14]-1.369306393762915*rhor[3]*uvar1r[14]-3.061862178478972*rhol[14]*uvar1l[14]-2.371708245126285*rhol[6]*uvar1l[14]-1.369306393762915*rhol[3]*uvar1l[14]+2.371708245126285*uvar1r[6]*rhor[14]-1.369306393762915*uvar1r[3]*rhor[14]-2.371708245126285*uvar1l[6]*rhol[14]-1.369306393762915*uvar1l[3]*rhol[14]-0.6123724356957945*rhor[13]*uvar1r[13]-0.6123724356957945*rhol[13]*uvar1l[13]-3.061862178478972*rhor[12]*uvar1r[12]+2.371708245126285*rhor[4]*uvar1r[12]-1.369306393762915*rhor[1]*uvar1r[12]-3.061862178478972*rhol[12]*uvar1l[12]-2.371708245126285*rhol[4]*uvar1l[12]-1.369306393762915*rhol[1]*uvar1l[12]+2.371708245126285*uvar1r[4]*rhor[12]-1.369306393762915*uvar1r[1]*rhor[12]-2.371708245126285*uvar1l[4]*rhol[12]-1.369306393762915*uvar1l[1]*rhol[12]-1.837117307087383*rhor[11]*uvar1r[11]+1.060660171779821*rhor[7]*uvar1r[11]-1.837117307087383*rhol[11]*uvar1l[11]-1.060660171779821*rhol[7]*uvar1l[11]+1.060660171779821*uvar1r[7]*rhor[11]-1.060660171779821*uvar1l[7]*rhol[11]-1.837117307087383*rhor[10]*uvar1r[10]+1.060660171779821*rhor[5]*uvar1r[10]-1.837117307087383*rhol[10]*uvar1l[10]-1.060660171779821*rhol[5]*uvar1l[10]+1.060660171779821*uvar1r[5]*rhor[10]-1.060660171779821*uvar1l[5]*rhol[10]-0.6123724356957945*rhor[9]*uvar1r[9]-0.6123724356957945*rhol[9]*uvar1l[9]-3.061862178478972*rhor[8]*uvar1r[8]+2.371708245126284*rhor[2]*uvar1r[8]-1.369306393762915*rhor[0]*uvar1r[8]-3.061862178478972*rhol[8]*uvar1l[8]-2.371708245126284*rhol[2]*uvar1l[8]-1.369306393762915*rhol[0]*uvar1l[8]+2.371708245126284*uvar1r[2]*rhor[8]-1.369306393762915*uvar1r[0]*rhor[8]-2.371708245126284*uvar1l[2]*rhol[8]-1.369306393762915*uvar1l[0]*rhol[8]-0.6123724356957945*rhor[7]*uvar1r[7]-0.6123724356957945*rhol[7]*uvar1l[7]-1.837117307087383*rhor[6]*uvar1r[6]+1.060660171779821*rhor[3]*uvar1r[6]-1.837117307087383*rhol[6]*uvar1l[6]-1.060660171779821*rhol[3]*uvar1l[6]+1.060660171779821*uvar1r[3]*rhor[6]-1.060660171779821*uvar1l[3]*rhol[6]-0.6123724356957945*rhor[5]*uvar1r[5]-0.6123724356957945*rhol[5]*uvar1l[5]-1.837117307087383*rhor[4]*uvar1r[4]+1.060660171779821*rhor[1]*uvar1r[4]-1.837117307087383*rhol[4]*uvar1l[4]-1.060660171779821*rhol[1]*uvar1l[4]+1.060660171779821*uvar1r[1]*rhor[4]-1.060660171779821*uvar1l[1]*rhol[4]-0.6123724356957945*rhor[3]*uvar1r[3]-0.6123724356957945*rhol[3]*uvar1l[3]-1.837117307087383*rhor[2]*uvar1r[2]+1.060660171779821*rhor[0]*uvar1r[2]-1.837117307087383*rhol[2]*uvar1l[2]-1.060660171779821*rhol[0]*uvar1l[2]+1.060660171779821*uvar1r[0]*rhor[2]-1.060660171779821*uvar1l[0]*rhol[2]-0.6123724356957945*rhor[1]*uvar1r[1]-0.6123724356957945*rhol[1]*uvar1l[1]-0.6123724356957945*rhor[0]*uvar1r[0]-0.6123724356957945*rhol[0]*uvar1l[0]; 
  incr[8] = (-11.85854122563142*rhor[8]*vthsq)+11.85854122563142*rhol[8]*vthsq+3.061862178478972*rhor[2]*vthsq+3.061862178478972*rhol[2]*vthsq+2.371708245126285*rhor[19]*uvar1r[19]-1.369306393762915*rhor[15]*uvar1r[19]+2.371708245126285*rhol[19]*uvar1l[19]+1.369306393762915*rhol[15]*uvar1l[19]-1.369306393762915*uvar1r[15]*rhor[19]+1.369306393762915*uvar1l[15]*rhol[19]+3.952847075210474*rhor[18]*uvar1r[18]-3.061862178478972*rhor[10]*uvar1r[18]+1.767766952966369*rhor[5]*uvar1r[18]+3.952847075210474*rhol[18]*uvar1l[18]+3.061862178478972*rhol[10]*uvar1l[18]+1.767766952966369*rhol[5]*uvar1l[18]-3.061862178478972*uvar1r[10]*rhor[18]+1.767766952966369*uvar1r[5]*rhor[18]+3.061862178478972*uvar1l[10]*rhol[18]+1.767766952966369*uvar1l[5]*rhol[18]+2.371708245126285*rhor[17]*uvar1r[17]-1.369306393762915*rhor[13]*uvar1r[17]+2.371708245126285*rhol[17]*uvar1l[17]+1.369306393762915*rhol[13]*uvar1l[17]-1.369306393762915*uvar1r[13]*rhor[17]+1.369306393762915*uvar1l[13]*rhol[17]+2.371708245126285*rhor[16]*uvar1r[16]-1.369306393762915*rhor[9]*uvar1r[16]+2.371708245126285*rhol[16]*uvar1l[16]+1.369306393762915*rhol[9]*uvar1l[16]-1.369306393762915*uvar1r[9]*rhor[16]+1.369306393762915*uvar1l[9]*rhol[16]+0.7905694150420948*rhor[15]*uvar1r[15]+0.7905694150420948*rhol[15]*uvar1l[15]+3.952847075210474*rhor[14]*uvar1r[14]-3.061862178478973*rhor[6]*uvar1r[14]+1.767766952966368*rhor[3]*uvar1r[14]+3.952847075210474*rhol[14]*uvar1l[14]+3.061862178478973*rhol[6]*uvar1l[14]+1.767766952966368*rhol[3]*uvar1l[14]-3.061862178478973*uvar1r[6]*rhor[14]+1.767766952966368*uvar1r[3]*rhor[14]+3.061862178478973*uvar1l[6]*rhol[14]+1.767766952966368*uvar1l[3]*rhol[14]+0.7905694150420948*rhor[13]*uvar1r[13]+0.7905694150420948*rhol[13]*uvar1l[13]+3.952847075210474*rhor[12]*uvar1r[12]-3.061862178478973*rhor[4]*uvar1r[12]+1.767766952966368*rhor[1]*uvar1r[12]+3.952847075210474*rhol[12]*uvar1l[12]+3.061862178478973*rhol[4]*uvar1l[12]+1.767766952966368*rhol[1]*uvar1l[12]-3.061862178478973*uvar1r[4]*rhor[12]+1.767766952966368*uvar1r[1]*rhor[12]+3.061862178478973*uvar1l[4]*rhol[12]+1.767766952966368*uvar1l[1]*rhol[12]+2.371708245126285*rhor[11]*uvar1r[11]-1.369306393762915*rhor[7]*uvar1r[11]+2.371708245126285*rhol[11]*uvar1l[11]+1.369306393762915*rhol[7]*uvar1l[11]-1.369306393762915*uvar1r[7]*rhor[11]+1.369306393762915*uvar1l[7]*rhol[11]+2.371708245126285*rhor[10]*uvar1r[10]-1.369306393762915*rhor[5]*uvar1r[10]+2.371708245126285*rhol[10]*uvar1l[10]+1.369306393762915*rhol[5]*uvar1l[10]-1.369306393762915*uvar1r[5]*rhor[10]+1.369306393762915*uvar1l[5]*rhol[10]+0.7905694150420948*rhor[9]*uvar1r[9]+0.7905694150420948*rhol[9]*uvar1l[9]+3.952847075210474*rhor[8]*uvar1r[8]-3.061862178478972*rhor[2]*uvar1r[8]+1.767766952966369*rhor[0]*uvar1r[8]+3.952847075210474*rhol[8]*uvar1l[8]+3.061862178478972*rhol[2]*uvar1l[8]+1.767766952966369*rhol[0]*uvar1l[8]-3.061862178478972*uvar1r[2]*rhor[8]+1.767766952966369*uvar1r[0]*rhor[8]+3.061862178478972*uvar1l[2]*rhol[8]+1.767766952966369*uvar1l[0]*rhol[8]+0.7905694150420948*rhor[7]*uvar1r[7]+0.7905694150420948*rhol[7]*uvar1l[7]+2.371708245126285*rhor[6]*uvar1r[6]-1.369306393762915*rhor[3]*uvar1r[6]+2.371708245126285*rhol[6]*uvar1l[6]+1.369306393762915*rhol[3]*uvar1l[6]-1.369306393762915*uvar1r[3]*rhor[6]+1.369306393762915*uvar1l[3]*rhol[6]+0.7905694150420948*rhor[5]*uvar1r[5]+0.7905694150420948*rhol[5]*uvar1l[5]+2.371708245126285*rhor[4]*uvar1r[4]-1.369306393762915*rhor[1]*uvar1r[4]+2.371708245126285*rhol[4]*uvar1l[4]+1.369306393762915*rhol[1]*uvar1l[4]-1.369306393762915*uvar1r[1]*rhor[4]+1.369306393762915*uvar1l[1]*rhol[4]+0.7905694150420948*rhor[3]*uvar1r[3]+0.7905694150420948*rhol[3]*uvar1l[3]+2.371708245126285*rhor[2]*uvar1r[2]-1.369306393762915*rhor[0]*uvar1r[2]+2.371708245126285*rhol[2]*uvar1l[2]+1.369306393762915*rhol[0]*uvar1l[2]-1.369306393762915*uvar1r[0]*rhor[2]+1.369306393762915*uvar1l[0]*rhol[2]+0.7905694150420948*rhor[1]*uvar1r[1]+0.7905694150420948*rhol[1]*uvar1l[1]+0.7905694150420948*rhor[0]*uvar1r[0]+0.7905694150420948*rhol[0]*uvar1l[0]; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 
  outrrhouy[8] += incr[8]*dxr1; 
  outrrhouy[9] += incr[9]*dxr1; 
  outrrhouy[10] += incr[10]*dxr1; 
  outrrhouy[11] += incr[11]*dxr1; 
  outrrhouy[12] += incr[12]*dxr1; 
  outrrhouy[13] += incr[13]*dxr1; 
  outrrhouy[14] += incr[14]*dxr1; 
  outrrhouy[15] += incr[15]*dxr1; 
  outrrhouy[16] += incr[16]*dxr1; 
  outrrhouy[17] += incr[17]*dxr1; 
  outrrhouy[18] += incr[18]*dxr1; 
  outrrhouy[19] += incr[19]*dxr1; 

  outlrhouy[0] += -1.0*incr[0]*dxl1; 
  outlrhouy[2] += incr[2]*dxl1; 
  outlrhouy[8] += -1.0*incr[8]*dxl1; 

 
  //FluxRhoUz; 

  outrrhouz[0] += incr[0]*dxr1; 
  outrrhouz[1] += incr[1]*dxr1; 
  outrrhouz[2] += incr[2]*dxr1; 
  outrrhouz[3] += incr[3]*dxr1; 
  outrrhouz[4] += incr[4]*dxr1; 
  outrrhouz[5] += incr[5]*dxr1; 
  outrrhouz[6] += incr[6]*dxr1; 
  outrrhouz[7] += incr[7]*dxr1; 
  outrrhouz[8] += incr[8]*dxr1; 
  outrrhouz[9] += incr[9]*dxr1; 
  outrrhouz[10] += incr[10]*dxr1; 
  outrrhouz[11] += incr[11]*dxr1; 
  outrrhouz[12] += incr[12]*dxr1; 
  outrrhouz[13] += incr[13]*dxr1; 
  outrrhouz[14] += incr[14]*dxr1; 
  outrrhouz[15] += incr[15]*dxr1; 
  outrrhouz[16] += incr[16]*dxr1; 
  outrrhouz[17] += incr[17]*dxr1; 
  outrrhouz[18] += incr[18]*dxr1; 
  outrrhouz[19] += incr[19]*dxr1; 


 
} 
GKYL_CU_DH void isoeuler_surfz_3x_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[2]; 
  const double dxr1 = 2.0/dxv[2]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[20]; 
  const double *rhou1l = &statevecl[40]; 
  const double *rhou2l = &statevecl[60]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[20]; 
  const double *rhou1r = &statevecr[40]; 
  const double *rhou2r = &statevecr[60]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[20]; 
  const double *uvar2l = &uvarl[40]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[20]; 
  const double *uvar2r = &uvarr[40]; 
  double *outlrho = &out[0]; 
  double *outlrhoux = &out[20]; 
  double *outlrhouy = &out[40]; 
  double *outlrhouz = &out[60]; 
  double *outrrho = &out[80]; 
  double *outrrhoux = &out[100]; 
  double *outrrhouy = &out[120]; 
  double *outrrhouz = &out[140]; 
  double incr[20]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  incr[0] = 0.5590169943749475*rhou2r[9]+0.5590169943749475*rhou2l[9]-0.4330127018922193*rhou2r[3]+0.4330127018922193*rhou2l[3]+0.25*rhou2r[0]+0.25*rhou2l[0]; 
  incr[1] = 0.5590169943749476*rhou2r[15]+0.5590169943749476*rhou2l[15]-0.4330127018922193*rhou2r[5]+0.4330127018922193*rhou2l[5]+0.25*rhou2r[1]+0.25*rhou2l[1]; 
  incr[2] = 0.5590169943749476*rhou2r[16]+0.5590169943749476*rhou2l[16]-0.4330127018922193*rhou2r[6]+0.4330127018922193*rhou2l[6]+0.25*rhou2r[2]+0.25*rhou2l[2]; 
  incr[3] = (-0.9682458365518543*rhou2r[9])-0.9682458365518543*rhou2l[9]+0.75*rhou2r[3]-0.75*rhou2l[3]-0.4330127018922193*rhou2r[0]-0.4330127018922193*rhou2l[0]; 
  incr[4] = 0.5590169943749475*rhou2r[19]+0.5590169943749475*rhou2l[19]-0.4330127018922194*rhou2r[10]+0.4330127018922194*rhou2l[10]+0.25*rhou2r[4]+0.25*rhou2l[4]; 
  incr[5] = (-0.9682458365518543*rhou2r[15])-0.9682458365518543*rhou2l[15]+0.75*rhou2r[5]-0.75*rhou2l[5]-0.4330127018922194*rhou2r[1]-0.4330127018922194*rhou2l[1]; 
  incr[6] = (-0.9682458365518543*rhou2r[16])-0.9682458365518543*rhou2l[16]+0.75*rhou2r[6]-0.75*rhou2l[6]-0.4330127018922194*rhou2r[2]-0.4330127018922194*rhou2l[2]; 
  incr[7] = (-0.4330127018922193*rhou2r[13])+0.4330127018922193*rhou2l[13]+0.25*rhou2r[7]+0.25*rhou2l[7]; 
  incr[8] = (-0.4330127018922193*rhou2r[14])+0.4330127018922193*rhou2l[14]+0.25*rhou2r[8]+0.25*rhou2l[8]; 
  incr[9] = 1.25*rhou2r[9]+1.25*rhou2l[9]-0.9682458365518543*rhou2r[3]+0.9682458365518543*rhou2l[3]+0.5590169943749475*rhou2r[0]+0.5590169943749475*rhou2l[0]; 
  incr[10] = (-0.9682458365518543*rhou2r[19])-0.9682458365518543*rhou2l[19]+0.75*rhou2r[10]-0.75*rhou2l[10]-0.4330127018922193*rhou2r[4]-0.4330127018922193*rhou2l[4]; 
  incr[11] = (-0.4330127018922194*rhou2r[17])+0.4330127018922194*rhou2l[17]+0.25*rhou2r[11]+0.25*rhou2l[11]; 
  incr[12] = (-0.4330127018922194*rhou2r[18])+0.4330127018922194*rhou2l[18]+0.25*rhou2r[12]+0.25*rhou2l[12]; 
  incr[13] = 0.75*rhou2r[13]-0.75*rhou2l[13]-0.4330127018922194*rhou2r[7]-0.4330127018922194*rhou2l[7]; 
  incr[14] = 0.75*rhou2r[14]-0.75*rhou2l[14]-0.4330127018922194*rhou2r[8]-0.4330127018922194*rhou2l[8]; 
  incr[15] = 1.25*rhou2r[15]+1.25*rhou2l[15]-0.9682458365518543*rhou2r[5]+0.9682458365518543*rhou2l[5]+0.5590169943749476*rhou2r[1]+0.5590169943749476*rhou2l[1]; 
  incr[16] = 1.25*rhou2r[16]+1.25*rhou2l[16]-0.9682458365518543*rhou2r[6]+0.9682458365518543*rhou2l[6]+0.5590169943749476*rhou2r[2]+0.5590169943749476*rhou2l[2]; 
  incr[17] = 0.75*rhou2r[17]-0.75*rhou2l[17]-0.4330127018922193*rhou2r[11]-0.4330127018922193*rhou2l[11]; 
  incr[18] = 0.75*rhou2r[18]-0.75*rhou2l[18]-0.4330127018922193*rhou2r[12]-0.4330127018922193*rhou2l[12]; 
  incr[19] = 1.25*rhou2r[19]+1.25*rhou2l[19]-0.9682458365518545*rhou2r[10]+0.9682458365518545*rhou2l[10]+0.5590169943749475*rhou2r[4]+0.5590169943749475*rhou2l[4]; 

  outrrho[0] += incr[0]*dxr1; 
  outrrho[1] += incr[1]*dxr1; 
  outrrho[2] += incr[2]*dxr1; 
  outrrho[3] += incr[3]*dxr1; 
  outrrho[4] += incr[4]*dxr1; 
  outrrho[5] += incr[5]*dxr1; 
  outrrho[6] += incr[6]*dxr1; 
  outrrho[7] += incr[7]*dxr1; 
  outrrho[8] += incr[8]*dxr1; 
  outrrho[9] += incr[9]*dxr1; 
  outrrho[10] += incr[10]*dxr1; 
  outrrho[11] += incr[11]*dxr1; 
  outrrho[12] += incr[12]*dxr1; 
  outrrho[13] += incr[13]*dxr1; 
  outrrho[14] += incr[14]*dxr1; 
  outrrho[15] += incr[15]*dxr1; 
  outrrho[16] += incr[16]*dxr1; 
  outrrho[17] += incr[17]*dxr1; 
  outrrho[18] += incr[18]*dxr1; 
  outrrho[19] += incr[19]*dxr1; 

  outlrho[0] += -1.0*incr[0]*dxl1; 
  outlrho[1] += -1.0*incr[1]*dxl1; 
  outlrho[2] += -1.0*incr[2]*dxl1; 
  outlrho[3] += incr[3]*dxl1; 
  outlrho[4] += -1.0*incr[4]*dxl1; 
  outlrho[5] += incr[5]*dxl1; 
  outlrho[6] += incr[6]*dxl1; 
  outlrho[7] += -1.0*incr[7]*dxl1; 
  outlrho[8] += -1.0*incr[8]*dxl1; 
  outlrho[9] += -1.0*incr[9]*dxl1; 
  outlrho[10] += incr[10]*dxl1; 
  outlrho[11] += -1.0*incr[11]*dxl1; 
  outlrho[12] += -1.0*incr[12]*dxl1; 
  outlrho[13] += incr[13]*dxl1; 
  outlrho[14] += incr[14]*dxl1; 
  outlrho[15] += -1.0*incr[15]*dxl1; 
  outlrho[16] += -1.0*incr[16]*dxl1; 
  outlrho[17] += incr[17]*dxl1; 
  outlrho[18] += incr[18]*dxl1; 
  outlrho[19] += -1.0*incr[19]*dxl1; 

 
  //FluxRhoUx; 

  outrrhoux[0] += incr[0]*dxr1; 
  outrrhoux[1] += incr[1]*dxr1; 
  outrrhoux[2] += incr[2]*dxr1; 
  outrrhoux[3] += incr[3]*dxr1; 
  outrrhoux[4] += incr[4]*dxr1; 
  outrrhoux[5] += incr[5]*dxr1; 
  outrrhoux[6] += incr[6]*dxr1; 
  outrrhoux[7] += incr[7]*dxr1; 
  outrrhoux[8] += incr[8]*dxr1; 
  outrrhoux[9] += incr[9]*dxr1; 
  outrrhoux[10] += incr[10]*dxr1; 
  outrrhoux[11] += incr[11]*dxr1; 
  outrrhoux[12] += incr[12]*dxr1; 
  outrrhoux[13] += incr[13]*dxr1; 
  outrrhoux[14] += incr[14]*dxr1; 
  outrrhoux[15] += incr[15]*dxr1; 
  outrrhoux[16] += incr[16]*dxr1; 
  outrrhoux[17] += incr[17]*dxr1; 
  outrrhoux[18] += incr[18]*dxr1; 
  outrrhoux[19] += incr[19]*dxr1; 


 
  //FluxRhoUy; 

  outrrhouy[0] += incr[0]*dxr1; 
  outrrhouy[1] += incr[1]*dxr1; 
  outrrhouy[2] += incr[2]*dxr1; 
  outrrhouy[3] += incr[3]*dxr1; 
  outrrhouy[4] += incr[4]*dxr1; 
  outrrhouy[5] += incr[5]*dxr1; 
  outrrhouy[6] += incr[6]*dxr1; 
  outrrhouy[7] += incr[7]*dxr1; 
  outrrhouy[8] += incr[8]*dxr1; 
  outrrhouy[9] += incr[9]*dxr1; 
  outrrhouy[10] += incr[10]*dxr1; 
  outrrhouy[11] += incr[11]*dxr1; 
  outrrhouy[12] += incr[12]*dxr1; 
  outrrhouy[13] += incr[13]*dxr1; 
  outrrhouy[14] += incr[14]*dxr1; 
  outrrhouy[15] += incr[15]*dxr1; 
  outrrhouy[16] += incr[16]*dxr1; 
  outrrhouy[17] += incr[17]*dxr1; 
  outrrhouy[18] += incr[18]*dxr1; 
  outrrhouy[19] += incr[19]*dxr1; 


 
  //FluxRhoUz; 
  incr[0] = (-2.371708245126285*rhor[9]*vthsq)+2.371708245126285*rhol[9]*vthsq+0.6123724356957945*rhor[3]*vthsq+0.6123724356957945*rhol[3]*vthsq+1.767766952966369*rhor[19]*uvar2r[19]-1.369306393762915*rhor[10]*uvar2r[19]+0.7905694150420948*rhor[4]*uvar2r[19]+1.767766952966369*rhol[19]*uvar2l[19]+1.369306393762915*rhol[10]*uvar2l[19]+0.7905694150420948*rhol[4]*uvar2l[19]-1.369306393762915*uvar2r[10]*rhor[19]+0.7905694150420948*uvar2r[4]*rhor[19]+1.369306393762915*uvar2l[10]*rhol[19]+0.7905694150420948*uvar2l[4]*rhol[19]+1.060660171779821*rhor[18]*uvar2r[18]-0.6123724356957945*rhor[12]*uvar2r[18]+1.060660171779821*rhol[18]*uvar2l[18]+0.6123724356957945*rhol[12]*uvar2l[18]-0.6123724356957945*uvar2r[12]*rhor[18]+0.6123724356957945*uvar2l[12]*rhol[18]+1.060660171779821*rhor[17]*uvar2r[17]-0.6123724356957945*rhor[11]*uvar2r[17]+1.060660171779821*rhol[17]*uvar2l[17]+0.6123724356957945*rhol[11]*uvar2l[17]-0.6123724356957945*uvar2r[11]*rhor[17]+0.6123724356957945*uvar2l[11]*rhol[17]+1.767766952966369*rhor[16]*uvar2r[16]-1.369306393762915*rhor[6]*uvar2r[16]+0.7905694150420947*rhor[2]*uvar2r[16]+1.767766952966369*rhol[16]*uvar2l[16]+1.369306393762915*rhol[6]*uvar2l[16]+0.7905694150420947*rhol[2]*uvar2l[16]-1.369306393762915*uvar2r[6]*rhor[16]+0.7905694150420947*uvar2r[2]*rhor[16]+1.369306393762915*uvar2l[6]*rhol[16]+0.7905694150420947*uvar2l[2]*rhol[16]+1.767766952966369*rhor[15]*uvar2r[15]-1.369306393762915*rhor[5]*uvar2r[15]+0.7905694150420947*rhor[1]*uvar2r[15]+1.767766952966369*rhol[15]*uvar2l[15]+1.369306393762915*rhol[5]*uvar2l[15]+0.7905694150420947*rhol[1]*uvar2l[15]-1.369306393762915*uvar2r[5]*rhor[15]+0.7905694150420947*uvar2r[1]*rhor[15]+1.369306393762915*uvar2l[5]*rhol[15]+0.7905694150420947*uvar2l[1]*rhol[15]+1.060660171779821*rhor[14]*uvar2r[14]-0.6123724356957945*rhor[8]*uvar2r[14]+1.060660171779821*rhol[14]*uvar2l[14]+0.6123724356957945*rhol[8]*uvar2l[14]-0.6123724356957945*uvar2r[8]*rhor[14]+0.6123724356957945*uvar2l[8]*rhol[14]+1.060660171779821*rhor[13]*uvar2r[13]-0.6123724356957945*rhor[7]*uvar2r[13]+1.060660171779821*rhol[13]*uvar2l[13]+0.6123724356957945*rhol[7]*uvar2l[13]-0.6123724356957945*uvar2r[7]*rhor[13]+0.6123724356957945*uvar2l[7]*rhol[13]+0.3535533905932737*rhor[12]*uvar2r[12]+0.3535533905932737*rhol[12]*uvar2l[12]+0.3535533905932737*rhor[11]*uvar2r[11]+0.3535533905932737*rhol[11]*uvar2l[11]+1.060660171779821*rhor[10]*uvar2r[10]-0.6123724356957945*rhor[4]*uvar2r[10]+1.060660171779821*rhol[10]*uvar2l[10]+0.6123724356957945*rhol[4]*uvar2l[10]-0.6123724356957945*uvar2r[4]*rhor[10]+0.6123724356957945*uvar2l[4]*rhol[10]+1.767766952966369*rhor[9]*uvar2r[9]-1.369306393762915*rhor[3]*uvar2r[9]+0.7905694150420948*rhor[0]*uvar2r[9]+1.767766952966369*rhol[9]*uvar2l[9]+1.369306393762915*rhol[3]*uvar2l[9]+0.7905694150420948*rhol[0]*uvar2l[9]-1.369306393762915*uvar2r[3]*rhor[9]+0.7905694150420948*uvar2r[0]*rhor[9]+1.369306393762915*uvar2l[3]*rhol[9]+0.7905694150420948*uvar2l[0]*rhol[9]+0.3535533905932737*rhor[8]*uvar2r[8]+0.3535533905932737*rhol[8]*uvar2l[8]+0.3535533905932737*rhor[7]*uvar2r[7]+0.3535533905932737*rhol[7]*uvar2l[7]+1.060660171779821*rhor[6]*uvar2r[6]-0.6123724356957945*rhor[2]*uvar2r[6]+1.060660171779821*rhol[6]*uvar2l[6]+0.6123724356957945*rhol[2]*uvar2l[6]-0.6123724356957945*uvar2r[2]*rhor[6]+0.6123724356957945*uvar2l[2]*rhol[6]+1.060660171779821*rhor[5]*uvar2r[5]-0.6123724356957945*rhor[1]*uvar2r[5]+1.060660171779821*rhol[5]*uvar2l[5]+0.6123724356957945*rhol[1]*uvar2l[5]-0.6123724356957945*uvar2r[1]*rhor[5]+0.6123724356957945*uvar2l[1]*rhol[5]+0.3535533905932737*rhor[4]*uvar2r[4]+0.3535533905932737*rhol[4]*uvar2l[4]+1.060660171779821*rhor[3]*uvar2r[3]-0.6123724356957945*rhor[0]*uvar2r[3]+1.060660171779821*rhol[3]*uvar2l[3]+0.6123724356957945*rhol[0]*uvar2l[3]-0.6123724356957945*uvar2r[0]*rhor[3]+0.6123724356957945*uvar2l[0]*rhol[3]+0.3535533905932737*rhor[2]*uvar2r[2]+0.3535533905932737*rhol[2]*uvar2l[2]+0.3535533905932737*rhor[1]*uvar2r[1]+0.3535533905932737*rhol[1]*uvar2l[1]+0.3535533905932737*rhor[0]*uvar2r[0]+0.3535533905932737*rhol[0]*uvar2l[0]; 
  incr[3] = (-7.115124735378854*rhor[9]*vthsq)-7.115124735378854*rhol[9]*vthsq+1.837117307087383*rhor[3]*vthsq-1.837117307087383*rhol[3]*vthsq-3.061862178478972*rhor[19]*uvar2r[19]+2.371708245126284*rhor[10]*uvar2r[19]-1.369306393762915*rhor[4]*uvar2r[19]-3.061862178478972*rhol[19]*uvar2l[19]-2.371708245126284*rhol[10]*uvar2l[19]-1.369306393762915*rhol[4]*uvar2l[19]+2.371708245126284*uvar2r[10]*rhor[19]-1.369306393762915*uvar2r[4]*rhor[19]-2.371708245126284*uvar2l[10]*rhol[19]-1.369306393762915*uvar2l[4]*rhol[19]-1.837117307087383*rhor[18]*uvar2r[18]+1.060660171779821*rhor[12]*uvar2r[18]-1.837117307087383*rhol[18]*uvar2l[18]-1.060660171779821*rhol[12]*uvar2l[18]+1.060660171779821*uvar2r[12]*rhor[18]-1.060660171779821*uvar2l[12]*rhol[18]-1.837117307087383*rhor[17]*uvar2r[17]+1.060660171779821*rhor[11]*uvar2r[17]-1.837117307087383*rhol[17]*uvar2l[17]-1.060660171779821*rhol[11]*uvar2l[17]+1.060660171779821*uvar2r[11]*rhor[17]-1.060660171779821*uvar2l[11]*rhol[17]-3.061862178478972*rhor[16]*uvar2r[16]+2.371708245126285*rhor[6]*uvar2r[16]-1.369306393762915*rhor[2]*uvar2r[16]-3.061862178478972*rhol[16]*uvar2l[16]-2.371708245126285*rhol[6]*uvar2l[16]-1.369306393762915*rhol[2]*uvar2l[16]+2.371708245126285*uvar2r[6]*rhor[16]-1.369306393762915*uvar2r[2]*rhor[16]-2.371708245126285*uvar2l[6]*rhol[16]-1.369306393762915*uvar2l[2]*rhol[16]-3.061862178478972*rhor[15]*uvar2r[15]+2.371708245126285*rhor[5]*uvar2r[15]-1.369306393762915*rhor[1]*uvar2r[15]-3.061862178478972*rhol[15]*uvar2l[15]-2.371708245126285*rhol[5]*uvar2l[15]-1.369306393762915*rhol[1]*uvar2l[15]+2.371708245126285*uvar2r[5]*rhor[15]-1.369306393762915*uvar2r[1]*rhor[15]-2.371708245126285*uvar2l[5]*rhol[15]-1.369306393762915*uvar2l[1]*rhol[15]-1.837117307087383*rhor[14]*uvar2r[14]+1.060660171779821*rhor[8]*uvar2r[14]-1.837117307087383*rhol[14]*uvar2l[14]-1.060660171779821*rhol[8]*uvar2l[14]+1.060660171779821*uvar2r[8]*rhor[14]-1.060660171779821*uvar2l[8]*rhol[14]-1.837117307087383*rhor[13]*uvar2r[13]+1.060660171779821*rhor[7]*uvar2r[13]-1.837117307087383*rhol[13]*uvar2l[13]-1.060660171779821*rhol[7]*uvar2l[13]+1.060660171779821*uvar2r[7]*rhor[13]-1.060660171779821*uvar2l[7]*rhol[13]-0.6123724356957945*rhor[12]*uvar2r[12]-0.6123724356957945*rhol[12]*uvar2l[12]-0.6123724356957945*rhor[11]*uvar2r[11]-0.6123724356957945*rhol[11]*uvar2l[11]-1.837117307087383*rhor[10]*uvar2r[10]+1.060660171779821*rhor[4]*uvar2r[10]-1.837117307087383*rhol[10]*uvar2l[10]-1.060660171779821*rhol[4]*uvar2l[10]+1.060660171779821*uvar2r[4]*rhor[10]-1.060660171779821*uvar2l[4]*rhol[10]-3.061862178478972*rhor[9]*uvar2r[9]+2.371708245126284*rhor[3]*uvar2r[9]-1.369306393762915*rhor[0]*uvar2r[9]-3.061862178478972*rhol[9]*uvar2l[9]-2.371708245126284*rhol[3]*uvar2l[9]-1.369306393762915*rhol[0]*uvar2l[9]+2.371708245126284*uvar2r[3]*rhor[9]-1.369306393762915*uvar2r[0]*rhor[9]-2.371708245126284*uvar2l[3]*rhol[9]-1.369306393762915*uvar2l[0]*rhol[9]-0.6123724356957945*rhor[8]*uvar2r[8]-0.6123724356957945*rhol[8]*uvar2l[8]-0.6123724356957945*rhor[7]*uvar2r[7]-0.6123724356957945*rhol[7]*uvar2l[7]-1.837117307087383*rhor[6]*uvar2r[6]+1.060660171779821*rhor[2]*uvar2r[6]-1.837117307087383*rhol[6]*uvar2l[6]-1.060660171779821*rhol[2]*uvar2l[6]+1.060660171779821*uvar2r[2]*rhor[6]-1.060660171779821*uvar2l[2]*rhol[6]-1.837117307087383*rhor[5]*uvar2r[5]+1.060660171779821*rhor[1]*uvar2r[5]-1.837117307087383*rhol[5]*uvar2l[5]-1.060660171779821*rhol[1]*uvar2l[5]+1.060660171779821*uvar2r[1]*rhor[5]-1.060660171779821*uvar2l[1]*rhol[5]-0.6123724356957945*rhor[4]*uvar2r[4]-0.6123724356957945*rhol[4]*uvar2l[4]-1.837117307087383*rhor[3]*uvar2r[3]+1.060660171779821*rhor[0]*uvar2r[3]-1.837117307087383*rhol[3]*uvar2l[3]-1.060660171779821*rhol[0]*uvar2l[3]+1.060660171779821*uvar2r[0]*rhor[3]-1.060660171779821*uvar2l[0]*rhol[3]-0.6123724356957945*rhor[2]*uvar2r[2]-0.6123724356957945*rhol[2]*uvar2l[2]-0.6123724356957945*rhor[1]*uvar2r[1]-0.6123724356957945*rhol[1]*uvar2l[1]-0.6123724356957945*rhor[0]*uvar2r[0]-0.6123724356957945*rhol[0]*uvar2l[0]; 
  incr[9] = (-11.85854122563142*rhor[9]*vthsq)+11.85854122563142*rhol[9]*vthsq+3.061862178478972*rhor[3]*vthsq+3.061862178478972*rhol[3]*vthsq+3.952847075210474*rhor[19]*uvar2r[19]-3.061862178478972*rhor[10]*uvar2r[19]+1.767766952966369*rhor[4]*uvar2r[19]+3.952847075210474*rhol[19]*uvar2l[19]+3.061862178478972*rhol[10]*uvar2l[19]+1.767766952966369*rhol[4]*uvar2l[19]-3.061862178478972*uvar2r[10]*rhor[19]+1.767766952966369*uvar2r[4]*rhor[19]+3.061862178478972*uvar2l[10]*rhol[19]+1.767766952966369*uvar2l[4]*rhol[19]+2.371708245126285*rhor[18]*uvar2r[18]-1.369306393762915*rhor[12]*uvar2r[18]+2.371708245126285*rhol[18]*uvar2l[18]+1.369306393762915*rhol[12]*uvar2l[18]-1.369306393762915*uvar2r[12]*rhor[18]+1.369306393762915*uvar2l[12]*rhol[18]+2.371708245126285*rhor[17]*uvar2r[17]-1.369306393762915*rhor[11]*uvar2r[17]+2.371708245126285*rhol[17]*uvar2l[17]+1.369306393762915*rhol[11]*uvar2l[17]-1.369306393762915*uvar2r[11]*rhor[17]+1.369306393762915*uvar2l[11]*rhol[17]+3.952847075210474*rhor[16]*uvar2r[16]-3.061862178478973*rhor[6]*uvar2r[16]+1.767766952966368*rhor[2]*uvar2r[16]+3.952847075210474*rhol[16]*uvar2l[16]+3.061862178478973*rhol[6]*uvar2l[16]+1.767766952966368*rhol[2]*uvar2l[16]-3.061862178478973*uvar2r[6]*rhor[16]+1.767766952966368*uvar2r[2]*rhor[16]+3.061862178478973*uvar2l[6]*rhol[16]+1.767766952966368*uvar2l[2]*rhol[16]+3.952847075210474*rhor[15]*uvar2r[15]-3.061862178478973*rhor[5]*uvar2r[15]+1.767766952966368*rhor[1]*uvar2r[15]+3.952847075210474*rhol[15]*uvar2l[15]+3.061862178478973*rhol[5]*uvar2l[15]+1.767766952966368*rhol[1]*uvar2l[15]-3.061862178478973*uvar2r[5]*rhor[15]+1.767766952966368*uvar2r[1]*rhor[15]+3.061862178478973*uvar2l[5]*rhol[15]+1.767766952966368*uvar2l[1]*rhol[15]+2.371708245126285*rhor[14]*uvar2r[14]-1.369306393762915*rhor[8]*uvar2r[14]+2.371708245126285*rhol[14]*uvar2l[14]+1.369306393762915*rhol[8]*uvar2l[14]-1.369306393762915*uvar2r[8]*rhor[14]+1.369306393762915*uvar2l[8]*rhol[14]+2.371708245126285*rhor[13]*uvar2r[13]-1.369306393762915*rhor[7]*uvar2r[13]+2.371708245126285*rhol[13]*uvar2l[13]+1.369306393762915*rhol[7]*uvar2l[13]-1.369306393762915*uvar2r[7]*rhor[13]+1.369306393762915*uvar2l[7]*rhol[13]+0.7905694150420948*rhor[12]*uvar2r[12]+0.7905694150420948*rhol[12]*uvar2l[12]+0.7905694150420948*rhor[11]*uvar2r[11]+0.7905694150420948*rhol[11]*uvar2l[11]+2.371708245126285*rhor[10]*uvar2r[10]-1.369306393762915*rhor[4]*uvar2r[10]+2.371708245126285*rhol[10]*uvar2l[10]+1.369306393762915*rhol[4]*uvar2l[10]-1.369306393762915*uvar2r[4]*rhor[10]+1.369306393762915*uvar2l[4]*rhol[10]+3.952847075210474*rhor[9]*uvar2r[9]-3.061862178478972*rhor[3]*uvar2r[9]+1.767766952966369*rhor[0]*uvar2r[9]+3.952847075210474*rhol[9]*uvar2l[9]+3.061862178478972*rhol[3]*uvar2l[9]+1.767766952966369*rhol[0]*uvar2l[9]-3.061862178478972*uvar2r[3]*rhor[9]+1.767766952966369*uvar2r[0]*rhor[9]+3.061862178478972*uvar2l[3]*rhol[9]+1.767766952966369*uvar2l[0]*rhol[9]+0.7905694150420948*rhor[8]*uvar2r[8]+0.7905694150420948*rhol[8]*uvar2l[8]+0.7905694150420948*rhor[7]*uvar2r[7]+0.7905694150420948*rhol[7]*uvar2l[7]+2.371708245126285*rhor[6]*uvar2r[6]-1.369306393762915*rhor[2]*uvar2r[6]+2.371708245126285*rhol[6]*uvar2l[6]+1.369306393762915*rhol[2]*uvar2l[6]-1.369306393762915*uvar2r[2]*rhor[6]+1.369306393762915*uvar2l[2]*rhol[6]+2.371708245126285*rhor[5]*uvar2r[5]-1.369306393762915*rhor[1]*uvar2r[5]+2.371708245126285*rhol[5]*uvar2l[5]+1.369306393762915*rhol[1]*uvar2l[5]-1.369306393762915*uvar2r[1]*rhor[5]+1.369306393762915*uvar2l[1]*rhol[5]+0.7905694150420948*rhor[4]*uvar2r[4]+0.7905694150420948*rhol[4]*uvar2l[4]+2.371708245126285*rhor[3]*uvar2r[3]-1.369306393762915*rhor[0]*uvar2r[3]+2.371708245126285*rhol[3]*uvar2l[3]+1.369306393762915*rhol[0]*uvar2l[3]-1.369306393762915*uvar2r[0]*rhor[3]+1.369306393762915*uvar2l[0]*rhol[3]+0.7905694150420948*rhor[2]*uvar2r[2]+0.7905694150420948*rhol[2]*uvar2l[2]+0.7905694150420948*rhor[1]*uvar2r[1]+0.7905694150420948*rhol[1]*uvar2l[1]+0.7905694150420948*rhor[0]*uvar2r[0]+0.7905694150420948*rhol[0]*uvar2l[0]; 

  outrrhouz[0] += incr[0]*dxr1; 
  outrrhouz[1] += incr[1]*dxr1; 
  outrrhouz[2] += incr[2]*dxr1; 
  outrrhouz[3] += incr[3]*dxr1; 
  outrrhouz[4] += incr[4]*dxr1; 
  outrrhouz[5] += incr[5]*dxr1; 
  outrrhouz[6] += incr[6]*dxr1; 
  outrrhouz[7] += incr[7]*dxr1; 
  outrrhouz[8] += incr[8]*dxr1; 
  outrrhouz[9] += incr[9]*dxr1; 
  outrrhouz[10] += incr[10]*dxr1; 
  outrrhouz[11] += incr[11]*dxr1; 
  outrrhouz[12] += incr[12]*dxr1; 
  outrrhouz[13] += incr[13]*dxr1; 
  outrrhouz[14] += incr[14]*dxr1; 
  outrrhouz[15] += incr[15]*dxr1; 
  outrrhouz[16] += incr[16]*dxr1; 
  outrrhouz[17] += incr[17]*dxr1; 
  outrrhouz[18] += incr[18]*dxr1; 
  outrrhouz[19] += incr[19]*dxr1; 

  outlrhouz[0] += -1.0*incr[0]*dxl1; 
  outlrhouz[3] += incr[3]*dxl1; 
  outlrhouz[9] += -1.0*incr[9]*dxl1; 

 
} 
