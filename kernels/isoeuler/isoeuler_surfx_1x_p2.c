#include <gkyl_isoeuler_kernels.h> 
GKYL_CU_DH void isoeuler_surfx_1x_ser_p2(const double *w, const double *dxv, double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out)  
{ 
  // w: Cell-center coordinates. dxv[NDIM]: Cell spacing. statevec(l/c/r): [rho, rho ux, rho uy, rho uz] in (left/center/right) cell, uvar(l/c/r): [ux, uv, uz]  in (left/center/right) cell
  // out: output in left/right cells 
  const double dxl1 = 2.0/dxv[0]; 
  const double dxr1 = 2.0/dxv[0]; 
  const double *rhol = &statevecl[0]; 
  const double *rhou0l = &statevecl[3]; 
  const double *rhou1l = &statevecl[6]; 
  const double *rhou2l = &statevecl[9]; 
  const double *rhoc = &statevecc[0]; 
  const double *rhou0c = &statevecc[3]; 
  const double *rhou1c = &statevecc[6]; 
  const double *rhou2c = &statevecc[9]; 
  const double *rhor = &statevecr[0]; 
  const double *rhou0r = &statevecr[3]; 
  const double *rhou1r = &statevecr[6]; 
  const double *rhou2r = &statevecr[9]; 
  const double *uvar0l = &uvarl[0]; 
  const double *uvar1l = &uvarl[3]; 
  const double *uvar2l = &uvarl[6]; 
  const double *uvar0c = &uvarc[0]; 
  const double *uvar1c = &uvarc[3]; 
  const double *uvar2c = &uvarc[6]; 
  const double *uvar0r = &uvarr[0]; 
  const double *uvar1r = &uvarr[3]; 
  const double *uvar2r = &uvarr[6]; 
  double *outrho = &out[0]; 
  double *outrhoux = &out[3]; 
  double *outrhouy = &out[6]; 
  double *outrhouz = &out[9]; 
  double incr[3]; 

  double vthsq = vth*vth; 
  //FluxRho; 
  outrho[0] += (-0.5590169943749475*rhou0r[2]*dxr1)-0.5590169943749475*rhou0c[2]*dxr1+0.4330127018922193*rhou0r[1]*dxr1-0.4330127018922193*rhou0c[1]*dxr1-0.25*rhou0r[0]*dxr1-0.25*rhou0c[0]*dxr1+0.5590169943749475*rhou0l[2]*dxl1+0.5590169943749475*rhou0c[2]*dxl1+0.4330127018922193*rhou0l[1]*dxl1-0.4330127018922193*rhou0c[1]*dxl1+0.25*rhou0l[0]*dxl1+0.25*rhou0c[0]*dxl1; 
  outrho[1] += (-0.9682458365518543*rhou0r[2]*dxr1)-0.9682458365518543*rhou0c[2]*dxr1+0.75*rhou0r[1]*dxr1-0.75*rhou0c[1]*dxr1-0.4330127018922193*rhou0r[0]*dxr1-0.4330127018922193*rhou0c[0]*dxr1-0.9682458365518543*rhou0l[2]*dxl1-0.9682458365518543*rhou0c[2]*dxl1-0.75*rhou0l[1]*dxl1+0.75*rhou0c[1]*dxl1-0.4330127018922193*rhou0l[0]*dxl1-0.4330127018922193*rhou0c[0]*dxl1; 
  outrho[2] += (-1.25*rhou0r[2]*dxr1)-1.25*rhou0c[2]*dxr1+0.9682458365518543*rhou0r[1]*dxr1-0.9682458365518543*rhou0c[1]*dxr1-0.5590169943749475*rhou0r[0]*dxr1-0.5590169943749475*rhou0c[0]*dxr1+1.25*rhou0l[2]*dxl1+1.25*rhou0c[2]*dxl1+0.9682458365518543*rhou0l[1]*dxl1-0.9682458365518543*rhou0c[1]*dxl1+0.5590169943749475*rhou0l[0]*dxl1+0.5590169943749475*rhou0c[0]*dxl1; 

 
  //FluxRhoUx; 
  outrhoux[0] += 1.185854122563142*rhor[2]*dxr1*vthsq-1.185854122563142*rhoc[2]*dxr1*vthsq-0.3061862178478972*rhor[1]*dxr1*vthsq-0.3061862178478972*rhoc[1]*dxr1*vthsq+1.185854122563142*rhol[2]*dxl1*vthsq-1.185854122563142*rhoc[2]*dxl1*vthsq+0.3061862178478972*rhol[1]*dxl1*vthsq+0.3061862178478972*rhoc[1]*dxl1*vthsq-0.8838834764831843*rhor[2]*uvar0r[2]*dxr1+0.6846531968814576*rhor[1]*uvar0r[2]*dxr1-0.3952847075210474*rhor[0]*uvar0r[2]*dxr1-0.8838834764831843*rhoc[2]*uvar0c[2]*dxr1-0.6846531968814576*rhoc[1]*uvar0c[2]*dxr1-0.3952847075210474*rhoc[0]*uvar0c[2]*dxr1+0.6846531968814576*uvar0r[1]*rhor[2]*dxr1-0.3952847075210474*uvar0r[0]*rhor[2]*dxr1-0.6846531968814576*uvar0c[1]*rhoc[2]*dxr1-0.3952847075210474*uvar0c[0]*rhoc[2]*dxr1-0.5303300858899106*rhor[1]*uvar0r[1]*dxr1+0.3061862178478972*rhor[0]*uvar0r[1]*dxr1-0.5303300858899106*rhoc[1]*uvar0c[1]*dxr1-0.3061862178478972*rhoc[0]*uvar0c[1]*dxr1+0.3061862178478972*uvar0r[0]*rhor[1]*dxr1-0.3061862178478972*uvar0c[0]*rhoc[1]*dxr1-0.1767766952966369*rhor[0]*uvar0r[0]*dxr1-0.1767766952966369*rhoc[0]*uvar0c[0]*dxr1+0.8838834764831843*rhol[2]*uvar0l[2]*dxl1+0.6846531968814576*rhol[1]*uvar0l[2]*dxl1+0.3952847075210474*rhol[0]*uvar0l[2]*dxl1+0.8838834764831843*rhoc[2]*uvar0c[2]*dxl1-0.6846531968814576*rhoc[1]*uvar0c[2]*dxl1+0.3952847075210474*rhoc[0]*uvar0c[2]*dxl1+0.6846531968814576*uvar0l[1]*rhol[2]*dxl1+0.3952847075210474*uvar0l[0]*rhol[2]*dxl1-0.6846531968814576*uvar0c[1]*rhoc[2]*dxl1+0.3952847075210474*uvar0c[0]*rhoc[2]*dxl1+0.5303300858899106*rhol[1]*uvar0l[1]*dxl1+0.3061862178478972*rhol[0]*uvar0l[1]*dxl1+0.5303300858899106*rhoc[1]*uvar0c[1]*dxl1-0.3061862178478972*rhoc[0]*uvar0c[1]*dxl1+0.3061862178478972*uvar0l[0]*rhol[1]*dxl1-0.3061862178478972*uvar0c[0]*rhoc[1]*dxl1+0.1767766952966369*rhol[0]*uvar0l[0]*dxl1+0.1767766952966369*rhoc[0]*uvar0c[0]*dxl1; 
  outrhoux[1] += (-3.557562367689427*rhor[2]*dxr1*vthsq)-3.557562367689427*rhoc[2]*dxr1*vthsq+0.9185586535436916*rhor[1]*dxr1*vthsq-0.9185586535436916*rhoc[1]*dxr1*vthsq-3.557562367689427*rhol[2]*dxl1*vthsq-3.557562367689427*rhoc[2]*dxl1*vthsq-0.9185586535436916*rhol[1]*dxl1*vthsq+0.9185586535436916*rhoc[1]*dxl1*vthsq-1.530931089239486*rhor[2]*uvar0r[2]*dxr1+1.185854122563142*rhor[1]*uvar0r[2]*dxr1-0.6846531968814576*rhor[0]*uvar0r[2]*dxr1-1.530931089239486*rhoc[2]*uvar0c[2]*dxr1-1.185854122563142*rhoc[1]*uvar0c[2]*dxr1-0.6846531968814576*rhoc[0]*uvar0c[2]*dxr1+1.185854122563142*uvar0r[1]*rhor[2]*dxr1-0.6846531968814576*uvar0r[0]*rhor[2]*dxr1-1.185854122563142*uvar0c[1]*rhoc[2]*dxr1-0.6846531968814576*uvar0c[0]*rhoc[2]*dxr1-0.9185586535436916*rhor[1]*uvar0r[1]*dxr1+0.5303300858899106*rhor[0]*uvar0r[1]*dxr1-0.9185586535436916*rhoc[1]*uvar0c[1]*dxr1-0.5303300858899106*rhoc[0]*uvar0c[1]*dxr1+0.5303300858899106*uvar0r[0]*rhor[1]*dxr1-0.5303300858899106*uvar0c[0]*rhoc[1]*dxr1-0.3061862178478972*rhor[0]*uvar0r[0]*dxr1-0.3061862178478972*rhoc[0]*uvar0c[0]*dxr1-1.530931089239486*rhol[2]*uvar0l[2]*dxl1-1.185854122563142*rhol[1]*uvar0l[2]*dxl1-0.6846531968814576*rhol[0]*uvar0l[2]*dxl1-1.530931089239486*rhoc[2]*uvar0c[2]*dxl1+1.185854122563142*rhoc[1]*uvar0c[2]*dxl1-0.6846531968814576*rhoc[0]*uvar0c[2]*dxl1-1.185854122563142*uvar0l[1]*rhol[2]*dxl1-0.6846531968814576*uvar0l[0]*rhol[2]*dxl1+1.185854122563142*uvar0c[1]*rhoc[2]*dxl1-0.6846531968814576*uvar0c[0]*rhoc[2]*dxl1-0.9185586535436916*rhol[1]*uvar0l[1]*dxl1-0.5303300858899106*rhol[0]*uvar0l[1]*dxl1-0.9185586535436916*rhoc[1]*uvar0c[1]*dxl1+0.5303300858899106*rhoc[0]*uvar0c[1]*dxl1-0.5303300858899106*uvar0l[0]*rhol[1]*dxl1+0.5303300858899106*uvar0c[0]*rhoc[1]*dxl1-0.3061862178478972*rhol[0]*uvar0l[0]*dxl1-0.3061862178478972*rhoc[0]*uvar0c[0]*dxl1; 
  outrhoux[2] += 5.929270612815712*rhor[2]*dxr1*vthsq-5.929270612815712*rhoc[2]*dxr1*vthsq-1.530931089239486*rhor[1]*dxr1*vthsq-1.530931089239486*rhoc[1]*dxr1*vthsq+5.929270612815712*rhol[2]*dxl1*vthsq-5.929270612815712*rhoc[2]*dxl1*vthsq+1.530931089239486*rhol[1]*dxl1*vthsq+1.530931089239486*rhoc[1]*dxl1*vthsq-1.976423537605237*rhor[2]*uvar0r[2]*dxr1+1.530931089239486*rhor[1]*uvar0r[2]*dxr1-0.8838834764831843*rhor[0]*uvar0r[2]*dxr1-1.976423537605237*rhoc[2]*uvar0c[2]*dxr1-1.530931089239486*rhoc[1]*uvar0c[2]*dxr1-0.8838834764831843*rhoc[0]*uvar0c[2]*dxr1+1.530931089239486*uvar0r[1]*rhor[2]*dxr1-0.8838834764831843*uvar0r[0]*rhor[2]*dxr1-1.530931089239486*uvar0c[1]*rhoc[2]*dxr1-0.8838834764831843*uvar0c[0]*rhoc[2]*dxr1-1.185854122563142*rhor[1]*uvar0r[1]*dxr1+0.6846531968814576*rhor[0]*uvar0r[1]*dxr1-1.185854122563142*rhoc[1]*uvar0c[1]*dxr1-0.6846531968814576*rhoc[0]*uvar0c[1]*dxr1+0.6846531968814576*uvar0r[0]*rhor[1]*dxr1-0.6846531968814576*uvar0c[0]*rhoc[1]*dxr1-0.3952847075210474*rhor[0]*uvar0r[0]*dxr1-0.3952847075210474*rhoc[0]*uvar0c[0]*dxr1+1.976423537605237*rhol[2]*uvar0l[2]*dxl1+1.530931089239486*rhol[1]*uvar0l[2]*dxl1+0.8838834764831843*rhol[0]*uvar0l[2]*dxl1+1.976423537605237*rhoc[2]*uvar0c[2]*dxl1-1.530931089239486*rhoc[1]*uvar0c[2]*dxl1+0.8838834764831843*rhoc[0]*uvar0c[2]*dxl1+1.530931089239486*uvar0l[1]*rhol[2]*dxl1+0.8838834764831843*uvar0l[0]*rhol[2]*dxl1-1.530931089239486*uvar0c[1]*rhoc[2]*dxl1+0.8838834764831843*uvar0c[0]*rhoc[2]*dxl1+1.185854122563142*rhol[1]*uvar0l[1]*dxl1+0.6846531968814576*rhol[0]*uvar0l[1]*dxl1+1.185854122563142*rhoc[1]*uvar0c[1]*dxl1-0.6846531968814576*rhoc[0]*uvar0c[1]*dxl1+0.6846531968814576*uvar0l[0]*rhol[1]*dxl1-0.6846531968814576*uvar0c[0]*rhoc[1]*dxl1+0.3952847075210474*rhol[0]*uvar0l[0]*dxl1+0.3952847075210474*rhoc[0]*uvar0c[0]*dxl1; 

 
} 
