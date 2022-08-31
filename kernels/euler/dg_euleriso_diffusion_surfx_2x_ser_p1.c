#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfx_2x_ser_p1(const double* w, const double* dx,
  const double* D_in,
  const double* uvarl, const double* uvarc, const double* uvarr,
  const double* statevecl, const double* statevecc, const double* statevecr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // uvarl: Input velocity in the left cell
  // uvarc: Input velocity in the center cell
  // uvarr: Input velocity in the right cell
  // statevecl: Input field in the left cell
  // statevecc: Input field in the center cell
  // statevecr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[0]/dx[0];
  double dxsqrd = dx[0]*dx[0];
  double dysqrd = dx[1]*dx[1];
  const double *D = &D_in[4]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[4]; 
  const double *uvarzl = &uvarl[8]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[4]; 
  const double *uvarzc = &uvarc[8]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[4]; 
  const double *uvarzr = &uvarr[8]; 
  out[0] += J*((-(0.875*D[3]*uvaryr[4])/dysqrd)-(1.082531754730548*D[1]*uvaryr[4])/dysqrd-(0.875*D[3]*uvaryl[4])/dysqrd+(1.082531754730548*D[1]*uvaryl[4])/dysqrd-(5.75*D[3]*uvaryc[4])/dysqrd-(0.875*D[2]*uvaryr[3])/dysqrd-(1.082531754730548*D[0]*uvaryr[3])/dysqrd-(0.875*D[2]*uvaryl[3])/dysqrd+(1.082531754730548*D[0]*uvaryl[3])/dysqrd-(5.75*D[2]*uvaryc[3])/dysqrd+(1.082531754730548*uvaryr[2]*D[3])/dysqrd-(1.082531754730548*uvaryl[2]*D[3])/dysqrd+(1.125*D[1]*uvaryr[2])/dysqrd+(1.125*D[1]*uvaryl[2])/dysqrd-(2.25*D[1]*uvaryc[2])/dysqrd+(1.082531754730548*uvaryr[1]*D[2])/dysqrd-(1.082531754730548*uvaryl[1]*D[2])/dysqrd+(1.125*D[0]*uvaryr[1])/dysqrd+(1.125*D[0]*uvaryl[1])/dysqrd-(2.25*D[0]*uvaryc[1])/dysqrd-(1.237436867076458*D[1]*uvarxr[2])/dxsqrd-(1.530931089239486*D[0]*uvarxr[2])/dxsqrd-(1.237436867076458*D[1]*uvarxl[2])/dxsqrd+(1.530931089239486*D[0]*uvarxl[2])/dxsqrd-(8.131727983645295*D[1]*uvarxc[2])/dxsqrd+(1.530931089239486*D[1]*uvarxr[1])/dxsqrd+(1.590990257669731*D[0]*uvarxr[1])/dxsqrd-(1.530931089239486*D[1]*uvarxl[1])/dxsqrd+(1.590990257669731*D[0]*uvarxl[1])/dxsqrd-(3.181980515339463*D[0]*uvarxc[1])/dxsqrd);
  out[1] += J*((-(0.875*D[2]*uvaryr[4])/dysqrd)-(1.082531754730548*D[0]*uvaryr[4])/dysqrd-(0.875*D[2]*uvaryl[4])/dysqrd+(1.082531754730548*D[0]*uvaryl[4])/dysqrd-(5.75*D[2]*uvaryc[4])/dysqrd-(0.875*D[3]*uvaryr[3])/dysqrd-(1.082531754730548*D[1]*uvaryr[3])/dysqrd-(0.875*D[3]*uvaryl[3])/dysqrd+(1.082531754730548*D[1]*uvaryl[3])/dysqrd-(5.75*D[3]*uvaryc[3])/dysqrd+(1.082531754730548*uvaryr[1]*D[3])/dysqrd-(1.082531754730548*uvaryl[1]*D[3])/dysqrd+(1.082531754730548*D[2]*uvaryr[2])/dysqrd+(1.125*D[0]*uvaryr[2])/dysqrd-(1.082531754730548*D[2]*uvaryl[2])/dysqrd+(1.125*D[0]*uvaryl[2])/dysqrd-(2.25*D[0]*uvaryc[2])/dysqrd+(1.125*D[1]*uvaryr[1])/dysqrd+(1.125*D[1]*uvaryl[1])/dysqrd-(2.25*D[1]*uvaryc[1])/dysqrd+(0.3061862178478971*D[1]*uvarxr[2])/dxsqrd-(1.237436867076458*D[0]*uvarxr[2])/dxsqrd-(0.3061862178478971*D[1]*uvarxl[2])/dxsqrd-(1.237436867076458*D[0]*uvarxl[2])/dxsqrd-(8.131727983645295*D[0]*uvarxc[2])/dxsqrd+(0.5303300858899105*D[1]*uvarxr[1])/dxsqrd+(1.530931089239486*D[0]*uvarxr[1])/dxsqrd+(0.5303300858899105*D[1]*uvarxl[1])/dxsqrd-(1.530931089239486*D[0]*uvarxl[1])/dxsqrd-(1.060660171779821*D[1]*uvarxc[1])/dxsqrd);
  out[2] += J*((0.2165063509461096*D[3]*uvaryr[4])/dysqrd-(0.875*D[1]*uvaryr[4])/dysqrd-(0.2165063509461096*D[3]*uvaryl[4])/dysqrd-(0.875*D[1]*uvaryl[4])/dysqrd-(5.75*D[1]*uvaryc[4])/dysqrd+(0.2165063509461096*D[2]*uvaryr[3])/dysqrd-(0.875*D[0]*uvaryr[3])/dysqrd-(0.2165063509461096*D[2]*uvaryl[3])/dysqrd-(0.875*D[0]*uvaryl[3])/dysqrd-(5.75*D[0]*uvaryc[3])/dysqrd+(0.375*uvaryr[2]*D[3])/dysqrd+(0.375*uvaryl[2]*D[3])/dysqrd-(0.75*uvaryc[2]*D[3])/dysqrd+(1.082531754730548*D[1]*uvaryr[2])/dysqrd-(1.082531754730548*D[1]*uvaryl[2])/dysqrd+(0.375*uvaryr[1]*D[2])/dysqrd+(0.375*uvaryl[1]*D[2])/dysqrd-(0.75*uvaryc[1]*D[2])/dysqrd+(1.082531754730548*D[0]*uvaryr[1])/dysqrd-(1.082531754730548*D[0]*uvaryl[1])/dysqrd-(1.237436867076458*uvarxr[2]*D[3])/dxsqrd-(1.237436867076458*uvarxl[2]*D[3])/dxsqrd-(8.131727983645295*uvarxc[2]*D[3])/dxsqrd+(1.530931089239486*uvarxr[1]*D[3])/dxsqrd-(1.530931089239486*uvarxl[1]*D[3])/dxsqrd-(1.530931089239486*D[2]*uvarxr[2])/dxsqrd+(1.530931089239486*D[2]*uvarxl[2])/dxsqrd+(1.590990257669731*uvarxr[1]*D[2])/dxsqrd+(1.590990257669731*uvarxl[1]*D[2])/dxsqrd-(3.181980515339463*uvarxc[1]*D[2])/dxsqrd);
  out[3] += J*((0.2165063509461096*D[2]*uvaryr[4])/dysqrd-(0.875*D[0]*uvaryr[4])/dysqrd-(0.2165063509461096*D[2]*uvaryl[4])/dysqrd-(0.875*D[0]*uvaryl[4])/dysqrd-(5.75*D[0]*uvaryc[4])/dysqrd+(0.2165063509461096*D[3]*uvaryr[3])/dysqrd-(0.875*D[1]*uvaryr[3])/dysqrd-(0.2165063509461096*D[3]*uvaryl[3])/dysqrd-(0.875*D[1]*uvaryl[3])/dysqrd-(5.75*D[1]*uvaryc[3])/dysqrd+(0.375*uvaryr[1]*D[3])/dysqrd+(0.375*uvaryl[1]*D[3])/dysqrd-(0.75*uvaryc[1]*D[3])/dysqrd+(0.375*D[2]*uvaryr[2])/dysqrd+(1.082531754730548*D[0]*uvaryr[2])/dysqrd+(0.375*D[2]*uvaryl[2])/dysqrd-(1.082531754730548*D[0]*uvaryl[2])/dysqrd-(0.75*D[2]*uvaryc[2])/dysqrd+(1.082531754730548*D[1]*uvaryr[1])/dysqrd-(1.082531754730548*D[1]*uvaryl[1])/dysqrd+(0.3061862178478971*uvarxr[2]*D[3])/dxsqrd-(0.3061862178478971*uvarxl[2]*D[3])/dxsqrd+(0.5303300858899105*uvarxr[1]*D[3])/dxsqrd+(0.5303300858899105*uvarxl[1]*D[3])/dxsqrd-(1.060660171779821*uvarxc[1]*D[3])/dxsqrd-(1.237436867076458*D[2]*uvarxr[2])/dxsqrd-(1.237436867076458*D[2]*uvarxl[2])/dxsqrd-(8.131727983645295*D[2]*uvarxc[2])/dxsqrd+(1.530931089239486*uvarxr[1]*D[2])/dxsqrd-(1.530931089239486*uvarxl[1]*D[2])/dxsqrd);
}
