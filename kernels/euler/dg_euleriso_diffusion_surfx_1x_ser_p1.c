#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfx_1x_ser_p1(const double* w, const double* dx,
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
  const double *D = &D_in[2]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[2]; 
  const double *uvarzl = &uvarl[4]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[2]; 
  const double *uvarzc = &uvarc[4]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[2]; 
  const double *uvarzr = &uvarr[4]; 
  out[0] += J*((-(1.237436867076458*D[1]*uvarxr[2])/dxsqrd)-(1.530931089239486*D[0]*uvarxr[2])/dxsqrd-(1.237436867076458*D[1]*uvarxl[2])/dxsqrd+(1.530931089239486*D[0]*uvarxl[2])/dxsqrd-(8.131727983645295*D[1]*uvarxc[2])/dxsqrd+(1.530931089239486*D[1]*uvarxr[1])/dxsqrd+(1.590990257669731*D[0]*uvarxr[1])/dxsqrd-(1.530931089239486*D[1]*uvarxl[1])/dxsqrd+(1.590990257669731*D[0]*uvarxl[1])/dxsqrd-(3.181980515339463*D[0]*uvarxc[1])/dxsqrd);
  out[1] += J*((0.3061862178478971*D[1]*uvarxr[2])/dxsqrd-(1.237436867076458*D[0]*uvarxr[2])/dxsqrd-(0.3061862178478971*D[1]*uvarxl[2])/dxsqrd-(1.237436867076458*D[0]*uvarxl[2])/dxsqrd-(8.131727983645295*D[0]*uvarxc[2])/dxsqrd+(0.5303300858899105*D[1]*uvarxr[1])/dxsqrd+(1.530931089239486*D[0]*uvarxr[1])/dxsqrd+(0.5303300858899105*D[1]*uvarxl[1])/dxsqrd-(1.530931089239486*D[0]*uvarxl[1])/dxsqrd-(1.060660171779821*D[1]*uvarxc[1])/dxsqrd);
}
