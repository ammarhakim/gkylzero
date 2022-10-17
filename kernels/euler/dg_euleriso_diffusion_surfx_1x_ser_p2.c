#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfx_1x_ser_p2(const double* w, const double* dx,
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
  const double *D = &D_in[3]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[3]; 
  const double *uvarzl = &uvarl[6]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[3]; 
  const double *uvarzc = &uvarc[6]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[3]; 
  const double *uvarzr = &uvarr[6]; 
  out[3] += J*((0.6708203932499369*uvarxr[2]*mu)/dxsqrd+(0.6708203932499369*uvarxl[2]*mu)/dxsqrd-(1.341640786499874*uvarxc[2]*mu)/dxsqrd-(1.190784930203603*uvarxr[1]*mu)/dxsqrd+(1.190784930203603*uvarxl[1]*mu)/dxsqrd+(0.9375*uvarxr[0]*mu)/dxsqrd+(0.9375*uvarxl[0]*mu)/dxsqrd-(1.875*uvarxc[0]*mu)/dxsqrd);
  out[4] += J*((0.7382874503707888*uvarxr[2]*mu)/dxsqrd-(0.7382874503707888*uvarxl[2]*mu)/dxsqrd-(1.453125*uvarxr[1]*mu)/dxsqrd-(1.453125*uvarxl[1]*mu)/dxsqrd-(5.34375*uvarxc[1]*mu)/dxsqrd+(1.190784930203603*uvarxr[0]*mu)/dxsqrd-(1.190784930203603*uvarxl[0]*mu)/dxsqrd);
  out[5] += J*((-(0.140625*uvarxr[2]*mu)/dxsqrd)-(0.140625*uvarxl[2]*mu)/dxsqrd-(6.28125*uvarxc[2]*mu)/dxsqrd-(0.3025768239224545*uvarxr[1]*mu)/dxsqrd+(0.3025768239224545*uvarxl[1]*mu)/dxsqrd+(0.4192627457812106*uvarxr[0]*mu)/dxsqrd+(0.4192627457812106*uvarxl[0]*mu)/dxsqrd-(0.8385254915624212*uvarxc[0]*mu)/dxsqrd);
}
