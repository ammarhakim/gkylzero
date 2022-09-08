#include <dg_euleriso_diffusion_kernels.h>

GKYL_CU_DH void
dg_euleriso_diffusion_surfx_2x_ser_p2(const double* w, const double* dx,
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
  const double *D = &D_in[8]; 
  double mu = D[0]; 
  const double *uvarxl = &uvarl[0]; 
  const double *uvaryl = &uvarl[8]; 
  const double *uvarzl = &uvarl[16]; 
  const double *uvarxc = &uvarc[0]; 
  const double *uvaryc = &uvarc[8]; 
  const double *uvarzc = &uvarc[16]; 
  const double *uvarxr = &uvarr[0]; 
  const double *uvaryr = &uvarr[8]; 
  const double *uvarzr = &uvarr[16]; 
  out[8] += J*((0.6708203932499369*uvaryr[5]*mu)/dysqrd+(0.6708203932499369*uvaryl[5]*mu)/dysqrd-(1.341640786499874*uvaryc[5]*mu)/dysqrd-(1.190784930203603*uvaryr[2]*mu)/dysqrd+(1.190784930203603*uvaryl[2]*mu)/dysqrd+(0.9375*uvaryr[0]*mu)/dysqrd+(0.9375*uvaryl[0]*mu)/dysqrd-(1.875*uvaryc[0]*mu)/dysqrd+(0.6708203932499369*uvarxr[4]*mu)/dxsqrd+(0.6708203932499369*uvarxl[4]*mu)/dxsqrd-(1.341640786499874*uvarxc[4]*mu)/dxsqrd-(1.190784930203603*uvarxr[1]*mu)/dxsqrd+(1.190784930203603*uvarxl[1]*mu)/dxsqrd+(0.9375*uvarxr[0]*mu)/dxsqrd+(0.9375*uvarxl[0]*mu)/dxsqrd-(1.875*uvarxc[0]*mu)/dxsqrd);
  out[9] += J*((0.6708203932499369*uvaryr[7]*mu)/dysqrd+(0.6708203932499369*uvaryl[7]*mu)/dysqrd-(1.341640786499874*uvaryc[7]*mu)/dysqrd-(1.190784930203603*uvaryr[3]*mu)/dysqrd+(1.190784930203603*uvaryl[3]*mu)/dysqrd+(0.9375*uvaryr[1]*mu)/dysqrd+(0.9375*uvaryl[1]*mu)/dysqrd-(1.875*uvaryc[1]*mu)/dysqrd+(0.7382874503707888*uvarxr[4]*mu)/dxsqrd-(0.7382874503707888*uvarxl[4]*mu)/dxsqrd-(1.453125*uvarxr[1]*mu)/dxsqrd-(1.453125*uvarxl[1]*mu)/dxsqrd-(5.34375*uvarxc[1]*mu)/dxsqrd+(1.190784930203603*uvarxr[0]*mu)/dxsqrd-(1.190784930203603*uvarxl[0]*mu)/dxsqrd);
  out[10] += J*((0.7382874503707888*uvaryr[5]*mu)/dysqrd-(0.7382874503707888*uvaryl[5]*mu)/dysqrd-(1.453125*uvaryr[2]*mu)/dysqrd-(1.453125*uvaryl[2]*mu)/dysqrd-(5.34375*uvaryc[2]*mu)/dysqrd+(1.190784930203603*uvaryr[0]*mu)/dysqrd-(1.190784930203603*uvaryl[0]*mu)/dysqrd+(0.6708203932499369*uvarxr[6]*mu)/dxsqrd+(0.6708203932499369*uvarxl[6]*mu)/dxsqrd-(1.341640786499874*uvarxc[6]*mu)/dxsqrd-(1.190784930203603*uvarxr[3]*mu)/dxsqrd+(1.190784930203603*uvarxl[3]*mu)/dxsqrd+(0.9375*uvarxr[2]*mu)/dxsqrd+(0.9375*uvarxl[2]*mu)/dxsqrd-(1.875*uvarxc[2]*mu)/dxsqrd);
  out[11] += J*((0.7382874503707888*uvaryr[7]*mu)/dysqrd-(0.7382874503707888*uvaryl[7]*mu)/dysqrd-(1.453125*uvaryr[3]*mu)/dysqrd-(1.453125*uvaryl[3]*mu)/dysqrd-(5.34375*uvaryc[3]*mu)/dysqrd+(1.190784930203603*uvaryr[1]*mu)/dysqrd-(1.190784930203603*uvaryl[1]*mu)/dysqrd+(0.7382874503707888*uvarxr[6]*mu)/dxsqrd-(0.7382874503707888*uvarxl[6]*mu)/dxsqrd-(1.453125*uvarxr[3]*mu)/dxsqrd-(1.453125*uvarxl[3]*mu)/dxsqrd-(5.34375*uvarxc[3]*mu)/dxsqrd+(1.190784930203603*uvarxr[2]*mu)/dxsqrd-(1.190784930203603*uvarxl[2]*mu)/dxsqrd);
  out[12] += J*((-(1.190784930203603*uvaryr[6]*mu)/dysqrd)+(1.190784930203603*uvaryl[6]*mu)/dysqrd+(0.9375*uvaryr[4]*mu)/dysqrd+(0.9375*uvaryl[4]*mu)/dysqrd-(1.875*uvaryc[4]*mu)/dysqrd-(0.140625*uvarxr[4]*mu)/dxsqrd-(0.140625*uvarxl[4]*mu)/dxsqrd-(6.28125*uvarxc[4]*mu)/dxsqrd-(0.3025768239224545*uvarxr[1]*mu)/dxsqrd+(0.3025768239224545*uvarxl[1]*mu)/dxsqrd+(0.4192627457812106*uvarxr[0]*mu)/dxsqrd+(0.4192627457812106*uvarxl[0]*mu)/dxsqrd-(0.8385254915624212*uvarxc[0]*mu)/dxsqrd);
  out[13] += J*((-(0.140625*uvaryr[5]*mu)/dysqrd)-(0.140625*uvaryl[5]*mu)/dysqrd-(6.28125*uvaryc[5]*mu)/dysqrd-(0.3025768239224545*uvaryr[2]*mu)/dysqrd+(0.3025768239224545*uvaryl[2]*mu)/dysqrd+(0.4192627457812106*uvaryr[0]*mu)/dysqrd+(0.4192627457812106*uvaryl[0]*mu)/dysqrd-(0.8385254915624212*uvaryc[0]*mu)/dysqrd-(1.190784930203603*uvarxr[7]*mu)/dxsqrd+(1.190784930203603*uvarxl[7]*mu)/dxsqrd+(0.9375*uvarxr[5]*mu)/dxsqrd+(0.9375*uvarxl[5]*mu)/dxsqrd-(1.875*uvarxc[5]*mu)/dxsqrd);
  out[14] += J*((-(1.453125*uvaryr[6]*mu)/dysqrd)-(1.453125*uvaryl[6]*mu)/dysqrd-(5.34375*uvaryc[6]*mu)/dysqrd+(1.190784930203603*uvaryr[4]*mu)/dysqrd-(1.190784930203603*uvaryl[4]*mu)/dysqrd-(0.140625*uvarxr[6]*mu)/dxsqrd-(0.140625*uvarxl[6]*mu)/dxsqrd-(6.28125*uvarxc[6]*mu)/dxsqrd-(0.3025768239224544*uvarxr[3]*mu)/dxsqrd+(0.3025768239224544*uvarxl[3]*mu)/dxsqrd+(0.4192627457812105*uvarxr[2]*mu)/dxsqrd+(0.4192627457812105*uvarxl[2]*mu)/dxsqrd-(0.8385254915624211*uvarxc[2]*mu)/dxsqrd);
  out[15] += J*((-(0.140625*uvaryr[7]*mu)/dysqrd)-(0.140625*uvaryl[7]*mu)/dysqrd-(6.28125*uvaryc[7]*mu)/dysqrd-(0.3025768239224544*uvaryr[3]*mu)/dysqrd+(0.3025768239224544*uvaryl[3]*mu)/dysqrd+(0.4192627457812105*uvaryr[1]*mu)/dysqrd+(0.4192627457812105*uvaryl[1]*mu)/dysqrd-(0.8385254915624211*uvaryc[1]*mu)/dysqrd-(1.453125*uvarxr[7]*mu)/dxsqrd-(1.453125*uvarxl[7]*mu)/dxsqrd-(5.34375*uvarxc[7]*mu)/dxsqrd+(1.190784930203603*uvarxr[5]*mu)/dxsqrd-(1.190784930203603*uvarxl[5]*mu)/dxsqrd);
  out[16] += J*((0.6708203932499369*uvaryr[5]*mu)/dysqrd+(0.6708203932499369*uvaryl[5]*mu)/dysqrd-(1.341640786499874*uvaryc[5]*mu)/dysqrd-(1.190784930203603*uvaryr[2]*mu)/dysqrd+(1.190784930203603*uvaryl[2]*mu)/dysqrd+(0.9375*uvaryr[0]*mu)/dysqrd+(0.9375*uvaryl[0]*mu)/dysqrd-(1.875*uvaryc[0]*mu)/dysqrd+(0.6708203932499369*uvarxr[4]*mu)/dxsqrd+(0.6708203932499369*uvarxl[4]*mu)/dxsqrd-(1.341640786499874*uvarxc[4]*mu)/dxsqrd-(1.190784930203603*uvarxr[1]*mu)/dxsqrd+(1.190784930203603*uvarxl[1]*mu)/dxsqrd+(0.9375*uvarxr[0]*mu)/dxsqrd+(0.9375*uvarxl[0]*mu)/dxsqrd-(1.875*uvarxc[0]*mu)/dxsqrd);
  out[17] += J*((0.6708203932499369*uvaryr[7]*mu)/dysqrd+(0.6708203932499369*uvaryl[7]*mu)/dysqrd-(1.341640786499874*uvaryc[7]*mu)/dysqrd-(1.190784930203603*uvaryr[3]*mu)/dysqrd+(1.190784930203603*uvaryl[3]*mu)/dysqrd+(0.9375*uvaryr[1]*mu)/dysqrd+(0.9375*uvaryl[1]*mu)/dysqrd-(1.875*uvaryc[1]*mu)/dysqrd+(0.7382874503707888*uvarxr[4]*mu)/dxsqrd-(0.7382874503707888*uvarxl[4]*mu)/dxsqrd-(1.453125*uvarxr[1]*mu)/dxsqrd-(1.453125*uvarxl[1]*mu)/dxsqrd-(5.34375*uvarxc[1]*mu)/dxsqrd+(1.190784930203603*uvarxr[0]*mu)/dxsqrd-(1.190784930203603*uvarxl[0]*mu)/dxsqrd);
  out[18] += J*((0.7382874503707888*uvaryr[5]*mu)/dysqrd-(0.7382874503707888*uvaryl[5]*mu)/dysqrd-(1.453125*uvaryr[2]*mu)/dysqrd-(1.453125*uvaryl[2]*mu)/dysqrd-(5.34375*uvaryc[2]*mu)/dysqrd+(1.190784930203603*uvaryr[0]*mu)/dysqrd-(1.190784930203603*uvaryl[0]*mu)/dysqrd+(0.6708203932499369*uvarxr[6]*mu)/dxsqrd+(0.6708203932499369*uvarxl[6]*mu)/dxsqrd-(1.341640786499874*uvarxc[6]*mu)/dxsqrd-(1.190784930203603*uvarxr[3]*mu)/dxsqrd+(1.190784930203603*uvarxl[3]*mu)/dxsqrd+(0.9375*uvarxr[2]*mu)/dxsqrd+(0.9375*uvarxl[2]*mu)/dxsqrd-(1.875*uvarxc[2]*mu)/dxsqrd);
  out[19] += J*((0.7382874503707888*uvaryr[7]*mu)/dysqrd-(0.7382874503707888*uvaryl[7]*mu)/dysqrd-(1.453125*uvaryr[3]*mu)/dysqrd-(1.453125*uvaryl[3]*mu)/dysqrd-(5.34375*uvaryc[3]*mu)/dysqrd+(1.190784930203603*uvaryr[1]*mu)/dysqrd-(1.190784930203603*uvaryl[1]*mu)/dysqrd+(0.7382874503707888*uvarxr[6]*mu)/dxsqrd-(0.7382874503707888*uvarxl[6]*mu)/dxsqrd-(1.453125*uvarxr[3]*mu)/dxsqrd-(1.453125*uvarxl[3]*mu)/dxsqrd-(5.34375*uvarxc[3]*mu)/dxsqrd+(1.190784930203603*uvarxr[2]*mu)/dxsqrd-(1.190784930203603*uvarxl[2]*mu)/dxsqrd);
  out[20] += J*((-(1.190784930203603*uvaryr[6]*mu)/dysqrd)+(1.190784930203603*uvaryl[6]*mu)/dysqrd+(0.9375*uvaryr[4]*mu)/dysqrd+(0.9375*uvaryl[4]*mu)/dysqrd-(1.875*uvaryc[4]*mu)/dysqrd-(0.140625*uvarxr[4]*mu)/dxsqrd-(0.140625*uvarxl[4]*mu)/dxsqrd-(6.28125*uvarxc[4]*mu)/dxsqrd-(0.3025768239224545*uvarxr[1]*mu)/dxsqrd+(0.3025768239224545*uvarxl[1]*mu)/dxsqrd+(0.4192627457812106*uvarxr[0]*mu)/dxsqrd+(0.4192627457812106*uvarxl[0]*mu)/dxsqrd-(0.8385254915624212*uvarxc[0]*mu)/dxsqrd);
  out[21] += J*((-(0.140625*uvaryr[5]*mu)/dysqrd)-(0.140625*uvaryl[5]*mu)/dysqrd-(6.28125*uvaryc[5]*mu)/dysqrd-(0.3025768239224545*uvaryr[2]*mu)/dysqrd+(0.3025768239224545*uvaryl[2]*mu)/dysqrd+(0.4192627457812106*uvaryr[0]*mu)/dysqrd+(0.4192627457812106*uvaryl[0]*mu)/dysqrd-(0.8385254915624212*uvaryc[0]*mu)/dysqrd-(1.190784930203603*uvarxr[7]*mu)/dxsqrd+(1.190784930203603*uvarxl[7]*mu)/dxsqrd+(0.9375*uvarxr[5]*mu)/dxsqrd+(0.9375*uvarxl[5]*mu)/dxsqrd-(1.875*uvarxc[5]*mu)/dxsqrd);
  out[22] += J*((-(1.453125*uvaryr[6]*mu)/dysqrd)-(1.453125*uvaryl[6]*mu)/dysqrd-(5.34375*uvaryc[6]*mu)/dysqrd+(1.190784930203603*uvaryr[4]*mu)/dysqrd-(1.190784930203603*uvaryl[4]*mu)/dysqrd-(0.140625*uvarxr[6]*mu)/dxsqrd-(0.140625*uvarxl[6]*mu)/dxsqrd-(6.28125*uvarxc[6]*mu)/dxsqrd-(0.3025768239224544*uvarxr[3]*mu)/dxsqrd+(0.3025768239224544*uvarxl[3]*mu)/dxsqrd+(0.4192627457812105*uvarxr[2]*mu)/dxsqrd+(0.4192627457812105*uvarxl[2]*mu)/dxsqrd-(0.8385254915624211*uvarxc[2]*mu)/dxsqrd);
  out[23] += J*((-(0.140625*uvaryr[7]*mu)/dysqrd)-(0.140625*uvaryl[7]*mu)/dysqrd-(6.28125*uvaryc[7]*mu)/dysqrd-(0.3025768239224544*uvaryr[3]*mu)/dysqrd+(0.3025768239224544*uvaryl[3]*mu)/dysqrd+(0.4192627457812105*uvaryr[1]*mu)/dysqrd+(0.4192627457812105*uvaryl[1]*mu)/dysqrd-(0.8385254915624211*uvaryc[1]*mu)/dysqrd-(1.453125*uvarxr[7]*mu)/dxsqrd-(1.453125*uvarxl[7]*mu)/dxsqrd-(5.34375*uvarxc[7]*mu)/dxsqrd+(1.190784930203603*uvarxr[5]*mu)/dxsqrd-(1.190784930203603*uvarxl[5]*mu)/dxsqrd);
}
