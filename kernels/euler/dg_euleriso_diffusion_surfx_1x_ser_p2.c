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


  out[0] += J*((-(0.3977475644174328*D[2]*uvarxr[3])/dxsqrd)+(2.088192250488444*D[1]*uvarxr[3])/dxsqrd+(1.897366596101028*D[0]*uvarxr[3])/dxsqrd-(0.3977475644174328*D[2]*uvarxl[3])/dxsqrd-(2.088192250488444*D[1]*uvarxl[3])/dxsqrd+(1.897366596101028*D[0]*uvarxl[3])/dxsqrd-(17.766057877312*D[2]*uvarxc[3])/dxsqrd-(3.794733192202056*D[0]*uvarxc[3])/dxsqrd-(0.8558164961018216*D[2]*uvarxr[2])/dxsqrd-(4.110058165646805*D[1]*uvarxr[2])/dxsqrd-(3.368048396326869*D[0]*uvarxr[2])/dxsqrd+(0.8558164961018216*D[2]*uvarxl[2])/dxsqrd-(4.110058165646805*D[1]*uvarxl[2])/dxsqrd+(3.368048396326869*D[0]*uvarxl[2])/dxsqrd-(15.11440744786245*D[1]*uvarxc[2])/dxsqrd+(1.185854122563142*uvarxr[1]*D[2])/dxsqrd+(1.185854122563142*uvarxl[1]*D[2])/dxsqrd-(2.371708245126284*uvarxc[1]*D[2])/dxsqrd+(3.368048396326869*D[1]*uvarxr[1])/dxsqrd+(2.651650429449552*D[0]*uvarxr[1])/dxsqrd-(3.368048396326869*D[1]*uvarxl[1])/dxsqrd+(2.651650429449552*D[0]*uvarxl[1])/dxsqrd-(5.303300858899105*D[0]*uvarxc[1])/dxsqrd);
  out[1] += J*((-(3.368048396326869*D[2]*uvarxr[3])/dxsqrd)+(1.541610359332084*D[1]*uvarxr[3])/dxsqrd+(2.088192250488444*D[0]*uvarxr[3])/dxsqrd+(3.368048396326869*D[2]*uvarxl[3])/dxsqrd+(1.541610359332084*D[1]*uvarxl[3])/dxsqrd-(2.088192250488444*D[0]*uvarxl[3])/dxsqrd-(19.68517843454815*D[1]*uvarxc[3])/dxsqrd+(2.371708245126284*D[2]*uvarxr[2])/dxsqrd-(4.133513940946609*D[1]*uvarxr[2])/dxsqrd-(4.110058165646805*D[0]*uvarxr[2])/dxsqrd+(2.371708245126284*D[2]*uvarxl[2])/dxsqrd+(4.133513940946609*D[1]*uvarxl[2])/dxsqrd-(4.110058165646805*D[0]*uvarxl[2])/dxsqrd-(15.11440744786245*D[0]*uvarxc[2])/dxsqrd-(0.6846531968814573*uvarxr[1]*D[2])/dxsqrd+(0.6846531968814573*uvarxl[1]*D[2])/dxsqrd+(3.712310601229373*D[1]*uvarxr[1])/dxsqrd+(3.368048396326869*D[0]*uvarxr[1])/dxsqrd+(3.712310601229373*D[1]*uvarxl[1])/dxsqrd-(3.368048396326869*D[0]*uvarxl[1])/dxsqrd-(7.424621202458747*D[1]*uvarxc[1])/dxsqrd);
  out[2] += J*((-(11.26561416434985*D[2]*uvarxr[3])/dxsqrd)-(3.368048396326869*D[1]*uvarxr[3])/dxsqrd-(0.3977475644174328*D[0]*uvarxr[3])/dxsqrd-(11.26561416434985*D[2]*uvarxl[3])/dxsqrd+(3.368048396326869*D[1]*uvarxl[3])/dxsqrd-(0.3977475644174328*D[0]*uvarxl[3])/dxsqrd+(53.36343551534139*D[2]*uvarxc[3])/dxsqrd-(17.766057877312*D[0]*uvarxc[3])/dxsqrd+(13.01291425853563*D[2]*uvarxr[2])/dxsqrd+(2.371708245126284*D[1]*uvarxr[2])/dxsqrd-(0.8558164961018216*D[0]*uvarxr[2])/dxsqrd-(13.01291425853563*D[2]*uvarxl[2])/dxsqrd+(2.371708245126284*D[1]*uvarxl[2])/dxsqrd+(0.8558164961018216*D[0]*uvarxl[2])/dxsqrd-(7.954951288348656*uvarxr[1]*D[2])/dxsqrd-(7.954951288348656*uvarxl[1]*D[2])/dxsqrd+(15.90990257669731*uvarxc[1]*D[2])/dxsqrd-(0.6846531968814573*D[1]*uvarxr[1])/dxsqrd+(1.185854122563142*D[0]*uvarxr[1])/dxsqrd+(0.6846531968814573*D[1]*uvarxl[1])/dxsqrd+(1.185854122563142*D[0]*uvarxl[1])/dxsqrd-(2.371708245126284*D[0]*uvarxc[1])/dxsqrd);
}
