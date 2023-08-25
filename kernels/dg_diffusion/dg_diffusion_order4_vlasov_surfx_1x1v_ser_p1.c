#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_vlasov_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[0]*qr[1]-25.98076211353316*coeff[0]*ql[1]-15.0*coeff[0]*qr[0]-15.0*coeff[0]*ql[0]+30.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.0625*(33.0*coeff[0]*qr[1]+33.0*coeff[0]*ql[1]+114.0*coeff[0]*qc[1]-25.98076211353316*coeff[0]*qr[0]+25.98076211353316*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0625*(25.98076211353316*coeff[0]*qr[3]-25.98076211353316*coeff[0]*ql[3]-15.0*coeff[0]*qr[2]-15.0*coeff[0]*ql[2]+30.0*coeff[0]*qc[2])*Jfac; 
  out[3] += -0.0625*(33.0*coeff[0]*qr[3]+33.0*coeff[0]*ql[3]+114.0*coeff[0]*qc[3]-25.98076211353316*coeff[0]*qr[2]+25.98076211353316*coeff[0]*ql[2])*Jfac; 
  out[4] += -0.0625*(25.98076211353316*coeff[0]*qr[5]-25.98076211353316*coeff[0]*ql[5]-15.0*coeff[0]*qr[4]-15.0*coeff[0]*ql[4]+30.0*coeff[0]*qc[4])*Jfac; 
  out[5] += -0.0625*(33.0*coeff[0]*qr[5]+33.0*coeff[0]*ql[5]+114.0*coeff[0]*qc[5]-25.98076211353316*coeff[0]*qr[4]+25.98076211353316*coeff[0]*ql[4])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_order4_vlasov_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -0.03125*((80.61017305526643*coeff[1]+36.74234614174767*coeff[0])*qr[1]+(80.61017305526643*coeff[1]-36.74234614174767*coeff[0])*ql[1]+93.33809511662429*coeff[1]*qc[1]+(36.74234614174767*ql[0]-36.74234614174767*qr[0])*coeff[1]-21.21320343559643*coeff[0]*qr[0]-21.21320343559643*coeff[0]*ql[0]+42.42640687119286*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.03125*((110.227038425243*coeff[1]+46.66904755831214*coeff[0])*qr[1]+(46.66904755831214*coeff[0]-110.227038425243*coeff[1])*ql[1]+161.2203461105329*coeff[0]*qc[1]+((-63.63961030678928*qr[0])-63.63961030678928*ql[0]+127.2792206135786*qc[0])*coeff[1]-36.74234614174767*coeff[0]*qr[0]+36.74234614174767*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.03125*((80.61017305526643*coeff[1]+36.74234614174767*coeff[0])*qr[3]+(80.61017305526643*coeff[1]-36.74234614174767*coeff[0])*ql[3]+93.33809511662429*coeff[1]*qc[3]+((-36.74234614174767*coeff[1])-21.21320343559643*coeff[0])*qr[2]+(36.74234614174767*coeff[1]-21.21320343559643*coeff[0])*ql[2]+42.42640687119286*coeff[0]*qc[2])*Jfac; 
  out[3] += -0.03125*((110.227038425243*coeff[1]+46.66904755831214*coeff[0])*qr[3]+(46.66904755831214*coeff[0]-110.227038425243*coeff[1])*ql[3]+161.2203461105329*coeff[0]*qc[3]+((-63.63961030678928*coeff[1])-36.74234614174767*coeff[0])*qr[2]+(36.74234614174767*coeff[0]-63.63961030678928*coeff[1])*ql[2]+127.2792206135786*coeff[1]*qc[2])*Jfac; 
  out[4] += -0.00625*((403.0508652763321*coeff[1]+183.7117307087384*coeff[0])*qr[5]+(403.0508652763321*coeff[1]-183.7117307087384*coeff[0])*ql[5]+466.6904755831214*coeff[1]*qc[5]+((-183.7117307087383*coeff[1])-106.0660171779821*coeff[0])*qr[4]+(183.7117307087383*coeff[1]-106.0660171779821*coeff[0])*ql[4]+212.1320343559643*coeff[0]*qc[4])*Jfac; 
  out[5] += -0.03125*((110.227038425243*coeff[1]+46.66904755831214*coeff[0])*qr[5]+(46.66904755831214*coeff[0]-110.227038425243*coeff[1])*ql[5]+161.2203461105329*coeff[0]*qc[5]+((-63.63961030678928*coeff[1])-36.74234614174768*coeff[0])*qr[4]+(36.74234614174768*coeff[0]-63.63961030678928*coeff[1])*ql[4]+127.2792206135786*coeff[1]*qc[4])*Jfac; 

  return 0.;

}

