#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += -0.0625*(25.98076211353316*coeff[1]*qr[2]-25.98076211353316*coeff[1]*ql[2]+((-15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[1])*Jfac; 
  out[1] += -0.0625*(25.98076211353316*coeff[1]*qr[3]-25.98076211353316*coeff[1]*ql[3]-15.0*coeff[1]*qr[1]-15.0*coeff[1]*ql[1]+30.0*coeff[1]*qc[1])*Jfac; 
  out[2] += -0.0625*(33.0*coeff[1]*qr[2]+33.0*coeff[1]*ql[2]+114.0*coeff[1]*qc[2]+(25.98076211353316*ql[0]-25.98076211353316*qr[0])*coeff[1])*Jfac; 
  out[3] += -0.0625*(33.0*coeff[1]*qr[3]+33.0*coeff[1]*ql[3]+114.0*coeff[1]*qc[3]-25.98076211353316*coeff[1]*qr[1]+25.98076211353316*coeff[1]*ql[1])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order4_surfy_2x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  out[0] += -0.03125*((57.0*qr[3]+57.0*ql[3]+66.0*qc[3]-25.98076211353316*qr[1]+25.98076211353316*ql[1])*coeff[7]+(57.0*qr[2]+57.0*ql[2]+66.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[6]+(25.98076211353316*qr[3]-25.98076211353316*ql[3]-15.0*qr[1]-15.0*ql[1]+30.0*qc[1])*coeff[5]+(25.98076211353316*qr[2]-25.98076211353316*ql[2]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[4])*Jfac; 
  out[1] += -0.03125*((57.0*qr[2]+57.0*ql[2]+66.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[7]+(57.0*qr[3]+57.0*ql[3]+66.0*qc[3]-25.98076211353316*qr[1]+25.98076211353316*ql[1])*coeff[6]+(25.98076211353316*qr[2]-25.98076211353316*ql[2]-15.0*qr[0]-15.0*ql[0]+30.0*qc[0])*coeff[5]+(25.98076211353316*qr[3]-25.98076211353316*ql[3]-15.0*qr[1]-15.0*ql[1]+30.0*qc[1])*coeff[4])*Jfac; 
  out[2] += -0.03125*((77.94228634059945*qr[3]-77.94228634059945*ql[3]-45.0*qr[1]-45.0*ql[1]+90.0*qc[1])*coeff[7]+(77.94228634059945*qr[2]-77.94228634059945*ql[2]-45.0*qr[0]-45.0*ql[0]+90.0*qc[0])*coeff[6]+(33.0*qr[3]+33.0*ql[3]+114.0*qc[3]-25.98076211353316*qr[1]+25.98076211353316*ql[1])*coeff[5]+(33.0*qr[2]+33.0*ql[2]+114.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[4])*Jfac; 
  out[3] += -0.03125*((77.94228634059945*qr[2]-77.94228634059945*ql[2]-45.0*qr[0]-45.0*ql[0]+90.0*qc[0])*coeff[7]+(77.94228634059945*qr[3]-77.94228634059945*ql[3]-45.0*qr[1]-45.0*ql[1]+90.0*qc[1])*coeff[6]+(33.0*qr[2]+33.0*ql[2]+114.0*qc[2]-25.98076211353316*qr[0]+25.98076211353316*ql[0])*coeff[5]+(33.0*qr[3]+33.0*ql[3]+114.0*qc[3]-25.98076211353316*qr[1]+25.98076211353316*ql[1])*coeff[4])*Jfac; 

  return 0.;

}

