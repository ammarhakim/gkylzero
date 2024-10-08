#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -0.0625*(8.660254037844386*coeff[0]*qr[1]-8.660254037844386*coeff[0]*ql[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.0625*(7.0*coeff[0]*qr[1]+7.0*coeff[0]*ql[1]+46.0*coeff[0]*qc[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_fluid_order2_surfx_1x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  out[0] += -0.0441941738241592*((15.0*coeff[1]+8.660254037844386*coeff[0])*qr[1]+(15.0*coeff[1]-8.660254037844386*coeff[0])*ql[1]+30.0*coeff[1]*qc[1]+(15.58845726811989*ql[0]-15.58845726811989*qr[0])*coeff[1]-9.0*coeff[0]*qr[0]-9.0*coeff[0]*ql[0]+18.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.0441941738241592*((12.12435565298214*coeff[1]+7.0*coeff[0])*qr[1]+(7.0*coeff[0]-12.12435565298214*coeff[1])*ql[1]+46.0*coeff[0]*qc[1]+((-15.0*qr[0])-15.0*ql[0]+30.0*qc[0])*coeff[1]-8.660254037844386*coeff[0]*qr[0]+8.660254037844386*coeff[0]*ql[0])*Jfac; 

  return 0.;

}

