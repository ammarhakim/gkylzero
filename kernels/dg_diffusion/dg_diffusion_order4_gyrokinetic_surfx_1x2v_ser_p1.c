#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_gyrokinetic_surfx_1x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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
  out[2] += -0.0625*(25.98076211353316*coeff[0]*qr[4]-25.98076211353316*coeff[0]*ql[4]-15.0*coeff[0]*qr[2]-15.0*coeff[0]*ql[2]+30.0*coeff[0]*qc[2])*Jfac; 
  out[3] += -0.0625*(25.98076211353316*coeff[0]*qr[5]-25.98076211353316*coeff[0]*ql[5]-15.0*coeff[0]*qr[3]-15.0*coeff[0]*ql[3]+30.0*coeff[0]*qc[3])*Jfac; 
  out[4] += -0.0625*(33.0*coeff[0]*qr[4]+33.0*coeff[0]*ql[4]+114.0*coeff[0]*qc[4]-25.98076211353316*coeff[0]*qr[2]+25.98076211353316*coeff[0]*ql[2])*Jfac; 
  out[5] += -0.0625*(33.0*coeff[0]*qr[5]+33.0*coeff[0]*ql[5]+114.0*coeff[0]*qc[5]-25.98076211353316*coeff[0]*qr[3]+25.98076211353316*coeff[0]*ql[3])*Jfac; 
  out[6] += -0.0625*(25.98076211353316*coeff[0]*qr[7]-25.98076211353316*coeff[0]*ql[7]-15.0*coeff[0]*qr[6]-15.0*coeff[0]*ql[6]+30.0*coeff[0]*qc[6])*Jfac; 
  out[7] += -0.0625*(33.0*coeff[0]*qr[7]+33.0*coeff[0]*ql[7]+114.0*coeff[0]*qc[7]-25.98076211353316*coeff[0]*qr[6]+25.98076211353316*coeff[0]*ql[6])*Jfac; 
  out[8] += -0.0625*(25.98076211353316*coeff[0]*qr[9]-25.98076211353316*coeff[0]*ql[9]-15.0*coeff[0]*qr[8]-15.0*coeff[0]*ql[8]+30.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0625*(33.0*coeff[0]*qr[9]+33.0*coeff[0]*ql[9]+114.0*coeff[0]*qc[9]-25.98076211353316*coeff[0]*qr[8]+25.98076211353316*coeff[0]*ql[8])*Jfac; 
  out[10] += -0.0625*(25.98076211353316*coeff[0]*qr[11]-25.98076211353316*coeff[0]*ql[11]-15.0*coeff[0]*qr[10]-15.0*coeff[0]*ql[10]+30.0*coeff[0]*qc[10])*Jfac; 
  out[11] += -0.0625*(33.0*coeff[0]*qr[11]+33.0*coeff[0]*ql[11]+114.0*coeff[0]*qc[11]-25.98076211353316*coeff[0]*qr[10]+25.98076211353316*coeff[0]*ql[10])*Jfac; 

  return 0.;

}

GKYL_CU_DH double dg_diffusion_order4_gyrokinetic_surfx_1x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  out[0] += -0.0441941738241592*((57.0*coeff[1]+25.98076211353316*coeff[0])*qr[1]+(57.0*coeff[1]-25.98076211353316*coeff[0])*ql[1]+66.0*coeff[1]*qc[1]+(25.98076211353316*ql[0]-25.98076211353316*qr[0])*coeff[1]-15.0*coeff[0]*qr[0]-15.0*coeff[0]*ql[0]+30.0*coeff[0]*qc[0])*Jfac; 
  out[1] += -0.0441941738241592*((77.94228634059945*coeff[1]+33.0*coeff[0])*qr[1]+(33.0*coeff[0]-77.94228634059945*coeff[1])*ql[1]+114.0*coeff[0]*qc[1]+((-45.0*qr[0])-45.0*ql[0]+90.0*qc[0])*coeff[1]-25.98076211353316*coeff[0]*qr[0]+25.98076211353316*coeff[0]*ql[0])*Jfac; 
  out[2] += -0.0441941738241592*((57.0*coeff[1]+25.98076211353316*coeff[0])*qr[4]+(57.0*coeff[1]-25.98076211353316*coeff[0])*ql[4]+66.0*coeff[1]*qc[4]+((-25.98076211353316*coeff[1])-15.0*coeff[0])*qr[2]+(25.98076211353316*coeff[1]-15.0*coeff[0])*ql[2]+30.0*coeff[0]*qc[2])*Jfac; 
  out[3] += -0.0441941738241592*((57.0*coeff[1]+25.98076211353316*coeff[0])*qr[5]+(57.0*coeff[1]-25.98076211353316*coeff[0])*ql[5]+66.0*coeff[1]*qc[5]+((-25.98076211353316*coeff[1])-15.0*coeff[0])*qr[3]+(25.98076211353316*coeff[1]-15.0*coeff[0])*ql[3]+30.0*coeff[0]*qc[3])*Jfac; 
  out[4] += -0.0441941738241592*((77.94228634059945*coeff[1]+33.0*coeff[0])*qr[4]+(33.0*coeff[0]-77.94228634059945*coeff[1])*ql[4]+114.0*coeff[0]*qc[4]+((-45.0*coeff[1])-25.98076211353316*coeff[0])*qr[2]+(25.98076211353316*coeff[0]-45.0*coeff[1])*ql[2]+90.0*coeff[1]*qc[2])*Jfac; 
  out[5] += -0.0441941738241592*((77.94228634059945*coeff[1]+33.0*coeff[0])*qr[5]+(33.0*coeff[0]-77.94228634059945*coeff[1])*ql[5]+114.0*coeff[0]*qc[5]+((-45.0*coeff[1])-25.98076211353316*coeff[0])*qr[3]+(25.98076211353316*coeff[0]-45.0*coeff[1])*ql[3]+90.0*coeff[1]*qc[3])*Jfac; 
  out[6] += -0.0441941738241592*((57.0*coeff[1]+25.98076211353316*coeff[0])*qr[7]+(57.0*coeff[1]-25.98076211353316*coeff[0])*ql[7]+66.0*coeff[1]*qc[7]+((-25.98076211353316*coeff[1])-15.0*coeff[0])*qr[6]+(25.98076211353316*coeff[1]-15.0*coeff[0])*ql[6]+30.0*coeff[0]*qc[6])*Jfac; 
  out[7] += -0.0441941738241592*((77.94228634059945*coeff[1]+33.0*coeff[0])*qr[7]+(33.0*coeff[0]-77.94228634059945*coeff[1])*ql[7]+114.0*coeff[0]*qc[7]+((-45.0*coeff[1])-25.98076211353316*coeff[0])*qr[6]+(25.98076211353316*coeff[0]-45.0*coeff[1])*ql[6]+90.0*coeff[1]*qc[6])*Jfac; 
  out[8] += -0.008838834764831839*((285.0*coeff[1]+129.9038105676658*coeff[0])*qr[9]+(285.0*coeff[1]-129.9038105676658*coeff[0])*ql[9]+330.0000000000001*coeff[1]*qc[9]+((-129.9038105676658*coeff[1])-75.0*coeff[0])*qr[8]+(129.9038105676658*coeff[1]-75.0*coeff[0])*ql[8]+150.0*coeff[0]*qc[8])*Jfac; 
  out[9] += -0.0441941738241592*((77.94228634059945*coeff[1]+33.0*coeff[0])*qr[9]+(33.0*coeff[0]-77.94228634059945*coeff[1])*ql[9]+114.0*coeff[0]*qc[9]+((-45.0*coeff[1])-25.98076211353316*coeff[0])*qr[8]+(25.98076211353316*coeff[0]-45.0*coeff[1])*ql[8]+90.0*coeff[1]*qc[8])*Jfac; 
  out[10] += -0.008838834764831839*((285.0*coeff[1]+129.9038105676658*coeff[0])*qr[11]+(285.0*coeff[1]-129.9038105676658*coeff[0])*ql[11]+330.0000000000001*coeff[1]*qc[11]+((-129.9038105676658*coeff[1])-75.0*coeff[0])*qr[10]+(129.9038105676658*coeff[1]-75.0*coeff[0])*ql[10]+150.0*coeff[0]*qc[10])*Jfac; 
  out[11] += -0.0441941738241592*((77.94228634059945*coeff[1]+33.0*coeff[0])*qr[11]+(33.0*coeff[0]-77.94228634059945*coeff[1])*ql[11]+114.0*coeff[0]*qc[11]+((-45.0*coeff[1])-25.98076211353316*coeff[0])*qr[10]+(25.98076211353316*coeff[0]-45.0*coeff[1])*ql[10]+90.0*coeff[1]*qc[10])*Jfac; 

  return 0.;

}

