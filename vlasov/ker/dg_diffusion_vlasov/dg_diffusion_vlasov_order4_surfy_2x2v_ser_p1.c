#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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
  out[1] += -0.0625*(25.98076211353316*coeff[1]*qr[5]-25.98076211353316*coeff[1]*ql[5]-15.0*coeff[1]*qr[1]-15.0*coeff[1]*ql[1]+30.0*coeff[1]*qc[1])*Jfac; 
  out[2] += -0.0625*(33.0*coeff[1]*qr[2]+33.0*coeff[1]*ql[2]+114.0*coeff[1]*qc[2]+(25.98076211353316*ql[0]-25.98076211353316*qr[0])*coeff[1])*Jfac; 
  out[3] += -0.0625*(25.98076211353316*coeff[1]*qr[7]-25.98076211353316*coeff[1]*ql[7]-15.0*coeff[1]*qr[3]-15.0*coeff[1]*ql[3]+30.0*coeff[1]*qc[3])*Jfac; 
  out[4] += -0.0625*(25.98076211353316*coeff[1]*qr[9]-25.98076211353316*coeff[1]*ql[9]-15.0*coeff[1]*qr[4]-15.0*coeff[1]*ql[4]+30.0*coeff[1]*qc[4])*Jfac; 
  out[5] += -0.0625*(33.0*coeff[1]*qr[5]+33.0*coeff[1]*ql[5]+114.0*coeff[1]*qc[5]-25.98076211353316*coeff[1]*qr[1]+25.98076211353316*coeff[1]*ql[1])*Jfac; 
  out[6] += -0.0625*(25.98076211353316*coeff[1]*qr[11]-25.98076211353316*coeff[1]*ql[11]-15.0*coeff[1]*qr[6]-15.0*coeff[1]*ql[6]+30.0*coeff[1]*qc[6])*Jfac; 
  out[7] += -0.0625*(33.0*coeff[1]*qr[7]+33.0*coeff[1]*ql[7]+114.0*coeff[1]*qc[7]-25.98076211353316*coeff[1]*qr[3]+25.98076211353316*coeff[1]*ql[3])*Jfac; 
  out[8] += -0.0625*(25.98076211353316*coeff[1]*qr[12]-25.98076211353316*coeff[1]*ql[12]-15.0*coeff[1]*qr[8]-15.0*coeff[1]*ql[8]+30.0*coeff[1]*qc[8])*Jfac; 
  out[9] += -0.0625*(33.0*coeff[1]*qr[9]+33.0*coeff[1]*ql[9]+114.0*coeff[1]*qc[9]-25.98076211353316*coeff[1]*qr[4]+25.98076211353316*coeff[1]*ql[4])*Jfac; 
  out[10] += -0.0625*(25.98076211353316*coeff[1]*qr[14]-25.98076211353316*coeff[1]*ql[14]-15.0*coeff[1]*qr[10]-15.0*coeff[1]*ql[10]+30.0*coeff[1]*qc[10])*Jfac; 
  out[11] += -0.0625*(33.0*coeff[1]*qr[11]+33.0*coeff[1]*ql[11]+114.0*coeff[1]*qc[11]-25.98076211353316*coeff[1]*qr[6]+25.98076211353316*coeff[1]*ql[6])*Jfac; 
  out[12] += -0.0625*(33.0*coeff[1]*qr[12]+33.0*coeff[1]*ql[12]+114.0*coeff[1]*qc[12]-25.98076211353316*coeff[1]*qr[8]+25.98076211353316*coeff[1]*ql[8])*Jfac; 
  out[13] += -0.0625*(25.98076211353316*coeff[1]*qr[15]-25.98076211353316*coeff[1]*ql[15]-15.0*coeff[1]*qr[13]-15.0*coeff[1]*ql[13]+30.0*coeff[1]*qc[13])*Jfac; 
  out[14] += -0.0625*(33.0*coeff[1]*qr[14]+33.0*coeff[1]*ql[14]+114.0*coeff[1]*qc[14]-25.98076211353316*coeff[1]*qr[10]+25.98076211353316*coeff[1]*ql[10])*Jfac; 
  out[15] += -0.0625*(33.0*coeff[1]*qr[15]+33.0*coeff[1]*ql[15]+114.0*coeff[1]*qc[15]-25.98076211353316*coeff[1]*qr[13]+25.98076211353316*coeff[1]*ql[13])*Jfac; 
  out[16] += -0.0625*(25.98076211353316*coeff[1]*qr[18]-25.98076211353316*coeff[1]*ql[18]-15.0*coeff[1]*qr[16]-15.0*coeff[1]*ql[16]+30.0*coeff[1]*qc[16])*Jfac; 
  out[17] += -0.0625*(25.98076211353316*coeff[1]*qr[20]-25.98076211353316*coeff[1]*ql[20]-15.0*coeff[1]*qr[17]-15.0*coeff[1]*ql[17]+30.0*coeff[1]*qc[17])*Jfac; 
  out[18] += -0.0625*(33.0*coeff[1]*qr[18]+33.0*coeff[1]*ql[18]+114.0*coeff[1]*qc[18]-25.98076211353316*coeff[1]*qr[16]+25.98076211353316*coeff[1]*ql[16])*Jfac; 
  out[19] += -0.0625*(25.98076211353316*coeff[1]*qr[22]-25.98076211353316*coeff[1]*ql[22]-15.0*coeff[1]*qr[19]-15.0*coeff[1]*ql[19]+30.0*coeff[1]*qc[19])*Jfac; 
  out[20] += -0.0625*(33.0*coeff[1]*qr[20]+33.0*coeff[1]*ql[20]+114.0*coeff[1]*qc[20]-25.98076211353316*coeff[1]*qr[17]+25.98076211353316*coeff[1]*ql[17])*Jfac; 
  out[21] += -0.0625*(25.98076211353316*coeff[1]*qr[23]-25.98076211353316*coeff[1]*ql[23]-15.0*coeff[1]*qr[21]-15.0*coeff[1]*ql[21]+30.0*coeff[1]*qc[21])*Jfac; 
  out[22] += -0.0625*(33.0*coeff[1]*qr[22]+33.0*coeff[1]*ql[22]+114.0*coeff[1]*qc[22]-25.98076211353316*coeff[1]*qr[19]+25.98076211353316*coeff[1]*ql[19])*Jfac; 
  out[23] += -0.0625*(33.0*coeff[1]*qr[23]+33.0*coeff[1]*ql[23]+114.0*coeff[1]*qc[23]-25.98076211353316*coeff[1]*qr[21]+25.98076211353316*coeff[1]*ql[21])*Jfac; 
  out[24] += -0.0625*(25.98076211353316*coeff[1]*qr[26]-25.98076211353316*coeff[1]*ql[26]-15.0*coeff[1]*qr[24]-15.0*coeff[1]*ql[24]+30.0*coeff[1]*qc[24])*Jfac; 
  out[25] += -0.0625*(25.98076211353316*coeff[1]*qr[28]-25.98076211353316*coeff[1]*ql[28]-15.0*coeff[1]*qr[25]-15.0*coeff[1]*ql[25]+30.0*coeff[1]*qc[25])*Jfac; 
  out[26] += -0.0625*(33.0*coeff[1]*qr[26]+33.0*coeff[1]*ql[26]+114.0*coeff[1]*qc[26]-25.98076211353316*coeff[1]*qr[24]+25.98076211353316*coeff[1]*ql[24])*Jfac; 
  out[27] += -0.0625*(25.98076211353316*coeff[1]*qr[30]-25.98076211353316*coeff[1]*ql[30]-15.0*coeff[1]*qr[27]-15.0*coeff[1]*ql[27]+30.0*coeff[1]*qc[27])*Jfac; 
  out[28] += -0.0625*(33.0*coeff[1]*qr[28]+33.0*coeff[1]*ql[28]+114.0*coeff[1]*qc[28]-25.98076211353316*coeff[1]*qr[25]+25.98076211353316*coeff[1]*ql[25])*Jfac; 
  out[29] += -0.0625*(25.98076211353316*coeff[1]*qr[31]-25.98076211353316*coeff[1]*ql[31]-15.0*coeff[1]*qr[29]-15.0*coeff[1]*ql[29]+30.0*coeff[1]*qc[29])*Jfac; 
  out[30] += -0.0625*(33.0*coeff[1]*qr[30]+33.0*coeff[1]*ql[30]+114.0*coeff[1]*qc[30]-25.98076211353316*coeff[1]*qr[27]+25.98076211353316*coeff[1]*ql[27])*Jfac; 
  out[31] += -0.0625*(33.0*coeff[1]*qr[31]+33.0*coeff[1]*ql[31]+114.0*coeff[1]*qc[31]-25.98076211353316*coeff[1]*qr[29]+25.98076211353316*coeff[1]*ql[29])*Jfac; 

  return 0.;

}

