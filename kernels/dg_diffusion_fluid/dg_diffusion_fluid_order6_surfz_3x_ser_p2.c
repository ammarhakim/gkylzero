#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_surfz_3x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],6.);

  out[0] += 0.0625*(563.489130329947*coeff[2]*qr[9]+563.489130329947*coeff[2]*ql[9]-1126.978260659894*coeff[2]*qc[9]-545.5960043841961*coeff[2]*qr[3]+545.5960043841961*coeff[2]*ql[3]+(315.0*qr[0]+315.0*ql[0]-630.0*qc[0])*coeff[2])*Jfac; 
  out[1] += 0.0625*(563.4891303299469*coeff[2]*qr[15]+563.4891303299469*coeff[2]*ql[15]-1126.978260659894*coeff[2]*qc[15]-545.5960043841961*coeff[2]*qr[5]+545.5960043841961*coeff[2]*ql[5]+(315.0*qr[1]+315.0*ql[1]-630.0*qc[1])*coeff[2])*Jfac; 
  out[2] += 0.0625*(563.4891303299469*coeff[2]*qr[16]+563.4891303299469*coeff[2]*ql[16]-1126.978260659894*coeff[2]*qc[16]-545.5960043841961*coeff[2]*qr[6]+545.5960043841961*coeff[2]*ql[6]+315.0*coeff[2]*qr[2]+315.0*coeff[2]*ql[2]-630.0*coeff[2]*qc[2])*Jfac; 
  out[3] += 0.0078125*(6587.944671898812*coeff[2]*qr[9]-6587.944671898812*coeff[2]*ql[9]-7245.0*coeff[2]*qr[3]-7245.0*coeff[2]*ql[3]-15750.0*coeff[2]*qc[3]+(4364.768035073569*qr[0]-4364.768035073569*ql[0])*coeff[2])*Jfac; 
  out[4] += 0.0625*(563.489130329947*coeff[2]*qr[19]+563.489130329947*coeff[2]*ql[19]-1126.978260659894*coeff[2]*qc[19]-545.5960043841961*coeff[2]*qr[10]+545.5960043841961*coeff[2]*ql[10]+315.0*coeff[2]*qr[4]+315.0*coeff[2]*ql[4]-630.0*coeff[2]*qc[4])*Jfac; 
  out[5] += 0.0078125*(6587.944671898817*coeff[2]*qr[15]-6587.944671898817*coeff[2]*ql[15]-7245.0*coeff[2]*qr[5]-7245.0*coeff[2]*ql[5]-15750.0*coeff[2]*qc[5]+(4364.768035073569*qr[1]-4364.768035073569*ql[1])*coeff[2])*Jfac; 
  out[6] += 0.0078125*(6587.944671898817*coeff[2]*qr[16]-6587.944671898817*coeff[2]*ql[16]-7245.0*coeff[2]*qr[6]-7245.0*coeff[2]*ql[6]-15750.0*coeff[2]*qc[6]+4364.768035073569*coeff[2]*qr[2]-4364.768035073569*coeff[2]*ql[2])*Jfac; 
  out[7] += -0.0625*(545.5960043841964*coeff[2]*qr[13]-545.5960043841964*coeff[2]*ql[13]-315.0*coeff[2]*qr[7]-315.0*coeff[2]*ql[7]+630.0*coeff[2]*qc[7])*Jfac; 
  out[8] += -0.0625*(545.5960043841964*coeff[2]*qr[14]-545.5960043841964*coeff[2]*ql[14]-315.0*coeff[2]*qr[8]-315.0*coeff[2]*ql[8]+630.0*coeff[2]*qc[8])*Jfac; 
  out[9] += -0.0078125*(405.0*coeff[2]*qr[9]+405.0*coeff[2]*ql[9]+18090.0*coeff[2]*qc[9]+1568.558255214003*coeff[2]*qr[3]-1568.558255214003*coeff[2]*ql[3]+((-1609.968943799849*qr[0])-1609.968943799849*ql[0]+3219.937887599698*qc[0])*coeff[2])*Jfac; 
  out[10] += 0.0078125*(6587.944671898812*coeff[2]*qr[19]-6587.944671898812*coeff[2]*ql[19]-7245.0*coeff[2]*qr[10]-7245.0*coeff[2]*ql[10]-15750.0*coeff[2]*qc[10]+4364.768035073569*coeff[2]*qr[4]-4364.768035073569*coeff[2]*ql[4])*Jfac; 
  out[11] += -0.0625*(545.5960043841964*coeff[2]*qr[17]-545.5960043841964*coeff[2]*ql[17]-315.0*coeff[2]*qr[11]-315.0*coeff[2]*ql[11]+630.0*coeff[2]*qc[11])*Jfac; 
  out[12] += -0.0625*(545.5960043841964*coeff[2]*qr[18]-545.5960043841964*coeff[2]*ql[18]-315.0*coeff[2]*qr[12]-315.0*coeff[2]*ql[12]+630.0*coeff[2]*qc[12])*Jfac; 
  out[13] += -0.0078125*(7245.0*coeff[2]*qr[13]+7245.0*coeff[2]*ql[13]+15750.0*coeff[2]*qc[13]-4364.768035073571*coeff[2]*qr[7]+4364.768035073571*coeff[2]*ql[7])*Jfac; 
  out[14] += -0.0078125*(7245.0*coeff[2]*qr[14]+7245.0*coeff[2]*ql[14]+15750.0*coeff[2]*qc[14]-4364.768035073571*coeff[2]*qr[8]+4364.768035073571*coeff[2]*ql[8])*Jfac; 
  out[15] += -0.0078125*(405.0*coeff[2]*qr[15]+405.0*coeff[2]*ql[15]+18090.0*coeff[2]*qc[15]+1568.558255214004*coeff[2]*qr[5]-1568.558255214004*coeff[2]*ql[5]+((-1609.968943799848*qr[1])-1609.968943799848*ql[1]+3219.937887599697*qc[1])*coeff[2])*Jfac; 
  out[16] += -0.0078125*(405.0*coeff[2]*qr[16]+405.0*coeff[2]*ql[16]+18090.0*coeff[2]*qc[16]+1568.558255214004*coeff[2]*qr[6]-1568.558255214004*coeff[2]*ql[6]-1609.968943799848*coeff[2]*qr[2]-1609.968943799848*coeff[2]*ql[2]+3219.937887599697*coeff[2]*qc[2])*Jfac; 
  out[17] += -0.0078125*(7245.0*coeff[2]*qr[17]+7245.0*coeff[2]*ql[17]+15750.0*coeff[2]*qc[17]-4364.768035073571*coeff[2]*qr[11]+4364.768035073571*coeff[2]*ql[11])*Jfac; 
  out[18] += -0.0078125*(7245.0*coeff[2]*qr[18]+7245.0*coeff[2]*ql[18]+15750.0*coeff[2]*qc[18]-4364.768035073571*coeff[2]*qr[12]+4364.768035073571*coeff[2]*ql[12])*Jfac; 
  out[19] += -0.0078125*(405.0*coeff[2]*qr[19]+405.0*coeff[2]*ql[19]+18090.0*coeff[2]*qc[19]+1568.558255214003*coeff[2]*qr[10]-1568.558255214003*coeff[2]*ql[10]-1609.968943799849*coeff[2]*qr[4]-1609.968943799849*coeff[2]*ql[4]+3219.937887599698*coeff[2]*qc[4])*Jfac; 

  return 0.;

}

