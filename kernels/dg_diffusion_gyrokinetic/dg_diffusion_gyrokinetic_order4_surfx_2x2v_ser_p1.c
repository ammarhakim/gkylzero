#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],4.);

  double fl[24];
  fl[0] = 0.5*(jacobgeo_inv[3]*ql[5]+jacobgeo_inv[2]*ql[2]+jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.5*(jacobgeo_inv[2]*ql[5]+ql[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1]); 
  fl[2] = 0.5*(jacobgeo_inv[1]*ql[5]+ql[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*ql[2]+ql[0]*jacobgeo_inv[2]); 
  fl[3] = 0.5*(jacobgeo_inv[3]*ql[11]+jacobgeo_inv[2]*ql[7]+jacobgeo_inv[1]*ql[6]+jacobgeo_inv[0]*ql[3]); 
  fl[4] = 0.5*(jacobgeo_inv[3]*ql[12]+jacobgeo_inv[2]*ql[9]+jacobgeo_inv[1]*ql[8]+jacobgeo_inv[0]*ql[4]); 
  fl[5] = 0.5*(jacobgeo_inv[0]*ql[5]+ql[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*ql[2]+ql[1]*jacobgeo_inv[2]); 
  fl[6] = 0.5*(jacobgeo_inv[2]*ql[11]+jacobgeo_inv[3]*ql[7]+jacobgeo_inv[0]*ql[6]+jacobgeo_inv[1]*ql[3]); 
  fl[7] = 0.5*(jacobgeo_inv[1]*ql[11]+jacobgeo_inv[0]*ql[7]+jacobgeo_inv[3]*ql[6]+jacobgeo_inv[2]*ql[3]); 
  fl[8] = 0.5*(jacobgeo_inv[2]*ql[12]+jacobgeo_inv[3]*ql[9]+jacobgeo_inv[0]*ql[8]+jacobgeo_inv[1]*ql[4]); 
  fl[9] = 0.5*(jacobgeo_inv[1]*ql[12]+jacobgeo_inv[0]*ql[9]+jacobgeo_inv[3]*ql[8]+jacobgeo_inv[2]*ql[4]); 
  fl[10] = 0.5*(jacobgeo_inv[3]*ql[15]+jacobgeo_inv[2]*ql[14]+jacobgeo_inv[1]*ql[13]+jacobgeo_inv[0]*ql[10]); 
  fl[11] = 0.5*(jacobgeo_inv[0]*ql[11]+jacobgeo_inv[1]*ql[7]+jacobgeo_inv[2]*ql[6]+jacobgeo_inv[3]*ql[3]); 
  fl[12] = 0.5*(jacobgeo_inv[0]*ql[12]+jacobgeo_inv[1]*ql[9]+jacobgeo_inv[2]*ql[8]+jacobgeo_inv[3]*ql[4]); 
  fl[13] = 0.5*(jacobgeo_inv[2]*ql[15]+jacobgeo_inv[3]*ql[14]+jacobgeo_inv[0]*ql[13]+jacobgeo_inv[1]*ql[10]); 
  fl[14] = 0.5*(jacobgeo_inv[1]*ql[15]+jacobgeo_inv[0]*ql[14]+jacobgeo_inv[3]*ql[13]+jacobgeo_inv[2]*ql[10]); 
  fl[15] = 0.5*(jacobgeo_inv[0]*ql[15]+jacobgeo_inv[1]*ql[14]+jacobgeo_inv[2]*ql[13]+jacobgeo_inv[3]*ql[10]); 
  fl[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*ql[20]+15.0*(jacobgeo_inv[2]*ql[18]+jacobgeo_inv[1]*ql[17])+15.0*jacobgeo_inv[0]*ql[16]); 
  fl[17] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*ql[20]+15.0*(jacobgeo_inv[3]*ql[18]+jacobgeo_inv[0]*ql[17])+15.0*jacobgeo_inv[1]*ql[16]); 
  fl[18] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*ql[20]+15.0*(jacobgeo_inv[0]*ql[18]+jacobgeo_inv[3]*ql[17])+15.0*jacobgeo_inv[2]*ql[16]); 
  fl[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*ql[23]+15.0*(jacobgeo_inv[2]*ql[22]+jacobgeo_inv[1]*ql[21])+15.0*jacobgeo_inv[0]*ql[19]); 
  fl[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*ql[20]+15.0*(jacobgeo_inv[1]*ql[18]+jacobgeo_inv[2]*ql[17])+15.0*jacobgeo_inv[3]*ql[16]); 
  fl[21] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*ql[23]+15.0*(jacobgeo_inv[3]*ql[22]+jacobgeo_inv[0]*ql[21])+15.0*jacobgeo_inv[1]*ql[19]); 
  fl[22] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*ql[23]+15.0*(jacobgeo_inv[0]*ql[22]+jacobgeo_inv[3]*ql[21])+15.0*jacobgeo_inv[2]*ql[19]); 
  fl[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*ql[23]+15.0*(jacobgeo_inv[1]*ql[22]+jacobgeo_inv[2]*ql[21])+15.0*jacobgeo_inv[3]*ql[19]); 

  double fc[24];
  fc[0] = 0.5*(jacobgeo_inv[3]*qc[5]+jacobgeo_inv[2]*qc[2]+jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.5*(jacobgeo_inv[2]*qc[5]+qc[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1]); 
  fc[2] = 0.5*(jacobgeo_inv[1]*qc[5]+qc[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qc[2]+qc[0]*jacobgeo_inv[2]); 
  fc[3] = 0.5*(jacobgeo_inv[3]*qc[11]+jacobgeo_inv[2]*qc[7]+jacobgeo_inv[1]*qc[6]+jacobgeo_inv[0]*qc[3]); 
  fc[4] = 0.5*(jacobgeo_inv[3]*qc[12]+jacobgeo_inv[2]*qc[9]+jacobgeo_inv[1]*qc[8]+jacobgeo_inv[0]*qc[4]); 
  fc[5] = 0.5*(jacobgeo_inv[0]*qc[5]+qc[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qc[2]+qc[1]*jacobgeo_inv[2]); 
  fc[6] = 0.5*(jacobgeo_inv[2]*qc[11]+jacobgeo_inv[3]*qc[7]+jacobgeo_inv[0]*qc[6]+jacobgeo_inv[1]*qc[3]); 
  fc[7] = 0.5*(jacobgeo_inv[1]*qc[11]+jacobgeo_inv[0]*qc[7]+jacobgeo_inv[3]*qc[6]+jacobgeo_inv[2]*qc[3]); 
  fc[8] = 0.5*(jacobgeo_inv[2]*qc[12]+jacobgeo_inv[3]*qc[9]+jacobgeo_inv[0]*qc[8]+jacobgeo_inv[1]*qc[4]); 
  fc[9] = 0.5*(jacobgeo_inv[1]*qc[12]+jacobgeo_inv[0]*qc[9]+jacobgeo_inv[3]*qc[8]+jacobgeo_inv[2]*qc[4]); 
  fc[10] = 0.5*(jacobgeo_inv[3]*qc[15]+jacobgeo_inv[2]*qc[14]+jacobgeo_inv[1]*qc[13]+jacobgeo_inv[0]*qc[10]); 
  fc[11] = 0.5*(jacobgeo_inv[0]*qc[11]+jacobgeo_inv[1]*qc[7]+jacobgeo_inv[2]*qc[6]+jacobgeo_inv[3]*qc[3]); 
  fc[12] = 0.5*(jacobgeo_inv[0]*qc[12]+jacobgeo_inv[1]*qc[9]+jacobgeo_inv[2]*qc[8]+jacobgeo_inv[3]*qc[4]); 
  fc[13] = 0.5*(jacobgeo_inv[2]*qc[15]+jacobgeo_inv[3]*qc[14]+jacobgeo_inv[0]*qc[13]+jacobgeo_inv[1]*qc[10]); 
  fc[14] = 0.5*(jacobgeo_inv[1]*qc[15]+jacobgeo_inv[0]*qc[14]+jacobgeo_inv[3]*qc[13]+jacobgeo_inv[2]*qc[10]); 
  fc[15] = 0.5*(jacobgeo_inv[0]*qc[15]+jacobgeo_inv[1]*qc[14]+jacobgeo_inv[2]*qc[13]+jacobgeo_inv[3]*qc[10]); 
  fc[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qc[20]+15.0*(jacobgeo_inv[2]*qc[18]+jacobgeo_inv[1]*qc[17])+15.0*jacobgeo_inv[0]*qc[16]); 
  fc[17] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qc[20]+15.0*(jacobgeo_inv[3]*qc[18]+jacobgeo_inv[0]*qc[17])+15.0*jacobgeo_inv[1]*qc[16]); 
  fc[18] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qc[20]+15.0*(jacobgeo_inv[0]*qc[18]+jacobgeo_inv[3]*qc[17])+15.0*jacobgeo_inv[2]*qc[16]); 
  fc[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qc[23]+15.0*(jacobgeo_inv[2]*qc[22]+jacobgeo_inv[1]*qc[21])+15.0*jacobgeo_inv[0]*qc[19]); 
  fc[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qc[20]+15.0*(jacobgeo_inv[1]*qc[18]+jacobgeo_inv[2]*qc[17])+15.0*jacobgeo_inv[3]*qc[16]); 
  fc[21] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qc[23]+15.0*(jacobgeo_inv[3]*qc[22]+jacobgeo_inv[0]*qc[21])+15.0*jacobgeo_inv[1]*qc[19]); 
  fc[22] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qc[23]+15.0*(jacobgeo_inv[0]*qc[22]+jacobgeo_inv[3]*qc[21])+15.0*jacobgeo_inv[2]*qc[19]); 
  fc[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qc[23]+15.0*(jacobgeo_inv[1]*qc[22]+jacobgeo_inv[2]*qc[21])+15.0*jacobgeo_inv[3]*qc[19]); 

  double fr[24];
  fr[0] = 0.5*(jacobgeo_inv[3]*qr[5]+jacobgeo_inv[2]*qr[2]+jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.5*(jacobgeo_inv[2]*qr[5]+qr[2]*jacobgeo_inv[3]+jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1]); 
  fr[2] = 0.5*(jacobgeo_inv[1]*qr[5]+qr[1]*jacobgeo_inv[3]+jacobgeo_inv[0]*qr[2]+qr[0]*jacobgeo_inv[2]); 
  fr[3] = 0.5*(jacobgeo_inv[3]*qr[11]+jacobgeo_inv[2]*qr[7]+jacobgeo_inv[1]*qr[6]+jacobgeo_inv[0]*qr[3]); 
  fr[4] = 0.5*(jacobgeo_inv[3]*qr[12]+jacobgeo_inv[2]*qr[9]+jacobgeo_inv[1]*qr[8]+jacobgeo_inv[0]*qr[4]); 
  fr[5] = 0.5*(jacobgeo_inv[0]*qr[5]+qr[0]*jacobgeo_inv[3]+jacobgeo_inv[1]*qr[2]+qr[1]*jacobgeo_inv[2]); 
  fr[6] = 0.5*(jacobgeo_inv[2]*qr[11]+jacobgeo_inv[3]*qr[7]+jacobgeo_inv[0]*qr[6]+jacobgeo_inv[1]*qr[3]); 
  fr[7] = 0.5*(jacobgeo_inv[1]*qr[11]+jacobgeo_inv[0]*qr[7]+jacobgeo_inv[3]*qr[6]+jacobgeo_inv[2]*qr[3]); 
  fr[8] = 0.5*(jacobgeo_inv[2]*qr[12]+jacobgeo_inv[3]*qr[9]+jacobgeo_inv[0]*qr[8]+jacobgeo_inv[1]*qr[4]); 
  fr[9] = 0.5*(jacobgeo_inv[1]*qr[12]+jacobgeo_inv[0]*qr[9]+jacobgeo_inv[3]*qr[8]+jacobgeo_inv[2]*qr[4]); 
  fr[10] = 0.5*(jacobgeo_inv[3]*qr[15]+jacobgeo_inv[2]*qr[14]+jacobgeo_inv[1]*qr[13]+jacobgeo_inv[0]*qr[10]); 
  fr[11] = 0.5*(jacobgeo_inv[0]*qr[11]+jacobgeo_inv[1]*qr[7]+jacobgeo_inv[2]*qr[6]+jacobgeo_inv[3]*qr[3]); 
  fr[12] = 0.5*(jacobgeo_inv[0]*qr[12]+jacobgeo_inv[1]*qr[9]+jacobgeo_inv[2]*qr[8]+jacobgeo_inv[3]*qr[4]); 
  fr[13] = 0.5*(jacobgeo_inv[2]*qr[15]+jacobgeo_inv[3]*qr[14]+jacobgeo_inv[0]*qr[13]+jacobgeo_inv[1]*qr[10]); 
  fr[14] = 0.5*(jacobgeo_inv[1]*qr[15]+jacobgeo_inv[0]*qr[14]+jacobgeo_inv[3]*qr[13]+jacobgeo_inv[2]*qr[10]); 
  fr[15] = 0.5*(jacobgeo_inv[0]*qr[15]+jacobgeo_inv[1]*qr[14]+jacobgeo_inv[2]*qr[13]+jacobgeo_inv[3]*qr[10]); 
  fr[16] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qr[20]+15.0*(jacobgeo_inv[2]*qr[18]+jacobgeo_inv[1]*qr[17])+15.0*jacobgeo_inv[0]*qr[16]); 
  fr[17] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qr[20]+15.0*(jacobgeo_inv[3]*qr[18]+jacobgeo_inv[0]*qr[17])+15.0*jacobgeo_inv[1]*qr[16]); 
  fr[18] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qr[20]+15.0*(jacobgeo_inv[0]*qr[18]+jacobgeo_inv[3]*qr[17])+15.0*jacobgeo_inv[2]*qr[16]); 
  fr[19] = 0.03333333333333333*(15.0*jacobgeo_inv[3]*qr[23]+15.0*(jacobgeo_inv[2]*qr[22]+jacobgeo_inv[1]*qr[21])+15.0*jacobgeo_inv[0]*qr[19]); 
  fr[20] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qr[20]+15.0*(jacobgeo_inv[1]*qr[18]+jacobgeo_inv[2]*qr[17])+15.0*jacobgeo_inv[3]*qr[16]); 
  fr[21] = 0.03333333333333333*(15.0*jacobgeo_inv[2]*qr[23]+15.0*(jacobgeo_inv[3]*qr[22]+jacobgeo_inv[0]*qr[21])+15.0*jacobgeo_inv[1]*qr[19]); 
  fr[22] = 0.03333333333333333*(15.0*jacobgeo_inv[1]*qr[23]+15.0*(jacobgeo_inv[0]*qr[22]+jacobgeo_inv[3]*qr[21])+15.0*jacobgeo_inv[2]*qr[19]); 
  fr[23] = 0.03333333333333333*(15.0*jacobgeo_inv[0]*qr[23]+15.0*(jacobgeo_inv[1]*qr[22]+jacobgeo_inv[2]*qr[21])+15.0*jacobgeo_inv[3]*qr[19]); 

  out[0] += -0.0625*(25.98076211353316*coeff[0]*fr[1]-25.98076211353316*coeff[0]*fl[1]-15.0*coeff[0]*fr[0]-15.0*coeff[0]*fl[0]+30.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += -0.0625*(33.0*coeff[0]*fr[1]+33.0*coeff[0]*fl[1]+114.0*coeff[0]*fc[1]-25.98076211353316*coeff[0]*fr[0]+25.98076211353316*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += -0.0625*(25.98076211353316*coeff[0]*fr[5]-25.98076211353316*coeff[0]*fl[5]-15.0*coeff[0]*fr[2]-15.0*coeff[0]*fl[2]+30.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += -0.0625*(25.98076211353316*coeff[0]*fr[6]-25.98076211353316*coeff[0]*fl[6]-15.0*coeff[0]*fr[3]-15.0*coeff[0]*fl[3]+30.0*coeff[0]*fc[3])*rdx2Sq; 
  out[4] += -0.0625*(25.98076211353316*coeff[0]*fr[8]-25.98076211353316*coeff[0]*fl[8]-15.0*coeff[0]*fr[4]-15.0*coeff[0]*fl[4]+30.0*coeff[0]*fc[4])*rdx2Sq; 
  out[5] += -0.0625*(33.0*coeff[0]*fr[5]+33.0*coeff[0]*fl[5]+114.0*coeff[0]*fc[5]-25.98076211353316*coeff[0]*fr[2]+25.98076211353316*coeff[0]*fl[2])*rdx2Sq; 
  out[6] += -0.0625*(33.0*coeff[0]*fr[6]+33.0*coeff[0]*fl[6]+114.0*coeff[0]*fc[6]-25.98076211353316*coeff[0]*fr[3]+25.98076211353316*coeff[0]*fl[3])*rdx2Sq; 
  out[7] += -0.0625*(25.98076211353316*coeff[0]*fr[11]-25.98076211353316*coeff[0]*fl[11]-15.0*coeff[0]*fr[7]-15.0*coeff[0]*fl[7]+30.0*coeff[0]*fc[7])*rdx2Sq; 
  out[8] += -0.0625*(33.0*coeff[0]*fr[8]+33.0*coeff[0]*fl[8]+114.0*coeff[0]*fc[8]-25.98076211353316*coeff[0]*fr[4]+25.98076211353316*coeff[0]*fl[4])*rdx2Sq; 
  out[9] += -0.0625*(25.98076211353316*coeff[0]*fr[12]-25.98076211353316*coeff[0]*fl[12]-15.0*coeff[0]*fr[9]-15.0*coeff[0]*fl[9]+30.0*coeff[0]*fc[9])*rdx2Sq; 
  out[10] += -0.0625*(25.98076211353316*coeff[0]*fr[13]-25.98076211353316*coeff[0]*fl[13]-15.0*coeff[0]*fr[10]-15.0*coeff[0]*fl[10]+30.0*coeff[0]*fc[10])*rdx2Sq; 
  out[11] += -0.0625*(33.0*coeff[0]*fr[11]+33.0*coeff[0]*fl[11]+114.0*coeff[0]*fc[11]-25.98076211353316*coeff[0]*fr[7]+25.98076211353316*coeff[0]*fl[7])*rdx2Sq; 
  out[12] += -0.0625*(33.0*coeff[0]*fr[12]+33.0*coeff[0]*fl[12]+114.0*coeff[0]*fc[12]-25.98076211353316*coeff[0]*fr[9]+25.98076211353316*coeff[0]*fl[9])*rdx2Sq; 
  out[13] += -0.0625*(33.0*coeff[0]*fr[13]+33.0*coeff[0]*fl[13]+114.0*coeff[0]*fc[13]-25.98076211353316*coeff[0]*fr[10]+25.98076211353316*coeff[0]*fl[10])*rdx2Sq; 
  out[14] += -0.0625*(25.98076211353316*coeff[0]*fr[15]-25.98076211353316*coeff[0]*fl[15]-15.0*coeff[0]*fr[14]-15.0*coeff[0]*fl[14]+30.0*coeff[0]*fc[14])*rdx2Sq; 
  out[15] += -0.0625*(33.0*coeff[0]*fr[15]+33.0*coeff[0]*fl[15]+114.0*coeff[0]*fc[15]-25.98076211353316*coeff[0]*fr[14]+25.98076211353316*coeff[0]*fl[14])*rdx2Sq; 
  out[16] += -0.0625*(25.98076211353316*coeff[0]*fr[17]-25.98076211353316*coeff[0]*fl[17]-15.0*coeff[0]*fr[16]-15.0*coeff[0]*fl[16]+30.0*coeff[0]*fc[16])*rdx2Sq; 
  out[17] += -0.0625*(33.0*coeff[0]*fr[17]+33.0*coeff[0]*fl[17]+114.0*coeff[0]*fc[17]-25.98076211353316*coeff[0]*fr[16]+25.98076211353316*coeff[0]*fl[16])*rdx2Sq; 
  out[18] += -0.0625*(25.98076211353316*coeff[0]*fr[20]-25.98076211353316*coeff[0]*fl[20]-15.0*coeff[0]*fr[18]-15.0*coeff[0]*fl[18]+30.0*coeff[0]*fc[18])*rdx2Sq; 
  out[19] += -0.0625*(25.98076211353316*coeff[0]*fr[21]-25.98076211353316*coeff[0]*fl[21]-15.0*coeff[0]*fr[19]-15.0*coeff[0]*fl[19]+30.0*coeff[0]*fc[19])*rdx2Sq; 
  out[20] += -0.0625*(33.0*coeff[0]*fr[20]+33.0*coeff[0]*fl[20]+114.0*coeff[0]*fc[20]-25.98076211353316*coeff[0]*fr[18]+25.98076211353316*coeff[0]*fl[18])*rdx2Sq; 
  out[21] += -0.0625*(33.0*coeff[0]*fr[21]+33.0*coeff[0]*fl[21]+114.0*coeff[0]*fc[21]-25.98076211353316*coeff[0]*fr[19]+25.98076211353316*coeff[0]*fl[19])*rdx2Sq; 
  out[22] += -0.0625*(25.98076211353316*coeff[0]*fr[23]-25.98076211353316*coeff[0]*fl[23]-15.0*coeff[0]*fr[22]-15.0*coeff[0]*fl[22]+30.0*coeff[0]*fc[22])*rdx2Sq; 
  out[23] += -0.0625*(33.0*coeff[0]*fr[23]+33.0*coeff[0]*fl[23]+114.0*coeff[0]*fc[23]-25.98076211353316*coeff[0]*fr[22]+25.98076211353316*coeff[0]*fl[22])*rdx2Sq; 

  return 0.;

}

