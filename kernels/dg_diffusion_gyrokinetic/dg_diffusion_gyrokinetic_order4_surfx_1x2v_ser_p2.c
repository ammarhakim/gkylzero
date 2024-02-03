#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  double fl[20];
  fl[0] = 0.7071067811865475*(jacobgeo_inv[2]*ql[7]+jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*ql[7]+ql[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1])); 
  fl[2] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*ql[11]+15.0*(jacobgeo_inv[1]*ql[4]+jacobgeo_inv[0]*ql[2])); 
  fl[3] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*ql[13]+15.0*(jacobgeo_inv[1]*ql[5]+jacobgeo_inv[0]*ql[3])); 
  fl[4] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*ql[11]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*ql[4]+15.0*jacobgeo_inv[1]*ql[2]); 
  fl[5] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*ql[13]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*ql[5]+15.0*jacobgeo_inv[1]*ql[3]); 
  fl[6] = 0.7071067811865475*(jacobgeo_inv[2]*ql[17]+jacobgeo_inv[1]*ql[10]+jacobgeo_inv[0]*ql[6]); 
  fl[7] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*ql[7]+35.0*ql[0]*jacobgeo_inv[2]+31.30495168499706*jacobgeo_inv[1]*ql[1]); 
  fl[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*ql[12]+15.0*jacobgeo_inv[0]*ql[8]); 
  fl[9] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*ql[15]+15.0*jacobgeo_inv[0]*ql[9]); 
  fl[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*ql[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*ql[10]+5.0*jacobgeo_inv[1]*ql[6]); 
  fl[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*ql[11]+93.91485505499116*jacobgeo_inv[1]*ql[4]+105.0*jacobgeo_inv[2]*ql[2]); 
  fl[12] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*ql[12]+15.0*jacobgeo_inv[1]*ql[8]); 
  fl[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*ql[13]+93.91485505499116*jacobgeo_inv[1]*ql[5]+105.0*jacobgeo_inv[2]*ql[3]); 
  fl[14] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*ql[18]+15.0*jacobgeo_inv[0]*ql[14]); 
  fl[15] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*ql[15]+15.0*jacobgeo_inv[1]*ql[9]); 
  fl[16] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*ql[19]+15.0*jacobgeo_inv[0]*ql[16]); 
  fl[17] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*ql[17]+31.30495168499706*jacobgeo_inv[1]*ql[10]+35.0*jacobgeo_inv[2]*ql[6]); 
  fl[18] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*ql[18]+15.0*jacobgeo_inv[1]*ql[14]); 
  fl[19] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*ql[19]+15.0*jacobgeo_inv[1]*ql[16]); 

  double fc[20];
  fc[0] = 0.7071067811865475*(jacobgeo_inv[2]*qc[7]+jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*qc[7]+qc[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1])); 
  fc[2] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qc[11]+15.0*(jacobgeo_inv[1]*qc[4]+jacobgeo_inv[0]*qc[2])); 
  fc[3] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qc[13]+15.0*(jacobgeo_inv[1]*qc[5]+jacobgeo_inv[0]*qc[3])); 
  fc[4] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qc[11]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qc[4]+15.0*jacobgeo_inv[1]*qc[2]); 
  fc[5] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qc[13]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qc[5]+15.0*jacobgeo_inv[1]*qc[3]); 
  fc[6] = 0.7071067811865475*(jacobgeo_inv[2]*qc[17]+jacobgeo_inv[1]*qc[10]+jacobgeo_inv[0]*qc[6]); 
  fc[7] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qc[7]+35.0*qc[0]*jacobgeo_inv[2]+31.30495168499706*jacobgeo_inv[1]*qc[1]); 
  fc[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qc[12]+15.0*jacobgeo_inv[0]*qc[8]); 
  fc[9] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qc[15]+15.0*jacobgeo_inv[0]*qc[9]); 
  fc[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*qc[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*qc[10]+5.0*jacobgeo_inv[1]*qc[6]); 
  fc[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qc[11]+93.91485505499116*jacobgeo_inv[1]*qc[4]+105.0*jacobgeo_inv[2]*qc[2]); 
  fc[12] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qc[12]+15.0*jacobgeo_inv[1]*qc[8]); 
  fc[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qc[13]+93.91485505499116*jacobgeo_inv[1]*qc[5]+105.0*jacobgeo_inv[2]*qc[3]); 
  fc[14] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qc[18]+15.0*jacobgeo_inv[0]*qc[14]); 
  fc[15] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qc[15]+15.0*jacobgeo_inv[1]*qc[9]); 
  fc[16] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qc[19]+15.0*jacobgeo_inv[0]*qc[16]); 
  fc[17] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qc[17]+31.30495168499706*jacobgeo_inv[1]*qc[10]+35.0*jacobgeo_inv[2]*qc[6]); 
  fc[18] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qc[18]+15.0*jacobgeo_inv[1]*qc[14]); 
  fc[19] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qc[19]+15.0*jacobgeo_inv[1]*qc[16]); 

  double fr[20];
  fr[0] = 0.7071067811865475*(jacobgeo_inv[2]*qr[7]+jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*qr[7]+qr[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1])); 
  fr[2] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qr[11]+15.0*(jacobgeo_inv[1]*qr[4]+jacobgeo_inv[0]*qr[2])); 
  fr[3] = 0.04714045207910316*(15.0*jacobgeo_inv[2]*qr[13]+15.0*(jacobgeo_inv[1]*qr[5]+jacobgeo_inv[0]*qr[3])); 
  fr[4] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qr[11]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qr[4]+15.0*jacobgeo_inv[1]*qr[2]); 
  fr[5] = 0.04714045207910316*(13.41640786499874*jacobgeo_inv[1]*qr[13]+(13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qr[5]+15.0*jacobgeo_inv[1]*qr[3]); 
  fr[6] = 0.7071067811865475*(jacobgeo_inv[2]*qr[17]+jacobgeo_inv[1]*qr[10]+jacobgeo_inv[0]*qr[6]); 
  fr[7] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qr[7]+35.0*qr[0]*jacobgeo_inv[2]+31.30495168499706*jacobgeo_inv[1]*qr[1]); 
  fr[8] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qr[12]+15.0*jacobgeo_inv[0]*qr[8]); 
  fr[9] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qr[15]+15.0*jacobgeo_inv[0]*qr[9]); 
  fr[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*qr[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*qr[10]+5.0*jacobgeo_inv[1]*qr[6]); 
  fr[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qr[11]+93.91485505499116*jacobgeo_inv[1]*qr[4]+105.0*jacobgeo_inv[2]*qr[2]); 
  fr[12] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qr[12]+15.0*jacobgeo_inv[1]*qr[8]); 
  fr[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qr[13]+93.91485505499116*jacobgeo_inv[1]*qr[5]+105.0*jacobgeo_inv[2]*qr[3]); 
  fr[14] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qr[18]+15.0*jacobgeo_inv[0]*qr[14]); 
  fr[15] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qr[15]+15.0*jacobgeo_inv[1]*qr[9]); 
  fr[16] = 0.04714045207910316*(15.0*jacobgeo_inv[1]*qr[19]+15.0*jacobgeo_inv[0]*qr[16]); 
  fr[17] = 0.02020305089104421*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qr[17]+31.30495168499706*jacobgeo_inv[1]*qr[10]+35.0*jacobgeo_inv[2]*qr[6]); 
  fr[18] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qr[18]+15.0*jacobgeo_inv[1]*qr[14]); 
  fr[19] = 0.04714045207910316*((13.41640786499874*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qr[19]+15.0*jacobgeo_inv[1]*qr[16]); 

  out[0] += 0.0625*(107.3312629199899*coeff[0]*fr[7]+107.3312629199899*coeff[0]*fl[7]-214.6625258399798*coeff[0]*fc[7]-129.9038105676658*coeff[0]*fr[1]+129.9038105676658*coeff[0]*fl[1]+75.0*coeff[0]*fr[0]+75.0*coeff[0]*fl[0]-150.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*fr[7]-290.4737509655563*coeff[0]*fl[7]-405.0*coeff[0]*fr[1]-405.0*coeff[0]*fl[1]-990.0*coeff[0]*fc[1]+259.8076211353315*coeff[0]*fr[0]-259.8076211353315*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += 0.0625*(107.3312629199899*coeff[0]*fr[11]+107.3312629199899*coeff[0]*fl[11]-214.6625258399798*coeff[0]*fc[11]-129.9038105676658*coeff[0]*fr[4]+129.9038105676658*coeff[0]*fl[4]+75.0*coeff[0]*fr[2]+75.0*coeff[0]*fl[2]-150.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += 0.0625*(107.3312629199899*coeff[0]*fr[13]+107.3312629199899*coeff[0]*fl[13]-214.6625258399798*coeff[0]*fc[13]-129.9038105676658*coeff[0]*fr[5]+129.9038105676658*coeff[0]*fl[5]+75.0*coeff[0]*fr[3]+75.0*coeff[0]*fl[3]-150.0*coeff[0]*fc[3])*rdx2Sq; 
  out[4] += 0.03125*(290.4737509655563*coeff[0]*fr[11]-290.4737509655563*coeff[0]*fl[11]-405.0*coeff[0]*fr[4]-405.0*coeff[0]*fl[4]-990.0*coeff[0]*fc[4]+259.8076211353315*coeff[0]*fr[2]-259.8076211353315*coeff[0]*fl[2])*rdx2Sq; 
  out[5] += 0.03125*(290.4737509655563*coeff[0]*fr[13]-290.4737509655563*coeff[0]*fl[13]-405.0*coeff[0]*fr[5]-405.0*coeff[0]*fl[5]-990.0*coeff[0]*fc[5]+259.8076211353315*coeff[0]*fr[3]-259.8076211353315*coeff[0]*fl[3])*rdx2Sq; 
  out[6] += 0.0625*(107.3312629199899*coeff[0]*fr[17]+107.3312629199899*coeff[0]*fl[17]-214.6625258399798*coeff[0]*fc[17]-129.9038105676658*coeff[0]*fr[10]+129.9038105676658*coeff[0]*fl[10]+75.0*coeff[0]*fr[6]+75.0*coeff[0]*fl[6]-150.0*coeff[0]*fc[6])*rdx2Sq; 
  out[7] += 0.03125*(21.0*coeff[0]*fr[7]+21.0*coeff[0]*fl[7]-1302.0*coeff[0]*fc[7]-151.0463505020892*coeff[0]*fr[1]+151.0463505020892*coeff[0]*fl[1]+134.1640786499874*coeff[0]*fr[0]+134.1640786499874*coeff[0]*fl[0]-268.3281572999748*coeff[0]*fc[0])*rdx2Sq; 
  out[8] += -0.0625*(129.9038105676658*coeff[0]*fr[12]-129.9038105676658*coeff[0]*fl[12]-75.0*coeff[0]*fr[8]-75.0*coeff[0]*fl[8]+150.0*coeff[0]*fc[8])*rdx2Sq; 
  out[9] += -0.0625*(129.9038105676658*coeff[0]*fr[15]-129.9038105676658*coeff[0]*fl[15]-75.0*coeff[0]*fr[9]-75.0*coeff[0]*fl[9]+150.0*coeff[0]*fc[9])*rdx2Sq; 
  out[10] += 0.03125*(290.4737509655563*coeff[0]*fr[17]-290.4737509655563*coeff[0]*fl[17]-405.0*coeff[0]*fr[10]-405.0*coeff[0]*fl[10]-990.0*coeff[0]*fc[10]+259.8076211353315*coeff[0]*fr[6]-259.8076211353315*coeff[0]*fl[6])*rdx2Sq; 
  out[11] += 0.03125*(21.0*coeff[0]*fr[11]+21.0*coeff[0]*fl[11]-1302.0*coeff[0]*fc[11]-151.0463505020893*coeff[0]*fr[4]+151.0463505020893*coeff[0]*fl[4]+134.1640786499874*coeff[0]*fr[2]+134.1640786499874*coeff[0]*fl[2]-268.3281572999747*coeff[0]*fc[2])*rdx2Sq; 
  out[12] += -0.03125*(405.0*coeff[0]*fr[12]+405.0*coeff[0]*fl[12]+990.0*coeff[0]*fc[12]-259.8076211353317*coeff[0]*fr[8]+259.8076211353317*coeff[0]*fl[8])*rdx2Sq; 
  out[13] += 0.03125*(21.0*coeff[0]*fr[13]+21.0*coeff[0]*fl[13]-1302.0*coeff[0]*fc[13]-151.0463505020893*coeff[0]*fr[5]+151.0463505020893*coeff[0]*fl[5]+134.1640786499874*coeff[0]*fr[3]+134.1640786499874*coeff[0]*fl[3]-268.3281572999747*coeff[0]*fc[3])*rdx2Sq; 
  out[14] += -0.0625*(129.9038105676658*coeff[0]*fr[18]-129.9038105676658*coeff[0]*fl[18]-75.0*coeff[0]*fr[14]-75.0*coeff[0]*fl[14]+150.0*coeff[0]*fc[14])*rdx2Sq; 
  out[15] += -0.03125*(405.0*coeff[0]*fr[15]+405.0*coeff[0]*fl[15]+990.0*coeff[0]*fc[15]-259.8076211353317*coeff[0]*fr[9]+259.8076211353317*coeff[0]*fl[9])*rdx2Sq; 
  out[16] += -0.0625*(129.9038105676658*coeff[0]*fr[19]-129.9038105676658*coeff[0]*fl[19]-75.0*coeff[0]*fr[16]-75.0*coeff[0]*fl[16]+150.0*coeff[0]*fc[16])*rdx2Sq; 
  out[17] += 0.03125*(21.0*coeff[0]*fr[17]+21.0*coeff[0]*fl[17]-1302.0*coeff[0]*fc[17]-151.0463505020892*coeff[0]*fr[10]+151.0463505020892*coeff[0]*fl[10]+134.1640786499874*coeff[0]*fr[6]+134.1640786499874*coeff[0]*fl[6]-268.3281572999748*coeff[0]*fc[6])*rdx2Sq; 
  out[18] += -0.03125*(405.0*coeff[0]*fr[18]+405.0*coeff[0]*fl[18]+990.0*coeff[0]*fc[18]-259.8076211353317*coeff[0]*fr[14]+259.8076211353317*coeff[0]*fl[14])*rdx2Sq; 
  out[19] += -0.03125*(405.0*coeff[0]*fr[19]+405.0*coeff[0]*fl[19]+990.0*coeff[0]*fc[19]-259.8076211353317*coeff[0]*fr[16]+259.8076211353317*coeff[0]*fl[16])*rdx2Sq; 

  return 0.;

}

