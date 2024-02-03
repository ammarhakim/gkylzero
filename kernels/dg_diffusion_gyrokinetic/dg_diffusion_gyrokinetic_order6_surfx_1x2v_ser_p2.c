#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // ql: Input field in the left cell.
  // qc: Input field in the center cell.
  // qr: Input field in the right cell.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],6.);

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

  out[0] += 0.0625*(563.489130329947*coeff[0]*fr[7]+563.489130329947*coeff[0]*fl[7]-1126.978260659894*coeff[0]*fc[7]-545.5960043841961*coeff[0]*fr[1]+545.5960043841961*coeff[0]*fl[1]+315.0*coeff[0]*fr[0]+315.0*coeff[0]*fl[0]-630.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*fr[7]-6587.944671898812*coeff[0]*fl[7]-7245.0*coeff[0]*fr[1]-7245.0*coeff[0]*fl[1]-15750.0*coeff[0]*fc[1]+4364.768035073569*coeff[0]*fr[0]-4364.768035073569*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*fr[11]+563.4891303299469*coeff[0]*fl[11]-1126.978260659894*coeff[0]*fc[11]-545.5960043841961*coeff[0]*fr[4]+545.5960043841961*coeff[0]*fl[4]+315.0*coeff[0]*fr[2]+315.0*coeff[0]*fl[2]-630.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += 0.0625*(563.4891303299469*coeff[0]*fr[13]+563.4891303299469*coeff[0]*fl[13]-1126.978260659894*coeff[0]*fc[13]-545.5960043841961*coeff[0]*fr[5]+545.5960043841961*coeff[0]*fl[5]+315.0*coeff[0]*fr[3]+315.0*coeff[0]*fl[3]-630.0*coeff[0]*fc[3])*rdx2Sq; 
  out[4] += 0.0078125*(6587.944671898817*coeff[0]*fr[11]-6587.944671898817*coeff[0]*fl[11]-7245.0*coeff[0]*fr[4]-7245.0*coeff[0]*fl[4]-15750.0*coeff[0]*fc[4]+4364.768035073569*coeff[0]*fr[2]-4364.768035073569*coeff[0]*fl[2])*rdx2Sq; 
  out[5] += 0.0078125*(6587.944671898817*coeff[0]*fr[13]-6587.944671898817*coeff[0]*fl[13]-7245.0*coeff[0]*fr[5]-7245.0*coeff[0]*fl[5]-15750.0*coeff[0]*fc[5]+4364.768035073569*coeff[0]*fr[3]-4364.768035073569*coeff[0]*fl[3])*rdx2Sq; 
  out[6] += 0.0625*(563.489130329947*coeff[0]*fr[17]+563.489130329947*coeff[0]*fl[17]-1126.978260659894*coeff[0]*fc[17]-545.5960043841961*coeff[0]*fr[10]+545.5960043841961*coeff[0]*fl[10]+315.0*coeff[0]*fr[6]+315.0*coeff[0]*fl[6]-630.0*coeff[0]*fc[6])*rdx2Sq; 
  out[7] += -0.0078125*(405.0*coeff[0]*fr[7]+405.0*coeff[0]*fl[7]+18090.0*coeff[0]*fc[7]+1568.558255214003*coeff[0]*fr[1]-1568.558255214003*coeff[0]*fl[1]-1609.968943799849*coeff[0]*fr[0]-1609.968943799849*coeff[0]*fl[0]+3219.937887599698*coeff[0]*fc[0])*rdx2Sq; 
  out[8] += -0.0625*(545.5960043841964*coeff[0]*fr[12]-545.5960043841964*coeff[0]*fl[12]-315.0*coeff[0]*fr[8]-315.0*coeff[0]*fl[8]+630.0*coeff[0]*fc[8])*rdx2Sq; 
  out[9] += -0.0625*(545.5960043841964*coeff[0]*fr[15]-545.5960043841964*coeff[0]*fl[15]-315.0*coeff[0]*fr[9]-315.0*coeff[0]*fl[9]+630.0*coeff[0]*fc[9])*rdx2Sq; 
  out[10] += 0.0078125*(6587.944671898812*coeff[0]*fr[17]-6587.944671898812*coeff[0]*fl[17]-7245.0*coeff[0]*fr[10]-7245.0*coeff[0]*fl[10]-15750.0*coeff[0]*fc[10]+4364.768035073569*coeff[0]*fr[6]-4364.768035073569*coeff[0]*fl[6])*rdx2Sq; 
  out[11] += -0.0078125*(405.0*coeff[0]*fr[11]+405.0*coeff[0]*fl[11]+18090.0*coeff[0]*fc[11]+1568.558255214004*coeff[0]*fr[4]-1568.558255214004*coeff[0]*fl[4]-1609.968943799848*coeff[0]*fr[2]-1609.968943799848*coeff[0]*fl[2]+3219.937887599697*coeff[0]*fc[2])*rdx2Sq; 
  out[12] += -0.0078125*(7245.0*coeff[0]*fr[12]+7245.0*coeff[0]*fl[12]+15750.0*coeff[0]*fc[12]-4364.768035073571*coeff[0]*fr[8]+4364.768035073571*coeff[0]*fl[8])*rdx2Sq; 
  out[13] += -0.0078125*(405.0*coeff[0]*fr[13]+405.0*coeff[0]*fl[13]+18090.0*coeff[0]*fc[13]+1568.558255214004*coeff[0]*fr[5]-1568.558255214004*coeff[0]*fl[5]-1609.968943799848*coeff[0]*fr[3]-1609.968943799848*coeff[0]*fl[3]+3219.937887599697*coeff[0]*fc[3])*rdx2Sq; 
  out[14] += -0.0625*(545.5960043841964*coeff[0]*fr[18]-545.5960043841964*coeff[0]*fl[18]-315.0*coeff[0]*fr[14]-315.0*coeff[0]*fl[14]+630.0*coeff[0]*fc[14])*rdx2Sq; 
  out[15] += -0.0078125*(7245.0*coeff[0]*fr[15]+7245.0*coeff[0]*fl[15]+15750.0*coeff[0]*fc[15]-4364.768035073571*coeff[0]*fr[9]+4364.768035073571*coeff[0]*fl[9])*rdx2Sq; 
  out[16] += -0.0625*(545.5960043841964*coeff[0]*fr[19]-545.5960043841964*coeff[0]*fl[19]-315.0*coeff[0]*fr[16]-315.0*coeff[0]*fl[16]+630.0*coeff[0]*fc[16])*rdx2Sq; 
  out[17] += -0.0078125*(405.0*coeff[0]*fr[17]+405.0*coeff[0]*fl[17]+18090.0*coeff[0]*fc[17]+1568.558255214003*coeff[0]*fr[10]-1568.558255214003*coeff[0]*fl[10]-1609.968943799849*coeff[0]*fr[6]-1609.968943799849*coeff[0]*fl[6]+3219.937887599698*coeff[0]*fc[6])*rdx2Sq; 
  out[18] += -0.0078125*(7245.0*coeff[0]*fr[18]+7245.0*coeff[0]*fl[18]+15750.0*coeff[0]*fc[18]-4364.768035073571*coeff[0]*fr[14]+4364.768035073571*coeff[0]*fl[14])*rdx2Sq; 
  out[19] += -0.0078125*(7245.0*coeff[0]*fr[19]+7245.0*coeff[0]*fl[19]+15750.0*coeff[0]*fc[19]-4364.768035073571*coeff[0]*fr[16]+4364.768035073571*coeff[0]*fl[16])*rdx2Sq; 

  return 0.;

}

