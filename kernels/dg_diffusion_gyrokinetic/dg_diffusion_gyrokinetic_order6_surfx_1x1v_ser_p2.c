#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  double fl[8];
  fl[0] = 0.7071067811865476*(jacobgeo_inv[2]*ql[4]+jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*ql[4]+ql[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1])); 
  fl[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*ql[6]+21.21320343559643*(jacobgeo_inv[1]*ql[3]+jacobgeo_inv[0]*ql[2])); 
  fl[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*ql[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*ql[3]+21.21320343559643*jacobgeo_inv[1]*ql[2]); 
  fl[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*ql[4]+49.49747468305833*ql[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*ql[1]); 
  fl[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*ql[7]+21.21320343559643*jacobgeo_inv[0]*ql[5]); 
  fl[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*ql[6]+132.815661727072*jacobgeo_inv[1]*ql[3]+148.492424049175*jacobgeo_inv[2]*ql[2]); 
  fl[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*ql[7]+21.21320343559643*jacobgeo_inv[1]*ql[5]); 

  double fc[8];
  fc[0] = 0.7071067811865476*(jacobgeo_inv[2]*qc[4]+jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qc[4]+qc[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1])); 
  fc[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qc[6]+21.21320343559643*(jacobgeo_inv[1]*qc[3]+jacobgeo_inv[0]*qc[2])); 
  fc[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qc[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qc[3]+21.21320343559643*jacobgeo_inv[1]*qc[2]); 
  fc[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qc[4]+49.49747468305833*qc[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qc[1]); 
  fc[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qc[7]+21.21320343559643*jacobgeo_inv[0]*qc[5]); 
  fc[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*qc[6]+132.815661727072*jacobgeo_inv[1]*qc[3]+148.492424049175*jacobgeo_inv[2]*qc[2]); 
  fc[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qc[7]+21.21320343559643*jacobgeo_inv[1]*qc[5]); 

  double fr[8];
  fr[0] = 0.7071067811865476*(jacobgeo_inv[2]*qr[4]+jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qr[4]+qr[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1])); 
  fr[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qr[6]+21.21320343559643*(jacobgeo_inv[1]*qr[3]+jacobgeo_inv[0]*qr[2])); 
  fr[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qr[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qr[3]+21.21320343559643*jacobgeo_inv[1]*qr[2]); 
  fr[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qr[4]+49.49747468305833*qr[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qr[1]); 
  fr[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qr[7]+21.21320343559643*jacobgeo_inv[0]*qr[5]); 
  fr[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*qr[6]+132.815661727072*jacobgeo_inv[1]*qr[3]+148.492424049175*jacobgeo_inv[2]*qr[2]); 
  fr[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qr[7]+21.21320343559643*jacobgeo_inv[1]*qr[5]); 

  out[0] += 0.0625*(563.489130329947*coeff[0]*fr[4]+563.489130329947*coeff[0]*fl[4]-1126.978260659894*coeff[0]*fc[4]-545.5960043841961*coeff[0]*fr[1]+545.5960043841961*coeff[0]*fl[1]+315.0*coeff[0]*fr[0]+315.0*coeff[0]*fl[0]-630.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += 0.0078125*(6587.944671898812*coeff[0]*fr[4]-6587.944671898812*coeff[0]*fl[4]-7245.0*coeff[0]*fr[1]-7245.0*coeff[0]*fl[1]-15750.0*coeff[0]*fc[1]+4364.768035073569*coeff[0]*fr[0]-4364.768035073569*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += 0.0625*(563.4891303299469*coeff[0]*fr[6]+563.4891303299469*coeff[0]*fl[6]-1126.978260659894*coeff[0]*fc[6]-545.5960043841961*coeff[0]*fr[3]+545.5960043841961*coeff[0]*fl[3]+315.0*coeff[0]*fr[2]+315.0*coeff[0]*fl[2]-630.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += 0.0078125*(6587.944671898817*coeff[0]*fr[6]-6587.944671898817*coeff[0]*fl[6]-7245.0*coeff[0]*fr[3]-7245.0*coeff[0]*fl[3]-15750.0*coeff[0]*fc[3]+4364.768035073569*coeff[0]*fr[2]-4364.768035073569*coeff[0]*fl[2])*rdx2Sq; 
  out[4] += -0.0078125*(405.0*coeff[0]*fr[4]+405.0*coeff[0]*fl[4]+18090.0*coeff[0]*fc[4]+1568.558255214003*coeff[0]*fr[1]-1568.558255214003*coeff[0]*fl[1]-1609.968943799849*coeff[0]*fr[0]-1609.968943799849*coeff[0]*fl[0]+3219.937887599698*coeff[0]*fc[0])*rdx2Sq; 
  out[5] += -0.0625*(545.5960043841964*coeff[0]*fr[7]-545.5960043841964*coeff[0]*fl[7]-315.0*coeff[0]*fr[5]-315.0*coeff[0]*fl[5]+630.0*coeff[0]*fc[5])*rdx2Sq; 
  out[6] += -0.0078125*(405.0*coeff[0]*fr[6]+405.0*coeff[0]*fl[6]+18090.0*coeff[0]*fc[6]+1568.558255214004*coeff[0]*fr[3]-1568.558255214004*coeff[0]*fl[3]-1609.968943799848*coeff[0]*fr[2]-1609.968943799848*coeff[0]*fl[2]+3219.937887599697*coeff[0]*fc[2])*rdx2Sq; 
  out[7] += -0.0078125*(7245.0*coeff[0]*fr[7]+7245.0*coeff[0]*fl[7]+15750.0*coeff[0]*fc[7]-4364.768035073571*coeff[0]*fr[5]+4364.768035073571*coeff[0]*fl[5])*rdx2Sq; 

  return 0.;

}

