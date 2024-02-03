#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  out[0] += 0.0625*(107.3312629199899*coeff[0]*fr[4]+107.3312629199899*coeff[0]*fl[4]-214.6625258399798*coeff[0]*fc[4]-129.9038105676658*coeff[0]*fr[1]+129.9038105676658*coeff[0]*fl[1]+75.0*coeff[0]*fr[0]+75.0*coeff[0]*fl[0]-150.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += 0.03125*(290.4737509655563*coeff[0]*fr[4]-290.4737509655563*coeff[0]*fl[4]-405.0*coeff[0]*fr[1]-405.0*coeff[0]*fl[1]-990.0*coeff[0]*fc[1]+259.8076211353315*coeff[0]*fr[0]-259.8076211353315*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += 0.0625*(107.3312629199899*coeff[0]*fr[6]+107.3312629199899*coeff[0]*fl[6]-214.6625258399798*coeff[0]*fc[6]-129.9038105676658*coeff[0]*fr[3]+129.9038105676658*coeff[0]*fl[3]+75.0*coeff[0]*fr[2]+75.0*coeff[0]*fl[2]-150.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += 0.03125*(290.4737509655563*coeff[0]*fr[6]-290.4737509655563*coeff[0]*fl[6]-405.0*coeff[0]*fr[3]-405.0*coeff[0]*fl[3]-990.0*coeff[0]*fc[3]+259.8076211353315*coeff[0]*fr[2]-259.8076211353315*coeff[0]*fl[2])*rdx2Sq; 
  out[4] += 0.03125*(21.0*coeff[0]*fr[4]+21.0*coeff[0]*fl[4]-1302.0*coeff[0]*fc[4]-151.0463505020892*coeff[0]*fr[1]+151.0463505020892*coeff[0]*fl[1]+134.1640786499874*coeff[0]*fr[0]+134.1640786499874*coeff[0]*fl[0]-268.3281572999748*coeff[0]*fc[0])*rdx2Sq; 
  out[5] += -0.0625*(129.9038105676658*coeff[0]*fr[7]-129.9038105676658*coeff[0]*fl[7]-75.0*coeff[0]*fr[5]-75.0*coeff[0]*fl[5]+150.0*coeff[0]*fc[5])*rdx2Sq; 
  out[6] += 0.03125*(21.0*coeff[0]*fr[6]+21.0*coeff[0]*fl[6]-1302.0*coeff[0]*fc[6]-151.0463505020893*coeff[0]*fr[3]+151.0463505020893*coeff[0]*fl[3]+134.1640786499874*coeff[0]*fr[2]+134.1640786499874*coeff[0]*fl[2]-268.3281572999747*coeff[0]*fc[2])*rdx2Sq; 
  out[7] += -0.03125*(405.0*coeff[0]*fr[7]+405.0*coeff[0]*fl[7]+990.0*coeff[0]*fc[7]-259.8076211353317*coeff[0]*fr[5]+259.8076211353317*coeff[0]*fl[5])*rdx2Sq; 

  return 0.;

}

