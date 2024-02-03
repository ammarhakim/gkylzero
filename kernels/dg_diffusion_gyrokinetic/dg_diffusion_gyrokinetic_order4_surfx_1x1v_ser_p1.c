#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
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

  double fl[6];
  fl[0] = 0.7071067811865476*(jacobgeo_inv[1]*ql[1]+jacobgeo_inv[0]*ql[0]); 
  fl[1] = 0.7071067811865476*(jacobgeo_inv[0]*ql[1]+ql[0]*jacobgeo_inv[1]); 
  fl[2] = 0.7071067811865476*(jacobgeo_inv[1]*ql[3]+jacobgeo_inv[0]*ql[2]); 
  fl[3] = 0.7071067811865476*(jacobgeo_inv[0]*ql[3]+jacobgeo_inv[1]*ql[2]); 
  fl[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*ql[5]+21.21320343559643*jacobgeo_inv[0]*ql[4]); 
  fl[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*ql[5]+21.21320343559643*jacobgeo_inv[1]*ql[4]); 

  double fc[6];
  fc[0] = 0.7071067811865476*(jacobgeo_inv[1]*qc[1]+jacobgeo_inv[0]*qc[0]); 
  fc[1] = 0.7071067811865476*(jacobgeo_inv[0]*qc[1]+qc[0]*jacobgeo_inv[1]); 
  fc[2] = 0.7071067811865476*(jacobgeo_inv[1]*qc[3]+jacobgeo_inv[0]*qc[2]); 
  fc[3] = 0.7071067811865476*(jacobgeo_inv[0]*qc[3]+jacobgeo_inv[1]*qc[2]); 
  fc[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qc[5]+21.21320343559643*jacobgeo_inv[0]*qc[4]); 
  fc[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*qc[5]+21.21320343559643*jacobgeo_inv[1]*qc[4]); 

  double fr[6];
  fr[0] = 0.7071067811865476*(jacobgeo_inv[1]*qr[1]+jacobgeo_inv[0]*qr[0]); 
  fr[1] = 0.7071067811865476*(jacobgeo_inv[0]*qr[1]+qr[0]*jacobgeo_inv[1]); 
  fr[2] = 0.7071067811865476*(jacobgeo_inv[1]*qr[3]+jacobgeo_inv[0]*qr[2]); 
  fr[3] = 0.7071067811865476*(jacobgeo_inv[0]*qr[3]+jacobgeo_inv[1]*qr[2]); 
  fr[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qr[5]+21.21320343559643*jacobgeo_inv[0]*qr[4]); 
  fr[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[0]*qr[5]+21.21320343559643*jacobgeo_inv[1]*qr[4]); 

  out[0] += -0.0625*(25.98076211353316*coeff[0]*fr[1]-25.98076211353316*coeff[0]*fl[1]-15.0*coeff[0]*fr[0]-15.0*coeff[0]*fl[0]+30.0*coeff[0]*fc[0])*rdx2Sq; 
  out[1] += -0.0625*(33.0*coeff[0]*fr[1]+33.0*coeff[0]*fl[1]+114.0*coeff[0]*fc[1]-25.98076211353316*coeff[0]*fr[0]+25.98076211353316*coeff[0]*fl[0])*rdx2Sq; 
  out[2] += -0.0625*(25.98076211353316*coeff[0]*fr[3]-25.98076211353316*coeff[0]*fl[3]-15.0*coeff[0]*fr[2]-15.0*coeff[0]*fl[2]+30.0*coeff[0]*fc[2])*rdx2Sq; 
  out[3] += -0.0625*(33.0*coeff[0]*fr[3]+33.0*coeff[0]*fl[3]+114.0*coeff[0]*fc[3]-25.98076211353316*coeff[0]*fr[2]+25.98076211353316*coeff[0]*fl[2])*rdx2Sq; 
  out[4] += -0.0625*(25.98076211353316*coeff[0]*fr[5]-25.98076211353316*coeff[0]*fl[5]-15.0*coeff[0]*fr[4]-15.0*coeff[0]*fl[4]+30.0*coeff[0]*fc[4])*rdx2Sq; 
  out[5] += -0.0625*(33.0*coeff[0]*fr[5]+33.0*coeff[0]*fl[5]+114.0*coeff[0]*fc[5]-25.98076211353316*coeff[0]*fr[4]+25.98076211353316*coeff[0]*fl[4])*rdx2Sq; 

  return 0.;

}

