#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_Bstar_Bmag_2x_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, double* GKYL_RESTRICT Bstar_Bmag) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species q_ and m_.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // Bstar_Bmag: output volume expansion of B*/Bmag.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double *BstarXdBmag = &Bstar_Bmag[0]; 
  double *BstarYdBmag = &Bstar_Bmag[8]; 

  BstarYdBmag[0] = (((-1.936491673103709*jacobtot_inv[1]*b_z[4])-0.8660254037844386*jacobtot_inv[0]*b_z[1])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[1] = ((b_z[4]*((-1.732050807568877*jacobtot_inv[4])-1.936491673103709*jacobtot_inv[0])-0.8660254037844386*b_z[1]*jacobtot_inv[1])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[2] = (((-1.118033988749895*jacobtot_inv[1]*b_z[4])-0.5*jacobtot_inv[0]*b_z[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[3] = ((b_z[4]*((-1.0*jacobtot_inv[4])-1.118033988749895*jacobtot_inv[0])-0.5*b_z[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[4] = (((-0.8660254037844386*b_z[1]*jacobtot_inv[4])-1.732050807568877*jacobtot_inv[1]*b_z[4])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[6] = (((-0.5000000000000001*b_z[1]*jacobtot_inv[4])-1.0*jacobtot_inv[1]*b_z[4])*m_*rdx2)/(q_*rdvpar2); 

} 
