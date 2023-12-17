#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_Bstar_Bmag_3x_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, double* GKYL_RESTRICT Bstar_Bmag) 
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double *BstarXdBmag = &Bstar_Bmag[0]; 
  double *BstarYdBmag = &Bstar_Bmag[8]; 
  double *BstarZdBmag = &Bstar_Bmag[16]; 
  BstarXdBmag[0] = -(0.6123724356957944*(jacobtot_inv[1]*b_y[5]+jacobtot_inv[0]*b_y[3])*m_*rdz2*wvpar)/q_; 
  BstarXdBmag[1] = -(0.6123724356957944*(jacobtot_inv[0]*b_y[5]+jacobtot_inv[1]*b_y[3])*m_*rdz2*wvpar)/q_; 
  BstarXdBmag[2] = -(0.6123724356957944*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])*m_*rdz2*wvpar)/q_; 
  BstarXdBmag[3] = -(0.3535533905932737*(jacobtot_inv[1]*b_y[5]+jacobtot_inv[0]*b_y[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[4] = -(0.6123724356957944*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])*m_*rdz2*wvpar)/q_; 
  BstarXdBmag[5] = -(0.3535533905932737*(jacobtot_inv[0]*b_y[5]+jacobtot_inv[1]*b_y[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[6] = -(0.3535533905932737*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[7] = -(0.3535533905932737*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])*m_*rdz2)/(q_*rdvpar2); 

  BstarYdBmag[0] = (m_*(0.6123724356957944*(jacobtot_inv[1]*b_x[5]+jacobtot_inv[0]*b_x[3])*rdz2-0.6123724356957944*(jacobtot_inv[3]*b_z[5]+jacobtot_inv[0]*b_z[1])*rdx2)*wvpar)/q_; 
  BstarYdBmag[1] = (m_*(0.6123724356957944*(jacobtot_inv[0]*b_x[5]+jacobtot_inv[1]*b_x[3])*rdz2-0.6123724356957944*(b_z[5]*jacobtot_inv[5]+b_z[1]*jacobtot_inv[1])*rdx2)*wvpar)/q_; 
  BstarYdBmag[2] = (m_*(0.6123724356957944*(b_x[5]*jacobtot_inv[5]+b_x[3]*jacobtot_inv[3])*rdz2-0.6123724356957944*(jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3])*rdx2)*wvpar)/q_; 
  BstarYdBmag[3] = (m_*(0.3535533905932737*(jacobtot_inv[1]*b_x[5]+jacobtot_inv[0]*b_x[3])*rdz2-0.3535533905932737*(jacobtot_inv[3]*b_z[5]+jacobtot_inv[0]*b_z[1])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[4] = (m_*(0.6123724356957944*(b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])*rdz2-0.6123724356957944*(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])*rdx2)*wvpar)/q_; 
  BstarYdBmag[5] = (m_*(0.3535533905932737*(jacobtot_inv[0]*b_x[5]+jacobtot_inv[1]*b_x[3])*rdz2-0.3535533905932737*(b_z[5]*jacobtot_inv[5]+b_z[1]*jacobtot_inv[1])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[6] = (m_*(0.3535533905932737*(b_x[5]*jacobtot_inv[5]+b_x[3]*jacobtot_inv[3])*rdz2-0.3535533905932737*(jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[7] = (m_*(0.3535533905932737*(b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])*rdz2-0.3535533905932737*(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])*rdx2))/(q_*rdvpar2); 

  BstarZdBmag[0] = (0.6123724356957944*(jacobtot_inv[3]*b_y[5]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar)/q_+0.3535533905932737*(cmag[5]*jacobtot_inv[5]+cmag[3]*jacobtot_inv[3]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0]); 
  BstarZdBmag[1] = (0.6123724356957944*(b_y[5]*jacobtot_inv[5]+b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar)/q_+0.3535533905932737*(cmag[3]*jacobtot_inv[5]+jacobtot_inv[3]*cmag[5]+cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1]); 
  BstarZdBmag[2] = (0.6123724356957944*(jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3])*m_*rdx2*wvpar)/q_+0.3535533905932737*(cmag[1]*jacobtot_inv[5]+jacobtot_inv[1]*cmag[5]+cmag[0]*jacobtot_inv[3]+jacobtot_inv[0]*cmag[3]); 
  BstarZdBmag[3] = (0.3535533905932737*(jacobtot_inv[3]*b_y[5]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (0.6123724356957944*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5])*m_*rdx2*wvpar)/q_+0.3535533905932737*(cmag[0]*jacobtot_inv[5]+jacobtot_inv[0]*cmag[5]+cmag[1]*jacobtot_inv[3]+jacobtot_inv[1]*cmag[3]); 
  BstarZdBmag[5] = (0.3535533905932737*(b_y[5]*jacobtot_inv[5]+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[6] = (0.3535533905932737*(jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.3535533905932737*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5])*m_*rdx2)/(q_*rdvpar2); 

} 
