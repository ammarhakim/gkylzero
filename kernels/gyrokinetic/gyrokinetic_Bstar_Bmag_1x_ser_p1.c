#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_Bstar_Bmag_1x_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, double* GKYL_RESTRICT Bstar_Bmag) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species q_ and m_.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // Bstar_Bmag: output volume expansion of B*/Bmag.


  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double *BstarZdBmag = &Bstar_Bmag[0]; 
  BstarZdBmag[0] = 0.7071067811865475*(cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0]); 
  BstarZdBmag[1] = 0.7071067811865475*(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1]); 

} 
