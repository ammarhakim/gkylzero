#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void gyrokinetic_prim_vars_vtSq_3x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: vtSq = 1/vdim(m2/m0 - upar^2). 
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[8]; 
  const double *m2 = &moms[16]; 
 
  double *vtSq = &prim_vars[0]; 
 
  double m0Sq[8] = {0.0}; 

  double m0Sq_inv[8] = {0.0}; 

  double m1Sq[8] = {0.0}; 
  double m2m0[8] = {0.0}; 

  binop_mul_3d_ser_p1(m0, m0, m0Sq); 
  binop_mul_3d_ser_p1(m1, m1, m1Sq); 
  binop_mul_3d_ser_p1(m0, m2, m2m0); 
  ser_3x_p1_inv(m0Sq, m0Sq_inv); 

  vtSq[0] = 0.1178511301977579*m0Sq_inv[7]*m2m0[7]-0.1178511301977579*m0Sq_inv[7]*m1Sq[7]+0.1178511301977579*m0Sq_inv[6]*m2m0[6]-0.1178511301977579*m0Sq_inv[6]*m1Sq[6]+0.1178511301977579*m0Sq_inv[5]*m2m0[5]-0.1178511301977579*m0Sq_inv[5]*m1Sq[5]+0.1178511301977579*m0Sq_inv[4]*m2m0[4]-0.1178511301977579*m0Sq_inv[4]*m1Sq[4]+0.1178511301977579*m0Sq_inv[3]*m2m0[3]-0.1178511301977579*m0Sq_inv[3]*m1Sq[3]+0.1178511301977579*m0Sq_inv[2]*m2m0[2]-0.1178511301977579*m0Sq_inv[2]*m1Sq[2]+0.1178511301977579*m0Sq_inv[1]*m2m0[1]-0.1178511301977579*m0Sq_inv[1]*m1Sq[1]+0.1178511301977579*m0Sq_inv[0]*m2m0[0]-0.1178511301977579*m0Sq_inv[0]*m1Sq[0]; 
  vtSq[1] = 0.1178511301977579*m0Sq_inv[6]*m2m0[7]-0.1178511301977579*m0Sq_inv[6]*m1Sq[7]+0.1178511301977579*m2m0[6]*m0Sq_inv[7]-0.1178511301977579*m1Sq[6]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[3]*m2m0[5]-0.1178511301977579*m0Sq_inv[3]*m1Sq[5]+0.1178511301977579*m2m0[3]*m0Sq_inv[5]-0.1178511301977579*m1Sq[3]*m0Sq_inv[5]+0.1178511301977579*m0Sq_inv[2]*m2m0[4]-0.1178511301977579*m0Sq_inv[2]*m1Sq[4]+0.1178511301977579*m2m0[2]*m0Sq_inv[4]-0.1178511301977579*m1Sq[2]*m0Sq_inv[4]+0.1178511301977579*m0Sq_inv[0]*m2m0[1]-0.1178511301977579*m0Sq_inv[0]*m1Sq[1]+0.1178511301977579*m2m0[0]*m0Sq_inv[1]-0.1178511301977579*m1Sq[0]*m0Sq_inv[1]; 
  vtSq[2] = 0.1178511301977579*m0Sq_inv[5]*m2m0[7]-0.1178511301977579*m0Sq_inv[5]*m1Sq[7]+0.1178511301977579*m2m0[5]*m0Sq_inv[7]-0.1178511301977579*m1Sq[5]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[3]*m2m0[6]-0.1178511301977579*m0Sq_inv[3]*m1Sq[6]+0.1178511301977579*m2m0[3]*m0Sq_inv[6]-0.1178511301977579*m1Sq[3]*m0Sq_inv[6]+0.1178511301977579*m0Sq_inv[1]*m2m0[4]-0.1178511301977579*m0Sq_inv[1]*m1Sq[4]+0.1178511301977579*m2m0[1]*m0Sq_inv[4]-0.1178511301977579*m1Sq[1]*m0Sq_inv[4]+0.1178511301977579*m0Sq_inv[0]*m2m0[2]-0.1178511301977579*m0Sq_inv[0]*m1Sq[2]+0.1178511301977579*m2m0[0]*m0Sq_inv[2]-0.1178511301977579*m1Sq[0]*m0Sq_inv[2]; 
  vtSq[3] = 0.1178511301977579*m0Sq_inv[4]*m2m0[7]-0.1178511301977579*m0Sq_inv[4]*m1Sq[7]+0.1178511301977579*m2m0[4]*m0Sq_inv[7]-0.1178511301977579*m1Sq[4]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[2]*m2m0[6]-0.1178511301977579*m0Sq_inv[2]*m1Sq[6]+0.1178511301977579*m2m0[2]*m0Sq_inv[6]-0.1178511301977579*m1Sq[2]*m0Sq_inv[6]+0.1178511301977579*m0Sq_inv[1]*m2m0[5]-0.1178511301977579*m0Sq_inv[1]*m1Sq[5]+0.1178511301977579*m2m0[1]*m0Sq_inv[5]-0.1178511301977579*m1Sq[1]*m0Sq_inv[5]+0.1178511301977579*m0Sq_inv[0]*m2m0[3]-0.1178511301977579*m0Sq_inv[0]*m1Sq[3]+0.1178511301977579*m2m0[0]*m0Sq_inv[3]-0.1178511301977579*m1Sq[0]*m0Sq_inv[3]; 
  vtSq[4] = 0.1178511301977579*m0Sq_inv[3]*m2m0[7]-0.1178511301977579*m0Sq_inv[3]*m1Sq[7]+0.1178511301977579*m2m0[3]*m0Sq_inv[7]-0.1178511301977579*m1Sq[3]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[5]*m2m0[6]-0.1178511301977579*m0Sq_inv[5]*m1Sq[6]+0.1178511301977579*m2m0[5]*m0Sq_inv[6]-0.1178511301977579*m1Sq[5]*m0Sq_inv[6]+0.1178511301977579*m0Sq_inv[0]*m2m0[4]-0.1178511301977579*m0Sq_inv[0]*m1Sq[4]+0.1178511301977579*m2m0[0]*m0Sq_inv[4]-0.1178511301977579*m1Sq[0]*m0Sq_inv[4]+0.1178511301977579*m0Sq_inv[1]*m2m0[2]-0.1178511301977579*m0Sq_inv[1]*m1Sq[2]+0.1178511301977579*m2m0[1]*m0Sq_inv[2]-0.1178511301977579*m1Sq[1]*m0Sq_inv[2]; 
  vtSq[5] = 0.1178511301977579*m0Sq_inv[2]*m2m0[7]-0.1178511301977579*m0Sq_inv[2]*m1Sq[7]+0.1178511301977579*m2m0[2]*m0Sq_inv[7]-0.1178511301977579*m1Sq[2]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[4]*m2m0[6]-0.1178511301977579*m0Sq_inv[4]*m1Sq[6]+0.1178511301977579*m2m0[4]*m0Sq_inv[6]-0.1178511301977579*m1Sq[4]*m0Sq_inv[6]+0.1178511301977579*m0Sq_inv[0]*m2m0[5]-0.1178511301977579*m0Sq_inv[0]*m1Sq[5]+0.1178511301977579*m2m0[0]*m0Sq_inv[5]-0.1178511301977579*m1Sq[0]*m0Sq_inv[5]+0.1178511301977579*m0Sq_inv[1]*m2m0[3]-0.1178511301977579*m0Sq_inv[1]*m1Sq[3]+0.1178511301977579*m2m0[1]*m0Sq_inv[3]-0.1178511301977579*m1Sq[1]*m0Sq_inv[3]; 
  vtSq[6] = 0.1178511301977579*m0Sq_inv[1]*m2m0[7]-0.1178511301977579*m0Sq_inv[1]*m1Sq[7]+0.1178511301977579*m2m0[1]*m0Sq_inv[7]-0.1178511301977579*m1Sq[1]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[0]*m2m0[6]-0.1178511301977579*m0Sq_inv[0]*m1Sq[6]+0.1178511301977579*m2m0[0]*m0Sq_inv[6]-0.1178511301977579*m1Sq[0]*m0Sq_inv[6]+0.1178511301977579*m0Sq_inv[4]*m2m0[5]-0.1178511301977579*m0Sq_inv[4]*m1Sq[5]+0.1178511301977579*m2m0[4]*m0Sq_inv[5]-0.1178511301977579*m1Sq[4]*m0Sq_inv[5]+0.1178511301977579*m0Sq_inv[2]*m2m0[3]-0.1178511301977579*m0Sq_inv[2]*m1Sq[3]+0.1178511301977579*m2m0[2]*m0Sq_inv[3]-0.1178511301977579*m1Sq[2]*m0Sq_inv[3]; 
  vtSq[7] = 0.1178511301977579*m0Sq_inv[0]*m2m0[7]-0.1178511301977579*m0Sq_inv[0]*m1Sq[7]+0.1178511301977579*m2m0[0]*m0Sq_inv[7]-0.1178511301977579*m1Sq[0]*m0Sq_inv[7]+0.1178511301977579*m0Sq_inv[1]*m2m0[6]-0.1178511301977579*m0Sq_inv[1]*m1Sq[6]+0.1178511301977579*m2m0[1]*m0Sq_inv[6]-0.1178511301977579*m1Sq[1]*m0Sq_inv[6]+0.1178511301977579*m0Sq_inv[2]*m2m0[5]-0.1178511301977579*m0Sq_inv[2]*m1Sq[5]+0.1178511301977579*m2m0[2]*m0Sq_inv[5]-0.1178511301977579*m1Sq[2]*m0Sq_inv[5]+0.1178511301977579*m0Sq_inv[3]*m2m0[4]-0.1178511301977579*m0Sq_inv[3]*m1Sq[4]+0.1178511301977579*m2m0[3]*m0Sq_inv[4]-0.1178511301977579*m1Sq[3]*m0Sq_inv[4]; 

} 
 