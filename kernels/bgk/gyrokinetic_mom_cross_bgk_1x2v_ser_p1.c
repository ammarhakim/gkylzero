#include <gkyl_mom_cross_bgk_gyrokinetic_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
 
GKYL_CU_DH void gyrokinetic_mom_cross_bgk_1x2v_ser_p1(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double *nu_sr, const double *nu_rs, double *moms_cross) 
{ 
  // m_:             mass. 
  // moms_:          moments of the distribution function. 
  // moms_cross:     cross moments. 
 
  const double ms = m_self; 
  const double mr = m_other; 
  const double *m0s = &moms_self[0]; 
  const double *m1s = &moms_self[2]; 
  const double *m2s = &moms_self[4]; 
  const double *m0r = &moms_other[0]; 
  const double *m1r = &moms_other[2]; 
  const double *m2r = &moms_other[4]; 
 
  double *m0sr = &moms_cross[0]; 
  double *m1sr = &moms_cross[2]; 
  double *m2sr = &moms_cross[4]; 
 
  double msM0sNusr[2] = {0.0}; 
  double mrM0rNurs[2] = {0.0}; 
  double mM0Nu[2] = {0.0}; 
  double mM0Nu_inv[2] = {0.0}; 
  double alphaE[2] = {0.0}; 

  double nu_sr_inv[2] = {0.0}; 
  double m0sm0r[2] = {0.0}; 
  double m0sm0r_inv[2] = {0.0}; 
  double m0sm1r[2] = {0.0}; 
  double m0rm1s[2] = {0.0}; 
  double cM1sr_temp[2] = {0.0}; 
  double cM1sr[2] = {0.0}; 

  double m0sm2r[2] = {0.0}; 
  double m0rm2s[2] = {0.0}; 
  double m1sm1r[2] = {0.0}; 
  double cM2sr_temp[2] = {0.0}; 
  double cM2sr[2] = {0.0}; 

  binop_mul_1d_ser_p1(m0s, nu_sr, msM0sNusr); 
  binop_mul_1d_ser_p1(m0r, nu_rs, mrM0rNurs); 
  msM0sNusr[0] = ms * msM0sNusr[0]; 
  mrM0rNurs[0] = mr * mrM0rNurs[0]; 
  mM0Nu[0] = msM0sNusr[0] + mrM0rNurs[0]; 
  msM0sNusr[1] = ms * msM0sNusr[1]; 
  mrM0rNurs[1] = mr * mrM0rNurs[1]; 
  mM0Nu[1] = msM0sNusr[1] + mrM0rNurs[1]; 
  ser_1x_p1_inv(mM0Nu, mM0Nu_inv); 
  binop_mul_1d_ser_p1(msM0sNusr, mrM0rNurs, alphaE); 
  binop_mul_1d_ser_p1(alphaE, mM0Nu_inv, alphaE); 
  alphaE[0] = alphaE[0] * 2 * (1+beta) / (ms+mr); 
  alphaE[1] = alphaE[1] * 2 * (1+beta) / (ms+mr); 

  m0sr[0] = m0s[0]; 
  m0sr[1] = m0s[1]; 
 
  ser_1x_p1_inv(nu_sr, nu_sr_inv); 
  binop_mul_1d_ser_p1(m0s, m0r, m0sm0r); 
  ser_1x_p1_inv(m0sm0r, m0sm0r_inv); 
  binop_mul_1d_ser_p1(m0s, m1r, m0sm1r); 
  binop_mul_1d_ser_p1(m0r, m1s, m0rm1s); 
  cM1sr_temp[0] = (0.3535533905932737*m0sm0r_inv[1]*m0sm1r[1]*mr)/ms-(0.3535533905932737*m0rm1s[1]*m0sm0r_inv[1]*mr)/ms+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[0]*mr)/ms-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[0]*mr)/ms+0.3535533905932737*m0sm0r_inv[1]*m0sm1r[1]-0.3535533905932737*m0rm1s[1]*m0sm0r_inv[1]+0.3535533905932737*m0sm0r_inv[0]*m0sm1r[0]-0.3535533905932737*m0rm1s[0]*m0sm0r_inv[0]; 
  cM1sr_temp[1] = (0.3535533905932737*m0sm0r_inv[0]*m0sm1r[1]*mr)/ms+(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[1]*mr)/ms-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[1]*mr)/ms-(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[1]*mr)/ms+0.3535533905932737*m0sm0r_inv[0]*m0sm1r[1]+0.3535533905932737*m0sm1r[0]*m0sm0r_inv[1]-0.3535533905932737*m0rm1s[0]*m0sm0r_inv[1]-0.3535533905932737*m0sm0r_inv[0]*m0rm1s[1]; 
  binop_mul_1d_ser_p1(nu_sr_inv, cM1sr_temp, cM1sr_temp); 
  binop_mul_1d_ser_p1(alphaE, cM1sr_temp, cM1sr); 
  m1sr[0] = m1s[0] + cM1sr[0]; 
  m1sr[1] = m1s[1] + cM1sr[1]; 
 
  binop_mul_1d_ser_p1(m0s, m2r, m0sm2r); 
  binop_mul_1d_ser_p1(m0r, m2s, m0rm2s); 
  binop_mul_1d_ser_p1(m1s, m1r, m1sm1r); 
  cM2sr_temp[0] = (-(0.7071067811865475*m0sm0r_inv[1]*m1sm1r[1]*mr)/ms)+(0.7071067811865475*m0sm0r_inv[1]*m0sm2r[1]*mr)/ms-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[0]*mr)/ms+(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[0]*mr)/ms+0.7071067811865475*m0sm0r_inv[1]*m1sm1r[1]-0.7071067811865475*m0rm2s[1]*m0sm0r_inv[1]+0.7071067811865475*m0sm0r_inv[0]*m1sm1r[0]-0.7071067811865475*m0rm2s[0]*m0sm0r_inv[0]; 
  cM2sr_temp[1] = (-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[1]*mr)/ms)+(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[1]*mr)/ms-(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[1]*mr)/ms+(0.7071067811865475*m0sm2r[0]*m0sm0r_inv[1]*mr)/ms+0.7071067811865475*m0sm0r_inv[0]*m1sm1r[1]+0.7071067811865475*m1sm1r[0]*m0sm0r_inv[1]-0.7071067811865475*m0rm2s[0]*m0sm0r_inv[1]-0.7071067811865475*m0sm0r_inv[0]*m0rm2s[1]; 
  binop_mul_1d_ser_p1(nu_sr_inv, cM2sr_temp, cM2sr_temp); 
  binop_mul_1d_ser_p1(alphaE, cM2sr_temp, cM2sr); 
  m2sr[0] = m2s[0] + cM2sr[0]; 
  m2sr[1] = m2s[1] + cM2sr[1]; 
} 
 