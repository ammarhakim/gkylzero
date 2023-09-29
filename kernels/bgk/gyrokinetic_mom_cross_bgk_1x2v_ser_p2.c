#include <gkyl_mom_cross_bgk_gyrokinetic_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
 
GKYL_CU_DH void gyrokinetic_mom_cross_bgk_1x2v_ser_p2(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross) 
{ 
  // m_:             mass. 
  // moms_:          moments of the distribution function. 
  // moms_cross:     cross moments. 
 
  const double ms = m_self; 
  const double mr = m_other; 
  const double *m0s = &moms_self[0]; 
  const double *m1s = &moms_self[3]; 
  const double *m2s = &moms_self[6]; 
  const double *m0r = &moms_other[0]; 
  const double *m1r = &moms_other[3]; 
  const double *m2r = &moms_other[6]; 
 
  double *m0sr = &moms_cross[0]; 
  double *m0rs = &moms_cross[3]; 
  double *m1sr = &moms_cross[6]; 
  double *m1rs = &moms_cross[9]; 
  double *m2sr = &moms_cross[12]; 
  double *m2rs = &moms_cross[15]; 
 
  double msM0sNusr[3] = {0.0}; 
  double mrM0rNurs[3] = {0.0}; 
  double mM0Nu[3] = {0.0}; 
  double mM0Nu_inv[3] = {0.0}; 
  double alphaE[3] = {0.0}; 

  double m0sm0r[3] = {0.0}; 
  double m0sm0r_inv[3] = {0.0}; 
  double m0sm1r[3] = {0.0}; 
  double m0rm1s[3] = {0.0}; 
  double cM1sr_temp[3] = {0.0}; 
  double cM1rs_temp[3] = {0.0}; 
  double cM1sr[3] = {0.0}; 
  double cM1rs[3] = {0.0}; 

  double m0sm2r[3] = {0.0}; 
  double m0rm2s[3] = {0.0}; 
  double m1sm1r[3] = {0.0}; 
  double cM2sr_temp[3] = {0.0}; 
  double cM2rs_temp[3] = {0.0}; 
  double cM2sr[3] = {0.0}; 
  double cM2rs[3] = {0.0}; 

  msM0sNusr[0] = ms * m0s[0] * nu_sr; 
  mrM0rNurs[0] = mr * m0r[0] * nu_rs; 
  mM0Nu[0] = msM0sNusr[0] + mrM0rNurs[0]; 
  msM0sNusr[1] = ms * m0s[1] * nu_sr; 
  mrM0rNurs[1] = mr * m0r[1] * nu_rs; 
  mM0Nu[1] = msM0sNusr[1] + mrM0rNurs[1]; 
  msM0sNusr[2] = ms * m0s[2] * nu_sr; 
  mrM0rNurs[2] = mr * m0r[2] * nu_rs; 
  mM0Nu[2] = msM0sNusr[2] + mrM0rNurs[2]; 
  ser_1x_p2_inv(mM0Nu, mM0Nu_inv); 
  binop_mul_1d_ser_p2(msM0sNusr, mrM0rNurs, alphaE); 
  binop_mul_1d_ser_p2(alphaE, mM0Nu_inv, alphaE); 
  alphaE[0] = alphaE[0] * 2 * (1+beta) / (ms+mr); 
  alphaE[1] = alphaE[1] * 2 * (1+beta) / (ms+mr); 
  alphaE[2] = alphaE[2] * 2 * (1+beta) / (ms+mr); 

  m0sr[0] = m0s[0]; 
  m0rs[0] = m0r[0]; 
  m0sr[1] = m0s[1]; 
  m0rs[1] = m0r[1]; 
  m0sr[2] = m0s[2]; 
  m0rs[2] = m0r[2]; 
 
  binop_mul_1d_ser_p2(m0s, m0r, m0sm0r); 
  ser_1x_p2_inv(m0sm0r, m0sm0r_inv); 
  binop_mul_1d_ser_p2(m0s, m1r, m0sm1r); 
  binop_mul_1d_ser_p2(m0r, m1s, m0rm1s); 
  cM1sr_temp[0] = (0.3535533905932737*m0sm0r_inv[2]*m0sm1r[2]*mr)/(ms*nu_sr)-(0.3535533905932737*m0rm1s[2]*m0sm0r_inv[2]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm0r_inv[1]*m0sm1r[1]*mr)/(ms*nu_sr)-(0.3535533905932737*m0rm1s[1]*m0sm0r_inv[1]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[0]*mr)/(ms*nu_sr)-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[0]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm0r_inv[2]*m0sm1r[2])/nu_sr-(0.3535533905932737*m0rm1s[2]*m0sm0r_inv[2])/nu_sr+(0.3535533905932737*m0sm0r_inv[1]*m0sm1r[1])/nu_sr-(0.3535533905932737*m0rm1s[1]*m0sm0r_inv[1])/nu_sr+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[0])/nu_sr-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[0])/nu_sr; 
  cM1sr_temp[1] = (0.3162277660168379*m0sm0r_inv[1]*m0sm1r[2]*mr)/(ms*nu_sr)+(0.3162277660168379*m0sm1r[1]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.3162277660168379*m0sm0r_inv[1]*m0rm1s[2]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[1]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[1]*mr)/(ms*nu_sr)-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[1]*mr)/(ms*nu_sr)-(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[1]*mr)/(ms*nu_sr)+(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[2])/nu_sr+(0.3162277660168379*m0sm1r[1]*m0sm0r_inv[2])/nu_sr-(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[2])/nu_sr-(0.3162277660168379*m0sm0r_inv[1]*m0rm1s[2])/nu_sr+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[1])/nu_sr+(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[1])/nu_sr-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[1])/nu_sr-(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[1])/nu_sr; 
  cM1sr_temp[2] = (0.2258769757263128*m0sm0r_inv[2]*m0sm1r[2]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[2]*mr)/(ms*nu_sr)-(0.2258769757263128*m0rm1s[2]*m0sm0r_inv[2]*mr)/(ms*nu_sr)+(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[2]*mr)/(ms*nu_sr)+(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[1]*mr)/(ms*nu_sr)-(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[1]*mr)/(ms*nu_sr)+(0.2258769757263128*m0sm0r_inv[2]*m0sm1r[2])/nu_sr+(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[2])/nu_sr-(0.2258769757263128*m0rm1s[2]*m0sm0r_inv[2])/nu_sr+(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[2])/nu_sr-(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[2])/nu_sr-(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[2])/nu_sr+(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[1])/nu_sr-(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[1])/nu_sr; 
  cM1rs_temp[0] = (-(0.3535533905932737*m0sm0r_inv[2]*m0sm1r[2]*ms)/(mr*nu_rs))+(0.3535533905932737*m0rm1s[2]*m0sm0r_inv[2]*ms)/(mr*nu_rs)-(0.3535533905932737*m0sm0r_inv[1]*m0sm1r[1]*ms)/(mr*nu_rs)+(0.3535533905932737*m0rm1s[1]*m0sm0r_inv[1]*ms)/(mr*nu_rs)-(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[0]*ms)/(mr*nu_rs)+(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[0]*ms)/(mr*nu_rs)-(0.3535533905932737*m0sm0r_inv[2]*m0sm1r[2])/nu_rs+(0.3535533905932737*m0rm1s[2]*m0sm0r_inv[2])/nu_rs-(0.3535533905932737*m0sm0r_inv[1]*m0sm1r[1])/nu_rs+(0.3535533905932737*m0rm1s[1]*m0sm0r_inv[1])/nu_rs-(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[0])/nu_rs+(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[0])/nu_rs; 
  cM1rs_temp[1] = (-(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[2]*ms)/(mr*nu_rs))-(0.3162277660168379*m0sm1r[1]*m0sm0r_inv[2]*ms)/(mr*nu_rs)+(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[2]*ms)/(mr*nu_rs)+(0.3162277660168379*m0sm0r_inv[1]*m0rm1s[2]*ms)/(mr*nu_rs)-(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[1]*ms)/(mr*nu_rs)-(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[1]*ms)/(mr*nu_rs)+(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[1]*ms)/(mr*nu_rs)+(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[1]*ms)/(mr*nu_rs)-(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[2])/nu_rs-(0.3162277660168379*m0sm1r[1]*m0sm0r_inv[2])/nu_rs+(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[2])/nu_rs+(0.3162277660168379*m0sm0r_inv[1]*m0rm1s[2])/nu_rs-(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[1])/nu_rs-(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[1])/nu_rs+(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[1])/nu_rs+(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[1])/nu_rs; 
  cM1rs_temp[2] = (-(0.2258769757263128*m0sm0r_inv[2]*m0sm1r[2]*ms)/(mr*nu_rs))-(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[2]*ms)/(mr*nu_rs)+(0.2258769757263128*m0rm1s[2]*m0sm0r_inv[2]*ms)/(mr*nu_rs)-(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[2]*ms)/(mr*nu_rs)+(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[2]*ms)/(mr*nu_rs)+(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[2]*ms)/(mr*nu_rs)-(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[1]*ms)/(mr*nu_rs)+(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[1]*ms)/(mr*nu_rs)-(0.2258769757263128*m0sm0r_inv[2]*m0sm1r[2])/nu_rs-(0.3535533905932737*m0sm0r_inv[0]*m0sm1r[2])/nu_rs+(0.2258769757263128*m0rm1s[2]*m0sm0r_inv[2])/nu_rs-(0.3535533905932737*m0sm1r[0]*m0sm0r_inv[2])/nu_rs+(0.3535533905932737*m0rm1s[0]*m0sm0r_inv[2])/nu_rs+(0.3535533905932737*m0sm0r_inv[0]*m0rm1s[2])/nu_rs-(0.3162277660168379*m0sm0r_inv[1]*m0sm1r[1])/nu_rs+(0.3162277660168379*m0rm1s[1]*m0sm0r_inv[1])/nu_rs; 
  binop_mul_1d_ser_p2(alphaE, cM1sr_temp, cM1sr); 
  binop_mul_1d_ser_p2(alphaE, cM1rs_temp, cM1rs); 
  m1sr[0] = m1s[0] + cM1sr[0]; 
  m1rs[0] = m1r[0] + cM1rs[0]; 
  m1sr[1] = m1s[1] + cM1sr[1]; 
  m1rs[1] = m1r[1] + cM1rs[1]; 
  m1sr[2] = m1s[2] + cM1sr[2]; 
  m1rs[2] = m1r[2] + cM1rs[2]; 
 
  binop_mul_1d_ser_p2(m0s, m2r, m0sm2r); 
  binop_mul_1d_ser_p2(m0r, m2s, m0rm2s); 
  binop_mul_1d_ser_p2(m1s, m1r, m1sm1r); 
  cM2sr_temp[0] = (-(0.7071067811865475*m0sm0r_inv[2]*m1sm1r[2]*mr)/(ms*nu_sr))+(0.7071067811865475*m0sm0r_inv[2]*m0sm2r[2]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[1]*m1sm1r[1]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[1]*m0sm2r[1]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[0]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[0]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[2]*m1sm1r[2])/nu_sr-(0.7071067811865475*m0rm2s[2]*m0sm0r_inv[2])/nu_sr+(0.7071067811865475*m0sm0r_inv[1]*m1sm1r[1])/nu_sr-(0.7071067811865475*m0rm2s[1]*m0sm0r_inv[1])/nu_sr+(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[0])/nu_sr-(0.7071067811865475*m0rm2s[0]*m0sm0r_inv[0])/nu_sr; 
  cM2sr_temp[1] = (-(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[2]*mr)/(ms*nu_sr))+(0.6324555320336759*m0sm0r_inv[1]*m0sm2r[2]*mr)/(ms*nu_sr)-(0.6324555320336759*m1sm1r[1]*m0sm0r_inv[2]*mr)/(ms*nu_sr)+(0.6324555320336759*m0sm2r[1]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[1]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[1]*mr)/(ms*nu_sr)-(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[1]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm2r[0]*m0sm0r_inv[1]*mr)/(ms*nu_sr)+(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[2])/nu_sr+(0.6324555320336759*m1sm1r[1]*m0sm0r_inv[2])/nu_sr-(0.6324555320336759*m0rm2s[1]*m0sm0r_inv[2])/nu_sr-(0.6324555320336759*m0sm0r_inv[1]*m0rm2s[2])/nu_sr+(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[1])/nu_sr+(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[1])/nu_sr-(0.7071067811865475*m0rm2s[0]*m0sm0r_inv[1])/nu_sr-(0.7071067811865475*m0sm0r_inv[0]*m0rm2s[1])/nu_sr; 
  cM2sr_temp[2] = (-(0.4517539514526256*m0sm0r_inv[2]*m1sm1r[2]*mr)/(ms*nu_sr))-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[2]*mr)/(ms*nu_sr)+(0.4517539514526256*m0sm0r_inv[2]*m0sm2r[2]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[2]*mr)/(ms*nu_sr)-(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[2]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm2r[0]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[1]*mr)/(ms*nu_sr)+(0.6324555320336759*m0sm0r_inv[1]*m0sm2r[1]*mr)/(ms*nu_sr)+(0.4517539514526256*m0sm0r_inv[2]*m1sm1r[2])/nu_sr+(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[2])/nu_sr-(0.4517539514526256*m0rm2s[2]*m0sm0r_inv[2])/nu_sr+(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[2])/nu_sr-(0.7071067811865475*m0rm2s[0]*m0sm0r_inv[2])/nu_sr-(0.7071067811865475*m0sm0r_inv[0]*m0rm2s[2])/nu_sr+(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[1])/nu_sr-(0.6324555320336759*m0rm2s[1]*m0sm0r_inv[1])/nu_sr; 
  cM2rs_temp[0] = (0.7071067811865475*m0sm0r_inv[2]*m1sm1r[2]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[2]*m0sm2r[2]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[1]*m1sm1r[1]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[1]*m0sm2r[1]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[0]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[0]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[2]*m1sm1r[2])/nu_sr+(0.7071067811865475*m0rm2s[2]*m0sm0r_inv[2])/nu_sr-(0.7071067811865475*m0sm0r_inv[1]*m1sm1r[1])/nu_sr+(0.7071067811865475*m0rm2s[1]*m0sm0r_inv[1])/nu_sr-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[0])/nu_sr+(0.7071067811865475*m0rm2s[0]*m0sm0r_inv[0])/nu_sr; 
  cM2rs_temp[1] = (0.6324555320336759*m0sm0r_inv[1]*m1sm1r[2]*mr)/(ms*nu_sr)-(0.6324555320336759*m0sm0r_inv[1]*m0sm2r[2]*mr)/(ms*nu_sr)+(0.6324555320336759*m1sm1r[1]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.6324555320336759*m0sm2r[1]*m0sm0r_inv[2]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[1]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[1]*mr)/(ms*nu_sr)+(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[1]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm2r[0]*m0sm0r_inv[1]*mr)/(ms*nu_sr)-(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[2])/nu_sr-(0.6324555320336759*m1sm1r[1]*m0sm0r_inv[2])/nu_sr+(0.6324555320336759*m0rm2s[1]*m0sm0r_inv[2])/nu_sr+(0.6324555320336759*m0sm0r_inv[1]*m0rm2s[2])/nu_sr-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[1])/nu_sr-(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[1])/nu_sr+(0.7071067811865475*m0rm2s[0]*m0sm0r_inv[1])/nu_sr+(0.7071067811865475*m0sm0r_inv[0]*m0rm2s[1])/nu_sr; 
  cM2rs_temp[2] = (0.4517539514526256*m0sm0r_inv[2]*m1sm1r[2]*mr)/(ms*nu_sr)+(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[2]*mr)/(ms*nu_sr)-(0.4517539514526256*m0sm0r_inv[2]*m0sm2r[2]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm0r_inv[0]*m0sm2r[2]*mr)/(ms*nu_sr)+(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[2]*mr)/(ms*nu_sr)-(0.7071067811865475*m0sm2r[0]*m0sm0r_inv[2]*mr)/(ms*nu_sr)+(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[1]*mr)/(ms*nu_sr)-(0.6324555320336759*m0sm0r_inv[1]*m0sm2r[1]*mr)/(ms*nu_sr)-(0.4517539514526256*m0sm0r_inv[2]*m1sm1r[2])/nu_sr-(0.7071067811865475*m0sm0r_inv[0]*m1sm1r[2])/nu_sr+(0.4517539514526256*m0rm2s[2]*m0sm0r_inv[2])/nu_sr-(0.7071067811865475*m1sm1r[0]*m0sm0r_inv[2])/nu_sr+(0.7071067811865475*m0rm2s[0]*m0sm0r_inv[2])/nu_sr+(0.7071067811865475*m0sm0r_inv[0]*m0rm2s[2])/nu_sr-(0.6324555320336759*m0sm0r_inv[1]*m1sm1r[1])/nu_sr+(0.6324555320336759*m0rm2s[1]*m0sm0r_inv[1])/nu_sr; 
  binop_mul_1d_ser_p2(alphaE, cM2sr_temp, cM2sr); 
  binop_mul_1d_ser_p2(alphaE, cM2rs_temp, cM2rs); 
  m2sr[0] = m2s[0] + cM2sr[0]; 
  m2rs[0] = m2r[0] + cM2rs[0]; 
  m2sr[1] = m2s[1] + cM2sr[1]; 
  m2rs[1] = m2r[1] + cM2rs[1]; 
  m2sr[2] = m2s[2] + cM2sr[2]; 
  m2rs[2] = m2r[2] + cM2rs[2]; 
} 
 
