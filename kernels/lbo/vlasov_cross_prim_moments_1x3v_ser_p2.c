#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greeneFac, const double m_self, const double *u_self, const double *vtsq_self, const double m_other, const double *u_other, const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections) 
{ 
  // greeneFac:            free parameter beta+1 multiplied by other factors. 
  // nu, m:                collisionality and mass. 
  // moms:                 moments of the distribution function. 
  // u,vtSq:               self primitive moments: mean flow velocity and thermal speed squared. 
  // boundary_corrections: corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:     cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*(2.23606797749979*moms[2]-1.732050807568877*moms[1]+moms[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*moms[2]+1.732050807568877*moms[1]+moms[0]) < 0) cellAvg = true; 
 
  double m0r[3] = {0.0}; 
  double m1r[9] = {0.0}; 
  double m2r[3] = {0.0}; 
  double cMr[9] = {0.0}; 
  double cEr[3] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = moms[3]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    m1r[3] = moms[6]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    m1r[6] = moms[9]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = 0.0; 
    cMr[8] = 0.0; 
    m2r[0] = moms[12]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    cEr[0] = boundary_corrections[9]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m1r[0] = moms[3]; 
    m1r[1] = moms[4]; 
    m1r[2] = moms[5]; 
    m1r[3] = moms[6]; 
    m1r[4] = moms[7]; 
    m1r[5] = moms[8]; 
    m1r[6] = moms[9]; 
    m1r[7] = moms[10]; 
    m1r[8] = moms[11]; 
    m2r[0] = moms[12]; 
    m2r[1] = moms[13]; 
    m2r[2] = moms[14]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = boundary_corrections[7]; 
    cMr[8] = boundary_corrections[8]; 
    cEr[0] = boundary_corrections[9]; 
    cEr[1] = boundary_corrections[10]; 
    cEr[2] = boundary_corrections[11]; 
  } 
 
  double mnu = 0.7071067811865475*nu[0]*m_self; 
  double momRHS[9] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,0,1,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,0,2,1.414213562373095*m0r[2]*mnu); 
  gkyl_mat_set(A,1,0,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,1,1,1.264911064067352*m0r[2]*mnu+1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,1,2,1.264911064067352*m0r[1]*mnu); 
  gkyl_mat_set(A,2,0,1.414213562373095*m0r[2]*mnu); 
  gkyl_mat_set(A,2,1,1.264911064067352*m0r[1]*mnu); 
  gkyl_mat_set(A,2,2,0.9035079029052515*m0r[2]*mnu+1.414213562373095*m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,9,-1.414213562373095*cMr[0]*mnu); 
  gkyl_mat_set(A,0,10,-1.414213562373095*cMr[1]*mnu); 
  gkyl_mat_set(A,0,11,-1.414213562373095*cMr[2]*mnu); 
  gkyl_mat_set(A,1,9,-1.414213562373095*cMr[1]*mnu); 
  gkyl_mat_set(A,1,10,(-1.264911064067352*cMr[2]*mnu)-1.414213562373095*cMr[0]*mnu); 
  gkyl_mat_set(A,1,11,-1.264911064067352*cMr[1]*mnu); 
  gkyl_mat_set(A,2,9,-1.414213562373095*cMr[2]*mnu); 
  gkyl_mat_set(A,2,10,-1.264911064067352*cMr[1]*mnu); 
  gkyl_mat_set(A,2,11,(-0.9035079029052515*cMr[2]*mnu)-1.414213562373095*cMr[0]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,9,0,(-0.5*m0r[2]*u_self[2]*mnu)-0.5*m0r[2]*u_other[2]*mnu-0.5*m0r[1]*u_self[1]*mnu-0.5*m0r[1]*u_other[1]*mnu-0.5*m0r[0]*u_self[0]*mnu-0.5*m0r[0]*u_other[0]*mnu+1.414213562373095*m1r[0]*mnu); 
  gkyl_mat_set(A,9,1,(-0.4472135954999579*m0r[1]*u_self[2]*mnu)-0.4472135954999579*m0r[1]*u_other[2]*mnu-0.4472135954999579*u_self[1]*m0r[2]*mnu-0.4472135954999579*u_other[1]*m0r[2]*mnu-0.5*m0r[0]*u_self[1]*mnu-0.5*m0r[0]*u_other[1]*mnu+1.414213562373095*m1r[1]*mnu-0.5*u_self[0]*m0r[1]*mnu-0.5*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,9,2,(-0.31943828249997*m0r[2]*u_self[2]*mnu)-0.5*m0r[0]*u_self[2]*mnu-0.31943828249997*m0r[2]*u_other[2]*mnu-0.5*m0r[0]*u_other[2]*mnu+1.414213562373095*m1r[2]*mnu-0.5*u_self[0]*m0r[2]*mnu-0.5*u_other[0]*m0r[2]*mnu-0.4472135954999579*m0r[1]*u_self[1]*mnu-0.4472135954999579*m0r[1]*u_other[1]*mnu); 
  gkyl_mat_set(A,10,0,(-0.4472135954999579*m0r[1]*u_self[2]*mnu)-0.4472135954999579*m0r[1]*u_other[2]*mnu-0.4472135954999579*u_self[1]*m0r[2]*mnu-0.4472135954999579*u_other[1]*m0r[2]*mnu-0.5*m0r[0]*u_self[1]*mnu-0.5*m0r[0]*u_other[1]*mnu+1.414213562373095*m1r[1]*mnu-0.5*u_self[0]*m0r[1]*mnu-0.5*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,10,1,(-0.7857142857142857*m0r[2]*u_self[2]*mnu)-0.4472135954999579*m0r[0]*u_self[2]*mnu-0.7857142857142857*m0r[2]*u_other[2]*mnu-0.4472135954999579*m0r[0]*u_other[2]*mnu+1.264911064067352*m1r[2]*mnu-0.4472135954999579*u_self[0]*m0r[2]*mnu-0.4472135954999579*u_other[0]*m0r[2]*mnu-0.9*m0r[1]*u_self[1]*mnu-0.9*m0r[1]*u_other[1]*mnu-0.5*m0r[0]*u_self[0]*mnu-0.5*m0r[0]*u_other[0]*mnu+1.414213562373095*m1r[0]*mnu); 
  gkyl_mat_set(A,10,2,(-0.7857142857142857*m0r[1]*u_self[2]*mnu)-0.7857142857142857*m0r[1]*u_other[2]*mnu-0.7857142857142857*u_self[1]*m0r[2]*mnu-0.7857142857142857*u_other[1]*m0r[2]*mnu-0.4472135954999579*m0r[0]*u_self[1]*mnu-0.4472135954999579*m0r[0]*u_other[1]*mnu+1.264911064067352*m1r[1]*mnu-0.4472135954999579*u_self[0]*m0r[1]*mnu-0.4472135954999579*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,11,0,(-0.31943828249997*m0r[2]*u_self[2]*mnu)-0.5*m0r[0]*u_self[2]*mnu-0.31943828249997*m0r[2]*u_other[2]*mnu-0.5*m0r[0]*u_other[2]*mnu+1.414213562373095*m1r[2]*mnu-0.5*u_self[0]*m0r[2]*mnu-0.5*u_other[0]*m0r[2]*mnu-0.4472135954999579*m0r[1]*u_self[1]*mnu-0.4472135954999579*m0r[1]*u_other[1]*mnu); 
  gkyl_mat_set(A,11,1,(-0.7857142857142857*m0r[1]*u_self[2]*mnu)-0.7857142857142857*m0r[1]*u_other[2]*mnu-0.7857142857142857*u_self[1]*m0r[2]*mnu-0.7857142857142857*u_other[1]*m0r[2]*mnu-0.4472135954999579*m0r[0]*u_self[1]*mnu-0.4472135954999579*m0r[0]*u_other[1]*mnu+1.264911064067352*m1r[1]*mnu-0.4472135954999579*u_self[0]*m0r[1]*mnu-0.4472135954999579*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,11,2,(-1.071428571428571*m0r[2]*u_self[2]*mnu)-0.31943828249997*m0r[0]*u_self[2]*mnu-1.071428571428571*m0r[2]*u_other[2]*mnu-0.31943828249997*m0r[0]*u_other[2]*mnu+0.9035079029052515*m1r[2]*mnu-0.31943828249997*u_self[0]*m0r[2]*mnu-0.31943828249997*u_other[0]*m0r[2]*mnu-0.7857142857142857*m0r[1]*u_self[1]*mnu-0.7857142857142857*m0r[1]*u_other[1]*mnu-0.5*m0r[0]*u_self[0]*mnu-0.5*m0r[0]*u_other[0]*mnu+1.414213562373095*m1r[0]*mnu); 
 
  momRHS[0] += -0.01428571428571429*(((22.3606797749979*greeneFac[2]+35.0*greeneFac[0])*m0r[2]+35.0*m0r[0]*greeneFac[2]+31.30495168499706*greeneFac[1]*m0r[1])*u_self[2]+(((-22.3606797749979*greeneFac[2])-35.0*greeneFac[0])*m0r[2]-35.0*m0r[0]*greeneFac[2]-31.30495168499706*greeneFac[1]*m0r[1])*u_other[2]+((35.0*u_self[0]-35.0*u_other[0])*greeneFac[2]+31.30495168499706*greeneFac[1]*u_self[1]-31.30495168499706*greeneFac[1]*u_other[1])*m0r[2]+(31.30495168499706*m0r[1]*u_self[1]-31.30495168499706*m0r[1]*u_other[1])*greeneFac[2]+(35.0*greeneFac[0]*m0r[1]+35.0*m0r[0]*greeneFac[1])*u_self[1]+((-35.0*greeneFac[0]*m0r[1])-35.0*m0r[0]*greeneFac[1])*u_other[1]+(35.0*u_self[0]-35.0*u_other[0])*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0]*u_self[0]-35.0*greeneFac[0]*m0r[0]*u_other[0]-140.0*m1r[0])*mnu; 
  momRHS[1] += -0.01428571428571429*((55.0*greeneFac[1]*m0r[2]+55.0*m0r[1]*greeneFac[2]+31.30495168499706*greeneFac[0]*m0r[1]+31.30495168499706*m0r[0]*greeneFac[1])*u_self[2]+((-55.0*greeneFac[1]*m0r[2])-55.0*m0r[1]*greeneFac[2]-31.30495168499706*greeneFac[0]*m0r[1]-31.30495168499706*m0r[0]*greeneFac[1])*u_other[2]+((55.0*u_self[1]-55.0*u_other[1])*greeneFac[2]+31.30495168499706*greeneFac[0]*u_self[1]-31.30495168499706*greeneFac[0]*u_other[1]+(31.30495168499706*u_self[0]-31.30495168499706*u_other[0])*greeneFac[1])*m0r[2]+(31.30495168499706*m0r[0]*u_self[1]-31.30495168499706*m0r[0]*u_other[1]+(31.30495168499706*u_self[0]-31.30495168499706*u_other[0])*m0r[1])*greeneFac[2]+(63.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[1]+((-63.0*greeneFac[1]*m0r[1])-35.0*greeneFac[0]*m0r[0])*u_other[1]-140.0*m1r[1]+(35.0*greeneFac[0]*u_self[0]-35.0*greeneFac[0]*u_other[0])*m0r[1]+(35.0*m0r[0]*u_self[0]-35.0*m0r[0]*u_other[0])*greeneFac[1])*mnu; 
  momRHS[2] += -0.01428571428571429*(((75.0*greeneFac[2]+22.3606797749979*greeneFac[0])*m0r[2]+22.3606797749979*m0r[0]*greeneFac[2]+55.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[2]+(((-75.0*greeneFac[2])-22.3606797749979*greeneFac[0])*m0r[2]-22.3606797749979*m0r[0]*greeneFac[2]-55.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[2]-140.0*m1r[2]+((22.3606797749979*u_self[0]-22.3606797749979*u_other[0])*greeneFac[2]+55.0*greeneFac[1]*u_self[1]-55.0*greeneFac[1]*u_other[1]+35.0*greeneFac[0]*u_self[0]-35.0*greeneFac[0]*u_other[0])*m0r[2]+(55.0*m0r[1]*u_self[1]-55.0*m0r[1]*u_other[1]+35.0*m0r[0]*u_self[0]-35.0*m0r[0]*u_other[0])*greeneFac[2]+(31.30495168499706*greeneFac[0]*m0r[1]+31.30495168499706*m0r[0]*greeneFac[1])*u_self[1]+((-31.30495168499706*greeneFac[0]*m0r[1])-31.30495168499706*m0r[0]*greeneFac[1])*u_other[1]+(31.30495168499706*u_self[0]-31.30495168499706*u_other[0])*greeneFac[1]*m0r[1])*mnu; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossY ... // 
  gkyl_mat_set(A,3,3,1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,3,4,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,3,5,1.414213562373095*m0r[2]*mnu); 
  gkyl_mat_set(A,4,3,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,4,4,1.264911064067352*m0r[2]*mnu+1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,4,5,1.264911064067352*m0r[1]*mnu); 
  gkyl_mat_set(A,5,3,1.414213562373095*m0r[2]*mnu); 
  gkyl_mat_set(A,5,4,1.264911064067352*m0r[1]*mnu); 
  gkyl_mat_set(A,5,5,0.9035079029052515*m0r[2]*mnu+1.414213562373095*m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,3,9,-1.414213562373095*cMr[3]*mnu); 
  gkyl_mat_set(A,3,10,-1.414213562373095*cMr[4]*mnu); 
  gkyl_mat_set(A,3,11,-1.414213562373095*cMr[5]*mnu); 
  gkyl_mat_set(A,4,9,-1.414213562373095*cMr[4]*mnu); 
  gkyl_mat_set(A,4,10,(-1.264911064067352*cMr[5]*mnu)-1.414213562373095*cMr[3]*mnu); 
  gkyl_mat_set(A,4,11,-1.264911064067352*cMr[4]*mnu); 
  gkyl_mat_set(A,5,9,-1.414213562373095*cMr[5]*mnu); 
  gkyl_mat_set(A,5,10,-1.264911064067352*cMr[4]*mnu); 
  gkyl_mat_set(A,5,11,(-0.9035079029052515*cMr[5]*mnu)-1.414213562373095*cMr[3]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1Y and uCrossY ... // 
  gkyl_mat_set(A,9,3,(-0.5*m0r[2]*u_self[5]*mnu)-0.5*m0r[2]*u_other[5]*mnu-0.5*m0r[1]*u_self[4]*mnu-0.5*m0r[1]*u_other[4]*mnu-0.5*m0r[0]*u_self[3]*mnu-0.5*m0r[0]*u_other[3]*mnu+1.414213562373095*m1r[3]*mnu); 
  gkyl_mat_set(A,9,4,(-0.4472135954999579*m0r[1]*u_self[5]*mnu)-0.4472135954999579*m0r[1]*u_other[5]*mnu-0.4472135954999579*m0r[2]*u_self[4]*mnu-0.5*m0r[0]*u_self[4]*mnu-0.4472135954999579*m0r[2]*u_other[4]*mnu-0.5*m0r[0]*u_other[4]*mnu+1.414213562373095*m1r[4]*mnu-0.5*m0r[1]*u_self[3]*mnu-0.5*m0r[1]*u_other[3]*mnu); 
  gkyl_mat_set(A,9,5,(-0.31943828249997*m0r[2]*u_self[5]*mnu)-0.5*m0r[0]*u_self[5]*mnu-0.31943828249997*m0r[2]*u_other[5]*mnu-0.5*m0r[0]*u_other[5]*mnu+1.414213562373095*m1r[5]*mnu-0.4472135954999579*m0r[1]*u_self[4]*mnu-0.4472135954999579*m0r[1]*u_other[4]*mnu-0.5*m0r[2]*u_self[3]*mnu-0.5*m0r[2]*u_other[3]*mnu); 
  gkyl_mat_set(A,10,3,(-0.4472135954999579*m0r[1]*u_self[5]*mnu)-0.4472135954999579*m0r[1]*u_other[5]*mnu-0.4472135954999579*m0r[2]*u_self[4]*mnu-0.5*m0r[0]*u_self[4]*mnu-0.4472135954999579*m0r[2]*u_other[4]*mnu-0.5*m0r[0]*u_other[4]*mnu+1.414213562373095*m1r[4]*mnu-0.5*m0r[1]*u_self[3]*mnu-0.5*m0r[1]*u_other[3]*mnu); 
  gkyl_mat_set(A,10,4,(-0.7857142857142857*m0r[2]*u_self[5]*mnu)-0.4472135954999579*m0r[0]*u_self[5]*mnu-0.7857142857142857*m0r[2]*u_other[5]*mnu-0.4472135954999579*m0r[0]*u_other[5]*mnu+1.264911064067352*m1r[5]*mnu-0.9*m0r[1]*u_self[4]*mnu-0.9*m0r[1]*u_other[4]*mnu-0.4472135954999579*m0r[2]*u_self[3]*mnu-0.5*m0r[0]*u_self[3]*mnu-0.4472135954999579*m0r[2]*u_other[3]*mnu-0.5*m0r[0]*u_other[3]*mnu+1.414213562373095*m1r[3]*mnu); 
  gkyl_mat_set(A,10,5,(-0.7857142857142857*m0r[1]*u_self[5]*mnu)-0.7857142857142857*m0r[1]*u_other[5]*mnu-0.7857142857142857*m0r[2]*u_self[4]*mnu-0.4472135954999579*m0r[0]*u_self[4]*mnu-0.7857142857142857*m0r[2]*u_other[4]*mnu-0.4472135954999579*m0r[0]*u_other[4]*mnu+1.264911064067352*m1r[4]*mnu-0.4472135954999579*m0r[1]*u_self[3]*mnu-0.4472135954999579*m0r[1]*u_other[3]*mnu); 
  gkyl_mat_set(A,11,3,(-0.31943828249997*m0r[2]*u_self[5]*mnu)-0.5*m0r[0]*u_self[5]*mnu-0.31943828249997*m0r[2]*u_other[5]*mnu-0.5*m0r[0]*u_other[5]*mnu+1.414213562373095*m1r[5]*mnu-0.4472135954999579*m0r[1]*u_self[4]*mnu-0.4472135954999579*m0r[1]*u_other[4]*mnu-0.5*m0r[2]*u_self[3]*mnu-0.5*m0r[2]*u_other[3]*mnu); 
  gkyl_mat_set(A,11,4,(-0.7857142857142857*m0r[1]*u_self[5]*mnu)-0.7857142857142857*m0r[1]*u_other[5]*mnu-0.7857142857142857*m0r[2]*u_self[4]*mnu-0.4472135954999579*m0r[0]*u_self[4]*mnu-0.7857142857142857*m0r[2]*u_other[4]*mnu-0.4472135954999579*m0r[0]*u_other[4]*mnu+1.264911064067352*m1r[4]*mnu-0.4472135954999579*m0r[1]*u_self[3]*mnu-0.4472135954999579*m0r[1]*u_other[3]*mnu); 
  gkyl_mat_set(A,11,5,(-1.071428571428571*m0r[2]*u_self[5]*mnu)-0.31943828249997*m0r[0]*u_self[5]*mnu-1.071428571428571*m0r[2]*u_other[5]*mnu-0.31943828249997*m0r[0]*u_other[5]*mnu+0.9035079029052515*m1r[5]*mnu-0.7857142857142857*m0r[1]*u_self[4]*mnu-0.7857142857142857*m0r[1]*u_other[4]*mnu-0.31943828249997*m0r[2]*u_self[3]*mnu-0.5*m0r[0]*u_self[3]*mnu-0.31943828249997*m0r[2]*u_other[3]*mnu-0.5*m0r[0]*u_other[3]*mnu+1.414213562373095*m1r[3]*mnu); 
 
  momRHS[3] += -0.01428571428571429*(((22.3606797749979*greeneFac[2]+35.0*greeneFac[0])*m0r[2]+35.0*m0r[0]*greeneFac[2]+31.30495168499706*greeneFac[1]*m0r[1])*u_self[5]+(((-22.3606797749979*greeneFac[2])-35.0*greeneFac[0])*m0r[2]-35.0*m0r[0]*greeneFac[2]-31.30495168499706*greeneFac[1]*m0r[1])*u_other[5]+(31.30495168499706*greeneFac[1]*m0r[2]+31.30495168499706*m0r[1]*greeneFac[2]+35.0*greeneFac[0]*m0r[1]+35.0*m0r[0]*greeneFac[1])*u_self[4]+((-31.30495168499706*greeneFac[1]*m0r[2])-31.30495168499706*m0r[1]*greeneFac[2]-35.0*greeneFac[0]*m0r[1]-35.0*m0r[0]*greeneFac[1])*u_other[4]+(35.0*greeneFac[2]*m0r[2]+35.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[3]+((-35.0*greeneFac[2]*m0r[2])-35.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[3]-140.0*m1r[3])*mnu; 
  momRHS[4] += -0.01428571428571429*((55.0*greeneFac[1]*m0r[2]+55.0*m0r[1]*greeneFac[2]+31.30495168499706*greeneFac[0]*m0r[1]+31.30495168499706*m0r[0]*greeneFac[1])*u_self[5]+((-55.0*greeneFac[1]*m0r[2])-55.0*m0r[1]*greeneFac[2]-31.30495168499706*greeneFac[0]*m0r[1]-31.30495168499706*m0r[0]*greeneFac[1])*u_other[5]+((55.0*greeneFac[2]+31.30495168499706*greeneFac[0])*m0r[2]+31.30495168499706*m0r[0]*greeneFac[2]+63.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[4]+(((-55.0*greeneFac[2])-31.30495168499706*greeneFac[0])*m0r[2]-31.30495168499706*m0r[0]*greeneFac[2]-63.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[4]-140.0*m1r[4]+(31.30495168499706*greeneFac[1]*m0r[2]+31.30495168499706*m0r[1]*greeneFac[2]+35.0*greeneFac[0]*m0r[1]+35.0*m0r[0]*greeneFac[1])*u_self[3]+((-31.30495168499706*greeneFac[1]*m0r[2])-31.30495168499706*m0r[1]*greeneFac[2]-35.0*greeneFac[0]*m0r[1]-35.0*m0r[0]*greeneFac[1])*u_other[3])*mnu; 
  momRHS[5] += -0.01428571428571429*(((75.0*greeneFac[2]+22.3606797749979*greeneFac[0])*m0r[2]+22.3606797749979*m0r[0]*greeneFac[2]+55.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[5]+(((-75.0*greeneFac[2])-22.3606797749979*greeneFac[0])*m0r[2]-22.3606797749979*m0r[0]*greeneFac[2]-55.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[5]-140.0*m1r[5]+(55.0*greeneFac[1]*m0r[2]+55.0*m0r[1]*greeneFac[2]+31.30495168499706*greeneFac[0]*m0r[1]+31.30495168499706*m0r[0]*greeneFac[1])*u_self[4]+((-55.0*greeneFac[1]*m0r[2])-55.0*m0r[1]*greeneFac[2]-31.30495168499706*greeneFac[0]*m0r[1]-31.30495168499706*m0r[0]*greeneFac[1])*u_other[4]+((22.3606797749979*greeneFac[2]+35.0*greeneFac[0])*m0r[2]+35.0*m0r[0]*greeneFac[2]+31.30495168499706*greeneFac[1]*m0r[1])*u_self[3]+(((-22.3606797749979*greeneFac[2])-35.0*greeneFac[0])*m0r[2]-35.0*m0r[0]*greeneFac[2]-31.30495168499706*greeneFac[1]*m0r[1])*u_other[3])*mnu; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossZ ... // 
  gkyl_mat_set(A,6,6,1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,6,7,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,6,8,1.414213562373095*m0r[2]*mnu); 
  gkyl_mat_set(A,7,6,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,7,7,1.264911064067352*m0r[2]*mnu+1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,7,8,1.264911064067352*m0r[1]*mnu); 
  gkyl_mat_set(A,8,6,1.414213562373095*m0r[2]*mnu); 
  gkyl_mat_set(A,8,7,1.264911064067352*m0r[1]*mnu); 
  gkyl_mat_set(A,8,8,0.9035079029052515*m0r[2]*mnu+1.414213562373095*m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,6,9,-1.414213562373095*cMr[6]*mnu); 
  gkyl_mat_set(A,6,10,-1.414213562373095*cMr[7]*mnu); 
  gkyl_mat_set(A,6,11,-1.414213562373095*cMr[8]*mnu); 
  gkyl_mat_set(A,7,9,-1.414213562373095*cMr[7]*mnu); 
  gkyl_mat_set(A,7,10,(-1.264911064067352*cMr[8]*mnu)-1.414213562373095*cMr[6]*mnu); 
  gkyl_mat_set(A,7,11,-1.264911064067352*cMr[7]*mnu); 
  gkyl_mat_set(A,8,9,-1.414213562373095*cMr[8]*mnu); 
  gkyl_mat_set(A,8,10,-1.264911064067352*cMr[7]*mnu); 
  gkyl_mat_set(A,8,11,(-0.9035079029052515*cMr[8]*mnu)-1.414213562373095*cMr[6]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1Z and uCrossZ ... // 
  gkyl_mat_set(A,9,6,(-0.5*m0r[2]*u_self[8]*mnu)-0.5*m0r[2]*u_other[8]*mnu-0.5*m0r[1]*u_self[7]*mnu-0.5*m0r[1]*u_other[7]*mnu-0.5*m0r[0]*u_self[6]*mnu-0.5*m0r[0]*u_other[6]*mnu+1.414213562373095*m1r[6]*mnu); 
  gkyl_mat_set(A,9,7,(-0.4472135954999579*m0r[1]*u_self[8]*mnu)-0.4472135954999579*m0r[1]*u_other[8]*mnu-0.4472135954999579*m0r[2]*u_self[7]*mnu-0.5*m0r[0]*u_self[7]*mnu-0.4472135954999579*m0r[2]*u_other[7]*mnu-0.5*m0r[0]*u_other[7]*mnu+1.414213562373095*m1r[7]*mnu-0.5*m0r[1]*u_self[6]*mnu-0.5*m0r[1]*u_other[6]*mnu); 
  gkyl_mat_set(A,9,8,(-0.31943828249997*m0r[2]*u_self[8]*mnu)-0.5*m0r[0]*u_self[8]*mnu-0.31943828249997*m0r[2]*u_other[8]*mnu-0.5*m0r[0]*u_other[8]*mnu+1.414213562373095*m1r[8]*mnu-0.4472135954999579*m0r[1]*u_self[7]*mnu-0.4472135954999579*m0r[1]*u_other[7]*mnu-0.5*m0r[2]*u_self[6]*mnu-0.5*m0r[2]*u_other[6]*mnu); 
  gkyl_mat_set(A,10,6,(-0.4472135954999579*m0r[1]*u_self[8]*mnu)-0.4472135954999579*m0r[1]*u_other[8]*mnu-0.4472135954999579*m0r[2]*u_self[7]*mnu-0.5*m0r[0]*u_self[7]*mnu-0.4472135954999579*m0r[2]*u_other[7]*mnu-0.5*m0r[0]*u_other[7]*mnu+1.414213562373095*m1r[7]*mnu-0.5*m0r[1]*u_self[6]*mnu-0.5*m0r[1]*u_other[6]*mnu); 
  gkyl_mat_set(A,10,7,(-0.7857142857142857*m0r[2]*u_self[8]*mnu)-0.4472135954999579*m0r[0]*u_self[8]*mnu-0.7857142857142857*m0r[2]*u_other[8]*mnu-0.4472135954999579*m0r[0]*u_other[8]*mnu+1.264911064067352*m1r[8]*mnu-0.9*m0r[1]*u_self[7]*mnu-0.9*m0r[1]*u_other[7]*mnu-0.4472135954999579*m0r[2]*u_self[6]*mnu-0.5*m0r[0]*u_self[6]*mnu-0.4472135954999579*m0r[2]*u_other[6]*mnu-0.5*m0r[0]*u_other[6]*mnu+1.414213562373095*m1r[6]*mnu); 
  gkyl_mat_set(A,10,8,(-0.7857142857142857*m0r[1]*u_self[8]*mnu)-0.7857142857142857*m0r[1]*u_other[8]*mnu-0.7857142857142857*m0r[2]*u_self[7]*mnu-0.4472135954999579*m0r[0]*u_self[7]*mnu-0.7857142857142857*m0r[2]*u_other[7]*mnu-0.4472135954999579*m0r[0]*u_other[7]*mnu+1.264911064067352*m1r[7]*mnu-0.4472135954999579*m0r[1]*u_self[6]*mnu-0.4472135954999579*m0r[1]*u_other[6]*mnu); 
  gkyl_mat_set(A,11,6,(-0.31943828249997*m0r[2]*u_self[8]*mnu)-0.5*m0r[0]*u_self[8]*mnu-0.31943828249997*m0r[2]*u_other[8]*mnu-0.5*m0r[0]*u_other[8]*mnu+1.414213562373095*m1r[8]*mnu-0.4472135954999579*m0r[1]*u_self[7]*mnu-0.4472135954999579*m0r[1]*u_other[7]*mnu-0.5*m0r[2]*u_self[6]*mnu-0.5*m0r[2]*u_other[6]*mnu); 
  gkyl_mat_set(A,11,7,(-0.7857142857142857*m0r[1]*u_self[8]*mnu)-0.7857142857142857*m0r[1]*u_other[8]*mnu-0.7857142857142857*m0r[2]*u_self[7]*mnu-0.4472135954999579*m0r[0]*u_self[7]*mnu-0.7857142857142857*m0r[2]*u_other[7]*mnu-0.4472135954999579*m0r[0]*u_other[7]*mnu+1.264911064067352*m1r[7]*mnu-0.4472135954999579*m0r[1]*u_self[6]*mnu-0.4472135954999579*m0r[1]*u_other[6]*mnu); 
  gkyl_mat_set(A,11,8,(-1.071428571428571*m0r[2]*u_self[8]*mnu)-0.31943828249997*m0r[0]*u_self[8]*mnu-1.071428571428571*m0r[2]*u_other[8]*mnu-0.31943828249997*m0r[0]*u_other[8]*mnu+0.9035079029052515*m1r[8]*mnu-0.7857142857142857*m0r[1]*u_self[7]*mnu-0.7857142857142857*m0r[1]*u_other[7]*mnu-0.31943828249997*m0r[2]*u_self[6]*mnu-0.5*m0r[0]*u_self[6]*mnu-0.31943828249997*m0r[2]*u_other[6]*mnu-0.5*m0r[0]*u_other[6]*mnu+1.414213562373095*m1r[6]*mnu); 
 
  momRHS[6] += -0.01428571428571429*(((22.3606797749979*greeneFac[2]+35.0*greeneFac[0])*m0r[2]+35.0*m0r[0]*greeneFac[2]+31.30495168499706*greeneFac[1]*m0r[1])*u_self[8]+(((-22.3606797749979*greeneFac[2])-35.0*greeneFac[0])*m0r[2]-35.0*m0r[0]*greeneFac[2]-31.30495168499706*greeneFac[1]*m0r[1])*u_other[8]+(31.30495168499706*greeneFac[1]*m0r[2]+31.30495168499706*m0r[1]*greeneFac[2]+35.0*greeneFac[0]*m0r[1]+35.0*m0r[0]*greeneFac[1])*u_self[7]+((-31.30495168499706*greeneFac[1]*m0r[2])-31.30495168499706*m0r[1]*greeneFac[2]-35.0*greeneFac[0]*m0r[1]-35.0*m0r[0]*greeneFac[1])*u_other[7]+(35.0*greeneFac[2]*m0r[2]+35.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[6]+((-35.0*greeneFac[2]*m0r[2])-35.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[6]-140.0*m1r[6])*mnu; 
  momRHS[7] += -0.01428571428571429*((55.0*greeneFac[1]*m0r[2]+55.0*m0r[1]*greeneFac[2]+31.30495168499706*greeneFac[0]*m0r[1]+31.30495168499706*m0r[0]*greeneFac[1])*u_self[8]+((-55.0*greeneFac[1]*m0r[2])-55.0*m0r[1]*greeneFac[2]-31.30495168499706*greeneFac[0]*m0r[1]-31.30495168499706*m0r[0]*greeneFac[1])*u_other[8]+((55.0*greeneFac[2]+31.30495168499706*greeneFac[0])*m0r[2]+31.30495168499706*m0r[0]*greeneFac[2]+63.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[7]+(((-55.0*greeneFac[2])-31.30495168499706*greeneFac[0])*m0r[2]-31.30495168499706*m0r[0]*greeneFac[2]-63.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[7]-140.0*m1r[7]+(31.30495168499706*greeneFac[1]*m0r[2]+31.30495168499706*m0r[1]*greeneFac[2]+35.0*greeneFac[0]*m0r[1]+35.0*m0r[0]*greeneFac[1])*u_self[6]+((-31.30495168499706*greeneFac[1]*m0r[2])-31.30495168499706*m0r[1]*greeneFac[2]-35.0*greeneFac[0]*m0r[1]-35.0*m0r[0]*greeneFac[1])*u_other[6])*mnu; 
  momRHS[8] += -0.01428571428571429*(((75.0*greeneFac[2]+22.3606797749979*greeneFac[0])*m0r[2]+22.3606797749979*m0r[0]*greeneFac[2]+55.0*greeneFac[1]*m0r[1]+35.0*greeneFac[0]*m0r[0])*u_self[8]+(((-75.0*greeneFac[2])-22.3606797749979*greeneFac[0])*m0r[2]-22.3606797749979*m0r[0]*greeneFac[2]-55.0*greeneFac[1]*m0r[1]-35.0*greeneFac[0]*m0r[0])*u_other[8]-140.0*m1r[8]+(55.0*greeneFac[1]*m0r[2]+55.0*m0r[1]*greeneFac[2]+31.30495168499706*greeneFac[0]*m0r[1]+31.30495168499706*m0r[0]*greeneFac[1])*u_self[7]+((-55.0*greeneFac[1]*m0r[2])-55.0*m0r[1]*greeneFac[2]-31.30495168499706*greeneFac[0]*m0r[1]-31.30495168499706*m0r[0]*greeneFac[1])*u_other[7]+((22.3606797749979*greeneFac[2]+35.0*greeneFac[0])*m0r[2]+35.0*m0r[0]*greeneFac[2]+31.30495168499706*greeneFac[1]*m0r[1])*u_self[6]+(((-22.3606797749979*greeneFac[2])-35.0*greeneFac[0])*m0r[2]-35.0*m0r[0]*greeneFac[2]-31.30495168499706*greeneFac[1]*m0r[1])*u_other[6])*mnu; 
 
  double ucMSelf[3] = {0.0}; 
  double ucMOther[3] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+2]*u_self[a0+2]+0.7071067811865475*cMr[a0+1]*u_self[a0+1]+0.7071067811865475*cMr[a0]*u_self[a0]; 
    ucMSelf[1] += 0.6324555320336759*cMr[a0+1]*u_self[a0+2]+0.6324555320336759*u_self[a0+1]*cMr[a0+2]+0.7071067811865475*cMr[a0]*u_self[a0+1]+0.7071067811865475*u_self[a0]*cMr[a0+1]; 
    ucMSelf[2] += 0.4517539514526256*cMr[a0+2]*u_self[a0+2]+0.7071067811865475*cMr[a0]*u_self[a0+2]+0.7071067811865475*u_self[a0]*cMr[a0+2]+0.6324555320336759*cMr[a0+1]*u_self[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+2]*u_other[a0+2]+0.7071067811865475*cMr[a0+1]*u_other[a0+1]+0.7071067811865475*cMr[a0]*u_other[a0]; 
    ucMOther[1] += 0.6324555320336759*cMr[a0+1]*u_other[a0+2]+0.6324555320336759*u_other[a0+1]*cMr[a0+2]+0.7071067811865475*cMr[a0]*u_other[a0+1]+0.7071067811865475*u_other[a0]*cMr[a0+1]; 
    ucMOther[2] += 0.4517539514526256*cMr[a0+2]*u_other[a0+2]+0.7071067811865475*cMr[a0]*u_other[a0+2]+0.7071067811865475*u_other[a0]*cMr[a0+2]+0.6324555320336759*cMr[a0+1]*u_other[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,21,9,0.7071067811865475*ucMSelf[0]*mnu+0.7071067811865475*ucMOther[0]*mnu+4.242640687119286*m0r[0]*mnu-1.414213562373095*cEr[0]*mnu); 
  gkyl_mat_set(A,21,10,0.7071067811865475*ucMSelf[1]*mnu+0.7071067811865475*ucMOther[1]*mnu+4.242640687119286*m0r[1]*mnu-1.414213562373095*cEr[1]*mnu); 
  gkyl_mat_set(A,21,11,0.7071067811865475*ucMSelf[2]*mnu+0.7071067811865475*ucMOther[2]*mnu+4.242640687119286*m0r[2]*mnu-1.414213562373095*cEr[2]*mnu); 
  gkyl_mat_set(A,22,9,0.7071067811865475*ucMSelf[1]*mnu+0.7071067811865475*ucMOther[1]*mnu+4.242640687119286*m0r[1]*mnu-1.414213562373095*cEr[1]*mnu); 
  gkyl_mat_set(A,22,10,0.6324555320336759*ucMSelf[2]*mnu+0.6324555320336759*ucMOther[2]*mnu+3.794733192202056*m0r[2]*mnu-1.264911064067352*cEr[2]*mnu+0.7071067811865475*ucMSelf[0]*mnu+0.7071067811865475*ucMOther[0]*mnu+4.242640687119286*m0r[0]*mnu-1.414213562373095*cEr[0]*mnu); 
  gkyl_mat_set(A,22,11,0.6324555320336759*ucMSelf[1]*mnu+0.6324555320336759*ucMOther[1]*mnu+3.794733192202056*m0r[1]*mnu-1.264911064067352*cEr[1]*mnu); 
  gkyl_mat_set(A,23,9,0.7071067811865475*ucMSelf[2]*mnu+0.7071067811865475*ucMOther[2]*mnu+4.242640687119286*m0r[2]*mnu-1.414213562373095*cEr[2]*mnu); 
  gkyl_mat_set(A,23,10,0.6324555320336759*ucMSelf[1]*mnu+0.6324555320336759*ucMOther[1]*mnu+3.794733192202056*m0r[1]*mnu-1.264911064067352*cEr[1]*mnu); 
  gkyl_mat_set(A,23,11,0.4517539514526256*ucMSelf[2]*mnu+0.4517539514526256*ucMOther[2]*mnu+2.710523708715754*m0r[2]*mnu-0.9035079029052515*cEr[2]*mnu+0.7071067811865475*ucMSelf[0]*mnu+0.7071067811865475*ucMOther[0]*mnu+4.242640687119286*m0r[0]*mnu-1.414213562373095*cEr[0]*mnu); 
 
  double uM1Self[3] = {0.0}; 
  double uM1Other[3] = {0.0}; 
  double uSumSq[3] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 3*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.7071067811865475*m1r[8]*u_self[a0+2]+0.7071067811865475*m1r[7]*u_self[a0+1]+0.7071067811865475*m1r[6]*u_self[a0]; 
    uM1Self[1] += 0.6324555320336759*m1r[7]*u_self[a0+2]+0.6324555320336759*m1r[8]*u_self[a0+1]+0.7071067811865475*m1r[6]*u_self[a0+1]+0.7071067811865475*m1r[7]*u_self[a0]; 
    uM1Self[2] += 0.4517539514526256*m1r[8]*u_self[a0+2]+0.7071067811865475*m1r[6]*u_self[a0+2]+0.6324555320336759*m1r[7]*u_self[a0+1]+0.7071067811865475*m1r[8]*u_self[a0]; 
    uM1Other[0] += 0.7071067811865475*m1r[8]*u_other[a0+2]+0.7071067811865475*m1r[7]*u_other[a0+1]+0.7071067811865475*m1r[6]*u_other[a0]; 
    uM1Other[1] += 0.6324555320336759*m1r[7]*u_other[a0+2]+0.6324555320336759*m1r[8]*u_other[a0+1]+0.7071067811865475*m1r[6]*u_other[a0+1]+0.7071067811865475*m1r[7]*u_other[a0]; 
    uM1Other[2] += 0.4517539514526256*m1r[8]*u_other[a0+2]+0.7071067811865475*m1r[6]*u_other[a0+2]+0.6324555320336759*m1r[7]*u_other[a0+1]+0.7071067811865475*m1r[8]*u_other[a0]; 
  const double u_self0R2 = pow(u_self[a0],2);
  const double u_self1R2 = pow(u_self[a0+1],2);
  const double u_self2R2 = pow(u_self[a0+2],2);
  const double u_other0R2 = pow(u_other[a0],2);
  const double u_other1R2 = pow(u_other[a0+1],2);
  const double u_other2R2 = pow(u_other[a0+2],2);

  uSumSq[0] += 0.7071067811865475*u_self2R2-1.414213562373095*u_other[a0+2]*u_self[a0+2]+0.7071067811865475*u_other2R2+0.7071067811865475*u_self1R2-1.414213562373095*u_other[a0+1]*u_self[a0+1]+0.7071067811865475*u_other1R2+0.7071067811865475*u_self0R2-1.414213562373095*u_other[a0]*u_self[a0]+0.7071067811865475*u_other0R2; 
  uSumSq[1] += 1.264911064067352*u_self[a0+1]*u_self[a0+2]-1.264911064067352*u_other[a0+1]*u_self[a0+2]-1.264911064067352*u_self[a0+1]*u_other[a0+2]+1.264911064067352*u_other[a0+1]*u_other[a0+2]+1.414213562373095*u_self[a0]*u_self[a0+1]-1.414213562373095*u_other[a0]*u_self[a0+1]-1.414213562373095*u_self[a0]*u_other[a0+1]+1.414213562373095*u_other[a0]*u_other[a0+1]; 
  uSumSq[2] += 0.4517539514526256*u_self2R2-0.9035079029052515*u_other[a0+2]*u_self[a0+2]+1.414213562373095*u_self[a0]*u_self[a0+2]-1.414213562373095*u_other[a0]*u_self[a0+2]+0.4517539514526256*u_other2R2-1.414213562373095*u_self[a0]*u_other[a0+2]+1.414213562373095*u_other[a0]*u_other[a0+2]+0.6324555320336759*u_self1R2-1.264911064067352*u_other[a0+1]*u_self[a0+1]+0.6324555320336759*u_other1R2; 
  } 
 
  double enRHS[3] = {0.0}; 
  enRHS[0] = (-(134.1640786499874*greeneFac[2]*m0r[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other))-(210.0*greeneFac[0]*m0r[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*m0r[0]*greeneFac[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*greeneFac[1]*m0r[1]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(22.3606797749979*greeneFac[2]*m0r[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*m0r[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*m0r[0]*greeneFac[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*greeneFac[1]*m0r[1]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*vtsq_self[0]*greeneFac[2]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*uSumSq[0]*greeneFac[2]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*greeneFac[1]*vtsq_self[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*greeneFac[1]*uSumSq[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*m0r[1]*vtsq_self[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*m0r[1]*uSumSq[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*greeneFac[0]*m0r[1]*vtsq_self[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*m0r[0]*greeneFac[1]*vtsq_self[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*m0r[1]*uSumSq[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*m0r[0]*greeneFac[1]*uSumSq[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*vtsq_self[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*uSumSq[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*greeneFac[0]*m0r[0]*vtsq_self[0]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*m0r[0]*uSumSq[0]*m_self*mnu)/(140.0*m_self+140.0*m_other)+(134.1640786499874*greeneFac[2]*m0r[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*m0r[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*m0r[0]*greeneFac[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*greeneFac[1]*m0r[1]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(22.3606797749979*greeneFac[2]*m0r[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*m0r[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*m0r[0]*greeneFac[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*greeneFac[1]*m0r[1]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*vtsq_other[0]*greeneFac[2]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*uSumSq[0]*greeneFac[2]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*greeneFac[1]*vtsq_other[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*greeneFac[1]*uSumSq[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*m0r[1]*vtsq_other[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*m0r[1]*uSumSq[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*m0r[1]*vtsq_other[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*m0r[0]*greeneFac[1]*vtsq_other[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*m0r[1]*uSumSq[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*m0r[0]*greeneFac[1]*uSumSq[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*vtsq_other[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*uSumSq[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*m0r[0]*vtsq_other[0]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*m0r[0]*uSumSq[0]*m_other*mnu)/(140.0*m_self+140.0*m_other)-1.0*uM1Self[0]*mnu-1.0*uM1Other[0]*mnu+2.0*m2r[0]*mnu; 
  enRHS[1] = (-(330.0*greeneFac[1]*m0r[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other))-(330.0*m0r[1]*greeneFac[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*greeneFac[0]*m0r[1]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*m0r[0]*greeneFac[1]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(55.0*greeneFac[1]*m0r[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(55.0*m0r[1]*greeneFac[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*greeneFac[0]*m0r[1]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*m0r[0]*greeneFac[1]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(330.0*vtsq_self[1]*greeneFac[2]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(55.0*uSumSq[1]*greeneFac[2]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*greeneFac[0]*vtsq_self[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*greeneFac[0]*uSumSq[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*vtsq_self[0]*greeneFac[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*uSumSq[0]*greeneFac[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*m0r[0]*vtsq_self[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*m0r[0]*uSumSq[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*vtsq_self[0]*m0r[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*uSumSq[0]*m0r[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(378.0*greeneFac[1]*m0r[1]*vtsq_self[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*greeneFac[0]*m0r[0]*vtsq_self[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(63.0*greeneFac[1]*m0r[1]*uSumSq[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*m0r[0]*uSumSq[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*greeneFac[0]*vtsq_self[0]*m0r[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*uSumSq[0]*m0r[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*m0r[0]*vtsq_self[0]*greeneFac[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*m0r[0]*uSumSq[0]*greeneFac[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)+(330.0*greeneFac[1]*m0r[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(330.0*m0r[1]*greeneFac[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*greeneFac[0]*m0r[1]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*m0r[0]*greeneFac[1]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(55.0*greeneFac[1]*m0r[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(55.0*m0r[1]*greeneFac[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*greeneFac[0]*m0r[1]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*m0r[0]*greeneFac[1]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(330.0*vtsq_other[1]*greeneFac[2]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(55.0*uSumSq[1]*greeneFac[2]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*greeneFac[0]*vtsq_other[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*greeneFac[0]*uSumSq[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*vtsq_other[0]*greeneFac[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*uSumSq[0]*greeneFac[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*m0r[0]*vtsq_other[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*m0r[0]*uSumSq[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*vtsq_other[0]*m0r[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*uSumSq[0]*m0r[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(378.0*greeneFac[1]*m0r[1]*vtsq_other[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*m0r[0]*vtsq_other[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(63.0*greeneFac[1]*m0r[1]*uSumSq[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*m0r[0]*uSumSq[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*vtsq_other[0]*m0r[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*uSumSq[0]*m0r[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*m0r[0]*vtsq_other[0]*greeneFac[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*m0r[0]*uSumSq[0]*greeneFac[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)-1.0*uM1Self[1]*mnu-1.0*uM1Other[1]*mnu+2.0*m2r[1]*mnu; 
  enRHS[2] = (-(450.0*greeneFac[2]*m0r[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other))-(134.1640786499874*greeneFac[0]*m0r[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(134.1640786499874*m0r[0]*greeneFac[2]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(330.0*greeneFac[1]*m0r[1]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*greeneFac[0]*m0r[0]*vtsq_self[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(75.0*greeneFac[2]*m0r[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(22.3606797749979*greeneFac[0]*m0r[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(22.3606797749979*m0r[0]*greeneFac[2]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(55.0*greeneFac[1]*m0r[1]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*m0r[0]*uSumSq[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(134.1640786499874*vtsq_self[0]*greeneFac[2]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(22.3606797749979*uSumSq[0]*greeneFac[2]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(330.0*greeneFac[1]*vtsq_self[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(55.0*greeneFac[1]*uSumSq[1]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*greeneFac[0]*vtsq_self[0]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*greeneFac[0]*uSumSq[0]*m0r[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(330.0*m0r[1]*vtsq_self[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(55.0*m0r[1]*uSumSq[1]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(210.0*m0r[0]*vtsq_self[0]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(35.0*m0r[0]*uSumSq[0]*greeneFac[2]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*greeneFac[0]*m0r[1]*vtsq_self[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*m0r[0]*greeneFac[1]*vtsq_self[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*greeneFac[0]*m0r[1]*uSumSq[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*m0r[0]*greeneFac[1]*uSumSq[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(187.8297101099823*vtsq_self[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)-(31.30495168499706*uSumSq[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(140.0*m_self+140.0*m_other)+(450.0*greeneFac[2]*m0r[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(134.1640786499874*greeneFac[0]*m0r[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(134.1640786499874*m0r[0]*greeneFac[2]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(330.0*greeneFac[1]*m0r[1]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*m0r[0]*vtsq_other[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(75.0*greeneFac[2]*m0r[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(22.3606797749979*greeneFac[0]*m0r[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(22.3606797749979*m0r[0]*greeneFac[2]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(55.0*greeneFac[1]*m0r[1]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*m0r[0]*uSumSq[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(134.1640786499874*vtsq_other[0]*greeneFac[2]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(22.3606797749979*uSumSq[0]*greeneFac[2]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(330.0*greeneFac[1]*vtsq_other[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(55.0*greeneFac[1]*uSumSq[1]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*greeneFac[0]*vtsq_other[0]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*greeneFac[0]*uSumSq[0]*m0r[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(330.0*m0r[1]*vtsq_other[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(55.0*m0r[1]*uSumSq[1]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(210.0*m0r[0]*vtsq_other[0]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(35.0*m0r[0]*uSumSq[0]*greeneFac[2]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*greeneFac[0]*m0r[1]*vtsq_other[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*m0r[0]*greeneFac[1]*vtsq_other[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*greeneFac[0]*m0r[1]*uSumSq[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*m0r[0]*greeneFac[1]*uSumSq[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(187.8297101099823*vtsq_other[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)+(31.30495168499706*uSumSq[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(140.0*m_self+140.0*m_other)-1.0*uM1Self[2]*mnu-1.0*uM1Other[2]*mnu+2.0*m2r[2]*mnu; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,momRHS[3]); 
  gkyl_mat_set(rhs,4,0,momRHS[4]); 
  gkyl_mat_set(rhs,5,0,momRHS[5]); 
  gkyl_mat_set(rhs,6,0,momRHS[6]); 
  gkyl_mat_set(rhs,7,0,momRHS[7]); 
  gkyl_mat_set(rhs,8,0,momRHS[8]); 
  gkyl_mat_set(rhs,9,0,enRHS[0]); 
  gkyl_mat_set(rhs,10,0,enRHS[1]); 
  gkyl_mat_set(rhs,11,0,enRHS[2]); 
} 
 
