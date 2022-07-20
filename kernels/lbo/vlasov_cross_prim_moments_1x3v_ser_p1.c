#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greeneFac, const double m_self, const double *u_self, const double *vtsq_self, const double m_other, const double *u_other, const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections) 
{ 
  // greeneFac:            free parameter beta+1 multiplied by other factors. 
  // nu, m:                collisionality and mass. 
  // moms:                 moments of the distribution function. 
  // u,vtSq:               self primitive moments: mean flow velocity and thermal speed squared. 
  // boundary_corrections: corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:     cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(2.449489742783178*moms[1]-1.414213562373095*moms[0]) < 0) cellAvg = true; 
  if (0.5*(2.449489742783178*moms[1]+1.414213562373095*moms[0]) < 0) cellAvg = true; 
 
  double m0r[2] = {0.0}; 
  double m1r[6] = {0.0}; 
  double m2r[2] = {0.0}; 
  double cMr[6] = {0.0}; 
  double cEr[2] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m1r[0] = moms[2]; 
    m1r[1] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    m1r[2] = moms[4]; 
    m1r[3] = 0.0; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = 0.0; 
    m1r[4] = moms[6]; 
    m1r[5] = 0.0; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = 0.0; 
    m2r[0] = moms[8]; 
    m2r[1] = 0.0; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = 0.0; 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m1r[0] = moms[2]; 
    m1r[1] = moms[3]; 
    m1r[2] = moms[4]; 
    m1r[3] = moms[5]; 
    m1r[4] = moms[6]; 
    m1r[5] = moms[7]; 
    m2r[0] = moms[8]; 
    m2r[1] = moms[9]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = boundary_corrections[7]; 
  } 
 
  double mnu = 0.7071067811865475*nu[0]*m_self; 
  double momRHS[6] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,0,1,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,1,0,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,1,1,1.414213562373095*m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,6,-1.414213562373095*cMr[0]*mnu); 
  gkyl_mat_set(A,0,7,-1.414213562373095*cMr[1]*mnu); 
  gkyl_mat_set(A,1,6,-1.414213562373095*cMr[1]*mnu); 
  gkyl_mat_set(A,1,7,-1.414213562373095*cMr[0]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,6,0,(-0.5*m0r[1]*u_self[1]*mnu)-0.5*m0r[1]*u_other[1]*mnu-0.5*m0r[0]*u_self[0]*mnu-0.5*m0r[0]*u_other[0]*mnu+1.414213562373095*m1r[0]*mnu); 
  gkyl_mat_set(A,6,1,(-0.5*m0r[0]*u_self[1]*mnu)-0.5*m0r[0]*u_other[1]*mnu+1.414213562373095*m1r[1]*mnu-0.5*u_self[0]*m0r[1]*mnu-0.5*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,7,0,(-0.5*m0r[0]*u_self[1]*mnu)-0.5*m0r[0]*u_other[1]*mnu+1.414213562373095*m1r[1]*mnu-0.5*u_self[0]*m0r[1]*mnu-0.5*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,7,1,(-0.9*m0r[1]*u_self[1]*mnu)-0.9*m0r[1]*u_other[1]*mnu-0.5*m0r[0]*u_self[0]*mnu-0.5*m0r[0]*u_other[0]*mnu+1.414213562373095*m1r[0]*mnu); 
 
  momRHS[0] += -0.5*((greeneFac[0]*m0r[1]+m0r[0]*greeneFac[1])*u_self[1]+((-1.0*greeneFac[0]*m0r[1])-1.0*m0r[0]*greeneFac[1])*u_other[1]+(u_self[0]-1.0*u_other[0])*greeneFac[1]*m0r[1]+greeneFac[0]*m0r[0]*u_self[0]-1.0*greeneFac[0]*m0r[0]*u_other[0]-4.0*m1r[0])*mnu; 
  momRHS[1] += -0.1*((9.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[1]+((-9.0*greeneFac[1]*m0r[1])-5.0*greeneFac[0]*m0r[0])*u_other[1]-20.0*m1r[1]+(5.0*greeneFac[0]*u_self[0]-5.0*greeneFac[0]*u_other[0])*m0r[1]+(5.0*m0r[0]*u_self[0]-5.0*m0r[0]*u_other[0])*greeneFac[1])*mnu; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossY ... // 
  gkyl_mat_set(A,2,2,1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,2,3,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,3,2,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,3,3,1.414213562373095*m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,2,6,-1.414213562373095*cMr[2]*mnu); 
  gkyl_mat_set(A,2,7,-1.414213562373095*cMr[3]*mnu); 
  gkyl_mat_set(A,3,6,-1.414213562373095*cMr[3]*mnu); 
  gkyl_mat_set(A,3,7,-1.414213562373095*cMr[2]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1Y and uCrossY ... // 
  gkyl_mat_set(A,6,2,(-0.5*m0r[1]*u_self[3]*mnu)-0.5*m0r[1]*u_other[3]*mnu-0.5*m0r[0]*u_self[2]*mnu-0.5*m0r[0]*u_other[2]*mnu+1.414213562373095*m1r[2]*mnu); 
  gkyl_mat_set(A,6,3,(-0.5*m0r[0]*u_self[3]*mnu)-0.5*m0r[0]*u_other[3]*mnu+1.414213562373095*m1r[3]*mnu-0.5*m0r[1]*u_self[2]*mnu-0.5*m0r[1]*u_other[2]*mnu); 
  gkyl_mat_set(A,7,2,(-0.5*m0r[0]*u_self[3]*mnu)-0.5*m0r[0]*u_other[3]*mnu+1.414213562373095*m1r[3]*mnu-0.5*m0r[1]*u_self[2]*mnu-0.5*m0r[1]*u_other[2]*mnu); 
  gkyl_mat_set(A,7,3,(-0.9*m0r[1]*u_self[3]*mnu)-0.9*m0r[1]*u_other[3]*mnu-0.5*m0r[0]*u_self[2]*mnu-0.5*m0r[0]*u_other[2]*mnu+1.414213562373095*m1r[2]*mnu); 
 
  momRHS[2] += -0.5*((greeneFac[0]*m0r[1]+m0r[0]*greeneFac[1])*u_self[3]+((-1.0*greeneFac[0]*m0r[1])-1.0*m0r[0]*greeneFac[1])*u_other[3]+(greeneFac[1]*m0r[1]+greeneFac[0]*m0r[0])*u_self[2]+((-1.0*greeneFac[1]*m0r[1])-1.0*greeneFac[0]*m0r[0])*u_other[2]-4.0*m1r[2])*mnu; 
  momRHS[3] += -0.1*((9.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[3]+((-9.0*greeneFac[1]*m0r[1])-5.0*greeneFac[0]*m0r[0])*u_other[3]-20.0*m1r[3]+(5.0*greeneFac[0]*m0r[1]+5.0*m0r[0]*greeneFac[1])*u_self[2]+((-5.0*greeneFac[0]*m0r[1])-5.0*m0r[0]*greeneFac[1])*u_other[2])*mnu; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossZ ... // 
  gkyl_mat_set(A,4,4,1.414213562373095*m0r[0]*mnu); 
  gkyl_mat_set(A,4,5,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,5,4,1.414213562373095*m0r[1]*mnu); 
  gkyl_mat_set(A,5,5,1.414213562373095*m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,4,6,-1.414213562373095*cMr[4]*mnu); 
  gkyl_mat_set(A,4,7,-1.414213562373095*cMr[5]*mnu); 
  gkyl_mat_set(A,5,6,-1.414213562373095*cMr[5]*mnu); 
  gkyl_mat_set(A,5,7,-1.414213562373095*cMr[4]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1Z and uCrossZ ... // 
  gkyl_mat_set(A,6,4,(-0.5*m0r[1]*u_self[5]*mnu)-0.5*m0r[1]*u_other[5]*mnu-0.5*m0r[0]*u_self[4]*mnu-0.5*m0r[0]*u_other[4]*mnu+1.414213562373095*m1r[4]*mnu); 
  gkyl_mat_set(A,6,5,(-0.5*m0r[0]*u_self[5]*mnu)-0.5*m0r[0]*u_other[5]*mnu+1.414213562373095*m1r[5]*mnu-0.5*m0r[1]*u_self[4]*mnu-0.5*m0r[1]*u_other[4]*mnu); 
  gkyl_mat_set(A,7,4,(-0.5*m0r[0]*u_self[5]*mnu)-0.5*m0r[0]*u_other[5]*mnu+1.414213562373095*m1r[5]*mnu-0.5*m0r[1]*u_self[4]*mnu-0.5*m0r[1]*u_other[4]*mnu); 
  gkyl_mat_set(A,7,5,(-0.9*m0r[1]*u_self[5]*mnu)-0.9*m0r[1]*u_other[5]*mnu-0.5*m0r[0]*u_self[4]*mnu-0.5*m0r[0]*u_other[4]*mnu+1.414213562373095*m1r[4]*mnu); 
 
  momRHS[4] += -0.5*((greeneFac[0]*m0r[1]+m0r[0]*greeneFac[1])*u_self[5]+((-1.0*greeneFac[0]*m0r[1])-1.0*m0r[0]*greeneFac[1])*u_other[5]+(greeneFac[1]*m0r[1]+greeneFac[0]*m0r[0])*u_self[4]+((-1.0*greeneFac[1]*m0r[1])-1.0*greeneFac[0]*m0r[0])*u_other[4]-4.0*m1r[4])*mnu; 
  momRHS[5] += -0.1*((9.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[5]+((-9.0*greeneFac[1]*m0r[1])-5.0*greeneFac[0]*m0r[0])*u_other[5]-20.0*m1r[5]+(5.0*greeneFac[0]*m0r[1]+5.0*m0r[0]*greeneFac[1])*u_self[4]+((-5.0*greeneFac[0]*m0r[1])-5.0*m0r[0]*greeneFac[1])*u_other[4])*mnu; 
 
  double ucMSelf[2] = {0.0}; 
  double ucMOther[2] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+1]*u_self[a0+1]+0.7071067811865475*cMr[a0]*u_self[a0]; 
    ucMSelf[1] += 0.7071067811865475*cMr[a0]*u_self[a0+1]+0.7071067811865475*u_self[a0]*cMr[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+1]*u_other[a0+1]+0.7071067811865475*cMr[a0]*u_other[a0]; 
    ucMOther[1] += 0.7071067811865475*cMr[a0]*u_other[a0+1]+0.7071067811865475*u_other[a0]*cMr[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,14,6,0.7071067811865475*ucMSelf[0]*mnu+0.7071067811865475*ucMOther[0]*mnu+4.242640687119286*m0r[0]*mnu-1.414213562373095*cEr[0]*mnu); 
  gkyl_mat_set(A,14,7,0.7071067811865475*ucMSelf[1]*mnu+0.7071067811865475*ucMOther[1]*mnu+4.242640687119286*m0r[1]*mnu-1.414213562373095*cEr[1]*mnu); 
  gkyl_mat_set(A,15,6,0.7071067811865475*ucMSelf[1]*mnu+0.7071067811865475*ucMOther[1]*mnu+4.242640687119286*m0r[1]*mnu-1.414213562373095*cEr[1]*mnu); 
  gkyl_mat_set(A,15,7,0.7071067811865475*ucMSelf[0]*mnu+0.7071067811865475*ucMOther[0]*mnu+4.242640687119286*m0r[0]*mnu-1.414213562373095*cEr[0]*mnu); 
 
  double uM1Self[2] = {0.0}; 
  double uM1Other[2] = {0.0}; 
  double uSumSq[2] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 2*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.7071067811865475*m1r[5]*u_self[a0+1]+0.7071067811865475*m1r[4]*u_self[a0]; 
    uM1Self[1] += 0.7071067811865475*m1r[4]*u_self[a0+1]+0.7071067811865475*m1r[5]*u_self[a0]; 
    uM1Other[0] += 0.7071067811865475*m1r[5]*u_other[a0+1]+0.7071067811865475*m1r[4]*u_other[a0]; 
    uM1Other[1] += 0.7071067811865475*m1r[4]*u_other[a0+1]+0.7071067811865475*m1r[5]*u_other[a0]; 
  const double u_self0R2 = pow(u_self[a0],2);
  const double u_self1R2 = pow(u_self[a0+1],2);
  const double u_other0R2 = pow(u_other[a0],2);
  const double u_other1R2 = pow(u_other[a0+1],2);

  uSumSq[0] += 0.7071067811865475*u_self1R2-1.414213562373095*u_other[a0+1]*u_self[a0+1]+0.7071067811865475*u_other1R2+0.7071067811865475*u_self0R2-1.414213562373095*u_other[a0]*u_self[a0]+0.7071067811865475*u_other0R2; 
  uSumSq[1] += 1.414213562373095*u_self[a0]*u_self[a0+1]-1.414213562373095*u_other[a0]*u_self[a0+1]-1.414213562373095*u_self[a0]*u_other[a0+1]+1.414213562373095*u_other[a0]*u_other[a0+1]; 
  } 
 
  double enRHS[2] = {0.0}; 
  enRHS[0] = (-(6.0*greeneFac[0]*m0r[1]*vtsq_self[1]*m_self*mnu)/(4.0*m_self+4.0*m_other))-(6.0*m0r[0]*greeneFac[1]*vtsq_self[1]*m_self*mnu)/(4.0*m_self+4.0*m_other)-(1.0*greeneFac[0]*m0r[1]*uSumSq[1]*m_self*mnu)/(4.0*m_self+4.0*m_other)-(1.0*m0r[0]*greeneFac[1]*uSumSq[1]*m_self*mnu)/(4.0*m_self+4.0*m_other)-(6.0*vtsq_self[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(4.0*m_self+4.0*m_other)-(1.0*uSumSq[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(4.0*m_self+4.0*m_other)-(6.0*greeneFac[0]*m0r[0]*vtsq_self[0]*m_self*mnu)/(4.0*m_self+4.0*m_other)-(1.0*greeneFac[0]*m0r[0]*uSumSq[0]*m_self*mnu)/(4.0*m_self+4.0*m_other)+(6.0*greeneFac[0]*m0r[1]*vtsq_other[1]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(6.0*m0r[0]*greeneFac[1]*vtsq_other[1]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(greeneFac[0]*m0r[1]*uSumSq[1]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(m0r[0]*greeneFac[1]*uSumSq[1]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(6.0*vtsq_other[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(uSumSq[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(6.0*greeneFac[0]*m0r[0]*vtsq_other[0]*m_other*mnu)/(4.0*m_self+4.0*m_other)+(greeneFac[0]*m0r[0]*uSumSq[0]*m_other*mnu)/(4.0*m_self+4.0*m_other)-1.0*uM1Self[0]*mnu-1.0*uM1Other[0]*mnu+2.0*m2r[0]*mnu; 
  enRHS[1] = (-(54.0*greeneFac[1]*m0r[1]*vtsq_self[1]*m_self*mnu)/(20.0*m_self+20.0*m_other))-(30.0*greeneFac[0]*m0r[0]*vtsq_self[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)-(9.0*greeneFac[1]*m0r[1]*uSumSq[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)-(5.0*greeneFac[0]*m0r[0]*uSumSq[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)-(30.0*greeneFac[0]*vtsq_self[0]*m0r[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)-(5.0*greeneFac[0]*uSumSq[0]*m0r[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)-(30.0*m0r[0]*vtsq_self[0]*greeneFac[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)-(5.0*m0r[0]*uSumSq[0]*greeneFac[1]*m_self*mnu)/(20.0*m_self+20.0*m_other)+(54.0*greeneFac[1]*m0r[1]*vtsq_other[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(30.0*greeneFac[0]*m0r[0]*vtsq_other[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(9.0*greeneFac[1]*m0r[1]*uSumSq[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(5.0*greeneFac[0]*m0r[0]*uSumSq[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(30.0*greeneFac[0]*vtsq_other[0]*m0r[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(5.0*greeneFac[0]*uSumSq[0]*m0r[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(30.0*m0r[0]*vtsq_other[0]*greeneFac[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)+(5.0*m0r[0]*uSumSq[0]*greeneFac[1]*m_other*mnu)/(20.0*m_self+20.0*m_other)-1.0*uM1Self[1]*mnu-1.0*uM1Other[1]*mnu+2.0*m2r[1]*mnu; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,momRHS[3]); 
  gkyl_mat_set(rhs,4,0,momRHS[4]); 
  gkyl_mat_set(rhs,5,0,momRHS[5]); 
  gkyl_mat_set(rhs,6,0,enRHS[0]); 
  gkyl_mat_set(rhs,7,0,enRHS[1]); 
} 
 
