#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *greeneFac, const double m_self, const double *u_self, const double *vtsq_self, const double m_other, const double *u_other, const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections) 
{ 
  // greeneFac:            free parameter beta+1 multiplied by other factors. 
  // nu, m:                collisionality and mass. 
  // moms:                 moments of the distribution function. 
  // u,vtSq:               self primitive moments: mean flow velocity and thermal speed squared. 
  // boundary_corrections: corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:     cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.5*(3.0*moms[3]-1.732050807568877*(moms[2]+moms[1])+moms[0]) < 0) cellAvg = true; 
  if (-0.5*(3.0*moms[3]+1.732050807568877*moms[2]-1.732050807568877*moms[1]-1.0*moms[0]) < 0) cellAvg = true; 
  if (-0.5*(3.0*moms[3]-1.732050807568877*moms[2]+1.732050807568877*moms[1]-1.0*moms[0]) < 0) cellAvg = true; 
  if (0.5*(3.0*moms[3]+1.732050807568877*(moms[2]+moms[1])+moms[0]) < 0) cellAvg = true; 
 
  double m0r[4] = {0.0}; 
  double m1r[8] = {0.0}; 
  double m2r[4] = {0.0}; 
  double cMr[8] = {0.0}; 
  double cEr[4] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = moms[4]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = moms[8]; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    m2r[0] = moms[12]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    cEr[0] = boundary_corrections[8]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m0r[3] = moms[3]; 
    m1r[0] = moms[4]; 
    m1r[1] = moms[5]; 
    m1r[2] = moms[6]; 
    m1r[3] = moms[7]; 
    m1r[4] = moms[8]; 
    m1r[5] = moms[9]; 
    m1r[6] = moms[10]; 
    m1r[7] = moms[11]; 
    m2r[0] = moms[12]; 
    m2r[1] = moms[13]; 
    m2r[2] = moms[14]; 
    m2r[3] = moms[15]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = boundary_corrections[7]; 
    cEr[0] = boundary_corrections[8]; 
    cEr[1] = boundary_corrections[9]; 
    cEr[2] = boundary_corrections[10]; 
    cEr[3] = boundary_corrections[11]; 
  } 
 
  double mnu = 0.5*nu[0]*m_self; 
  double momRHS[8] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,m0r[0]*mnu); 
  gkyl_mat_set(A,0,1,m0r[1]*mnu); 
  gkyl_mat_set(A,0,2,m0r[2]*mnu); 
  gkyl_mat_set(A,0,3,m0r[3]*mnu); 
  gkyl_mat_set(A,1,0,m0r[1]*mnu); 
  gkyl_mat_set(A,1,1,m0r[0]*mnu); 
  gkyl_mat_set(A,1,2,m0r[3]*mnu); 
  gkyl_mat_set(A,1,3,m0r[2]*mnu); 
  gkyl_mat_set(A,2,0,m0r[2]*mnu); 
  gkyl_mat_set(A,2,1,m0r[3]*mnu); 
  gkyl_mat_set(A,2,2,m0r[0]*mnu); 
  gkyl_mat_set(A,2,3,m0r[1]*mnu); 
  gkyl_mat_set(A,3,0,m0r[3]*mnu); 
  gkyl_mat_set(A,3,1,m0r[2]*mnu); 
  gkyl_mat_set(A,3,2,m0r[1]*mnu); 
  gkyl_mat_set(A,3,3,m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,8,-1.0*cMr[0]*mnu); 
  gkyl_mat_set(A,0,9,-1.0*cMr[1]*mnu); 
  gkyl_mat_set(A,0,10,-1.0*cMr[2]*mnu); 
  gkyl_mat_set(A,0,11,-1.0*cMr[3]*mnu); 
  gkyl_mat_set(A,1,8,-1.0*cMr[1]*mnu); 
  gkyl_mat_set(A,1,9,-1.0*cMr[0]*mnu); 
  gkyl_mat_set(A,1,10,-1.0*cMr[3]*mnu); 
  gkyl_mat_set(A,1,11,-1.0*cMr[2]*mnu); 
  gkyl_mat_set(A,2,8,-1.0*cMr[2]*mnu); 
  gkyl_mat_set(A,2,9,-1.0*cMr[3]*mnu); 
  gkyl_mat_set(A,2,10,-1.0*cMr[0]*mnu); 
  gkyl_mat_set(A,2,11,-1.0*cMr[1]*mnu); 
  gkyl_mat_set(A,3,8,-1.0*cMr[3]*mnu); 
  gkyl_mat_set(A,3,9,-1.0*cMr[2]*mnu); 
  gkyl_mat_set(A,3,10,-1.0*cMr[1]*mnu); 
  gkyl_mat_set(A,3,11,-1.0*cMr[0]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,8,0,(-0.25*m0r[3]*u_self[3]*mnu)-0.25*m0r[3]*u_other[3]*mnu-0.25*m0r[2]*u_self[2]*mnu-0.25*m0r[2]*u_other[2]*mnu-0.25*m0r[1]*u_self[1]*mnu-0.25*m0r[1]*u_other[1]*mnu-0.25*m0r[0]*u_self[0]*mnu-0.25*m0r[0]*u_other[0]*mnu+m1r[0]*mnu); 
  gkyl_mat_set(A,8,1,(-0.25*m0r[2]*u_self[3]*mnu)-0.25*m0r[2]*u_other[3]*mnu-0.25*u_self[2]*m0r[3]*mnu-0.25*u_other[2]*m0r[3]*mnu-0.25*m0r[0]*u_self[1]*mnu-0.25*m0r[0]*u_other[1]*mnu+m1r[1]*mnu-0.25*u_self[0]*m0r[1]*mnu-0.25*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,8,2,(-0.25*m0r[1]*u_self[3]*mnu)-0.25*m0r[1]*u_other[3]*mnu-0.25*u_self[1]*m0r[3]*mnu-0.25*u_other[1]*m0r[3]*mnu-0.25*m0r[0]*u_self[2]*mnu-0.25*m0r[0]*u_other[2]*mnu+m1r[2]*mnu-0.25*u_self[0]*m0r[2]*mnu-0.25*u_other[0]*m0r[2]*mnu); 
  gkyl_mat_set(A,8,3,(-0.25*m0r[0]*u_self[3]*mnu)-0.25*m0r[0]*u_other[3]*mnu+m1r[3]*mnu-0.25*u_self[0]*m0r[3]*mnu-0.25*u_other[0]*m0r[3]*mnu-0.25*m0r[1]*u_self[2]*mnu-0.25*m0r[1]*u_other[2]*mnu-0.25*u_self[1]*m0r[2]*mnu-0.25*u_other[1]*m0r[2]*mnu); 
  gkyl_mat_set(A,9,0,(-0.25*m0r[2]*u_self[3]*mnu)-0.25*m0r[2]*u_other[3]*mnu-0.25*u_self[2]*m0r[3]*mnu-0.25*u_other[2]*m0r[3]*mnu-0.25*m0r[0]*u_self[1]*mnu-0.25*m0r[0]*u_other[1]*mnu+m1r[1]*mnu-0.25*u_self[0]*m0r[1]*mnu-0.25*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,9,1,(-0.45*m0r[3]*u_self[3]*mnu)-0.45*m0r[3]*u_other[3]*mnu-0.25*m0r[2]*u_self[2]*mnu-0.25*m0r[2]*u_other[2]*mnu-0.45*m0r[1]*u_self[1]*mnu-0.45*m0r[1]*u_other[1]*mnu-0.25*m0r[0]*u_self[0]*mnu-0.25*m0r[0]*u_other[0]*mnu+m1r[0]*mnu); 
  gkyl_mat_set(A,9,2,(-0.25*m0r[0]*u_self[3]*mnu)-0.25*m0r[0]*u_other[3]*mnu+m1r[3]*mnu-0.25*u_self[0]*m0r[3]*mnu-0.25*u_other[0]*m0r[3]*mnu-0.25*m0r[1]*u_self[2]*mnu-0.25*m0r[1]*u_other[2]*mnu-0.25*u_self[1]*m0r[2]*mnu-0.25*u_other[1]*m0r[2]*mnu); 
  gkyl_mat_set(A,9,3,(-0.45*m0r[1]*u_self[3]*mnu)-0.45*m0r[1]*u_other[3]*mnu-0.45*u_self[1]*m0r[3]*mnu-0.45*u_other[1]*m0r[3]*mnu-0.25*m0r[0]*u_self[2]*mnu-0.25*m0r[0]*u_other[2]*mnu+m1r[2]*mnu-0.25*u_self[0]*m0r[2]*mnu-0.25*u_other[0]*m0r[2]*mnu); 
  gkyl_mat_set(A,10,0,(-0.25*m0r[1]*u_self[3]*mnu)-0.25*m0r[1]*u_other[3]*mnu-0.25*u_self[1]*m0r[3]*mnu-0.25*u_other[1]*m0r[3]*mnu-0.25*m0r[0]*u_self[2]*mnu-0.25*m0r[0]*u_other[2]*mnu+m1r[2]*mnu-0.25*u_self[0]*m0r[2]*mnu-0.25*u_other[0]*m0r[2]*mnu); 
  gkyl_mat_set(A,10,1,(-0.25*m0r[0]*u_self[3]*mnu)-0.25*m0r[0]*u_other[3]*mnu+m1r[3]*mnu-0.25*u_self[0]*m0r[3]*mnu-0.25*u_other[0]*m0r[3]*mnu-0.25*m0r[1]*u_self[2]*mnu-0.25*m0r[1]*u_other[2]*mnu-0.25*u_self[1]*m0r[2]*mnu-0.25*u_other[1]*m0r[2]*mnu); 
  gkyl_mat_set(A,10,2,(-0.45*m0r[3]*u_self[3]*mnu)-0.45*m0r[3]*u_other[3]*mnu-0.45*m0r[2]*u_self[2]*mnu-0.45*m0r[2]*u_other[2]*mnu-0.25*m0r[1]*u_self[1]*mnu-0.25*m0r[1]*u_other[1]*mnu-0.25*m0r[0]*u_self[0]*mnu-0.25*m0r[0]*u_other[0]*mnu+m1r[0]*mnu); 
  gkyl_mat_set(A,10,3,(-0.45*m0r[2]*u_self[3]*mnu)-0.45*m0r[2]*u_other[3]*mnu-0.45*u_self[2]*m0r[3]*mnu-0.45*u_other[2]*m0r[3]*mnu-0.25*m0r[0]*u_self[1]*mnu-0.25*m0r[0]*u_other[1]*mnu+m1r[1]*mnu-0.25*u_self[0]*m0r[1]*mnu-0.25*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,11,0,(-0.25*m0r[0]*u_self[3]*mnu)-0.25*m0r[0]*u_other[3]*mnu+m1r[3]*mnu-0.25*u_self[0]*m0r[3]*mnu-0.25*u_other[0]*m0r[3]*mnu-0.25*m0r[1]*u_self[2]*mnu-0.25*m0r[1]*u_other[2]*mnu-0.25*u_self[1]*m0r[2]*mnu-0.25*u_other[1]*m0r[2]*mnu); 
  gkyl_mat_set(A,11,1,(-0.45*m0r[1]*u_self[3]*mnu)-0.45*m0r[1]*u_other[3]*mnu-0.45*u_self[1]*m0r[3]*mnu-0.45*u_other[1]*m0r[3]*mnu-0.25*m0r[0]*u_self[2]*mnu-0.25*m0r[0]*u_other[2]*mnu+m1r[2]*mnu-0.25*u_self[0]*m0r[2]*mnu-0.25*u_other[0]*m0r[2]*mnu); 
  gkyl_mat_set(A,11,2,(-0.45*m0r[2]*u_self[3]*mnu)-0.45*m0r[2]*u_other[3]*mnu-0.45*u_self[2]*m0r[3]*mnu-0.45*u_other[2]*m0r[3]*mnu-0.25*m0r[0]*u_self[1]*mnu-0.25*m0r[0]*u_other[1]*mnu+m1r[1]*mnu-0.25*u_self[0]*m0r[1]*mnu-0.25*u_other[0]*m0r[1]*mnu); 
  gkyl_mat_set(A,11,3,(-0.81*m0r[3]*u_self[3]*mnu)-0.81*m0r[3]*u_other[3]*mnu-0.45*m0r[2]*u_self[2]*mnu-0.45*m0r[2]*u_other[2]*mnu-0.45*m0r[1]*u_self[1]*mnu-0.45*m0r[1]*u_other[1]*mnu-0.25*m0r[0]*u_self[0]*mnu-0.25*m0r[0]*u_other[0]*mnu+m1r[0]*mnu); 
 
  momRHS[0] += -0.25*((greeneFac[0]*m0r[3]+m0r[0]*greeneFac[3]+greeneFac[1]*m0r[2]+m0r[1]*greeneFac[2])*u_self[3]+((-1.0*greeneFac[0]*m0r[3])-1.0*m0r[0]*greeneFac[3]-1.0*greeneFac[1]*m0r[2]-1.0*m0r[1]*greeneFac[2])*u_other[3]+((u_self[0]-1.0*u_other[0])*greeneFac[3]+greeneFac[1]*u_self[2]-1.0*greeneFac[1]*u_other[2]+(u_self[1]-1.0*u_other[1])*greeneFac[2])*m0r[3]+(m0r[1]*u_self[2]-1.0*m0r[1]*u_other[2]+(u_self[1]-1.0*u_other[1])*m0r[2])*greeneFac[3]+(greeneFac[0]*m0r[2]+m0r[0]*greeneFac[2])*u_self[2]+((-1.0*greeneFac[0]*m0r[2])-1.0*m0r[0]*greeneFac[2])*u_other[2]+(u_self[0]-1.0*u_other[0])*greeneFac[2]*m0r[2]+(greeneFac[0]*m0r[1]+m0r[0]*greeneFac[1])*u_self[1]+((-1.0*greeneFac[0]*m0r[1])-1.0*m0r[0]*greeneFac[1])*u_other[1]+(u_self[0]-1.0*u_other[0])*greeneFac[1]*m0r[1]+greeneFac[0]*m0r[0]*u_self[0]-1.0*greeneFac[0]*m0r[0]*u_other[0]-8.0*m1r[0])*mnu; 
  momRHS[1] += -0.05*((9.0*greeneFac[1]*m0r[3]+9.0*m0r[1]*greeneFac[3]+5.0*greeneFac[0]*m0r[2]+5.0*m0r[0]*greeneFac[2])*u_self[3]+((-9.0*greeneFac[1]*m0r[3])-9.0*m0r[1]*greeneFac[3]-5.0*greeneFac[0]*m0r[2]-5.0*m0r[0]*greeneFac[2])*u_other[3]+((9.0*u_self[1]-9.0*u_other[1])*greeneFac[3]+5.0*greeneFac[0]*u_self[2]-5.0*greeneFac[0]*u_other[2]+(5.0*u_self[0]-5.0*u_other[0])*greeneFac[2])*m0r[3]+(5.0*m0r[0]*u_self[2]-5.0*m0r[0]*u_other[2]+(5.0*u_self[0]-5.0*u_other[0])*m0r[2])*greeneFac[3]+(5.0*greeneFac[1]*m0r[2]+5.0*m0r[1]*greeneFac[2])*u_self[2]+((-5.0*greeneFac[1]*m0r[2])-5.0*m0r[1]*greeneFac[2])*u_other[2]+(5.0*u_self[1]-5.0*u_other[1])*greeneFac[2]*m0r[2]+(9.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[1]+((-9.0*greeneFac[1]*m0r[1])-5.0*greeneFac[0]*m0r[0])*u_other[1]-40.0*m1r[1]+(5.0*greeneFac[0]*u_self[0]-5.0*greeneFac[0]*u_other[0])*m0r[1]+(5.0*m0r[0]*u_self[0]-5.0*m0r[0]*u_other[0])*greeneFac[1])*mnu; 
  momRHS[2] += -0.05*((9.0*greeneFac[2]*m0r[3]+9.0*m0r[2]*greeneFac[3]+5.0*greeneFac[0]*m0r[1]+5.0*m0r[0]*greeneFac[1])*u_self[3]+((-9.0*greeneFac[2]*m0r[3])-9.0*m0r[2]*greeneFac[3]-5.0*greeneFac[0]*m0r[1]-5.0*m0r[0]*greeneFac[1])*u_other[3]+((9.0*u_self[2]-9.0*u_other[2])*greeneFac[3]+5.0*greeneFac[0]*u_self[1]-5.0*greeneFac[0]*u_other[1]+(5.0*u_self[0]-5.0*u_other[0])*greeneFac[1])*m0r[3]+(5.0*m0r[0]*u_self[1]-5.0*m0r[0]*u_other[1]+(5.0*u_self[0]-5.0*u_other[0])*m0r[1])*greeneFac[3]+(9.0*greeneFac[2]*m0r[2]+5.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[2]+((-9.0*greeneFac[2]*m0r[2])-5.0*greeneFac[1]*m0r[1]-5.0*greeneFac[0]*m0r[0])*u_other[2]-40.0*m1r[2]+(5.0*greeneFac[1]*u_self[1]-5.0*greeneFac[1]*u_other[1]+5.0*greeneFac[0]*u_self[0]-5.0*greeneFac[0]*u_other[0])*m0r[2]+(5.0*m0r[1]*u_self[1]-5.0*m0r[1]*u_other[1]+5.0*m0r[0]*u_self[0]-5.0*m0r[0]*u_other[0])*greeneFac[2])*mnu; 
  momRHS[3] += -0.01*((81.0*greeneFac[3]*m0r[3]+45.0*greeneFac[2]*m0r[2]+45.0*greeneFac[1]*m0r[1]+25.0*greeneFac[0]*m0r[0])*u_self[3]+((-81.0*greeneFac[3]*m0r[3])-45.0*greeneFac[2]*m0r[2]-45.0*greeneFac[1]*m0r[1]-25.0*greeneFac[0]*m0r[0])*u_other[3]-200.0*m1r[3]+(45.0*greeneFac[2]*u_self[2]-45.0*greeneFac[2]*u_other[2]+45.0*greeneFac[1]*u_self[1]-45.0*greeneFac[1]*u_other[1]+25.0*greeneFac[0]*u_self[0]-25.0*greeneFac[0]*u_other[0])*m0r[3]+(45.0*m0r[2]*u_self[2]-45.0*m0r[2]*u_other[2]+45.0*m0r[1]*u_self[1]-45.0*m0r[1]*u_other[1]+25.0*m0r[0]*u_self[0]-25.0*m0r[0]*u_other[0])*greeneFac[3]+(25.0*greeneFac[0]*m0r[1]+25.0*m0r[0]*greeneFac[1])*u_self[2]+((-25.0*greeneFac[0]*m0r[1])-25.0*m0r[0]*greeneFac[1])*u_other[2]+(25.0*greeneFac[0]*u_self[1]-25.0*greeneFac[0]*u_other[1]+(25.0*u_self[0]-25.0*u_other[0])*greeneFac[1])*m0r[2]+(25.0*m0r[0]*u_self[1]-25.0*m0r[0]*u_other[1]+(25.0*u_self[0]-25.0*u_other[0])*m0r[1])*greeneFac[2])*mnu; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossY ... // 
  gkyl_mat_set(A,4,4,m0r[0]*mnu); 
  gkyl_mat_set(A,4,5,m0r[1]*mnu); 
  gkyl_mat_set(A,4,6,m0r[2]*mnu); 
  gkyl_mat_set(A,4,7,m0r[3]*mnu); 
  gkyl_mat_set(A,5,4,m0r[1]*mnu); 
  gkyl_mat_set(A,5,5,m0r[0]*mnu); 
  gkyl_mat_set(A,5,6,m0r[3]*mnu); 
  gkyl_mat_set(A,5,7,m0r[2]*mnu); 
  gkyl_mat_set(A,6,4,m0r[2]*mnu); 
  gkyl_mat_set(A,6,5,m0r[3]*mnu); 
  gkyl_mat_set(A,6,6,m0r[0]*mnu); 
  gkyl_mat_set(A,6,7,m0r[1]*mnu); 
  gkyl_mat_set(A,7,4,m0r[3]*mnu); 
  gkyl_mat_set(A,7,5,m0r[2]*mnu); 
  gkyl_mat_set(A,7,6,m0r[1]*mnu); 
  gkyl_mat_set(A,7,7,m0r[0]*mnu); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,4,8,-1.0*cMr[4]*mnu); 
  gkyl_mat_set(A,4,9,-1.0*cMr[5]*mnu); 
  gkyl_mat_set(A,4,10,-1.0*cMr[6]*mnu); 
  gkyl_mat_set(A,4,11,-1.0*cMr[7]*mnu); 
  gkyl_mat_set(A,5,8,-1.0*cMr[5]*mnu); 
  gkyl_mat_set(A,5,9,-1.0*cMr[4]*mnu); 
  gkyl_mat_set(A,5,10,-1.0*cMr[7]*mnu); 
  gkyl_mat_set(A,5,11,-1.0*cMr[6]*mnu); 
  gkyl_mat_set(A,6,8,-1.0*cMr[6]*mnu); 
  gkyl_mat_set(A,6,9,-1.0*cMr[7]*mnu); 
  gkyl_mat_set(A,6,10,-1.0*cMr[4]*mnu); 
  gkyl_mat_set(A,6,11,-1.0*cMr[5]*mnu); 
  gkyl_mat_set(A,7,8,-1.0*cMr[7]*mnu); 
  gkyl_mat_set(A,7,9,-1.0*cMr[6]*mnu); 
  gkyl_mat_set(A,7,10,-1.0*cMr[5]*mnu); 
  gkyl_mat_set(A,7,11,-1.0*cMr[4]*mnu); 
 
  // ... Block from weak multiply of m_self, nu, m1Y and uCrossY ... // 
  gkyl_mat_set(A,8,4,(-0.25*m0r[3]*u_self[7]*mnu)-0.25*m0r[3]*u_other[7]*mnu-0.25*m0r[2]*u_self[6]*mnu-0.25*m0r[2]*u_other[6]*mnu-0.25*m0r[1]*u_self[5]*mnu-0.25*m0r[1]*u_other[5]*mnu-0.25*m0r[0]*u_self[4]*mnu-0.25*m0r[0]*u_other[4]*mnu+m1r[4]*mnu); 
  gkyl_mat_set(A,8,5,(-0.25*m0r[2]*u_self[7]*mnu)-0.25*m0r[2]*u_other[7]*mnu-0.25*m0r[3]*u_self[6]*mnu-0.25*m0r[3]*u_other[6]*mnu-0.25*m0r[0]*u_self[5]*mnu-0.25*m0r[0]*u_other[5]*mnu+m1r[5]*mnu-0.25*m0r[1]*u_self[4]*mnu-0.25*m0r[1]*u_other[4]*mnu); 
  gkyl_mat_set(A,8,6,(-0.25*m0r[1]*u_self[7]*mnu)-0.25*m0r[1]*u_other[7]*mnu-0.25*m0r[0]*u_self[6]*mnu-0.25*m0r[0]*u_other[6]*mnu+m1r[6]*mnu-0.25*m0r[3]*u_self[5]*mnu-0.25*m0r[3]*u_other[5]*mnu-0.25*m0r[2]*u_self[4]*mnu-0.25*m0r[2]*u_other[4]*mnu); 
  gkyl_mat_set(A,8,7,(-0.25*m0r[0]*u_self[7]*mnu)-0.25*m0r[0]*u_other[7]*mnu+m1r[7]*mnu-0.25*m0r[1]*u_self[6]*mnu-0.25*m0r[1]*u_other[6]*mnu-0.25*m0r[2]*u_self[5]*mnu-0.25*m0r[2]*u_other[5]*mnu-0.25*m0r[3]*u_self[4]*mnu-0.25*m0r[3]*u_other[4]*mnu); 
  gkyl_mat_set(A,9,4,(-0.25*m0r[2]*u_self[7]*mnu)-0.25*m0r[2]*u_other[7]*mnu-0.25*m0r[3]*u_self[6]*mnu-0.25*m0r[3]*u_other[6]*mnu-0.25*m0r[0]*u_self[5]*mnu-0.25*m0r[0]*u_other[5]*mnu+m1r[5]*mnu-0.25*m0r[1]*u_self[4]*mnu-0.25*m0r[1]*u_other[4]*mnu); 
  gkyl_mat_set(A,9,5,(-0.45*m0r[3]*u_self[7]*mnu)-0.45*m0r[3]*u_other[7]*mnu-0.25*m0r[2]*u_self[6]*mnu-0.25*m0r[2]*u_other[6]*mnu-0.45*m0r[1]*u_self[5]*mnu-0.45*m0r[1]*u_other[5]*mnu-0.25*m0r[0]*u_self[4]*mnu-0.25*m0r[0]*u_other[4]*mnu+m1r[4]*mnu); 
  gkyl_mat_set(A,9,6,(-0.25*m0r[0]*u_self[7]*mnu)-0.25*m0r[0]*u_other[7]*mnu+m1r[7]*mnu-0.25*m0r[1]*u_self[6]*mnu-0.25*m0r[1]*u_other[6]*mnu-0.25*m0r[2]*u_self[5]*mnu-0.25*m0r[2]*u_other[5]*mnu-0.25*m0r[3]*u_self[4]*mnu-0.25*m0r[3]*u_other[4]*mnu); 
  gkyl_mat_set(A,9,7,(-0.45*m0r[1]*u_self[7]*mnu)-0.45*m0r[1]*u_other[7]*mnu-0.25*m0r[0]*u_self[6]*mnu-0.25*m0r[0]*u_other[6]*mnu+m1r[6]*mnu-0.45*m0r[3]*u_self[5]*mnu-0.45*m0r[3]*u_other[5]*mnu-0.25*m0r[2]*u_self[4]*mnu-0.25*m0r[2]*u_other[4]*mnu); 
  gkyl_mat_set(A,10,4,(-0.25*m0r[1]*u_self[7]*mnu)-0.25*m0r[1]*u_other[7]*mnu-0.25*m0r[0]*u_self[6]*mnu-0.25*m0r[0]*u_other[6]*mnu+m1r[6]*mnu-0.25*m0r[3]*u_self[5]*mnu-0.25*m0r[3]*u_other[5]*mnu-0.25*m0r[2]*u_self[4]*mnu-0.25*m0r[2]*u_other[4]*mnu); 
  gkyl_mat_set(A,10,5,(-0.25*m0r[0]*u_self[7]*mnu)-0.25*m0r[0]*u_other[7]*mnu+m1r[7]*mnu-0.25*m0r[1]*u_self[6]*mnu-0.25*m0r[1]*u_other[6]*mnu-0.25*m0r[2]*u_self[5]*mnu-0.25*m0r[2]*u_other[5]*mnu-0.25*m0r[3]*u_self[4]*mnu-0.25*m0r[3]*u_other[4]*mnu); 
  gkyl_mat_set(A,10,6,(-0.45*m0r[3]*u_self[7]*mnu)-0.45*m0r[3]*u_other[7]*mnu-0.45*m0r[2]*u_self[6]*mnu-0.45*m0r[2]*u_other[6]*mnu-0.25*m0r[1]*u_self[5]*mnu-0.25*m0r[1]*u_other[5]*mnu-0.25*m0r[0]*u_self[4]*mnu-0.25*m0r[0]*u_other[4]*mnu+m1r[4]*mnu); 
  gkyl_mat_set(A,10,7,(-0.45*m0r[2]*u_self[7]*mnu)-0.45*m0r[2]*u_other[7]*mnu-0.45*m0r[3]*u_self[6]*mnu-0.45*m0r[3]*u_other[6]*mnu-0.25*m0r[0]*u_self[5]*mnu-0.25*m0r[0]*u_other[5]*mnu+m1r[5]*mnu-0.25*m0r[1]*u_self[4]*mnu-0.25*m0r[1]*u_other[4]*mnu); 
  gkyl_mat_set(A,11,4,(-0.25*m0r[0]*u_self[7]*mnu)-0.25*m0r[0]*u_other[7]*mnu+m1r[7]*mnu-0.25*m0r[1]*u_self[6]*mnu-0.25*m0r[1]*u_other[6]*mnu-0.25*m0r[2]*u_self[5]*mnu-0.25*m0r[2]*u_other[5]*mnu-0.25*m0r[3]*u_self[4]*mnu-0.25*m0r[3]*u_other[4]*mnu); 
  gkyl_mat_set(A,11,5,(-0.45*m0r[1]*u_self[7]*mnu)-0.45*m0r[1]*u_other[7]*mnu-0.25*m0r[0]*u_self[6]*mnu-0.25*m0r[0]*u_other[6]*mnu+m1r[6]*mnu-0.45*m0r[3]*u_self[5]*mnu-0.45*m0r[3]*u_other[5]*mnu-0.25*m0r[2]*u_self[4]*mnu-0.25*m0r[2]*u_other[4]*mnu); 
  gkyl_mat_set(A,11,6,(-0.45*m0r[2]*u_self[7]*mnu)-0.45*m0r[2]*u_other[7]*mnu-0.45*m0r[3]*u_self[6]*mnu-0.45*m0r[3]*u_other[6]*mnu-0.25*m0r[0]*u_self[5]*mnu-0.25*m0r[0]*u_other[5]*mnu+m1r[5]*mnu-0.25*m0r[1]*u_self[4]*mnu-0.25*m0r[1]*u_other[4]*mnu); 
  gkyl_mat_set(A,11,7,(-0.81*m0r[3]*u_self[7]*mnu)-0.81*m0r[3]*u_other[7]*mnu-0.45*m0r[2]*u_self[6]*mnu-0.45*m0r[2]*u_other[6]*mnu-0.45*m0r[1]*u_self[5]*mnu-0.45*m0r[1]*u_other[5]*mnu-0.25*m0r[0]*u_self[4]*mnu-0.25*m0r[0]*u_other[4]*mnu+m1r[4]*mnu); 
 
  momRHS[4] += -0.25*((greeneFac[0]*m0r[3]+m0r[0]*greeneFac[3]+greeneFac[1]*m0r[2]+m0r[1]*greeneFac[2])*u_self[7]+((-1.0*greeneFac[0]*m0r[3])-1.0*m0r[0]*greeneFac[3]-1.0*greeneFac[1]*m0r[2]-1.0*m0r[1]*greeneFac[2])*u_other[7]+(greeneFac[1]*m0r[3]+m0r[1]*greeneFac[3]+greeneFac[0]*m0r[2]+m0r[0]*greeneFac[2])*u_self[6]+((-1.0*greeneFac[1]*m0r[3])-1.0*m0r[1]*greeneFac[3]-1.0*greeneFac[0]*m0r[2]-1.0*m0r[0]*greeneFac[2])*u_other[6]+(greeneFac[2]*m0r[3]+m0r[2]*greeneFac[3]+greeneFac[0]*m0r[1]+m0r[0]*greeneFac[1])*u_self[5]+((-1.0*greeneFac[2]*m0r[3])-1.0*m0r[2]*greeneFac[3]-1.0*greeneFac[0]*m0r[1]-1.0*m0r[0]*greeneFac[1])*u_other[5]+(greeneFac[3]*m0r[3]+greeneFac[2]*m0r[2]+greeneFac[1]*m0r[1]+greeneFac[0]*m0r[0])*u_self[4]+((-1.0*greeneFac[3]*m0r[3])-1.0*greeneFac[2]*m0r[2]-1.0*greeneFac[1]*m0r[1]-1.0*greeneFac[0]*m0r[0])*u_other[4]-8.0*m1r[4])*mnu; 
  momRHS[5] += -0.05*((9.0*greeneFac[1]*m0r[3]+9.0*m0r[1]*greeneFac[3]+5.0*greeneFac[0]*m0r[2]+5.0*m0r[0]*greeneFac[2])*u_self[7]+((-9.0*greeneFac[1]*m0r[3])-9.0*m0r[1]*greeneFac[3]-5.0*greeneFac[0]*m0r[2]-5.0*m0r[0]*greeneFac[2])*u_other[7]+(5.0*greeneFac[0]*m0r[3]+5.0*m0r[0]*greeneFac[3]+5.0*greeneFac[1]*m0r[2]+5.0*m0r[1]*greeneFac[2])*u_self[6]+((-5.0*greeneFac[0]*m0r[3])-5.0*m0r[0]*greeneFac[3]-5.0*greeneFac[1]*m0r[2]-5.0*m0r[1]*greeneFac[2])*u_other[6]+(9.0*greeneFac[3]*m0r[3]+5.0*greeneFac[2]*m0r[2]+9.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[5]+((-9.0*greeneFac[3]*m0r[3])-5.0*greeneFac[2]*m0r[2]-9.0*greeneFac[1]*m0r[1]-5.0*greeneFac[0]*m0r[0])*u_other[5]-40.0*m1r[5]+(5.0*greeneFac[2]*m0r[3]+5.0*m0r[2]*greeneFac[3]+5.0*greeneFac[0]*m0r[1]+5.0*m0r[0]*greeneFac[1])*u_self[4]+((-5.0*greeneFac[2]*m0r[3])-5.0*m0r[2]*greeneFac[3]-5.0*greeneFac[0]*m0r[1]-5.0*m0r[0]*greeneFac[1])*u_other[4])*mnu; 
  momRHS[6] += -0.05*((9.0*greeneFac[2]*m0r[3]+9.0*m0r[2]*greeneFac[3]+5.0*greeneFac[0]*m0r[1]+5.0*m0r[0]*greeneFac[1])*u_self[7]+((-9.0*greeneFac[2]*m0r[3])-9.0*m0r[2]*greeneFac[3]-5.0*greeneFac[0]*m0r[1]-5.0*m0r[0]*greeneFac[1])*u_other[7]+(9.0*greeneFac[3]*m0r[3]+9.0*greeneFac[2]*m0r[2]+5.0*greeneFac[1]*m0r[1]+5.0*greeneFac[0]*m0r[0])*u_self[6]+((-9.0*greeneFac[3]*m0r[3])-9.0*greeneFac[2]*m0r[2]-5.0*greeneFac[1]*m0r[1]-5.0*greeneFac[0]*m0r[0])*u_other[6]-40.0*m1r[6]+(5.0*greeneFac[0]*m0r[3]+5.0*m0r[0]*greeneFac[3]+5.0*greeneFac[1]*m0r[2]+5.0*m0r[1]*greeneFac[2])*u_self[5]+((-5.0*greeneFac[0]*m0r[3])-5.0*m0r[0]*greeneFac[3]-5.0*greeneFac[1]*m0r[2]-5.0*m0r[1]*greeneFac[2])*u_other[5]+(5.0*greeneFac[1]*m0r[3]+5.0*m0r[1]*greeneFac[3]+5.0*greeneFac[0]*m0r[2]+5.0*m0r[0]*greeneFac[2])*u_self[4]+((-5.0*greeneFac[1]*m0r[3])-5.0*m0r[1]*greeneFac[3]-5.0*greeneFac[0]*m0r[2]-5.0*m0r[0]*greeneFac[2])*u_other[4])*mnu; 
  momRHS[7] += -0.01*((81.0*greeneFac[3]*m0r[3]+45.0*greeneFac[2]*m0r[2]+45.0*greeneFac[1]*m0r[1]+25.0*greeneFac[0]*m0r[0])*u_self[7]+((-81.0*greeneFac[3]*m0r[3])-45.0*greeneFac[2]*m0r[2]-45.0*greeneFac[1]*m0r[1]-25.0*greeneFac[0]*m0r[0])*u_other[7]-200.0*m1r[7]+(45.0*greeneFac[2]*m0r[3]+45.0*m0r[2]*greeneFac[3]+25.0*greeneFac[0]*m0r[1]+25.0*m0r[0]*greeneFac[1])*u_self[6]+((-45.0*greeneFac[2]*m0r[3])-45.0*m0r[2]*greeneFac[3]-25.0*greeneFac[0]*m0r[1]-25.0*m0r[0]*greeneFac[1])*u_other[6]+(45.0*greeneFac[1]*m0r[3]+45.0*m0r[1]*greeneFac[3]+25.0*greeneFac[0]*m0r[2]+25.0*m0r[0]*greeneFac[2])*u_self[5]+((-45.0*greeneFac[1]*m0r[3])-45.0*m0r[1]*greeneFac[3]-25.0*greeneFac[0]*m0r[2]-25.0*m0r[0]*greeneFac[2])*u_other[5]+(25.0*greeneFac[0]*m0r[3]+25.0*m0r[0]*greeneFac[3]+25.0*greeneFac[1]*m0r[2]+25.0*m0r[1]*greeneFac[2])*u_self[4]+((-25.0*greeneFac[0]*m0r[3])-25.0*m0r[0]*greeneFac[3]-25.0*greeneFac[1]*m0r[2]-25.0*m0r[1]*greeneFac[2])*u_other[4])*mnu; 
 
  double ucMSelf[4] = {0.0}; 
  double ucMOther[4] = {0.0}; 
  for (int vd=0; vd<2; vd++) 
  { 
    int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMr[a0+3]*u_self[a0+3]+0.5*cMr[a0+2]*u_self[a0+2]+0.5*cMr[a0+1]*u_self[a0+1]+0.5*cMr[a0]*u_self[a0]; 
    ucMSelf[1] += 0.5*cMr[a0+2]*u_self[a0+3]+0.5*u_self[a0+2]*cMr[a0+3]+0.5*cMr[a0]*u_self[a0+1]+0.5*u_self[a0]*cMr[a0+1]; 
    ucMSelf[2] += 0.5*cMr[a0+1]*u_self[a0+3]+0.5*u_self[a0+1]*cMr[a0+3]+0.5*cMr[a0]*u_self[a0+2]+0.5*u_self[a0]*cMr[a0+2]; 
    ucMSelf[3] += 0.5*cMr[a0]*u_self[a0+3]+0.5*u_self[a0]*cMr[a0+3]+0.5*cMr[a0+1]*u_self[a0+2]+0.5*u_self[a0+1]*cMr[a0+2]; 
    ucMOther[0] += 0.5*cMr[a0+3]*u_other[a0+3]+0.5*cMr[a0+2]*u_other[a0+2]+0.5*cMr[a0+1]*u_other[a0+1]+0.5*cMr[a0]*u_other[a0]; 
    ucMOther[1] += 0.5*cMr[a0+2]*u_other[a0+3]+0.5*u_other[a0+2]*cMr[a0+3]+0.5*cMr[a0]*u_other[a0+1]+0.5*u_other[a0]*cMr[a0+1]; 
    ucMOther[2] += 0.5*cMr[a0+1]*u_other[a0+3]+0.5*u_other[a0+1]*cMr[a0+3]+0.5*cMr[a0]*u_other[a0+2]+0.5*u_other[a0]*cMr[a0+2]; 
    ucMOther[3] += 0.5*cMr[a0]*u_other[a0+3]+0.5*u_other[a0]*cMr[a0+3]+0.5*cMr[a0+1]*u_other[a0+2]+0.5*u_other[a0+1]*cMr[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,20,8,0.5*ucMSelf[0]*mnu+0.5*ucMOther[0]*mnu+2.0*m0r[0]*mnu-1.0*cEr[0]*mnu); 
  gkyl_mat_set(A,20,9,0.5*ucMSelf[1]*mnu+0.5*ucMOther[1]*mnu+2.0*m0r[1]*mnu-1.0*cEr[1]*mnu); 
  gkyl_mat_set(A,20,10,0.5*ucMSelf[2]*mnu+0.5*ucMOther[2]*mnu+2.0*m0r[2]*mnu-1.0*cEr[2]*mnu); 
  gkyl_mat_set(A,20,11,0.5*ucMSelf[3]*mnu+0.5*ucMOther[3]*mnu+2.0*m0r[3]*mnu-1.0*cEr[3]*mnu); 
  gkyl_mat_set(A,21,8,0.5*ucMSelf[1]*mnu+0.5*ucMOther[1]*mnu+2.0*m0r[1]*mnu-1.0*cEr[1]*mnu); 
  gkyl_mat_set(A,21,9,0.5*ucMSelf[0]*mnu+0.5*ucMOther[0]*mnu+2.0*m0r[0]*mnu-1.0*cEr[0]*mnu); 
  gkyl_mat_set(A,21,10,0.5*ucMSelf[3]*mnu+0.5*ucMOther[3]*mnu+2.0*m0r[3]*mnu-1.0*cEr[3]*mnu); 
  gkyl_mat_set(A,21,11,0.5*ucMSelf[2]*mnu+0.5*ucMOther[2]*mnu+2.0*m0r[2]*mnu-1.0*cEr[2]*mnu); 
  gkyl_mat_set(A,22,8,0.5*ucMSelf[2]*mnu+0.5*ucMOther[2]*mnu+2.0*m0r[2]*mnu-1.0*cEr[2]*mnu); 
  gkyl_mat_set(A,22,9,0.5*ucMSelf[3]*mnu+0.5*ucMOther[3]*mnu+2.0*m0r[3]*mnu-1.0*cEr[3]*mnu); 
  gkyl_mat_set(A,22,10,0.5*ucMSelf[0]*mnu+0.5*ucMOther[0]*mnu+2.0*m0r[0]*mnu-1.0*cEr[0]*mnu); 
  gkyl_mat_set(A,22,11,0.5*ucMSelf[1]*mnu+0.5*ucMOther[1]*mnu+2.0*m0r[1]*mnu-1.0*cEr[1]*mnu); 
  gkyl_mat_set(A,23,8,0.5*ucMSelf[3]*mnu+0.5*ucMOther[3]*mnu+2.0*m0r[3]*mnu-1.0*cEr[3]*mnu); 
  gkyl_mat_set(A,23,9,0.5*ucMSelf[2]*mnu+0.5*ucMOther[2]*mnu+2.0*m0r[2]*mnu-1.0*cEr[2]*mnu); 
  gkyl_mat_set(A,23,10,0.5*ucMSelf[1]*mnu+0.5*ucMOther[1]*mnu+2.0*m0r[1]*mnu-1.0*cEr[1]*mnu); 
  gkyl_mat_set(A,23,11,0.5*ucMSelf[0]*mnu+0.5*ucMOther[0]*mnu+2.0*m0r[0]*mnu-1.0*cEr[0]*mnu); 
 
  double uM1Self[4] = {0.0}; 
  double uM1Other[4] = {0.0}; 
  double uSumSq[4] = {0.0}; 
  for (int vd=0; vd<2; vd++) 
  { 
    int a0 = 4*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.5*m1r[7]*u_self[a0+3]+0.5*m1r[6]*u_self[a0+2]+0.5*m1r[5]*u_self[a0+1]+0.5*m1r[4]*u_self[a0]; 
    uM1Self[1] += 0.5*m1r[6]*u_self[a0+3]+0.5*m1r[7]*u_self[a0+2]+0.5*m1r[4]*u_self[a0+1]+0.5*m1r[5]*u_self[a0]; 
    uM1Self[2] += 0.5*m1r[5]*u_self[a0+3]+0.5*m1r[4]*u_self[a0+2]+0.5*m1r[7]*u_self[a0+1]+0.5*m1r[6]*u_self[a0]; 
    uM1Self[3] += 0.5*m1r[4]*u_self[a0+3]+0.5*m1r[5]*u_self[a0+2]+0.5*m1r[6]*u_self[a0+1]+0.5*m1r[7]*u_self[a0]; 
    uM1Other[0] += 0.5*m1r[7]*u_other[a0+3]+0.5*m1r[6]*u_other[a0+2]+0.5*m1r[5]*u_other[a0+1]+0.5*m1r[4]*u_other[a0]; 
    uM1Other[1] += 0.5*m1r[6]*u_other[a0+3]+0.5*m1r[7]*u_other[a0+2]+0.5*m1r[4]*u_other[a0+1]+0.5*m1r[5]*u_other[a0]; 
    uM1Other[2] += 0.5*m1r[5]*u_other[a0+3]+0.5*m1r[4]*u_other[a0+2]+0.5*m1r[7]*u_other[a0+1]+0.5*m1r[6]*u_other[a0]; 
    uM1Other[3] += 0.5*m1r[4]*u_other[a0+3]+0.5*m1r[5]*u_other[a0+2]+0.5*m1r[6]*u_other[a0+1]+0.5*m1r[7]*u_other[a0]; 
  const double u_self0R2 = pow(u_self[a0],2);
  const double u_self1R2 = pow(u_self[a0+1],2);
  const double u_self2R2 = pow(u_self[a0+2],2);
  const double u_self3R2 = pow(u_self[a0+3],2);
  const double u_other0R2 = pow(u_other[a0],2);
  const double u_other1R2 = pow(u_other[a0+1],2);
  const double u_other2R2 = pow(u_other[a0+2],2);
  const double u_other3R2 = pow(u_other[a0+3],2);

  uSumSq[0] += 0.5*u_self3R2-1.0*u_other[a0+3]*u_self[a0+3]+0.5*u_other3R2+0.5*u_self2R2-1.0*u_other[a0+2]*u_self[a0+2]+0.5*u_other2R2+0.5*u_self1R2-1.0*u_other[a0+1]*u_self[a0+1]+0.5*u_other1R2+0.5*u_self0R2-1.0*u_other[a0]*u_self[a0]+0.5*u_other0R2; 
  uSumSq[1] += u_self[a0+2]*u_self[a0+3]-1.0*u_other[a0+2]*u_self[a0+3]-1.0*u_self[a0+2]*u_other[a0+3]+u_other[a0+2]*u_other[a0+3]+u_self[a0]*u_self[a0+1]-1.0*u_other[a0]*u_self[a0+1]-1.0*u_self[a0]*u_other[a0+1]+u_other[a0]*u_other[a0+1]; 
  uSumSq[2] += u_self[a0+1]*u_self[a0+3]-1.0*u_other[a0+1]*u_self[a0+3]-1.0*u_self[a0+1]*u_other[a0+3]+u_other[a0+1]*u_other[a0+3]+u_self[a0]*u_self[a0+2]-1.0*u_other[a0]*u_self[a0+2]-1.0*u_self[a0]*u_other[a0+2]+u_other[a0]*u_other[a0+2]; 
  uSumSq[3] += u_self[a0]*u_self[a0+3]-1.0*u_other[a0]*u_self[a0+3]-1.0*u_self[a0]*u_other[a0+3]+u_other[a0]*u_other[a0+3]+u_self[a0+1]*u_self[a0+2]-1.0*u_other[a0+1]*u_self[a0+2]-1.0*u_self[a0+1]*u_other[a0+2]+u_other[a0+1]*u_other[a0+2]; 
  } 
 
  double enRHS[4] = {0.0}; 
  enRHS[0] = (-(4.0*greeneFac[0]*m0r[3]*vtsq_self[3]*m_self*mnu)/(8.0*m_self+8.0*m_other))-(4.0*m0r[0]*greeneFac[3]*vtsq_self[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*greeneFac[1]*m0r[2]*vtsq_self[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*m0r[1]*greeneFac[2]*vtsq_self[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*greeneFac[0]*m0r[3]*uSumSq[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*m0r[0]*greeneFac[3]*uSumSq[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*greeneFac[1]*m0r[2]*uSumSq[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*m0r[1]*greeneFac[2]*uSumSq[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*vtsq_self[0]*greeneFac[3]*m0r[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*uSumSq[0]*greeneFac[3]*m0r[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*greeneFac[1]*vtsq_self[2]*m0r[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*greeneFac[1]*uSumSq[2]*m0r[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*vtsq_self[1]*greeneFac[2]*m0r[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*uSumSq[1]*greeneFac[2]*m0r[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*m0r[1]*vtsq_self[2]*greeneFac[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*m0r[1]*uSumSq[2]*greeneFac[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*vtsq_self[1]*m0r[2]*greeneFac[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*uSumSq[1]*m0r[2]*greeneFac[3]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*greeneFac[0]*m0r[2]*vtsq_self[2]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*m0r[0]*greeneFac[2]*vtsq_self[2]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*greeneFac[0]*m0r[2]*uSumSq[2]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*m0r[0]*greeneFac[2]*uSumSq[2]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*vtsq_self[0]*greeneFac[2]*m0r[2]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*uSumSq[0]*greeneFac[2]*m0r[2]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*greeneFac[0]*m0r[1]*vtsq_self[1]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*m0r[0]*greeneFac[1]*vtsq_self[1]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*greeneFac[0]*m0r[1]*uSumSq[1]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*m0r[0]*greeneFac[1]*uSumSq[1]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*vtsq_self[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*uSumSq[0]*greeneFac[1]*m0r[1]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(4.0*greeneFac[0]*m0r[0]*vtsq_self[0]*m_self*mnu)/(8.0*m_self+8.0*m_other)-(1.0*greeneFac[0]*m0r[0]*uSumSq[0]*m_self*mnu)/(8.0*m_self+8.0*m_other)+(4.0*greeneFac[0]*m0r[3]*vtsq_other[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*m0r[0]*greeneFac[3]*vtsq_other[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*greeneFac[1]*m0r[2]*vtsq_other[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*m0r[1]*greeneFac[2]*vtsq_other[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(greeneFac[0]*m0r[3]*uSumSq[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(m0r[0]*greeneFac[3]*uSumSq[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(greeneFac[1]*m0r[2]*uSumSq[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(m0r[1]*greeneFac[2]*uSumSq[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*vtsq_other[0]*greeneFac[3]*m0r[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(uSumSq[0]*greeneFac[3]*m0r[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*greeneFac[1]*vtsq_other[2]*m0r[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(greeneFac[1]*uSumSq[2]*m0r[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*vtsq_other[1]*greeneFac[2]*m0r[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(uSumSq[1]*greeneFac[2]*m0r[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*m0r[1]*vtsq_other[2]*greeneFac[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(m0r[1]*uSumSq[2]*greeneFac[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*vtsq_other[1]*m0r[2]*greeneFac[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(uSumSq[1]*m0r[2]*greeneFac[3]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*greeneFac[0]*m0r[2]*vtsq_other[2]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*m0r[0]*greeneFac[2]*vtsq_other[2]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(greeneFac[0]*m0r[2]*uSumSq[2]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(m0r[0]*greeneFac[2]*uSumSq[2]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*vtsq_other[0]*greeneFac[2]*m0r[2]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(uSumSq[0]*greeneFac[2]*m0r[2]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*greeneFac[0]*m0r[1]*vtsq_other[1]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*m0r[0]*greeneFac[1]*vtsq_other[1]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(greeneFac[0]*m0r[1]*uSumSq[1]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(m0r[0]*greeneFac[1]*uSumSq[1]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*vtsq_other[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(uSumSq[0]*greeneFac[1]*m0r[1]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(4.0*greeneFac[0]*m0r[0]*vtsq_other[0]*m_other*mnu)/(8.0*m_self+8.0*m_other)+(greeneFac[0]*m0r[0]*uSumSq[0]*m_other*mnu)/(8.0*m_self+8.0*m_other)-1.0*uM1Self[0]*mnu-1.0*uM1Other[0]*mnu+2.0*m2r[0]*mnu; 
  enRHS[1] = (-(36.0*greeneFac[1]*m0r[3]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other))-(36.0*m0r[1]*greeneFac[3]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*m0r[2]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[0]*greeneFac[2]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*greeneFac[1]*m0r[3]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*m0r[1]*greeneFac[3]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*m0r[2]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[0]*greeneFac[2]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(36.0*vtsq_self[1]*greeneFac[3]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*uSumSq[1]*greeneFac[3]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*vtsq_self[2]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*uSumSq[2]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*vtsq_self[0]*greeneFac[2]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*uSumSq[0]*greeneFac[2]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[0]*vtsq_self[2]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[0]*uSumSq[2]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*vtsq_self[0]*m0r[2]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*uSumSq[0]*m0r[2]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[1]*m0r[2]*vtsq_self[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[1]*greeneFac[2]*vtsq_self[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[1]*m0r[2]*uSumSq[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[1]*greeneFac[2]*uSumSq[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*vtsq_self[1]*greeneFac[2]*m0r[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*uSumSq[1]*greeneFac[2]*m0r[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(36.0*greeneFac[1]*m0r[1]*vtsq_self[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*m0r[0]*vtsq_self[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*greeneFac[1]*m0r[1]*uSumSq[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*m0r[0]*uSumSq[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*vtsq_self[0]*m0r[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*uSumSq[0]*m0r[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[0]*vtsq_self[0]*greeneFac[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[0]*uSumSq[0]*greeneFac[1]*m_self*mnu)/(40.0*m_self+40.0*m_other)+(36.0*greeneFac[1]*m0r[3]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(36.0*m0r[1]*greeneFac[3]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*m0r[2]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[0]*greeneFac[2]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*greeneFac[1]*m0r[3]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*m0r[1]*greeneFac[3]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*m0r[2]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[0]*greeneFac[2]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(36.0*vtsq_other[1]*greeneFac[3]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*uSumSq[1]*greeneFac[3]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*vtsq_other[2]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*uSumSq[2]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*vtsq_other[0]*greeneFac[2]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*uSumSq[0]*greeneFac[2]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[0]*vtsq_other[2]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[0]*uSumSq[2]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*vtsq_other[0]*m0r[2]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*uSumSq[0]*m0r[2]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[1]*m0r[2]*vtsq_other[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[1]*greeneFac[2]*vtsq_other[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[1]*m0r[2]*uSumSq[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[1]*greeneFac[2]*uSumSq[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*vtsq_other[1]*greeneFac[2]*m0r[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*uSumSq[1]*greeneFac[2]*m0r[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(36.0*greeneFac[1]*m0r[1]*vtsq_other[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*m0r[0]*vtsq_other[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*greeneFac[1]*m0r[1]*uSumSq[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*m0r[0]*uSumSq[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*vtsq_other[0]*m0r[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*uSumSq[0]*m0r[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[0]*vtsq_other[0]*greeneFac[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[0]*uSumSq[0]*greeneFac[1]*m_other*mnu)/(40.0*m_self+40.0*m_other)-1.0*uM1Self[1]*mnu-1.0*uM1Other[1]*mnu+2.0*m2r[1]*mnu; 
  enRHS[2] = (-(36.0*greeneFac[2]*m0r[3]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other))-(36.0*m0r[2]*greeneFac[3]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*m0r[1]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[0]*greeneFac[1]*vtsq_self[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*greeneFac[2]*m0r[3]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*m0r[2]*greeneFac[3]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*m0r[1]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[0]*greeneFac[1]*uSumSq[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(36.0*vtsq_self[2]*greeneFac[3]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*uSumSq[2]*greeneFac[3]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*vtsq_self[1]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*uSumSq[1]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*vtsq_self[0]*greeneFac[1]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*uSumSq[0]*greeneFac[1]*m0r[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[0]*vtsq_self[1]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[0]*uSumSq[1]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*vtsq_self[0]*m0r[1]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*uSumSq[0]*m0r[1]*greeneFac[3]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(36.0*greeneFac[2]*m0r[2]*vtsq_self[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[1]*m0r[1]*vtsq_self[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*m0r[0]*vtsq_self[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(9.0*greeneFac[2]*m0r[2]*uSumSq[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[1]*m0r[1]*uSumSq[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*m0r[0]*uSumSq[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[1]*vtsq_self[1]*m0r[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[1]*uSumSq[1]*m0r[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*greeneFac[0]*vtsq_self[0]*m0r[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*greeneFac[0]*uSumSq[0]*m0r[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[1]*vtsq_self[1]*greeneFac[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[1]*uSumSq[1]*greeneFac[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(20.0*m0r[0]*vtsq_self[0]*greeneFac[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)-(5.0*m0r[0]*uSumSq[0]*greeneFac[2]*m_self*mnu)/(40.0*m_self+40.0*m_other)+(36.0*greeneFac[2]*m0r[3]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(36.0*m0r[2]*greeneFac[3]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*m0r[1]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[0]*greeneFac[1]*vtsq_other[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*greeneFac[2]*m0r[3]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*m0r[2]*greeneFac[3]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*m0r[1]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[0]*greeneFac[1]*uSumSq[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(36.0*vtsq_other[2]*greeneFac[3]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*uSumSq[2]*greeneFac[3]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*vtsq_other[1]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*uSumSq[1]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*vtsq_other[0]*greeneFac[1]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*uSumSq[0]*greeneFac[1]*m0r[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[0]*vtsq_other[1]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[0]*uSumSq[1]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*vtsq_other[0]*m0r[1]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*uSumSq[0]*m0r[1]*greeneFac[3]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(36.0*greeneFac[2]*m0r[2]*vtsq_other[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[1]*m0r[1]*vtsq_other[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*m0r[0]*vtsq_other[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(9.0*greeneFac[2]*m0r[2]*uSumSq[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[1]*m0r[1]*uSumSq[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*m0r[0]*uSumSq[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[1]*vtsq_other[1]*m0r[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[1]*uSumSq[1]*m0r[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*greeneFac[0]*vtsq_other[0]*m0r[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*greeneFac[0]*uSumSq[0]*m0r[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[1]*vtsq_other[1]*greeneFac[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[1]*uSumSq[1]*greeneFac[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(20.0*m0r[0]*vtsq_other[0]*greeneFac[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)+(5.0*m0r[0]*uSumSq[0]*greeneFac[2]*m_other*mnu)/(40.0*m_self+40.0*m_other)-1.0*uM1Self[2]*mnu-1.0*uM1Other[2]*mnu+2.0*m2r[2]*mnu; 
  enRHS[3] = (-(324.0*greeneFac[3]*m0r[3]*vtsq_self[3]*m_self*mnu)/(200.0*m_self+200.0*m_other))-(180.0*greeneFac[2]*m0r[2]*vtsq_self[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(180.0*greeneFac[1]*m0r[1]*vtsq_self[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*greeneFac[0]*m0r[0]*vtsq_self[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(81.0*greeneFac[3]*m0r[3]*uSumSq[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(45.0*greeneFac[2]*m0r[2]*uSumSq[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(45.0*greeneFac[1]*m0r[1]*uSumSq[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*greeneFac[0]*m0r[0]*uSumSq[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(180.0*greeneFac[2]*vtsq_self[2]*m0r[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(45.0*greeneFac[2]*uSumSq[2]*m0r[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(180.0*greeneFac[1]*vtsq_self[1]*m0r[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(45.0*greeneFac[1]*uSumSq[1]*m0r[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*greeneFac[0]*vtsq_self[0]*m0r[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*greeneFac[0]*uSumSq[0]*m0r[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(180.0*m0r[2]*vtsq_self[2]*greeneFac[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(45.0*m0r[2]*uSumSq[2]*greeneFac[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(180.0*m0r[1]*vtsq_self[1]*greeneFac[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(45.0*m0r[1]*uSumSq[1]*greeneFac[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*m0r[0]*vtsq_self[0]*greeneFac[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*m0r[0]*uSumSq[0]*greeneFac[3]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*greeneFac[0]*m0r[1]*vtsq_self[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*m0r[0]*greeneFac[1]*vtsq_self[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*greeneFac[0]*m0r[1]*uSumSq[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*m0r[0]*greeneFac[1]*uSumSq[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*greeneFac[0]*vtsq_self[1]*m0r[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*greeneFac[0]*uSumSq[1]*m0r[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*vtsq_self[0]*greeneFac[1]*m0r[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*uSumSq[0]*greeneFac[1]*m0r[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*m0r[0]*vtsq_self[1]*greeneFac[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*m0r[0]*uSumSq[1]*greeneFac[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(100.0*vtsq_self[0]*m0r[1]*greeneFac[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)-(25.0*uSumSq[0]*m0r[1]*greeneFac[2]*m_self*mnu)/(200.0*m_self+200.0*m_other)+(324.0*greeneFac[3]*m0r[3]*vtsq_other[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(180.0*greeneFac[2]*m0r[2]*vtsq_other[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(180.0*greeneFac[1]*m0r[1]*vtsq_other[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*greeneFac[0]*m0r[0]*vtsq_other[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(81.0*greeneFac[3]*m0r[3]*uSumSq[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(45.0*greeneFac[2]*m0r[2]*uSumSq[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(45.0*greeneFac[1]*m0r[1]*uSumSq[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*greeneFac[0]*m0r[0]*uSumSq[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(180.0*greeneFac[2]*vtsq_other[2]*m0r[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(45.0*greeneFac[2]*uSumSq[2]*m0r[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(180.0*greeneFac[1]*vtsq_other[1]*m0r[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(45.0*greeneFac[1]*uSumSq[1]*m0r[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*greeneFac[0]*vtsq_other[0]*m0r[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*greeneFac[0]*uSumSq[0]*m0r[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(180.0*m0r[2]*vtsq_other[2]*greeneFac[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(45.0*m0r[2]*uSumSq[2]*greeneFac[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(180.0*m0r[1]*vtsq_other[1]*greeneFac[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(45.0*m0r[1]*uSumSq[1]*greeneFac[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*m0r[0]*vtsq_other[0]*greeneFac[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*m0r[0]*uSumSq[0]*greeneFac[3]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*greeneFac[0]*m0r[1]*vtsq_other[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*m0r[0]*greeneFac[1]*vtsq_other[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*greeneFac[0]*m0r[1]*uSumSq[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*m0r[0]*greeneFac[1]*uSumSq[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*greeneFac[0]*vtsq_other[1]*m0r[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*greeneFac[0]*uSumSq[1]*m0r[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*vtsq_other[0]*greeneFac[1]*m0r[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*uSumSq[0]*greeneFac[1]*m0r[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*m0r[0]*vtsq_other[1]*greeneFac[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*m0r[0]*uSumSq[1]*greeneFac[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(100.0*vtsq_other[0]*m0r[1]*greeneFac[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)+(25.0*uSumSq[0]*m0r[1]*greeneFac[2]*m_other*mnu)/(200.0*m_self+200.0*m_other)-1.0*uM1Self[3]*mnu-1.0*uM1Other[3]*mnu+2.0*m2r[3]*mnu; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,momRHS[3]); 
  gkyl_mat_set(rhs,4,0,momRHS[4]); 
  gkyl_mat_set(rhs,5,0,momRHS[5]); 
  gkyl_mat_set(rhs,6,0,momRHS[6]); 
  gkyl_mat_set(rhs,7,0,momRHS[7]); 
  gkyl_mat_set(rhs,8,0,enRHS[0]); 
  gkyl_mat_set(rhs,9,0,enRHS[1]); 
  gkyl_mat_set(rhs,10,0,enRHS[2]); 
  gkyl_mat_set(rhs,11,0,enRHS[3]); 
} 
 
