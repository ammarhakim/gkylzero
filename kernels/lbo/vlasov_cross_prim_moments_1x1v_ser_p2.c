#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf, const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther, const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *uCross, double *vtSqCross) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.581138830084189*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.581138830084189*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[3]; 
  double m1r[3]; 
  double m2r[3]; 
  double cMr[3]; 
  double cEr[3]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
  } 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double momRHS[3]; 
  // zero out momentum RHS array. 
  for (int vd=0; vd<3; vd++) 
  { 
    momRHS[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.414213562373095*m0r[0]*mnuSelf); 
  gkyl_mat_set(A,0,1,1.414213562373095*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,0,2,1.414213562373095*m0r[2]*mnuSelf); 
  gkyl_mat_set(A,1,0,1.414213562373095*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,1,1,1.264911064067352*m0r[2]*mnuSelf+1.414213562373095*m0r[0]*mnuSelf); 
  gkyl_mat_set(A,1,2,1.264911064067352*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,2,0,1.414213562373095*m0r[2]*mnuSelf); 
  gkyl_mat_set(A,2,1,1.264911064067352*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,2,2,0.9035079029052515*m0r[2]*mnuSelf+1.414213562373095*m0r[0]*mnuSelf); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,3,-1.414213562373095*cMr[0]*mnuSelf); 
  gkyl_mat_set(A,0,4,-1.414213562373095*cMr[1]*mnuSelf); 
  gkyl_mat_set(A,0,5,-1.414213562373095*cMr[2]*mnuSelf); 
  gkyl_mat_set(A,1,3,-1.414213562373095*cMr[1]*mnuSelf); 
  gkyl_mat_set(A,1,4,(-1.264911064067352*cMr[2]*mnuSelf)-1.414213562373095*cMr[0]*mnuSelf); 
  gkyl_mat_set(A,1,5,-1.264911064067352*cMr[1]*mnuSelf); 
  gkyl_mat_set(A,2,3,-1.414213562373095*cMr[2]*mnuSelf); 
  gkyl_mat_set(A,2,4,-1.264911064067352*cMr[1]*mnuSelf); 
  gkyl_mat_set(A,2,5,(-0.9035079029052515*cMr[2]*mnuSelf)-1.414213562373095*cMr[0]*mnuSelf); 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1X and uCrossX ... // 
  gkyl_mat_set(A,3,0,(-0.5*m0r[2]*uSelf[2]*mnuSelf)-0.5*m0r[2]*uOther[2]*mnuSelf-0.5*m0r[1]*uSelf[1]*mnuSelf-0.5*m0r[1]*uOther[1]*mnuSelf-0.5*m0r[0]*uSelf[0]*mnuSelf-0.5*m0r[0]*uOther[0]*mnuSelf+1.414213562373095*m1r[0]*mnuSelf); 
  gkyl_mat_set(A,3,1,(-0.4472135954999579*m0r[1]*uSelf[2]*mnuSelf)-0.4472135954999579*m0r[1]*uOther[2]*mnuSelf-0.4472135954999579*uSelf[1]*m0r[2]*mnuSelf-0.4472135954999579*uOther[1]*m0r[2]*mnuSelf-0.5*m0r[0]*uSelf[1]*mnuSelf-0.5*m0r[0]*uOther[1]*mnuSelf+1.414213562373095*m1r[1]*mnuSelf-0.5*uSelf[0]*m0r[1]*mnuSelf-0.5*uOther[0]*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,3,2,(-0.3194382824999699*m0r[2]*uSelf[2]*mnuSelf)-0.5*m0r[0]*uSelf[2]*mnuSelf-0.3194382824999699*m0r[2]*uOther[2]*mnuSelf-0.5*m0r[0]*uOther[2]*mnuSelf+1.414213562373095*m1r[2]*mnuSelf-0.5*uSelf[0]*m0r[2]*mnuSelf-0.5*uOther[0]*m0r[2]*mnuSelf-0.4472135954999579*m0r[1]*uSelf[1]*mnuSelf-0.4472135954999579*m0r[1]*uOther[1]*mnuSelf); 
  gkyl_mat_set(A,4,0,(-0.4472135954999579*m0r[1]*uSelf[2]*mnuSelf)-0.4472135954999579*m0r[1]*uOther[2]*mnuSelf-0.4472135954999579*uSelf[1]*m0r[2]*mnuSelf-0.4472135954999579*uOther[1]*m0r[2]*mnuSelf-0.5*m0r[0]*uSelf[1]*mnuSelf-0.5*m0r[0]*uOther[1]*mnuSelf+1.414213562373095*m1r[1]*mnuSelf-0.5*uSelf[0]*m0r[1]*mnuSelf-0.5*uOther[0]*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,4,1,(-0.7857142857142857*m0r[2]*uSelf[2]*mnuSelf)-0.4472135954999579*m0r[0]*uSelf[2]*mnuSelf-0.7857142857142857*m0r[2]*uOther[2]*mnuSelf-0.4472135954999579*m0r[0]*uOther[2]*mnuSelf+1.264911064067352*m1r[2]*mnuSelf-0.4472135954999579*uSelf[0]*m0r[2]*mnuSelf-0.4472135954999579*uOther[0]*m0r[2]*mnuSelf-0.9*m0r[1]*uSelf[1]*mnuSelf-0.9*m0r[1]*uOther[1]*mnuSelf-0.5*m0r[0]*uSelf[0]*mnuSelf-0.5*m0r[0]*uOther[0]*mnuSelf+1.414213562373095*m1r[0]*mnuSelf); 
  gkyl_mat_set(A,4,2,(-0.7857142857142857*m0r[1]*uSelf[2]*mnuSelf)-0.7857142857142857*m0r[1]*uOther[2]*mnuSelf-0.7857142857142857*uSelf[1]*m0r[2]*mnuSelf-0.7857142857142857*uOther[1]*m0r[2]*mnuSelf-0.4472135954999579*m0r[0]*uSelf[1]*mnuSelf-0.4472135954999579*m0r[0]*uOther[1]*mnuSelf+1.264911064067352*m1r[1]*mnuSelf-0.4472135954999579*uSelf[0]*m0r[1]*mnuSelf-0.4472135954999579*uOther[0]*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,5,0,(-0.3194382824999699*m0r[2]*uSelf[2]*mnuSelf)-0.5*m0r[0]*uSelf[2]*mnuSelf-0.3194382824999699*m0r[2]*uOther[2]*mnuSelf-0.5*m0r[0]*uOther[2]*mnuSelf+1.414213562373095*m1r[2]*mnuSelf-0.5*uSelf[0]*m0r[2]*mnuSelf-0.5*uOther[0]*m0r[2]*mnuSelf-0.4472135954999579*m0r[1]*uSelf[1]*mnuSelf-0.4472135954999579*m0r[1]*uOther[1]*mnuSelf); 
  gkyl_mat_set(A,5,1,(-0.7857142857142857*m0r[1]*uSelf[2]*mnuSelf)-0.7857142857142857*m0r[1]*uOther[2]*mnuSelf-0.7857142857142857*uSelf[1]*m0r[2]*mnuSelf-0.7857142857142857*uOther[1]*m0r[2]*mnuSelf-0.4472135954999579*m0r[0]*uSelf[1]*mnuSelf-0.4472135954999579*m0r[0]*uOther[1]*mnuSelf+1.264911064067352*m1r[1]*mnuSelf-0.4472135954999579*uSelf[0]*m0r[1]*mnuSelf-0.4472135954999579*uOther[0]*m0r[1]*mnuSelf); 
  gkyl_mat_set(A,5,2,(-1.071428571428571*m0r[2]*uSelf[2]*mnuSelf)-0.3194382824999699*m0r[0]*uSelf[2]*mnuSelf-1.071428571428571*m0r[2]*uOther[2]*mnuSelf-0.3194382824999699*m0r[0]*uOther[2]*mnuSelf+0.9035079029052515*m1r[2]*mnuSelf-0.3194382824999699*uSelf[0]*m0r[2]*mnuSelf-0.3194382824999699*uOther[0]*m0r[2]*mnuSelf-0.7857142857142857*m0r[1]*uSelf[1]*mnuSelf-0.7857142857142857*m0r[1]*uOther[1]*mnuSelf-0.5*m0r[0]*uSelf[0]*mnuSelf-0.5*m0r[0]*uOther[0]*mnuSelf+1.414213562373095*m1r[0]*mnuSelf); 
 
  momRHS[0] += -0.7071067811865475*((m0r[2]*uSelf[2]-1.0*m0r[2]*uOther[2]+m0r[1]*uSelf[1]-1.0*m0r[1]*uOther[1]+m0r[0]*uSelf[0]-1.0*m0r[0]*uOther[0])*betaGreenep1-2.82842712474619*m1r[0])*mnuSelf; 
  momRHS[1] += -0.1414213562373095*((4.47213595499958*m0r[1]*uSelf[2]-4.47213595499958*m0r[1]*uOther[2]+(4.47213595499958*uSelf[1]-4.47213595499958*uOther[1])*m0r[2]+5.0*m0r[0]*uSelf[1]-5.0*m0r[0]*uOther[1]+(5.0*uSelf[0]-5.0*uOther[0])*m0r[1])*betaGreenep1-14.14213562373095*m1r[1])*mnuSelf; 
  momRHS[2] += -0.02020305089104421*(((22.3606797749979*m0r[2]+35.0*m0r[0])*uSelf[2]+((-22.3606797749979*m0r[2])-35.0*m0r[0])*uOther[2]+(35.0*uSelf[0]-35.0*uOther[0])*m0r[2]+31.30495168499705*m0r[1]*uSelf[1]-31.30495168499705*m0r[1]*uOther[1])*betaGreenep1-98.99494936611667*m1r[2])*mnuSelf; 
 
  double ucMSelf[3]; 
  double ucMOther[3]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (int vd=0; vd<3; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
    ucMOther[vd] = 0.0; 
  } 
  for (int vd=0; vd<1; vd++) 
  { 
    int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+2]*uSelf[a0+2]+0.7071067811865475*cMr[a0+1]*uSelf[a0+1]+0.7071067811865475*cMr[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.6324555320336759*cMr[a0+1]*uSelf[a0+2]+0.6324555320336759*uSelf[a0+1]*cMr[a0+2]+0.7071067811865475*cMr[a0]*uSelf[a0+1]+0.7071067811865475*uSelf[a0]*cMr[a0+1]; 
    ucMSelf[2] += 0.4517539514526256*cMr[a0+2]*uSelf[a0+2]+0.7071067811865475*cMr[a0]*uSelf[a0+2]+0.7071067811865475*uSelf[a0]*cMr[a0+2]+0.6324555320336759*cMr[a0+1]*uSelf[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+2]*uOther[a0+2]+0.7071067811865475*cMr[a0+1]*uOther[a0+1]+0.7071067811865475*cMr[a0]*uOther[a0]; 
    ucMOther[1] += 0.6324555320336759*cMr[a0+1]*uOther[a0+2]+0.6324555320336759*uOther[a0+1]*cMr[a0+2]+0.7071067811865475*cMr[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*cMr[a0+1]; 
    ucMOther[2] += 0.4517539514526256*cMr[a0+2]*uOther[a0+2]+0.7071067811865475*cMr[a0]*uOther[a0+2]+0.7071067811865475*uOther[a0]*cMr[a0+2]+0.6324555320336759*cMr[a0+1]*uOther[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,9,3,0.7071067811865475*ucMSelf[0]*mnuSelf+0.7071067811865475*ucMOther[0]*mnuSelf+1.414213562373095*m0r[0]*mnuSelf-1.414213562373095*cEr[0]*mnuSelf); 
  gkyl_mat_set(A,9,4,0.7071067811865475*ucMSelf[1]*mnuSelf+0.7071067811865475*ucMOther[1]*mnuSelf+1.414213562373095*m0r[1]*mnuSelf-1.414213562373095*cEr[1]*mnuSelf); 
  gkyl_mat_set(A,9,5,0.7071067811865475*ucMSelf[2]*mnuSelf+0.7071067811865475*ucMOther[2]*mnuSelf+1.414213562373095*m0r[2]*mnuSelf-1.414213562373095*cEr[2]*mnuSelf); 
  gkyl_mat_set(A,10,3,0.7071067811865475*ucMSelf[1]*mnuSelf+0.7071067811865475*ucMOther[1]*mnuSelf+1.414213562373095*m0r[1]*mnuSelf-1.414213562373095*cEr[1]*mnuSelf); 
  gkyl_mat_set(A,10,4,0.6324555320336759*ucMSelf[2]*mnuSelf+0.6324555320336759*ucMOther[2]*mnuSelf+1.264911064067352*m0r[2]*mnuSelf-1.264911064067352*cEr[2]*mnuSelf+0.7071067811865475*ucMSelf[0]*mnuSelf+0.7071067811865475*ucMOther[0]*mnuSelf+1.414213562373095*m0r[0]*mnuSelf-1.414213562373095*cEr[0]*mnuSelf); 
  gkyl_mat_set(A,10,5,0.6324555320336759*ucMSelf[1]*mnuSelf+0.6324555320336759*ucMOther[1]*mnuSelf+1.264911064067352*m0r[1]*mnuSelf-1.264911064067352*cEr[1]*mnuSelf); 
  gkyl_mat_set(A,11,3,0.7071067811865475*ucMSelf[2]*mnuSelf+0.7071067811865475*ucMOther[2]*mnuSelf+1.414213562373095*m0r[2]*mnuSelf-1.414213562373095*cEr[2]*mnuSelf); 
  gkyl_mat_set(A,11,4,0.6324555320336759*ucMSelf[1]*mnuSelf+0.6324555320336759*ucMOther[1]*mnuSelf+1.264911064067352*m0r[1]*mnuSelf-1.264911064067352*cEr[1]*mnuSelf); 
  gkyl_mat_set(A,11,5,0.4517539514526256*ucMSelf[2]*mnuSelf+0.4517539514526256*ucMOther[2]*mnuSelf+0.9035079029052515*m0r[2]*mnuSelf-0.9035079029052515*cEr[2]*mnuSelf+0.7071067811865475*ucMSelf[0]*mnuSelf+0.7071067811865475*ucMOther[0]*mnuSelf+1.414213562373095*m0r[0]*mnuSelf-1.414213562373095*cEr[0]*mnuSelf); 
 
  double uM1Self[3]; 
  double uM1Other[3]; 
  double uSumSq[3]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (int vd=0; vd<3; vd++) 
  { 
    uM1Self[vd] = 0.0; 
    uM1Other[vd] = 0.0; 
    uSumSq[vd] = 0.0; 
  } 
  for (int vd=0; vd<1; vd++) 
  { 
    int a0 = 3*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.7071067811865475*m1r[2]*uSelf[a0+2]+0.7071067811865475*m1r[1]*uSelf[a0+1]+0.7071067811865475*m1r[0]*uSelf[a0]; 
    uM1Self[1] += 0.6324555320336759*m1r[1]*uSelf[a0+2]+0.6324555320336759*m1r[2]*uSelf[a0+1]+0.7071067811865475*m1r[0]*uSelf[a0+1]+0.7071067811865475*m1r[1]*uSelf[a0]; 
    uM1Self[2] += 0.4517539514526256*m1r[2]*uSelf[a0+2]+0.7071067811865475*m1r[0]*uSelf[a0+2]+0.6324555320336759*m1r[1]*uSelf[a0+1]+0.7071067811865475*m1r[2]*uSelf[a0]; 
    uM1Other[0] += 0.7071067811865475*m1r[2]*uOther[a0+2]+0.7071067811865475*m1r[1]*uOther[a0+1]+0.7071067811865475*m1r[0]*uOther[a0]; 
    uM1Other[1] += 0.6324555320336759*m1r[1]*uOther[a0+2]+0.6324555320336759*m1r[2]*uOther[a0+1]+0.7071067811865475*m1r[0]*uOther[a0+1]+0.7071067811865475*m1r[1]*uOther[a0]; 
    uM1Other[2] += 0.4517539514526256*m1r[2]*uOther[a0+2]+0.7071067811865475*m1r[0]*uOther[a0+2]+0.6324555320336759*m1r[1]*uOther[a0+1]+0.7071067811865475*m1r[2]*uOther[a0]; 
  const double uSelf0R2 = pow(uSelf[a0],2);
  const double uSelf1R2 = pow(uSelf[a0+1],2);
  const double uSelf2R2 = pow(uSelf[a0+2],2);
  const double uOther0R2 = pow(uOther[a0],2);
  const double uOther1R2 = pow(uOther[a0+1],2);
  const double uOther2R2 = pow(uOther[a0+2],2);

  uSumSq[0] += 0.7071067811865475*uSelf2R2-1.414213562373095*uOther[a0+2]*uSelf[a0+2]+0.7071067811865475*uOther2R2+0.7071067811865475*uSelf1R2-1.414213562373095*uOther[a0+1]*uSelf[a0+1]+0.7071067811865475*uOther1R2+0.7071067811865475*uSelf0R2-1.414213562373095*uOther[a0]*uSelf[a0]+0.7071067811865475*uOther0R2; 
  uSumSq[1] += 1.264911064067352*uSelf[a0+1]*uSelf[a0+2]-1.264911064067352*uOther[a0+1]*uSelf[a0+2]-1.264911064067352*uSelf[a0+1]*uOther[a0+2]+1.264911064067352*uOther[a0+1]*uOther[a0+2]+1.414213562373095*uSelf[a0]*uSelf[a0+1]-1.414213562373095*uOther[a0]*uSelf[a0+1]-1.414213562373095*uSelf[a0]*uOther[a0+1]+1.414213562373095*uOther[a0]*uOther[a0+1]; 
  uSumSq[2] += 0.4517539514526256*uSelf2R2-0.9035079029052515*uOther[a0+2]*uSelf[a0+2]+1.414213562373095*uSelf[a0]*uSelf[a0+2]-1.414213562373095*uOther[a0]*uSelf[a0+2]+0.4517539514526256*uOther2R2-1.414213562373095*uSelf[a0]*uOther[a0+2]+1.414213562373095*uOther[a0]*uOther[a0+2]+0.6324555320336759*uSelf1R2-1.264911064067352*uOther[a0+1]*uSelf[a0+1]+0.6324555320336759*uOther1R2; 
  } 
 
  double enRHS[3]; 
  enRHS[0] = (-(2.0*m0r[2]*vtSqSelf[2]*betaGreenep1*mSelf*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther))-(1.0*m0r[2]*uSumSq[2]*betaGreenep1*mSelf*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)-(2.0*m0r[1]*vtSqSelf[1]*betaGreenep1*mSelf*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)-(1.0*m0r[1]*uSumSq[1]*betaGreenep1*mSelf*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)-(2.0*m0r[0]*vtSqSelf[0]*betaGreenep1*mSelf*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)-(1.0*m0r[0]*uSumSq[0]*betaGreenep1*mSelf*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)+(2.0*m0r[2]*vtSqOther[2]*betaGreenep1*mOther*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)+(m0r[2]*uSumSq[2]*betaGreenep1*mOther*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)+(2.0*m0r[1]*vtSqOther[1]*betaGreenep1*mOther*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)+(m0r[1]*uSumSq[1]*betaGreenep1*mOther*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)+(2.0*m0r[0]*vtSqOther[0]*betaGreenep1*mOther*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)+(m0r[0]*uSumSq[0]*betaGreenep1*mOther*mnuSelf)/(2.82842712474619*mSelf+2.82842712474619*mOther)-1.0*uM1Self[0]*mnuSelf-1.0*uM1Other[0]*mnuSelf+2.0*m2r[0]*mnuSelf; 
  enRHS[1] = (-(8.94427190999916*m0r[1]*vtSqSelf[2]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther))-(4.47213595499958*m0r[1]*uSumSq[2]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-(8.94427190999916*vtSqSelf[1]*m0r[2]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-(4.47213595499958*uSumSq[1]*m0r[2]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-(10.0*m0r[0]*vtSqSelf[1]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-(5.0*m0r[0]*uSumSq[1]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-(10.0*vtSqSelf[0]*m0r[1]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-(5.0*uSumSq[0]*m0r[1]*betaGreenep1*mSelf*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(8.94427190999916*m0r[1]*vtSqOther[2]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(4.47213595499958*m0r[1]*uSumSq[2]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(8.94427190999916*vtSqOther[1]*m0r[2]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(4.47213595499958*uSumSq[1]*m0r[2]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(10.0*m0r[0]*vtSqOther[1]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(5.0*m0r[0]*uSumSq[1]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(10.0*vtSqOther[0]*m0r[1]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)+(5.0*uSumSq[0]*m0r[1]*betaGreenep1*mOther*mnuSelf)/(14.14213562373095*mSelf+14.14213562373095*mOther)-1.0*uM1Self[1]*mnuSelf-1.0*uM1Other[1]*mnuSelf+2.0*m2r[1]*mnuSelf; 
  enRHS[2] = (-(44.7213595499958*m0r[2]*vtSqSelf[2]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther))-(70.0*m0r[0]*vtSqSelf[2]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-(22.3606797749979*m0r[2]*uSumSq[2]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-(35.0*m0r[0]*uSumSq[2]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-(70.0*vtSqSelf[0]*m0r[2]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-(35.0*uSumSq[0]*m0r[2]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-(62.60990336999411*m0r[1]*vtSqSelf[1]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-(31.30495168499705*m0r[1]*uSumSq[1]*betaGreenep1*mSelf*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(44.7213595499958*m0r[2]*vtSqOther[2]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(70.0*m0r[0]*vtSqOther[2]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(22.3606797749979*m0r[2]*uSumSq[2]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(35.0*m0r[0]*uSumSq[2]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(70.0*vtSqOther[0]*m0r[2]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(35.0*uSumSq[0]*m0r[2]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(62.60990336999411*m0r[1]*vtSqOther[1]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)+(31.30495168499705*m0r[1]*uSumSq[1]*betaGreenep1*mOther*mnuSelf)/(98.99494936611667*mSelf+98.99494936611667*mOther)-1.0*uM1Self[2]*mnuSelf-1.0*uM1Other[2]*mnuSelf+2.0*m2r[2]*mnuSelf; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,enRHS[0]); 
  gkyl_mat_set(rhs,4,0,enRHS[1]); 
  gkyl_mat_set(rhs,5,0,enRHS[2]); 
} 
 
