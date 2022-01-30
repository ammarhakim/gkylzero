#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*(2.23606797749979*m0[2]-1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m0[2]+1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m2[2]-1.732050807568877*m2[1]+m2[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m2[2]+1.732050807568877*m2[1]+m2[0]) < 0) cellAvg = true; 
 
  double m0r[3] = {0.0}; 
  double m1r[3] = {0.0}; 
  double cMr[3] = {0.0}; 
  double cEr[3] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    gkyl_mat_set(rhs,0,0,m1[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,3,0,m2[0]); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    gkyl_mat_set(rhs,0,0,m1[0]); 
    gkyl_mat_set(rhs,1,0,m1[1]); 
    gkyl_mat_set(rhs,2,0,m1[2]); 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    gkyl_mat_set(rhs,3,0,m2[0]); 
    gkyl_mat_set(rhs,4,0,m2[1]); 
    gkyl_mat_set(rhs,5,0,m2[2]); 
  } 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,1,2,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,2,1,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,2,0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to uX .......... // 
  gkyl_mat_set(A,0,3,-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,0,4,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,0,5,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,1,3,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,1,4,(-0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,5,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,3,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,2,4,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,5,(-0.4517539514526256*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  gkyl_mat_set(A,3,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,3,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,4,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,4,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,4,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,5,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,2,0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,3,3,0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,3,4,0.7071067811865475*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,5,0.7071067811865475*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,4,3,0.7071067811865475*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,4,4,0.6324555320336759*m0r[2]-0.6324555320336759*cEr[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,4,5,0.6324555320336759*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,3,0.7071067811865475*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,5,4,0.6324555320336759*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,5,0.4517539514526256*m0r[2]-0.4517539514526256*cEr[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
