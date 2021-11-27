#include <gkyl_prim_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq) 
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
  double m1r[6] = {0.0}; 
  double cMr[6] = {0.0}; 
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
    m1r[3] = m1[3]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    gkyl_mat_set(rhs,0,0,m1[3]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[3] = cM[3]; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,6,0,m2[0]); 
    gkyl_mat_set(rhs,7,0,0.0); 
    gkyl_mat_set(rhs,8,0,0.0); 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    gkyl_mat_set(rhs,0,0,m1[0]); 
    gkyl_mat_set(rhs,1,0,m1[1]); 
    gkyl_mat_set(rhs,2,0,m1[2]); 
    gkyl_mat_set(rhs,3,0,m1[3]); 
    gkyl_mat_set(rhs,4,0,m1[4]); 
    gkyl_mat_set(rhs,5,0,m1[5]); 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cMr[4] = cM[4]; 
    cMr[5] = cM[5]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    gkyl_mat_set(rhs,6,0,m2[0]); 
    gkyl_mat_set(rhs,7,0,m2[1]); 
    gkyl_mat_set(rhs,8,0,m2[2]); 
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
  gkyl_mat_set(A,0,6,-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,0,7,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,0,8,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,1,6,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,1,7,(-0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,8,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,6,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,2,7,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,8,(-0.4517539514526256*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  gkyl_mat_set(A,6,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,6,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,6,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,7,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,7,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,7,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,8,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,8,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,8,2,0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]); 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  gkyl_mat_set(A,3,3,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,3,4,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,3,5,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,4,3,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,4,4,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,4,5,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,5,3,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,5,4,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,5,5,0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to uY .......... // 
  gkyl_mat_set(A,3,6,-0.7071067811865475*cMr[3]); 
  gkyl_mat_set(A,3,7,-0.7071067811865475*cMr[4]); 
  gkyl_mat_set(A,3,8,-0.7071067811865475*cMr[5]); 
  gkyl_mat_set(A,4,6,-0.7071067811865475*cMr[4]); 
  gkyl_mat_set(A,4,7,(-0.6324555320336759*cMr[5])-0.7071067811865475*cMr[3]); 
  gkyl_mat_set(A,4,8,-0.6324555320336759*cMr[4]); 
  gkyl_mat_set(A,5,6,-0.7071067811865475*cMr[5]); 
  gkyl_mat_set(A,5,7,-0.6324555320336759*cMr[4]); 
  gkyl_mat_set(A,5,8,(-0.4517539514526256*cMr[5])-0.7071067811865475*cMr[3]); 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  gkyl_mat_set(A,6,3,0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,6,4,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,6,5,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,7,3,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,7,4,0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,7,5,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,8,3,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,8,4,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,8,5,0.4517539514526256*m1r[5]+0.7071067811865475*m1r[3]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,6,6,2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,6,7,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,6,8,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,7,6,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,7,7,1.897366596101028*m0r[2]-0.6324555320336759*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,7,8,1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,8,6,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,8,7,1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,8,8,1.355261854357877*m0r[2]-0.4517539514526256*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
 
  long ipiv[9] = {0.0}; 
  gkyl_mat_linsolve_lu(A,rhs,ipiv); 
  for(size_t i=0; i<9; i++) { 
    if (i<3) { 
      vtSq[i] = gkyl_mat_get(rhs,i+6,0); 
    } 
    u[i] = gkyl_mat_get(rhs,i,0); 
  } 
 
} 
 
