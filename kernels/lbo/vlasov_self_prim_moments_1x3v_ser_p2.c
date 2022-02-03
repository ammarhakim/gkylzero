#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE) 
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
  double m1r[9] = {0.0}; 
  double cMr[9] = {0.0}; 
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
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    gkyl_mat_set(rhs,0,0,m1[6]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[6] = cM[6]; 
    cMr[7] = 0.0; 
    cMr[8] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,9,0,m2[0]); 
    gkyl_mat_set(rhs,10,0,0.0); 
    gkyl_mat_set(rhs,11,0,0.0); 
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
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    m1r[8] = m1[8]; 
    gkyl_mat_set(rhs,0,0,m1[0]); 
    gkyl_mat_set(rhs,1,0,m1[1]); 
    gkyl_mat_set(rhs,2,0,m1[2]); 
    gkyl_mat_set(rhs,3,0,m1[3]); 
    gkyl_mat_set(rhs,4,0,m1[4]); 
    gkyl_mat_set(rhs,5,0,m1[5]); 
    gkyl_mat_set(rhs,6,0,m1[6]); 
    gkyl_mat_set(rhs,7,0,m1[7]); 
    gkyl_mat_set(rhs,8,0,m1[8]); 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cMr[4] = cM[4]; 
    cMr[5] = cM[5]; 
    cMr[6] = cM[6]; 
    cMr[7] = cM[7]; 
    cMr[8] = cM[8]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    gkyl_mat_set(rhs,9,0,m2[0]); 
    gkyl_mat_set(rhs,10,0,m2[1]); 
    gkyl_mat_set(rhs,11,0,m2[2]); 
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
  gkyl_mat_set(A,0,9,-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,0,10,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,0,11,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,1,9,-0.7071067811865475*cMr[1]); 
  gkyl_mat_set(A,1,10,(-0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,11,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,9,-0.7071067811865475*cMr[2]); 
  gkyl_mat_set(A,2,10,-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,11,(-0.4517539514526256*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  gkyl_mat_set(A,9,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,9,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,9,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,10,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,10,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,10,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,11,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,11,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,11,2,0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]); 
 
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
  gkyl_mat_set(A,3,9,-0.7071067811865475*cMr[3]); 
  gkyl_mat_set(A,3,10,-0.7071067811865475*cMr[4]); 
  gkyl_mat_set(A,3,11,-0.7071067811865475*cMr[5]); 
  gkyl_mat_set(A,4,9,-0.7071067811865475*cMr[4]); 
  gkyl_mat_set(A,4,10,(-0.6324555320336759*cMr[5])-0.7071067811865475*cMr[3]); 
  gkyl_mat_set(A,4,11,-0.6324555320336759*cMr[4]); 
  gkyl_mat_set(A,5,9,-0.7071067811865475*cMr[5]); 
  gkyl_mat_set(A,5,10,-0.6324555320336759*cMr[4]); 
  gkyl_mat_set(A,5,11,(-0.4517539514526256*cMr[5])-0.7071067811865475*cMr[3]); 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  gkyl_mat_set(A,9,3,0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,9,4,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,9,5,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,10,3,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,10,4,0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,10,5,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,11,3,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,11,4,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,11,5,0.4517539514526256*m1r[5]+0.7071067811865475*m1r[3]); 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  gkyl_mat_set(A,6,6,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,6,7,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,6,8,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,7,6,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,7,7,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,7,8,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,8,6,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,8,7,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,8,8,0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to uZ .......... // 
  gkyl_mat_set(A,6,9,-0.7071067811865475*cMr[6]); 
  gkyl_mat_set(A,6,10,-0.7071067811865475*cMr[7]); 
  gkyl_mat_set(A,6,11,-0.7071067811865475*cMr[8]); 
  gkyl_mat_set(A,7,9,-0.7071067811865475*cMr[7]); 
  gkyl_mat_set(A,7,10,(-0.6324555320336759*cMr[8])-0.7071067811865475*cMr[6]); 
  gkyl_mat_set(A,7,11,-0.6324555320336759*cMr[7]); 
  gkyl_mat_set(A,8,9,-0.7071067811865475*cMr[8]); 
  gkyl_mat_set(A,8,10,-0.6324555320336759*cMr[7]); 
  gkyl_mat_set(A,8,11,(-0.4517539514526256*cMr[8])-0.7071067811865475*cMr[6]); 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  gkyl_mat_set(A,9,6,0.7071067811865475*m1r[6]); 
  gkyl_mat_set(A,9,7,0.7071067811865475*m1r[7]); 
  gkyl_mat_set(A,9,8,0.7071067811865475*m1r[8]); 
  gkyl_mat_set(A,10,6,0.7071067811865475*m1r[7]); 
  gkyl_mat_set(A,10,7,0.6324555320336759*m1r[8]+0.7071067811865475*m1r[6]); 
  gkyl_mat_set(A,10,8,0.6324555320336759*m1r[7]); 
  gkyl_mat_set(A,11,6,0.7071067811865475*m1r[8]); 
  gkyl_mat_set(A,11,7,0.6324555320336759*m1r[7]); 
  gkyl_mat_set(A,11,8,0.4517539514526256*m1r[8]+0.7071067811865475*m1r[6]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,9,9,2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,9,10,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,9,11,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,10,9,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,10,10,1.897366596101028*m0r[2]-0.6324555320336759*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,10,11,1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,11,9,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,11,10,1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,11,11,1.355261854357877*m0r[2]-0.4517539514526256*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
