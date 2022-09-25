#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_pkpm_self_prim_moments_1x1v_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *p_ij, const double *boundary_corrections) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, vpar^2, vpar^3, only need Zeroth). 
  // p_ij:                 Fluid pressure tensor (trace of pressure tensor used to obtain vt^2). 
  // boundary_corrections: boundary corrections to vtSq. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(3.741657386773942*moms[3]-3.16227766016838*moms[2]+2.449489742783178*moms[1]-1.414213562373095*moms[0]) < 0) cellAvg = true; 
  if (0.5*(3.741657386773942*moms[3]+3.16227766016838*moms[2]+2.449489742783178*moms[1]+1.414213562373095*moms[0]) < 0) cellAvg = true; 
 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pzz = &p_ij[20]; 
  double m0r[4] = {0.0}; 
  double cEr[4] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    cEr[0] = boundary_corrections[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    gkyl_mat_set(rhs,0,0,Pzz[0]+Pyy[0]+Pxx[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m0r[3] = moms[3]; 
    cEr[0] = boundary_corrections[0]; 
    cEr[1] = boundary_corrections[1]; 
    cEr[2] = boundary_corrections[2]; 
    cEr[3] = boundary_corrections[3]; 
    gkyl_mat_set(rhs,0,0,Pzz[0]+Pyy[0]+Pxx[0]); 
    gkyl_mat_set(rhs,1,0,Pzz[1]+Pyy[1]+Pxx[1]); 
    gkyl_mat_set(rhs,2,0,Pzz[2]+Pyy[2]+Pxx[2]); 
    gkyl_mat_set(rhs,3,0,Pzz[3]+Pyy[3]+Pxx[3]); 
  } 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,0,0,2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,0,1,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,0,2,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,0,3,2.121320343559642*m0r[3]-0.7071067811865475*cEr[3]); 
  gkyl_mat_set(A,1,0,2.121320343559642*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,1,1,1.897366596101028*m0r[2]-0.6324555320336759*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,1,2,1.863177010224355*m0r[3]-0.6210590034081186*cEr[3]+1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,1,3,1.863177010224355*m0r[2]-0.6210590034081186*cEr[2]); 
  gkyl_mat_set(A,2,0,2.121320343559642*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,2,1,1.863177010224355*m0r[3]-0.6210590034081186*cEr[3]+1.897366596101028*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,2,2,1.355261854357877*m0r[2]-0.4517539514526256*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,2,3,1.264911064067352*m0r[3]-0.421637021355784*cEr[3]+1.863177010224355*m0r[1]-0.6210590034081186*cEr[1]); 
  gkyl_mat_set(A,3,0,2.121320343559642*m0r[3]-0.7071067811865475*cEr[3]); 
  gkyl_mat_set(A,3,1,1.863177010224355*m0r[2]-0.6210590034081186*cEr[2]); 
  gkyl_mat_set(A,3,2,1.264911064067352*m0r[3]-0.421637021355784*cEr[3]+1.863177010224355*m0r[1]-0.6210590034081186*cEr[1]); 
  gkyl_mat_set(A,3,3,1.264911064067352*m0r[2]-0.421637021355784*cEr[2]+2.121320343559642*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
