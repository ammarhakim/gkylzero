#include <gkyl_prim_lbo_vlasov_pkpm_kernels.h> 
 
GKYL_CU_DH void vlasov_pkpm_self_prim_moments_3x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *vlasov_pkpm_moms, const double *boundary_corrections) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // vlasov_pkpm_moms:     [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model. 
  // boundary_corrections: boundary corrections to vtSq. 
 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *p_parallel = &vlasov_pkpm_moms[8]; 
  const double *p_perp = &vlasov_pkpm_moms[16]; 
  // If a corner value is below zero, use cell average.
  bool cellAvg = false;
  if (-0.25*(7.348469228349534*rho[7]-4.242640687119286*(rho[6]+rho[5]+rho[4])+2.449489742783178*(rho[3]+rho[2]+rho[1])-1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_parallel[7]-4.242640687119286*(p_parallel[6]+p_parallel[5]+p_parallel[4])+2.449489742783178*(p_parallel[3]+p_parallel[2]+p_parallel[1])-1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_perp[7]-4.242640687119286*(p_perp[6]+p_perp[5]+p_perp[4])+2.449489742783178*(p_perp[3]+p_perp[2]+p_perp[1])-1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*rho[7]+4.242640687119286*rho[6]-4.242640687119286*(rho[5]+rho[4])-2.449489742783178*(rho[3]+rho[2])+2.449489742783178*rho[1]+1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_parallel[7]+4.242640687119286*p_parallel[6]-4.242640687119286*(p_parallel[5]+p_parallel[4])-2.449489742783178*(p_parallel[3]+p_parallel[2])+2.449489742783178*p_parallel[1]+1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_perp[7]+4.242640687119286*p_perp[6]-4.242640687119286*(p_perp[5]+p_perp[4])-2.449489742783178*(p_perp[3]+p_perp[2])+2.449489742783178*p_perp[1]+1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*rho[7]-4.242640687119286*rho[6]+4.242640687119286*rho[5]-4.242640687119286*rho[4]-2.449489742783178*rho[3]+2.449489742783178*rho[2]-2.449489742783178*rho[1]+1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_parallel[7]-4.242640687119286*p_parallel[6]+4.242640687119286*p_parallel[5]-4.242640687119286*p_parallel[4]-2.449489742783178*p_parallel[3]+2.449489742783178*p_parallel[2]-2.449489742783178*p_parallel[1]+1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_perp[7]-4.242640687119286*p_perp[6]+4.242640687119286*p_perp[5]-4.242640687119286*p_perp[4]-2.449489742783178*p_perp[3]+2.449489742783178*p_perp[2]-2.449489742783178*p_perp[1]+1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*rho[7]+4.242640687119286*(rho[6]+rho[5])-4.242640687119286*rho[4]+2.449489742783178*rho[3]-2.449489742783178*(rho[2]+rho[1])-1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_parallel[7]+4.242640687119286*(p_parallel[6]+p_parallel[5])-4.242640687119286*p_parallel[4]+2.449489742783178*p_parallel[3]-2.449489742783178*(p_parallel[2]+p_parallel[1])-1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_perp[7]+4.242640687119286*(p_perp[6]+p_perp[5])-4.242640687119286*p_perp[4]+2.449489742783178*p_perp[3]-2.449489742783178*(p_perp[2]+p_perp[1])-1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*rho[7]-4.242640687119286*(rho[6]+rho[5])+4.242640687119286*rho[4]+2.449489742783178*rho[3]-2.449489742783178*(rho[2]+rho[1])+1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_parallel[7]-4.242640687119286*(p_parallel[6]+p_parallel[5])+4.242640687119286*p_parallel[4]+2.449489742783178*p_parallel[3]-2.449489742783178*(p_parallel[2]+p_parallel[1])+1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_perp[7]-4.242640687119286*(p_perp[6]+p_perp[5])+4.242640687119286*p_perp[4]+2.449489742783178*p_perp[3]-2.449489742783178*(p_perp[2]+p_perp[1])+1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*rho[7]+4.242640687119286*rho[6]-4.242640687119286*rho[5]+4.242640687119286*rho[4]-2.449489742783178*rho[3]+2.449489742783178*rho[2]-2.449489742783178*rho[1]-1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_parallel[7]+4.242640687119286*p_parallel[6]-4.242640687119286*p_parallel[5]+4.242640687119286*p_parallel[4]-2.449489742783178*p_parallel[3]+2.449489742783178*p_parallel[2]-2.449489742783178*p_parallel[1]-1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_perp[7]+4.242640687119286*p_perp[6]-4.242640687119286*p_perp[5]+4.242640687119286*p_perp[4]-2.449489742783178*p_perp[3]+2.449489742783178*p_perp[2]-2.449489742783178*p_perp[1]-1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*rho[7]-4.242640687119286*rho[6]+4.242640687119286*(rho[5]+rho[4])-2.449489742783178*(rho[3]+rho[2])+2.449489742783178*rho[1]-1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_parallel[7]-4.242640687119286*p_parallel[6]+4.242640687119286*(p_parallel[5]+p_parallel[4])-2.449489742783178*(p_parallel[3]+p_parallel[2])+2.449489742783178*p_parallel[1]-1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (-0.25*(7.348469228349534*p_perp[7]-4.242640687119286*p_perp[6]+4.242640687119286*(p_perp[5]+p_perp[4])-2.449489742783178*(p_perp[3]+p_perp[2])+2.449489742783178*p_perp[1]-1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*rho[7]+4.242640687119286*(rho[6]+rho[5]+rho[4])+2.449489742783178*(rho[3]+rho[2]+rho[1])+1.414213562373095*rho[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_parallel[7]+4.242640687119286*(p_parallel[6]+p_parallel[5]+p_parallel[4])+2.449489742783178*(p_parallel[3]+p_parallel[2]+p_parallel[1])+1.414213562373095*p_parallel[0]) < 0) cellAvg = true; 
  if (0.25*(7.348469228349534*p_perp[7]+4.242640687119286*(p_perp[6]+p_perp[5]+p_perp[4])+2.449489742783178*(p_perp[3]+p_perp[2]+p_perp[1])+1.414213562373095*p_perp[0]) < 0) cellAvg = true; 
 
  double m0r[8] = {0.0}; 
  double cEr[8] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = rho[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    cEr[0] = boundary_corrections[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    cEr[4] = 0.0; 
    cEr[5] = 0.0; 
    cEr[6] = 0.0; 
    cEr[7] = 0.0; 
    gkyl_mat_set(rhs,0,0,2.0*p_perp[0]+p_parallel[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
  } else { 
    m0r[0] = rho[0]; 
    m0r[1] = rho[1]; 
    m0r[2] = rho[2]; 
    m0r[3] = rho[3]; 
    m0r[4] = rho[4]; 
    m0r[5] = rho[5]; 
    m0r[6] = rho[6]; 
    m0r[7] = rho[7]; 
    cEr[0] = boundary_corrections[0]; 
    cEr[1] = boundary_corrections[1]; 
    cEr[2] = boundary_corrections[2]; 
    cEr[3] = boundary_corrections[3]; 
    cEr[4] = boundary_corrections[4]; 
    cEr[5] = boundary_corrections[5]; 
    cEr[6] = boundary_corrections[6]; 
    cEr[7] = boundary_corrections[7]; 
    gkyl_mat_set(rhs,0,0,2.0*p_perp[0]+p_parallel[0]); 
    gkyl_mat_set(rhs,1,0,2.0*p_perp[1]+p_parallel[1]); 
    gkyl_mat_set(rhs,2,0,2.0*p_perp[2]+p_parallel[2]); 
    gkyl_mat_set(rhs,3,0,2.0*p_perp[3]+p_parallel[3]); 
    gkyl_mat_set(rhs,4,0,2.0*p_perp[4]+p_parallel[4]); 
    gkyl_mat_set(rhs,5,0,2.0*p_perp[5]+p_parallel[5]); 
    gkyl_mat_set(rhs,6,0,2.0*p_perp[6]+p_parallel[6]); 
    gkyl_mat_set(rhs,7,0,2.0*p_perp[7]+p_parallel[7]); 
  } 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,0,0,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,0,1,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,0,2,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,0,3,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,0,4,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,0,5,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,0,6,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,0,7,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,1,0,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,1,1,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,1,2,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,1,3,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,1,4,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,1,5,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,1,6,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,1,7,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,2,0,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,2,1,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,2,2,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,2,3,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,2,4,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,2,5,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,2,6,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,2,7,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,3,0,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,3,1,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,3,2,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,3,3,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,3,4,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,3,5,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,3,6,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,3,7,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,4,0,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,4,1,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,4,2,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,4,3,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,4,4,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,4,5,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,4,6,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,4,7,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,5,0,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,5,1,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,5,2,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,5,3,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,5,4,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,5,5,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,5,6,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,5,7,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,6,0,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,6,1,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,6,2,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,6,3,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,6,4,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,6,5,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,6,6,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
  gkyl_mat_set(A,6,7,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,7,0,1.060660171779821*m0r[7]-0.3535533905932737*cEr[7]); 
  gkyl_mat_set(A,7,1,1.060660171779821*m0r[6]-0.3535533905932737*cEr[6]); 
  gkyl_mat_set(A,7,2,1.060660171779821*m0r[5]-0.3535533905932737*cEr[5]); 
  gkyl_mat_set(A,7,3,1.060660171779821*m0r[4]-0.3535533905932737*cEr[4]); 
  gkyl_mat_set(A,7,4,1.060660171779821*m0r[3]-0.3535533905932737*cEr[3]); 
  gkyl_mat_set(A,7,5,1.060660171779821*m0r[2]-0.3535533905932737*cEr[2]); 
  gkyl_mat_set(A,7,6,1.060660171779821*m0r[1]-0.3535533905932737*cEr[1]); 
  gkyl_mat_set(A,7,7,1.060660171779821*m0r[0]-0.3535533905932737*cEr[0]); 
 
} 
 
