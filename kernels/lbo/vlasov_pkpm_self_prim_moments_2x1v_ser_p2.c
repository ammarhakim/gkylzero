#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_pkpm_self_prim_moments_2x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *p_ij, const double *boundary_corrections) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, vpar^2, vpar^3, only need Zeroth). 
  // p_ij:                 Fluid pressure tensor (trace of pressure tensor used to obtain vt^2). 
  // boundary_corrections: boundary corrections to vtSq. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(3.872983346207417*(moms[7]+moms[6])-2.23606797749979*(moms[5]+moms[4])-3.0*moms[3]+1.732050807568877*(moms[2]+moms[1])-1.0*moms[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*moms[7]-3.872983346207417*moms[6]+2.23606797749979*(moms[5]+moms[4])-3.0*moms[3]-1.732050807568877*moms[2]+1.732050807568877*moms[1]+moms[0]) < 0) cellAvg = true; 
  if (-0.5*(3.872983346207417*moms[7]-3.872983346207417*moms[6]-2.23606797749979*(moms[5]+moms[4])+3.0*moms[3]-1.732050807568877*moms[2]+1.732050807568877*moms[1]-1.0*moms[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*(moms[7]+moms[6])+2.23606797749979*(moms[5]+moms[4])+3.0*moms[3]+1.732050807568877*(moms[2]+moms[1])+moms[0]) < 0) cellAvg = true; 
 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[24]; 
  const double *Pzz = &p_ij[40]; 
  double m0r[8] = {0.0}; 
  double cEr[8] = {0.0}; 
  if (cellAvg) { 
    m0r[0] = moms[0]; 
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
    gkyl_mat_set(rhs,0,0,Pzz[0]+Pyy[0]+Pxx[0]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    gkyl_mat_set(rhs,4,0,0.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m0r[3] = moms[3]; 
    m0r[4] = moms[4]; 
    m0r[5] = moms[5]; 
    m0r[6] = moms[6]; 
    m0r[7] = moms[7]; 
    cEr[0] = boundary_corrections[0]; 
    cEr[1] = boundary_corrections[1]; 
    cEr[2] = boundary_corrections[2]; 
    cEr[3] = boundary_corrections[3]; 
    cEr[4] = boundary_corrections[4]; 
    cEr[5] = boundary_corrections[5]; 
    cEr[6] = boundary_corrections[6]; 
    cEr[7] = boundary_corrections[7]; 
    gkyl_mat_set(rhs,0,0,Pzz[0]+Pyy[0]+Pxx[0]); 
    gkyl_mat_set(rhs,1,0,Pzz[1]+Pyy[1]+Pxx[1]); 
    gkyl_mat_set(rhs,2,0,Pzz[2]+Pyy[2]+Pxx[2]); 
    gkyl_mat_set(rhs,3,0,Pzz[3]+Pyy[3]+Pxx[3]); 
    gkyl_mat_set(rhs,4,0,Pzz[4]+Pyy[4]+Pxx[4]); 
    gkyl_mat_set(rhs,5,0,Pzz[5]+Pyy[5]+Pxx[5]); 
    gkyl_mat_set(rhs,6,0,Pzz[6]+Pyy[6]+Pxx[6]); 
    gkyl_mat_set(rhs,7,0,Pzz[7]+Pyy[7]+Pxx[7]); 
  } 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,0,0,1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,0,1,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,0,2,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,0,3,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,0,4,1.5*m0r[4]-0.5*cEr[4]); 
  gkyl_mat_set(A,0,5,1.5*m0r[5]-0.5*cEr[5]); 
  gkyl_mat_set(A,0,6,1.5*m0r[6]-0.5*cEr[6]); 
  gkyl_mat_set(A,0,7,1.5*m0r[7]-0.5*cEr[7]); 
  gkyl_mat_set(A,1,0,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,1,1,1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,1,2,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,1,3,1.341640786499874*m0r[6]-0.447213595499958*cEr[6]+1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,1,4,1.341640786499874*m0r[1]-0.4472135954999579*cEr[1]); 
  gkyl_mat_set(A,1,5,1.5*m0r[7]-0.5000000000000001*cEr[7]); 
  gkyl_mat_set(A,1,6,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,1,7,1.5*m0r[5]-0.5000000000000001*cEr[5]); 
  gkyl_mat_set(A,2,0,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,2,1,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,2,2,1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,2,3,1.341640786499874*m0r[7]-0.447213595499958*cEr[7]+1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,2,4,1.5*m0r[6]-0.5000000000000001*cEr[6]); 
  gkyl_mat_set(A,2,5,1.341640786499874*m0r[2]-0.4472135954999579*cEr[2]); 
  gkyl_mat_set(A,2,6,1.5*m0r[4]-0.5000000000000001*cEr[4]); 
  gkyl_mat_set(A,2,7,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,3,0,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,3,1,1.341640786499874*m0r[6]-0.447213595499958*cEr[6]+1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,3,2,1.341640786499874*m0r[7]-0.447213595499958*cEr[7]+1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,3,3,1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,3,4,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,3,5,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,3,6,1.2*m0r[7]-0.4*cEr[7]+1.341640786499874*m0r[1]-0.447213595499958*cEr[1]); 
  gkyl_mat_set(A,3,7,1.2*m0r[6]-0.4*cEr[6]+1.341640786499874*m0r[2]-0.447213595499958*cEr[2]); 
  gkyl_mat_set(A,4,0,1.5*m0r[4]-0.5*cEr[4]); 
  gkyl_mat_set(A,4,1,1.341640786499874*m0r[1]-0.4472135954999579*cEr[1]); 
  gkyl_mat_set(A,4,2,1.5*m0r[6]-0.5000000000000001*cEr[6]); 
  gkyl_mat_set(A,4,3,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,4,4,0.9583148474999099*m0r[4]-0.31943828249997*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,4,6,0.9583148474999099*m0r[6]-0.31943828249997*cEr[6]+1.5*m0r[2]-0.5000000000000001*cEr[2]); 
  gkyl_mat_set(A,4,7,1.341640786499874*m0r[7]-0.4472135954999579*cEr[7]); 
  gkyl_mat_set(A,5,0,1.5*m0r[5]-0.5*cEr[5]); 
  gkyl_mat_set(A,5,1,1.5*m0r[7]-0.5000000000000001*cEr[7]); 
  gkyl_mat_set(A,5,2,1.341640786499874*m0r[2]-0.4472135954999579*cEr[2]); 
  gkyl_mat_set(A,5,3,1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]); 
  gkyl_mat_set(A,5,5,0.9583148474999099*m0r[5]-0.31943828249997*cEr[5]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,5,6,1.341640786499874*m0r[6]-0.4472135954999579*cEr[6]); 
  gkyl_mat_set(A,5,7,0.9583148474999099*m0r[7]-0.31943828249997*cEr[7]+1.5*m0r[1]-0.5000000000000001*cEr[1]); 
  gkyl_mat_set(A,6,0,1.5*m0r[6]-0.5*cEr[6]); 
  gkyl_mat_set(A,6,1,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,6,2,1.5*m0r[4]-0.5000000000000001*cEr[4]); 
  gkyl_mat_set(A,6,3,1.2*m0r[7]-0.4*cEr[7]+1.341640786499874*m0r[1]-0.447213595499958*cEr[1]); 
  gkyl_mat_set(A,6,4,0.9583148474999099*m0r[6]-0.31943828249997*cEr[6]+1.5*m0r[2]-0.5000000000000001*cEr[2]); 
  gkyl_mat_set(A,6,5,1.341640786499874*m0r[6]-0.4472135954999579*cEr[6]); 
  gkyl_mat_set(A,6,6,1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+0.9583148474999099*m0r[4]-0.31943828249997*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,6,7,1.2*m0r[3]-0.4*cEr[3]); 
  gkyl_mat_set(A,7,0,1.5*m0r[7]-0.5*cEr[7]); 
  gkyl_mat_set(A,7,1,1.5*m0r[5]-0.5000000000000001*cEr[5]); 
  gkyl_mat_set(A,7,2,1.341640786499874*m0r[3]-0.447213595499958*cEr[3]); 
  gkyl_mat_set(A,7,3,1.2*m0r[6]-0.4*cEr[6]+1.341640786499874*m0r[2]-0.447213595499958*cEr[2]); 
  gkyl_mat_set(A,7,4,1.341640786499874*m0r[7]-0.4472135954999579*cEr[7]); 
  gkyl_mat_set(A,7,5,0.9583148474999099*m0r[7]-0.31943828249997*cEr[7]+1.5*m0r[1]-0.5000000000000001*cEr[1]); 
  gkyl_mat_set(A,7,6,1.2*m0r[3]-0.4*cEr[3]); 
  gkyl_mat_set(A,7,7,0.9583148474999099*m0r[5]-0.31943828249997*cEr[5]+1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]); 
 
} 
 
