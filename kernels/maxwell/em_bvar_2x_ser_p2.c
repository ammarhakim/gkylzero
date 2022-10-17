#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p2_exp_sq.h> 
#include <gkyl_basis_ser_2x_p2_inv.h> 
#include <gkyl_basis_ser_2x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_bvar_2x_ser_p2(const double *em, double* GKYL_RESTRICT bvar) 
{ 
  // em:   Input electromagnetic fields. 
  // bvar: b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
 
  const double *B_x = &em[24]; 
  const double *B_y = &em[32]; 
  const double *B_z = &em[40]; 
 
  double *bx = &bvar[0]; 
  double *by = &bvar[8]; 
  double *bz = &bvar[16]; 
  double *bxbx = &bvar[24]; 
  double *bxby = &bvar[32]; 
  double *bxbz = &bvar[40]; 
  double *byby = &bvar[48]; 
  double *bybz = &bvar[56]; 
  double *bzbz = &bvar[64]; 
 
  // Calculate |B|^2 and get expansion of 1/|B|^2. 
  double B_x_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(B_x, B_x_sq); 
 
  double B_y_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(B_y, B_y_sq); 
 
  double B_z_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(B_z, B_z_sq); 
 
  double magB_sq[8] = {0.0}; 

  double B_x_B_y[8] = {0.0}; 
  B_x_B_y[0] = 0.5*B_x[7]*B_y[7]+0.5*B_x[6]*B_y[6]+0.5*B_x[5]*B_y[5]+0.5*B_x[4]*B_y[4]+0.5*B_x[3]*B_y[3]+0.5*B_x[2]*B_y[2]+0.5*B_x[1]*B_y[1]+0.5*B_x[0]*B_y[0]; 
  B_x_B_y[1] = 0.5000000000000001*B_x[5]*B_y[7]+0.5000000000000001*B_y[5]*B_x[7]+0.447213595499958*B_x[3]*B_y[6]+0.447213595499958*B_y[3]*B_x[6]+0.4472135954999579*B_x[1]*B_y[4]+0.4472135954999579*B_y[1]*B_x[4]+0.5*B_x[2]*B_y[3]+0.5*B_y[2]*B_x[3]+0.5*B_x[0]*B_y[1]+0.5*B_y[0]*B_x[1]; 
  B_x_B_y[2] = 0.447213595499958*B_x[3]*B_y[7]+0.447213595499958*B_y[3]*B_x[7]+0.5000000000000001*B_x[4]*B_y[6]+0.5000000000000001*B_y[4]*B_x[6]+0.4472135954999579*B_x[2]*B_y[5]+0.4472135954999579*B_y[2]*B_x[5]+0.5*B_x[1]*B_y[3]+0.5*B_y[1]*B_x[3]+0.5*B_x[0]*B_y[2]+0.5*B_y[0]*B_x[2]; 
  B_x_B_y[3] = 0.4*B_x[6]*B_y[7]+0.447213595499958*B_x[2]*B_y[7]+0.4*B_y[6]*B_x[7]+0.447213595499958*B_y[2]*B_x[7]+0.447213595499958*B_x[1]*B_y[6]+0.447213595499958*B_y[1]*B_x[6]+0.4472135954999579*B_x[3]*B_y[5]+0.4472135954999579*B_y[3]*B_x[5]+0.4472135954999579*B_x[3]*B_y[4]+0.4472135954999579*B_y[3]*B_x[4]+0.5*B_x[0]*B_y[3]+0.5*B_y[0]*B_x[3]+0.5*B_x[1]*B_y[2]+0.5*B_y[1]*B_x[2]; 
  B_x_B_y[4] = 0.4472135954999579*B_x[7]*B_y[7]+0.31943828249997*B_x[6]*B_y[6]+0.5000000000000001*B_x[2]*B_y[6]+0.5000000000000001*B_y[2]*B_x[6]+0.31943828249997*B_x[4]*B_y[4]+0.5*B_x[0]*B_y[4]+0.5*B_y[0]*B_x[4]+0.4472135954999579*B_x[3]*B_y[3]+0.4472135954999579*B_x[1]*B_y[1]; 
  B_x_B_y[5] = 0.31943828249997*B_x[7]*B_y[7]+0.5000000000000001*B_x[1]*B_y[7]+0.5000000000000001*B_y[1]*B_x[7]+0.4472135954999579*B_x[6]*B_y[6]+0.31943828249997*B_x[5]*B_y[5]+0.5*B_x[0]*B_y[5]+0.5*B_y[0]*B_x[5]+0.4472135954999579*B_x[3]*B_y[3]+0.4472135954999579*B_x[2]*B_y[2]; 
  B_x_B_y[6] = 0.4*B_x[3]*B_y[7]+0.4*B_y[3]*B_x[7]+0.4472135954999579*B_x[5]*B_y[6]+0.31943828249997*B_x[4]*B_y[6]+0.5*B_x[0]*B_y[6]+0.4472135954999579*B_y[5]*B_x[6]+0.31943828249997*B_y[4]*B_x[6]+0.5*B_y[0]*B_x[6]+0.5000000000000001*B_x[2]*B_y[4]+0.5000000000000001*B_y[2]*B_x[4]+0.447213595499958*B_x[1]*B_y[3]+0.447213595499958*B_y[1]*B_x[3]; 
  B_x_B_y[7] = 0.31943828249997*B_x[5]*B_y[7]+0.4472135954999579*B_x[4]*B_y[7]+0.5*B_x[0]*B_y[7]+0.31943828249997*B_y[5]*B_x[7]+0.4472135954999579*B_y[4]*B_x[7]+0.5*B_y[0]*B_x[7]+0.4*B_x[3]*B_y[6]+0.4*B_y[3]*B_x[6]+0.5000000000000001*B_x[1]*B_y[5]+0.5000000000000001*B_y[1]*B_x[5]+0.447213595499958*B_x[2]*B_y[3]+0.447213595499958*B_y[2]*B_x[3]; 

  double B_x_B_z[8] = {0.0}; 
  B_x_B_z[0] = 0.5*B_x[7]*B_z[7]+0.5*B_x[6]*B_z[6]+0.5*B_x[5]*B_z[5]+0.5*B_x[4]*B_z[4]+0.5*B_x[3]*B_z[3]+0.5*B_x[2]*B_z[2]+0.5*B_x[1]*B_z[1]+0.5*B_x[0]*B_z[0]; 
  B_x_B_z[1] = 0.5000000000000001*B_x[5]*B_z[7]+0.5000000000000001*B_z[5]*B_x[7]+0.447213595499958*B_x[3]*B_z[6]+0.447213595499958*B_z[3]*B_x[6]+0.4472135954999579*B_x[1]*B_z[4]+0.4472135954999579*B_z[1]*B_x[4]+0.5*B_x[2]*B_z[3]+0.5*B_z[2]*B_x[3]+0.5*B_x[0]*B_z[1]+0.5*B_z[0]*B_x[1]; 
  B_x_B_z[2] = 0.447213595499958*B_x[3]*B_z[7]+0.447213595499958*B_z[3]*B_x[7]+0.5000000000000001*B_x[4]*B_z[6]+0.5000000000000001*B_z[4]*B_x[6]+0.4472135954999579*B_x[2]*B_z[5]+0.4472135954999579*B_z[2]*B_x[5]+0.5*B_x[1]*B_z[3]+0.5*B_z[1]*B_x[3]+0.5*B_x[0]*B_z[2]+0.5*B_z[0]*B_x[2]; 
  B_x_B_z[3] = 0.4*B_x[6]*B_z[7]+0.447213595499958*B_x[2]*B_z[7]+0.4*B_z[6]*B_x[7]+0.447213595499958*B_z[2]*B_x[7]+0.447213595499958*B_x[1]*B_z[6]+0.447213595499958*B_z[1]*B_x[6]+0.4472135954999579*B_x[3]*B_z[5]+0.4472135954999579*B_z[3]*B_x[5]+0.4472135954999579*B_x[3]*B_z[4]+0.4472135954999579*B_z[3]*B_x[4]+0.5*B_x[0]*B_z[3]+0.5*B_z[0]*B_x[3]+0.5*B_x[1]*B_z[2]+0.5*B_z[1]*B_x[2]; 
  B_x_B_z[4] = 0.4472135954999579*B_x[7]*B_z[7]+0.31943828249997*B_x[6]*B_z[6]+0.5000000000000001*B_x[2]*B_z[6]+0.5000000000000001*B_z[2]*B_x[6]+0.31943828249997*B_x[4]*B_z[4]+0.5*B_x[0]*B_z[4]+0.5*B_z[0]*B_x[4]+0.4472135954999579*B_x[3]*B_z[3]+0.4472135954999579*B_x[1]*B_z[1]; 
  B_x_B_z[5] = 0.31943828249997*B_x[7]*B_z[7]+0.5000000000000001*B_x[1]*B_z[7]+0.5000000000000001*B_z[1]*B_x[7]+0.4472135954999579*B_x[6]*B_z[6]+0.31943828249997*B_x[5]*B_z[5]+0.5*B_x[0]*B_z[5]+0.5*B_z[0]*B_x[5]+0.4472135954999579*B_x[3]*B_z[3]+0.4472135954999579*B_x[2]*B_z[2]; 
  B_x_B_z[6] = 0.4*B_x[3]*B_z[7]+0.4*B_z[3]*B_x[7]+0.4472135954999579*B_x[5]*B_z[6]+0.31943828249997*B_x[4]*B_z[6]+0.5*B_x[0]*B_z[6]+0.4472135954999579*B_z[5]*B_x[6]+0.31943828249997*B_z[4]*B_x[6]+0.5*B_z[0]*B_x[6]+0.5000000000000001*B_x[2]*B_z[4]+0.5000000000000001*B_z[2]*B_x[4]+0.447213595499958*B_x[1]*B_z[3]+0.447213595499958*B_z[1]*B_x[3]; 
  B_x_B_z[7] = 0.31943828249997*B_x[5]*B_z[7]+0.4472135954999579*B_x[4]*B_z[7]+0.5*B_x[0]*B_z[7]+0.31943828249997*B_z[5]*B_x[7]+0.4472135954999579*B_z[4]*B_x[7]+0.5*B_z[0]*B_x[7]+0.4*B_x[3]*B_z[6]+0.4*B_z[3]*B_x[6]+0.5000000000000001*B_x[1]*B_z[5]+0.5000000000000001*B_z[1]*B_x[5]+0.447213595499958*B_x[2]*B_z[3]+0.447213595499958*B_z[2]*B_x[3]; 

  double B_y_B_z[8] = {0.0}; 
  B_y_B_z[0] = 0.5*B_y[7]*B_z[7]+0.5*B_y[6]*B_z[6]+0.5*B_y[5]*B_z[5]+0.5*B_y[4]*B_z[4]+0.5*B_y[3]*B_z[3]+0.5*B_y[2]*B_z[2]+0.5*B_y[1]*B_z[1]+0.5*B_y[0]*B_z[0]; 
  B_y_B_z[1] = 0.5000000000000001*B_y[5]*B_z[7]+0.5000000000000001*B_z[5]*B_y[7]+0.447213595499958*B_y[3]*B_z[6]+0.447213595499958*B_z[3]*B_y[6]+0.4472135954999579*B_y[1]*B_z[4]+0.4472135954999579*B_z[1]*B_y[4]+0.5*B_y[2]*B_z[3]+0.5*B_z[2]*B_y[3]+0.5*B_y[0]*B_z[1]+0.5*B_z[0]*B_y[1]; 
  B_y_B_z[2] = 0.447213595499958*B_y[3]*B_z[7]+0.447213595499958*B_z[3]*B_y[7]+0.5000000000000001*B_y[4]*B_z[6]+0.5000000000000001*B_z[4]*B_y[6]+0.4472135954999579*B_y[2]*B_z[5]+0.4472135954999579*B_z[2]*B_y[5]+0.5*B_y[1]*B_z[3]+0.5*B_z[1]*B_y[3]+0.5*B_y[0]*B_z[2]+0.5*B_z[0]*B_y[2]; 
  B_y_B_z[3] = 0.4*B_y[6]*B_z[7]+0.447213595499958*B_y[2]*B_z[7]+0.4*B_z[6]*B_y[7]+0.447213595499958*B_z[2]*B_y[7]+0.447213595499958*B_y[1]*B_z[6]+0.447213595499958*B_z[1]*B_y[6]+0.4472135954999579*B_y[3]*B_z[5]+0.4472135954999579*B_z[3]*B_y[5]+0.4472135954999579*B_y[3]*B_z[4]+0.4472135954999579*B_z[3]*B_y[4]+0.5*B_y[0]*B_z[3]+0.5*B_z[0]*B_y[3]+0.5*B_y[1]*B_z[2]+0.5*B_z[1]*B_y[2]; 
  B_y_B_z[4] = 0.4472135954999579*B_y[7]*B_z[7]+0.31943828249997*B_y[6]*B_z[6]+0.5000000000000001*B_y[2]*B_z[6]+0.5000000000000001*B_z[2]*B_y[6]+0.31943828249997*B_y[4]*B_z[4]+0.5*B_y[0]*B_z[4]+0.5*B_z[0]*B_y[4]+0.4472135954999579*B_y[3]*B_z[3]+0.4472135954999579*B_y[1]*B_z[1]; 
  B_y_B_z[5] = 0.31943828249997*B_y[7]*B_z[7]+0.5000000000000001*B_y[1]*B_z[7]+0.5000000000000001*B_z[1]*B_y[7]+0.4472135954999579*B_y[6]*B_z[6]+0.31943828249997*B_y[5]*B_z[5]+0.5*B_y[0]*B_z[5]+0.5*B_z[0]*B_y[5]+0.4472135954999579*B_y[3]*B_z[3]+0.4472135954999579*B_y[2]*B_z[2]; 
  B_y_B_z[6] = 0.4*B_y[3]*B_z[7]+0.4*B_z[3]*B_y[7]+0.4472135954999579*B_y[5]*B_z[6]+0.31943828249997*B_y[4]*B_z[6]+0.5*B_y[0]*B_z[6]+0.4472135954999579*B_z[5]*B_y[6]+0.31943828249997*B_z[4]*B_y[6]+0.5*B_z[0]*B_y[6]+0.5000000000000001*B_y[2]*B_z[4]+0.5000000000000001*B_z[2]*B_y[4]+0.447213595499958*B_y[1]*B_z[3]+0.447213595499958*B_z[1]*B_y[3]; 
  B_y_B_z[7] = 0.31943828249997*B_y[5]*B_z[7]+0.4472135954999579*B_y[4]*B_z[7]+0.5*B_y[0]*B_z[7]+0.31943828249997*B_z[5]*B_y[7]+0.4472135954999579*B_z[4]*B_y[7]+0.5*B_z[0]*B_y[7]+0.4*B_y[3]*B_z[6]+0.4*B_z[3]*B_y[6]+0.5000000000000001*B_y[1]*B_z[5]+0.5000000000000001*B_z[1]*B_y[5]+0.447213595499958*B_y[2]*B_z[3]+0.447213595499958*B_z[2]*B_y[3]; 

  magB_sq[0] = B_z_sq[0]+B_y_sq[0]+B_x_sq[0]; 
  magB_sq[1] = B_z_sq[1]+B_y_sq[1]+B_x_sq[1]; 
  magB_sq[2] = B_z_sq[2]+B_y_sq[2]+B_x_sq[2]; 
  magB_sq[3] = B_z_sq[3]+B_y_sq[3]+B_x_sq[3]; 
  magB_sq[4] = B_z_sq[4]+B_y_sq[4]+B_x_sq[4]; 
  magB_sq[5] = B_z_sq[5]+B_y_sq[5]+B_x_sq[5]; 
  magB_sq[6] = B_z_sq[6]+B_y_sq[6]+B_x_sq[6]; 
  magB_sq[7] = B_z_sq[7]+B_y_sq[7]+B_x_sq[7]; 

  bool notCellAvg = true;
  if (notCellAvg && ((-1.936491673103709*magB_sq[7])-1.936491673103709*magB_sq[6]+1.118033988749895*magB_sq[5]+1.118033988749895*magB_sq[4]+1.5*magB_sq[3]-0.8660254037844386*magB_sq[2]-0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.936491673103709*magB_sq[7]-1.936491673103709*magB_sq[6]+1.118033988749895*magB_sq[5]+1.118033988749895*magB_sq[4]-1.5*magB_sq[3]-0.8660254037844386*magB_sq[2]+0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && ((-1.936491673103709*magB_sq[7])+1.936491673103709*magB_sq[6]+1.118033988749895*magB_sq[5]+1.118033988749895*magB_sq[4]-1.5*magB_sq[3]+0.8660254037844386*magB_sq[2]-0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  if (notCellAvg && (1.936491673103709*magB_sq[7]+1.936491673103709*magB_sq[6]+1.118033988749895*magB_sq[5]+1.118033988749895*magB_sq[4]+1.5*magB_sq[3]+0.8660254037844386*magB_sq[2]+0.8660254037844386*magB_sq[1]+0.5*magB_sq[0] < 0)) notCellAvg = false; 
  double magB_sq_inv[8] = {0.0}; 

  // Calculate expansions of B_i B_j/|B|^2, which can be calculated free of aliasing errors. 
  if (notCellAvg) { 
  ser_2x_p2_inv(magB_sq, magB_sq_inv); 
  bxbx[0] = 0.5*B_x_sq[7]*magB_sq_inv[7]+0.5*B_x_sq[6]*magB_sq_inv[6]+0.5*B_x_sq[5]*magB_sq_inv[5]+0.5*B_x_sq[4]*magB_sq_inv[4]+0.5*B_x_sq[3]*magB_sq_inv[3]+0.5*B_x_sq[2]*magB_sq_inv[2]+0.5*B_x_sq[1]*magB_sq_inv[1]+0.5*B_x_sq[0]*magB_sq_inv[0]; 
  bxbx[1] = 0.5000000000000001*B_x_sq[5]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[5]*B_x_sq[7]+0.447213595499958*B_x_sq[3]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[3]*B_x_sq[6]+0.4472135954999579*B_x_sq[1]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[1]*B_x_sq[4]+0.5*B_x_sq[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_x_sq[3]+0.5*B_x_sq[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x_sq[1]; 
  bxbx[2] = 0.447213595499958*B_x_sq[3]*magB_sq_inv[7]+0.447213595499958*magB_sq_inv[3]*B_x_sq[7]+0.5000000000000001*B_x_sq[4]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[4]*B_x_sq[6]+0.4472135954999579*B_x_sq[2]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[2]*B_x_sq[5]+0.5*B_x_sq[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_x_sq[3]+0.5*B_x_sq[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_x_sq[2]; 
  bxbx[3] = 0.4*B_x_sq[6]*magB_sq_inv[7]+0.447213595499958*B_x_sq[2]*magB_sq_inv[7]+0.4*magB_sq_inv[6]*B_x_sq[7]+0.447213595499958*magB_sq_inv[2]*B_x_sq[7]+0.447213595499958*B_x_sq[1]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[1]*B_x_sq[6]+0.4472135954999579*B_x_sq[3]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[3]*B_x_sq[5]+0.4472135954999579*B_x_sq[3]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[3]*B_x_sq[4]+0.5*B_x_sq[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_x_sq[3]+0.5*B_x_sq[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_x_sq[2]; 
  bxbx[4] = 0.4472135954999579*B_x_sq[7]*magB_sq_inv[7]+0.31943828249997*B_x_sq[6]*magB_sq_inv[6]+0.5000000000000001*B_x_sq[2]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[2]*B_x_sq[6]+0.31943828249997*B_x_sq[4]*magB_sq_inv[4]+0.5*B_x_sq[0]*magB_sq_inv[4]+0.5*magB_sq_inv[0]*B_x_sq[4]+0.4472135954999579*B_x_sq[3]*magB_sq_inv[3]+0.4472135954999579*B_x_sq[1]*magB_sq_inv[1]; 
  bxbx[5] = 0.31943828249997*B_x_sq[7]*magB_sq_inv[7]+0.5000000000000001*B_x_sq[1]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[1]*B_x_sq[7]+0.4472135954999579*B_x_sq[6]*magB_sq_inv[6]+0.31943828249997*B_x_sq[5]*magB_sq_inv[5]+0.5*B_x_sq[0]*magB_sq_inv[5]+0.5*magB_sq_inv[0]*B_x_sq[5]+0.4472135954999579*B_x_sq[3]*magB_sq_inv[3]+0.4472135954999579*B_x_sq[2]*magB_sq_inv[2]; 
  bxbx[6] = 0.4*B_x_sq[3]*magB_sq_inv[7]+0.4*magB_sq_inv[3]*B_x_sq[7]+0.4472135954999579*B_x_sq[5]*magB_sq_inv[6]+0.31943828249997*B_x_sq[4]*magB_sq_inv[6]+0.5*B_x_sq[0]*magB_sq_inv[6]+0.4472135954999579*magB_sq_inv[5]*B_x_sq[6]+0.31943828249997*magB_sq_inv[4]*B_x_sq[6]+0.5*magB_sq_inv[0]*B_x_sq[6]+0.5000000000000001*B_x_sq[2]*magB_sq_inv[4]+0.5000000000000001*magB_sq_inv[2]*B_x_sq[4]+0.447213595499958*B_x_sq[1]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[1]*B_x_sq[3]; 
  bxbx[7] = 0.31943828249997*B_x_sq[5]*magB_sq_inv[7]+0.4472135954999579*B_x_sq[4]*magB_sq_inv[7]+0.5*B_x_sq[0]*magB_sq_inv[7]+0.31943828249997*magB_sq_inv[5]*B_x_sq[7]+0.4472135954999579*magB_sq_inv[4]*B_x_sq[7]+0.5*magB_sq_inv[0]*B_x_sq[7]+0.4*B_x_sq[3]*magB_sq_inv[6]+0.4*magB_sq_inv[3]*B_x_sq[6]+0.5000000000000001*B_x_sq[1]*magB_sq_inv[5]+0.5000000000000001*magB_sq_inv[1]*B_x_sq[5]+0.447213595499958*B_x_sq[2]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[2]*B_x_sq[3]; 

  bxby[0] = 0.5*B_x_B_y[7]*magB_sq_inv[7]+0.5*B_x_B_y[6]*magB_sq_inv[6]+0.5*B_x_B_y[5]*magB_sq_inv[5]+0.5*B_x_B_y[4]*magB_sq_inv[4]+0.5*B_x_B_y[3]*magB_sq_inv[3]+0.5*B_x_B_y[2]*magB_sq_inv[2]+0.5*B_x_B_y[1]*magB_sq_inv[1]+0.5*B_x_B_y[0]*magB_sq_inv[0]; 
  bxby[1] = 0.5000000000000001*B_x_B_y[5]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[5]*B_x_B_y[7]+0.447213595499958*B_x_B_y[3]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[3]*B_x_B_y[6]+0.4472135954999579*B_x_B_y[1]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[1]*B_x_B_y[4]+0.5*B_x_B_y[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_x_B_y[3]+0.5*B_x_B_y[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x_B_y[1]; 
  bxby[2] = 0.447213595499958*B_x_B_y[3]*magB_sq_inv[7]+0.447213595499958*magB_sq_inv[3]*B_x_B_y[7]+0.5000000000000001*B_x_B_y[4]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[4]*B_x_B_y[6]+0.4472135954999579*B_x_B_y[2]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[2]*B_x_B_y[5]+0.5*B_x_B_y[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_x_B_y[3]+0.5*B_x_B_y[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_x_B_y[2]; 
  bxby[3] = 0.4*B_x_B_y[6]*magB_sq_inv[7]+0.447213595499958*B_x_B_y[2]*magB_sq_inv[7]+0.4*magB_sq_inv[6]*B_x_B_y[7]+0.447213595499958*magB_sq_inv[2]*B_x_B_y[7]+0.447213595499958*B_x_B_y[1]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[1]*B_x_B_y[6]+0.4472135954999579*B_x_B_y[3]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[3]*B_x_B_y[5]+0.4472135954999579*B_x_B_y[3]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[3]*B_x_B_y[4]+0.5*B_x_B_y[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_x_B_y[3]+0.5*B_x_B_y[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_x_B_y[2]; 
  bxby[4] = 0.4472135954999579*B_x_B_y[7]*magB_sq_inv[7]+0.31943828249997*B_x_B_y[6]*magB_sq_inv[6]+0.5000000000000001*B_x_B_y[2]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[2]*B_x_B_y[6]+0.31943828249997*B_x_B_y[4]*magB_sq_inv[4]+0.5*B_x_B_y[0]*magB_sq_inv[4]+0.5*magB_sq_inv[0]*B_x_B_y[4]+0.4472135954999579*B_x_B_y[3]*magB_sq_inv[3]+0.4472135954999579*B_x_B_y[1]*magB_sq_inv[1]; 
  bxby[5] = 0.31943828249997*B_x_B_y[7]*magB_sq_inv[7]+0.5000000000000001*B_x_B_y[1]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[1]*B_x_B_y[7]+0.4472135954999579*B_x_B_y[6]*magB_sq_inv[6]+0.31943828249997*B_x_B_y[5]*magB_sq_inv[5]+0.5*B_x_B_y[0]*magB_sq_inv[5]+0.5*magB_sq_inv[0]*B_x_B_y[5]+0.4472135954999579*B_x_B_y[3]*magB_sq_inv[3]+0.4472135954999579*B_x_B_y[2]*magB_sq_inv[2]; 
  bxby[6] = 0.4*B_x_B_y[3]*magB_sq_inv[7]+0.4*magB_sq_inv[3]*B_x_B_y[7]+0.4472135954999579*B_x_B_y[5]*magB_sq_inv[6]+0.31943828249997*B_x_B_y[4]*magB_sq_inv[6]+0.5*B_x_B_y[0]*magB_sq_inv[6]+0.4472135954999579*magB_sq_inv[5]*B_x_B_y[6]+0.31943828249997*magB_sq_inv[4]*B_x_B_y[6]+0.5*magB_sq_inv[0]*B_x_B_y[6]+0.5000000000000001*B_x_B_y[2]*magB_sq_inv[4]+0.5000000000000001*magB_sq_inv[2]*B_x_B_y[4]+0.447213595499958*B_x_B_y[1]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[1]*B_x_B_y[3]; 
  bxby[7] = 0.31943828249997*B_x_B_y[5]*magB_sq_inv[7]+0.4472135954999579*B_x_B_y[4]*magB_sq_inv[7]+0.5*B_x_B_y[0]*magB_sq_inv[7]+0.31943828249997*magB_sq_inv[5]*B_x_B_y[7]+0.4472135954999579*magB_sq_inv[4]*B_x_B_y[7]+0.5*magB_sq_inv[0]*B_x_B_y[7]+0.4*B_x_B_y[3]*magB_sq_inv[6]+0.4*magB_sq_inv[3]*B_x_B_y[6]+0.5000000000000001*B_x_B_y[1]*magB_sq_inv[5]+0.5000000000000001*magB_sq_inv[1]*B_x_B_y[5]+0.447213595499958*B_x_B_y[2]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[2]*B_x_B_y[3]; 

  bxbz[0] = 0.5*B_x_B_z[7]*magB_sq_inv[7]+0.5*B_x_B_z[6]*magB_sq_inv[6]+0.5*B_x_B_z[5]*magB_sq_inv[5]+0.5*B_x_B_z[4]*magB_sq_inv[4]+0.5*B_x_B_z[3]*magB_sq_inv[3]+0.5*B_x_B_z[2]*magB_sq_inv[2]+0.5*B_x_B_z[1]*magB_sq_inv[1]+0.5*B_x_B_z[0]*magB_sq_inv[0]; 
  bxbz[1] = 0.5000000000000001*B_x_B_z[5]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[5]*B_x_B_z[7]+0.447213595499958*B_x_B_z[3]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[3]*B_x_B_z[6]+0.4472135954999579*B_x_B_z[1]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[1]*B_x_B_z[4]+0.5*B_x_B_z[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_x_B_z[3]+0.5*B_x_B_z[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_x_B_z[1]; 
  bxbz[2] = 0.447213595499958*B_x_B_z[3]*magB_sq_inv[7]+0.447213595499958*magB_sq_inv[3]*B_x_B_z[7]+0.5000000000000001*B_x_B_z[4]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[4]*B_x_B_z[6]+0.4472135954999579*B_x_B_z[2]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[2]*B_x_B_z[5]+0.5*B_x_B_z[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_x_B_z[3]+0.5*B_x_B_z[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_x_B_z[2]; 
  bxbz[3] = 0.4*B_x_B_z[6]*magB_sq_inv[7]+0.447213595499958*B_x_B_z[2]*magB_sq_inv[7]+0.4*magB_sq_inv[6]*B_x_B_z[7]+0.447213595499958*magB_sq_inv[2]*B_x_B_z[7]+0.447213595499958*B_x_B_z[1]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[1]*B_x_B_z[6]+0.4472135954999579*B_x_B_z[3]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[3]*B_x_B_z[5]+0.4472135954999579*B_x_B_z[3]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[3]*B_x_B_z[4]+0.5*B_x_B_z[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_x_B_z[3]+0.5*B_x_B_z[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_x_B_z[2]; 
  bxbz[4] = 0.4472135954999579*B_x_B_z[7]*magB_sq_inv[7]+0.31943828249997*B_x_B_z[6]*magB_sq_inv[6]+0.5000000000000001*B_x_B_z[2]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[2]*B_x_B_z[6]+0.31943828249997*B_x_B_z[4]*magB_sq_inv[4]+0.5*B_x_B_z[0]*magB_sq_inv[4]+0.5*magB_sq_inv[0]*B_x_B_z[4]+0.4472135954999579*B_x_B_z[3]*magB_sq_inv[3]+0.4472135954999579*B_x_B_z[1]*magB_sq_inv[1]; 
  bxbz[5] = 0.31943828249997*B_x_B_z[7]*magB_sq_inv[7]+0.5000000000000001*B_x_B_z[1]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[1]*B_x_B_z[7]+0.4472135954999579*B_x_B_z[6]*magB_sq_inv[6]+0.31943828249997*B_x_B_z[5]*magB_sq_inv[5]+0.5*B_x_B_z[0]*magB_sq_inv[5]+0.5*magB_sq_inv[0]*B_x_B_z[5]+0.4472135954999579*B_x_B_z[3]*magB_sq_inv[3]+0.4472135954999579*B_x_B_z[2]*magB_sq_inv[2]; 
  bxbz[6] = 0.4*B_x_B_z[3]*magB_sq_inv[7]+0.4*magB_sq_inv[3]*B_x_B_z[7]+0.4472135954999579*B_x_B_z[5]*magB_sq_inv[6]+0.31943828249997*B_x_B_z[4]*magB_sq_inv[6]+0.5*B_x_B_z[0]*magB_sq_inv[6]+0.4472135954999579*magB_sq_inv[5]*B_x_B_z[6]+0.31943828249997*magB_sq_inv[4]*B_x_B_z[6]+0.5*magB_sq_inv[0]*B_x_B_z[6]+0.5000000000000001*B_x_B_z[2]*magB_sq_inv[4]+0.5000000000000001*magB_sq_inv[2]*B_x_B_z[4]+0.447213595499958*B_x_B_z[1]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[1]*B_x_B_z[3]; 
  bxbz[7] = 0.31943828249997*B_x_B_z[5]*magB_sq_inv[7]+0.4472135954999579*B_x_B_z[4]*magB_sq_inv[7]+0.5*B_x_B_z[0]*magB_sq_inv[7]+0.31943828249997*magB_sq_inv[5]*B_x_B_z[7]+0.4472135954999579*magB_sq_inv[4]*B_x_B_z[7]+0.5*magB_sq_inv[0]*B_x_B_z[7]+0.4*B_x_B_z[3]*magB_sq_inv[6]+0.4*magB_sq_inv[3]*B_x_B_z[6]+0.5000000000000001*B_x_B_z[1]*magB_sq_inv[5]+0.5000000000000001*magB_sq_inv[1]*B_x_B_z[5]+0.447213595499958*B_x_B_z[2]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[2]*B_x_B_z[3]; 

  byby[0] = 0.5*B_y_sq[7]*magB_sq_inv[7]+0.5*B_y_sq[6]*magB_sq_inv[6]+0.5*B_y_sq[5]*magB_sq_inv[5]+0.5*B_y_sq[4]*magB_sq_inv[4]+0.5*B_y_sq[3]*magB_sq_inv[3]+0.5*B_y_sq[2]*magB_sq_inv[2]+0.5*B_y_sq[1]*magB_sq_inv[1]+0.5*B_y_sq[0]*magB_sq_inv[0]; 
  byby[1] = 0.5000000000000001*B_y_sq[5]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[5]*B_y_sq[7]+0.447213595499958*B_y_sq[3]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[3]*B_y_sq[6]+0.4472135954999579*B_y_sq[1]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[1]*B_y_sq[4]+0.5*B_y_sq[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_y_sq[3]+0.5*B_y_sq[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_y_sq[1]; 
  byby[2] = 0.447213595499958*B_y_sq[3]*magB_sq_inv[7]+0.447213595499958*magB_sq_inv[3]*B_y_sq[7]+0.5000000000000001*B_y_sq[4]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[4]*B_y_sq[6]+0.4472135954999579*B_y_sq[2]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[2]*B_y_sq[5]+0.5*B_y_sq[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_y_sq[3]+0.5*B_y_sq[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_y_sq[2]; 
  byby[3] = 0.4*B_y_sq[6]*magB_sq_inv[7]+0.447213595499958*B_y_sq[2]*magB_sq_inv[7]+0.4*magB_sq_inv[6]*B_y_sq[7]+0.447213595499958*magB_sq_inv[2]*B_y_sq[7]+0.447213595499958*B_y_sq[1]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[1]*B_y_sq[6]+0.4472135954999579*B_y_sq[3]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[3]*B_y_sq[5]+0.4472135954999579*B_y_sq[3]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[3]*B_y_sq[4]+0.5*B_y_sq[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_y_sq[3]+0.5*B_y_sq[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_y_sq[2]; 
  byby[4] = 0.4472135954999579*B_y_sq[7]*magB_sq_inv[7]+0.31943828249997*B_y_sq[6]*magB_sq_inv[6]+0.5000000000000001*B_y_sq[2]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[2]*B_y_sq[6]+0.31943828249997*B_y_sq[4]*magB_sq_inv[4]+0.5*B_y_sq[0]*magB_sq_inv[4]+0.5*magB_sq_inv[0]*B_y_sq[4]+0.4472135954999579*B_y_sq[3]*magB_sq_inv[3]+0.4472135954999579*B_y_sq[1]*magB_sq_inv[1]; 
  byby[5] = 0.31943828249997*B_y_sq[7]*magB_sq_inv[7]+0.5000000000000001*B_y_sq[1]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[1]*B_y_sq[7]+0.4472135954999579*B_y_sq[6]*magB_sq_inv[6]+0.31943828249997*B_y_sq[5]*magB_sq_inv[5]+0.5*B_y_sq[0]*magB_sq_inv[5]+0.5*magB_sq_inv[0]*B_y_sq[5]+0.4472135954999579*B_y_sq[3]*magB_sq_inv[3]+0.4472135954999579*B_y_sq[2]*magB_sq_inv[2]; 
  byby[6] = 0.4*B_y_sq[3]*magB_sq_inv[7]+0.4*magB_sq_inv[3]*B_y_sq[7]+0.4472135954999579*B_y_sq[5]*magB_sq_inv[6]+0.31943828249997*B_y_sq[4]*magB_sq_inv[6]+0.5*B_y_sq[0]*magB_sq_inv[6]+0.4472135954999579*magB_sq_inv[5]*B_y_sq[6]+0.31943828249997*magB_sq_inv[4]*B_y_sq[6]+0.5*magB_sq_inv[0]*B_y_sq[6]+0.5000000000000001*B_y_sq[2]*magB_sq_inv[4]+0.5000000000000001*magB_sq_inv[2]*B_y_sq[4]+0.447213595499958*B_y_sq[1]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[1]*B_y_sq[3]; 
  byby[7] = 0.31943828249997*B_y_sq[5]*magB_sq_inv[7]+0.4472135954999579*B_y_sq[4]*magB_sq_inv[7]+0.5*B_y_sq[0]*magB_sq_inv[7]+0.31943828249997*magB_sq_inv[5]*B_y_sq[7]+0.4472135954999579*magB_sq_inv[4]*B_y_sq[7]+0.5*magB_sq_inv[0]*B_y_sq[7]+0.4*B_y_sq[3]*magB_sq_inv[6]+0.4*magB_sq_inv[3]*B_y_sq[6]+0.5000000000000001*B_y_sq[1]*magB_sq_inv[5]+0.5000000000000001*magB_sq_inv[1]*B_y_sq[5]+0.447213595499958*B_y_sq[2]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[2]*B_y_sq[3]; 

  bybz[0] = 0.5*B_y_B_z[7]*magB_sq_inv[7]+0.5*B_y_B_z[6]*magB_sq_inv[6]+0.5*B_y_B_z[5]*magB_sq_inv[5]+0.5*B_y_B_z[4]*magB_sq_inv[4]+0.5*B_y_B_z[3]*magB_sq_inv[3]+0.5*B_y_B_z[2]*magB_sq_inv[2]+0.5*B_y_B_z[1]*magB_sq_inv[1]+0.5*B_y_B_z[0]*magB_sq_inv[0]; 
  bybz[1] = 0.5000000000000001*B_y_B_z[5]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[5]*B_y_B_z[7]+0.447213595499958*B_y_B_z[3]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[3]*B_y_B_z[6]+0.4472135954999579*B_y_B_z[1]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[1]*B_y_B_z[4]+0.5*B_y_B_z[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_y_B_z[3]+0.5*B_y_B_z[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_y_B_z[1]; 
  bybz[2] = 0.447213595499958*B_y_B_z[3]*magB_sq_inv[7]+0.447213595499958*magB_sq_inv[3]*B_y_B_z[7]+0.5000000000000001*B_y_B_z[4]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[4]*B_y_B_z[6]+0.4472135954999579*B_y_B_z[2]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[2]*B_y_B_z[5]+0.5*B_y_B_z[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_y_B_z[3]+0.5*B_y_B_z[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_y_B_z[2]; 
  bybz[3] = 0.4*B_y_B_z[6]*magB_sq_inv[7]+0.447213595499958*B_y_B_z[2]*magB_sq_inv[7]+0.4*magB_sq_inv[6]*B_y_B_z[7]+0.447213595499958*magB_sq_inv[2]*B_y_B_z[7]+0.447213595499958*B_y_B_z[1]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[1]*B_y_B_z[6]+0.4472135954999579*B_y_B_z[3]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[3]*B_y_B_z[5]+0.4472135954999579*B_y_B_z[3]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[3]*B_y_B_z[4]+0.5*B_y_B_z[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_y_B_z[3]+0.5*B_y_B_z[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_y_B_z[2]; 
  bybz[4] = 0.4472135954999579*B_y_B_z[7]*magB_sq_inv[7]+0.31943828249997*B_y_B_z[6]*magB_sq_inv[6]+0.5000000000000001*B_y_B_z[2]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[2]*B_y_B_z[6]+0.31943828249997*B_y_B_z[4]*magB_sq_inv[4]+0.5*B_y_B_z[0]*magB_sq_inv[4]+0.5*magB_sq_inv[0]*B_y_B_z[4]+0.4472135954999579*B_y_B_z[3]*magB_sq_inv[3]+0.4472135954999579*B_y_B_z[1]*magB_sq_inv[1]; 
  bybz[5] = 0.31943828249997*B_y_B_z[7]*magB_sq_inv[7]+0.5000000000000001*B_y_B_z[1]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[1]*B_y_B_z[7]+0.4472135954999579*B_y_B_z[6]*magB_sq_inv[6]+0.31943828249997*B_y_B_z[5]*magB_sq_inv[5]+0.5*B_y_B_z[0]*magB_sq_inv[5]+0.5*magB_sq_inv[0]*B_y_B_z[5]+0.4472135954999579*B_y_B_z[3]*magB_sq_inv[3]+0.4472135954999579*B_y_B_z[2]*magB_sq_inv[2]; 
  bybz[6] = 0.4*B_y_B_z[3]*magB_sq_inv[7]+0.4*magB_sq_inv[3]*B_y_B_z[7]+0.4472135954999579*B_y_B_z[5]*magB_sq_inv[6]+0.31943828249997*B_y_B_z[4]*magB_sq_inv[6]+0.5*B_y_B_z[0]*magB_sq_inv[6]+0.4472135954999579*magB_sq_inv[5]*B_y_B_z[6]+0.31943828249997*magB_sq_inv[4]*B_y_B_z[6]+0.5*magB_sq_inv[0]*B_y_B_z[6]+0.5000000000000001*B_y_B_z[2]*magB_sq_inv[4]+0.5000000000000001*magB_sq_inv[2]*B_y_B_z[4]+0.447213595499958*B_y_B_z[1]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[1]*B_y_B_z[3]; 
  bybz[7] = 0.31943828249997*B_y_B_z[5]*magB_sq_inv[7]+0.4472135954999579*B_y_B_z[4]*magB_sq_inv[7]+0.5*B_y_B_z[0]*magB_sq_inv[7]+0.31943828249997*magB_sq_inv[5]*B_y_B_z[7]+0.4472135954999579*magB_sq_inv[4]*B_y_B_z[7]+0.5*magB_sq_inv[0]*B_y_B_z[7]+0.4*B_y_B_z[3]*magB_sq_inv[6]+0.4*magB_sq_inv[3]*B_y_B_z[6]+0.5000000000000001*B_y_B_z[1]*magB_sq_inv[5]+0.5000000000000001*magB_sq_inv[1]*B_y_B_z[5]+0.447213595499958*B_y_B_z[2]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[2]*B_y_B_z[3]; 

  bzbz[0] = 0.5*B_z_sq[7]*magB_sq_inv[7]+0.5*B_z_sq[6]*magB_sq_inv[6]+0.5*B_z_sq[5]*magB_sq_inv[5]+0.5*B_z_sq[4]*magB_sq_inv[4]+0.5*B_z_sq[3]*magB_sq_inv[3]+0.5*B_z_sq[2]*magB_sq_inv[2]+0.5*B_z_sq[1]*magB_sq_inv[1]+0.5*B_z_sq[0]*magB_sq_inv[0]; 
  bzbz[1] = 0.5000000000000001*B_z_sq[5]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[5]*B_z_sq[7]+0.447213595499958*B_z_sq[3]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[3]*B_z_sq[6]+0.4472135954999579*B_z_sq[1]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[1]*B_z_sq[4]+0.5*B_z_sq[2]*magB_sq_inv[3]+0.5*magB_sq_inv[2]*B_z_sq[3]+0.5*B_z_sq[0]*magB_sq_inv[1]+0.5*magB_sq_inv[0]*B_z_sq[1]; 
  bzbz[2] = 0.447213595499958*B_z_sq[3]*magB_sq_inv[7]+0.447213595499958*magB_sq_inv[3]*B_z_sq[7]+0.5000000000000001*B_z_sq[4]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[4]*B_z_sq[6]+0.4472135954999579*B_z_sq[2]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[2]*B_z_sq[5]+0.5*B_z_sq[1]*magB_sq_inv[3]+0.5*magB_sq_inv[1]*B_z_sq[3]+0.5*B_z_sq[0]*magB_sq_inv[2]+0.5*magB_sq_inv[0]*B_z_sq[2]; 
  bzbz[3] = 0.4*B_z_sq[6]*magB_sq_inv[7]+0.447213595499958*B_z_sq[2]*magB_sq_inv[7]+0.4*magB_sq_inv[6]*B_z_sq[7]+0.447213595499958*magB_sq_inv[2]*B_z_sq[7]+0.447213595499958*B_z_sq[1]*magB_sq_inv[6]+0.447213595499958*magB_sq_inv[1]*B_z_sq[6]+0.4472135954999579*B_z_sq[3]*magB_sq_inv[5]+0.4472135954999579*magB_sq_inv[3]*B_z_sq[5]+0.4472135954999579*B_z_sq[3]*magB_sq_inv[4]+0.4472135954999579*magB_sq_inv[3]*B_z_sq[4]+0.5*B_z_sq[0]*magB_sq_inv[3]+0.5*magB_sq_inv[0]*B_z_sq[3]+0.5*B_z_sq[1]*magB_sq_inv[2]+0.5*magB_sq_inv[1]*B_z_sq[2]; 
  bzbz[4] = 0.4472135954999579*B_z_sq[7]*magB_sq_inv[7]+0.31943828249997*B_z_sq[6]*magB_sq_inv[6]+0.5000000000000001*B_z_sq[2]*magB_sq_inv[6]+0.5000000000000001*magB_sq_inv[2]*B_z_sq[6]+0.31943828249997*B_z_sq[4]*magB_sq_inv[4]+0.5*B_z_sq[0]*magB_sq_inv[4]+0.5*magB_sq_inv[0]*B_z_sq[4]+0.4472135954999579*B_z_sq[3]*magB_sq_inv[3]+0.4472135954999579*B_z_sq[1]*magB_sq_inv[1]; 
  bzbz[5] = 0.31943828249997*B_z_sq[7]*magB_sq_inv[7]+0.5000000000000001*B_z_sq[1]*magB_sq_inv[7]+0.5000000000000001*magB_sq_inv[1]*B_z_sq[7]+0.4472135954999579*B_z_sq[6]*magB_sq_inv[6]+0.31943828249997*B_z_sq[5]*magB_sq_inv[5]+0.5*B_z_sq[0]*magB_sq_inv[5]+0.5*magB_sq_inv[0]*B_z_sq[5]+0.4472135954999579*B_z_sq[3]*magB_sq_inv[3]+0.4472135954999579*B_z_sq[2]*magB_sq_inv[2]; 
  bzbz[6] = 0.4*B_z_sq[3]*magB_sq_inv[7]+0.4*magB_sq_inv[3]*B_z_sq[7]+0.4472135954999579*B_z_sq[5]*magB_sq_inv[6]+0.31943828249997*B_z_sq[4]*magB_sq_inv[6]+0.5*B_z_sq[0]*magB_sq_inv[6]+0.4472135954999579*magB_sq_inv[5]*B_z_sq[6]+0.31943828249997*magB_sq_inv[4]*B_z_sq[6]+0.5*magB_sq_inv[0]*B_z_sq[6]+0.5000000000000001*B_z_sq[2]*magB_sq_inv[4]+0.5000000000000001*magB_sq_inv[2]*B_z_sq[4]+0.447213595499958*B_z_sq[1]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[1]*B_z_sq[3]; 
  bzbz[7] = 0.31943828249997*B_z_sq[5]*magB_sq_inv[7]+0.4472135954999579*B_z_sq[4]*magB_sq_inv[7]+0.5*B_z_sq[0]*magB_sq_inv[7]+0.31943828249997*magB_sq_inv[5]*B_z_sq[7]+0.4472135954999579*magB_sq_inv[4]*B_z_sq[7]+0.5*magB_sq_inv[0]*B_z_sq[7]+0.4*B_z_sq[3]*magB_sq_inv[6]+0.4*magB_sq_inv[3]*B_z_sq[6]+0.5000000000000001*B_z_sq[1]*magB_sq_inv[5]+0.5000000000000001*magB_sq_inv[1]*B_z_sq[5]+0.447213595499958*B_z_sq[2]*magB_sq_inv[3]+0.447213595499958*magB_sq_inv[2]*B_z_sq[3]; 

  } else { 
  magB_sq_inv[0] = 4.0/magB_sq[0]; 
  bxbx[0] = 0.5*B_x_sq[7]*magB_sq_inv[7]+0.5*B_x_sq[6]*magB_sq_inv[6]+0.5*B_x_sq[5]*magB_sq_inv[5]+0.5*B_x_sq[4]*magB_sq_inv[4]+0.5*B_x_sq[3]*magB_sq_inv[3]+0.5*B_x_sq[2]*magB_sq_inv[2]+0.5*B_x_sq[1]*magB_sq_inv[1]+0.5*B_x_sq[0]*magB_sq_inv[0]; 
  bxby[0] = 0.5*B_x_B_y[7]*magB_sq_inv[7]+0.5*B_x_B_y[6]*magB_sq_inv[6]+0.5*B_x_B_y[5]*magB_sq_inv[5]+0.5*B_x_B_y[4]*magB_sq_inv[4]+0.5*B_x_B_y[3]*magB_sq_inv[3]+0.5*B_x_B_y[2]*magB_sq_inv[2]+0.5*B_x_B_y[1]*magB_sq_inv[1]+0.5*B_x_B_y[0]*magB_sq_inv[0]; 
  bxbz[0] = 0.5*B_x_B_z[7]*magB_sq_inv[7]+0.5*B_x_B_z[6]*magB_sq_inv[6]+0.5*B_x_B_z[5]*magB_sq_inv[5]+0.5*B_x_B_z[4]*magB_sq_inv[4]+0.5*B_x_B_z[3]*magB_sq_inv[3]+0.5*B_x_B_z[2]*magB_sq_inv[2]+0.5*B_x_B_z[1]*magB_sq_inv[1]+0.5*B_x_B_z[0]*magB_sq_inv[0]; 
  byby[0] = 0.5*B_y_sq[7]*magB_sq_inv[7]+0.5*B_y_sq[6]*magB_sq_inv[6]+0.5*B_y_sq[5]*magB_sq_inv[5]+0.5*B_y_sq[4]*magB_sq_inv[4]+0.5*B_y_sq[3]*magB_sq_inv[3]+0.5*B_y_sq[2]*magB_sq_inv[2]+0.5*B_y_sq[1]*magB_sq_inv[1]+0.5*B_y_sq[0]*magB_sq_inv[0]; 
  bybz[0] = 0.5*B_y_B_z[7]*magB_sq_inv[7]+0.5*B_y_B_z[6]*magB_sq_inv[6]+0.5*B_y_B_z[5]*magB_sq_inv[5]+0.5*B_y_B_z[4]*magB_sq_inv[4]+0.5*B_y_B_z[3]*magB_sq_inv[3]+0.5*B_y_B_z[2]*magB_sq_inv[2]+0.5*B_y_B_z[1]*magB_sq_inv[1]+0.5*B_y_B_z[0]*magB_sq_inv[0]; 
  bzbz[0] = 0.5*B_z_sq[7]*magB_sq_inv[7]+0.5*B_z_sq[6]*magB_sq_inv[6]+0.5*B_z_sq[5]*magB_sq_inv[5]+0.5*B_z_sq[4]*magB_sq_inv[4]+0.5*B_z_sq[3]*magB_sq_inv[3]+0.5*B_z_sq[2]*magB_sq_inv[2]+0.5*B_z_sq[1]*magB_sq_inv[1]+0.5*B_z_sq[0]*magB_sq_inv[0]; 

  bxbx[1] = 0.0; 
  bxby[1] = 0.0; 
  bxbz[1] = 0.0; 
  byby[1] = 0.0; 
  bybz[1] = 0.0; 
  bzbz[1] = 0.0; 

  bxbx[2] = 0.0; 
  bxby[2] = 0.0; 
  bxbz[2] = 0.0; 
  byby[2] = 0.0; 
  bybz[2] = 0.0; 
  bzbz[2] = 0.0; 

  bxbx[3] = 0.0; 
  bxby[3] = 0.0; 
  bxbz[3] = 0.0; 
  byby[3] = 0.0; 
  bybz[3] = 0.0; 
  bzbz[3] = 0.0; 

  bxbx[4] = 0.0; 
  bxby[4] = 0.0; 
  bxbz[4] = 0.0; 
  byby[4] = 0.0; 
  bybz[4] = 0.0; 
  bzbz[4] = 0.0; 

  bxbx[5] = 0.0; 
  bxby[5] = 0.0; 
  bxbz[5] = 0.0; 
  byby[5] = 0.0; 
  bybz[5] = 0.0; 
  bzbz[5] = 0.0; 

  bxbx[6] = 0.0; 
  bxby[6] = 0.0; 
  bxbz[6] = 0.0; 
  byby[6] = 0.0; 
  bybz[6] = 0.0; 
  bzbz[6] = 0.0; 

  bxbx[7] = 0.0; 
  bxby[7] = 0.0; 
  bxbz[7] = 0.0; 
  byby[7] = 0.0; 
  bybz[7] = 0.0; 
  bzbz[7] = 0.0; 

  } 
  // Calculate b_i = B_i/|B| by taking square root of B_i^2/|B|^2 at quadrature points. 
  ser_2x_p2_sqrt_with_sign(bxbx, bx); 
  ser_2x_p2_sqrt_with_sign(byby, by); 
  ser_2x_p2_sqrt_with_sign(bzbz, bz); 
} 
 
