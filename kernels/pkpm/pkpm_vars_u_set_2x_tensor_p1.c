#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH int pkpm_vars_u_set_2x_tensor_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       [rho ux, rho uy, rho uz], Fluid input state vector.

  struct gkyl_mat A_ux = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_uy = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_uz = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  // Clear matrix and rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&A_ux, 0.0); gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&A_uy, 0.0); gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&A_uz, 0.0); gkyl_mat_clear(&rhs_uz, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[9]; 
  const double *rhouz = &euler_pkpm[18]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  int cell_avg = 0;
  // Check if rho < 0 at control points for diagnostic purposes. 
  if (2.5*rho[8]-1.936491673103709*rho[7]-1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]+1.5*rho[3]-0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]+1.936491673103709*rho[7]-1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]-1.5*rho[3]-0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]-1.936491673103709*rho[7]+1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]-1.5*rho[3]+0.8660254037844386*rho[2]-0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
  if (2.5*rho[8]+1.936491673103709*rho[7]+1.936491673103709*rho[6]+1.118033988749895*rho[5]+1.118033988749895*rho[4]+1.5*rho[3]+0.8660254037844386*rho[2]+0.8660254037844386*rho[1]+0.5*rho[0] < 0.0) cell_avg = 1; 
 
  gkyl_mat_set(&rhs_ux,0,0,rhoux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,rhouy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,rhouz[0]); 
  gkyl_mat_set(&rhs_ux,1,0,rhoux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,rhouy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,rhouz[1]); 
  gkyl_mat_set(&rhs_ux,2,0,rhoux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,rhouy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,rhouz[2]); 
  gkyl_mat_set(&rhs_ux,3,0,rhoux[3]); 
  gkyl_mat_set(&rhs_uy,3,0,rhouy[3]); 
  gkyl_mat_set(&rhs_uz,3,0,rhouz[3]); 
 
  double temp_rho = 0.0; 
  temp_rho = 0.5*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
 
  temp_rho = 0.5*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
 
  temp_rho = 0.5*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,0,3,temp_rho); 
  gkyl_mat_set(&A_uy,0,3,temp_rho); 
  gkyl_mat_set(&A_uz,0,3,temp_rho); 
 
  temp_rho = 0.5*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[6]+0.5*rho[2]; 
  gkyl_mat_set(&A_ux,1,3,temp_rho); 
  gkyl_mat_set(&A_uy,1,3,temp_rho); 
  gkyl_mat_set(&A_uz,1,3,temp_rho); 
 
  temp_rho = 0.5*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
 
  temp_rho = 0.4472135954999579*rho[5]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[7]+0.5*rho[1]; 
  gkyl_mat_set(&A_ux,2,3,temp_rho); 
  gkyl_mat_set(&A_uy,2,3,temp_rho); 
  gkyl_mat_set(&A_uz,2,3,temp_rho); 
 
  temp_rho = 0.5*rho[3]; 
  gkyl_mat_set(&A_ux,3,0,temp_rho); 
  gkyl_mat_set(&A_uy,3,0,temp_rho); 
  gkyl_mat_set(&A_uz,3,0,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[6]+0.5*rho[2]; 
  gkyl_mat_set(&A_ux,3,1,temp_rho); 
  gkyl_mat_set(&A_uy,3,1,temp_rho); 
  gkyl_mat_set(&A_uz,3,1,temp_rho); 
 
  temp_rho = 0.447213595499958*rho[7]+0.5*rho[1]; 
  gkyl_mat_set(&A_ux,3,2,temp_rho); 
  gkyl_mat_set(&A_uy,3,2,temp_rho); 
  gkyl_mat_set(&A_uz,3,2,temp_rho); 
 
  temp_rho = 0.4*rho[8]+0.4472135954999579*rho[5]+0.4472135954999579*rho[4]+0.5*rho[0]; 
  gkyl_mat_set(&A_ux,3,3,temp_rho); 
  gkyl_mat_set(&A_uy,3,3,temp_rho); 
  gkyl_mat_set(&A_uz,3,3,temp_rho); 
 
  return cell_avg;
} 
