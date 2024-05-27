#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH int pkpm_vars_u_set_3x_tensor_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  const double *rhouy = &euler_pkpm[8]; 
  const double *rhouz = &euler_pkpm[16]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  int cell_avg = 0;
  // Check if rho < 0 at control points for diagnostic purposes. 
  if (3.952847075210474*rho[26]-3.061862178478972*rho[25]-3.061862178478972*rho[24]-3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]+2.371708245126284*rho[19]+2.371708245126284*rho[18]+2.371708245126284*rho[17]-1.369306393762915*rho[16]-1.369306393762915*rho[15]-1.369306393762915*rho[14]-1.369306393762915*rho[13]-1.369306393762915*rho[12]-1.369306393762915*rho[11]-1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]+1.060660171779821*rho[6]+1.060660171779821*rho[5]+1.060660171779821*rho[4]-0.6123724356957944*rho[3]-0.6123724356957944*rho[2]-0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]+3.061862178478972*rho[25]-3.061862178478972*rho[24]-3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]-2.371708245126284*rho[19]-2.371708245126284*rho[18]+2.371708245126284*rho[17]-1.369306393762915*rho[16]+1.369306393762915*rho[15]-1.369306393762915*rho[14]-1.369306393762915*rho[13]+1.369306393762915*rho[12]-1.369306393762915*rho[11]+1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]+1.060660171779821*rho[6]-1.060660171779821*rho[5]-1.060660171779821*rho[4]-0.6123724356957944*rho[3]-0.6123724356957944*rho[2]+0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]-3.061862178478972*rho[25]+3.061862178478972*rho[24]-3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]-2.371708245126284*rho[19]+2.371708245126284*rho[18]-2.371708245126284*rho[17]+1.369306393762915*rho[16]-1.369306393762915*rho[15]-1.369306393762915*rho[14]-1.369306393762915*rho[13]-1.369306393762915*rho[12]+1.369306393762915*rho[11]+1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]-1.060660171779821*rho[6]+1.060660171779821*rho[5]-1.060660171779821*rho[4]-0.6123724356957944*rho[3]+0.6123724356957944*rho[2]-0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]+3.061862178478972*rho[25]+3.061862178478972*rho[24]-3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]+2.371708245126284*rho[19]-2.371708245126284*rho[18]-2.371708245126284*rho[17]+1.369306393762915*rho[16]+1.369306393762915*rho[15]-1.369306393762915*rho[14]-1.369306393762915*rho[13]+1.369306393762915*rho[12]+1.369306393762915*rho[11]-1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]-1.060660171779821*rho[6]-1.060660171779821*rho[5]+1.060660171779821*rho[4]-0.6123724356957944*rho[3]+0.6123724356957944*rho[2]+0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]-3.061862178478972*rho[25]-3.061862178478972*rho[24]+3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]+2.371708245126284*rho[19]-2.371708245126284*rho[18]-2.371708245126284*rho[17]-1.369306393762915*rho[16]-1.369306393762915*rho[15]+1.369306393762915*rho[14]+1.369306393762915*rho[13]-1.369306393762915*rho[12]-1.369306393762915*rho[11]+1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]-1.060660171779821*rho[6]-1.060660171779821*rho[5]+1.060660171779821*rho[4]+0.6123724356957944*rho[3]-0.6123724356957944*rho[2]-0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]+3.061862178478972*rho[25]-3.061862178478972*rho[24]+3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]-2.371708245126284*rho[19]+2.371708245126284*rho[18]-2.371708245126284*rho[17]-1.369306393762915*rho[16]+1.369306393762915*rho[15]+1.369306393762915*rho[14]+1.369306393762915*rho[13]+1.369306393762915*rho[12]-1.369306393762915*rho[11]-1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]-1.060660171779821*rho[6]+1.060660171779821*rho[5]-1.060660171779821*rho[4]+0.6123724356957944*rho[3]-0.6123724356957944*rho[2]+0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]-3.061862178478972*rho[25]+3.061862178478972*rho[24]+3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]-2.371708245126284*rho[19]-2.371708245126284*rho[18]+2.371708245126284*rho[17]+1.369306393762915*rho[16]-1.369306393762915*rho[15]+1.369306393762915*rho[14]+1.369306393762915*rho[13]-1.369306393762915*rho[12]+1.369306393762915*rho[11]-1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]+1.060660171779821*rho[6]-1.060660171779821*rho[5]-1.060660171779821*rho[4]+0.6123724356957944*rho[3]+0.6123724356957944*rho[2]-0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
  if (3.952847075210474*rho[26]+3.061862178478972*rho[25]+3.061862178478972*rho[24]+3.061862178478972*rho[23]+1.767766952966368*rho[22]+1.767766952966368*rho[21]+1.767766952966368*rho[20]+2.371708245126284*rho[19]+2.371708245126284*rho[18]+2.371708245126284*rho[17]+1.369306393762915*rho[16]+1.369306393762915*rho[15]+1.369306393762915*rho[14]+1.369306393762915*rho[13]+1.369306393762915*rho[12]+1.369306393762915*rho[11]+1.837117307087383*rho[10]+0.7905694150420947*rho[9]+0.7905694150420947*rho[8]+0.7905694150420947*rho[7]+1.060660171779821*rho[6]+1.060660171779821*rho[5]+1.060660171779821*rho[4]+0.6123724356957944*rho[3]+0.6123724356957944*rho[2]+0.6123724356957944*rho[1]+0.3535533905932737*rho[0] < 0.0) cell_avg = 1; 
 
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
  gkyl_mat_set(&rhs_ux,4,0,rhoux[4]); 
  gkyl_mat_set(&rhs_uy,4,0,rhouy[4]); 
  gkyl_mat_set(&rhs_uz,4,0,rhouz[4]); 
  gkyl_mat_set(&rhs_ux,5,0,rhoux[5]); 
  gkyl_mat_set(&rhs_uy,5,0,rhouy[5]); 
  gkyl_mat_set(&rhs_uz,5,0,rhouz[5]); 
  gkyl_mat_set(&rhs_ux,6,0,rhoux[6]); 
  gkyl_mat_set(&rhs_uy,6,0,rhouy[6]); 
  gkyl_mat_set(&rhs_uz,6,0,rhouz[6]); 
  gkyl_mat_set(&rhs_ux,7,0,rhoux[7]); 
  gkyl_mat_set(&rhs_uy,7,0,rhouy[7]); 
  gkyl_mat_set(&rhs_uz,7,0,rhouz[7]); 
 
  double temp_rho = 0.0; 
  temp_rho = 0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,0,0,temp_rho); 
  gkyl_mat_set(&A_uy,0,0,temp_rho); 
  gkyl_mat_set(&A_uz,0,0,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,0,1,temp_rho); 
  gkyl_mat_set(&A_uy,0,1,temp_rho); 
  gkyl_mat_set(&A_uz,0,1,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,0,2,temp_rho); 
  gkyl_mat_set(&A_uy,0,2,temp_rho); 
  gkyl_mat_set(&A_uz,0,2,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,0,3,temp_rho); 
  gkyl_mat_set(&A_uy,0,3,temp_rho); 
  gkyl_mat_set(&A_uz,0,3,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,0,4,temp_rho); 
  gkyl_mat_set(&A_uy,0,4,temp_rho); 
  gkyl_mat_set(&A_uz,0,4,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,0,5,temp_rho); 
  gkyl_mat_set(&A_uy,0,5,temp_rho); 
  gkyl_mat_set(&A_uz,0,5,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,0,6,temp_rho); 
  gkyl_mat_set(&A_uy,0,6,temp_rho); 
  gkyl_mat_set(&A_uz,0,6,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,0,7,temp_rho); 
  gkyl_mat_set(&A_uy,0,7,temp_rho); 
  gkyl_mat_set(&A_uz,0,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,1,0,temp_rho); 
  gkyl_mat_set(&A_uy,1,0,temp_rho); 
  gkyl_mat_set(&A_uz,1,0,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[7]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,1,1,temp_rho); 
  gkyl_mat_set(&A_uy,1,1,temp_rho); 
  gkyl_mat_set(&A_uz,1,1,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,1,2,temp_rho); 
  gkyl_mat_set(&A_uy,1,2,temp_rho); 
  gkyl_mat_set(&A_uz,1,2,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,1,3,temp_rho); 
  gkyl_mat_set(&A_uy,1,3,temp_rho); 
  gkyl_mat_set(&A_uz,1,3,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[11]+0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,1,4,temp_rho); 
  gkyl_mat_set(&A_uy,1,4,temp_rho); 
  gkyl_mat_set(&A_uz,1,4,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[13]+0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,1,5,temp_rho); 
  gkyl_mat_set(&A_uy,1,5,temp_rho); 
  gkyl_mat_set(&A_uz,1,5,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,1,6,temp_rho); 
  gkyl_mat_set(&A_uy,1,6,temp_rho); 
  gkyl_mat_set(&A_uz,1,6,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[17]+0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,1,7,temp_rho); 
  gkyl_mat_set(&A_uy,1,7,temp_rho); 
  gkyl_mat_set(&A_uz,1,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,2,0,temp_rho); 
  gkyl_mat_set(&A_uy,2,0,temp_rho); 
  gkyl_mat_set(&A_uz,2,0,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,2,1,temp_rho); 
  gkyl_mat_set(&A_uy,2,1,temp_rho); 
  gkyl_mat_set(&A_uz,2,1,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[8]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,2,2,temp_rho); 
  gkyl_mat_set(&A_uy,2,2,temp_rho); 
  gkyl_mat_set(&A_uz,2,2,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,2,3,temp_rho); 
  gkyl_mat_set(&A_uy,2,3,temp_rho); 
  gkyl_mat_set(&A_uz,2,3,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[12]+0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,2,4,temp_rho); 
  gkyl_mat_set(&A_uy,2,4,temp_rho); 
  gkyl_mat_set(&A_uz,2,4,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,2,5,temp_rho); 
  gkyl_mat_set(&A_uy,2,5,temp_rho); 
  gkyl_mat_set(&A_uz,2,5,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[14]+0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,2,6,temp_rho); 
  gkyl_mat_set(&A_uy,2,6,temp_rho); 
  gkyl_mat_set(&A_uz,2,6,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[18]+0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,2,7,temp_rho); 
  gkyl_mat_set(&A_uy,2,7,temp_rho); 
  gkyl_mat_set(&A_uz,2,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,3,0,temp_rho); 
  gkyl_mat_set(&A_uy,3,0,temp_rho); 
  gkyl_mat_set(&A_uz,3,0,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,3,1,temp_rho); 
  gkyl_mat_set(&A_uy,3,1,temp_rho); 
  gkyl_mat_set(&A_uz,3,1,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,3,2,temp_rho); 
  gkyl_mat_set(&A_uy,3,2,temp_rho); 
  gkyl_mat_set(&A_uz,3,2,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[9]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,3,3,temp_rho); 
  gkyl_mat_set(&A_uy,3,3,temp_rho); 
  gkyl_mat_set(&A_uz,3,3,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,3,4,temp_rho); 
  gkyl_mat_set(&A_uy,3,4,temp_rho); 
  gkyl_mat_set(&A_uz,3,4,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[15]+0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,3,5,temp_rho); 
  gkyl_mat_set(&A_uy,3,5,temp_rho); 
  gkyl_mat_set(&A_uz,3,5,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[16]+0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,3,6,temp_rho); 
  gkyl_mat_set(&A_uy,3,6,temp_rho); 
  gkyl_mat_set(&A_uz,3,6,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[19]+0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,3,7,temp_rho); 
  gkyl_mat_set(&A_uy,3,7,temp_rho); 
  gkyl_mat_set(&A_uz,3,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,4,0,temp_rho); 
  gkyl_mat_set(&A_uy,4,0,temp_rho); 
  gkyl_mat_set(&A_uz,4,0,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[11]+0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,4,1,temp_rho); 
  gkyl_mat_set(&A_uy,4,1,temp_rho); 
  gkyl_mat_set(&A_uz,4,1,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[12]+0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,4,2,temp_rho); 
  gkyl_mat_set(&A_uy,4,2,temp_rho); 
  gkyl_mat_set(&A_uz,4,2,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,4,3,temp_rho); 
  gkyl_mat_set(&A_uy,4,3,temp_rho); 
  gkyl_mat_set(&A_uz,4,3,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[20]+0.3162277660168379*rho[8]+0.3162277660168379*rho[7]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,4,4,temp_rho); 
  gkyl_mat_set(&A_uy,4,4,temp_rho); 
  gkyl_mat_set(&A_uz,4,4,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[17]+0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,4,5,temp_rho); 
  gkyl_mat_set(&A_uy,4,5,temp_rho); 
  gkyl_mat_set(&A_uz,4,5,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[18]+0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,4,6,temp_rho); 
  gkyl_mat_set(&A_uy,4,6,temp_rho); 
  gkyl_mat_set(&A_uz,4,6,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[23]+0.3162277660168379*rho[14]+0.3162277660168379*rho[13]+0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,4,7,temp_rho); 
  gkyl_mat_set(&A_uy,4,7,temp_rho); 
  gkyl_mat_set(&A_uz,4,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,5,0,temp_rho); 
  gkyl_mat_set(&A_uy,5,0,temp_rho); 
  gkyl_mat_set(&A_uz,5,0,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[13]+0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,5,1,temp_rho); 
  gkyl_mat_set(&A_uy,5,1,temp_rho); 
  gkyl_mat_set(&A_uz,5,1,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,5,2,temp_rho); 
  gkyl_mat_set(&A_uy,5,2,temp_rho); 
  gkyl_mat_set(&A_uz,5,2,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[15]+0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,5,3,temp_rho); 
  gkyl_mat_set(&A_uy,5,3,temp_rho); 
  gkyl_mat_set(&A_uz,5,3,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[17]+0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,5,4,temp_rho); 
  gkyl_mat_set(&A_uy,5,4,temp_rho); 
  gkyl_mat_set(&A_uz,5,4,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[21]+0.3162277660168379*rho[9]+0.3162277660168379*rho[7]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,5,5,temp_rho); 
  gkyl_mat_set(&A_uy,5,5,temp_rho); 
  gkyl_mat_set(&A_uz,5,5,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[19]+0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,5,6,temp_rho); 
  gkyl_mat_set(&A_uy,5,6,temp_rho); 
  gkyl_mat_set(&A_uz,5,6,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[24]+0.3162277660168379*rho[16]+0.3162277660168379*rho[11]+0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,5,7,temp_rho); 
  gkyl_mat_set(&A_uy,5,7,temp_rho); 
  gkyl_mat_set(&A_uz,5,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,6,0,temp_rho); 
  gkyl_mat_set(&A_uy,6,0,temp_rho); 
  gkyl_mat_set(&A_uz,6,0,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,6,1,temp_rho); 
  gkyl_mat_set(&A_uy,6,1,temp_rho); 
  gkyl_mat_set(&A_uz,6,1,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[14]+0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,6,2,temp_rho); 
  gkyl_mat_set(&A_uy,6,2,temp_rho); 
  gkyl_mat_set(&A_uz,6,2,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[16]+0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,6,3,temp_rho); 
  gkyl_mat_set(&A_uy,6,3,temp_rho); 
  gkyl_mat_set(&A_uz,6,3,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[18]+0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,6,4,temp_rho); 
  gkyl_mat_set(&A_uy,6,4,temp_rho); 
  gkyl_mat_set(&A_uz,6,4,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[19]+0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,6,5,temp_rho); 
  gkyl_mat_set(&A_uy,6,5,temp_rho); 
  gkyl_mat_set(&A_uz,6,5,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[22]+0.3162277660168379*rho[9]+0.3162277660168379*rho[8]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,6,6,temp_rho); 
  gkyl_mat_set(&A_uy,6,6,temp_rho); 
  gkyl_mat_set(&A_uz,6,6,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[25]+0.3162277660168379*rho[15]+0.3162277660168379*rho[12]+0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,6,7,temp_rho); 
  gkyl_mat_set(&A_uy,6,7,temp_rho); 
  gkyl_mat_set(&A_uz,6,7,temp_rho); 
 
  temp_rho = 0.3535533905932737*rho[10]; 
  gkyl_mat_set(&A_ux,7,0,temp_rho); 
  gkyl_mat_set(&A_uy,7,0,temp_rho); 
  gkyl_mat_set(&A_uz,7,0,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[17]+0.3535533905932737*rho[6]; 
  gkyl_mat_set(&A_ux,7,1,temp_rho); 
  gkyl_mat_set(&A_uy,7,1,temp_rho); 
  gkyl_mat_set(&A_uz,7,1,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[18]+0.3535533905932737*rho[5]; 
  gkyl_mat_set(&A_ux,7,2,temp_rho); 
  gkyl_mat_set(&A_uy,7,2,temp_rho); 
  gkyl_mat_set(&A_uz,7,2,temp_rho); 
 
  temp_rho = 0.3162277660168379*rho[19]+0.3535533905932737*rho[4]; 
  gkyl_mat_set(&A_ux,7,3,temp_rho); 
  gkyl_mat_set(&A_uy,7,3,temp_rho); 
  gkyl_mat_set(&A_uz,7,3,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[23]+0.3162277660168379*rho[14]+0.3162277660168379*rho[13]+0.3535533905932737*rho[3]; 
  gkyl_mat_set(&A_ux,7,4,temp_rho); 
  gkyl_mat_set(&A_uy,7,4,temp_rho); 
  gkyl_mat_set(&A_uz,7,4,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[24]+0.3162277660168379*rho[16]+0.3162277660168379*rho[11]+0.3535533905932737*rho[2]; 
  gkyl_mat_set(&A_ux,7,5,temp_rho); 
  gkyl_mat_set(&A_uy,7,5,temp_rho); 
  gkyl_mat_set(&A_uz,7,5,temp_rho); 
 
  temp_rho = 0.2828427124746191*rho[25]+0.3162277660168379*rho[15]+0.3162277660168379*rho[12]+0.3535533905932737*rho[1]; 
  gkyl_mat_set(&A_ux,7,6,temp_rho); 
  gkyl_mat_set(&A_uy,7,6,temp_rho); 
  gkyl_mat_set(&A_uz,7,6,temp_rho); 
 
  temp_rho = 0.2529822128134704*rho[26]+0.2828427124746191*rho[22]+0.2828427124746191*rho[21]+0.2828427124746191*rho[20]+0.3162277660168379*rho[9]+0.3162277660168379*rho[8]+0.3162277660168379*rho[7]+0.3535533905932737*rho[0]; 
  gkyl_mat_set(&A_ux,7,7,temp_rho); 
  gkyl_mat_set(&A_uy,7,7,temp_rho); 
  gkyl_mat_set(&A_uz,7,7,temp_rho); 
 
  return cell_avg;
} 
