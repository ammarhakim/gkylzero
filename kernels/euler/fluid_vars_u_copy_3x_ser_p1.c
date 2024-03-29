#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_u_copy_3x_ser_p1(int count, struct gkyl_nmat *x, 
    double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf) 
{ 
  // count:  integer to indicate which matrix being fetched. 
  // x:      Input solution vector. 
  // u:      Output volume expansion of flow velocity: 
  //         [ux, uy, uz]. 
  // u_surf: Output surface expansion of flow velocity 
  //         [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //          ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr,  
  //          ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
 
  struct gkyl_mat x_ux = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(x, count+2); 
  double *ux = &u[0]; 
  double *uy = &u[8]; 
  double *uz = &u[16]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 
  ux[4] = gkyl_mat_get(&x_ux,4,0); 
  uy[4] = gkyl_mat_get(&x_uy,4,0); 
  uz[4] = gkyl_mat_get(&x_uz,4,0); 
  ux[5] = gkyl_mat_get(&x_ux,5,0); 
  uy[5] = gkyl_mat_get(&x_uy,5,0); 
  uz[5] = gkyl_mat_get(&x_uz,5,0); 
  ux[6] = gkyl_mat_get(&x_ux,6,0); 
  uy[6] = gkyl_mat_get(&x_uy,6,0); 
  uz[6] = gkyl_mat_get(&x_uz,6,0); 
  ux[7] = gkyl_mat_get(&x_ux,7,0); 
  uy[7] = gkyl_mat_get(&x_uy,7,0); 
  uz[7] = gkyl_mat_get(&x_uz,7,0); 

  double *ux_xl = &u_surf[0]; 
  double *ux_xr = &u_surf[4]; 
  double *uy_xl = &u_surf[8]; 
  double *uy_xr = &u_surf[12]; 
  double *uz_xl = &u_surf[16]; 
  double *uz_xr = &u_surf[20]; 
 
  ux_xl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[1]; 
  ux_xl[1] = 0.7071067811865475*ux[2]-1.224744871391589*ux[4]; 
  ux_xl[2] = 0.7071067811865475*ux[3]-1.224744871391589*ux[5]; 
  ux_xl[3] = 0.7071067811865475*ux[6]-1.224744871391589*ux[7]; 
  uy_xl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[1]; 
  uy_xl[1] = 0.7071067811865475*uy[2]-1.224744871391589*uy[4]; 
  uy_xl[2] = 0.7071067811865475*uy[3]-1.224744871391589*uy[5]; 
  uy_xl[3] = 0.7071067811865475*uy[6]-1.224744871391589*uy[7]; 
  uz_xl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[1]; 
  uz_xl[1] = 0.7071067811865475*uz[2]-1.224744871391589*uz[4]; 
  uz_xl[2] = 0.7071067811865475*uz[3]-1.224744871391589*uz[5]; 
  uz_xl[3] = 0.7071067811865475*uz[6]-1.224744871391589*uz[7]; 
 
  ux_xr[0] = 1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  ux_xr[1] = 1.224744871391589*ux[4]+0.7071067811865475*ux[2]; 
  ux_xr[2] = 1.224744871391589*ux[5]+0.7071067811865475*ux[3]; 
  ux_xr[3] = 1.224744871391589*ux[7]+0.7071067811865475*ux[6]; 
  uy_xr[0] = 1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uy_xr[1] = 1.224744871391589*uy[4]+0.7071067811865475*uy[2]; 
  uy_xr[2] = 1.224744871391589*uy[5]+0.7071067811865475*uy[3]; 
  uy_xr[3] = 1.224744871391589*uy[7]+0.7071067811865475*uy[6]; 
  uz_xr[0] = 1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  uz_xr[1] = 1.224744871391589*uz[4]+0.7071067811865475*uz[2]; 
  uz_xr[2] = 1.224744871391589*uz[5]+0.7071067811865475*uz[3]; 
  uz_xr[3] = 1.224744871391589*uz[7]+0.7071067811865475*uz[6]; 
 
  double *ux_yl = &u_surf[24]; 
  double *ux_yr = &u_surf[28]; 
  double *uy_yl = &u_surf[32]; 
  double *uy_yr = &u_surf[36]; 
  double *uz_yl = &u_surf[40]; 
  double *uz_yr = &u_surf[44]; 
 
  ux_yl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[2]; 
  ux_yl[1] = 0.7071067811865475*ux[1]-1.224744871391589*ux[4]; 
  ux_yl[2] = 0.7071067811865475*ux[3]-1.224744871391589*ux[6]; 
  ux_yl[3] = 0.7071067811865475*ux[5]-1.224744871391589*ux[7]; 
  uy_yl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[2]; 
  uy_yl[1] = 0.7071067811865475*uy[1]-1.224744871391589*uy[4]; 
  uy_yl[2] = 0.7071067811865475*uy[3]-1.224744871391589*uy[6]; 
  uy_yl[3] = 0.7071067811865475*uy[5]-1.224744871391589*uy[7]; 
  uz_yl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[2]; 
  uz_yl[1] = 0.7071067811865475*uz[1]-1.224744871391589*uz[4]; 
  uz_yl[2] = 0.7071067811865475*uz[3]-1.224744871391589*uz[6]; 
  uz_yl[3] = 0.7071067811865475*uz[5]-1.224744871391589*uz[7]; 
 
  ux_yr[0] = 1.224744871391589*ux[2]+0.7071067811865475*ux[0]; 
  ux_yr[1] = 1.224744871391589*ux[4]+0.7071067811865475*ux[1]; 
  ux_yr[2] = 1.224744871391589*ux[6]+0.7071067811865475*ux[3]; 
  ux_yr[3] = 1.224744871391589*ux[7]+0.7071067811865475*ux[5]; 
  uy_yr[0] = 1.224744871391589*uy[2]+0.7071067811865475*uy[0]; 
  uy_yr[1] = 1.224744871391589*uy[4]+0.7071067811865475*uy[1]; 
  uy_yr[2] = 1.224744871391589*uy[6]+0.7071067811865475*uy[3]; 
  uy_yr[3] = 1.224744871391589*uy[7]+0.7071067811865475*uy[5]; 
  uz_yr[0] = 1.224744871391589*uz[2]+0.7071067811865475*uz[0]; 
  uz_yr[1] = 1.224744871391589*uz[4]+0.7071067811865475*uz[1]; 
  uz_yr[2] = 1.224744871391589*uz[6]+0.7071067811865475*uz[3]; 
  uz_yr[3] = 1.224744871391589*uz[7]+0.7071067811865475*uz[5]; 
 
  double *ux_zl = &u_surf[48]; 
  double *ux_zr = &u_surf[52]; 
  double *uy_zl = &u_surf[56]; 
  double *uy_zr = &u_surf[60]; 
  double *uz_zl = &u_surf[64]; 
  double *uz_zr = &u_surf[68]; 
 
  ux_zl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[3]; 
  ux_zl[1] = 0.7071067811865475*ux[1]-1.224744871391589*ux[5]; 
  ux_zl[2] = 0.7071067811865475*ux[2]-1.224744871391589*ux[6]; 
  ux_zl[3] = 0.7071067811865475*ux[4]-1.224744871391589*ux[7]; 
  uy_zl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[3]; 
  uy_zl[1] = 0.7071067811865475*uy[1]-1.224744871391589*uy[5]; 
  uy_zl[2] = 0.7071067811865475*uy[2]-1.224744871391589*uy[6]; 
  uy_zl[3] = 0.7071067811865475*uy[4]-1.224744871391589*uy[7]; 
  uz_zl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[3]; 
  uz_zl[1] = 0.7071067811865475*uz[1]-1.224744871391589*uz[5]; 
  uz_zl[2] = 0.7071067811865475*uz[2]-1.224744871391589*uz[6]; 
  uz_zl[3] = 0.7071067811865475*uz[4]-1.224744871391589*uz[7]; 
 
  ux_zr[0] = 1.224744871391589*ux[3]+0.7071067811865475*ux[0]; 
  ux_zr[1] = 1.224744871391589*ux[5]+0.7071067811865475*ux[1]; 
  ux_zr[2] = 1.224744871391589*ux[6]+0.7071067811865475*ux[2]; 
  ux_zr[3] = 1.224744871391589*ux[7]+0.7071067811865475*ux[4]; 
  uy_zr[0] = 1.224744871391589*uy[3]+0.7071067811865475*uy[0]; 
  uy_zr[1] = 1.224744871391589*uy[5]+0.7071067811865475*uy[1]; 
  uy_zr[2] = 1.224744871391589*uy[6]+0.7071067811865475*uy[2]; 
  uy_zr[3] = 1.224744871391589*uy[7]+0.7071067811865475*uy[4]; 
  uz_zr[0] = 1.224744871391589*uz[3]+0.7071067811865475*uz[0]; 
  uz_zr[1] = 1.224744871391589*uz[5]+0.7071067811865475*uz[1]; 
  uz_zr[2] = 1.224744871391589*uz[6]+0.7071067811865475*uz[2]; 
  uz_zr[3] = 1.224744871391589*uz[7]+0.7071067811865475*uz[4]; 
 
} 
 
