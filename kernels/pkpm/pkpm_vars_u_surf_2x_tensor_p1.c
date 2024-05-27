#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_u_surf_2x_tensor_p1(const double *pkpm_u, double* GKYL_RESTRICT pkpm_u_surf) 
{ 
  // pkpm_u:  Input volume expansion of [ux, uy, uz]. 
  // pkpm_u_surf: Output surface expansion of flow velocity. 
  //              [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //               ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
  //               ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
 
  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[4]; 
  const double *uz = &pkpm_u[8]; 
  double *ux_xl = &pkpm_u_surf[0]; 
  double *ux_xr = &pkpm_u_surf[2]; 
  double *uy_xl = &pkpm_u_surf[4]; 
  double *uy_xr = &pkpm_u_surf[6]; 
  double *uz_xl = &pkpm_u_surf[8]; 
  double *uz_xr = &pkpm_u_surf[10]; 
 
  ux_xl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[1]; 
  ux_xl[1] = 0.7071067811865475*ux[2]-1.224744871391589*ux[3]; 
  uy_xl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[1]; 
  uy_xl[1] = 0.7071067811865475*uy[2]-1.224744871391589*uy[3]; 
  uz_xl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[1]; 
  uz_xl[1] = 0.7071067811865475*uz[2]-1.224744871391589*uz[3]; 
 
  ux_xr[0] = 1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  ux_xr[1] = 1.224744871391589*ux[3]+0.7071067811865475*ux[2]; 
  uy_xr[0] = 1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uy_xr[1] = 1.224744871391589*uy[3]+0.7071067811865475*uy[2]; 
  uz_xr[0] = 1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  uz_xr[1] = 1.224744871391589*uz[3]+0.7071067811865475*uz[2]; 
 
  double *ux_yl = &pkpm_u_surf[12]; 
  double *ux_yr = &pkpm_u_surf[14]; 
  double *uy_yl = &pkpm_u_surf[16]; 
  double *uy_yr = &pkpm_u_surf[18]; 
  double *uz_yl = &pkpm_u_surf[20]; 
  double *uz_yr = &pkpm_u_surf[22]; 
 
  ux_yl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[2]; 
  ux_yl[1] = 0.7071067811865475*ux[1]-1.224744871391589*ux[3]; 
  uy_yl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[2]; 
  uy_yl[1] = 0.7071067811865475*uy[1]-1.224744871391589*uy[3]; 
  uz_yl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[2]; 
  uz_yl[1] = 0.7071067811865475*uz[1]-1.224744871391589*uz[3]; 
 
  ux_yr[0] = 1.224744871391589*ux[2]+0.7071067811865475*ux[0]; 
  ux_yr[1] = 1.224744871391589*ux[3]+0.7071067811865475*ux[1]; 
  uy_yr[0] = 1.224744871391589*uy[2]+0.7071067811865475*uy[0]; 
  uy_yr[1] = 1.224744871391589*uy[3]+0.7071067811865475*uy[1]; 
  uz_yr[0] = 1.224744871391589*uz[2]+0.7071067811865475*uz[0]; 
  uz_yr[1] = 1.224744871391589*uz[3]+0.7071067811865475*uz[1]; 
 
} 
 
