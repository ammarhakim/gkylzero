#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_u_surf_1x_tensor_p1(const double *pkpm_u, double* GKYL_RESTRICT pkpm_u_surf) 
{ 
  // pkpm_u:  Input volume expansion of [ux, uy, uz]. 
  // pkpm_u_surf: Output surface expansion of flow velocity. 
  //              [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //               ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr, 
  //               ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
 
  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[2]; 
  const double *uz = &pkpm_u[4]; 
  double *ux_xl = &pkpm_u_surf[0]; 
  double *ux_xr = &pkpm_u_surf[1]; 
  double *uy_xl = &pkpm_u_surf[2]; 
  double *uy_xr = &pkpm_u_surf[3]; 
  double *uz_xl = &pkpm_u_surf[4]; 
  double *uz_xr = &pkpm_u_surf[5]; 
 
  ux_xl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[1]; 
  ux_xr[0] = 1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  uy_xl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[1]; 
  uy_xr[0] = 1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uz_xl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[1]; 
  uz_xr[0] = 1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
 
} 
 
