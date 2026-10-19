#include <gkyl_flash_coupling_kernels.h> 
void flash_coupling_ener_1x1v_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *pkpm_prim, 
  const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r,
  double *flux_p_gradu) 
{ 
  // dxv: Cell spacing 
  // vlasov_pkpm_moms: Moments of the distribution functions, [rho, p_par, p_perp] 
  // pkpm_prim: Volume expansions of primitive variables: 
  //            [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m]. 
  // pkpm_prim_surf_l/c/r: Surface expansion of primitive variables [u_i, 3*T_ii/m] 
  // flux_p_gradu: Output volume expansion of P:grad(u) for contribution to energy flux. 

  const double dx1 = 2.0/dxv[0]; 

  const double *ux = &pkpm_prim[0]; 
  const double *p_par = &vlasov_pkpm_moms[2]; 
  const double *p_perp = &vlasov_pkpm_moms[4]; 

  double pi_xx[2]; 
  pi_xx[0] = 0.6666666666666666*p_par[0]-0.6666666666666666*p_perp[0]; 
  pi_xx[1] = 0.6666666666666666*p_par[1]-0.6666666666666666*p_perp[1]; 
  const double *ux_surf_lr = &pkpm_prim_surf_l[1]; 
  const double *ux_surf_cl = &pkpm_prim_surf_c[0]; 
  const double *ux_surf_cr = &pkpm_prim_surf_c[1]; 
  const double *ux_surf_rl = &pkpm_prim_surf_r[0]; 

  double grad_u_x[2] = {0.0}; 
  grad_u_x[0] = (0.3535533905932737*ux_surf_rl[0]-0.3535533905932737*ux_surf_lr[0]+0.3535533905932737*ux_surf_cr[0]-0.3535533905932737*ux_surf_cl[0])*dx1; 
  grad_u_x[1] = (0.6123724356957944*(ux_surf_rl[0]+ux_surf_lr[0]+ux_surf_cr[0]+ux_surf_cl[0])-1.7320508075688772*ux[0])*dx1; 

  flux_p_gradu[0] = 0.7071067811865475*grad_u_x[1]*pi_xx[1]+0.7071067811865475*grad_u_x[0]*pi_xx[0]; 
  flux_p_gradu[1] = 0.7071067811865475*grad_u_x[0]*pi_xx[1]+0.7071067811865475*pi_xx[0]*grad_u_x[1]; 
}
