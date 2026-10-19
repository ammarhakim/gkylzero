#include <gkyl_flash_coupling_kernels.h> 
void flash_coupling_mom_1x1v_ser_p1(const double *dxv, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
  const double *pkpm_div_ppar, double *flux_div_p) 
{ 
  // dxv: Cell spacing 
  // vlasov_pkpm_moms_l/c/r: Moments of the distribution functions, [rho, p_par, p_perp] 
  // pkpm_div_ppar: Divergence of pressure tensor, computed from kinetic flux for consistent pressure force. 
  // flux_div_p: Output volume expansion of momentum flux contributions, div(P). 

  const double dx1 = 2.0/dxv[0]; 

  const double *p_par_l = &vlasov_pkpm_moms_l[2]; 
  const double *p_par_c = &vlasov_pkpm_moms_c[2]; 
  const double *p_par_r = &vlasov_pkpm_moms_r[2]; 
  const double *p_perp_l = &vlasov_pkpm_moms_l[4]; 
  const double *p_perp_c = &vlasov_pkpm_moms_c[4]; 
  const double *p_perp_r = &vlasov_pkpm_moms_r[4]; 

  double p_l[2]; 
  double p_c[2]; 
  double p_r[2]; 
  p_l[0] = 0.6666666666666666*p_perp_l[0]+0.3333333333333333*p_par_l[0]; 
  p_l[1] = 0.6666666666666666*p_perp_l[1]+0.3333333333333333*p_par_l[1]; 

  p_c[0] = 0.6666666666666666*p_perp_c[0]+0.3333333333333333*p_par_c[0]; 
  p_c[1] = 0.6666666666666666*p_perp_c[1]+0.3333333333333333*p_par_c[1]; 

  p_r[0] = 0.6666666666666666*p_perp_r[0]+0.3333333333333333*p_par_r[0]; 
  p_r[1] = 0.6666666666666666*p_perp_r[1]+0.3333333333333333*p_par_r[1]; 

  double p_surf_lr[1]; 
  double p_surf_cl[1]; 
  double p_surf_cr[1]; 
  double p_surf_rl[1]; 

  p_surf_lr[0] = 1.224744871391589*p_l[1]+0.7071067811865475*p_l[0]; 
  p_surf_cl[0] = 0.7071067811865475*p_c[0]-1.224744871391589*p_c[1]; 
  p_surf_cr[0] = 1.224744871391589*p_c[1]+0.7071067811865475*p_c[0]; 
  p_surf_rl[0] = 0.7071067811865475*p_r[0]-1.224744871391589*p_r[1]; 

  double gradx_p[2]; 
  gradx_p[0] = (0.3535533905932737*p_surf_rl[0]-0.3535533905932737*p_surf_lr[0]+0.3535533905932737*p_surf_cr[0]-0.3535533905932737*p_surf_cl[0])*dx1; 
  gradx_p[1] = (0.6123724356957944*(p_surf_rl[0]+p_surf_lr[0]+p_surf_cr[0]+p_surf_cl[0])-1.7320508075688772*p_c[0])*dx1; 

  flux_div_p[0] = pkpm_div_ppar[0]-1.0*gradx_p[0]; 
  flux_div_p[1] = pkpm_div_ppar[1]-1.0*gradx_p[1]; 
}
