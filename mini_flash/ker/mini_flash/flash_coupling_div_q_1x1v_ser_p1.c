#include <gkyl_flash_coupling_kernels.h> 
void flash_coupling_div_q_1x1v_ser_p1(const double *w, const double *dxv, double mass, 
  const double *fl, const double *fc, const double *fr,
  const double *bvar_c, const double *pkpm_max_b,
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  double *flux_div_q) 
{ 
  // dxv: Cell spacing 

  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double wvpar_sq = wvpar*wvpar; 
  const double dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar; 
  const double dvpar_cu = dvpar*dvpar*dvpar; 
  const double volFact = dxv[1]/2.0; 

  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  const double *G_1l = &fl[6]; 
  const double *G_1c = &fc[6]; 
  const double *G_1r = &fr[6]; 

  double F_0_lr[3] = {0.0}; 
  double F_0_cl[3] = {0.0}; 
  double F_0_cr[3] = {0.0}; 
  double F_0_rl[3] = {0.0}; 
  double G_1_lr[3] = {0.0}; 
  double G_1_cl[3] = {0.0}; 
  double G_1_cr[3] = {0.0}; 
  double G_1_rl[3] = {0.0}; 

  double Ghat_F_0_vpar_l[3] = {0.0}; 
  double Ghat_F_0_vpar_r[3] = {0.0}; 
  double Ghat_G_1_vpar_l[3] = {0.0}; 
  double Ghat_G_1_vpar_r[3] = {0.0}; 

  F_0_lr[0] = 1.224744871391589*F_0l[1]+0.7071067811865475*F_0l[0]; 
  F_0_lr[1] = 1.224744871391589*F_0l[3]+0.7071067811865475*F_0l[2]; 
  F_0_lr[2] = 1.224744871391589*F_0l[5]+0.7071067811865475*F_0l[4]; 

  F_0_cl[0] = 0.7071067811865475*F_0c[0]-1.224744871391589*F_0c[1]; 
  F_0_cl[1] = 0.7071067811865475*F_0c[2]-1.224744871391589*F_0c[3]; 
  F_0_cl[2] = 0.7071067811865475*F_0c[4]-1.224744871391589*F_0c[5]; 

  F_0_cr[0] = 1.224744871391589*F_0c[1]+0.7071067811865475*F_0c[0]; 
  F_0_cr[1] = 1.224744871391589*F_0c[3]+0.7071067811865475*F_0c[2]; 
  F_0_cr[2] = 1.224744871391589*F_0c[5]+0.7071067811865475*F_0c[4]; 

  F_0_rl[0] = 0.7071067811865475*F_0r[0]-1.224744871391589*F_0r[1]; 
  F_0_rl[1] = 0.7071067811865475*F_0r[2]-1.224744871391589*F_0r[3]; 
  F_0_rl[2] = 0.7071067811865475*F_0r[4]-1.224744871391589*F_0r[5]; 


  G_1_lr[0] = 1.224744871391589*G_1l[1]+0.7071067811865475*G_1l[0]; 
  G_1_lr[1] = 1.224744871391589*G_1l[3]+0.7071067811865475*G_1l[2]; 
  G_1_lr[2] = 1.224744871391589*G_1l[5]+0.7071067811865475*G_1l[4]; 

  G_1_cl[0] = 0.7071067811865475*G_1c[0]-1.224744871391589*G_1c[1]; 
  G_1_cl[1] = 0.7071067811865475*G_1c[2]-1.224744871391589*G_1c[3]; 
  G_1_cl[2] = 0.7071067811865475*G_1c[4]-1.224744871391589*G_1c[5]; 

  G_1_cr[0] = 1.224744871391589*G_1c[1]+0.7071067811865475*G_1c[0]; 
  G_1_cr[1] = 1.224744871391589*G_1c[3]+0.7071067811865475*G_1c[2]; 
  G_1_cr[2] = 1.224744871391589*G_1c[5]+0.7071067811865475*G_1c[4]; 

  G_1_rl[0] = 0.7071067811865475*G_1r[0]-1.224744871391589*G_1r[1]; 
  G_1_rl[1] = 0.7071067811865475*G_1r[2]-1.224744871391589*G_1r[3]; 
  G_1_rl[2] = 0.7071067811865475*G_1r[4]-1.224744871391589*G_1r[5]; 

  const double *b_surf_lr = &bvar_surf_l[1]; 
  const double *b_surf_cl = &bvar_surf_c[0]; 
  const double *b_surf_cr = &bvar_surf_c[1]; 
  const double *b_surf_rl = &bvar_surf_r[0]; 

  const double *pkpm_max_b_l = &pkpm_max_b[0]; 
  const double *pkpm_max_b_r = &pkpm_max_b[1]; 

  const double *b_c = &bvar_c[0]; 
  double alpha_c[6] = {0.0}; 
  alpha_c[0] = 1.4142135623730951*b_c[0]*wvpar; 
  alpha_c[1] = 1.4142135623730951*b_c[1]*wvpar; 
  alpha_c[2] = 0.408248290463863*b_c[0]*dvpar; 
  alpha_c[3] = 0.408248290463863*b_c[1]*dvpar; 

  double bl_r = b_surf_lr[0]; 
  double bc_l = b_surf_cl[0]; 
  double bc_r = b_surf_cr[0]; 
  double br_l = b_surf_rl[0]; 
  double avg_b_l = 0.5*(bl_r + bc_l); 
  double avg_b_r = 0.5*(bc_r + br_l); 

  double max_b_l = pkpm_max_b_l[0]; 
  double max_b_r = pkpm_max_b_r[0]; 

  if (wvpar>0) { 

  Ghat_F_0_vpar_l[0] = ((0.5*F_0_lr[0]-0.5*F_0_cl[0])*max_b_l+0.5*(F_0_lr[0]+F_0_cl[0])*avg_b_l)*wvpar+dvpar*((0.14433756729740646*F_0_lr[1]-0.14433756729740646*F_0_cl[1])*max_b_l+0.14433756729740646*(F_0_lr[1]+F_0_cl[1])*avg_b_l); 
  Ghat_F_0_vpar_l[1] = ((0.5*F_0_lr[1]-0.5*F_0_cl[1])*max_b_l+0.5*(F_0_lr[1]+F_0_cl[1])*avg_b_l)*wvpar+dvpar*((0.12909944487358058*F_0_lr[2]-0.12909944487358058*F_0_cl[2]+0.14433756729740646*F_0_lr[0]-0.14433756729740646*F_0_cl[0])*max_b_l+(0.12909944487358058*(F_0_lr[2]+F_0_cl[2])+0.14433756729740646*(F_0_lr[0]+F_0_cl[0]))*avg_b_l); 
  Ghat_F_0_vpar_l[2] = ((0.5*F_0_lr[2]-0.5*F_0_cl[2])*max_b_l+0.5*(F_0_lr[2]+F_0_cl[2])*avg_b_l)*wvpar+dvpar*((0.12909944487358058*F_0_lr[1]-0.12909944487358058*F_0_cl[1])*max_b_l+0.12909944487358058*(F_0_lr[1]+F_0_cl[1])*avg_b_l); 

  Ghat_F_0_vpar_r[0] = (0.5*(F_0_cr[0]*max_b_r+(F_0_rl[0]+F_0_cr[0])*avg_b_r)-0.5*F_0_rl[0]*max_b_r)*wvpar+dvpar*(0.14433756729740646*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.14433756729740646*F_0_rl[1]*max_b_r); 
  Ghat_F_0_vpar_r[1] = (0.5*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.5*F_0_rl[1]*max_b_r)*wvpar+dvpar*((-(0.12909944487358058*F_0_rl[2])+0.12909944487358058*F_0_cr[2]-0.14433756729740646*F_0_rl[0]+0.14433756729740646*F_0_cr[0])*max_b_r+(0.12909944487358058*(F_0_rl[2]+F_0_cr[2])+0.14433756729740646*(F_0_rl[0]+F_0_cr[0]))*avg_b_r); 
  Ghat_F_0_vpar_r[2] = (0.5*(F_0_cr[2]*max_b_r+(F_0_rl[2]+F_0_cr[2])*avg_b_r)-0.5*F_0_rl[2]*max_b_r)*wvpar+dvpar*(0.12909944487358058*(F_0_cr[1]*max_b_r+(F_0_rl[1]+F_0_cr[1])*avg_b_r)-0.12909944487358058*F_0_rl[1]*max_b_r); 

  Ghat_G_1_vpar_l[0] = ((0.5*G_1_lr[0]-0.5*G_1_cl[0])*max_b_l+0.5*(G_1_lr[0]+G_1_cl[0])*avg_b_l)*wvpar+dvpar*((0.14433756729740646*G_1_lr[1]-0.14433756729740646*G_1_cl[1])*max_b_l+0.14433756729740646*(G_1_lr[1]+G_1_cl[1])*avg_b_l); 
  Ghat_G_1_vpar_l[1] = ((0.5*G_1_lr[1]-0.5*G_1_cl[1])*max_b_l+0.5*(G_1_lr[1]+G_1_cl[1])*avg_b_l)*wvpar+dvpar*((0.12909944487358058*G_1_lr[2]-0.12909944487358058*G_1_cl[2]+0.14433756729740646*G_1_lr[0]-0.14433756729740646*G_1_cl[0])*max_b_l+(0.12909944487358058*(G_1_lr[2]+G_1_cl[2])+0.14433756729740646*(G_1_lr[0]+G_1_cl[0]))*avg_b_l); 
  Ghat_G_1_vpar_l[2] = ((0.5*G_1_lr[2]-0.5*G_1_cl[2])*max_b_l+0.5*(G_1_lr[2]+G_1_cl[2])*avg_b_l)*wvpar+dvpar*((0.12909944487358058*G_1_lr[1]-0.12909944487358058*G_1_cl[1])*max_b_l+0.12909944487358058*(G_1_lr[1]+G_1_cl[1])*avg_b_l); 

  Ghat_G_1_vpar_r[0] = (0.5*(G_1_cr[0]*max_b_r+(G_1_rl[0]+G_1_cr[0])*avg_b_r)-0.5*G_1_rl[0]*max_b_r)*wvpar+dvpar*(0.14433756729740646*(G_1_cr[1]*max_b_r+(G_1_rl[1]+G_1_cr[1])*avg_b_r)-0.14433756729740646*G_1_rl[1]*max_b_r); 
  Ghat_G_1_vpar_r[1] = (0.5*(G_1_cr[1]*max_b_r+(G_1_rl[1]+G_1_cr[1])*avg_b_r)-0.5*G_1_rl[1]*max_b_r)*wvpar+dvpar*((-(0.12909944487358058*G_1_rl[2])+0.12909944487358058*G_1_cr[2]-0.14433756729740646*G_1_rl[0]+0.14433756729740646*G_1_cr[0])*max_b_r+(0.12909944487358058*(G_1_rl[2]+G_1_cr[2])+0.14433756729740646*(G_1_rl[0]+G_1_cr[0]))*avg_b_r); 
  Ghat_G_1_vpar_r[2] = (0.5*(G_1_cr[2]*max_b_r+(G_1_rl[2]+G_1_cr[2])*avg_b_r)-0.5*G_1_rl[2]*max_b_r)*wvpar+dvpar*(0.12909944487358058*(G_1_cr[1]*max_b_r+(G_1_rl[1]+G_1_cr[1])*avg_b_r)-0.12909944487358058*G_1_rl[1]*max_b_r); 

  } else { 

  Ghat_F_0_vpar_l[0] = (0.5*(F_0_cl[0]*max_b_l+(F_0_lr[0]+F_0_cl[0])*avg_b_l)-0.5*F_0_lr[0]*max_b_l)*wvpar+dvpar*(0.14433756729740646*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.14433756729740646*F_0_lr[1]*max_b_l); 
  Ghat_F_0_vpar_l[1] = (0.5*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.5*F_0_lr[1]*max_b_l)*wvpar+dvpar*((-(0.12909944487358058*F_0_lr[2])+0.12909944487358058*F_0_cl[2]-0.14433756729740646*F_0_lr[0]+0.14433756729740646*F_0_cl[0])*max_b_l+(0.12909944487358058*(F_0_lr[2]+F_0_cl[2])+0.14433756729740646*(F_0_lr[0]+F_0_cl[0]))*avg_b_l); 
  Ghat_F_0_vpar_l[2] = (0.5*(F_0_cl[2]*max_b_l+(F_0_lr[2]+F_0_cl[2])*avg_b_l)-0.5*F_0_lr[2]*max_b_l)*wvpar+dvpar*(0.12909944487358058*(F_0_cl[1]*max_b_l+(F_0_lr[1]+F_0_cl[1])*avg_b_l)-0.12909944487358058*F_0_lr[1]*max_b_l); 

  Ghat_F_0_vpar_r[0] = ((0.5*F_0_rl[0]-0.5*F_0_cr[0])*max_b_r+0.5*(F_0_rl[0]+F_0_cr[0])*avg_b_r)*wvpar+dvpar*((0.14433756729740646*F_0_rl[1]-0.14433756729740646*F_0_cr[1])*max_b_r+0.14433756729740646*(F_0_rl[1]+F_0_cr[1])*avg_b_r); 
  Ghat_F_0_vpar_r[1] = ((0.5*F_0_rl[1]-0.5*F_0_cr[1])*max_b_r+0.5*(F_0_rl[1]+F_0_cr[1])*avg_b_r)*wvpar+dvpar*((0.12909944487358058*F_0_rl[2]-0.12909944487358058*F_0_cr[2]+0.14433756729740646*F_0_rl[0]-0.14433756729740646*F_0_cr[0])*max_b_r+(0.12909944487358058*(F_0_rl[2]+F_0_cr[2])+0.14433756729740646*(F_0_rl[0]+F_0_cr[0]))*avg_b_r); 
  Ghat_F_0_vpar_r[2] = ((0.5*F_0_rl[2]-0.5*F_0_cr[2])*max_b_r+0.5*(F_0_rl[2]+F_0_cr[2])*avg_b_r)*wvpar+dvpar*((0.12909944487358058*F_0_rl[1]-0.12909944487358058*F_0_cr[1])*max_b_r+0.12909944487358058*(F_0_rl[1]+F_0_cr[1])*avg_b_r); 

  Ghat_G_1_vpar_l[0] = (0.5*(G_1_cl[0]*max_b_l+(G_1_lr[0]+G_1_cl[0])*avg_b_l)-0.5*G_1_lr[0]*max_b_l)*wvpar+dvpar*(0.14433756729740646*(G_1_cl[1]*max_b_l+(G_1_lr[1]+G_1_cl[1])*avg_b_l)-0.14433756729740646*G_1_lr[1]*max_b_l); 
  Ghat_G_1_vpar_l[1] = (0.5*(G_1_cl[1]*max_b_l+(G_1_lr[1]+G_1_cl[1])*avg_b_l)-0.5*G_1_lr[1]*max_b_l)*wvpar+dvpar*((-(0.12909944487358058*G_1_lr[2])+0.12909944487358058*G_1_cl[2]-0.14433756729740646*G_1_lr[0]+0.14433756729740646*G_1_cl[0])*max_b_l+(0.12909944487358058*(G_1_lr[2]+G_1_cl[2])+0.14433756729740646*(G_1_lr[0]+G_1_cl[0]))*avg_b_l); 
  Ghat_G_1_vpar_l[2] = (0.5*(G_1_cl[2]*max_b_l+(G_1_lr[2]+G_1_cl[2])*avg_b_l)-0.5*G_1_lr[2]*max_b_l)*wvpar+dvpar*(0.12909944487358058*(G_1_cl[1]*max_b_l+(G_1_lr[1]+G_1_cl[1])*avg_b_l)-0.12909944487358058*G_1_lr[1]*max_b_l); 

  Ghat_G_1_vpar_r[0] = ((0.5*G_1_rl[0]-0.5*G_1_cr[0])*max_b_r+0.5*(G_1_rl[0]+G_1_cr[0])*avg_b_r)*wvpar+dvpar*((0.14433756729740646*G_1_rl[1]-0.14433756729740646*G_1_cr[1])*max_b_r+0.14433756729740646*(G_1_rl[1]+G_1_cr[1])*avg_b_r); 
  Ghat_G_1_vpar_r[1] = ((0.5*G_1_rl[1]-0.5*G_1_cr[1])*max_b_r+0.5*(G_1_rl[1]+G_1_cr[1])*avg_b_r)*wvpar+dvpar*((0.12909944487358058*G_1_rl[2]-0.12909944487358058*G_1_cr[2]+0.14433756729740646*G_1_rl[0]-0.14433756729740646*G_1_cr[0])*max_b_r+(0.12909944487358058*(G_1_rl[2]+G_1_cr[2])+0.14433756729740646*(G_1_rl[0]+G_1_cr[0]))*avg_b_r); 
  Ghat_G_1_vpar_r[2] = ((0.5*G_1_rl[2]-0.5*G_1_cr[2])*max_b_r+0.5*(G_1_rl[2]+G_1_cr[2])*avg_b_r)*wvpar+dvpar*((0.12909944487358058*G_1_rl[1]-0.12909944487358058*G_1_cr[1])*max_b_r+0.12909944487358058*(G_1_rl[1]+G_1_cr[1])*avg_b_r); 

  } 

  flux_div_q[0] += dx1*mass*volFact*((0.5*Ghat_F_0_vpar_r[0]-0.5*Ghat_F_0_vpar_l[0])*wvpar_sq+(0.2886751345948129*Ghat_F_0_vpar_r[1]-0.2886751345948129*Ghat_F_0_vpar_l[1])*dvpar*wvpar+(0.03726779962499649*Ghat_F_0_vpar_r[2]-0.03726779962499649*Ghat_F_0_vpar_l[2]+0.041666666666666664*Ghat_F_0_vpar_r[0]-0.041666666666666664*Ghat_F_0_vpar_l[0])*dvpar_sq); 
  flux_div_q[1] += dx1*mass*volFact*((0.8660254037844386*(Ghat_F_0_vpar_r[0]+Ghat_F_0_vpar_l[0])-0.6123724356957944*(F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar_sq+(-(0.31622776601683794*alpha_c[3]*F_0c[5])-0.3162277660168379*alpha_c[2]*F_0c[4]-0.3535533905932737*(F_0c[1]*alpha_c[3]+alpha_c[1]*F_0c[3]+F_0c[0]*alpha_c[2]+alpha_c[0]*F_0c[2])+0.5*(Ghat_F_0_vpar_r[1]+Ghat_F_0_vpar_l[1]))*dvpar*wvpar+(-(0.04564354645876382*alpha_c[1]*F_0c[5])-0.04564354645876383*alpha_c[0]*F_0c[4]-0.09185586535436913*(F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2])+0.06454972243679029*(Ghat_F_0_vpar_r[2]+Ghat_F_0_vpar_l[2])-0.05103103630798286*(F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0])+0.07216878364870323*(Ghat_F_0_vpar_r[0]+Ghat_F_0_vpar_l[0]))*dvpar_sq); 

  flux_div_q[0] += (Ghat_G_1_vpar_r[0]-1.0*Ghat_G_1_vpar_l[0])*dx1*mass*volFact; 
  flux_div_q[1] += (1.7320508075688772*(Ghat_G_1_vpar_r[0]+Ghat_G_1_vpar_l[0])-1.224744871391589*(G_1c[3]*alpha_c[3]+G_1c[2]*alpha_c[2]+G_1c[1]*alpha_c[1]+G_1c[0]*alpha_c[0]))*dx1*mass*volFact; 

} 
