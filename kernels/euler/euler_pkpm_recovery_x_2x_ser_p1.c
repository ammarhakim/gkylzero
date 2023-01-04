#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_ser_p1(const double *dxv, double nuHyp, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
  const double *statevecl, const double *statevecc, const double *statevecr, 
  const double *T_perp_over_m, const double *nu, const double *nu_vthsq, 
  double* div_b, double* bb_grad_u, double* div_p, double* p_force, double* p_perp_source, double* p_perp_div_b) 
{ 
  // dxv[NDIM]:             Cell spacing.
  // nuHyp:                 Hyper-diffusion coefficient.
  // bvarl/c/r:             Input magnetic field unit vector in left/center/right cells.
  // u_il/c/r:              Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // p_ijl/c/r:             Input pressure tensor in left/center/right cells.
  // vlasov_pkpm_momsl/c/r: Input pkpm moments (rho, p_parallel, p_perp) in left/center/right cells.
  // statevecl/c/r:         [rho ux, rho uy, rho uz], Fluid input state vector in center cell.
  // T_perp_over_m:         Input p_perp/rho = T_perp/m in center cell.
  // nu:                    Input collisionality in center cell.
  // nu_vthsq:              Input nu*vth^2 in center cell.
  // div_b:                 Increment to volume expansion of div(b) in one direction.
  // bb_grad_u:             Increment to volume expansion of bb : grad(u) in one direction.
  // div_p:                 Increment to volume expansion of div(p) in one direction.
  // p_force:               Increment to volume expansion of p_force = [1/rho * div(p_parallel b_hat),  1/rho * div(p_parallel b_hat) + 3 T_perp/m*div(b)] in one direction.
  // p_perp_source:         Increment to volume expansion of source for G_1 equation bb : grad(u) - div(u) - nu + nu m vth^2/T_perp .
  // p_perp_div_b:          Increment to volume expansion of p_perp/rho*div(b) = T_perp/m*div(b) in one direction.

  const double dx1 = 2.0/dxv[0]; 
  const double dx14 = dx1*dx1*dx1*dx1; 

  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[4]; 
  const double *uz_l = &u_il[8]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[4]; 
  const double *uz_c = &u_ic[8]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[4]; 
  const double *uz_r = &u_ir[8]; 

  const double *bxbx = &bvarc[12]; 
  const double *bxby = &bvarc[16]; 
  const double *bxbz = &bvarc[20]; 
  const double *byby = &bvarc[24]; 
  const double *bybz = &bvarc[28]; 
  const double *bzbz = &bvarc[32]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[4]; 
  const double *Pxz_l = &p_ijl[8]; 
  const double *Pyy_l = &p_ijl[12]; 
  const double *Pyz_l = &p_ijl[16]; 
  const double *Pzz_l = &p_ijl[20]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[4]; 
  const double *Pxz_c = &p_ijc[8]; 
  const double *Pyy_c = &p_ijc[12]; 
  const double *Pyz_c = &p_ijc[16]; 
  const double *Pzz_c = &p_ijc[20]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[4]; 
  const double *Pxz_r = &p_ijr[8]; 
  const double *Pyy_r = &p_ijr[12]; 
  const double *Pyz_r = &p_ijr[16]; 
  const double *Pzz_r = &p_ijr[20]; 

  const double *ppar_l = &vlasov_pkpm_momsl[4]; 
  const double *ppar_c = &vlasov_pkpm_momsc[4]; 
  const double *ppar_r = &vlasov_pkpm_momsr[4]; 
  const double *rho = &vlasov_pkpm_momsc[0]; 
  const double *p_perp = &vlasov_pkpm_momsc[8]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[4]; 
  const double *rhouz_l = &statevecl[8]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[4]; 
  const double *rhouz_c = &statevecc[8]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[4]; 
  const double *rhouz_r = &statevecr[8]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[4]; 
  double *div_p_z = &div_p[8]; 

  double *p_force_F_0 = &p_force[0]; 
  double *p_force_G_1 = &p_force[4]; 

  double grad_u_x[4] = {0.0}; 
  double grad_u_y[4] = {0.0}; 
  double grad_u_z[4] = {0.0}; 
  grad_u_x[0] = (-0.2886751345948129*ux_r[1]*dx1)-0.2886751345948129*ux_l[1]*dx1+0.5773502691896258*ux_c[1]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_u_x[1] = (-0.5*ux_r[1]*dx1)+0.5*ux_l[1]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844386*ux_c[0]*dx1; 
  grad_u_x[2] = (-0.2886751345948129*ux_r[3]*dx1)-0.2886751345948129*ux_l[3]*dx1+0.5773502691896258*ux_c[3]*dx1+0.25*ux_r[2]*dx1-0.25*ux_l[2]*dx1; 
  grad_u_x[3] = (-0.5*ux_r[3]*dx1)+0.5*ux_l[3]*dx1+0.4330127018922193*ux_r[2]*dx1+0.4330127018922193*ux_l[2]*dx1-0.8660254037844386*ux_c[2]*dx1; 

  grad_u_y[0] = (-0.2886751345948129*uy_r[1]*dx1)-0.2886751345948129*uy_l[1]*dx1+0.5773502691896258*uy_c[1]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_u_y[1] = (-0.5*uy_r[1]*dx1)+0.5*uy_l[1]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844386*uy_c[0]*dx1; 
  grad_u_y[2] = (-0.2886751345948129*uy_r[3]*dx1)-0.2886751345948129*uy_l[3]*dx1+0.5773502691896258*uy_c[3]*dx1+0.25*uy_r[2]*dx1-0.25*uy_l[2]*dx1; 
  grad_u_y[3] = (-0.5*uy_r[3]*dx1)+0.5*uy_l[3]*dx1+0.4330127018922193*uy_r[2]*dx1+0.4330127018922193*uy_l[2]*dx1-0.8660254037844386*uy_c[2]*dx1; 

  grad_u_z[0] = (-0.2886751345948129*uz_r[1]*dx1)-0.2886751345948129*uz_l[1]*dx1+0.5773502691896258*uz_c[1]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_u_z[1] = (-0.5*uz_r[1]*dx1)+0.5*uz_l[1]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844386*uz_c[0]*dx1; 
  grad_u_z[2] = (-0.2886751345948129*uz_r[3]*dx1)-0.2886751345948129*uz_l[3]*dx1+0.5773502691896258*uz_c[3]*dx1+0.25*uz_r[2]*dx1-0.25*uz_l[2]*dx1; 
  grad_u_z[3] = (-0.5*uz_r[3]*dx1)+0.5*uz_l[3]*dx1+0.4330127018922193*uz_r[2]*dx1+0.4330127018922193*uz_l[2]*dx1-0.8660254037844386*uz_c[2]*dx1; 

  double ppar_b_l[4] = {0.0}; 
  double ppar_b_c[4] = {0.0}; 
  double ppar_b_r[4] = {0.0}; 
  double div_b_comp[4] = {0.0}; 
  double bb_grad_u_comp[4] = {0.0}; 
  ppar_b_l[0] = 0.5*b_l[3]*ppar_l[3]+0.5*b_l[2]*ppar_l[2]+0.5*b_l[1]*ppar_l[1]+0.5*b_l[0]*ppar_l[0]; 
  ppar_b_l[1] = 0.5*b_l[2]*ppar_l[3]+0.5*ppar_l[2]*b_l[3]+0.5*b_l[0]*ppar_l[1]+0.5*ppar_l[0]*b_l[1]; 
  ppar_b_l[2] = 0.5*b_l[1]*ppar_l[3]+0.5*ppar_l[1]*b_l[3]+0.5*b_l[0]*ppar_l[2]+0.5*ppar_l[0]*b_l[2]; 
  ppar_b_l[3] = 0.5*b_l[0]*ppar_l[3]+0.5*ppar_l[0]*b_l[3]+0.5*b_l[1]*ppar_l[2]+0.5*ppar_l[1]*b_l[2]; 

  ppar_b_c[0] = 0.5*b_c[3]*ppar_c[3]+0.5*b_c[2]*ppar_c[2]+0.5*b_c[1]*ppar_c[1]+0.5*b_c[0]*ppar_c[0]; 
  ppar_b_c[1] = 0.5*b_c[2]*ppar_c[3]+0.5*ppar_c[2]*b_c[3]+0.5*b_c[0]*ppar_c[1]+0.5*ppar_c[0]*b_c[1]; 
  ppar_b_c[2] = 0.5*b_c[1]*ppar_c[3]+0.5*ppar_c[1]*b_c[3]+0.5*b_c[0]*ppar_c[2]+0.5*ppar_c[0]*b_c[2]; 
  ppar_b_c[3] = 0.5*b_c[0]*ppar_c[3]+0.5*ppar_c[0]*b_c[3]+0.5*b_c[1]*ppar_c[2]+0.5*ppar_c[1]*b_c[2]; 

  ppar_b_r[0] = 0.5*b_r[3]*ppar_r[3]+0.5*b_r[2]*ppar_r[2]+0.5*b_r[1]*ppar_r[1]+0.5*b_r[0]*ppar_r[0]; 
  ppar_b_r[1] = 0.5*b_r[2]*ppar_r[3]+0.5*ppar_r[2]*b_r[3]+0.5*b_r[0]*ppar_r[1]+0.5*ppar_r[0]*b_r[1]; 
  ppar_b_r[2] = 0.5*b_r[1]*ppar_r[3]+0.5*ppar_r[1]*b_r[3]+0.5*b_r[0]*ppar_r[2]+0.5*ppar_r[0]*b_r[2]; 
  ppar_b_r[3] = 0.5*b_r[0]*ppar_r[3]+0.5*ppar_r[0]*b_r[3]+0.5*b_r[1]*ppar_r[2]+0.5*ppar_r[1]*b_r[2]; 

  div_b_comp[0] = (-0.2886751345948129*b_r[1]*dx1)-0.2886751345948129*b_l[1]*dx1+0.5773502691896258*b_c[1]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.5*b_r[1]*dx1)+0.5*b_l[1]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 
  div_b_comp[2] = (-0.2886751345948129*b_r[3]*dx1)-0.2886751345948129*b_l[3]*dx1+0.5773502691896258*b_c[3]*dx1+0.25*b_r[2]*dx1-0.25*b_l[2]*dx1; 
  div_b_comp[3] = (-0.5*b_r[3]*dx1)+0.5*b_l[3]*dx1+0.4330127018922193*b_r[2]*dx1+0.4330127018922193*b_l[2]*dx1-0.8660254037844386*b_c[2]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 
  div_b[3] += div_b_comp[3]; 

  bb_grad_u_comp[0] = 0.5*bxbz[3]*grad_u_z[3]+0.5*bxby[3]*grad_u_y[3]+0.5*bxbx[3]*grad_u_x[3]+0.5*bxbz[2]*grad_u_z[2]+0.5*bxby[2]*grad_u_y[2]+0.5*bxbx[2]*grad_u_x[2]+0.5*bxbz[1]*grad_u_z[1]+0.5*bxby[1]*grad_u_y[1]+0.5*bxbx[1]*grad_u_x[1]+0.5*bxbz[0]*grad_u_z[0]+0.5*bxby[0]*grad_u_y[0]+0.5*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.5*bxbz[2]*grad_u_z[3]+0.5*bxby[2]*grad_u_y[3]+0.5*bxbx[2]*grad_u_x[3]+0.5*grad_u_z[2]*bxbz[3]+0.5*grad_u_y[2]*bxby[3]+0.5*grad_u_x[2]*bxbx[3]+0.5*bxbz[0]*grad_u_z[1]+0.5*bxby[0]*grad_u_y[1]+0.5*bxbx[0]*grad_u_x[1]+0.5*grad_u_z[0]*bxbz[1]+0.5*grad_u_y[0]*bxby[1]+0.5*grad_u_x[0]*bxbx[1]; 
  bb_grad_u_comp[2] = 0.5*bxbz[1]*grad_u_z[3]+0.5*bxby[1]*grad_u_y[3]+0.5*bxbx[1]*grad_u_x[3]+0.5*grad_u_z[1]*bxbz[3]+0.5*grad_u_y[1]*bxby[3]+0.5*grad_u_x[1]*bxbx[3]+0.5*bxbz[0]*grad_u_z[2]+0.5*bxby[0]*grad_u_y[2]+0.5*bxbx[0]*grad_u_x[2]+0.5*grad_u_z[0]*bxbz[2]+0.5*grad_u_y[0]*bxby[2]+0.5*grad_u_x[0]*bxbx[2]; 
  bb_grad_u_comp[3] = 0.5*bxbz[0]*grad_u_z[3]+0.5*bxby[0]*grad_u_y[3]+0.5*bxbx[0]*grad_u_x[3]+0.5*grad_u_z[0]*bxbz[3]+0.5*grad_u_y[0]*bxby[3]+0.5*grad_u_x[0]*bxbx[3]+0.5*bxbz[1]*grad_u_z[2]+0.5*bxby[1]*grad_u_y[2]+0.5*bxbx[1]*grad_u_x[2]+0.5*grad_u_z[1]*bxbz[2]+0.5*grad_u_y[1]*bxby[2]+0.5*grad_u_x[1]*bxbx[2]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 

  div_p_x[0] += (4.871392896287466*rhoux_r[1]-4.871392896287466*rhoux_l[1]-2.8125*(rhoux_r[0]+rhoux_l[0])+5.625*rhoux_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxx_r[1]+Pxx_l[1]))+0.5773502691896258*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_x[1] += (72.1875*(rhoux_r[1]+rhoux_l[1])+249.375*rhoux_c[1]-56.83291712335378*rhoux_r[0]+56.83291712335378*rhoux_l[0])*dx14*nuHyp+((-0.5*Pxx_r[1])+0.5*Pxx_l[1]+0.4330127018922193*(Pxx_r[0]+Pxx_l[0])-0.8660254037844386*Pxx_c[0])*dx1; 
  div_p_x[2] += (4.871392896287466*rhoux_r[3]-4.871392896287466*rhoux_l[3]-2.8125*(rhoux_r[2]+rhoux_l[2])+5.625*rhoux_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxx_r[3]+Pxx_l[3]))+0.5773502691896258*Pxx_c[3]+0.25*Pxx_r[2]-0.25*Pxx_l[2])*dx1; 
  div_p_x[3] += (72.1875*(rhoux_r[3]+rhoux_l[3])+249.375*rhoux_c[3]-56.83291712335378*rhoux_r[2]+56.83291712335378*rhoux_l[2])*dx14*nuHyp+((-0.5*Pxx_r[3])+0.5*Pxx_l[3]+0.4330127018922193*(Pxx_r[2]+Pxx_l[2])-0.8660254037844386*Pxx_c[2])*dx1; 

  div_p_y[0] += (4.871392896287466*rhouy_r[1]-4.871392896287466*rhouy_l[1]-2.8125*(rhouy_r[0]+rhouy_l[0])+5.625*rhouy_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxy_r[1]+Pxy_l[1]))+0.5773502691896258*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_y[1] += (72.1875*(rhouy_r[1]+rhouy_l[1])+249.375*rhouy_c[1]-56.83291712335378*rhouy_r[0]+56.83291712335378*rhouy_l[0])*dx14*nuHyp+((-0.5*Pxy_r[1])+0.5*Pxy_l[1]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_y[2] += (4.871392896287466*rhouy_r[3]-4.871392896287466*rhouy_l[3]-2.8125*(rhouy_r[2]+rhouy_l[2])+5.625*rhouy_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxy_r[3]+Pxy_l[3]))+0.5773502691896258*Pxy_c[3]+0.25*Pxy_r[2]-0.25*Pxy_l[2])*dx1; 
  div_p_y[3] += (72.1875*(rhouy_r[3]+rhouy_l[3])+249.375*rhouy_c[3]-56.83291712335378*rhouy_r[2]+56.83291712335378*rhouy_l[2])*dx14*nuHyp+((-0.5*Pxy_r[3])+0.5*Pxy_l[3]+0.4330127018922193*(Pxy_r[2]+Pxy_l[2])-0.8660254037844386*Pxy_c[2])*dx1; 

  div_p_z[0] += (4.871392896287466*rhouz_r[1]-4.871392896287466*rhouz_l[1]-2.8125*(rhouz_r[0]+rhouz_l[0])+5.625*rhouz_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[1]+Pxz_l[1]))+0.5773502691896258*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_z[1] += (72.1875*(rhouz_r[1]+rhouz_l[1])+249.375*rhouz_c[1]-56.83291712335378*rhouz_r[0]+56.83291712335378*rhouz_l[0])*dx14*nuHyp+((-0.5*Pxz_r[1])+0.5*Pxz_l[1]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 
  div_p_z[2] += (4.871392896287466*rhouz_r[3]-4.871392896287466*rhouz_l[3]-2.8125*(rhouz_r[2]+rhouz_l[2])+5.625*rhouz_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[3]+Pxz_l[3]))+0.5773502691896258*Pxz_c[3]+0.25*Pxz_r[2]-0.25*Pxz_l[2])*dx1; 
  div_p_z[3] += (72.1875*(rhouz_r[3]+rhouz_l[3])+249.375*rhouz_c[3]-56.83291712335378*rhouz_r[2]+56.83291712335378*rhouz_l[2])*dx14*nuHyp+((-0.5*Pxz_r[3])+0.5*Pxz_l[3]+0.4330127018922193*(Pxz_r[2]+Pxz_l[2])-0.8660254037844386*Pxz_c[2])*dx1; 

  double div_ppar_b[4] = {0.0}; 
  double rho_inv[4] = {0.0}; 
  double T_perp_over_m_inv[4] = {0.0}; 
  ser_2x_p1_inv(rho, rho_inv); 
  ser_2x_p1_inv(T_perp_over_m,T_perp_over_m_inv); 
  div_ppar_b[0] = (-0.2886751345948129*ppar_b_r[1]*dx1)-0.2886751345948129*ppar_b_l[1]*dx1+0.5773502691896258*ppar_b_c[1]*dx1+0.25*ppar_b_r[0]*dx1-0.25*ppar_b_l[0]*dx1; 
  div_ppar_b[1] = (-0.5*ppar_b_r[1]*dx1)+0.5*ppar_b_l[1]*dx1+0.4330127018922193*ppar_b_r[0]*dx1+0.4330127018922193*ppar_b_l[0]*dx1-0.8660254037844386*ppar_b_c[0]*dx1; 
  div_ppar_b[2] = (-0.2886751345948129*ppar_b_r[3]*dx1)-0.2886751345948129*ppar_b_l[3]*dx1+0.5773502691896258*ppar_b_c[3]*dx1+0.25*ppar_b_r[2]*dx1-0.25*ppar_b_l[2]*dx1; 
  div_ppar_b[3] = (-0.5*ppar_b_r[3]*dx1)+0.5*ppar_b_l[3]*dx1+0.4330127018922193*ppar_b_r[2]*dx1+0.4330127018922193*ppar_b_l[2]*dx1-0.8660254037844386*ppar_b_c[2]*dx1; 

  p_force_F_0[0] += 0.5*div_ppar_b[3]*rho_inv[3]-0.5*T_perp_over_m[3]*div_b_comp[3]+0.5*div_ppar_b[2]*rho_inv[2]-0.5*T_perp_over_m[2]*div_b_comp[2]+0.5*div_ppar_b[1]*rho_inv[1]-0.5*T_perp_over_m[1]*div_b_comp[1]+0.5*div_ppar_b[0]*rho_inv[0]-0.5*T_perp_over_m[0]*div_b_comp[0]; 
  p_force_F_0[1] += 0.5*(div_ppar_b[2]*rho_inv[3]+rho_inv[2]*div_ppar_b[3])-0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3])+0.5*(div_ppar_b[0]*rho_inv[1]+rho_inv[0]*div_ppar_b[1])-0.5*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force_F_0[2] += 0.5*(div_ppar_b[1]*rho_inv[3]+rho_inv[1]*div_ppar_b[3])-0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3])+0.5*(div_ppar_b[0]*rho_inv[2]+rho_inv[0]*div_ppar_b[2])-0.5*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force_F_0[3] += 0.5*(div_ppar_b[0]*rho_inv[3]+rho_inv[0]*div_ppar_b[3])-0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3])+0.5*(div_ppar_b[1]*rho_inv[2]+rho_inv[1]*div_ppar_b[2])-0.5*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 
  p_force_G_1[0] += 0.5*(div_ppar_b[3]*rho_inv[3]+T_perp_over_m[3]*div_b_comp[3]+div_ppar_b[2]*rho_inv[2]+T_perp_over_m[2]*div_b_comp[2]+div_ppar_b[1]*rho_inv[1]+T_perp_over_m[1]*div_b_comp[1]+div_ppar_b[0]*rho_inv[0]+T_perp_over_m[0]*div_b_comp[0]); 
  p_force_G_1[1] += 0.5*(div_ppar_b[2]*rho_inv[3]+rho_inv[2]*div_ppar_b[3]+T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]+div_ppar_b[0]*rho_inv[1]+rho_inv[0]*div_ppar_b[1]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force_G_1[2] += 0.5*(div_ppar_b[1]*rho_inv[3]+rho_inv[1]*div_ppar_b[3]+T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]+div_ppar_b[0]*rho_inv[2]+rho_inv[0]*div_ppar_b[2]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force_G_1[3] += 0.5*(div_ppar_b[0]*rho_inv[3]+rho_inv[0]*div_ppar_b[3]+T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+div_ppar_b[1]*rho_inv[2]+rho_inv[1]*div_ppar_b[2]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

  p_perp_source[0] += 0.25*(T_perp_over_m_inv[3]*nu_vthsq[3]+T_perp_over_m_inv[2]*nu_vthsq[2]+T_perp_over_m_inv[1]*nu_vthsq[1]+T_perp_over_m_inv[0]*nu_vthsq[0])-0.5*nu[0]-1.0*grad_u_x[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += 0.25*(T_perp_over_m_inv[2]*nu_vthsq[3]+nu_vthsq[2]*T_perp_over_m_inv[3]+T_perp_over_m_inv[0]*nu_vthsq[1])-0.5*nu[1]-1.0*grad_u_x[1]+bb_grad_u_comp[1]+0.25*nu_vthsq[0]*T_perp_over_m_inv[1]; 
  p_perp_source[2] += 0.25*(T_perp_over_m_inv[1]*nu_vthsq[3]+nu_vthsq[1]*T_perp_over_m_inv[3]+T_perp_over_m_inv[0]*nu_vthsq[2])-0.5*nu[2]-1.0*grad_u_x[2]+bb_grad_u_comp[2]+0.25*nu_vthsq[0]*T_perp_over_m_inv[2]; 
  p_perp_source[3] += 0.25*T_perp_over_m_inv[0]*nu_vthsq[3]-0.5*nu[3]-1.0*grad_u_x[3]+bb_grad_u_comp[3]+0.25*(nu_vthsq[0]*T_perp_over_m_inv[3]+T_perp_over_m_inv[1]*nu_vthsq[2]+nu_vthsq[1]*T_perp_over_m_inv[2]); 

  p_perp_div_b[0] += 0.5*(T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

} 
