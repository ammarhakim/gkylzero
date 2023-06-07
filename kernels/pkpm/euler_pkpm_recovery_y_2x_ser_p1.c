#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p1(const double *dxv, double nuHyp, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *div_p_cell_avgl, const double *div_p_cell_avgc, const double *div_p_cell_avgr, 
  const double *statevecl, const double *statevecc, const double *statevecr, 
  const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
  const double *T_perp_over_m_inv, const double *nu, 
  double* div_p, double* pkpm_accel_vars) 
{ 
  // dxv[NDIM]:             Cell spacing.
  // nuHyp:                 Hyper-diffusion coefficient.
  // bvarl/c/r:             Input magnetic field unit vector in left/center/right cells.
  // u_il/c/r:              Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // p_ijl/c/r:             Input pressure tensor in left/center/right cells.
  // div_p_cell_avgl/c/r:   Input cell average d/dx_i p_ij for limiting div(p).
  // statevecl/c/r:         [rho ux, rho uy, rho uz], Fluid input state vector in center cell.
  // pkpm_div_ppar:         Input div(p_parallel b_hat) in center cell.
  // rho_inv:               Input 1/rho in center cell.
  // T_perp_over_m:         Input p_perp/rho = T_perp/m in center cell.
  // T_perp_over_m_inv:     Input (T_perp/m)^-1 in center cell.
  // nu:                    Input collisionality in center cell.
  // div_p:                 Increment to volume expansion of div(p) in one direction; includes hyper-diffusion for momentum.
  // pkpm_accel_vars:       Increment to volume expansion of pkpm acceleration variables.

  const double dx1 = 2.0/dxv[1]; 
  const double *b_l = &bvarl[4]; 
  const double *b_c = &bvarc[4]; 
  const double *b_r = &bvarr[4]; 

  const double *pkpm_div_ppar_c = &pkpm_div_ppar[4]; 

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

  const double *div_p_x_cell_avg_l = &div_p_cell_avgl[3]; 
  const double *div_p_x_cell_avg_c = &div_p_cell_avgc[3]; 
  const double *div_p_x_cell_avg_r = &div_p_cell_avgr[3]; 
  const double *div_p_y_cell_avg_l = &div_p_cell_avgl[4]; 
  const double *div_p_y_cell_avg_c = &div_p_cell_avgc[4]; 
  const double *div_p_y_cell_avg_r = &div_p_cell_avgr[4]; 
  const double *div_p_z_cell_avg_l = &div_p_cell_avgl[5]; 
  const double *div_p_z_cell_avg_c = &div_p_cell_avgc[5]; 
  const double *div_p_z_cell_avg_r = &div_p_cell_avgr[5]; 

  const double *grad_u_x_cell_avg_l = &div_p_cell_avgl[12]; 
  const double *grad_u_x_cell_avg_c = &div_p_cell_avgc[12]; 
  const double *grad_u_x_cell_avg_r = &div_p_cell_avgr[12]; 
  const double *grad_u_y_cell_avg_l = &div_p_cell_avgl[13]; 
  const double *grad_u_y_cell_avg_c = &div_p_cell_avgc[13]; 
  const double *grad_u_y_cell_avg_r = &div_p_cell_avgr[13]; 
  const double *grad_u_z_cell_avg_l = &div_p_cell_avgl[14]; 
  const double *grad_u_z_cell_avg_c = &div_p_cell_avgc[14]; 
  const double *grad_u_z_cell_avg_r = &div_p_cell_avgr[14]; 

  double r_p_x_l = div_p_x_cell_avg_l[0]/div_p_x_cell_avg_c[0]; 
  double r_p_x_r = div_p_x_cell_avg_c[0]/div_p_x_cell_avg_r[0]; 
  double r_p_y_l = div_p_y_cell_avg_l[0]/div_p_y_cell_avg_c[0]; 
  double r_p_y_r = div_p_y_cell_avg_c[0]/div_p_y_cell_avg_r[0]; 
  double r_p_z_l = div_p_z_cell_avg_l[0]/div_p_z_cell_avg_c[0]; 
  double r_p_z_r = div_p_z_cell_avg_c[0]/div_p_z_cell_avg_r[0]; 

  double r_u_x_l = grad_u_x_cell_avg_l[0]/grad_u_x_cell_avg_c[0]; 
  double r_u_x_r = grad_u_x_cell_avg_c[0]/grad_u_x_cell_avg_r[0]; 
  double r_u_y_l = grad_u_y_cell_avg_l[0]/grad_u_y_cell_avg_c[0]; 
  double r_u_y_r = grad_u_y_cell_avg_c[0]/grad_u_y_cell_avg_r[0]; 
  double r_u_z_l = grad_u_z_cell_avg_l[0]/grad_u_z_cell_avg_c[0]; 
  double r_u_z_r = grad_u_z_cell_avg_c[0]/grad_u_z_cell_avg_r[0]; 

  double minmod_x_l = fmax(0.0, fmin(fmin(1.0, r_p_x_l), fmin(1.0, r_u_x_l))); 
  double minmod_x_r = fmax(0.0, fmin(fmin(1.0, r_p_x_r), fmin(1.0, r_u_x_r))); 
  double minmod_y_l = fmax(0.0, fmin(fmin(1.0, r_p_y_l), fmin(1.0, r_u_y_l))); 
  double minmod_y_r = fmax(0.0, fmin(fmin(1.0, r_p_y_r), fmin(1.0, r_u_y_r))); 
  double minmod_z_l = fmax(0.0, fmin(fmin(1.0, r_p_z_l), fmin(1.0, r_u_z_l))); 
  double minmod_z_r = fmax(0.0, fmin(fmin(1.0, r_p_z_r), fmin(1.0, r_u_z_r))); 

  double phi_x_l = minmod_x_l == 0.0 ? minmod_x_l : 1.0; 
  double phi_x_r = minmod_x_r == 0.0 ? minmod_x_r : 1.0; 
  double phi_y_l = minmod_y_l == 0.0 ? minmod_y_l : 1.0; 
  double phi_y_r = minmod_y_r == 0.0 ? minmod_y_r : 1.0; 
  double phi_z_l = minmod_z_l == 0.0 ? minmod_z_l : 1.0; 
  double phi_z_r = minmod_z_r == 0.0 ? minmod_z_r : 1.0; 

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

  double *div_b = &pkpm_accel_vars[0]; 
  double *bb_grad_u = &pkpm_accel_vars[4]; 
  double *p_force = &pkpm_accel_vars[8]; 
  double *p_perp_source = &pkpm_accel_vars[12]; 
  double *p_perp_div_b = &pkpm_accel_vars[16]; 

  const double dxHyp = dx1*dx1*dx1*dx1; 

  double grad_u_x[4] = {0.0}; 
  double grad_u_y[4] = {0.0}; 
  double grad_u_z[4] = {0.0}; 
  double div_p_x_comp[4] = {0.0}; 
  double div_p_y_comp[4] = {0.0}; 
  double div_p_z_comp[4] = {0.0}; 
  grad_u_x[2] += -1.732050807568877*ux_c[0]*dx1; 
  grad_u_x[3] += -1.732050807568877*ux_c[1]*dx1; 
  div_p_x_comp[2] += -1.732050807568877*Pxy_c[0]*dx1; 
  div_p_x_comp[3] += -1.732050807568877*Pxy_c[1]*dx1; 
  if (phi_x_l == 0.0) { 
  grad_u_x[0] += ((-0.4330127018922193*ux_l[2])+0.4330127018922193*ux_c[2]-0.25*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += ((-0.4330127018922193*ux_l[3])+0.4330127018922193*ux_c[3]-0.25*(ux_l[1]+ux_c[1]))*dx1; 
  grad_u_x[2] += (0.75*ux_l[2]-0.75*ux_c[2]+0.4330127018922193*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[3] += (0.75*ux_l[3]-0.75*ux_c[3]+0.4330127018922193*(ux_l[1]+ux_c[1]))*dx1; 
  div_p_x_comp[0] += ((-0.4330127018922193*Pxy_l[2])+0.4330127018922193*Pxy_c[2]-0.25*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[1] += ((-0.4330127018922193*Pxy_l[3])+0.4330127018922193*Pxy_c[3]-0.25*(Pxy_l[1]+Pxy_c[1]))*dx1; 
  div_p_x_comp[2] += (0.75*Pxy_l[2]-0.75*Pxy_c[2]+0.4330127018922193*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[3] += (0.75*Pxy_l[3]-0.75*Pxy_c[3]+0.4330127018922193*(Pxy_l[1]+Pxy_c[1]))*dx1; 
  } else { 
  grad_u_x[0] += ((-0.2886751345948129*ux_l[2])+0.2886751345948129*ux_c[2]-0.25*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += ((-0.2886751345948129*ux_l[3])+0.2886751345948129*ux_c[3]-0.25*(ux_l[1]+ux_c[1]))*dx1; 
  grad_u_x[2] += (0.5*ux_l[2]-0.5*ux_c[2]+0.4330127018922193*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[3] += (0.5*ux_l[3]-0.5*ux_c[3]+0.4330127018922193*(ux_l[1]+ux_c[1]))*dx1; 
  div_p_x_comp[0] += ((-0.2886751345948129*Pxy_l[2])+0.2886751345948129*Pxy_c[2]-0.25*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[1] += ((-0.2886751345948129*Pxy_l[3])+0.2886751345948129*Pxy_c[3]-0.25*(Pxy_l[1]+Pxy_c[1]))*dx1; 
  div_p_x_comp[2] += (0.5*Pxy_l[2]-0.5*Pxy_c[2]+0.4330127018922193*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[3] += (0.5*Pxy_l[3]-0.5*Pxy_c[3]+0.4330127018922193*(Pxy_l[1]+Pxy_c[1]))*dx1; 
  } 
  if (phi_x_r == 0.0) { 
  grad_u_x[0] += ((-0.4330127018922193*ux_r[2])+0.4330127018922193*ux_c[2]+0.25*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += ((-0.4330127018922193*ux_r[3])+0.4330127018922193*ux_c[3]+0.25*(ux_r[1]+ux_c[1]))*dx1; 
  grad_u_x[2] += ((-0.75*ux_r[2])+0.75*ux_c[2]+0.4330127018922193*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[3] += ((-0.75*ux_r[3])+0.75*ux_c[3]+0.4330127018922193*(ux_r[1]+ux_c[1]))*dx1; 
  div_p_x_comp[0] += ((-0.4330127018922193*Pxy_r[2])+0.4330127018922193*Pxy_c[2]+0.25*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[1] += ((-0.4330127018922193*Pxy_r[3])+0.4330127018922193*Pxy_c[3]+0.25*(Pxy_r[1]+Pxy_c[1]))*dx1; 
  div_p_x_comp[2] += ((-0.75*Pxy_r[2])+0.75*Pxy_c[2]+0.4330127018922193*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[3] += ((-0.75*Pxy_r[3])+0.75*Pxy_c[3]+0.4330127018922193*(Pxy_r[1]+Pxy_c[1]))*dx1; 
  } else { 
  grad_u_x[0] += ((-0.2886751345948129*ux_r[2])+0.2886751345948129*ux_c[2]+0.25*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += ((-0.2886751345948129*ux_r[3])+0.2886751345948129*ux_c[3]+0.25*(ux_r[1]+ux_c[1]))*dx1; 
  grad_u_x[2] += ((-0.5*ux_r[2])+0.5*ux_c[2]+0.4330127018922193*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[3] += ((-0.5*ux_r[3])+0.5*ux_c[3]+0.4330127018922193*(ux_r[1]+ux_c[1]))*dx1; 
  div_p_x_comp[0] += ((-0.2886751345948129*Pxy_r[2])+0.2886751345948129*Pxy_c[2]+0.25*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[1] += ((-0.2886751345948129*Pxy_r[3])+0.2886751345948129*Pxy_c[3]+0.25*(Pxy_r[1]+Pxy_c[1]))*dx1; 
  div_p_x_comp[2] += ((-0.5*Pxy_r[2])+0.5*Pxy_c[2]+0.4330127018922193*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_x_comp[3] += ((-0.5*Pxy_r[3])+0.5*Pxy_c[3]+0.4330127018922193*(Pxy_r[1]+Pxy_c[1]))*dx1; 
  } 

  grad_u_y[2] += -1.732050807568877*uy_c[0]*dx1; 
  grad_u_y[3] += -1.732050807568877*uy_c[1]*dx1; 
  div_p_y_comp[2] += -1.732050807568877*Pyy_c[0]*dx1; 
  div_p_y_comp[3] += -1.732050807568877*Pyy_c[1]*dx1; 
  if (phi_y_l == 0.0) { 
  grad_u_y[0] += ((-0.4330127018922193*uy_l[2])+0.4330127018922193*uy_c[2]-0.25*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += ((-0.4330127018922193*uy_l[3])+0.4330127018922193*uy_c[3]-0.25*(uy_l[1]+uy_c[1]))*dx1; 
  grad_u_y[2] += (0.75*uy_l[2]-0.75*uy_c[2]+0.4330127018922193*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[3] += (0.75*uy_l[3]-0.75*uy_c[3]+0.4330127018922193*(uy_l[1]+uy_c[1]))*dx1; 
  div_p_y_comp[0] += ((-0.4330127018922193*Pyy_l[2])+0.4330127018922193*Pyy_c[2]-0.25*(Pyy_l[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[1] += ((-0.4330127018922193*Pyy_l[3])+0.4330127018922193*Pyy_c[3]-0.25*(Pyy_l[1]+Pyy_c[1]))*dx1; 
  div_p_y_comp[2] += (0.75*Pyy_l[2]-0.75*Pyy_c[2]+0.4330127018922193*(Pyy_l[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[3] += (0.75*Pyy_l[3]-0.75*Pyy_c[3]+0.4330127018922193*(Pyy_l[1]+Pyy_c[1]))*dx1; 
  } else { 
  grad_u_y[0] += ((-0.2886751345948129*uy_l[2])+0.2886751345948129*uy_c[2]-0.25*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += ((-0.2886751345948129*uy_l[3])+0.2886751345948129*uy_c[3]-0.25*(uy_l[1]+uy_c[1]))*dx1; 
  grad_u_y[2] += (0.5*uy_l[2]-0.5*uy_c[2]+0.4330127018922193*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[3] += (0.5*uy_l[3]-0.5*uy_c[3]+0.4330127018922193*(uy_l[1]+uy_c[1]))*dx1; 
  div_p_y_comp[0] += ((-0.2886751345948129*Pyy_l[2])+0.2886751345948129*Pyy_c[2]-0.25*(Pyy_l[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[1] += ((-0.2886751345948129*Pyy_l[3])+0.2886751345948129*Pyy_c[3]-0.25*(Pyy_l[1]+Pyy_c[1]))*dx1; 
  div_p_y_comp[2] += (0.5*Pyy_l[2]-0.5*Pyy_c[2]+0.4330127018922193*(Pyy_l[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[3] += (0.5*Pyy_l[3]-0.5*Pyy_c[3]+0.4330127018922193*(Pyy_l[1]+Pyy_c[1]))*dx1; 
  } 
  if (phi_y_r == 0.0) { 
  grad_u_y[0] += ((-0.4330127018922193*uy_r[2])+0.4330127018922193*uy_c[2]+0.25*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += ((-0.4330127018922193*uy_r[3])+0.4330127018922193*uy_c[3]+0.25*(uy_r[1]+uy_c[1]))*dx1; 
  grad_u_y[2] += ((-0.75*uy_r[2])+0.75*uy_c[2]+0.4330127018922193*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[3] += ((-0.75*uy_r[3])+0.75*uy_c[3]+0.4330127018922193*(uy_r[1]+uy_c[1]))*dx1; 
  div_p_y_comp[0] += ((-0.4330127018922193*Pyy_r[2])+0.4330127018922193*Pyy_c[2]+0.25*(Pyy_r[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[1] += ((-0.4330127018922193*Pyy_r[3])+0.4330127018922193*Pyy_c[3]+0.25*(Pyy_r[1]+Pyy_c[1]))*dx1; 
  div_p_y_comp[2] += ((-0.75*Pyy_r[2])+0.75*Pyy_c[2]+0.4330127018922193*(Pyy_r[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[3] += ((-0.75*Pyy_r[3])+0.75*Pyy_c[3]+0.4330127018922193*(Pyy_r[1]+Pyy_c[1]))*dx1; 
  } else { 
  grad_u_y[0] += ((-0.2886751345948129*uy_r[2])+0.2886751345948129*uy_c[2]+0.25*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += ((-0.2886751345948129*uy_r[3])+0.2886751345948129*uy_c[3]+0.25*(uy_r[1]+uy_c[1]))*dx1; 
  grad_u_y[2] += ((-0.5*uy_r[2])+0.5*uy_c[2]+0.4330127018922193*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[3] += ((-0.5*uy_r[3])+0.5*uy_c[3]+0.4330127018922193*(uy_r[1]+uy_c[1]))*dx1; 
  div_p_y_comp[0] += ((-0.2886751345948129*Pyy_r[2])+0.2886751345948129*Pyy_c[2]+0.25*(Pyy_r[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[1] += ((-0.2886751345948129*Pyy_r[3])+0.2886751345948129*Pyy_c[3]+0.25*(Pyy_r[1]+Pyy_c[1]))*dx1; 
  div_p_y_comp[2] += ((-0.5*Pyy_r[2])+0.5*Pyy_c[2]+0.4330127018922193*(Pyy_r[0]+Pyy_c[0]))*dx1; 
  div_p_y_comp[3] += ((-0.5*Pyy_r[3])+0.5*Pyy_c[3]+0.4330127018922193*(Pyy_r[1]+Pyy_c[1]))*dx1; 
  } 

  grad_u_z[2] += -1.732050807568877*uz_c[0]*dx1; 
  grad_u_z[3] += -1.732050807568877*uz_c[1]*dx1; 
  div_p_z_comp[2] += -1.732050807568877*Pyz_c[0]*dx1; 
  div_p_z_comp[3] += -1.732050807568877*Pyz_c[1]*dx1; 
  if (phi_z_l == 0.0) { 
  grad_u_z[0] += ((-0.4330127018922193*uz_l[2])+0.4330127018922193*uz_c[2]-0.25*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += ((-0.4330127018922193*uz_l[3])+0.4330127018922193*uz_c[3]-0.25*(uz_l[1]+uz_c[1]))*dx1; 
  grad_u_z[2] += (0.75*uz_l[2]-0.75*uz_c[2]+0.4330127018922193*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[3] += (0.75*uz_l[3]-0.75*uz_c[3]+0.4330127018922193*(uz_l[1]+uz_c[1]))*dx1; 
  div_p_z_comp[0] += ((-0.4330127018922193*Pyz_l[2])+0.4330127018922193*Pyz_c[2]-0.25*(Pyz_l[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[1] += ((-0.4330127018922193*Pyz_l[3])+0.4330127018922193*Pyz_c[3]-0.25*(Pyz_l[1]+Pyz_c[1]))*dx1; 
  div_p_z_comp[2] += (0.75*Pyz_l[2]-0.75*Pyz_c[2]+0.4330127018922193*(Pyz_l[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[3] += (0.75*Pyz_l[3]-0.75*Pyz_c[3]+0.4330127018922193*(Pyz_l[1]+Pyz_c[1]))*dx1; 
  } else { 
  grad_u_z[0] += ((-0.2886751345948129*uz_l[2])+0.2886751345948129*uz_c[2]-0.25*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += ((-0.2886751345948129*uz_l[3])+0.2886751345948129*uz_c[3]-0.25*(uz_l[1]+uz_c[1]))*dx1; 
  grad_u_z[2] += (0.5*uz_l[2]-0.5*uz_c[2]+0.4330127018922193*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[3] += (0.5*uz_l[3]-0.5*uz_c[3]+0.4330127018922193*(uz_l[1]+uz_c[1]))*dx1; 
  div_p_z_comp[0] += ((-0.2886751345948129*Pyz_l[2])+0.2886751345948129*Pyz_c[2]-0.25*(Pyz_l[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[1] += ((-0.2886751345948129*Pyz_l[3])+0.2886751345948129*Pyz_c[3]-0.25*(Pyz_l[1]+Pyz_c[1]))*dx1; 
  div_p_z_comp[2] += (0.5*Pyz_l[2]-0.5*Pyz_c[2]+0.4330127018922193*(Pyz_l[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[3] += (0.5*Pyz_l[3]-0.5*Pyz_c[3]+0.4330127018922193*(Pyz_l[1]+Pyz_c[1]))*dx1; 
  } 
  if (phi_z_r == 0.0) { 
  grad_u_z[0] += ((-0.4330127018922193*uz_r[2])+0.4330127018922193*uz_c[2]+0.25*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += ((-0.4330127018922193*uz_r[3])+0.4330127018922193*uz_c[3]+0.25*(uz_r[1]+uz_c[1]))*dx1; 
  grad_u_z[2] += ((-0.75*uz_r[2])+0.75*uz_c[2]+0.4330127018922193*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[3] += ((-0.75*uz_r[3])+0.75*uz_c[3]+0.4330127018922193*(uz_r[1]+uz_c[1]))*dx1; 
  div_p_z_comp[0] += ((-0.4330127018922193*Pyz_r[2])+0.4330127018922193*Pyz_c[2]+0.25*(Pyz_r[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[1] += ((-0.4330127018922193*Pyz_r[3])+0.4330127018922193*Pyz_c[3]+0.25*(Pyz_r[1]+Pyz_c[1]))*dx1; 
  div_p_z_comp[2] += ((-0.75*Pyz_r[2])+0.75*Pyz_c[2]+0.4330127018922193*(Pyz_r[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[3] += ((-0.75*Pyz_r[3])+0.75*Pyz_c[3]+0.4330127018922193*(Pyz_r[1]+Pyz_c[1]))*dx1; 
  } else { 
  grad_u_z[0] += ((-0.2886751345948129*uz_r[2])+0.2886751345948129*uz_c[2]+0.25*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += ((-0.2886751345948129*uz_r[3])+0.2886751345948129*uz_c[3]+0.25*(uz_r[1]+uz_c[1]))*dx1; 
  grad_u_z[2] += ((-0.5*uz_r[2])+0.5*uz_c[2]+0.4330127018922193*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[3] += ((-0.5*uz_r[3])+0.5*uz_c[3]+0.4330127018922193*(uz_r[1]+uz_c[1]))*dx1; 
  div_p_z_comp[0] += ((-0.2886751345948129*Pyz_r[2])+0.2886751345948129*Pyz_c[2]+0.25*(Pyz_r[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[1] += ((-0.2886751345948129*Pyz_r[3])+0.2886751345948129*Pyz_c[3]+0.25*(Pyz_r[1]+Pyz_c[1]))*dx1; 
  div_p_z_comp[2] += ((-0.5*Pyz_r[2])+0.5*Pyz_c[2]+0.4330127018922193*(Pyz_r[0]+Pyz_c[0]))*dx1; 
  div_p_z_comp[3] += ((-0.5*Pyz_r[3])+0.5*Pyz_c[3]+0.4330127018922193*(Pyz_r[1]+Pyz_c[1]))*dx1; 
  } 

  div_p_x[0] += (4.871392896287466*rhoux_r[2]-4.871392896287466*rhoux_l[2]-2.8125*(rhoux_r[0]+rhoux_l[0])+5.625*rhoux_c[0])*dxHyp*nuHyp+div_p_x_comp[0]; 
  div_p_x[1] += (4.871392896287466*rhoux_r[3]-4.871392896287466*rhoux_l[3]-2.8125*(rhoux_r[1]+rhoux_l[1])+5.625*rhoux_c[1])*dxHyp*nuHyp+div_p_x_comp[1]; 
  div_p_x[2] += (72.1875*(rhoux_r[2]+rhoux_l[2])+249.375*rhoux_c[2]-56.83291712335378*rhoux_r[0]+56.83291712335378*rhoux_l[0])*dxHyp*nuHyp+div_p_x_comp[2]; 
  div_p_x[3] += (72.1875*(rhoux_r[3]+rhoux_l[3])+249.375*rhoux_c[3]-56.83291712335378*rhoux_r[1]+56.83291712335378*rhoux_l[1])*dxHyp*nuHyp+div_p_x_comp[3]; 

  div_p_y[0] += (4.871392896287466*rhouy_r[2]-4.871392896287466*rhouy_l[2]-2.8125*(rhouy_r[0]+rhouy_l[0])+5.625*rhouy_c[0])*dxHyp*nuHyp+div_p_y_comp[0]; 
  div_p_y[1] += (4.871392896287466*rhouy_r[3]-4.871392896287466*rhouy_l[3]-2.8125*(rhouy_r[1]+rhouy_l[1])+5.625*rhouy_c[1])*dxHyp*nuHyp+div_p_y_comp[1]; 
  div_p_y[2] += (72.1875*(rhouy_r[2]+rhouy_l[2])+249.375*rhouy_c[2]-56.83291712335378*rhouy_r[0]+56.83291712335378*rhouy_l[0])*dxHyp*nuHyp+div_p_y_comp[2]; 
  div_p_y[3] += (72.1875*(rhouy_r[3]+rhouy_l[3])+249.375*rhouy_c[3]-56.83291712335378*rhouy_r[1]+56.83291712335378*rhouy_l[1])*dxHyp*nuHyp+div_p_y_comp[3]; 

  div_p_z[0] += (4.871392896287466*rhouz_r[2]-4.871392896287466*rhouz_l[2]-2.8125*(rhouz_r[0]+rhouz_l[0])+5.625*rhouz_c[0])*dxHyp*nuHyp+div_p_z_comp[0]; 
  div_p_z[1] += (4.871392896287466*rhouz_r[3]-4.871392896287466*rhouz_l[3]-2.8125*(rhouz_r[1]+rhouz_l[1])+5.625*rhouz_c[1])*dxHyp*nuHyp+div_p_z_comp[1]; 
  div_p_z[2] += (72.1875*(rhouz_r[2]+rhouz_l[2])+249.375*rhouz_c[2]-56.83291712335378*rhouz_r[0]+56.83291712335378*rhouz_l[0])*dxHyp*nuHyp+div_p_z_comp[2]; 
  div_p_z[3] += (72.1875*(rhouz_r[3]+rhouz_l[3])+249.375*rhouz_c[3]-56.83291712335378*rhouz_r[1]+56.83291712335378*rhouz_l[1])*dxHyp*nuHyp+div_p_z_comp[3]; 

  double div_b_comp[4] = {0.0}; 
  double bb_grad_u_comp[4] = {0.0}; 
  div_b_comp[0] = (-0.2886751345948129*b_r[2]*dx1)-0.2886751345948129*b_l[2]*dx1+0.5773502691896258*b_c[2]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.2886751345948129*b_r[3]*dx1)-0.2886751345948129*b_l[3]*dx1+0.5773502691896258*b_c[3]*dx1+0.25*b_r[1]*dx1-0.25*b_l[1]*dx1; 
  div_b_comp[2] = (-0.5*b_r[2]*dx1)+0.5*b_l[2]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 
  div_b_comp[3] = (-0.5*b_r[3]*dx1)+0.5*b_l[3]*dx1+0.4330127018922193*b_r[1]*dx1+0.4330127018922193*b_l[1]*dx1-0.8660254037844386*b_c[1]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 
  div_b[3] += div_b_comp[3]; 

  bb_grad_u_comp[0] = 0.5*bybz[3]*grad_u_z[3]+0.5*byby[3]*grad_u_y[3]+0.5*bxby[3]*grad_u_x[3]+0.5*bybz[2]*grad_u_z[2]+0.5*byby[2]*grad_u_y[2]+0.5*bxby[2]*grad_u_x[2]+0.5*bybz[1]*grad_u_z[1]+0.5*byby[1]*grad_u_y[1]+0.5*bxby[1]*grad_u_x[1]+0.5*bybz[0]*grad_u_z[0]+0.5*byby[0]*grad_u_y[0]+0.5*bxby[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.5*bybz[2]*grad_u_z[3]+0.5*byby[2]*grad_u_y[3]+0.5*bxby[2]*grad_u_x[3]+0.5*grad_u_z[2]*bybz[3]+0.5*grad_u_y[2]*byby[3]+0.5*grad_u_x[2]*bxby[3]+0.5*bybz[0]*grad_u_z[1]+0.5*byby[0]*grad_u_y[1]+0.5*bxby[0]*grad_u_x[1]+0.5*grad_u_z[0]*bybz[1]+0.5*grad_u_y[0]*byby[1]+0.5*grad_u_x[0]*bxby[1]; 
  bb_grad_u_comp[2] = 0.5*bybz[1]*grad_u_z[3]+0.5*byby[1]*grad_u_y[3]+0.5*bxby[1]*grad_u_x[3]+0.5*grad_u_z[1]*bybz[3]+0.5*grad_u_y[1]*byby[3]+0.5*grad_u_x[1]*bxby[3]+0.5*bybz[0]*grad_u_z[2]+0.5*byby[0]*grad_u_y[2]+0.5*bxby[0]*grad_u_x[2]+0.5*grad_u_z[0]*bybz[2]+0.5*grad_u_y[0]*byby[2]+0.5*grad_u_x[0]*bxby[2]; 
  bb_grad_u_comp[3] = 0.5*bybz[0]*grad_u_z[3]+0.5*byby[0]*grad_u_y[3]+0.5*bxby[0]*grad_u_x[3]+0.5*grad_u_z[0]*bybz[3]+0.5*grad_u_y[0]*byby[3]+0.5*grad_u_x[0]*bxby[3]+0.5*bybz[1]*grad_u_z[2]+0.5*byby[1]*grad_u_y[2]+0.5*bxby[1]*grad_u_x[2]+0.5*grad_u_z[1]*bybz[2]+0.5*grad_u_y[1]*byby[2]+0.5*grad_u_x[1]*bxby[2]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 
  bb_grad_u[3] += bb_grad_u_comp[3]; 

  p_force[0] += 0.5*pkpm_div_ppar_c[3]*rho_inv[3]-0.5*T_perp_over_m[3]*div_b_comp[3]+0.5*pkpm_div_ppar_c[2]*rho_inv[2]-0.5*T_perp_over_m[2]*div_b_comp[2]+0.5*pkpm_div_ppar_c[1]*rho_inv[1]-0.5*T_perp_over_m[1]*div_b_comp[1]+0.5*pkpm_div_ppar_c[0]*rho_inv[0]-0.5*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += 0.5*(pkpm_div_ppar_c[2]*rho_inv[3]+rho_inv[2]*pkpm_div_ppar_c[3])-0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3])+0.5*(pkpm_div_ppar_c[0]*rho_inv[1]+rho_inv[0]*pkpm_div_ppar_c[1])-0.5*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force[2] += 0.5*(pkpm_div_ppar_c[1]*rho_inv[3]+rho_inv[1]*pkpm_div_ppar_c[3])-0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3])+0.5*(pkpm_div_ppar_c[0]*rho_inv[2]+rho_inv[0]*pkpm_div_ppar_c[2])-0.5*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_force[3] += 0.5*(pkpm_div_ppar_c[0]*rho_inv[3]+rho_inv[0]*pkpm_div_ppar_c[3])-0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3])+0.5*(pkpm_div_ppar_c[1]*rho_inv[2]+rho_inv[1]*pkpm_div_ppar_c[2])-0.5*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

  p_perp_source[0] += bb_grad_u_comp[0]-1.0*(nu[0]+grad_u_y[0]); 
  p_perp_source[1] += bb_grad_u_comp[1]-1.0*(nu[1]+grad_u_y[1]); 
  p_perp_source[2] += bb_grad_u_comp[2]-1.0*(nu[2]+grad_u_y[2]); 
  p_perp_source[3] += bb_grad_u_comp[3]-1.0*(nu[3]+grad_u_y[3]); 

  p_perp_div_b[0] += 0.5*(T_perp_over_m[3]*div_b_comp[3]+T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.5*(T_perp_over_m[2]*div_b_comp[3]+div_b_comp[2]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.5*(T_perp_over_m[1]*div_b_comp[3]+div_b_comp[1]*T_perp_over_m[3]+T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2]); 
  p_perp_div_b[3] += 0.5*(T_perp_over_m[0]*div_b_comp[3]+div_b_comp[0]*T_perp_over_m[3]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2]); 

} 
