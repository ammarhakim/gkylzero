#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p1(const double *dxv, double nuHyp, 
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

  const double dx1 = 2.0/dxv[0]; 
  const double *b_l = &bvarl[0]; 
  const double *b_c = &bvarc[0]; 
  const double *b_r = &bvarr[0]; 

  const double *pkpm_div_ppar_c = &pkpm_div_ppar[0]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[2]; 
  const double *uz_l = &u_il[4]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[2]; 
  const double *uz_c = &u_ic[4]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[2]; 
  const double *uz_r = &u_ir[4]; 

  const double *bxbx = &bvarc[6]; 
  const double *bxby = &bvarc[8]; 
  const double *bxbz = &bvarc[10]; 
  const double *byby = &bvarc[12]; 
  const double *bybz = &bvarc[14]; 
  const double *bzbz = &bvarc[16]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[2]; 
  const double *Pxz_l = &p_ijl[4]; 
  const double *Pyy_l = &p_ijl[6]; 
  const double *Pyz_l = &p_ijl[8]; 
  const double *Pzz_l = &p_ijl[10]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[2]; 
  const double *Pxz_c = &p_ijc[4]; 
  const double *Pyy_c = &p_ijc[6]; 
  const double *Pyz_c = &p_ijc[8]; 
  const double *Pzz_c = &p_ijc[10]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[2]; 
  const double *Pxz_r = &p_ijr[4]; 
  const double *Pyy_r = &p_ijr[6]; 
  const double *Pyz_r = &p_ijr[8]; 
  const double *Pzz_r = &p_ijr[10]; 

  const double *div_p_x_cell_avg_l = &div_p_cell_avgl[0]; 
  const double *div_p_x_cell_avg_c = &div_p_cell_avgc[0]; 
  const double *div_p_x_cell_avg_r = &div_p_cell_avgr[0]; 
  const double *div_p_y_cell_avg_l = &div_p_cell_avgl[1]; 
  const double *div_p_y_cell_avg_c = &div_p_cell_avgc[1]; 
  const double *div_p_y_cell_avg_r = &div_p_cell_avgr[1]; 
  const double *div_p_z_cell_avg_l = &div_p_cell_avgl[2]; 
  const double *div_p_z_cell_avg_c = &div_p_cell_avgc[2]; 
  const double *div_p_z_cell_avg_r = &div_p_cell_avgr[2]; 

  const double *grad_u_x_cell_avg_l = &div_p_cell_avgl[9]; 
  const double *grad_u_x_cell_avg_c = &div_p_cell_avgc[9]; 
  const double *grad_u_x_cell_avg_r = &div_p_cell_avgr[9]; 
  const double *grad_u_y_cell_avg_l = &div_p_cell_avgl[10]; 
  const double *grad_u_y_cell_avg_c = &div_p_cell_avgc[10]; 
  const double *grad_u_y_cell_avg_r = &div_p_cell_avgr[10]; 
  const double *grad_u_z_cell_avg_l = &div_p_cell_avgl[11]; 
  const double *grad_u_z_cell_avg_c = &div_p_cell_avgc[11]; 
  const double *grad_u_z_cell_avg_r = &div_p_cell_avgr[11]; 

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
  const double *rhouy_l = &statevecl[2]; 
  const double *rhouz_l = &statevecl[4]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[2]; 
  const double *rhouz_c = &statevecc[4]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[2]; 
  const double *rhouz_r = &statevecr[4]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[2]; 
  double *div_p_z = &div_p[4]; 

  double *div_b = &pkpm_accel_vars[0]; 
  double *bb_grad_u = &pkpm_accel_vars[2]; 
  double *p_force = &pkpm_accel_vars[4]; 
  double *p_perp_source = &pkpm_accel_vars[6]; 
  double *p_perp_div_b = &pkpm_accel_vars[8]; 

  const double dxHyp = dx1*dx1*dx1*dx1; 

  double grad_u_x[2] = {0.0}; 
  double grad_u_y[2] = {0.0}; 
  double grad_u_z[2] = {0.0}; 
  double div_p_x_comp[2] = {0.0}; 
  double div_p_y_comp[2] = {0.0}; 
  double div_p_z_comp[2] = {0.0}; 
  grad_u_x[1] += -1.732050807568877*ux_c[0]*dx1; 
  div_p_x_comp[1] += -1.732050807568877*Pxx_c[0]*dx1; 
  if (phi_x_l == 0.0) { 
  grad_u_x[0] += ((-0.4330127018922193*ux_l[1])+0.4330127018922193*ux_c[1]-0.25*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += (0.75*ux_l[1]-0.75*ux_c[1]+0.4330127018922193*(ux_l[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += ((-0.4330127018922193*Pxx_l[1])+0.4330127018922193*Pxx_c[1]-0.25*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += (0.75*Pxx_l[1]-0.75*Pxx_c[1]+0.4330127018922193*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  } else { 
  grad_u_x[0] += ((-0.2886751345948129*ux_l[1])+0.2886751345948129*ux_c[1]-0.25*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += (0.5*ux_l[1]-0.5*ux_c[1]+0.4330127018922193*(ux_l[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += ((-0.2886751345948129*Pxx_l[1])+0.2886751345948129*Pxx_c[1]-0.25*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += (0.5*Pxx_l[1]-0.5*Pxx_c[1]+0.4330127018922193*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  } 
  if (phi_x_r == 0.0) { 
  grad_u_x[0] += ((-0.4330127018922193*ux_r[1])+0.4330127018922193*ux_c[1]+0.25*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += ((-0.75*ux_r[1])+0.75*ux_c[1]+0.4330127018922193*(ux_r[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += ((-0.4330127018922193*Pxx_r[1])+0.4330127018922193*Pxx_c[1]+0.25*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += ((-0.75*Pxx_r[1])+0.75*Pxx_c[1]+0.4330127018922193*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  } else { 
  grad_u_x[0] += ((-0.2886751345948129*ux_r[1])+0.2886751345948129*ux_c[1]+0.25*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += ((-0.5*ux_r[1])+0.5*ux_c[1]+0.4330127018922193*(ux_r[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += ((-0.2886751345948129*Pxx_r[1])+0.2886751345948129*Pxx_c[1]+0.25*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += ((-0.5*Pxx_r[1])+0.5*Pxx_c[1]+0.4330127018922193*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  } 

  grad_u_y[1] += -1.732050807568877*uy_c[0]*dx1; 
  div_p_y_comp[1] += -1.732050807568877*Pxy_c[0]*dx1; 
  if (phi_y_l == 0.0) { 
  grad_u_y[0] += ((-0.4330127018922193*uy_l[1])+0.4330127018922193*uy_c[1]-0.25*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += (0.75*uy_l[1]-0.75*uy_c[1]+0.4330127018922193*(uy_l[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += ((-0.4330127018922193*Pxy_l[1])+0.4330127018922193*Pxy_c[1]-0.25*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += (0.75*Pxy_l[1]-0.75*Pxy_c[1]+0.4330127018922193*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  } else { 
  grad_u_y[0] += ((-0.2886751345948129*uy_l[1])+0.2886751345948129*uy_c[1]-0.25*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += (0.5*uy_l[1]-0.5*uy_c[1]+0.4330127018922193*(uy_l[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += ((-0.2886751345948129*Pxy_l[1])+0.2886751345948129*Pxy_c[1]-0.25*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += (0.5*Pxy_l[1]-0.5*Pxy_c[1]+0.4330127018922193*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  } 
  if (phi_y_r == 0.0) { 
  grad_u_y[0] += ((-0.4330127018922193*uy_r[1])+0.4330127018922193*uy_c[1]+0.25*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += ((-0.75*uy_r[1])+0.75*uy_c[1]+0.4330127018922193*(uy_r[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += ((-0.4330127018922193*Pxy_r[1])+0.4330127018922193*Pxy_c[1]+0.25*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += ((-0.75*Pxy_r[1])+0.75*Pxy_c[1]+0.4330127018922193*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  } else { 
  grad_u_y[0] += ((-0.2886751345948129*uy_r[1])+0.2886751345948129*uy_c[1]+0.25*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += ((-0.5*uy_r[1])+0.5*uy_c[1]+0.4330127018922193*(uy_r[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += ((-0.2886751345948129*Pxy_r[1])+0.2886751345948129*Pxy_c[1]+0.25*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += ((-0.5*Pxy_r[1])+0.5*Pxy_c[1]+0.4330127018922193*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  } 

  grad_u_z[1] += -1.732050807568877*uz_c[0]*dx1; 
  div_p_z_comp[1] += -1.732050807568877*Pxz_c[0]*dx1; 
  if (phi_z_l == 0.0) { 
  grad_u_z[0] += ((-0.4330127018922193*uz_l[1])+0.4330127018922193*uz_c[1]-0.25*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += (0.75*uz_l[1]-0.75*uz_c[1]+0.4330127018922193*(uz_l[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += ((-0.4330127018922193*Pxz_l[1])+0.4330127018922193*Pxz_c[1]-0.25*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += (0.75*Pxz_l[1]-0.75*Pxz_c[1]+0.4330127018922193*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  } else { 
  grad_u_z[0] += ((-0.2886751345948129*uz_l[1])+0.2886751345948129*uz_c[1]-0.25*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += (0.5*uz_l[1]-0.5*uz_c[1]+0.4330127018922193*(uz_l[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += ((-0.2886751345948129*Pxz_l[1])+0.2886751345948129*Pxz_c[1]-0.25*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += (0.5*Pxz_l[1]-0.5*Pxz_c[1]+0.4330127018922193*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  } 
  if (phi_z_r == 0.0) { 
  grad_u_z[0] += ((-0.4330127018922193*uz_r[1])+0.4330127018922193*uz_c[1]+0.25*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += ((-0.75*uz_r[1])+0.75*uz_c[1]+0.4330127018922193*(uz_r[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += ((-0.4330127018922193*Pxz_r[1])+0.4330127018922193*Pxz_c[1]+0.25*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += ((-0.75*Pxz_r[1])+0.75*Pxz_c[1]+0.4330127018922193*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  } else { 
  grad_u_z[0] += ((-0.2886751345948129*uz_r[1])+0.2886751345948129*uz_c[1]+0.25*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += ((-0.5*uz_r[1])+0.5*uz_c[1]+0.4330127018922193*(uz_r[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += ((-0.2886751345948129*Pxz_r[1])+0.2886751345948129*Pxz_c[1]+0.25*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += ((-0.5*Pxz_r[1])+0.5*Pxz_c[1]+0.4330127018922193*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  } 

  div_p_x[0] += (4.871392896287466*rhoux_r[1]-4.871392896287466*rhoux_l[1]-2.8125*(rhoux_r[0]+rhoux_l[0])+5.625*rhoux_c[0])*dxHyp*nuHyp+div_p_x_comp[0]; 
  div_p_x[1] += (72.1875*(rhoux_r[1]+rhoux_l[1])+249.375*rhoux_c[1]-56.83291712335378*rhoux_r[0]+56.83291712335378*rhoux_l[0])*dxHyp*nuHyp+div_p_x_comp[1]; 

  div_p_y[0] += (4.871392896287466*rhouy_r[1]-4.871392896287466*rhouy_l[1]-2.8125*(rhouy_r[0]+rhouy_l[0])+5.625*rhouy_c[0])*dxHyp*nuHyp+div_p_y_comp[0]; 
  div_p_y[1] += (72.1875*(rhouy_r[1]+rhouy_l[1])+249.375*rhouy_c[1]-56.83291712335378*rhouy_r[0]+56.83291712335378*rhouy_l[0])*dxHyp*nuHyp+div_p_y_comp[1]; 

  div_p_z[0] += (4.871392896287466*rhouz_r[1]-4.871392896287466*rhouz_l[1]-2.8125*(rhouz_r[0]+rhouz_l[0])+5.625*rhouz_c[0])*dxHyp*nuHyp+div_p_z_comp[0]; 
  div_p_z[1] += (72.1875*(rhouz_r[1]+rhouz_l[1])+249.375*rhouz_c[1]-56.83291712335378*rhouz_r[0]+56.83291712335378*rhouz_l[0])*dxHyp*nuHyp+div_p_z_comp[1]; 

  double div_b_comp[2] = {0.0}; 
  double bb_grad_u_comp[2] = {0.0}; 
  div_b_comp[0] = (-0.2886751345948129*b_r[1]*dx1)-0.2886751345948129*b_l[1]*dx1+0.5773502691896258*b_c[1]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = (-0.5*b_r[1]*dx1)+0.5*b_l[1]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 

  bb_grad_u_comp[0] = 0.7071067811865475*bxbz[1]*grad_u_z[1]+0.7071067811865475*bxby[1]*grad_u_y[1]+0.7071067811865475*bxbx[1]*grad_u_x[1]+0.7071067811865475*bxbz[0]*grad_u_z[0]+0.7071067811865475*bxby[0]*grad_u_y[0]+0.7071067811865475*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.7071067811865475*bxbz[0]*grad_u_z[1]+0.7071067811865475*bxby[0]*grad_u_y[1]+0.7071067811865475*bxbx[0]*grad_u_x[1]+0.7071067811865475*grad_u_z[0]*bxbz[1]+0.7071067811865475*grad_u_y[0]*bxby[1]+0.7071067811865475*grad_u_x[0]*bxbx[1]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 

  p_force[0] += 0.7071067811865475*pkpm_div_ppar_c[1]*rho_inv[1]-0.7071067811865475*T_perp_over_m[1]*div_b_comp[1]+0.7071067811865475*pkpm_div_ppar_c[0]*rho_inv[0]-0.7071067811865475*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += 0.7071067811865475*(pkpm_div_ppar_c[0]*rho_inv[1]+rho_inv[0]*pkpm_div_ppar_c[1])-0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 

  p_perp_source[0] += (-2.0*nu[0])-1.0*grad_u_x[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-2.0*nu[1])-1.0*grad_u_x[1]+bb_grad_u_comp[1]; 

  p_perp_div_b[0] += 0.7071067811865475*(T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 

} 
