#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p2(const double *dxv, double nuHyp, 
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
  const double *uy_l = &u_il[3]; 
  const double *uz_l = &u_il[6]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[3]; 
  const double *uz_c = &u_ic[6]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[3]; 
  const double *uz_r = &u_ir[6]; 

  const double *bxbx = &bvarc[9]; 
  const double *bxby = &bvarc[12]; 
  const double *bxbz = &bvarc[15]; 
  const double *byby = &bvarc[18]; 
  const double *bybz = &bvarc[21]; 
  const double *bzbz = &bvarc[24]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[3]; 
  const double *Pxz_l = &p_ijl[6]; 
  const double *Pyy_l = &p_ijl[9]; 
  const double *Pyz_l = &p_ijl[12]; 
  const double *Pzz_l = &p_ijl[15]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[3]; 
  const double *Pxz_c = &p_ijc[6]; 
  const double *Pyy_c = &p_ijc[9]; 
  const double *Pyz_c = &p_ijc[12]; 
  const double *Pzz_c = &p_ijc[15]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[3]; 
  const double *Pxz_r = &p_ijr[6]; 
  const double *Pyy_r = &p_ijr[9]; 
  const double *Pyz_r = &p_ijr[12]; 
  const double *Pzz_r = &p_ijr[15]; 

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
  const double *rhouy_l = &statevecl[3]; 
  const double *rhouz_l = &statevecl[6]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[3]; 
  const double *rhouz_c = &statevecc[6]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[3]; 
  const double *rhouz_r = &statevecr[6]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[3]; 
  double *div_p_z = &div_p[6]; 

  double *div_b = &pkpm_accel_vars[0]; 
  double *bb_grad_u = &pkpm_accel_vars[3]; 
  double *p_force = &pkpm_accel_vars[6]; 
  double *p_perp_source = &pkpm_accel_vars[9]; 
  double *p_perp_div_b = &pkpm_accel_vars[12]; 

  const double dxHyp = dx1*dx1*dx1*dx1; 

  double grad_u_x[3] = {0.0}; 
  double grad_u_y[3] = {0.0}; 
  double grad_u_z[3] = {0.0}; 
  double div_p_x_comp[3] = {0.0}; 
  double div_p_y_comp[3] = {0.0}; 
  double div_p_z_comp[3] = {0.0}; 
  grad_u_x[1] += -1.732050807568877*ux_c[0]*dx1; 
  grad_u_x[2] += -3.872983346207417*ux_c[1]*dx1; 
  div_p_x_comp[1] += -1.732050807568877*Pxx_c[0]*dx1; 
  div_p_x_comp[2] += -3.872983346207417*Pxx_c[1]*dx1; 
  if (phi_x_l == 0.0) { 
  grad_u_x[0] += ((-0.5590169943749475*(ux_l[2]+ux_c[2]))-0.4330127018922193*ux_l[1]+0.4330127018922193*ux_c[1]-0.25*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += (0.9682458365518543*(ux_l[2]+ux_c[2])+0.75*ux_l[1]-0.75*ux_c[1]+0.4330127018922193*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[2] += ((-1.25*(ux_l[2]+ux_c[2]))-0.9682458365518543*ux_l[1]+0.9682458365518543*ux_c[1]-0.5590169943749475*(ux_l[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += ((-0.5590169943749475*(Pxx_l[2]+Pxx_c[2]))-0.4330127018922193*Pxx_l[1]+0.4330127018922193*Pxx_c[1]-0.25*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += (0.9682458365518543*(Pxx_l[2]+Pxx_c[2])+0.75*Pxx_l[1]-0.75*Pxx_c[1]+0.4330127018922193*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[2] += ((-1.25*(Pxx_l[2]+Pxx_c[2]))-0.9682458365518543*Pxx_l[1]+0.9682458365518543*Pxx_c[1]-0.5590169943749475*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  } else { 
  grad_u_x[0] += ((-0.2445699350390395*(ux_l[2]+ux_c[2]))-0.3518228202874282*ux_l[1]+0.3518228202874282*ux_c[1]-0.25*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += (0.4236075534914363*(ux_l[2]+ux_c[2])+0.609375*ux_l[1]-0.609375*ux_c[1]+0.4330127018922193*(ux_l[0]+ux_c[0]))*dx1; 
  grad_u_x[2] += ((-0.546875*(ux_l[2]+ux_c[2]))-0.7866997421983816*ux_l[1]+0.7866997421983816*ux_c[1]-0.5590169943749475*(ux_l[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += ((-0.2445699350390395*(Pxx_l[2]+Pxx_c[2]))-0.3518228202874282*Pxx_l[1]+0.3518228202874282*Pxx_c[1]-0.25*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += (0.4236075534914363*(Pxx_l[2]+Pxx_c[2])+0.609375*Pxx_l[1]-0.609375*Pxx_c[1]+0.4330127018922193*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[2] += ((-0.546875*(Pxx_l[2]+Pxx_c[2]))-0.7866997421983816*Pxx_l[1]+0.7866997421983816*Pxx_c[1]-0.5590169943749475*(Pxx_l[0]+Pxx_c[0]))*dx1; 
  } 
  if (phi_x_r == 0.0) { 
  grad_u_x[0] += (0.5590169943749475*(ux_r[2]+ux_c[2])-0.4330127018922193*ux_r[1]+0.4330127018922193*ux_c[1]+0.25*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += (0.9682458365518543*(ux_r[2]+ux_c[2])-0.75*ux_r[1]+0.75*ux_c[1]+0.4330127018922193*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[2] += (1.25*(ux_r[2]+ux_c[2])-0.9682458365518543*ux_r[1]+0.9682458365518543*ux_c[1]+0.5590169943749475*(ux_r[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += (0.5590169943749475*(Pxx_r[2]+Pxx_c[2])-0.4330127018922193*Pxx_r[1]+0.4330127018922193*Pxx_c[1]+0.25*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += (0.9682458365518543*(Pxx_r[2]+Pxx_c[2])-0.75*Pxx_r[1]+0.75*Pxx_c[1]+0.4330127018922193*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[2] += (1.25*(Pxx_r[2]+Pxx_c[2])-0.9682458365518543*Pxx_r[1]+0.9682458365518543*Pxx_c[1]+0.5590169943749475*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  } else { 
  grad_u_x[0] += (0.2445699350390395*(ux_r[2]+ux_c[2])-0.3518228202874282*ux_r[1]+0.3518228202874282*ux_c[1]+0.25*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[1] += (0.4236075534914363*(ux_r[2]+ux_c[2])-0.609375*ux_r[1]+0.609375*ux_c[1]+0.4330127018922193*(ux_r[0]+ux_c[0]))*dx1; 
  grad_u_x[2] += (0.546875*(ux_r[2]+ux_c[2])-0.7866997421983816*ux_r[1]+0.7866997421983816*ux_c[1]+0.5590169943749475*(ux_r[0]+ux_c[0]))*dx1; 
  div_p_x_comp[0] += (0.2445699350390395*(Pxx_r[2]+Pxx_c[2])-0.3518228202874282*Pxx_r[1]+0.3518228202874282*Pxx_c[1]+0.25*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[1] += (0.4236075534914363*(Pxx_r[2]+Pxx_c[2])-0.609375*Pxx_r[1]+0.609375*Pxx_c[1]+0.4330127018922193*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  div_p_x_comp[2] += (0.546875*(Pxx_r[2]+Pxx_c[2])-0.7866997421983816*Pxx_r[1]+0.7866997421983816*Pxx_c[1]+0.5590169943749475*(Pxx_r[0]+Pxx_c[0]))*dx1; 
  } 

  grad_u_y[1] += -1.732050807568877*uy_c[0]*dx1; 
  grad_u_y[2] += -3.872983346207417*uy_c[1]*dx1; 
  div_p_y_comp[1] += -1.732050807568877*Pxy_c[0]*dx1; 
  div_p_y_comp[2] += -3.872983346207417*Pxy_c[1]*dx1; 
  if (phi_y_l == 0.0) { 
  grad_u_y[0] += ((-0.5590169943749475*(uy_l[2]+uy_c[2]))-0.4330127018922193*uy_l[1]+0.4330127018922193*uy_c[1]-0.25*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += (0.9682458365518543*(uy_l[2]+uy_c[2])+0.75*uy_l[1]-0.75*uy_c[1]+0.4330127018922193*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[2] += ((-1.25*(uy_l[2]+uy_c[2]))-0.9682458365518543*uy_l[1]+0.9682458365518543*uy_c[1]-0.5590169943749475*(uy_l[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += ((-0.5590169943749475*(Pxy_l[2]+Pxy_c[2]))-0.4330127018922193*Pxy_l[1]+0.4330127018922193*Pxy_c[1]-0.25*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += (0.9682458365518543*(Pxy_l[2]+Pxy_c[2])+0.75*Pxy_l[1]-0.75*Pxy_c[1]+0.4330127018922193*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[2] += ((-1.25*(Pxy_l[2]+Pxy_c[2]))-0.9682458365518543*Pxy_l[1]+0.9682458365518543*Pxy_c[1]-0.5590169943749475*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  } else { 
  grad_u_y[0] += ((-0.2445699350390395*(uy_l[2]+uy_c[2]))-0.3518228202874282*uy_l[1]+0.3518228202874282*uy_c[1]-0.25*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += (0.4236075534914363*(uy_l[2]+uy_c[2])+0.609375*uy_l[1]-0.609375*uy_c[1]+0.4330127018922193*(uy_l[0]+uy_c[0]))*dx1; 
  grad_u_y[2] += ((-0.546875*(uy_l[2]+uy_c[2]))-0.7866997421983816*uy_l[1]+0.7866997421983816*uy_c[1]-0.5590169943749475*(uy_l[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += ((-0.2445699350390395*(Pxy_l[2]+Pxy_c[2]))-0.3518228202874282*Pxy_l[1]+0.3518228202874282*Pxy_c[1]-0.25*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += (0.4236075534914363*(Pxy_l[2]+Pxy_c[2])+0.609375*Pxy_l[1]-0.609375*Pxy_c[1]+0.4330127018922193*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[2] += ((-0.546875*(Pxy_l[2]+Pxy_c[2]))-0.7866997421983816*Pxy_l[1]+0.7866997421983816*Pxy_c[1]-0.5590169943749475*(Pxy_l[0]+Pxy_c[0]))*dx1; 
  } 
  if (phi_y_r == 0.0) { 
  grad_u_y[0] += (0.5590169943749475*(uy_r[2]+uy_c[2])-0.4330127018922193*uy_r[1]+0.4330127018922193*uy_c[1]+0.25*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += (0.9682458365518543*(uy_r[2]+uy_c[2])-0.75*uy_r[1]+0.75*uy_c[1]+0.4330127018922193*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[2] += (1.25*(uy_r[2]+uy_c[2])-0.9682458365518543*uy_r[1]+0.9682458365518543*uy_c[1]+0.5590169943749475*(uy_r[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += (0.5590169943749475*(Pxy_r[2]+Pxy_c[2])-0.4330127018922193*Pxy_r[1]+0.4330127018922193*Pxy_c[1]+0.25*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += (0.9682458365518543*(Pxy_r[2]+Pxy_c[2])-0.75*Pxy_r[1]+0.75*Pxy_c[1]+0.4330127018922193*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[2] += (1.25*(Pxy_r[2]+Pxy_c[2])-0.9682458365518543*Pxy_r[1]+0.9682458365518543*Pxy_c[1]+0.5590169943749475*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  } else { 
  grad_u_y[0] += (0.2445699350390395*(uy_r[2]+uy_c[2])-0.3518228202874282*uy_r[1]+0.3518228202874282*uy_c[1]+0.25*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[1] += (0.4236075534914363*(uy_r[2]+uy_c[2])-0.609375*uy_r[1]+0.609375*uy_c[1]+0.4330127018922193*(uy_r[0]+uy_c[0]))*dx1; 
  grad_u_y[2] += (0.546875*(uy_r[2]+uy_c[2])-0.7866997421983816*uy_r[1]+0.7866997421983816*uy_c[1]+0.5590169943749475*(uy_r[0]+uy_c[0]))*dx1; 
  div_p_y_comp[0] += (0.2445699350390395*(Pxy_r[2]+Pxy_c[2])-0.3518228202874282*Pxy_r[1]+0.3518228202874282*Pxy_c[1]+0.25*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[1] += (0.4236075534914363*(Pxy_r[2]+Pxy_c[2])-0.609375*Pxy_r[1]+0.609375*Pxy_c[1]+0.4330127018922193*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  div_p_y_comp[2] += (0.546875*(Pxy_r[2]+Pxy_c[2])-0.7866997421983816*Pxy_r[1]+0.7866997421983816*Pxy_c[1]+0.5590169943749475*(Pxy_r[0]+Pxy_c[0]))*dx1; 
  } 

  grad_u_z[1] += -1.732050807568877*uz_c[0]*dx1; 
  grad_u_z[2] += -3.872983346207417*uz_c[1]*dx1; 
  div_p_z_comp[1] += -1.732050807568877*Pxz_c[0]*dx1; 
  div_p_z_comp[2] += -3.872983346207417*Pxz_c[1]*dx1; 
  if (phi_z_l == 0.0) { 
  grad_u_z[0] += ((-0.5590169943749475*(uz_l[2]+uz_c[2]))-0.4330127018922193*uz_l[1]+0.4330127018922193*uz_c[1]-0.25*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += (0.9682458365518543*(uz_l[2]+uz_c[2])+0.75*uz_l[1]-0.75*uz_c[1]+0.4330127018922193*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[2] += ((-1.25*(uz_l[2]+uz_c[2]))-0.9682458365518543*uz_l[1]+0.9682458365518543*uz_c[1]-0.5590169943749475*(uz_l[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += ((-0.5590169943749475*(Pxz_l[2]+Pxz_c[2]))-0.4330127018922193*Pxz_l[1]+0.4330127018922193*Pxz_c[1]-0.25*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += (0.9682458365518543*(Pxz_l[2]+Pxz_c[2])+0.75*Pxz_l[1]-0.75*Pxz_c[1]+0.4330127018922193*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[2] += ((-1.25*(Pxz_l[2]+Pxz_c[2]))-0.9682458365518543*Pxz_l[1]+0.9682458365518543*Pxz_c[1]-0.5590169943749475*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  } else { 
  grad_u_z[0] += ((-0.2445699350390395*(uz_l[2]+uz_c[2]))-0.3518228202874282*uz_l[1]+0.3518228202874282*uz_c[1]-0.25*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += (0.4236075534914363*(uz_l[2]+uz_c[2])+0.609375*uz_l[1]-0.609375*uz_c[1]+0.4330127018922193*(uz_l[0]+uz_c[0]))*dx1; 
  grad_u_z[2] += ((-0.546875*(uz_l[2]+uz_c[2]))-0.7866997421983816*uz_l[1]+0.7866997421983816*uz_c[1]-0.5590169943749475*(uz_l[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += ((-0.2445699350390395*(Pxz_l[2]+Pxz_c[2]))-0.3518228202874282*Pxz_l[1]+0.3518228202874282*Pxz_c[1]-0.25*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += (0.4236075534914363*(Pxz_l[2]+Pxz_c[2])+0.609375*Pxz_l[1]-0.609375*Pxz_c[1]+0.4330127018922193*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[2] += ((-0.546875*(Pxz_l[2]+Pxz_c[2]))-0.7866997421983816*Pxz_l[1]+0.7866997421983816*Pxz_c[1]-0.5590169943749475*(Pxz_l[0]+Pxz_c[0]))*dx1; 
  } 
  if (phi_z_r == 0.0) { 
  grad_u_z[0] += (0.5590169943749475*(uz_r[2]+uz_c[2])-0.4330127018922193*uz_r[1]+0.4330127018922193*uz_c[1]+0.25*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += (0.9682458365518543*(uz_r[2]+uz_c[2])-0.75*uz_r[1]+0.75*uz_c[1]+0.4330127018922193*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[2] += (1.25*(uz_r[2]+uz_c[2])-0.9682458365518543*uz_r[1]+0.9682458365518543*uz_c[1]+0.5590169943749475*(uz_r[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += (0.5590169943749475*(Pxz_r[2]+Pxz_c[2])-0.4330127018922193*Pxz_r[1]+0.4330127018922193*Pxz_c[1]+0.25*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += (0.9682458365518543*(Pxz_r[2]+Pxz_c[2])-0.75*Pxz_r[1]+0.75*Pxz_c[1]+0.4330127018922193*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[2] += (1.25*(Pxz_r[2]+Pxz_c[2])-0.9682458365518543*Pxz_r[1]+0.9682458365518543*Pxz_c[1]+0.5590169943749475*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  } else { 
  grad_u_z[0] += (0.2445699350390395*(uz_r[2]+uz_c[2])-0.3518228202874282*uz_r[1]+0.3518228202874282*uz_c[1]+0.25*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[1] += (0.4236075534914363*(uz_r[2]+uz_c[2])-0.609375*uz_r[1]+0.609375*uz_c[1]+0.4330127018922193*(uz_r[0]+uz_c[0]))*dx1; 
  grad_u_z[2] += (0.546875*(uz_r[2]+uz_c[2])-0.7866997421983816*uz_r[1]+0.7866997421983816*uz_c[1]+0.5590169943749475*(uz_r[0]+uz_c[0]))*dx1; 
  div_p_z_comp[0] += (0.2445699350390395*(Pxz_r[2]+Pxz_c[2])-0.3518228202874282*Pxz_r[1]+0.3518228202874282*Pxz_c[1]+0.25*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[1] += (0.4236075534914363*(Pxz_r[2]+Pxz_c[2])-0.609375*Pxz_r[1]+0.609375*Pxz_c[1]+0.4330127018922193*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  div_p_z_comp[2] += (0.546875*(Pxz_r[2]+Pxz_c[2])-0.7866997421983816*Pxz_r[1]+0.7866997421983816*Pxz_c[1]+0.5590169943749475*(Pxz_r[0]+Pxz_c[0]))*dx1; 
  } 

  div_p_x[0] += ((-160.6824473206489*(rhoux_r[2]+rhoux_l[2]))+805.6133660185961*rhoux_c[2]+207.4401475002413*rhoux_r[1]-207.4401475002413*rhoux_l[1]-137.8125*(rhoux_r[0]+rhoux_l[0])+275.625*rhoux_c[0])*dxHyp*nuHyp+div_p_x_comp[0]; 
  div_p_x[1] += ((-48.29126109802372*rhoux_r[2])+48.29126109802372*rhoux_l[2]+55.78125*(rhoux_r[1]+rhoux_l[1])+124.6875*rhoux_c[1]-34.09975027401226*rhoux_r[0]+34.09975027401226*rhoux_l[0])*dxHyp*nuHyp+div_p_x_comp[1]; 
  div_p_x[2] += ((-250.03125*(rhoux_r[2]+rhoux_l[2]))+1326.9375*rhoux_c[2]+327.8722464023716*rhoux_r[1]-327.8722464023716*rhoux_l[1]-220.1129415351356*(rhoux_r[0]+rhoux_l[0])+440.2258830702712*rhoux_c[0])*dxHyp*nuHyp+div_p_x_comp[2]; 

  div_p_y[0] += ((-160.6824473206489*(rhouy_r[2]+rhouy_l[2]))+805.6133660185961*rhouy_c[2]+207.4401475002413*rhouy_r[1]-207.4401475002413*rhouy_l[1]-137.8125*(rhouy_r[0]+rhouy_l[0])+275.625*rhouy_c[0])*dxHyp*nuHyp+div_p_y_comp[0]; 
  div_p_y[1] += ((-48.29126109802372*rhouy_r[2])+48.29126109802372*rhouy_l[2]+55.78125*(rhouy_r[1]+rhouy_l[1])+124.6875*rhouy_c[1]-34.09975027401226*rhouy_r[0]+34.09975027401226*rhouy_l[0])*dxHyp*nuHyp+div_p_y_comp[1]; 
  div_p_y[2] += ((-250.03125*(rhouy_r[2]+rhouy_l[2]))+1326.9375*rhouy_c[2]+327.8722464023716*rhouy_r[1]-327.8722464023716*rhouy_l[1]-220.1129415351356*(rhouy_r[0]+rhouy_l[0])+440.2258830702712*rhouy_c[0])*dxHyp*nuHyp+div_p_y_comp[2]; 

  div_p_z[0] += ((-160.6824473206489*(rhouz_r[2]+rhouz_l[2]))+805.6133660185961*rhouz_c[2]+207.4401475002413*rhouz_r[1]-207.4401475002413*rhouz_l[1]-137.8125*(rhouz_r[0]+rhouz_l[0])+275.625*rhouz_c[0])*dxHyp*nuHyp+div_p_z_comp[0]; 
  div_p_z[1] += ((-48.29126109802372*rhouz_r[2])+48.29126109802372*rhouz_l[2]+55.78125*(rhouz_r[1]+rhouz_l[1])+124.6875*rhouz_c[1]-34.09975027401226*rhouz_r[0]+34.09975027401226*rhouz_l[0])*dxHyp*nuHyp+div_p_z_comp[1]; 
  div_p_z[2] += ((-250.03125*(rhouz_r[2]+rhouz_l[2]))+1326.9375*rhouz_c[2]+327.8722464023716*rhouz_r[1]-327.8722464023716*rhouz_l[1]-220.1129415351356*(rhouz_r[0]+rhouz_l[0])+440.2258830702712*rhouz_c[0])*dxHyp*nuHyp+div_p_z_comp[2]; 

  double div_b_comp[3] = {0.0}; 
  double bb_grad_u_comp[3] = {0.0}; 
  div_b_comp[0] = 0.2445699350390395*b_r[2]*dx1-0.2445699350390395*b_l[2]*dx1-0.3518228202874282*b_r[1]*dx1-0.3518228202874282*b_l[1]*dx1+0.7036456405748563*b_c[1]*dx1+0.25*b_r[0]*dx1-0.25*b_l[0]*dx1; 
  div_b_comp[1] = 0.4236075534914363*b_r[2]*dx1+0.4236075534914363*b_l[2]*dx1+0.8472151069828725*b_c[2]*dx1-0.609375*b_r[1]*dx1+0.609375*b_l[1]*dx1+0.4330127018922193*b_r[0]*dx1+0.4330127018922193*b_l[0]*dx1-0.8660254037844386*b_c[0]*dx1; 
  div_b_comp[2] = 0.546875*b_r[2]*dx1-0.546875*b_l[2]*dx1-0.7866997421983816*b_r[1]*dx1-0.7866997421983816*b_l[1]*dx1-2.299583861810654*b_c[1]*dx1+0.5590169943749475*b_r[0]*dx1-0.5590169943749475*b_l[0]*dx1; 

  div_b[0] += div_b_comp[0]; 
  div_b[1] += div_b_comp[1]; 
  div_b[2] += div_b_comp[2]; 

  bb_grad_u_comp[0] = 0.7071067811865475*bxbz[2]*grad_u_z[2]+0.7071067811865475*bxby[2]*grad_u_y[2]+0.7071067811865475*bxbx[2]*grad_u_x[2]+0.7071067811865475*bxbz[1]*grad_u_z[1]+0.7071067811865475*bxby[1]*grad_u_y[1]+0.7071067811865475*bxbx[1]*grad_u_x[1]+0.7071067811865475*bxbz[0]*grad_u_z[0]+0.7071067811865475*bxby[0]*grad_u_y[0]+0.7071067811865475*bxbx[0]*grad_u_x[0]; 
  bb_grad_u_comp[1] = 0.6324555320336759*bxbz[1]*grad_u_z[2]+0.6324555320336759*bxby[1]*grad_u_y[2]+0.6324555320336759*bxbx[1]*grad_u_x[2]+0.6324555320336759*grad_u_z[1]*bxbz[2]+0.6324555320336759*grad_u_y[1]*bxby[2]+0.6324555320336759*grad_u_x[1]*bxbx[2]+0.7071067811865475*bxbz[0]*grad_u_z[1]+0.7071067811865475*bxby[0]*grad_u_y[1]+0.7071067811865475*bxbx[0]*grad_u_x[1]+0.7071067811865475*grad_u_z[0]*bxbz[1]+0.7071067811865475*grad_u_y[0]*bxby[1]+0.7071067811865475*grad_u_x[0]*bxbx[1]; 
  bb_grad_u_comp[2] = 0.4517539514526256*bxbz[2]*grad_u_z[2]+0.7071067811865475*bxbz[0]*grad_u_z[2]+0.4517539514526256*bxby[2]*grad_u_y[2]+0.7071067811865475*bxby[0]*grad_u_y[2]+0.4517539514526256*bxbx[2]*grad_u_x[2]+0.7071067811865475*bxbx[0]*grad_u_x[2]+0.7071067811865475*grad_u_z[0]*bxbz[2]+0.7071067811865475*grad_u_y[0]*bxby[2]+0.7071067811865475*grad_u_x[0]*bxbx[2]+0.6324555320336759*bxbz[1]*grad_u_z[1]+0.6324555320336759*bxby[1]*grad_u_y[1]+0.6324555320336759*bxbx[1]*grad_u_x[1]; 

  bb_grad_u[0] += bb_grad_u_comp[0]; 
  bb_grad_u[1] += bb_grad_u_comp[1]; 
  bb_grad_u[2] += bb_grad_u_comp[2]; 

  p_force[0] += 0.7071067811865475*pkpm_div_ppar_c[2]*rho_inv[2]-0.7071067811865475*T_perp_over_m[2]*div_b_comp[2]+0.7071067811865475*pkpm_div_ppar_c[1]*rho_inv[1]-0.7071067811865475*T_perp_over_m[1]*div_b_comp[1]+0.7071067811865475*pkpm_div_ppar_c[0]*rho_inv[0]-0.7071067811865475*T_perp_over_m[0]*div_b_comp[0]; 
  p_force[1] += 0.6324555320336759*(pkpm_div_ppar_c[1]*rho_inv[2]+rho_inv[1]*pkpm_div_ppar_c[2])-0.6324555320336759*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2])+0.7071067811865475*(pkpm_div_ppar_c[0]*rho_inv[1]+rho_inv[0]*pkpm_div_ppar_c[1])-0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force[2] += 0.4517539514526256*pkpm_div_ppar_c[2]*rho_inv[2]+0.7071067811865475*(pkpm_div_ppar_c[0]*rho_inv[2]+rho_inv[0]*pkpm_div_ppar_c[2])-0.4517539514526256*T_perp_over_m[2]*div_b_comp[2]-0.7071067811865475*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2])+0.6324555320336759*pkpm_div_ppar_c[1]*rho_inv[1]-0.6324555320336759*T_perp_over_m[1]*div_b_comp[1]; 

  p_perp_source[0] += (-2.0*nu[0])-1.0*grad_u_x[0]+bb_grad_u_comp[0]; 
  p_perp_source[1] += (-2.0*nu[1])-1.0*grad_u_x[1]+bb_grad_u_comp[1]; 
  p_perp_source[2] += (-2.0*nu[2])-1.0*grad_u_x[2]+bb_grad_u_comp[2]; 

  p_perp_div_b[0] += 0.7071067811865475*(T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.6324555320336759*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2])+0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.4517539514526256*T_perp_over_m[2]*div_b_comp[2]+0.7071067811865475*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2])+0.6324555320336759*T_perp_over_m[1]*div_b_comp[1]; 

} 
