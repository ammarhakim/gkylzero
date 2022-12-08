#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void euler_pkpm_recovery_z_3x_ser_p1(const double *dxv, double nuHyp, 
  const double *bvarl, const double *bvarc, const double *bvarr, 
  const double *u_il, const double *u_ic, const double *u_ir, 
  const double *p_ijl, const double *p_ijc, const double *p_ijr, 
  const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
  const double *statevecl, const double *statevecc, const double *statevecr, 
  double* div_b, double* bb_grad_u, double* div_p, double* p_force, double* p_perp_source) 
{ 
  // dxv[NDIM]: Cell spacing.
  // nuHyp: Hyper-diffusion coefficient.
  // bvarl/bvarc/bvarr:  Input magnetic field unit vector in left/center/right cells.
  // u_il/u_ic/u_ir:     Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // p_ijl/p_ijc/p_ijr:  Input pressure tensor in left/center/right cells.
  // vlasov_pkpm_momsl/vlasov_pkpm_momsc/vlasov_pkpm_momsr: Input pkpm moments (rho, p_parallel, q_parallel) in left/center/right cells.
  // statevecl/statevecc/statevecr:          [rho ux, rho uy, rho uz, p_perp], Fluid input state vector in center cell.
  // div_b:              Increment to volume expansion of div(b) in one direction.
  // bb_grad_u:          Increment to volume expansion of bb : grad(u) in one direction.
  // div_p:              Increment to volume expansion of div(p) in one direction.
  // p_force:            Increment to volume expansion of p_force = 1/rho * div(p_parallel b_hat) in one direction.
  // p_perp_source:      Increment to volume expansion of perpendicular pressure compression source (p_perp div(u) - p_perp bb : grad(u)).

  const double dx1 = 2.0/dxv[2]; 
  const double dx14 = dx1*dx1*dx1*dx1; 

  const double *b_l = &bvarl[16]; 
  const double *b_c = &bvarc[16]; 
  const double *b_r = &bvarr[16]; 

  const double *ux_l = &u_il[0]; 
  const double *uy_l = &u_il[8]; 
  const double *uz_l = &u_il[16]; 

  const double *ux_c = &u_ic[0]; 
  const double *uy_c = &u_ic[8]; 
  const double *uz_c = &u_ic[16]; 

  const double *ux_r = &u_ir[0]; 
  const double *uy_r = &u_ir[8]; 
  const double *uz_r = &u_ir[16]; 

  const double *bxbx = &bvarc[24]; 
  const double *bxby = &bvarc[32]; 
  const double *bxbz = &bvarc[40]; 
  const double *byby = &bvarc[48]; 
  const double *bybz = &bvarc[56]; 
  const double *bzbz = &bvarc[64]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[8]; 
  const double *Pxz_l = &p_ijl[16]; 
  const double *Pyy_l = &p_ijl[24]; 
  const double *Pyz_l = &p_ijl[32]; 
  const double *Pzz_l = &p_ijl[40]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[8]; 
  const double *Pxz_c = &p_ijc[16]; 
  const double *Pyy_c = &p_ijc[24]; 
  const double *Pyz_c = &p_ijc[32]; 
  const double *Pzz_c = &p_ijc[40]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[8]; 
  const double *Pxz_r = &p_ijr[16]; 
  const double *Pyy_r = &p_ijr[24]; 
  const double *Pyz_r = &p_ijr[32]; 
  const double *Pzz_r = &p_ijr[40]; 

  const double *ppar_l = &vlasov_pkpm_momsl[8]; 
  const double *ppar_c = &vlasov_pkpm_momsc[8]; 
  const double *ppar_r = &vlasov_pkpm_momsr[8]; 
  const double *rho = &vlasov_pkpm_momsc[0]; 

  const double *rhoux_l = &statevecl[0]; 
  const double *rhouy_l = &statevecl[8]; 
  const double *rhouz_l = &statevecl[16]; 

  const double *rhoux_c = &statevecc[0]; 
  const double *rhouy_c = &statevecc[8]; 
  const double *rhouz_c = &statevecc[16]; 

  const double *rhoux_r = &statevecr[0]; 
  const double *rhouy_r = &statevecr[8]; 
  const double *rhouz_r = &statevecr[16]; 

  const double *p_perp = &statevecc[24]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[8]; 
  double *div_p_z = &div_p[16]; 

  double grad_u_x[8] = {0.0}; 
  double grad_u_y[8] = {0.0}; 
  double grad_u_z[8] = {0.0}; 
  grad_u_x[0] = (-0.2886751345948129*ux_r[3])-0.2886751345948129*ux_l[3]+0.5773502691896258*ux_c[3]+0.25*ux_r[0]-0.25*ux_l[0]; 
  grad_u_x[1] = (-0.2886751345948129*ux_r[5])-0.2886751345948129*ux_l[5]+0.5773502691896258*ux_c[5]+0.25*ux_r[1]-0.25*ux_l[1]; 
  grad_u_x[2] = (-0.2886751345948129*ux_r[6])-0.2886751345948129*ux_l[6]+0.5773502691896258*ux_c[6]+0.25*ux_r[2]-0.25*ux_l[2]; 
  grad_u_x[3] = (-0.5*ux_r[3])+0.5*ux_l[3]+0.4330127018922193*ux_r[0]+0.4330127018922193*ux_l[0]-0.8660254037844386*ux_c[0]; 
  grad_u_x[4] = (-0.2886751345948129*ux_r[7])-0.2886751345948129*ux_l[7]+0.5773502691896258*ux_c[7]+0.25*ux_r[4]-0.25*ux_l[4]; 
  grad_u_x[5] = (-0.5*ux_r[5])+0.5*ux_l[5]+0.4330127018922193*ux_r[1]+0.4330127018922193*ux_l[1]-0.8660254037844386*ux_c[1]; 
  grad_u_x[6] = (-0.5*ux_r[6])+0.5*ux_l[6]+0.4330127018922193*ux_r[2]+0.4330127018922193*ux_l[2]-0.8660254037844386*ux_c[2]; 
  grad_u_x[7] = (-0.5*ux_r[7])+0.5*ux_l[7]+0.4330127018922193*ux_r[4]+0.4330127018922193*ux_l[4]-0.8660254037844386*ux_c[4]; 

  grad_u_y[0] = (-0.2886751345948129*uy_r[3])-0.2886751345948129*uy_l[3]+0.5773502691896258*uy_c[3]+0.25*uy_r[0]-0.25*uy_l[0]; 
  grad_u_y[1] = (-0.2886751345948129*uy_r[5])-0.2886751345948129*uy_l[5]+0.5773502691896258*uy_c[5]+0.25*uy_r[1]-0.25*uy_l[1]; 
  grad_u_y[2] = (-0.2886751345948129*uy_r[6])-0.2886751345948129*uy_l[6]+0.5773502691896258*uy_c[6]+0.25*uy_r[2]-0.25*uy_l[2]; 
  grad_u_y[3] = (-0.5*uy_r[3])+0.5*uy_l[3]+0.4330127018922193*uy_r[0]+0.4330127018922193*uy_l[0]-0.8660254037844386*uy_c[0]; 
  grad_u_y[4] = (-0.2886751345948129*uy_r[7])-0.2886751345948129*uy_l[7]+0.5773502691896258*uy_c[7]+0.25*uy_r[4]-0.25*uy_l[4]; 
  grad_u_y[5] = (-0.5*uy_r[5])+0.5*uy_l[5]+0.4330127018922193*uy_r[1]+0.4330127018922193*uy_l[1]-0.8660254037844386*uy_c[1]; 
  grad_u_y[6] = (-0.5*uy_r[6])+0.5*uy_l[6]+0.4330127018922193*uy_r[2]+0.4330127018922193*uy_l[2]-0.8660254037844386*uy_c[2]; 
  grad_u_y[7] = (-0.5*uy_r[7])+0.5*uy_l[7]+0.4330127018922193*uy_r[4]+0.4330127018922193*uy_l[4]-0.8660254037844386*uy_c[4]; 

  grad_u_z[0] = (-0.2886751345948129*uz_r[3])-0.2886751345948129*uz_l[3]+0.5773502691896258*uz_c[3]+0.25*uz_r[0]-0.25*uz_l[0]; 
  grad_u_z[1] = (-0.2886751345948129*uz_r[5])-0.2886751345948129*uz_l[5]+0.5773502691896258*uz_c[5]+0.25*uz_r[1]-0.25*uz_l[1]; 
  grad_u_z[2] = (-0.2886751345948129*uz_r[6])-0.2886751345948129*uz_l[6]+0.5773502691896258*uz_c[6]+0.25*uz_r[2]-0.25*uz_l[2]; 
  grad_u_z[3] = (-0.5*uz_r[3])+0.5*uz_l[3]+0.4330127018922193*uz_r[0]+0.4330127018922193*uz_l[0]-0.8660254037844386*uz_c[0]; 
  grad_u_z[4] = (-0.2886751345948129*uz_r[7])-0.2886751345948129*uz_l[7]+0.5773502691896258*uz_c[7]+0.25*uz_r[4]-0.25*uz_l[4]; 
  grad_u_z[5] = (-0.5*uz_r[5])+0.5*uz_l[5]+0.4330127018922193*uz_r[1]+0.4330127018922193*uz_l[1]-0.8660254037844386*uz_c[1]; 
  grad_u_z[6] = (-0.5*uz_r[6])+0.5*uz_l[6]+0.4330127018922193*uz_r[2]+0.4330127018922193*uz_l[2]-0.8660254037844386*uz_c[2]; 
  grad_u_z[7] = (-0.5*uz_r[7])+0.5*uz_l[7]+0.4330127018922193*uz_r[4]+0.4330127018922193*uz_l[4]-0.8660254037844386*uz_c[4]; 

  double ppar_b_l[8] = {0.0}; 
  double ppar_b_c[8] = {0.0}; 
  double ppar_b_r[8] = {0.0}; 
  ppar_b_l[0] = 0.3535533905932737*b_l[7]*ppar_l[7]+0.3535533905932737*b_l[6]*ppar_l[6]+0.3535533905932737*b_l[5]*ppar_l[5]+0.3535533905932737*b_l[4]*ppar_l[4]+0.3535533905932737*b_l[3]*ppar_l[3]+0.3535533905932737*b_l[2]*ppar_l[2]+0.3535533905932737*b_l[1]*ppar_l[1]+0.3535533905932737*b_l[0]*ppar_l[0]; 
  ppar_b_l[1] = 0.3535533905932737*b_l[6]*ppar_l[7]+0.3535533905932737*ppar_l[6]*b_l[7]+0.3535533905932737*b_l[3]*ppar_l[5]+0.3535533905932737*ppar_l[3]*b_l[5]+0.3535533905932737*b_l[2]*ppar_l[4]+0.3535533905932737*ppar_l[2]*b_l[4]+0.3535533905932737*b_l[0]*ppar_l[1]+0.3535533905932737*ppar_l[0]*b_l[1]; 
  ppar_b_l[2] = 0.3535533905932737*b_l[5]*ppar_l[7]+0.3535533905932737*ppar_l[5]*b_l[7]+0.3535533905932737*b_l[3]*ppar_l[6]+0.3535533905932737*ppar_l[3]*b_l[6]+0.3535533905932737*b_l[1]*ppar_l[4]+0.3535533905932737*ppar_l[1]*b_l[4]+0.3535533905932737*b_l[0]*ppar_l[2]+0.3535533905932737*ppar_l[0]*b_l[2]; 
  ppar_b_l[3] = 0.3535533905932737*b_l[4]*ppar_l[7]+0.3535533905932737*ppar_l[4]*b_l[7]+0.3535533905932737*b_l[2]*ppar_l[6]+0.3535533905932737*ppar_l[2]*b_l[6]+0.3535533905932737*b_l[1]*ppar_l[5]+0.3535533905932737*ppar_l[1]*b_l[5]+0.3535533905932737*b_l[0]*ppar_l[3]+0.3535533905932737*ppar_l[0]*b_l[3]; 
  ppar_b_l[4] = 0.3535533905932737*b_l[3]*ppar_l[7]+0.3535533905932737*ppar_l[3]*b_l[7]+0.3535533905932737*b_l[5]*ppar_l[6]+0.3535533905932737*ppar_l[5]*b_l[6]+0.3535533905932737*b_l[0]*ppar_l[4]+0.3535533905932737*ppar_l[0]*b_l[4]+0.3535533905932737*b_l[1]*ppar_l[2]+0.3535533905932737*ppar_l[1]*b_l[2]; 
  ppar_b_l[5] = 0.3535533905932737*b_l[2]*ppar_l[7]+0.3535533905932737*ppar_l[2]*b_l[7]+0.3535533905932737*b_l[4]*ppar_l[6]+0.3535533905932737*ppar_l[4]*b_l[6]+0.3535533905932737*b_l[0]*ppar_l[5]+0.3535533905932737*ppar_l[0]*b_l[5]+0.3535533905932737*b_l[1]*ppar_l[3]+0.3535533905932737*ppar_l[1]*b_l[3]; 
  ppar_b_l[6] = 0.3535533905932737*b_l[1]*ppar_l[7]+0.3535533905932737*ppar_l[1]*b_l[7]+0.3535533905932737*b_l[0]*ppar_l[6]+0.3535533905932737*ppar_l[0]*b_l[6]+0.3535533905932737*b_l[4]*ppar_l[5]+0.3535533905932737*ppar_l[4]*b_l[5]+0.3535533905932737*b_l[2]*ppar_l[3]+0.3535533905932737*ppar_l[2]*b_l[3]; 
  ppar_b_l[7] = 0.3535533905932737*b_l[0]*ppar_l[7]+0.3535533905932737*ppar_l[0]*b_l[7]+0.3535533905932737*b_l[1]*ppar_l[6]+0.3535533905932737*ppar_l[1]*b_l[6]+0.3535533905932737*b_l[2]*ppar_l[5]+0.3535533905932737*ppar_l[2]*b_l[5]+0.3535533905932737*b_l[3]*ppar_l[4]+0.3535533905932737*ppar_l[3]*b_l[4]; 

  ppar_b_c[0] = 0.3535533905932737*b_c[7]*ppar_c[7]+0.3535533905932737*b_c[6]*ppar_c[6]+0.3535533905932737*b_c[5]*ppar_c[5]+0.3535533905932737*b_c[4]*ppar_c[4]+0.3535533905932737*b_c[3]*ppar_c[3]+0.3535533905932737*b_c[2]*ppar_c[2]+0.3535533905932737*b_c[1]*ppar_c[1]+0.3535533905932737*b_c[0]*ppar_c[0]; 
  ppar_b_c[1] = 0.3535533905932737*b_c[6]*ppar_c[7]+0.3535533905932737*ppar_c[6]*b_c[7]+0.3535533905932737*b_c[3]*ppar_c[5]+0.3535533905932737*ppar_c[3]*b_c[5]+0.3535533905932737*b_c[2]*ppar_c[4]+0.3535533905932737*ppar_c[2]*b_c[4]+0.3535533905932737*b_c[0]*ppar_c[1]+0.3535533905932737*ppar_c[0]*b_c[1]; 
  ppar_b_c[2] = 0.3535533905932737*b_c[5]*ppar_c[7]+0.3535533905932737*ppar_c[5]*b_c[7]+0.3535533905932737*b_c[3]*ppar_c[6]+0.3535533905932737*ppar_c[3]*b_c[6]+0.3535533905932737*b_c[1]*ppar_c[4]+0.3535533905932737*ppar_c[1]*b_c[4]+0.3535533905932737*b_c[0]*ppar_c[2]+0.3535533905932737*ppar_c[0]*b_c[2]; 
  ppar_b_c[3] = 0.3535533905932737*b_c[4]*ppar_c[7]+0.3535533905932737*ppar_c[4]*b_c[7]+0.3535533905932737*b_c[2]*ppar_c[6]+0.3535533905932737*ppar_c[2]*b_c[6]+0.3535533905932737*b_c[1]*ppar_c[5]+0.3535533905932737*ppar_c[1]*b_c[5]+0.3535533905932737*b_c[0]*ppar_c[3]+0.3535533905932737*ppar_c[0]*b_c[3]; 
  ppar_b_c[4] = 0.3535533905932737*b_c[3]*ppar_c[7]+0.3535533905932737*ppar_c[3]*b_c[7]+0.3535533905932737*b_c[5]*ppar_c[6]+0.3535533905932737*ppar_c[5]*b_c[6]+0.3535533905932737*b_c[0]*ppar_c[4]+0.3535533905932737*ppar_c[0]*b_c[4]+0.3535533905932737*b_c[1]*ppar_c[2]+0.3535533905932737*ppar_c[1]*b_c[2]; 
  ppar_b_c[5] = 0.3535533905932737*b_c[2]*ppar_c[7]+0.3535533905932737*ppar_c[2]*b_c[7]+0.3535533905932737*b_c[4]*ppar_c[6]+0.3535533905932737*ppar_c[4]*b_c[6]+0.3535533905932737*b_c[0]*ppar_c[5]+0.3535533905932737*ppar_c[0]*b_c[5]+0.3535533905932737*b_c[1]*ppar_c[3]+0.3535533905932737*ppar_c[1]*b_c[3]; 
  ppar_b_c[6] = 0.3535533905932737*b_c[1]*ppar_c[7]+0.3535533905932737*ppar_c[1]*b_c[7]+0.3535533905932737*b_c[0]*ppar_c[6]+0.3535533905932737*ppar_c[0]*b_c[6]+0.3535533905932737*b_c[4]*ppar_c[5]+0.3535533905932737*ppar_c[4]*b_c[5]+0.3535533905932737*b_c[2]*ppar_c[3]+0.3535533905932737*ppar_c[2]*b_c[3]; 
  ppar_b_c[7] = 0.3535533905932737*b_c[0]*ppar_c[7]+0.3535533905932737*ppar_c[0]*b_c[7]+0.3535533905932737*b_c[1]*ppar_c[6]+0.3535533905932737*ppar_c[1]*b_c[6]+0.3535533905932737*b_c[2]*ppar_c[5]+0.3535533905932737*ppar_c[2]*b_c[5]+0.3535533905932737*b_c[3]*ppar_c[4]+0.3535533905932737*ppar_c[3]*b_c[4]; 

  ppar_b_r[0] = 0.3535533905932737*b_r[7]*ppar_r[7]+0.3535533905932737*b_r[6]*ppar_r[6]+0.3535533905932737*b_r[5]*ppar_r[5]+0.3535533905932737*b_r[4]*ppar_r[4]+0.3535533905932737*b_r[3]*ppar_r[3]+0.3535533905932737*b_r[2]*ppar_r[2]+0.3535533905932737*b_r[1]*ppar_r[1]+0.3535533905932737*b_r[0]*ppar_r[0]; 
  ppar_b_r[1] = 0.3535533905932737*b_r[6]*ppar_r[7]+0.3535533905932737*ppar_r[6]*b_r[7]+0.3535533905932737*b_r[3]*ppar_r[5]+0.3535533905932737*ppar_r[3]*b_r[5]+0.3535533905932737*b_r[2]*ppar_r[4]+0.3535533905932737*ppar_r[2]*b_r[4]+0.3535533905932737*b_r[0]*ppar_r[1]+0.3535533905932737*ppar_r[0]*b_r[1]; 
  ppar_b_r[2] = 0.3535533905932737*b_r[5]*ppar_r[7]+0.3535533905932737*ppar_r[5]*b_r[7]+0.3535533905932737*b_r[3]*ppar_r[6]+0.3535533905932737*ppar_r[3]*b_r[6]+0.3535533905932737*b_r[1]*ppar_r[4]+0.3535533905932737*ppar_r[1]*b_r[4]+0.3535533905932737*b_r[0]*ppar_r[2]+0.3535533905932737*ppar_r[0]*b_r[2]; 
  ppar_b_r[3] = 0.3535533905932737*b_r[4]*ppar_r[7]+0.3535533905932737*ppar_r[4]*b_r[7]+0.3535533905932737*b_r[2]*ppar_r[6]+0.3535533905932737*ppar_r[2]*b_r[6]+0.3535533905932737*b_r[1]*ppar_r[5]+0.3535533905932737*ppar_r[1]*b_r[5]+0.3535533905932737*b_r[0]*ppar_r[3]+0.3535533905932737*ppar_r[0]*b_r[3]; 
  ppar_b_r[4] = 0.3535533905932737*b_r[3]*ppar_r[7]+0.3535533905932737*ppar_r[3]*b_r[7]+0.3535533905932737*b_r[5]*ppar_r[6]+0.3535533905932737*ppar_r[5]*b_r[6]+0.3535533905932737*b_r[0]*ppar_r[4]+0.3535533905932737*ppar_r[0]*b_r[4]+0.3535533905932737*b_r[1]*ppar_r[2]+0.3535533905932737*ppar_r[1]*b_r[2]; 
  ppar_b_r[5] = 0.3535533905932737*b_r[2]*ppar_r[7]+0.3535533905932737*ppar_r[2]*b_r[7]+0.3535533905932737*b_r[4]*ppar_r[6]+0.3535533905932737*ppar_r[4]*b_r[6]+0.3535533905932737*b_r[0]*ppar_r[5]+0.3535533905932737*ppar_r[0]*b_r[5]+0.3535533905932737*b_r[1]*ppar_r[3]+0.3535533905932737*ppar_r[1]*b_r[3]; 
  ppar_b_r[6] = 0.3535533905932737*b_r[1]*ppar_r[7]+0.3535533905932737*ppar_r[1]*b_r[7]+0.3535533905932737*b_r[0]*ppar_r[6]+0.3535533905932737*ppar_r[0]*b_r[6]+0.3535533905932737*b_r[4]*ppar_r[5]+0.3535533905932737*ppar_r[4]*b_r[5]+0.3535533905932737*b_r[2]*ppar_r[3]+0.3535533905932737*ppar_r[2]*b_r[3]; 
  ppar_b_r[7] = 0.3535533905932737*b_r[0]*ppar_r[7]+0.3535533905932737*ppar_r[0]*b_r[7]+0.3535533905932737*b_r[1]*ppar_r[6]+0.3535533905932737*ppar_r[1]*b_r[6]+0.3535533905932737*b_r[2]*ppar_r[5]+0.3535533905932737*ppar_r[2]*b_r[5]+0.3535533905932737*b_r[3]*ppar_r[4]+0.3535533905932737*ppar_r[3]*b_r[4]; 

  div_b[0] += ((-0.2886751345948129*(b_r[3]+b_l[3]))+0.5773502691896258*b_c[3]+0.25*b_r[0]-0.25*b_l[0])*dx1; 
  div_b[1] += ((-0.2886751345948129*(b_r[5]+b_l[5]))+0.5773502691896258*b_c[5]+0.25*b_r[1]-0.25*b_l[1])*dx1; 
  div_b[2] += ((-0.2886751345948129*(b_r[6]+b_l[6]))+0.5773502691896258*b_c[6]+0.25*b_r[2]-0.25*b_l[2])*dx1; 
  div_b[3] += ((-0.5*b_r[3])+0.5*b_l[3]+0.4330127018922193*(b_r[0]+b_l[0])-0.8660254037844386*b_c[0])*dx1; 
  div_b[4] += ((-0.2886751345948129*(b_r[7]+b_l[7]))+0.5773502691896258*b_c[7]+0.25*b_r[4]-0.25*b_l[4])*dx1; 
  div_b[5] += ((-0.5*b_r[5])+0.5*b_l[5]+0.4330127018922193*(b_r[1]+b_l[1])-0.8660254037844386*b_c[1])*dx1; 
  div_b[6] += ((-0.5*b_r[6])+0.5*b_l[6]+0.4330127018922193*(b_r[2]+b_l[2])-0.8660254037844386*b_c[2])*dx1; 
  div_b[7] += ((-0.5*b_r[7])+0.5*b_l[7]+0.4330127018922193*(b_r[4]+b_l[4])-0.8660254037844386*b_c[4])*dx1; 

  bb_grad_u[0] += 0.3535533905932737*(bzbz[7]*grad_u_z[7]+bybz[7]*grad_u_y[7]+bxbz[7]*grad_u_x[7]+bzbz[6]*grad_u_z[6]+bybz[6]*grad_u_y[6]+bxbz[6]*grad_u_x[6]+bzbz[5]*grad_u_z[5]+bybz[5]*grad_u_y[5]+bxbz[5]*grad_u_x[5]+bzbz[4]*grad_u_z[4]+bybz[4]*grad_u_y[4]+bxbz[4]*grad_u_x[4]+bzbz[3]*grad_u_z[3]+bybz[3]*grad_u_y[3]+bxbz[3]*grad_u_x[3]+bzbz[2]*grad_u_z[2]+bybz[2]*grad_u_y[2]+bxbz[2]*grad_u_x[2]+bzbz[1]*grad_u_z[1]+bybz[1]*grad_u_y[1]+bxbz[1]*grad_u_x[1]+bzbz[0]*grad_u_z[0]+bybz[0]*grad_u_y[0]+bxbz[0]*grad_u_x[0])*dx1; 
  bb_grad_u[1] += 0.3535533905932737*(bzbz[6]*grad_u_z[7]+bybz[6]*grad_u_y[7]+bxbz[6]*grad_u_x[7]+grad_u_z[6]*bzbz[7]+grad_u_y[6]*bybz[7]+grad_u_x[6]*bxbz[7]+bzbz[3]*grad_u_z[5]+bybz[3]*grad_u_y[5]+bxbz[3]*grad_u_x[5]+grad_u_z[3]*bzbz[5]+grad_u_y[3]*bybz[5]+grad_u_x[3]*bxbz[5]+bzbz[2]*grad_u_z[4]+bybz[2]*grad_u_y[4]+bxbz[2]*grad_u_x[4]+grad_u_z[2]*bzbz[4]+grad_u_y[2]*bybz[4]+grad_u_x[2]*bxbz[4]+bzbz[0]*grad_u_z[1]+bybz[0]*grad_u_y[1]+bxbz[0]*grad_u_x[1]+grad_u_z[0]*bzbz[1]+grad_u_y[0]*bybz[1]+grad_u_x[0]*bxbz[1])*dx1; 
  bb_grad_u[2] += 0.3535533905932737*(bzbz[5]*grad_u_z[7]+bybz[5]*grad_u_y[7]+bxbz[5]*grad_u_x[7]+grad_u_z[5]*bzbz[7]+grad_u_y[5]*bybz[7]+grad_u_x[5]*bxbz[7]+bzbz[3]*grad_u_z[6]+bybz[3]*grad_u_y[6]+bxbz[3]*grad_u_x[6]+grad_u_z[3]*bzbz[6]+grad_u_y[3]*bybz[6]+grad_u_x[3]*bxbz[6]+bzbz[1]*grad_u_z[4]+bybz[1]*grad_u_y[4]+bxbz[1]*grad_u_x[4]+grad_u_z[1]*bzbz[4]+grad_u_y[1]*bybz[4]+grad_u_x[1]*bxbz[4]+bzbz[0]*grad_u_z[2]+bybz[0]*grad_u_y[2]+bxbz[0]*grad_u_x[2]+grad_u_z[0]*bzbz[2]+grad_u_y[0]*bybz[2]+grad_u_x[0]*bxbz[2])*dx1; 
  bb_grad_u[3] += 0.3535533905932737*(bzbz[4]*grad_u_z[7]+bybz[4]*grad_u_y[7]+bxbz[4]*grad_u_x[7]+grad_u_z[4]*bzbz[7]+grad_u_y[4]*bybz[7]+grad_u_x[4]*bxbz[7]+bzbz[2]*grad_u_z[6]+bybz[2]*grad_u_y[6]+bxbz[2]*grad_u_x[6]+grad_u_z[2]*bzbz[6]+grad_u_y[2]*bybz[6]+grad_u_x[2]*bxbz[6]+bzbz[1]*grad_u_z[5]+bybz[1]*grad_u_y[5]+bxbz[1]*grad_u_x[5]+grad_u_z[1]*bzbz[5]+grad_u_y[1]*bybz[5]+grad_u_x[1]*bxbz[5]+bzbz[0]*grad_u_z[3]+bybz[0]*grad_u_y[3]+bxbz[0]*grad_u_x[3]+grad_u_z[0]*bzbz[3]+grad_u_y[0]*bybz[3]+grad_u_x[0]*bxbz[3])*dx1; 
  bb_grad_u[4] += 0.3535533905932737*(bzbz[3]*grad_u_z[7]+bybz[3]*grad_u_y[7]+bxbz[3]*grad_u_x[7]+grad_u_z[3]*bzbz[7]+grad_u_y[3]*bybz[7]+grad_u_x[3]*bxbz[7]+bzbz[5]*grad_u_z[6]+bybz[5]*grad_u_y[6]+bxbz[5]*grad_u_x[6]+grad_u_z[5]*bzbz[6]+grad_u_y[5]*bybz[6]+grad_u_x[5]*bxbz[6]+bzbz[0]*grad_u_z[4]+bybz[0]*grad_u_y[4]+bxbz[0]*grad_u_x[4]+grad_u_z[0]*bzbz[4]+grad_u_y[0]*bybz[4]+grad_u_x[0]*bxbz[4]+bzbz[1]*grad_u_z[2]+bybz[1]*grad_u_y[2]+bxbz[1]*grad_u_x[2]+grad_u_z[1]*bzbz[2]+grad_u_y[1]*bybz[2]+grad_u_x[1]*bxbz[2])*dx1; 
  bb_grad_u[5] += 0.3535533905932737*(bzbz[2]*grad_u_z[7]+bybz[2]*grad_u_y[7]+bxbz[2]*grad_u_x[7]+grad_u_z[2]*bzbz[7]+grad_u_y[2]*bybz[7]+grad_u_x[2]*bxbz[7]+bzbz[4]*grad_u_z[6]+bybz[4]*grad_u_y[6]+bxbz[4]*grad_u_x[6]+grad_u_z[4]*bzbz[6]+grad_u_y[4]*bybz[6]+grad_u_x[4]*bxbz[6]+bzbz[0]*grad_u_z[5]+bybz[0]*grad_u_y[5]+bxbz[0]*grad_u_x[5]+grad_u_z[0]*bzbz[5]+grad_u_y[0]*bybz[5]+grad_u_x[0]*bxbz[5]+bzbz[1]*grad_u_z[3]+bybz[1]*grad_u_y[3]+bxbz[1]*grad_u_x[3]+grad_u_z[1]*bzbz[3]+grad_u_y[1]*bybz[3]+grad_u_x[1]*bxbz[3])*dx1; 
  bb_grad_u[6] += 0.3535533905932737*(bzbz[1]*grad_u_z[7]+bybz[1]*grad_u_y[7]+bxbz[1]*grad_u_x[7]+grad_u_z[1]*bzbz[7]+grad_u_y[1]*bybz[7]+grad_u_x[1]*bxbz[7]+bzbz[0]*grad_u_z[6]+bybz[0]*grad_u_y[6]+bxbz[0]*grad_u_x[6]+grad_u_z[0]*bzbz[6]+grad_u_y[0]*bybz[6]+grad_u_x[0]*bxbz[6]+bzbz[4]*grad_u_z[5]+bybz[4]*grad_u_y[5]+bxbz[4]*grad_u_x[5]+grad_u_z[4]*bzbz[5]+grad_u_y[4]*bybz[5]+grad_u_x[4]*bxbz[5]+bzbz[2]*grad_u_z[3]+bybz[2]*grad_u_y[3]+bxbz[2]*grad_u_x[3]+grad_u_z[2]*bzbz[3]+grad_u_y[2]*bybz[3]+grad_u_x[2]*bxbz[3])*dx1; 
  bb_grad_u[7] += 0.3535533905932737*(bzbz[0]*grad_u_z[7]+bybz[0]*grad_u_y[7]+bxbz[0]*grad_u_x[7]+grad_u_z[0]*bzbz[7]+grad_u_y[0]*bybz[7]+grad_u_x[0]*bxbz[7]+bzbz[1]*grad_u_z[6]+bybz[1]*grad_u_y[6]+bxbz[1]*grad_u_x[6]+grad_u_z[1]*bzbz[6]+grad_u_y[1]*bybz[6]+grad_u_x[1]*bxbz[6]+bzbz[2]*grad_u_z[5]+bybz[2]*grad_u_y[5]+bxbz[2]*grad_u_x[5]+grad_u_z[2]*bzbz[5]+grad_u_y[2]*bybz[5]+grad_u_x[2]*bxbz[5]+bzbz[3]*grad_u_z[4]+bybz[3]*grad_u_y[4]+bxbz[3]*grad_u_x[4]+grad_u_z[3]*bzbz[4]+grad_u_y[3]*bybz[4]+grad_u_x[3]*bxbz[4])*dx1; 

  div_p_x[0] += (4.871392896287466*rhoux_r[3]-4.871392896287466*rhoux_l[3]-2.8125*(rhoux_r[0]+rhoux_l[0])+5.625*rhoux_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[3]+Pxz_l[3]))+0.5773502691896258*Pxz_c[3]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_x[1] += (4.871392896287466*rhoux_r[5]-4.871392896287466*rhoux_l[5]-2.8125*(rhoux_r[1]+rhoux_l[1])+5.625*rhoux_c[1])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[5]+Pxz_l[5]))+0.5773502691896258*Pxz_c[5]+0.25*Pxz_r[1]-0.25*Pxz_l[1])*dx1; 
  div_p_x[2] += (4.871392896287466*rhoux_r[6]-4.871392896287466*rhoux_l[6]-2.8125*(rhoux_r[2]+rhoux_l[2])+5.625*rhoux_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[6]+Pxz_l[6]))+0.5773502691896258*Pxz_c[6]+0.25*Pxz_r[2]-0.25*Pxz_l[2])*dx1; 
  div_p_x[3] += (72.1875*(rhoux_r[3]+rhoux_l[3])+249.375*rhoux_c[3]-56.83291712335378*rhoux_r[0]+56.83291712335378*rhoux_l[0])*dx14*nuHyp+((-0.5*Pxz_r[3])+0.5*Pxz_l[3]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 
  div_p_x[4] += (4.871392896287466*rhoux_r[7]-4.871392896287466*rhoux_l[7]-2.8125*(rhoux_r[4]+rhoux_l[4])+5.625*rhoux_c[4])*dx14*nuHyp+((-0.2886751345948129*(Pxz_r[7]+Pxz_l[7]))+0.5773502691896258*Pxz_c[7]+0.25*Pxz_r[4]-0.25*Pxz_l[4])*dx1; 
  div_p_x[5] += (72.1875*(rhoux_r[5]+rhoux_l[5])+249.375*rhoux_c[5]-56.83291712335378*rhoux_r[1]+56.83291712335378*rhoux_l[1])*dx14*nuHyp+((-0.5*Pxz_r[5])+0.5*Pxz_l[5]+0.4330127018922193*(Pxz_r[1]+Pxz_l[1])-0.8660254037844386*Pxz_c[1])*dx1; 
  div_p_x[6] += (72.1875*(rhoux_r[6]+rhoux_l[6])+249.375*rhoux_c[6]-56.83291712335378*rhoux_r[2]+56.83291712335378*rhoux_l[2])*dx14*nuHyp+((-0.5*Pxz_r[6])+0.5*Pxz_l[6]+0.4330127018922193*(Pxz_r[2]+Pxz_l[2])-0.8660254037844386*Pxz_c[2])*dx1; 
  div_p_x[7] += (72.1875*(rhoux_r[7]+rhoux_l[7])+249.375*rhoux_c[7]-56.83291712335378*rhoux_r[4]+56.83291712335378*rhoux_l[4])*dx14*nuHyp+((-0.5*Pxz_r[7])+0.5*Pxz_l[7]+0.4330127018922193*(Pxz_r[4]+Pxz_l[4])-0.8660254037844386*Pxz_c[4])*dx1; 

  div_p_y[0] += (4.871392896287466*rhouy_r[3]-4.871392896287466*rhouy_l[3]-2.8125*(rhouy_r[0]+rhouy_l[0])+5.625*rhouy_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pyz_r[3]+Pyz_l[3]))+0.5773502691896258*Pyz_c[3]+0.25*Pyz_r[0]-0.25*Pyz_l[0])*dx1; 
  div_p_y[1] += (4.871392896287466*rhouy_r[5]-4.871392896287466*rhouy_l[5]-2.8125*(rhouy_r[1]+rhouy_l[1])+5.625*rhouy_c[1])*dx14*nuHyp+((-0.2886751345948129*(Pyz_r[5]+Pyz_l[5]))+0.5773502691896258*Pyz_c[5]+0.25*Pyz_r[1]-0.25*Pyz_l[1])*dx1; 
  div_p_y[2] += (4.871392896287466*rhouy_r[6]-4.871392896287466*rhouy_l[6]-2.8125*(rhouy_r[2]+rhouy_l[2])+5.625*rhouy_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pyz_r[6]+Pyz_l[6]))+0.5773502691896258*Pyz_c[6]+0.25*Pyz_r[2]-0.25*Pyz_l[2])*dx1; 
  div_p_y[3] += (72.1875*(rhouy_r[3]+rhouy_l[3])+249.375*rhouy_c[3]-56.83291712335378*rhouy_r[0]+56.83291712335378*rhouy_l[0])*dx14*nuHyp+((-0.5*Pyz_r[3])+0.5*Pyz_l[3]+0.4330127018922193*(Pyz_r[0]+Pyz_l[0])-0.8660254037844386*Pyz_c[0])*dx1; 
  div_p_y[4] += (4.871392896287466*rhouy_r[7]-4.871392896287466*rhouy_l[7]-2.8125*(rhouy_r[4]+rhouy_l[4])+5.625*rhouy_c[4])*dx14*nuHyp+((-0.2886751345948129*(Pyz_r[7]+Pyz_l[7]))+0.5773502691896258*Pyz_c[7]+0.25*Pyz_r[4]-0.25*Pyz_l[4])*dx1; 
  div_p_y[5] += (72.1875*(rhouy_r[5]+rhouy_l[5])+249.375*rhouy_c[5]-56.83291712335378*rhouy_r[1]+56.83291712335378*rhouy_l[1])*dx14*nuHyp+((-0.5*Pyz_r[5])+0.5*Pyz_l[5]+0.4330127018922193*(Pyz_r[1]+Pyz_l[1])-0.8660254037844386*Pyz_c[1])*dx1; 
  div_p_y[6] += (72.1875*(rhouy_r[6]+rhouy_l[6])+249.375*rhouy_c[6]-56.83291712335378*rhouy_r[2]+56.83291712335378*rhouy_l[2])*dx14*nuHyp+((-0.5*Pyz_r[6])+0.5*Pyz_l[6]+0.4330127018922193*(Pyz_r[2]+Pyz_l[2])-0.8660254037844386*Pyz_c[2])*dx1; 
  div_p_y[7] += (72.1875*(rhouy_r[7]+rhouy_l[7])+249.375*rhouy_c[7]-56.83291712335378*rhouy_r[4]+56.83291712335378*rhouy_l[4])*dx14*nuHyp+((-0.5*Pyz_r[7])+0.5*Pyz_l[7]+0.4330127018922193*(Pyz_r[4]+Pyz_l[4])-0.8660254037844386*Pyz_c[4])*dx1; 

  div_p_z[0] += (4.871392896287466*rhouy_r[3]-4.871392896287466*rhouy_l[3]-2.8125*(rhouy_r[0]+rhouy_l[0])+5.625*rhouy_c[0])*dx14*nuHyp+((-0.2886751345948129*(Pzz_r[3]+Pzz_l[3]))+0.5773502691896258*Pzz_c[3]+0.25*Pzz_r[0]-0.25*Pzz_l[0])*dx1; 
  div_p_z[1] += (4.871392896287466*rhouy_r[5]-4.871392896287466*rhouy_l[5]-2.8125*(rhouy_r[1]+rhouy_l[1])+5.625*rhouy_c[1])*dx14*nuHyp+((-0.2886751345948129*(Pzz_r[5]+Pzz_l[5]))+0.5773502691896258*Pzz_c[5]+0.25*Pzz_r[1]-0.25*Pzz_l[1])*dx1; 
  div_p_z[2] += (4.871392896287466*rhouy_r[6]-4.871392896287466*rhouy_l[6]-2.8125*(rhouy_r[2]+rhouy_l[2])+5.625*rhouy_c[2])*dx14*nuHyp+((-0.2886751345948129*(Pzz_r[6]+Pzz_l[6]))+0.5773502691896258*Pzz_c[6]+0.25*Pzz_r[2]-0.25*Pzz_l[2])*dx1; 
  div_p_z[3] += (72.1875*(rhouy_r[3]+rhouy_l[3])+249.375*rhouy_c[3]-56.83291712335378*rhouy_r[0]+56.83291712335378*rhouy_l[0])*dx14*nuHyp+((-0.5*Pzz_r[3])+0.5*Pzz_l[3]+0.4330127018922193*(Pzz_r[0]+Pzz_l[0])-0.8660254037844386*Pzz_c[0])*dx1; 
  div_p_z[4] += (4.871392896287466*rhouy_r[7]-4.871392896287466*rhouy_l[7]-2.8125*(rhouy_r[4]+rhouy_l[4])+5.625*rhouy_c[4])*dx14*nuHyp+((-0.2886751345948129*(Pzz_r[7]+Pzz_l[7]))+0.5773502691896258*Pzz_c[7]+0.25*Pzz_r[4]-0.25*Pzz_l[4])*dx1; 
  div_p_z[5] += (72.1875*(rhouy_r[5]+rhouy_l[5])+249.375*rhouy_c[5]-56.83291712335378*rhouy_r[1]+56.83291712335378*rhouy_l[1])*dx14*nuHyp+((-0.5*Pzz_r[5])+0.5*Pzz_l[5]+0.4330127018922193*(Pzz_r[1]+Pzz_l[1])-0.8660254037844386*Pzz_c[1])*dx1; 
  div_p_z[6] += (72.1875*(rhouy_r[6]+rhouy_l[6])+249.375*rhouy_c[6]-56.83291712335378*rhouy_r[2]+56.83291712335378*rhouy_l[2])*dx14*nuHyp+((-0.5*Pzz_r[6])+0.5*Pzz_l[6]+0.4330127018922193*(Pzz_r[2]+Pzz_l[2])-0.8660254037844386*Pzz_c[2])*dx1; 
  div_p_z[7] += (72.1875*(rhouy_r[7]+rhouy_l[7])+249.375*rhouy_c[7]-56.83291712335378*rhouy_r[4]+56.83291712335378*rhouy_l[4])*dx14*nuHyp+((-0.5*Pzz_r[7])+0.5*Pzz_l[7]+0.4330127018922193*(Pzz_r[4]+Pzz_l[4])-0.8660254037844386*Pzz_c[4])*dx1; 

  double div_ppar_b[8] = {0.0}; 
  double rho_inv[8] = {0.0}; 
  ser_3x_p1_inv(rho, rho_inv); 
  div_ppar_b[0] = (-0.2886751345948129*ppar_b_r[3])-0.2886751345948129*ppar_b_l[3]+0.5773502691896258*ppar_b_c[3]+0.25*ppar_b_r[0]-0.25*ppar_b_l[0]; 
  div_ppar_b[1] = (-0.2886751345948129*ppar_b_r[5])-0.2886751345948129*ppar_b_l[5]+0.5773502691896258*ppar_b_c[5]+0.25*ppar_b_r[1]-0.25*ppar_b_l[1]; 
  div_ppar_b[2] = (-0.2886751345948129*ppar_b_r[6])-0.2886751345948129*ppar_b_l[6]+0.5773502691896258*ppar_b_c[6]+0.25*ppar_b_r[2]-0.25*ppar_b_l[2]; 
  div_ppar_b[3] = (-0.5*ppar_b_r[3])+0.5*ppar_b_l[3]+0.4330127018922193*ppar_b_r[0]+0.4330127018922193*ppar_b_l[0]-0.8660254037844386*ppar_b_c[0]; 
  div_ppar_b[4] = (-0.2886751345948129*ppar_b_r[7])-0.2886751345948129*ppar_b_l[7]+0.5773502691896258*ppar_b_c[7]+0.25*ppar_b_r[4]-0.25*ppar_b_l[4]; 
  div_ppar_b[5] = (-0.5*ppar_b_r[5])+0.5*ppar_b_l[5]+0.4330127018922193*ppar_b_r[1]+0.4330127018922193*ppar_b_l[1]-0.8660254037844386*ppar_b_c[1]; 
  div_ppar_b[6] = (-0.5*ppar_b_r[6])+0.5*ppar_b_l[6]+0.4330127018922193*ppar_b_r[2]+0.4330127018922193*ppar_b_l[2]-0.8660254037844386*ppar_b_c[2]; 
  div_ppar_b[7] = (-0.5*ppar_b_r[7])+0.5*ppar_b_l[7]+0.4330127018922193*ppar_b_r[4]+0.4330127018922193*ppar_b_l[4]-0.8660254037844386*ppar_b_c[4]; 

  p_force[0] += 0.3535533905932737*(div_ppar_b[7]*rho_inv[7]+div_ppar_b[6]*rho_inv[6]+div_ppar_b[5]*rho_inv[5]+div_ppar_b[4]*rho_inv[4]+div_ppar_b[3]*rho_inv[3]+div_ppar_b[2]*rho_inv[2]+div_ppar_b[1]*rho_inv[1]+div_ppar_b[0]*rho_inv[0])*dx1; 
  p_force[1] += 0.3535533905932737*(div_ppar_b[6]*rho_inv[7]+rho_inv[6]*div_ppar_b[7]+div_ppar_b[3]*rho_inv[5]+rho_inv[3]*div_ppar_b[5]+div_ppar_b[2]*rho_inv[4]+rho_inv[2]*div_ppar_b[4]+div_ppar_b[0]*rho_inv[1]+rho_inv[0]*div_ppar_b[1])*dx1; 
  p_force[2] += 0.3535533905932737*(div_ppar_b[5]*rho_inv[7]+rho_inv[5]*div_ppar_b[7]+div_ppar_b[3]*rho_inv[6]+rho_inv[3]*div_ppar_b[6]+div_ppar_b[1]*rho_inv[4]+rho_inv[1]*div_ppar_b[4]+div_ppar_b[0]*rho_inv[2]+rho_inv[0]*div_ppar_b[2])*dx1; 
  p_force[3] += 0.3535533905932737*(div_ppar_b[4]*rho_inv[7]+rho_inv[4]*div_ppar_b[7]+div_ppar_b[2]*rho_inv[6]+rho_inv[2]*div_ppar_b[6]+div_ppar_b[1]*rho_inv[5]+rho_inv[1]*div_ppar_b[5]+div_ppar_b[0]*rho_inv[3]+rho_inv[0]*div_ppar_b[3])*dx1; 
  p_force[4] += 0.3535533905932737*(div_ppar_b[3]*rho_inv[7]+rho_inv[3]*div_ppar_b[7]+div_ppar_b[5]*rho_inv[6]+rho_inv[5]*div_ppar_b[6]+div_ppar_b[0]*rho_inv[4]+rho_inv[0]*div_ppar_b[4]+div_ppar_b[1]*rho_inv[2]+rho_inv[1]*div_ppar_b[2])*dx1; 
  p_force[5] += 0.3535533905932737*(div_ppar_b[2]*rho_inv[7]+rho_inv[2]*div_ppar_b[7]+div_ppar_b[4]*rho_inv[6]+rho_inv[4]*div_ppar_b[6]+div_ppar_b[0]*rho_inv[5]+rho_inv[0]*div_ppar_b[5]+div_ppar_b[1]*rho_inv[3]+rho_inv[1]*div_ppar_b[3])*dx1; 
  p_force[6] += 0.3535533905932737*(div_ppar_b[1]*rho_inv[7]+rho_inv[1]*div_ppar_b[7]+div_ppar_b[0]*rho_inv[6]+rho_inv[0]*div_ppar_b[6]+div_ppar_b[4]*rho_inv[5]+rho_inv[4]*div_ppar_b[5]+div_ppar_b[2]*rho_inv[3]+rho_inv[2]*div_ppar_b[3])*dx1; 
  p_force[7] += 0.3535533905932737*(div_ppar_b[0]*rho_inv[7]+rho_inv[0]*div_ppar_b[7]+div_ppar_b[1]*rho_inv[6]+rho_inv[1]*div_ppar_b[6]+div_ppar_b[2]*rho_inv[5]+rho_inv[2]*div_ppar_b[5]+div_ppar_b[3]*rho_inv[4]+rho_inv[3]*div_ppar_b[4])*dx1; 

  p_perp_source[0] += 0.3535533905932737*(bb_grad_u[7]*p_perp[7]+bb_grad_u[6]*p_perp[6]+bb_grad_u[5]*p_perp[5]+bb_grad_u[4]*p_perp[4]+bb_grad_u[3]*p_perp[3]+bb_grad_u[2]*p_perp[2]+bb_grad_u[1]*p_perp[1]+bb_grad_u[0]*p_perp[0])-0.3535533905932737*(grad_u_z[7]*p_perp[7]+grad_u_z[6]*p_perp[6]+grad_u_z[5]*p_perp[5]+grad_u_z[4]*p_perp[4]+grad_u_z[3]*p_perp[3]+grad_u_z[2]*p_perp[2]+grad_u_z[1]*p_perp[1]+grad_u_z[0]*p_perp[0])*dx1; 
  p_perp_source[1] += 0.3535533905932737*(bb_grad_u[6]*p_perp[7]+p_perp[6]*bb_grad_u[7]+bb_grad_u[3]*p_perp[5]+p_perp[3]*bb_grad_u[5]+bb_grad_u[2]*p_perp[4]+p_perp[2]*bb_grad_u[4]+bb_grad_u[0]*p_perp[1]+p_perp[0]*bb_grad_u[1])-0.3535533905932737*(grad_u_z[6]*p_perp[7]+p_perp[6]*grad_u_z[7]+grad_u_z[3]*p_perp[5]+p_perp[3]*grad_u_z[5]+grad_u_z[2]*p_perp[4]+p_perp[2]*grad_u_z[4]+grad_u_z[0]*p_perp[1]+p_perp[0]*grad_u_z[1])*dx1; 
  p_perp_source[2] += 0.3535533905932737*(bb_grad_u[5]*p_perp[7]+p_perp[5]*bb_grad_u[7]+bb_grad_u[3]*p_perp[6]+p_perp[3]*bb_grad_u[6]+bb_grad_u[1]*p_perp[4]+p_perp[1]*bb_grad_u[4]+bb_grad_u[0]*p_perp[2]+p_perp[0]*bb_grad_u[2])-0.3535533905932737*(grad_u_z[5]*p_perp[7]+p_perp[5]*grad_u_z[7]+grad_u_z[3]*p_perp[6]+p_perp[3]*grad_u_z[6]+grad_u_z[1]*p_perp[4]+p_perp[1]*grad_u_z[4]+grad_u_z[0]*p_perp[2]+p_perp[0]*grad_u_z[2])*dx1; 
  p_perp_source[3] += 0.3535533905932737*(bb_grad_u[4]*p_perp[7]+p_perp[4]*bb_grad_u[7]+bb_grad_u[2]*p_perp[6]+p_perp[2]*bb_grad_u[6]+bb_grad_u[1]*p_perp[5]+p_perp[1]*bb_grad_u[5]+bb_grad_u[0]*p_perp[3]+p_perp[0]*bb_grad_u[3])-0.3535533905932737*(grad_u_z[4]*p_perp[7]+p_perp[4]*grad_u_z[7]+grad_u_z[2]*p_perp[6]+p_perp[2]*grad_u_z[6]+grad_u_z[1]*p_perp[5]+p_perp[1]*grad_u_z[5]+grad_u_z[0]*p_perp[3]+p_perp[0]*grad_u_z[3])*dx1; 
  p_perp_source[4] += 0.3535533905932737*(bb_grad_u[3]*p_perp[7]+p_perp[3]*bb_grad_u[7]+bb_grad_u[5]*p_perp[6]+p_perp[5]*bb_grad_u[6]+bb_grad_u[0]*p_perp[4]+p_perp[0]*bb_grad_u[4]+bb_grad_u[1]*p_perp[2]+p_perp[1]*bb_grad_u[2])-0.3535533905932737*(grad_u_z[3]*p_perp[7]+p_perp[3]*grad_u_z[7]+grad_u_z[5]*p_perp[6]+p_perp[5]*grad_u_z[6]+grad_u_z[0]*p_perp[4]+p_perp[0]*grad_u_z[4]+grad_u_z[1]*p_perp[2]+p_perp[1]*grad_u_z[2])*dx1; 
  p_perp_source[5] += 0.3535533905932737*(bb_grad_u[2]*p_perp[7]+p_perp[2]*bb_grad_u[7]+bb_grad_u[4]*p_perp[6]+p_perp[4]*bb_grad_u[6]+bb_grad_u[0]*p_perp[5]+p_perp[0]*bb_grad_u[5]+bb_grad_u[1]*p_perp[3]+p_perp[1]*bb_grad_u[3])-0.3535533905932737*(grad_u_z[2]*p_perp[7]+p_perp[2]*grad_u_z[7]+grad_u_z[4]*p_perp[6]+p_perp[4]*grad_u_z[6]+grad_u_z[0]*p_perp[5]+p_perp[0]*grad_u_z[5]+grad_u_z[1]*p_perp[3]+p_perp[1]*grad_u_z[3])*dx1; 
  p_perp_source[6] += 0.3535533905932737*(bb_grad_u[1]*p_perp[7]+p_perp[1]*bb_grad_u[7]+bb_grad_u[0]*p_perp[6]+p_perp[0]*bb_grad_u[6]+bb_grad_u[4]*p_perp[5]+p_perp[4]*bb_grad_u[5]+bb_grad_u[2]*p_perp[3]+p_perp[2]*bb_grad_u[3])-0.3535533905932737*(grad_u_z[1]*p_perp[7]+p_perp[1]*grad_u_z[7]+grad_u_z[0]*p_perp[6]+p_perp[0]*grad_u_z[6]+grad_u_z[4]*p_perp[5]+p_perp[4]*grad_u_z[5]+grad_u_z[2]*p_perp[3]+p_perp[2]*grad_u_z[3])*dx1; 
  p_perp_source[7] += 0.3535533905932737*(bb_grad_u[0]*p_perp[7]+p_perp[0]*bb_grad_u[7]+bb_grad_u[1]*p_perp[6]+p_perp[1]*bb_grad_u[6]+bb_grad_u[2]*p_perp[5]+p_perp[2]*bb_grad_u[5]+bb_grad_u[3]*p_perp[4]+p_perp[3]*bb_grad_u[4])-0.3535533905932737*(grad_u_z[0]*p_perp[7]+p_perp[0]*grad_u_z[7]+grad_u_z[1]*p_perp[6]+p_perp[1]*grad_u_z[6]+grad_u_z[2]*p_perp[5]+p_perp[2]*grad_u_z[5]+grad_u_z[3]*p_perp[4]+p_perp[3]*grad_u_z[4])*dx1; 

} 
