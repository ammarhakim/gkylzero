#include <gkyl_euler_kernels.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p2(const double *dxv, double nuHyp, 
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

  const double *ppar_l = &vlasov_pkpm_momsl[3]; 
  const double *ppar_c = &vlasov_pkpm_momsc[3]; 
  const double *ppar_r = &vlasov_pkpm_momsr[3]; 
  const double *rho = &vlasov_pkpm_momsc[0]; 
  const double *p_perp = &vlasov_pkpm_momsc[6]; 

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

  double *p_force_F_0 = &p_force[0]; 
  double *p_force_G_1 = &p_force[3]; 

  double grad_u_x[3] = {0.0}; 
  double grad_u_y[3] = {0.0}; 
  double grad_u_z[3] = {0.0}; 
  grad_u_x[0] = 0.2445699350390395*ux_r[2]*dx1-0.2445699350390395*ux_l[2]*dx1-0.3518228202874282*ux_r[1]*dx1-0.3518228202874282*ux_l[1]*dx1+0.7036456405748563*ux_c[1]*dx1+0.25*ux_r[0]*dx1-0.25*ux_l[0]*dx1; 
  grad_u_x[1] = 0.4236075534914363*ux_r[2]*dx1+0.4236075534914363*ux_l[2]*dx1+0.8472151069828725*ux_c[2]*dx1-0.609375*ux_r[1]*dx1+0.609375*ux_l[1]*dx1+0.4330127018922193*ux_r[0]*dx1+0.4330127018922193*ux_l[0]*dx1-0.8660254037844386*ux_c[0]*dx1; 
  grad_u_x[2] = 0.546875*ux_r[2]*dx1-0.546875*ux_l[2]*dx1-0.7866997421983816*ux_r[1]*dx1-0.7866997421983816*ux_l[1]*dx1-2.299583861810654*ux_c[1]*dx1+0.5590169943749475*ux_r[0]*dx1-0.5590169943749475*ux_l[0]*dx1; 

  grad_u_y[0] = 0.2445699350390395*uy_r[2]*dx1-0.2445699350390395*uy_l[2]*dx1-0.3518228202874282*uy_r[1]*dx1-0.3518228202874282*uy_l[1]*dx1+0.7036456405748563*uy_c[1]*dx1+0.25*uy_r[0]*dx1-0.25*uy_l[0]*dx1; 
  grad_u_y[1] = 0.4236075534914363*uy_r[2]*dx1+0.4236075534914363*uy_l[2]*dx1+0.8472151069828725*uy_c[2]*dx1-0.609375*uy_r[1]*dx1+0.609375*uy_l[1]*dx1+0.4330127018922193*uy_r[0]*dx1+0.4330127018922193*uy_l[0]*dx1-0.8660254037844386*uy_c[0]*dx1; 
  grad_u_y[2] = 0.546875*uy_r[2]*dx1-0.546875*uy_l[2]*dx1-0.7866997421983816*uy_r[1]*dx1-0.7866997421983816*uy_l[1]*dx1-2.299583861810654*uy_c[1]*dx1+0.5590169943749475*uy_r[0]*dx1-0.5590169943749475*uy_l[0]*dx1; 

  grad_u_z[0] = 0.2445699350390395*uz_r[2]*dx1-0.2445699350390395*uz_l[2]*dx1-0.3518228202874282*uz_r[1]*dx1-0.3518228202874282*uz_l[1]*dx1+0.7036456405748563*uz_c[1]*dx1+0.25*uz_r[0]*dx1-0.25*uz_l[0]*dx1; 
  grad_u_z[1] = 0.4236075534914363*uz_r[2]*dx1+0.4236075534914363*uz_l[2]*dx1+0.8472151069828725*uz_c[2]*dx1-0.609375*uz_r[1]*dx1+0.609375*uz_l[1]*dx1+0.4330127018922193*uz_r[0]*dx1+0.4330127018922193*uz_l[0]*dx1-0.8660254037844386*uz_c[0]*dx1; 
  grad_u_z[2] = 0.546875*uz_r[2]*dx1-0.546875*uz_l[2]*dx1-0.7866997421983816*uz_r[1]*dx1-0.7866997421983816*uz_l[1]*dx1-2.299583861810654*uz_c[1]*dx1+0.5590169943749475*uz_r[0]*dx1-0.5590169943749475*uz_l[0]*dx1; 

  double ppar_b_l[3] = {0.0}; 
  double ppar_b_c[3] = {0.0}; 
  double ppar_b_r[3] = {0.0}; 
  double div_b_comp[3] = {0.0}; 
  double bb_grad_u_comp[3] = {0.0}; 
  ppar_b_l[0] = 0.7071067811865475*b_l[2]*ppar_l[2]+0.7071067811865475*b_l[1]*ppar_l[1]+0.7071067811865475*b_l[0]*ppar_l[0]; 
  ppar_b_l[1] = 0.6324555320336759*b_l[1]*ppar_l[2]+0.6324555320336759*ppar_l[1]*b_l[2]+0.7071067811865475*b_l[0]*ppar_l[1]+0.7071067811865475*ppar_l[0]*b_l[1]; 
  ppar_b_l[2] = 0.4517539514526256*b_l[2]*ppar_l[2]+0.7071067811865475*b_l[0]*ppar_l[2]+0.7071067811865475*ppar_l[0]*b_l[2]+0.6324555320336759*b_l[1]*ppar_l[1]; 

  ppar_b_c[0] = 0.7071067811865475*b_c[2]*ppar_c[2]+0.7071067811865475*b_c[1]*ppar_c[1]+0.7071067811865475*b_c[0]*ppar_c[0]; 
  ppar_b_c[1] = 0.6324555320336759*b_c[1]*ppar_c[2]+0.6324555320336759*ppar_c[1]*b_c[2]+0.7071067811865475*b_c[0]*ppar_c[1]+0.7071067811865475*ppar_c[0]*b_c[1]; 
  ppar_b_c[2] = 0.4517539514526256*b_c[2]*ppar_c[2]+0.7071067811865475*b_c[0]*ppar_c[2]+0.7071067811865475*ppar_c[0]*b_c[2]+0.6324555320336759*b_c[1]*ppar_c[1]; 

  ppar_b_r[0] = 0.7071067811865475*b_r[2]*ppar_r[2]+0.7071067811865475*b_r[1]*ppar_r[1]+0.7071067811865475*b_r[0]*ppar_r[0]; 
  ppar_b_r[1] = 0.6324555320336759*b_r[1]*ppar_r[2]+0.6324555320336759*ppar_r[1]*b_r[2]+0.7071067811865475*b_r[0]*ppar_r[1]+0.7071067811865475*ppar_r[0]*b_r[1]; 
  ppar_b_r[2] = 0.4517539514526256*b_r[2]*ppar_r[2]+0.7071067811865475*b_r[0]*ppar_r[2]+0.7071067811865475*ppar_r[0]*b_r[2]+0.6324555320336759*b_r[1]*ppar_r[1]; 

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

  div_p_x[0] += ((-160.6824473206489*(rhoux_r[2]+rhoux_l[2]))+805.6133660185961*rhoux_c[2]+207.4401475002413*rhoux_r[1]-207.4401475002413*rhoux_l[1]-137.8125*(rhoux_r[0]+rhoux_l[0])+275.625*rhoux_c[0])*dx14*nuHyp+(0.2445699350390395*Pxx_r[2]-0.2445699350390395*Pxx_l[2]-0.3518228202874282*(Pxx_r[1]+Pxx_l[1])+0.7036456405748563*Pxx_c[1]+0.25*Pxx_r[0]-0.25*Pxx_l[0])*dx1; 
  div_p_x[1] += ((-48.29126109802372*rhoux_r[2])+48.29126109802372*rhoux_l[2]+55.78125*(rhoux_r[1]+rhoux_l[1])+124.6875*rhoux_c[1]-34.09975027401226*rhoux_r[0]+34.09975027401226*rhoux_l[0])*dx14*nuHyp+(0.4236075534914363*(Pxx_r[2]+Pxx_l[2])+0.8472151069828725*Pxx_c[2]-0.609375*Pxx_r[1]+0.609375*Pxx_l[1]+0.4330127018922193*(Pxx_r[0]+Pxx_l[0])-0.8660254037844386*Pxx_c[0])*dx1; 
  div_p_x[2] += ((-250.03125*(rhoux_r[2]+rhoux_l[2]))+1326.9375*rhoux_c[2]+327.8722464023716*rhoux_r[1]-327.8722464023716*rhoux_l[1]-220.1129415351356*(rhoux_r[0]+rhoux_l[0])+440.2258830702712*rhoux_c[0])*dx14*nuHyp+(0.546875*Pxx_r[2]-0.546875*Pxx_l[2]-0.7866997421983816*(Pxx_r[1]+Pxx_l[1])-2.299583861810654*Pxx_c[1]+0.5590169943749475*Pxx_r[0]-0.5590169943749475*Pxx_l[0])*dx1; 

  div_p_y[0] += ((-160.6824473206489*(rhouy_r[2]+rhouy_l[2]))+805.6133660185961*rhouy_c[2]+207.4401475002413*rhouy_r[1]-207.4401475002413*rhouy_l[1]-137.8125*(rhouy_r[0]+rhouy_l[0])+275.625*rhouy_c[0])*dx14*nuHyp+(0.2445699350390395*Pxy_r[2]-0.2445699350390395*Pxy_l[2]-0.3518228202874282*(Pxy_r[1]+Pxy_l[1])+0.7036456405748563*Pxy_c[1]+0.25*Pxy_r[0]-0.25*Pxy_l[0])*dx1; 
  div_p_y[1] += ((-48.29126109802372*rhouy_r[2])+48.29126109802372*rhouy_l[2]+55.78125*(rhouy_r[1]+rhouy_l[1])+124.6875*rhouy_c[1]-34.09975027401226*rhouy_r[0]+34.09975027401226*rhouy_l[0])*dx14*nuHyp+(0.4236075534914363*(Pxy_r[2]+Pxy_l[2])+0.8472151069828725*Pxy_c[2]-0.609375*Pxy_r[1]+0.609375*Pxy_l[1]+0.4330127018922193*(Pxy_r[0]+Pxy_l[0])-0.8660254037844386*Pxy_c[0])*dx1; 
  div_p_y[2] += ((-250.03125*(rhouy_r[2]+rhouy_l[2]))+1326.9375*rhouy_c[2]+327.8722464023716*rhouy_r[1]-327.8722464023716*rhouy_l[1]-220.1129415351356*(rhouy_r[0]+rhouy_l[0])+440.2258830702712*rhouy_c[0])*dx14*nuHyp+(0.546875*Pxy_r[2]-0.546875*Pxy_l[2]-0.7866997421983816*(Pxy_r[1]+Pxy_l[1])-2.299583861810654*Pxy_c[1]+0.5590169943749475*Pxy_r[0]-0.5590169943749475*Pxy_l[0])*dx1; 

  div_p_z[0] += ((-160.6824473206489*(rhouz_r[2]+rhouz_l[2]))+805.6133660185961*rhouz_c[2]+207.4401475002413*rhouz_r[1]-207.4401475002413*rhouz_l[1]-137.8125*(rhouz_r[0]+rhouz_l[0])+275.625*rhouz_c[0])*dx14*nuHyp+(0.2445699350390395*Pxz_r[2]-0.2445699350390395*Pxz_l[2]-0.3518228202874282*(Pxz_r[1]+Pxz_l[1])+0.7036456405748563*Pxz_c[1]+0.25*Pxz_r[0]-0.25*Pxz_l[0])*dx1; 
  div_p_z[1] += ((-48.29126109802372*rhouz_r[2])+48.29126109802372*rhouz_l[2]+55.78125*(rhouz_r[1]+rhouz_l[1])+124.6875*rhouz_c[1]-34.09975027401226*rhouz_r[0]+34.09975027401226*rhouz_l[0])*dx14*nuHyp+(0.4236075534914363*(Pxz_r[2]+Pxz_l[2])+0.8472151069828725*Pxz_c[2]-0.609375*Pxz_r[1]+0.609375*Pxz_l[1]+0.4330127018922193*(Pxz_r[0]+Pxz_l[0])-0.8660254037844386*Pxz_c[0])*dx1; 
  div_p_z[2] += ((-250.03125*(rhouz_r[2]+rhouz_l[2]))+1326.9375*rhouz_c[2]+327.8722464023716*rhouz_r[1]-327.8722464023716*rhouz_l[1]-220.1129415351356*(rhouz_r[0]+rhouz_l[0])+440.2258830702712*rhouz_c[0])*dx14*nuHyp+(0.546875*Pxz_r[2]-0.546875*Pxz_l[2]-0.7866997421983816*(Pxz_r[1]+Pxz_l[1])-2.299583861810654*Pxz_c[1]+0.5590169943749475*Pxz_r[0]-0.5590169943749475*Pxz_l[0])*dx1; 

  double div_ppar_b[3] = {0.0}; 
  double rho_inv[3] = {0.0}; 
  double T_perp_over_m_inv[3] = {0.0}; 
  ser_1x_p2_inv(rho, rho_inv); 
  ser_1x_p2_inv(T_perp_over_m,T_perp_over_m_inv); 
  div_ppar_b[0] = 0.2445699350390395*ppar_b_r[2]*dx1-0.2445699350390395*ppar_b_l[2]*dx1-0.3518228202874282*ppar_b_r[1]*dx1-0.3518228202874282*ppar_b_l[1]*dx1+0.7036456405748563*ppar_b_c[1]*dx1+0.25*ppar_b_r[0]*dx1-0.25*ppar_b_l[0]*dx1; 
  div_ppar_b[1] = 0.4236075534914363*ppar_b_r[2]*dx1+0.4236075534914363*ppar_b_l[2]*dx1+0.8472151069828725*ppar_b_c[2]*dx1-0.609375*ppar_b_r[1]*dx1+0.609375*ppar_b_l[1]*dx1+0.4330127018922193*ppar_b_r[0]*dx1+0.4330127018922193*ppar_b_l[0]*dx1-0.8660254037844386*ppar_b_c[0]*dx1; 
  div_ppar_b[2] = 0.546875*ppar_b_r[2]*dx1-0.546875*ppar_b_l[2]*dx1-0.7866997421983816*ppar_b_r[1]*dx1-0.7866997421983816*ppar_b_l[1]*dx1-2.299583861810654*ppar_b_c[1]*dx1+0.5590169943749475*ppar_b_r[0]*dx1-0.5590169943749475*ppar_b_l[0]*dx1; 

  p_force_F_0[0] += 0.7071067811865475*div_ppar_b[2]*rho_inv[2]-0.7071067811865475*T_perp_over_m[2]*div_b_comp[2]+0.7071067811865475*div_ppar_b[1]*rho_inv[1]-0.7071067811865475*T_perp_over_m[1]*div_b_comp[1]+0.7071067811865475*div_ppar_b[0]*rho_inv[0]-0.7071067811865475*T_perp_over_m[0]*div_b_comp[0]; 
  p_force_F_0[1] += 0.6324555320336759*(div_ppar_b[1]*rho_inv[2]+rho_inv[1]*div_ppar_b[2])-0.6324555320336759*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2])+0.7071067811865475*(div_ppar_b[0]*rho_inv[1]+rho_inv[0]*div_ppar_b[1])-0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force_F_0[2] += 0.4517539514526256*div_ppar_b[2]*rho_inv[2]+0.7071067811865475*(div_ppar_b[0]*rho_inv[2]+rho_inv[0]*div_ppar_b[2])-0.4517539514526256*T_perp_over_m[2]*div_b_comp[2]-0.7071067811865475*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2])+0.6324555320336759*div_ppar_b[1]*rho_inv[1]-0.6324555320336759*T_perp_over_m[1]*div_b_comp[1]; 
  p_force_G_1[0] += 0.7071067811865475*(div_ppar_b[2]*rho_inv[2]+T_perp_over_m[2]*div_b_comp[2]+div_ppar_b[1]*rho_inv[1]+T_perp_over_m[1]*div_b_comp[1]+div_ppar_b[0]*rho_inv[0]+T_perp_over_m[0]*div_b_comp[0]); 
  p_force_G_1[1] += 0.6324555320336759*(div_ppar_b[1]*rho_inv[2]+rho_inv[1]*div_ppar_b[2]+T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2])+0.7071067811865475*(div_ppar_b[0]*rho_inv[1]+rho_inv[0]*div_ppar_b[1]+T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_force_G_1[2] += 0.4517539514526256*div_ppar_b[2]*rho_inv[2]+0.7071067811865475*(div_ppar_b[0]*rho_inv[2]+rho_inv[0]*div_ppar_b[2])+0.4517539514526256*T_perp_over_m[2]*div_b_comp[2]+0.7071067811865475*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2])+0.6324555320336759*(div_ppar_b[1]*rho_inv[1]+T_perp_over_m[1]*div_b_comp[1]); 

  p_perp_source[0] += 0.7071067811865475*(T_perp_over_m_inv[2]*nu_vthsq[2]+T_perp_over_m_inv[1]*nu_vthsq[1]+T_perp_over_m_inv[0]*nu_vthsq[0])-1.0*(nu[0]+grad_u_x[0])+bb_grad_u_comp[0]; 
  p_perp_source[1] += 0.6324555320336759*T_perp_over_m_inv[1]*nu_vthsq[2]+nu_vthsq[1]*(0.6324555320336759*T_perp_over_m_inv[2]+0.7071067811865475*T_perp_over_m_inv[0])-1.0*(nu[1]+grad_u_x[1])+bb_grad_u_comp[1]+0.7071067811865475*nu_vthsq[0]*T_perp_over_m_inv[1]; 
  p_perp_source[2] += (0.4517539514526256*T_perp_over_m_inv[2]+0.7071067811865475*T_perp_over_m_inv[0])*nu_vthsq[2]-1.0*(nu[2]+grad_u_x[2])+bb_grad_u_comp[2]+0.7071067811865475*nu_vthsq[0]*T_perp_over_m_inv[2]+0.6324555320336759*T_perp_over_m_inv[1]*nu_vthsq[1]; 

  p_perp_div_b[0] += 0.7071067811865475*(T_perp_over_m[2]*div_b_comp[2]+T_perp_over_m[1]*div_b_comp[1]+T_perp_over_m[0]*div_b_comp[0]); 
  p_perp_div_b[1] += 0.6324555320336759*(T_perp_over_m[1]*div_b_comp[2]+div_b_comp[1]*T_perp_over_m[2])+0.7071067811865475*(T_perp_over_m[0]*div_b_comp[1]+div_b_comp[0]*T_perp_over_m[1]); 
  p_perp_div_b[2] += 0.4517539514526256*T_perp_over_m[2]*div_b_comp[2]+0.7071067811865475*(T_perp_over_m[0]*div_b_comp[2]+div_b_comp[0]*T_perp_over_m[2])+0.6324555320336759*T_perp_over_m[1]*div_b_comp[1]; 

} 
