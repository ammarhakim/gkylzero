#include <gkyl_vlasov_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_dist_div_ppar_y_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
     const double *fl, const double *fc, const double *fr, 
     const double *bvar_c, const double *pkpm_max_b, double* GKYL_RESTRICT pkpm_div_ppar) 
{ 
  // w[NDIM]:         Cell-center coordinates.
  // dxv[NDIM]:       Cell spacing.
  // bvar_surf_l/c/r: Input surface magnetic field unit vector and tensor in left/center/right cells in each direction.
  // fl/fc/fr:        Input distribution functions [F_0, T_perp/m G = T_perp/m (F_0 - F_1)] in left/center/right cells.
  // bvar_c:          Input volume expansion of magnetic field unit vector and tensor in center cell.
  // pkpm_max_b:      Input surface expansion of max |b| for Lax penalization of streaming: lambda_i = |b_i|.
  // pkpm_div_ppar:   Increment to volume expansion of div(p_par b).
  const double dx1 = 2.0/dxv[1]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double volFact = dxv[2]/2.0; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  double F_0_lr[6] = {0.0}; 
  double F_0_cl[6] = {0.0}; 
  double F_0_cr[6] = {0.0}; 
  double F_0_rl[6] = {0.0}; 
  double Ghat_F_0_vpar_l[4] = {0.0}; 
  double Ghat_F_0_vpar_r[4] = {0.0}; 

  F_0_lr[0] = 1.224744871391589*F_0l[2]+0.7071067811865475*F_0l[0]; 
  F_0_lr[1] = 1.224744871391589*F_0l[4]+0.7071067811865475*F_0l[1]; 
  F_0_lr[2] = 1.224744871391589*F_0l[6]+0.7071067811865475*F_0l[3]; 
  F_0_lr[3] = 1.224744871391589*F_0l[7]+0.7071067811865475*F_0l[5]; 
  F_0_lr[4] = 1.224744871391589*F_0l[10]+0.7071067811865475*F_0l[8]; 
  F_0_lr[5] = 1.224744871391589*F_0l[11]+0.7071067811865475*F_0l[9]; 

  F_0_cl[0] = 0.7071067811865475*F_0c[0]-1.224744871391589*F_0c[2]; 
  F_0_cl[1] = 0.7071067811865475*F_0c[1]-1.224744871391589*F_0c[4]; 
  F_0_cl[2] = 0.7071067811865475*F_0c[3]-1.224744871391589*F_0c[6]; 
  F_0_cl[3] = 0.7071067811865475*F_0c[5]-1.224744871391589*F_0c[7]; 
  F_0_cl[4] = 0.7071067811865475*F_0c[8]-1.224744871391589*F_0c[10]; 
  F_0_cl[5] = 0.7071067811865475*F_0c[9]-1.224744871391589*F_0c[11]; 

  F_0_cr[0] = 1.224744871391589*F_0c[2]+0.7071067811865475*F_0c[0]; 
  F_0_cr[1] = 1.224744871391589*F_0c[4]+0.7071067811865475*F_0c[1]; 
  F_0_cr[2] = 1.224744871391589*F_0c[6]+0.7071067811865475*F_0c[3]; 
  F_0_cr[3] = 1.224744871391589*F_0c[7]+0.7071067811865475*F_0c[5]; 
  F_0_cr[4] = 1.224744871391589*F_0c[10]+0.7071067811865475*F_0c[8]; 
  F_0_cr[5] = 1.224744871391589*F_0c[11]+0.7071067811865475*F_0c[9]; 

  F_0_rl[0] = 0.7071067811865475*F_0r[0]-1.224744871391589*F_0r[2]; 
  F_0_rl[1] = 0.7071067811865475*F_0r[1]-1.224744871391589*F_0r[4]; 
  F_0_rl[2] = 0.7071067811865475*F_0r[3]-1.224744871391589*F_0r[6]; 
  F_0_rl[3] = 0.7071067811865475*F_0r[5]-1.224744871391589*F_0r[7]; 
  F_0_rl[4] = 0.7071067811865475*F_0r[8]-1.224744871391589*F_0r[10]; 
  F_0_rl[5] = 0.7071067811865475*F_0r[9]-1.224744871391589*F_0r[11]; 

  const double *b_surf_lr = &bvar_surf_l[18]; 
  const double *b_surf_cl = &bvar_surf_c[16]; 
  const double *b_surf_cr = &bvar_surf_c[18]; 
  const double *b_surf_rl = &bvar_surf_r[16]; 

  const double *pkpm_max_b_l = &pkpm_max_b[4]; 
  const double *pkpm_max_b_r = &pkpm_max_b[6]; 

  const double *b_c = &bvar_c[4]; 
  double alpha_c[12] = {0.0}; 
  alpha_c[0] = 1.414213562373095*b_c[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*b_c[1]*wvpar; 
  alpha_c[2] = 1.414213562373095*b_c[2]*wvpar; 
  alpha_c[3] = 0.408248290463863*b_c[0]*dvpar; 
  alpha_c[4] = 1.414213562373095*b_c[3]*wvpar; 
  alpha_c[5] = 0.408248290463863*b_c[1]*dvpar; 
  alpha_c[6] = 0.408248290463863*b_c[2]*dvpar; 
  alpha_c[7] = 0.408248290463863*b_c[3]*dvpar; 

  if (wvpar>0) { 

  Ghat_F_0_vpar_l[0] = ((0.3535533905932737*F_0_lr[1]-0.3535533905932737*F_0_cl[1])*pkpm_max_b_l[1]+0.1767766952966368*((F_0_lr[1]+F_0_cl[1])*b_surf_lr[1]+(F_0_lr[1]+F_0_cl[1])*b_surf_cl[1])+(0.3535533905932737*F_0_lr[0]-0.3535533905932737*F_0_cl[0])*pkpm_max_b_l[0]+0.1767766952966368*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[0]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[0]))*wvpar+((0.1020620726159657*pkpm_max_b_l[1]+0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1]))*F_0_lr[3]+(0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1])-0.1020620726159657*pkpm_max_b_l[1])*F_0_cl[3]+(0.1020620726159657*pkpm_max_b_l[0]+0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[2]+(0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0])-0.1020620726159657*pkpm_max_b_l[0])*F_0_cl[2])*dvpar; 
  Ghat_F_0_vpar_l[1] = ((0.3535533905932737*F_0_lr[0]-0.3535533905932737*F_0_cl[0])*pkpm_max_b_l[1]+0.1767766952966368*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[1]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[1])+(0.3535533905932737*pkpm_max_b_l[0]+0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[1]+(0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0])-0.3535533905932737*pkpm_max_b_l[0])*F_0_cl[1])*wvpar+((0.1020620726159657*pkpm_max_b_l[0]+0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[3]+(0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0])-0.1020620726159657*pkpm_max_b_l[0])*F_0_cl[3]+(0.1020620726159657*pkpm_max_b_l[1]+0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1]))*F_0_lr[2]+(0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1])-0.1020620726159657*pkpm_max_b_l[1])*F_0_cl[2])*dvpar; 
  Ghat_F_0_vpar_l[2] = ((0.3535533905932737*pkpm_max_b_l[1]+0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1]))*F_0_lr[3]+(0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1])-0.3535533905932737*pkpm_max_b_l[1])*F_0_cl[3]+(0.3535533905932737*pkpm_max_b_l[0]+0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[2]+(0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0])-0.3535533905932737*pkpm_max_b_l[0])*F_0_cl[2])*wvpar+((0.09128709291752765*pkpm_max_b_l[1]+0.04564354645876382*(b_surf_lr[1]+b_surf_cl[1]))*F_0_lr[5]+(0.04564354645876382*(b_surf_lr[1]+b_surf_cl[1])-0.09128709291752765*pkpm_max_b_l[1])*F_0_cl[5]+(0.09128709291752767*pkpm_max_b_l[0]+0.04564354645876383*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[4]+(0.04564354645876383*(b_surf_lr[0]+b_surf_cl[0])-0.09128709291752767*pkpm_max_b_l[0])*F_0_cl[4]+(0.1020620726159657*F_0_lr[1]-0.1020620726159657*F_0_cl[1])*pkpm_max_b_l[1]+0.05103103630798286*((F_0_lr[1]+F_0_cl[1])*b_surf_lr[1]+(F_0_lr[1]+F_0_cl[1])*b_surf_cl[1])+(0.1020620726159657*F_0_lr[0]-0.1020620726159657*F_0_cl[0])*pkpm_max_b_l[0]+0.05103103630798286*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[0]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[0]))*dvpar; 
  Ghat_F_0_vpar_l[3] = ((0.3535533905932737*pkpm_max_b_l[0]+0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[3]+(0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0])-0.3535533905932737*pkpm_max_b_l[0])*F_0_cl[3]+(0.3535533905932737*pkpm_max_b_l[1]+0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1]))*F_0_lr[2]+(0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1])-0.3535533905932737*pkpm_max_b_l[1])*F_0_cl[2])*wvpar+((0.09128709291752765*pkpm_max_b_l[0]+0.04564354645876382*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[5]+(0.04564354645876382*(b_surf_lr[0]+b_surf_cl[0])-0.09128709291752765*pkpm_max_b_l[0])*F_0_cl[5]+(0.09128709291752767*pkpm_max_b_l[1]+0.04564354645876383*(b_surf_lr[1]+b_surf_cl[1]))*F_0_lr[4]+(0.04564354645876383*(b_surf_lr[1]+b_surf_cl[1])-0.09128709291752767*pkpm_max_b_l[1])*F_0_cl[4]+(0.1020620726159657*F_0_lr[0]-0.1020620726159657*F_0_cl[0])*pkpm_max_b_l[1]+0.05103103630798286*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[1]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[1])+(0.1020620726159657*pkpm_max_b_l[0]+0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0]))*F_0_lr[1]+(0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0])-0.1020620726159657*pkpm_max_b_l[0])*F_0_cl[1])*dvpar; 

  Ghat_F_0_vpar_r[0] = ((0.3535533905932737*F_0_cr[1]-0.3535533905932737*F_0_rl[1])*pkpm_max_b_r[1]+0.1767766952966368*((F_0_rl[1]+F_0_cr[1])*b_surf_rl[1]+(F_0_rl[1]+F_0_cr[1])*b_surf_cr[1])+(0.3535533905932737*F_0_cr[0]-0.3535533905932737*F_0_rl[0])*pkpm_max_b_r[0]+0.1767766952966368*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[0]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[0]))*wvpar+((0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1])-0.1020620726159657*pkpm_max_b_r[1])*F_0_rl[3]+(0.1020620726159657*pkpm_max_b_r[1]+0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1]))*F_0_cr[3]+(0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0])-0.1020620726159657*pkpm_max_b_r[0])*F_0_rl[2]+(0.1020620726159657*pkpm_max_b_r[0]+0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[2])*dvpar; 
  Ghat_F_0_vpar_r[1] = ((0.3535533905932737*F_0_cr[0]-0.3535533905932737*F_0_rl[0])*pkpm_max_b_r[1]+0.1767766952966368*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[1]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[1])+(0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0])-0.3535533905932737*pkpm_max_b_r[0])*F_0_rl[1]+(0.3535533905932737*pkpm_max_b_r[0]+0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[1])*wvpar+((0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0])-0.1020620726159657*pkpm_max_b_r[0])*F_0_rl[3]+(0.1020620726159657*pkpm_max_b_r[0]+0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[3]+(0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1])-0.1020620726159657*pkpm_max_b_r[1])*F_0_rl[2]+(0.1020620726159657*pkpm_max_b_r[1]+0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1]))*F_0_cr[2])*dvpar; 
  Ghat_F_0_vpar_r[2] = ((0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1])-0.3535533905932737*pkpm_max_b_r[1])*F_0_rl[3]+(0.3535533905932737*pkpm_max_b_r[1]+0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1]))*F_0_cr[3]+(0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0])-0.3535533905932737*pkpm_max_b_r[0])*F_0_rl[2]+(0.3535533905932737*pkpm_max_b_r[0]+0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[2])*wvpar+((0.04564354645876382*(b_surf_rl[1]+b_surf_cr[1])-0.09128709291752765*pkpm_max_b_r[1])*F_0_rl[5]+(0.09128709291752765*pkpm_max_b_r[1]+0.04564354645876382*(b_surf_rl[1]+b_surf_cr[1]))*F_0_cr[5]+(0.04564354645876383*(b_surf_rl[0]+b_surf_cr[0])-0.09128709291752767*pkpm_max_b_r[0])*F_0_rl[4]+(0.09128709291752767*pkpm_max_b_r[0]+0.04564354645876383*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[4]+(0.1020620726159657*F_0_cr[1]-0.1020620726159657*F_0_rl[1])*pkpm_max_b_r[1]+0.05103103630798286*((F_0_rl[1]+F_0_cr[1])*b_surf_rl[1]+(F_0_rl[1]+F_0_cr[1])*b_surf_cr[1])+(0.1020620726159657*F_0_cr[0]-0.1020620726159657*F_0_rl[0])*pkpm_max_b_r[0]+0.05103103630798286*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[0]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[0]))*dvpar; 
  Ghat_F_0_vpar_r[3] = ((0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0])-0.3535533905932737*pkpm_max_b_r[0])*F_0_rl[3]+(0.3535533905932737*pkpm_max_b_r[0]+0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[3]+(0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1])-0.3535533905932737*pkpm_max_b_r[1])*F_0_rl[2]+(0.3535533905932737*pkpm_max_b_r[1]+0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1]))*F_0_cr[2])*wvpar+((0.04564354645876382*(b_surf_rl[0]+b_surf_cr[0])-0.09128709291752765*pkpm_max_b_r[0])*F_0_rl[5]+(0.09128709291752765*pkpm_max_b_r[0]+0.04564354645876382*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[5]+(0.04564354645876383*(b_surf_rl[1]+b_surf_cr[1])-0.09128709291752767*pkpm_max_b_r[1])*F_0_rl[4]+(0.09128709291752767*pkpm_max_b_r[1]+0.04564354645876383*(b_surf_rl[1]+b_surf_cr[1]))*F_0_cr[4]+(0.1020620726159657*F_0_cr[0]-0.1020620726159657*F_0_rl[0])*pkpm_max_b_r[1]+0.05103103630798286*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[1]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[1])+(0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0])-0.1020620726159657*pkpm_max_b_r[0])*F_0_rl[1]+(0.1020620726159657*pkpm_max_b_r[0]+0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0]))*F_0_cr[1])*dvpar; 

  } else { 

  Ghat_F_0_vpar_l[0] = ((0.3535533905932737*F_0_cl[1]-0.3535533905932737*F_0_lr[1])*pkpm_max_b_l[1]+0.1767766952966368*((F_0_lr[1]+F_0_cl[1])*b_surf_lr[1]+(F_0_lr[1]+F_0_cl[1])*b_surf_cl[1])+(0.3535533905932737*F_0_cl[0]-0.3535533905932737*F_0_lr[0])*pkpm_max_b_l[0]+0.1767766952966368*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[0]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[0]))*wvpar+((0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1])-0.1020620726159657*pkpm_max_b_l[1])*F_0_lr[3]+(0.1020620726159657*pkpm_max_b_l[1]+0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1]))*F_0_cl[3]+(0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0])-0.1020620726159657*pkpm_max_b_l[0])*F_0_lr[2]+(0.1020620726159657*pkpm_max_b_l[0]+0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[2])*dvpar; 
  Ghat_F_0_vpar_l[1] = ((0.3535533905932737*F_0_cl[0]-0.3535533905932737*F_0_lr[0])*pkpm_max_b_l[1]+0.1767766952966368*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[1]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[1])+(0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0])-0.3535533905932737*pkpm_max_b_l[0])*F_0_lr[1]+(0.3535533905932737*pkpm_max_b_l[0]+0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[1])*wvpar+((0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0])-0.1020620726159657*pkpm_max_b_l[0])*F_0_lr[3]+(0.1020620726159657*pkpm_max_b_l[0]+0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[3]+(0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1])-0.1020620726159657*pkpm_max_b_l[1])*F_0_lr[2]+(0.1020620726159657*pkpm_max_b_l[1]+0.05103103630798286*(b_surf_lr[1]+b_surf_cl[1]))*F_0_cl[2])*dvpar; 
  Ghat_F_0_vpar_l[2] = ((0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1])-0.3535533905932737*pkpm_max_b_l[1])*F_0_lr[3]+(0.3535533905932737*pkpm_max_b_l[1]+0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1]))*F_0_cl[3]+(0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0])-0.3535533905932737*pkpm_max_b_l[0])*F_0_lr[2]+(0.3535533905932737*pkpm_max_b_l[0]+0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[2])*wvpar+((0.04564354645876382*(b_surf_lr[1]+b_surf_cl[1])-0.09128709291752765*pkpm_max_b_l[1])*F_0_lr[5]+(0.09128709291752765*pkpm_max_b_l[1]+0.04564354645876382*(b_surf_lr[1]+b_surf_cl[1]))*F_0_cl[5]+(0.04564354645876383*(b_surf_lr[0]+b_surf_cl[0])-0.09128709291752767*pkpm_max_b_l[0])*F_0_lr[4]+(0.09128709291752767*pkpm_max_b_l[0]+0.04564354645876383*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[4]+(0.1020620726159657*F_0_cl[1]-0.1020620726159657*F_0_lr[1])*pkpm_max_b_l[1]+0.05103103630798286*((F_0_lr[1]+F_0_cl[1])*b_surf_lr[1]+(F_0_lr[1]+F_0_cl[1])*b_surf_cl[1])+(0.1020620726159657*F_0_cl[0]-0.1020620726159657*F_0_lr[0])*pkpm_max_b_l[0]+0.05103103630798286*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[0]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[0]))*dvpar; 
  Ghat_F_0_vpar_l[3] = ((0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0])-0.3535533905932737*pkpm_max_b_l[0])*F_0_lr[3]+(0.3535533905932737*pkpm_max_b_l[0]+0.1767766952966368*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[3]+(0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1])-0.3535533905932737*pkpm_max_b_l[1])*F_0_lr[2]+(0.3535533905932737*pkpm_max_b_l[1]+0.1767766952966368*(b_surf_lr[1]+b_surf_cl[1]))*F_0_cl[2])*wvpar+((0.04564354645876382*(b_surf_lr[0]+b_surf_cl[0])-0.09128709291752765*pkpm_max_b_l[0])*F_0_lr[5]+(0.09128709291752765*pkpm_max_b_l[0]+0.04564354645876382*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[5]+(0.04564354645876383*(b_surf_lr[1]+b_surf_cl[1])-0.09128709291752767*pkpm_max_b_l[1])*F_0_lr[4]+(0.09128709291752767*pkpm_max_b_l[1]+0.04564354645876383*(b_surf_lr[1]+b_surf_cl[1]))*F_0_cl[4]+(0.1020620726159657*F_0_cl[0]-0.1020620726159657*F_0_lr[0])*pkpm_max_b_l[1]+0.05103103630798286*((F_0_lr[0]+F_0_cl[0])*b_surf_lr[1]+(F_0_lr[0]+F_0_cl[0])*b_surf_cl[1])+(0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0])-0.1020620726159657*pkpm_max_b_l[0])*F_0_lr[1]+(0.1020620726159657*pkpm_max_b_l[0]+0.05103103630798286*(b_surf_lr[0]+b_surf_cl[0]))*F_0_cl[1])*dvpar; 

  Ghat_F_0_vpar_r[0] = ((0.3535533905932737*F_0_rl[1]-0.3535533905932737*F_0_cr[1])*pkpm_max_b_r[1]+0.1767766952966368*((F_0_rl[1]+F_0_cr[1])*b_surf_rl[1]+(F_0_rl[1]+F_0_cr[1])*b_surf_cr[1])+(0.3535533905932737*F_0_rl[0]-0.3535533905932737*F_0_cr[0])*pkpm_max_b_r[0]+0.1767766952966368*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[0]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[0]))*wvpar+((0.1020620726159657*pkpm_max_b_r[1]+0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1]))*F_0_rl[3]+(0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1])-0.1020620726159657*pkpm_max_b_r[1])*F_0_cr[3]+(0.1020620726159657*pkpm_max_b_r[0]+0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[2]+(0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0])-0.1020620726159657*pkpm_max_b_r[0])*F_0_cr[2])*dvpar; 
  Ghat_F_0_vpar_r[1] = ((0.3535533905932737*F_0_rl[0]-0.3535533905932737*F_0_cr[0])*pkpm_max_b_r[1]+0.1767766952966368*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[1]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[1])+(0.3535533905932737*pkpm_max_b_r[0]+0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[1]+(0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0])-0.3535533905932737*pkpm_max_b_r[0])*F_0_cr[1])*wvpar+((0.1020620726159657*pkpm_max_b_r[0]+0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[3]+(0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0])-0.1020620726159657*pkpm_max_b_r[0])*F_0_cr[3]+(0.1020620726159657*pkpm_max_b_r[1]+0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1]))*F_0_rl[2]+(0.05103103630798286*(b_surf_rl[1]+b_surf_cr[1])-0.1020620726159657*pkpm_max_b_r[1])*F_0_cr[2])*dvpar; 
  Ghat_F_0_vpar_r[2] = ((0.3535533905932737*pkpm_max_b_r[1]+0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1]))*F_0_rl[3]+(0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1])-0.3535533905932737*pkpm_max_b_r[1])*F_0_cr[3]+(0.3535533905932737*pkpm_max_b_r[0]+0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[2]+(0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0])-0.3535533905932737*pkpm_max_b_r[0])*F_0_cr[2])*wvpar+((0.09128709291752765*pkpm_max_b_r[1]+0.04564354645876382*(b_surf_rl[1]+b_surf_cr[1]))*F_0_rl[5]+(0.04564354645876382*(b_surf_rl[1]+b_surf_cr[1])-0.09128709291752765*pkpm_max_b_r[1])*F_0_cr[5]+(0.09128709291752767*pkpm_max_b_r[0]+0.04564354645876383*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[4]+(0.04564354645876383*(b_surf_rl[0]+b_surf_cr[0])-0.09128709291752767*pkpm_max_b_r[0])*F_0_cr[4]+(0.1020620726159657*F_0_rl[1]-0.1020620726159657*F_0_cr[1])*pkpm_max_b_r[1]+0.05103103630798286*((F_0_rl[1]+F_0_cr[1])*b_surf_rl[1]+(F_0_rl[1]+F_0_cr[1])*b_surf_cr[1])+(0.1020620726159657*F_0_rl[0]-0.1020620726159657*F_0_cr[0])*pkpm_max_b_r[0]+0.05103103630798286*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[0]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[0]))*dvpar; 
  Ghat_F_0_vpar_r[3] = ((0.3535533905932737*pkpm_max_b_r[0]+0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[3]+(0.1767766952966368*(b_surf_rl[0]+b_surf_cr[0])-0.3535533905932737*pkpm_max_b_r[0])*F_0_cr[3]+(0.3535533905932737*pkpm_max_b_r[1]+0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1]))*F_0_rl[2]+(0.1767766952966368*(b_surf_rl[1]+b_surf_cr[1])-0.3535533905932737*pkpm_max_b_r[1])*F_0_cr[2])*wvpar+((0.09128709291752765*pkpm_max_b_r[0]+0.04564354645876382*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[5]+(0.04564354645876382*(b_surf_rl[0]+b_surf_cr[0])-0.09128709291752765*pkpm_max_b_r[0])*F_0_cr[5]+(0.09128709291752767*pkpm_max_b_r[1]+0.04564354645876383*(b_surf_rl[1]+b_surf_cr[1]))*F_0_rl[4]+(0.04564354645876383*(b_surf_rl[1]+b_surf_cr[1])-0.09128709291752767*pkpm_max_b_r[1])*F_0_cr[4]+(0.1020620726159657*F_0_rl[0]-0.1020620726159657*F_0_cr[0])*pkpm_max_b_r[1]+0.05103103630798286*((F_0_rl[0]+F_0_cr[0])*b_surf_rl[1]+(F_0_rl[0]+F_0_cr[0])*b_surf_cr[1])+(0.1020620726159657*pkpm_max_b_r[0]+0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0]))*F_0_rl[1]+(0.05103103630798286*(b_surf_rl[0]+b_surf_cr[0])-0.1020620726159657*pkpm_max_b_r[0])*F_0_cr[1])*dvpar; 

  } 
  pkpm_div_ppar[0] += dx1*volFact*((Ghat_F_0_vpar_r[0]-1.0*Ghat_F_0_vpar_l[0])*wvpar+(0.2886751345948129*Ghat_F_0_vpar_r[2]-0.2886751345948129*Ghat_F_0_vpar_l[2])*dvpar); 
  pkpm_div_ppar[1] += dx1*volFact*((Ghat_F_0_vpar_r[1]-1.0*Ghat_F_0_vpar_l[1])*wvpar+(0.2886751345948129*Ghat_F_0_vpar_r[3]-0.2886751345948129*Ghat_F_0_vpar_l[3])*dvpar); 
  pkpm_div_ppar[2] += dx1*volFact*((1.732050807568877*(Ghat_F_0_vpar_r[0]+Ghat_F_0_vpar_l[0])-0.8660254037844386*(F_0c[7]*alpha_c[7]+F_0c[6]*alpha_c[6]+F_0c[5]*alpha_c[5]+F_0c[4]*alpha_c[4]+F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar+((-0.223606797749979*alpha_c[7]*F_0c[11])-0.223606797749979*(alpha_c[6]*F_0c[10]+alpha_c[5]*F_0c[9])-0.223606797749979*alpha_c[3]*F_0c[8]-0.25*(F_0c[4]*alpha_c[7]+alpha_c[4]*F_0c[7]+F_0c[2]*alpha_c[6]+alpha_c[2]*F_0c[6]+F_0c[1]*alpha_c[5]+alpha_c[1]*F_0c[5]+F_0c[0]*alpha_c[3]+alpha_c[0]*F_0c[3])+0.5*(Ghat_F_0_vpar_r[2]+Ghat_F_0_vpar_l[2]))*dvpar); 
  pkpm_div_ppar[3] += dx1*volFact*(((-0.8660254037844386*(F_0c[6]*alpha_c[7]+alpha_c[6]*F_0c[7]+F_0c[3]*alpha_c[5]+alpha_c[3]*F_0c[5]+F_0c[2]*alpha_c[4]+alpha_c[2]*F_0c[4]+F_0c[0]*alpha_c[1]))+1.732050807568877*(Ghat_F_0_vpar_r[1]+Ghat_F_0_vpar_l[1])-0.8660254037844386*alpha_c[0]*F_0c[1])*wvpar+((-0.223606797749979*alpha_c[6]*F_0c[11])-0.223606797749979*(alpha_c[7]*F_0c[10]+alpha_c[3]*F_0c[9])-0.223606797749979*alpha_c[5]*F_0c[8]-0.25*(F_0c[2]*alpha_c[7]+alpha_c[2]*F_0c[7]+F_0c[4]*alpha_c[6]+alpha_c[4]*F_0c[6]+F_0c[0]*alpha_c[5]+alpha_c[0]*F_0c[5]+F_0c[1]*alpha_c[3])+0.5*(Ghat_F_0_vpar_r[3]+Ghat_F_0_vpar_l[3])-0.25*alpha_c[1]*F_0c[3])*dvpar); 

} 
