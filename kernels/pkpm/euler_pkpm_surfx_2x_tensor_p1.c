#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double euler_pkpm_surfx_2x_tensor_p1(const double *w, const double *dxv,
    const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
    const double *pkpm_u_surf_l, const double *pkpm_u_surf_c, const double *pkpm_u_surf_r,
    const double *p_ij_l, const double *p_ij_c, const double *p_ij_r,
    const double *euler_pkpm_l, const double *euler_pkpm_c, const double *euler_pkpm_r, 
    const double *pkpm_lax, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                Cell-center coordinates.
  // dxv[NDIM]:              Cell spacing.
  // wv_eqn:                 Wave equation for computing fluctuations at the interface for upwinding.
  // geom_l:                 Geometry for the left surface update.
  // geom_r:                 Geometry for the right surface update.
  // vlasov_pkpm_moms_l/c/r: Input pkpm moments in left/center/right cells.
  // pkpm_u_surf_l/c/r:      Input surface flow velocity in left/center/right cells in each direction.
  // p_ij_l/c/r:             Input volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij in left/center/right cells.
  // euler_pkpm_l/c/r:       Input [rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.
  // pkpm_lax:               Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // out: Incremented output.

  const double dx1 = 2.0/dxv[0]; 

  const double *rho_l = &vlasov_pkpm_moms_l[0]; 
  const double *rho_c = &vlasov_pkpm_moms_c[0]; 
  const double *rho_r = &vlasov_pkpm_moms_r[0]; 

  const double *Pxx_l = &p_ij_l[0]; 
  const double *Pxy_l = &p_ij_l[9]; 
  const double *Pxz_l = &p_ij_l[18]; 
  const double *Pyy_l = &p_ij_l[27]; 
  const double *Pyz_l = &p_ij_l[36]; 
  const double *Pzz_l = &p_ij_l[45]; 

  const double *Pxx_c = &p_ij_c[0]; 
  const double *Pxy_c = &p_ij_c[9]; 
  const double *Pxz_c = &p_ij_c[18]; 
  const double *Pyy_c = &p_ij_c[27]; 
  const double *Pyz_c = &p_ij_c[36]; 
  const double *Pzz_c = &p_ij_c[45]; 

  const double *Pxx_r = &p_ij_r[0]; 
  const double *Pxy_r = &p_ij_r[9]; 
  const double *Pxz_r = &p_ij_r[18]; 
  const double *Pyy_r = &p_ij_r[27]; 
  const double *Pyz_r = &p_ij_r[36]; 
  const double *Pzz_r = &p_ij_r[45]; 

  const double *ux_surf_lr = &pkpm_u_surf_l[2]; 
  const double *uy_surf_lr = &pkpm_u_surf_l[6]; 
  const double *uz_surf_lr = &pkpm_u_surf_l[10]; 

  const double *ux_surf_cl = &pkpm_u_surf_c[0]; 
  const double *uy_surf_cl = &pkpm_u_surf_c[4]; 
  const double *uz_surf_cl = &pkpm_u_surf_c[8]; 

  const double *ux_surf_cr = &pkpm_u_surf_c[2]; 
  const double *uy_surf_cr = &pkpm_u_surf_c[6]; 
  const double *uz_surf_cr = &pkpm_u_surf_c[10]; 

  const double *ux_surf_rl = &pkpm_u_surf_r[0]; 
  const double *uy_surf_rl = &pkpm_u_surf_r[4]; 
  const double *uz_surf_rl = &pkpm_u_surf_r[8]; 

  const double *pkpm_lax_l = &pkpm_lax[0]; 
  const double *pkpm_lax_r = &pkpm_lax[3]; 

  double *outrhou0 = &out[0]; 
  double *outrhou1 = &out[9]; 
  double *outrhou2 = &out[18]; 

  double flux_rho_l[3] = {0.0}; 
  double flux_rho_r[3] = {0.0}; 

  double avg_p_ij_x_l[3] = {0.0}; 
  double avg_p_ij_x_r[3] = {0.0}; 
  double avg_p_ij_y_l[3] = {0.0}; 
  double avg_p_ij_y_r[3] = {0.0}; 
  double avg_p_ij_z_l[3] = {0.0}; 
  double avg_p_ij_z_r[3] = {0.0}; 

  flux_rho_l[0] = 0.5590169943749475*pkpm_lax_l[2]*rho_l[8]-0.5590169943749475*pkpm_lax_l[2]*rho_c[8]+0.4330127018922194*pkpm_lax_l[2]*rho_l[7]+0.4330127018922194*pkpm_lax_l[2]*rho_c[7]+0.2795084971874738*ux_surf_lr[1]*rho_l[6]+0.2795084971874738*ux_surf_cl[1]*rho_l[6]+0.5590169943749476*pkpm_lax_l[1]*rho_l[6]+0.2795084971874738*ux_surf_lr[1]*rho_c[6]+0.2795084971874738*ux_surf_cl[1]*rho_c[6]-0.5590169943749476*pkpm_lax_l[1]*rho_c[6]+0.25*pkpm_lax_l[2]*rho_l[5]-0.25*pkpm_lax_l[2]*rho_c[5]+0.2795084971874737*ux_surf_lr[0]*rho_l[4]+0.2795084971874737*ux_surf_cl[0]*rho_l[4]+0.5590169943749475*pkpm_lax_l[0]*rho_l[4]+0.2795084971874737*ux_surf_lr[0]*rho_c[4]+0.2795084971874737*ux_surf_cl[0]*rho_c[4]-0.5590169943749475*pkpm_lax_l[0]*rho_c[4]+0.2165063509461096*ux_surf_lr[1]*rho_l[3]+0.2165063509461096*ux_surf_cl[1]*rho_l[3]+0.4330127018922193*pkpm_lax_l[1]*rho_l[3]-0.2165063509461096*ux_surf_lr[1]*rho_c[3]-0.2165063509461096*ux_surf_cl[1]*rho_c[3]+0.4330127018922193*pkpm_lax_l[1]*rho_c[3]+0.125*ux_surf_lr[1]*rho_l[2]+0.125*ux_surf_cl[1]*rho_l[2]+0.25*pkpm_lax_l[1]*rho_l[2]+0.125*ux_surf_lr[1]*rho_c[2]+0.125*ux_surf_cl[1]*rho_c[2]-0.25*pkpm_lax_l[1]*rho_c[2]+0.2165063509461096*ux_surf_lr[0]*rho_l[1]+0.2165063509461096*ux_surf_cl[0]*rho_l[1]+0.4330127018922193*pkpm_lax_l[0]*rho_l[1]-0.2165063509461096*ux_surf_lr[0]*rho_c[1]-0.2165063509461096*ux_surf_cl[0]*rho_c[1]+0.4330127018922193*pkpm_lax_l[0]*rho_c[1]+0.125*rho_l[0]*ux_surf_lr[0]+0.125*rho_c[0]*ux_surf_lr[0]+0.125*rho_l[0]*ux_surf_cl[0]+0.125*rho_c[0]*ux_surf_cl[0]+0.25*pkpm_lax_l[0]*rho_l[0]-0.25*pkpm_lax_l[0]*rho_c[0]; 
  flux_rho_l[1] = 0.25*ux_surf_lr[1]*rho_l[8]+0.25*ux_surf_cl[1]*rho_l[8]+0.5*pkpm_lax_l[1]*rho_l[8]+0.25*ux_surf_lr[1]*rho_c[8]+0.25*ux_surf_cl[1]*rho_c[8]-0.5*pkpm_lax_l[1]*rho_c[8]+0.1936491673103709*ux_surf_lr[1]*rho_l[7]+0.1936491673103709*ux_surf_cl[1]*rho_l[7]+0.3872983346207417*pkpm_lax_l[1]*rho_l[7]-0.1936491673103709*ux_surf_lr[1]*rho_c[7]-0.1936491673103709*ux_surf_cl[1]*rho_c[7]+0.3872983346207417*pkpm_lax_l[1]*rho_c[7]+0.5000000000000001*pkpm_lax_l[2]*rho_l[6]+0.2795084971874738*ux_surf_lr[0]*rho_l[6]+0.2795084971874738*ux_surf_cl[0]*rho_l[6]+0.5590169943749476*pkpm_lax_l[0]*rho_l[6]-0.5000000000000001*pkpm_lax_l[2]*rho_c[6]+0.2795084971874738*ux_surf_lr[0]*rho_c[6]+0.2795084971874738*ux_surf_cl[0]*rho_c[6]-0.5590169943749476*pkpm_lax_l[0]*rho_c[6]+0.1118033988749895*ux_surf_lr[1]*rho_l[5]+0.1118033988749895*ux_surf_cl[1]*rho_l[5]+0.223606797749979*pkpm_lax_l[1]*rho_l[5]+0.1118033988749895*ux_surf_lr[1]*rho_c[5]+0.1118033988749895*ux_surf_cl[1]*rho_c[5]-0.223606797749979*pkpm_lax_l[1]*rho_c[5]+0.2795084971874737*ux_surf_lr[1]*rho_l[4]+0.2795084971874737*ux_surf_cl[1]*rho_l[4]+0.5590169943749475*pkpm_lax_l[1]*rho_l[4]+0.2795084971874737*ux_surf_lr[1]*rho_c[4]+0.2795084971874737*ux_surf_cl[1]*rho_c[4]-0.5590169943749475*pkpm_lax_l[1]*rho_c[4]+0.3872983346207416*pkpm_lax_l[2]*rho_l[3]+0.2165063509461096*ux_surf_lr[0]*rho_l[3]+0.2165063509461096*ux_surf_cl[0]*rho_l[3]+0.4330127018922193*pkpm_lax_l[0]*rho_l[3]+0.3872983346207416*pkpm_lax_l[2]*rho_c[3]-0.2165063509461096*ux_surf_lr[0]*rho_c[3]-0.2165063509461096*ux_surf_cl[0]*rho_c[3]+0.4330127018922193*pkpm_lax_l[0]*rho_c[3]+0.223606797749979*pkpm_lax_l[2]*rho_l[2]+0.125*ux_surf_lr[0]*rho_l[2]+0.125*ux_surf_cl[0]*rho_l[2]+0.25*pkpm_lax_l[0]*rho_l[2]-0.223606797749979*pkpm_lax_l[2]*rho_c[2]+0.125*ux_surf_lr[0]*rho_c[2]+0.125*ux_surf_cl[0]*rho_c[2]-0.25*pkpm_lax_l[0]*rho_c[2]+0.2165063509461096*rho_l[1]*ux_surf_lr[1]-0.2165063509461096*rho_c[1]*ux_surf_lr[1]+0.125*rho_l[0]*ux_surf_lr[1]+0.125*rho_c[0]*ux_surf_lr[1]+0.2165063509461096*rho_l[1]*ux_surf_cl[1]-0.2165063509461096*rho_c[1]*ux_surf_cl[1]+0.125*rho_l[0]*ux_surf_cl[1]+0.125*rho_c[0]*ux_surf_cl[1]+0.4330127018922193*pkpm_lax_l[1]*rho_l[1]+0.4330127018922193*pkpm_lax_l[1]*rho_c[1]+0.25*rho_l[0]*pkpm_lax_l[1]-0.25*rho_c[0]*pkpm_lax_l[1]; 
  flux_rho_l[2] = 0.3571428571428572*pkpm_lax_l[2]*rho_l[8]+0.2795084971874737*ux_surf_lr[0]*rho_l[8]+0.2795084971874737*ux_surf_cl[0]*rho_l[8]+0.5590169943749475*pkpm_lax_l[0]*rho_l[8]-0.3571428571428572*pkpm_lax_l[2]*rho_c[8]+0.2795084971874737*ux_surf_lr[0]*rho_c[8]+0.2795084971874737*ux_surf_cl[0]*rho_c[8]-0.5590169943749475*pkpm_lax_l[0]*rho_c[8]+0.276641667586244*pkpm_lax_l[2]*rho_l[7]+0.2165063509461097*ux_surf_lr[0]*rho_l[7]+0.2165063509461097*ux_surf_cl[0]*rho_l[7]+0.4330127018922194*pkpm_lax_l[0]*rho_l[7]+0.276641667586244*pkpm_lax_l[2]*rho_c[7]-0.2165063509461097*ux_surf_lr[0]*rho_c[7]-0.2165063509461097*ux_surf_cl[0]*rho_c[7]+0.4330127018922194*pkpm_lax_l[0]*rho_c[7]+0.2500000000000001*ux_surf_lr[1]*rho_l[6]+0.2500000000000001*ux_surf_cl[1]*rho_l[6]+0.5000000000000001*pkpm_lax_l[1]*rho_l[6]+0.2500000000000001*ux_surf_lr[1]*rho_c[6]+0.2500000000000001*ux_surf_cl[1]*rho_c[6]-0.5000000000000001*pkpm_lax_l[1]*rho_c[6]+0.159719141249985*pkpm_lax_l[2]*rho_l[5]+0.125*ux_surf_lr[0]*rho_l[5]+0.125*ux_surf_cl[0]*rho_l[5]+0.25*pkpm_lax_l[0]*rho_l[5]-0.159719141249985*pkpm_lax_l[2]*rho_c[5]+0.125*ux_surf_lr[0]*rho_c[5]+0.125*ux_surf_cl[0]*rho_c[5]-0.25*pkpm_lax_l[0]*rho_c[5]+0.5590169943749475*pkpm_lax_l[2]*rho_l[4]-0.5590169943749475*pkpm_lax_l[2]*rho_c[4]+0.1936491673103708*ux_surf_lr[1]*rho_l[3]+0.1936491673103708*ux_surf_cl[1]*rho_l[3]+0.3872983346207416*pkpm_lax_l[1]*rho_l[3]-0.1936491673103708*ux_surf_lr[1]*rho_c[3]-0.1936491673103708*ux_surf_cl[1]*rho_c[3]+0.3872983346207416*pkpm_lax_l[1]*rho_c[3]+0.1118033988749895*ux_surf_lr[1]*rho_l[2]+0.1118033988749895*ux_surf_cl[1]*rho_l[2]+0.223606797749979*pkpm_lax_l[1]*rho_l[2]+0.1118033988749895*ux_surf_lr[1]*rho_c[2]+0.1118033988749895*ux_surf_cl[1]*rho_c[2]-0.223606797749979*pkpm_lax_l[1]*rho_c[2]+0.4330127018922193*rho_l[1]*pkpm_lax_l[2]+0.4330127018922193*rho_c[1]*pkpm_lax_l[2]+0.25*rho_l[0]*pkpm_lax_l[2]-0.25*rho_c[0]*pkpm_lax_l[2]; 

  flux_rho_r[0] = (-0.5590169943749475*pkpm_lax_r[2]*rho_r[8])+0.5590169943749475*pkpm_lax_r[2]*rho_c[8]+0.4330127018922194*pkpm_lax_r[2]*rho_r[7]+0.4330127018922194*pkpm_lax_r[2]*rho_c[7]+0.2795084971874738*ux_surf_rl[1]*rho_r[6]+0.2795084971874738*ux_surf_cr[1]*rho_r[6]-0.5590169943749476*pkpm_lax_r[1]*rho_r[6]+0.2795084971874738*ux_surf_rl[1]*rho_c[6]+0.2795084971874738*ux_surf_cr[1]*rho_c[6]+0.5590169943749476*pkpm_lax_r[1]*rho_c[6]-0.25*pkpm_lax_r[2]*rho_r[5]+0.25*pkpm_lax_r[2]*rho_c[5]+0.2795084971874737*ux_surf_rl[0]*rho_r[4]+0.2795084971874737*ux_surf_cr[0]*rho_r[4]-0.5590169943749475*pkpm_lax_r[0]*rho_r[4]+0.2795084971874737*ux_surf_rl[0]*rho_c[4]+0.2795084971874737*ux_surf_cr[0]*rho_c[4]+0.5590169943749475*pkpm_lax_r[0]*rho_c[4]-0.2165063509461096*ux_surf_rl[1]*rho_r[3]-0.2165063509461096*ux_surf_cr[1]*rho_r[3]+0.4330127018922193*pkpm_lax_r[1]*rho_r[3]+0.2165063509461096*ux_surf_rl[1]*rho_c[3]+0.2165063509461096*ux_surf_cr[1]*rho_c[3]+0.4330127018922193*pkpm_lax_r[1]*rho_c[3]+0.125*ux_surf_rl[1]*rho_r[2]+0.125*ux_surf_cr[1]*rho_r[2]-0.25*pkpm_lax_r[1]*rho_r[2]+0.125*ux_surf_rl[1]*rho_c[2]+0.125*ux_surf_cr[1]*rho_c[2]+0.25*pkpm_lax_r[1]*rho_c[2]-0.2165063509461096*ux_surf_rl[0]*rho_r[1]-0.2165063509461096*ux_surf_cr[0]*rho_r[1]+0.4330127018922193*pkpm_lax_r[0]*rho_r[1]+0.2165063509461096*ux_surf_rl[0]*rho_c[1]+0.2165063509461096*ux_surf_cr[0]*rho_c[1]+0.4330127018922193*pkpm_lax_r[0]*rho_c[1]+0.125*rho_r[0]*ux_surf_rl[0]+0.125*rho_c[0]*ux_surf_rl[0]+0.125*rho_r[0]*ux_surf_cr[0]+0.125*rho_c[0]*ux_surf_cr[0]-0.25*pkpm_lax_r[0]*rho_r[0]+0.25*pkpm_lax_r[0]*rho_c[0]; 
  flux_rho_r[1] = 0.25*ux_surf_rl[1]*rho_r[8]+0.25*ux_surf_cr[1]*rho_r[8]-0.5*pkpm_lax_r[1]*rho_r[8]+0.25*ux_surf_rl[1]*rho_c[8]+0.25*ux_surf_cr[1]*rho_c[8]+0.5*pkpm_lax_r[1]*rho_c[8]-0.1936491673103709*ux_surf_rl[1]*rho_r[7]-0.1936491673103709*ux_surf_cr[1]*rho_r[7]+0.3872983346207417*pkpm_lax_r[1]*rho_r[7]+0.1936491673103709*ux_surf_rl[1]*rho_c[7]+0.1936491673103709*ux_surf_cr[1]*rho_c[7]+0.3872983346207417*pkpm_lax_r[1]*rho_c[7]-0.5000000000000001*pkpm_lax_r[2]*rho_r[6]+0.2795084971874738*ux_surf_rl[0]*rho_r[6]+0.2795084971874738*ux_surf_cr[0]*rho_r[6]-0.5590169943749476*pkpm_lax_r[0]*rho_r[6]+0.5000000000000001*pkpm_lax_r[2]*rho_c[6]+0.2795084971874738*ux_surf_rl[0]*rho_c[6]+0.2795084971874738*ux_surf_cr[0]*rho_c[6]+0.5590169943749476*pkpm_lax_r[0]*rho_c[6]+0.1118033988749895*ux_surf_rl[1]*rho_r[5]+0.1118033988749895*ux_surf_cr[1]*rho_r[5]-0.223606797749979*pkpm_lax_r[1]*rho_r[5]+0.1118033988749895*ux_surf_rl[1]*rho_c[5]+0.1118033988749895*ux_surf_cr[1]*rho_c[5]+0.223606797749979*pkpm_lax_r[1]*rho_c[5]+0.2795084971874737*ux_surf_rl[1]*rho_r[4]+0.2795084971874737*ux_surf_cr[1]*rho_r[4]-0.5590169943749475*pkpm_lax_r[1]*rho_r[4]+0.2795084971874737*ux_surf_rl[1]*rho_c[4]+0.2795084971874737*ux_surf_cr[1]*rho_c[4]+0.5590169943749475*pkpm_lax_r[1]*rho_c[4]+0.3872983346207416*pkpm_lax_r[2]*rho_r[3]-0.2165063509461096*ux_surf_rl[0]*rho_r[3]-0.2165063509461096*ux_surf_cr[0]*rho_r[3]+0.4330127018922193*pkpm_lax_r[0]*rho_r[3]+0.3872983346207416*pkpm_lax_r[2]*rho_c[3]+0.2165063509461096*ux_surf_rl[0]*rho_c[3]+0.2165063509461096*ux_surf_cr[0]*rho_c[3]+0.4330127018922193*pkpm_lax_r[0]*rho_c[3]-0.223606797749979*pkpm_lax_r[2]*rho_r[2]+0.125*ux_surf_rl[0]*rho_r[2]+0.125*ux_surf_cr[0]*rho_r[2]-0.25*pkpm_lax_r[0]*rho_r[2]+0.223606797749979*pkpm_lax_r[2]*rho_c[2]+0.125*ux_surf_rl[0]*rho_c[2]+0.125*ux_surf_cr[0]*rho_c[2]+0.25*pkpm_lax_r[0]*rho_c[2]-0.2165063509461096*rho_r[1]*ux_surf_rl[1]+0.2165063509461096*rho_c[1]*ux_surf_rl[1]+0.125*rho_r[0]*ux_surf_rl[1]+0.125*rho_c[0]*ux_surf_rl[1]-0.2165063509461096*rho_r[1]*ux_surf_cr[1]+0.2165063509461096*rho_c[1]*ux_surf_cr[1]+0.125*rho_r[0]*ux_surf_cr[1]+0.125*rho_c[0]*ux_surf_cr[1]+0.4330127018922193*pkpm_lax_r[1]*rho_r[1]+0.4330127018922193*pkpm_lax_r[1]*rho_c[1]-0.25*rho_r[0]*pkpm_lax_r[1]+0.25*rho_c[0]*pkpm_lax_r[1]; 
  flux_rho_r[2] = (-0.3571428571428572*pkpm_lax_r[2]*rho_r[8])+0.2795084971874737*ux_surf_rl[0]*rho_r[8]+0.2795084971874737*ux_surf_cr[0]*rho_r[8]-0.5590169943749475*pkpm_lax_r[0]*rho_r[8]+0.3571428571428572*pkpm_lax_r[2]*rho_c[8]+0.2795084971874737*ux_surf_rl[0]*rho_c[8]+0.2795084971874737*ux_surf_cr[0]*rho_c[8]+0.5590169943749475*pkpm_lax_r[0]*rho_c[8]+0.276641667586244*pkpm_lax_r[2]*rho_r[7]-0.2165063509461097*ux_surf_rl[0]*rho_r[7]-0.2165063509461097*ux_surf_cr[0]*rho_r[7]+0.4330127018922194*pkpm_lax_r[0]*rho_r[7]+0.276641667586244*pkpm_lax_r[2]*rho_c[7]+0.2165063509461097*ux_surf_rl[0]*rho_c[7]+0.2165063509461097*ux_surf_cr[0]*rho_c[7]+0.4330127018922194*pkpm_lax_r[0]*rho_c[7]+0.2500000000000001*ux_surf_rl[1]*rho_r[6]+0.2500000000000001*ux_surf_cr[1]*rho_r[6]-0.5000000000000001*pkpm_lax_r[1]*rho_r[6]+0.2500000000000001*ux_surf_rl[1]*rho_c[6]+0.2500000000000001*ux_surf_cr[1]*rho_c[6]+0.5000000000000001*pkpm_lax_r[1]*rho_c[6]-0.159719141249985*pkpm_lax_r[2]*rho_r[5]+0.125*ux_surf_rl[0]*rho_r[5]+0.125*ux_surf_cr[0]*rho_r[5]-0.25*pkpm_lax_r[0]*rho_r[5]+0.159719141249985*pkpm_lax_r[2]*rho_c[5]+0.125*ux_surf_rl[0]*rho_c[5]+0.125*ux_surf_cr[0]*rho_c[5]+0.25*pkpm_lax_r[0]*rho_c[5]-0.5590169943749475*pkpm_lax_r[2]*rho_r[4]+0.5590169943749475*pkpm_lax_r[2]*rho_c[4]-0.1936491673103708*ux_surf_rl[1]*rho_r[3]-0.1936491673103708*ux_surf_cr[1]*rho_r[3]+0.3872983346207416*pkpm_lax_r[1]*rho_r[3]+0.1936491673103708*ux_surf_rl[1]*rho_c[3]+0.1936491673103708*ux_surf_cr[1]*rho_c[3]+0.3872983346207416*pkpm_lax_r[1]*rho_c[3]+0.1118033988749895*ux_surf_rl[1]*rho_r[2]+0.1118033988749895*ux_surf_cr[1]*rho_r[2]-0.223606797749979*pkpm_lax_r[1]*rho_r[2]+0.1118033988749895*ux_surf_rl[1]*rho_c[2]+0.1118033988749895*ux_surf_cr[1]*rho_c[2]+0.223606797749979*pkpm_lax_r[1]*rho_c[2]+0.4330127018922193*rho_r[1]*pkpm_lax_r[2]+0.4330127018922193*rho_c[1]*pkpm_lax_r[2]-0.25*rho_r[0]*pkpm_lax_r[2]+0.25*rho_c[0]*pkpm_lax_r[2]; 

  avg_p_ij_x_l[0] = 0.7905694150420947*Pxx_l[4]+0.7905694150420947*Pxx_c[4]+0.6123724356957944*Pxx_l[1]-0.6123724356957944*Pxx_c[1]+0.3535533905932737*Pxx_l[0]+0.3535533905932737*Pxx_c[0]; 
  avg_p_ij_x_l[1] = 0.7905694150420948*Pxx_l[6]+0.7905694150420948*Pxx_c[6]+0.6123724356957944*Pxx_l[3]-0.6123724356957944*Pxx_c[3]+0.3535533905932737*Pxx_l[2]+0.3535533905932737*Pxx_c[2]; 
  avg_p_ij_x_l[2] = 0.7905694150420947*Pxx_l[8]+0.7905694150420947*Pxx_c[8]+0.6123724356957944*Pxx_l[7]-0.6123724356957944*Pxx_c[7]+0.3535533905932737*Pxx_l[5]+0.3535533905932737*Pxx_c[5]; 

  avg_p_ij_x_r[0] = 0.7905694150420947*Pxx_r[4]+0.7905694150420947*Pxx_c[4]-0.6123724356957944*Pxx_r[1]+0.6123724356957944*Pxx_c[1]+0.3535533905932737*Pxx_r[0]+0.3535533905932737*Pxx_c[0]; 
  avg_p_ij_x_r[1] = 0.7905694150420948*Pxx_r[6]+0.7905694150420948*Pxx_c[6]-0.6123724356957944*Pxx_r[3]+0.6123724356957944*Pxx_c[3]+0.3535533905932737*Pxx_r[2]+0.3535533905932737*Pxx_c[2]; 
  avg_p_ij_x_r[2] = 0.7905694150420947*Pxx_r[8]+0.7905694150420947*Pxx_c[8]-0.6123724356957944*Pxx_r[7]+0.6123724356957944*Pxx_c[7]+0.3535533905932737*Pxx_r[5]+0.3535533905932737*Pxx_c[5]; 

  avg_p_ij_y_l[0] = 0.7905694150420947*Pxy_l[4]+0.7905694150420947*Pxy_c[4]+0.6123724356957944*Pxy_l[1]-0.6123724356957944*Pxy_c[1]+0.3535533905932737*Pxy_l[0]+0.3535533905932737*Pxy_c[0]; 
  avg_p_ij_y_l[1] = 0.7905694150420948*Pxy_l[6]+0.7905694150420948*Pxy_c[6]+0.6123724356957944*Pxy_l[3]-0.6123724356957944*Pxy_c[3]+0.3535533905932737*Pxy_l[2]+0.3535533905932737*Pxy_c[2]; 
  avg_p_ij_y_l[2] = 0.7905694150420947*Pxy_l[8]+0.7905694150420947*Pxy_c[8]+0.6123724356957944*Pxy_l[7]-0.6123724356957944*Pxy_c[7]+0.3535533905932737*Pxy_l[5]+0.3535533905932737*Pxy_c[5]; 

  avg_p_ij_y_r[0] = 0.7905694150420947*Pxy_r[4]+0.7905694150420947*Pxy_c[4]-0.6123724356957944*Pxy_r[1]+0.6123724356957944*Pxy_c[1]+0.3535533905932737*Pxy_r[0]+0.3535533905932737*Pxy_c[0]; 
  avg_p_ij_y_r[1] = 0.7905694150420948*Pxy_r[6]+0.7905694150420948*Pxy_c[6]-0.6123724356957944*Pxy_r[3]+0.6123724356957944*Pxy_c[3]+0.3535533905932737*Pxy_r[2]+0.3535533905932737*Pxy_c[2]; 
  avg_p_ij_y_r[2] = 0.7905694150420947*Pxy_r[8]+0.7905694150420947*Pxy_c[8]-0.6123724356957944*Pxy_r[7]+0.6123724356957944*Pxy_c[7]+0.3535533905932737*Pxy_r[5]+0.3535533905932737*Pxy_c[5]; 

  avg_p_ij_z_l[0] = 0.7905694150420947*Pxz_l[4]+0.7905694150420947*Pxz_c[4]+0.6123724356957944*Pxz_l[1]-0.6123724356957944*Pxz_c[1]+0.3535533905932737*Pxz_l[0]+0.3535533905932737*Pxz_c[0]; 
  avg_p_ij_z_l[1] = 0.7905694150420948*Pxz_l[6]+0.7905694150420948*Pxz_c[6]+0.6123724356957944*Pxz_l[3]-0.6123724356957944*Pxz_c[3]+0.3535533905932737*Pxz_l[2]+0.3535533905932737*Pxz_c[2]; 
  avg_p_ij_z_l[2] = 0.7905694150420947*Pxz_l[8]+0.7905694150420947*Pxz_c[8]+0.6123724356957944*Pxz_l[7]-0.6123724356957944*Pxz_c[7]+0.3535533905932737*Pxz_l[5]+0.3535533905932737*Pxz_c[5]; 

  avg_p_ij_z_r[0] = 0.7905694150420947*Pxz_r[4]+0.7905694150420947*Pxz_c[4]-0.6123724356957944*Pxz_r[1]+0.6123724356957944*Pxz_c[1]+0.3535533905932737*Pxz_r[0]+0.3535533905932737*Pxz_c[0]; 
  avg_p_ij_z_r[1] = 0.7905694150420948*Pxz_r[6]+0.7905694150420948*Pxz_c[6]-0.6123724356957944*Pxz_r[3]+0.6123724356957944*Pxz_c[3]+0.3535533905932737*Pxz_r[2]+0.3535533905932737*Pxz_c[2]; 
  avg_p_ij_z_r[2] = 0.7905694150420947*Pxz_r[8]+0.7905694150420947*Pxz_c[8]-0.6123724356957944*Pxz_r[7]+0.6123724356957944*Pxz_c[7]+0.3535533905932737*Pxz_r[5]+0.3535533905932737*Pxz_c[5]; 

  double amdq_rhoux_l[3] = {0.0}; 
  double apdq_rhoux_l[3] = {0.0}; 
  double amdq_rhouy_l[3] = {0.0}; 
  double apdq_rhouy_l[3] = {0.0}; 
  double amdq_rhouz_l[3] = {0.0}; 
  double apdq_rhouz_l[3] = {0.0}; 

  double amdq_rhoux_r[3] = {0.0}; 
  double apdq_rhoux_r[3] = {0.0}; 
  double amdq_rhouy_r[3] = {0.0}; 
  double apdq_rhouy_r[3] = {0.0}; 
  double amdq_rhouz_r[3] = {0.0}; 
  double apdq_rhouz_r[3] = {0.0}; 

  double amdq_rhoux_quad_l[3] = {0.0}; 
  double apdq_rhoux_quad_l[3] = {0.0}; 
  double amdq_rhouy_quad_l[3] = {0.0}; 
  double apdq_rhouy_quad_l[3] = {0.0}; 
  double amdq_rhouz_quad_l[3] = {0.0}; 
  double apdq_rhouz_quad_l[3] = {0.0}; 

  double amdq_rhoux_quad_r[3] = {0.0}; 
  double apdq_rhoux_quad_r[3] = {0.0}; 
  double amdq_rhouy_quad_r[3] = {0.0}; 
  double apdq_rhouy_quad_r[3] = {0.0}; 
  double amdq_rhouz_quad_r[3] = {0.0}; 
  double apdq_rhouz_quad_r[3] = {0.0}; 

  double q_lr[10] = {0.0}; 
  double q_cl[10] = {0.0}; 
  double q_cr[10] = {0.0}; 
  double q_rl[10] = {0.0}; 
  double q_lr_local[10] = {0.0}; 
  double q_cl_local[10] = {0.0}; 
  double q_cr_local[10] = {0.0}; 
  double q_rl_local[10] = {0.0}; 
  double delta_l[10] = {0.0}; 
  double delta_r[10] = {0.0}; 
  double my_max_speed_l = 0.0; 
  double my_max_speed_r = 0.0; 
  double lenr_l = 0.0; 
  double lenr_r = 0.0; 
  double waves_l[50] = {0.0}; 
  double waves_r[50] = {0.0}; 
  double speeds_l[5] = {0.0}; 
  double speeds_r[5] = {0.0}; 
  double amdq_l_local[10] = {0.0}; 
  double apdq_l_local[10] = {0.0}; 
  double amdq_r_local[10] = {0.0}; 
  double apdq_r_local[10] = {0.0}; 
  double amdq_l[10] = {0.0}; 
  double apdq_l[10] = {0.0}; 
  double amdq_r[10] = {0.0}; 
  double apdq_r[10] = {0.0}; 

  q_lr[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rho_l); 
  q_lr[1] = q_lr[0]*(0.7071067811865475*ux_surf_lr[0]-0.9486832980505137*ux_surf_lr[1]); 
  q_lr[2] = q_lr[0]*(0.7071067811865475*uy_surf_lr[0]-0.9486832980505137*uy_surf_lr[1]); 
  q_lr[3] = q_lr[0]*(0.7071067811865475*uz_surf_lr[0]-0.9486832980505137*uz_surf_lr[1]); 
  q_lr[4] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rho_c); 
  q_cl[1] = q_cl[0]*(0.7071067811865475*ux_surf_cl[0]-0.9486832980505137*ux_surf_cl[1]); 
  q_cl[2] = q_cl[0]*(0.7071067811865475*uy_surf_cl[0]-0.9486832980505137*uy_surf_cl[1]); 
  q_cl[3] = q_cl[0]*(0.7071067811865475*uz_surf_cl[0]-0.9486832980505137*uz_surf_cl[1]); 
  q_cl[4] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(rho_c); 
  q_cr[1] = q_cr[0]*(0.7071067811865475*ux_surf_cr[0]-0.9486832980505137*ux_surf_cr[1]); 
  q_cr[2] = q_cr[0]*(0.7071067811865475*uy_surf_cr[0]-0.9486832980505137*uy_surf_cr[1]); 
  q_cr[3] = q_cr[0]*(0.7071067811865475*uz_surf_cr[0]-0.9486832980505137*uz_surf_cr[1]); 
  q_cr[4] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = tensor_2x_p2_surfx1_eval_quad_node_0_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(rho_r); 
  q_rl[1] = q_rl[0]*(0.7071067811865475*ux_surf_rl[0]-0.9486832980505137*ux_surf_rl[1]); 
  q_rl[2] = q_rl[0]*(0.7071067811865475*uy_surf_rl[0]-0.9486832980505137*uy_surf_rl[1]); 
  q_rl[3] = q_rl[0]*(0.7071067811865475*uz_surf_rl[0]-0.9486832980505137*uz_surf_rl[1]); 
  q_rl[4] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = tensor_2x_p2_surfx1_eval_quad_node_0_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[0]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[0]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[0] = amdq_l[1]; 
  apdq_rhoux_quad_l[0] = apdq_l[1]; 
  amdq_rhouy_quad_l[0] = amdq_l[2]; 
  apdq_rhouy_quad_l[0] = apdq_l[2]; 
  amdq_rhouz_quad_l[0] = amdq_l[3]; 
  apdq_rhouz_quad_l[0] = apdq_l[3]; 

  amdq_rhoux_quad_r[0] = amdq_r[1]; 
  apdq_rhoux_quad_r[0] = apdq_r[1]; 
  amdq_rhouy_quad_r[0] = amdq_r[2]; 
  apdq_rhouy_quad_r[0] = apdq_r[2]; 
  amdq_rhouz_quad_r[0] = amdq_r[3]; 
  apdq_rhouz_quad_r[0] = apdq_r[3]; 

  q_lr[0] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rho_l); 
  q_lr[1] = q_lr[0]*(0.7071067811865475*ux_surf_lr[0]); 
  q_lr[2] = q_lr[0]*(0.7071067811865475*uy_surf_lr[0]); 
  q_lr[3] = q_lr[0]*(0.7071067811865475*uz_surf_lr[0]); 
  q_lr[4] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rho_c); 
  q_cl[1] = q_cl[0]*(0.7071067811865475*ux_surf_cl[0]); 
  q_cl[2] = q_cl[0]*(0.7071067811865475*uy_surf_cl[0]); 
  q_cl[3] = q_cl[0]*(0.7071067811865475*uz_surf_cl[0]); 
  q_cl[4] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = tensor_2x_p2_surfx1_eval_quad_node_1_r(rho_c); 
  q_cr[1] = q_cr[0]*(0.7071067811865475*ux_surf_cr[0]); 
  q_cr[2] = q_cr[0]*(0.7071067811865475*uy_surf_cr[0]); 
  q_cr[3] = q_cr[0]*(0.7071067811865475*uz_surf_cr[0]); 
  q_cr[4] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = tensor_2x_p2_surfx1_eval_quad_node_1_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = tensor_2x_p2_surfx1_eval_quad_node_1_l(rho_r); 
  q_rl[1] = q_rl[0]*(0.7071067811865475*ux_surf_rl[0]); 
  q_rl[2] = q_rl[0]*(0.7071067811865475*uy_surf_rl[0]); 
  q_rl[3] = q_rl[0]*(0.7071067811865475*uz_surf_rl[0]); 
  q_rl[4] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = tensor_2x_p2_surfx1_eval_quad_node_1_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[0]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[0]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[1] = amdq_l[1]; 
  apdq_rhoux_quad_l[1] = apdq_l[1]; 
  amdq_rhouy_quad_l[1] = amdq_l[2]; 
  apdq_rhouy_quad_l[1] = apdq_l[2]; 
  amdq_rhouz_quad_l[1] = amdq_l[3]; 
  apdq_rhouz_quad_l[1] = apdq_l[3]; 

  amdq_rhoux_quad_r[1] = amdq_r[1]; 
  apdq_rhoux_quad_r[1] = apdq_r[1]; 
  amdq_rhouy_quad_r[1] = amdq_r[2]; 
  apdq_rhouy_quad_r[1] = apdq_r[2]; 
  amdq_rhouz_quad_r[1] = amdq_r[3]; 
  apdq_rhouz_quad_r[1] = apdq_r[3]; 

  q_lr[0] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rho_l); 
  q_lr[1] = q_lr[0]*(0.9486832980505137*ux_surf_lr[1]+0.7071067811865475*ux_surf_lr[0]); 
  q_lr[2] = q_lr[0]*(0.9486832980505137*uy_surf_lr[1]+0.7071067811865475*uy_surf_lr[0]); 
  q_lr[3] = q_lr[0]*(0.9486832980505137*uz_surf_lr[1]+0.7071067811865475*uz_surf_lr[0]); 
  q_lr[4] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pxx_l) + q_lr[1]*q_lr[1]/q_lr[0]; 
  q_lr[5] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pxy_l) + q_lr[1]*q_lr[2]/q_lr[0]; 
  q_lr[6] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pxz_l) + q_lr[1]*q_lr[3]/q_lr[0]; 
  q_lr[7] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pyy_l) + q_lr[2]*q_lr[2]/q_lr[0]; 
  q_lr[8] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pyz_l) + q_lr[2]*q_lr[3]/q_lr[0]; 
  q_lr[9] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pzz_l) + q_lr[3]*q_lr[3]/q_lr[0]; 
  q_cl[0] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rho_c); 
  q_cl[1] = q_cl[0]*(0.9486832980505137*ux_surf_cl[1]+0.7071067811865475*ux_surf_cl[0]); 
  q_cl[2] = q_cl[0]*(0.9486832980505137*uy_surf_cl[1]+0.7071067811865475*uy_surf_cl[0]); 
  q_cl[3] = q_cl[0]*(0.9486832980505137*uz_surf_cl[1]+0.7071067811865475*uz_surf_cl[0]); 
  q_cl[4] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pxx_c) + q_cl[1]*q_cl[1]/q_cl[0]; 
  q_cl[5] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pxy_c) + q_cl[1]*q_cl[2]/q_cl[0]; 
  q_cl[6] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pxz_c) + q_cl[1]*q_cl[3]/q_cl[0]; 
  q_cl[7] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pyy_c) + q_cl[2]*q_cl[2]/q_cl[0]; 
  q_cl[8] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pyz_c) + q_cl[2]*q_cl[3]/q_cl[0]; 
  q_cl[9] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pzz_c) + q_cl[3]*q_cl[3]/q_cl[0]; 
  q_cr[0] = tensor_2x_p2_surfx1_eval_quad_node_2_r(rho_c); 
  q_cr[1] = q_cr[0]*(0.9486832980505137*ux_surf_cr[1]+0.7071067811865475*ux_surf_cr[0]); 
  q_cr[2] = q_cr[0]*(0.9486832980505137*uy_surf_cr[1]+0.7071067811865475*uy_surf_cr[0]); 
  q_cr[3] = q_cr[0]*(0.9486832980505137*uz_surf_cr[1]+0.7071067811865475*uz_surf_cr[0]); 
  q_cr[4] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pxx_c) + q_cr[1]*q_cr[1]/q_cr[0]; 
  q_cr[5] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pxy_c) + q_cr[1]*q_cr[2]/q_cr[0]; 
  q_cr[6] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pxz_c) + q_cr[1]*q_cr[3]/q_cr[0]; 
  q_cr[7] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pyy_c) + q_cr[2]*q_cr[2]/q_cr[0]; 
  q_cr[8] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pyz_c) + q_cr[2]*q_cr[3]/q_cr[0]; 
  q_cr[9] = tensor_2x_p2_surfx1_eval_quad_node_2_r(Pzz_c) + q_cr[3]*q_cr[3]/q_cr[0]; 
  q_rl[0] = tensor_2x_p2_surfx1_eval_quad_node_2_l(rho_r); 
  q_rl[1] = q_rl[0]*(0.9486832980505137*ux_surf_rl[1]+0.7071067811865475*ux_surf_rl[0]); 
  q_rl[2] = q_rl[0]*(0.9486832980505137*uy_surf_rl[1]+0.7071067811865475*uy_surf_rl[0]); 
  q_rl[3] = q_rl[0]*(0.9486832980505137*uz_surf_rl[1]+0.7071067811865475*uz_surf_rl[0]); 
  q_rl[4] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pxx_r) + q_rl[1]*q_rl[1]/q_rl[0]; 
  q_rl[5] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pxy_r) + q_rl[1]*q_rl[2]/q_rl[0]; 
  q_rl[6] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pxz_r) + q_rl[1]*q_rl[3]/q_rl[0]; 
  q_rl[7] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pyy_r) + q_rl[2]*q_rl[2]/q_rl[0]; 
  q_rl[8] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pyz_r) + q_rl[2]*q_rl[3]/q_rl[0]; 
  q_rl[9] = tensor_2x_p2_surfx1_eval_quad_node_2_l(Pzz_r) + q_rl[3]*q_rl[3]/q_rl[0]; 

  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_lr, q_lr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], q_cl, q_cl_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_cr, q_cr_local); 
  gkyl_wv_eqn_rotate_to_local(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], q_rl, q_rl_local); 

  delta_l[0] = q_cl_local[0] - q_lr_local[0]; 
  delta_l[1] = q_cl_local[1] - q_lr_local[1]; 
  delta_l[2] = q_cl_local[2] - q_lr_local[2]; 
  delta_l[3] = q_cl_local[3] - q_lr_local[3]; 
  delta_l[4] = q_cl_local[4] - q_lr_local[4]; 
  delta_l[5] = q_cl_local[5] - q_lr_local[5]; 
  delta_l[6] = q_cl_local[6] - q_lr_local[6]; 
  delta_l[7] = q_cl_local[7] - q_lr_local[7]; 
  delta_l[8] = q_cl_local[8] - q_lr_local[8]; 
  delta_l[9] = q_cl_local[9] - q_lr_local[9]; 
  delta_r[0] = q_rl_local[0] - q_cr_local[0]; 
  delta_r[1] = q_rl_local[1] - q_cr_local[1]; 
  delta_r[2] = q_rl_local[2] - q_cr_local[2]; 
  delta_r[3] = q_rl_local[3] - q_cr_local[3]; 
  delta_r[4] = q_rl_local[4] - q_cr_local[4]; 
  delta_r[5] = q_rl_local[5] - q_cr_local[5]; 
  delta_r[6] = q_rl_local[6] - q_cr_local[6]; 
  delta_r[7] = q_rl_local[7] - q_cr_local[7]; 
  delta_r[8] = q_rl_local[8] - q_cr_local[8]; 
  delta_r[9] = q_rl_local[9] - q_cr_local[9]; 

  my_max_speed_l = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_lr_local, q_cl_local, waves_l, speeds_l); 
  my_max_speed_r = gkyl_wv_eqn_waves(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_cr_local, q_rl_local, waves_r, speeds_r); 
  lenr_l = geom_l->lenr[0]; 
  speeds_l[0] *= lenr_l; 
  speeds_l[1] *= lenr_l; 
  speeds_l[2] *= lenr_l; 
  speeds_l[3] *= lenr_l; 
  speeds_l[4] *= lenr_l; 
  lenr_r = geom_r->lenr[0]; 
  speeds_r[0] *= lenr_r; 
  speeds_r[1] *= lenr_r; 
  speeds_r[2] *= lenr_r; 
  speeds_r[3] *= lenr_r; 
  speeds_r[4] *= lenr_r; 

  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_lr_local, q_cl_local, waves_l, speeds_l, amdq_l_local, apdq_l_local); 
  gkyl_wv_eqn_qfluct(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, q_cr_local, q_rl_local, waves_r, speeds_r, amdq_r_local, apdq_r_local); 

  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], amdq_l_local, amdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_l->tau1[0], geom_l->tau2[0], geom_l->norm[0], apdq_l_local, apdq_l); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], amdq_r_local, amdq_r); 
  gkyl_wv_eqn_rotate_to_global(wv_eqn, geom_r->tau1[0], geom_r->tau2[0], geom_r->norm[0], apdq_r_local, apdq_r); 

  amdq_rhoux_quad_l[2] = amdq_l[1]; 
  apdq_rhoux_quad_l[2] = apdq_l[1]; 
  amdq_rhouy_quad_l[2] = amdq_l[2]; 
  apdq_rhouy_quad_l[2] = apdq_l[2]; 
  amdq_rhouz_quad_l[2] = amdq_l[3]; 
  apdq_rhouz_quad_l[2] = apdq_l[3]; 

  amdq_rhoux_quad_r[2] = amdq_r[1]; 
  apdq_rhoux_quad_r[2] = apdq_r[1]; 
  amdq_rhouy_quad_r[2] = amdq_r[2]; 
  apdq_rhouy_quad_r[2] = apdq_r[2]; 
  amdq_rhouz_quad_r[2] = amdq_r[3]; 
  apdq_rhouz_quad_r[2] = apdq_r[3]; 

  tensor_2x_p2_upwind_quad_to_modal(amdq_rhoux_quad_l, amdq_rhoux_l); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouy_quad_l, amdq_rhouy_l); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouz_quad_l, amdq_rhouz_l); 

  tensor_2x_p2_upwind_quad_to_modal(apdq_rhoux_quad_l, apdq_rhoux_l); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouy_quad_l, apdq_rhouy_l); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouz_quad_l, apdq_rhouz_l); 

  tensor_2x_p2_upwind_quad_to_modal(amdq_rhoux_quad_r, amdq_rhoux_r); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouy_quad_r, amdq_rhouy_r); 
  tensor_2x_p2_upwind_quad_to_modal(amdq_rhouz_quad_r, amdq_rhouz_r); 

  tensor_2x_p2_upwind_quad_to_modal(apdq_rhoux_quad_r, apdq_rhoux_r); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouy_quad_r, apdq_rhouy_r); 
  tensor_2x_p2_upwind_quad_to_modal(apdq_rhouz_quad_r, apdq_rhouz_r); 

  outrhou0[0] += ((-0.25*flux_rho_r[1]*ux_surf_rl[1])+0.25*flux_rho_l[1]*ux_surf_lr[1]-0.25*flux_rho_r[1]*ux_surf_cr[1]+0.25*flux_rho_l[1]*ux_surf_cl[1]-0.25*flux_rho_r[0]*ux_surf_rl[0]+0.25*flux_rho_l[0]*ux_surf_lr[0]-0.25*flux_rho_r[0]*ux_surf_cr[0]+0.25*flux_rho_l[0]*ux_surf_cl[0]-0.7071067811865475*avg_p_ij_x_r[0]+0.7071067811865475*avg_p_ij_x_l[0]+0.3535533905932737*apdq_rhoux_r[0]-0.3535533905932737*(apdq_rhoux_l[0]+amdq_rhoux_r[0])+0.3535533905932737*amdq_rhoux_l[0])*dx1; 
  outrhou0[1] += ((-0.4330127018922193*(flux_rho_r[1]*ux_surf_rl[1]+flux_rho_l[1]*ux_surf_lr[1]+flux_rho_r[1]*ux_surf_cr[1]+flux_rho_l[1]*ux_surf_cl[1]+flux_rho_r[0]*ux_surf_rl[0]+flux_rho_l[0]*ux_surf_lr[0]+flux_rho_r[0]*ux_surf_cr[0]+flux_rho_l[0]*ux_surf_cl[0]))-1.224744871391589*(avg_p_ij_x_r[0]+avg_p_ij_x_l[0])+0.6123724356957944*(apdq_rhoux_r[0]+apdq_rhoux_l[0])-0.6123724356957944*(amdq_rhoux_r[0]+amdq_rhoux_l[0]))*dx1; 
  outrhou0[2] += ((-0.223606797749979*(ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[2])+0.223606797749979*(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[2]-0.25*flux_rho_r[0]*ux_surf_rl[1]+0.25*flux_rho_l[0]*ux_surf_lr[1]-0.25*flux_rho_r[0]*ux_surf_cr[1]+0.25*flux_rho_l[0]*ux_surf_cl[1]-0.25*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+0.25*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_x_r[1]+0.7071067811865475*avg_p_ij_x_l[1]+0.3535533905932737*apdq_rhoux_r[1]-0.3535533905932737*(apdq_rhoux_l[1]+amdq_rhoux_r[1])+0.3535533905932737*amdq_rhoux_l[1])*dx1; 
  outrhou0[3] += ((-0.3872983346207416*((ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[2]+(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[2]))-0.4330127018922193*(flux_rho_r[0]*ux_surf_rl[1]+flux_rho_l[0]*ux_surf_lr[1]+flux_rho_r[0]*ux_surf_cr[1]+flux_rho_l[0]*ux_surf_cl[1]+(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_x_r[1]+avg_p_ij_x_l[1])+0.6123724356957944*(apdq_rhoux_r[1]+apdq_rhoux_l[1])-0.6123724356957944*(amdq_rhoux_r[1]+amdq_rhoux_l[1]))*dx1; 
  outrhou0[4] += ((-0.5590169943749475*flux_rho_r[1]*ux_surf_rl[1])+0.5590169943749475*flux_rho_l[1]*ux_surf_lr[1]-0.5590169943749475*flux_rho_r[1]*ux_surf_cr[1]+0.5590169943749475*flux_rho_l[1]*ux_surf_cl[1]-0.5590169943749475*flux_rho_r[0]*ux_surf_rl[0]+0.5590169943749475*flux_rho_l[0]*ux_surf_lr[0]-0.5590169943749475*flux_rho_r[0]*ux_surf_cr[0]+0.5590169943749475*flux_rho_l[0]*ux_surf_cl[0]-1.58113883008419*avg_p_ij_x_r[0]+1.58113883008419*avg_p_ij_x_l[0]+0.7905694150420947*apdq_rhoux_r[0]-0.7905694150420947*(apdq_rhoux_l[0]+amdq_rhoux_r[0])+0.7905694150420947*amdq_rhoux_l[0])*dx1; 
  outrhou0[5] += ((-0.25*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[2])+0.25*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[2]-0.7071067811865475*avg_p_ij_x_r[2]+0.7071067811865475*avg_p_ij_x_l[2]+0.3535533905932737*apdq_rhoux_r[2]-0.3535533905932737*(apdq_rhoux_l[2]+amdq_rhoux_r[2])+0.3535533905932737*amdq_rhoux_l[2]-0.223606797749979*flux_rho_r[1]*ux_surf_rl[1]+0.223606797749979*flux_rho_l[1]*ux_surf_lr[1]-0.223606797749979*flux_rho_r[1]*ux_surf_cr[1]+0.223606797749979*flux_rho_l[1]*ux_surf_cl[1])*dx1; 
  outrhou0[6] += ((-0.5000000000000001*(ux_surf_rl[1]+ux_surf_cr[1])*flux_rho_r[2])+0.5000000000000001*(ux_surf_lr[1]+ux_surf_cl[1])*flux_rho_l[2]-0.5590169943749476*flux_rho_r[0]*ux_surf_rl[1]+0.5590169943749476*flux_rho_l[0]*ux_surf_lr[1]-0.5590169943749476*flux_rho_r[0]*ux_surf_cr[1]+0.5590169943749476*flux_rho_l[0]*ux_surf_cl[1]-0.5590169943749476*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[1]+0.5590169943749476*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[1]-1.58113883008419*avg_p_ij_x_r[1]+1.58113883008419*avg_p_ij_x_l[1]+0.7905694150420948*apdq_rhoux_r[1]-0.7905694150420948*(apdq_rhoux_l[1]+amdq_rhoux_r[1])+0.7905694150420948*amdq_rhoux_l[1])*dx1; 
  outrhou0[7] += ((-0.4330127018922194*((ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[2]+(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[2]))-1.224744871391589*(avg_p_ij_x_r[2]+avg_p_ij_x_l[2])+0.6123724356957944*(apdq_rhoux_r[2]+apdq_rhoux_l[2])-0.6123724356957944*(amdq_rhoux_r[2]+amdq_rhoux_l[2])-0.3872983346207417*(flux_rho_r[1]*ux_surf_rl[1]+flux_rho_l[1]*ux_surf_lr[1]+flux_rho_r[1]*ux_surf_cr[1]+flux_rho_l[1]*ux_surf_cl[1]))*dx1; 
  outrhou0[8] += ((-0.5590169943749475*(ux_surf_rl[0]+ux_surf_cr[0])*flux_rho_r[2])+0.5590169943749475*(ux_surf_lr[0]+ux_surf_cl[0])*flux_rho_l[2]-1.58113883008419*avg_p_ij_x_r[2]+1.58113883008419*avg_p_ij_x_l[2]+0.7905694150420947*apdq_rhoux_r[2]-0.7905694150420947*(apdq_rhoux_l[2]+amdq_rhoux_r[2])+0.7905694150420947*amdq_rhoux_l[2]-0.5*flux_rho_r[1]*ux_surf_rl[1]+0.5*flux_rho_l[1]*ux_surf_lr[1]-0.5*flux_rho_r[1]*ux_surf_cr[1]+0.5*flux_rho_l[1]*ux_surf_cl[1])*dx1; 

  outrhou1[0] += ((-0.25*flux_rho_r[1]*uy_surf_rl[1])+0.25*flux_rho_l[1]*uy_surf_lr[1]-0.25*flux_rho_r[1]*uy_surf_cr[1]+0.25*flux_rho_l[1]*uy_surf_cl[1]-0.25*flux_rho_r[0]*uy_surf_rl[0]+0.25*flux_rho_l[0]*uy_surf_lr[0]-0.25*flux_rho_r[0]*uy_surf_cr[0]+0.25*flux_rho_l[0]*uy_surf_cl[0]-0.7071067811865475*avg_p_ij_y_r[0]+0.7071067811865475*avg_p_ij_y_l[0]+0.3535533905932737*apdq_rhouy_r[0]-0.3535533905932737*(apdq_rhouy_l[0]+amdq_rhouy_r[0])+0.3535533905932737*amdq_rhouy_l[0])*dx1; 
  outrhou1[1] += ((-0.4330127018922193*(flux_rho_r[1]*uy_surf_rl[1]+flux_rho_l[1]*uy_surf_lr[1]+flux_rho_r[1]*uy_surf_cr[1]+flux_rho_l[1]*uy_surf_cl[1]+flux_rho_r[0]*uy_surf_rl[0]+flux_rho_l[0]*uy_surf_lr[0]+flux_rho_r[0]*uy_surf_cr[0]+flux_rho_l[0]*uy_surf_cl[0]))-1.224744871391589*(avg_p_ij_y_r[0]+avg_p_ij_y_l[0])+0.6123724356957944*(apdq_rhouy_r[0]+apdq_rhouy_l[0])-0.6123724356957944*(amdq_rhouy_r[0]+amdq_rhouy_l[0]))*dx1; 
  outrhou1[2] += ((-0.223606797749979*(uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[2])+0.223606797749979*(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[2]-0.25*flux_rho_r[0]*uy_surf_rl[1]+0.25*flux_rho_l[0]*uy_surf_lr[1]-0.25*flux_rho_r[0]*uy_surf_cr[1]+0.25*flux_rho_l[0]*uy_surf_cl[1]-0.25*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+0.25*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_y_r[1]+0.7071067811865475*avg_p_ij_y_l[1]+0.3535533905932737*apdq_rhouy_r[1]-0.3535533905932737*(apdq_rhouy_l[1]+amdq_rhouy_r[1])+0.3535533905932737*amdq_rhouy_l[1])*dx1; 
  outrhou1[3] += ((-0.3872983346207416*((uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[2]+(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[2]))-0.4330127018922193*(flux_rho_r[0]*uy_surf_rl[1]+flux_rho_l[0]*uy_surf_lr[1]+flux_rho_r[0]*uy_surf_cr[1]+flux_rho_l[0]*uy_surf_cl[1]+(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_y_r[1]+avg_p_ij_y_l[1])+0.6123724356957944*(apdq_rhouy_r[1]+apdq_rhouy_l[1])-0.6123724356957944*(amdq_rhouy_r[1]+amdq_rhouy_l[1]))*dx1; 
  outrhou1[4] += ((-0.5590169943749475*flux_rho_r[1]*uy_surf_rl[1])+0.5590169943749475*flux_rho_l[1]*uy_surf_lr[1]-0.5590169943749475*flux_rho_r[1]*uy_surf_cr[1]+0.5590169943749475*flux_rho_l[1]*uy_surf_cl[1]-0.5590169943749475*flux_rho_r[0]*uy_surf_rl[0]+0.5590169943749475*flux_rho_l[0]*uy_surf_lr[0]-0.5590169943749475*flux_rho_r[0]*uy_surf_cr[0]+0.5590169943749475*flux_rho_l[0]*uy_surf_cl[0]-1.58113883008419*avg_p_ij_y_r[0]+1.58113883008419*avg_p_ij_y_l[0]+0.7905694150420947*apdq_rhouy_r[0]-0.7905694150420947*(apdq_rhouy_l[0]+amdq_rhouy_r[0])+0.7905694150420947*amdq_rhouy_l[0])*dx1; 
  outrhou1[5] += ((-0.25*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[2])+0.25*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[2]-0.7071067811865475*avg_p_ij_y_r[2]+0.7071067811865475*avg_p_ij_y_l[2]+0.3535533905932737*apdq_rhouy_r[2]-0.3535533905932737*(apdq_rhouy_l[2]+amdq_rhouy_r[2])+0.3535533905932737*amdq_rhouy_l[2]-0.223606797749979*flux_rho_r[1]*uy_surf_rl[1]+0.223606797749979*flux_rho_l[1]*uy_surf_lr[1]-0.223606797749979*flux_rho_r[1]*uy_surf_cr[1]+0.223606797749979*flux_rho_l[1]*uy_surf_cl[1])*dx1; 
  outrhou1[6] += ((-0.5000000000000001*(uy_surf_rl[1]+uy_surf_cr[1])*flux_rho_r[2])+0.5000000000000001*(uy_surf_lr[1]+uy_surf_cl[1])*flux_rho_l[2]-0.5590169943749476*flux_rho_r[0]*uy_surf_rl[1]+0.5590169943749476*flux_rho_l[0]*uy_surf_lr[1]-0.5590169943749476*flux_rho_r[0]*uy_surf_cr[1]+0.5590169943749476*flux_rho_l[0]*uy_surf_cl[1]-0.5590169943749476*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[1]+0.5590169943749476*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[1]-1.58113883008419*avg_p_ij_y_r[1]+1.58113883008419*avg_p_ij_y_l[1]+0.7905694150420948*apdq_rhouy_r[1]-0.7905694150420948*(apdq_rhouy_l[1]+amdq_rhouy_r[1])+0.7905694150420948*amdq_rhouy_l[1])*dx1; 
  outrhou1[7] += ((-0.4330127018922194*((uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[2]+(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[2]))-1.224744871391589*(avg_p_ij_y_r[2]+avg_p_ij_y_l[2])+0.6123724356957944*(apdq_rhouy_r[2]+apdq_rhouy_l[2])-0.6123724356957944*(amdq_rhouy_r[2]+amdq_rhouy_l[2])-0.3872983346207417*(flux_rho_r[1]*uy_surf_rl[1]+flux_rho_l[1]*uy_surf_lr[1]+flux_rho_r[1]*uy_surf_cr[1]+flux_rho_l[1]*uy_surf_cl[1]))*dx1; 
  outrhou1[8] += ((-0.5590169943749475*(uy_surf_rl[0]+uy_surf_cr[0])*flux_rho_r[2])+0.5590169943749475*(uy_surf_lr[0]+uy_surf_cl[0])*flux_rho_l[2]-1.58113883008419*avg_p_ij_y_r[2]+1.58113883008419*avg_p_ij_y_l[2]+0.7905694150420947*apdq_rhouy_r[2]-0.7905694150420947*(apdq_rhouy_l[2]+amdq_rhouy_r[2])+0.7905694150420947*amdq_rhouy_l[2]-0.5*flux_rho_r[1]*uy_surf_rl[1]+0.5*flux_rho_l[1]*uy_surf_lr[1]-0.5*flux_rho_r[1]*uy_surf_cr[1]+0.5*flux_rho_l[1]*uy_surf_cl[1])*dx1; 

  outrhou2[0] += ((-0.25*flux_rho_r[1]*uz_surf_rl[1])+0.25*flux_rho_l[1]*uz_surf_lr[1]-0.25*flux_rho_r[1]*uz_surf_cr[1]+0.25*flux_rho_l[1]*uz_surf_cl[1]-0.25*flux_rho_r[0]*uz_surf_rl[0]+0.25*flux_rho_l[0]*uz_surf_lr[0]-0.25*flux_rho_r[0]*uz_surf_cr[0]+0.25*flux_rho_l[0]*uz_surf_cl[0]-0.7071067811865475*avg_p_ij_z_r[0]+0.7071067811865475*avg_p_ij_z_l[0]+0.3535533905932737*apdq_rhouz_r[0]-0.3535533905932737*(apdq_rhouz_l[0]+amdq_rhouz_r[0])+0.3535533905932737*amdq_rhouz_l[0])*dx1; 
  outrhou2[1] += ((-0.4330127018922193*(flux_rho_r[1]*uz_surf_rl[1]+flux_rho_l[1]*uz_surf_lr[1]+flux_rho_r[1]*uz_surf_cr[1]+flux_rho_l[1]*uz_surf_cl[1]+flux_rho_r[0]*uz_surf_rl[0]+flux_rho_l[0]*uz_surf_lr[0]+flux_rho_r[0]*uz_surf_cr[0]+flux_rho_l[0]*uz_surf_cl[0]))-1.224744871391589*(avg_p_ij_z_r[0]+avg_p_ij_z_l[0])+0.6123724356957944*(apdq_rhouz_r[0]+apdq_rhouz_l[0])-0.6123724356957944*(amdq_rhouz_r[0]+amdq_rhouz_l[0]))*dx1; 
  outrhou2[2] += ((-0.223606797749979*(uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[2])+0.223606797749979*(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[2]-0.25*flux_rho_r[0]*uz_surf_rl[1]+0.25*flux_rho_l[0]*uz_surf_lr[1]-0.25*flux_rho_r[0]*uz_surf_cr[1]+0.25*flux_rho_l[0]*uz_surf_cl[1]-0.25*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+0.25*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1]-0.7071067811865475*avg_p_ij_z_r[1]+0.7071067811865475*avg_p_ij_z_l[1]+0.3535533905932737*apdq_rhouz_r[1]-0.3535533905932737*(apdq_rhouz_l[1]+amdq_rhouz_r[1])+0.3535533905932737*amdq_rhouz_l[1])*dx1; 
  outrhou2[3] += ((-0.3872983346207416*((uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[2]+(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[2]))-0.4330127018922193*(flux_rho_r[0]*uz_surf_rl[1]+flux_rho_l[0]*uz_surf_lr[1]+flux_rho_r[0]*uz_surf_cr[1]+flux_rho_l[0]*uz_surf_cl[1]+(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1])-1.224744871391589*(avg_p_ij_z_r[1]+avg_p_ij_z_l[1])+0.6123724356957944*(apdq_rhouz_r[1]+apdq_rhouz_l[1])-0.6123724356957944*(amdq_rhouz_r[1]+amdq_rhouz_l[1]))*dx1; 
  outrhou2[4] += ((-0.5590169943749475*flux_rho_r[1]*uz_surf_rl[1])+0.5590169943749475*flux_rho_l[1]*uz_surf_lr[1]-0.5590169943749475*flux_rho_r[1]*uz_surf_cr[1]+0.5590169943749475*flux_rho_l[1]*uz_surf_cl[1]-0.5590169943749475*flux_rho_r[0]*uz_surf_rl[0]+0.5590169943749475*flux_rho_l[0]*uz_surf_lr[0]-0.5590169943749475*flux_rho_r[0]*uz_surf_cr[0]+0.5590169943749475*flux_rho_l[0]*uz_surf_cl[0]-1.58113883008419*avg_p_ij_z_r[0]+1.58113883008419*avg_p_ij_z_l[0]+0.7905694150420947*apdq_rhouz_r[0]-0.7905694150420947*(apdq_rhouz_l[0]+amdq_rhouz_r[0])+0.7905694150420947*amdq_rhouz_l[0])*dx1; 
  outrhou2[5] += ((-0.25*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[2])+0.25*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[2]-0.7071067811865475*avg_p_ij_z_r[2]+0.7071067811865475*avg_p_ij_z_l[2]+0.3535533905932737*apdq_rhouz_r[2]-0.3535533905932737*(apdq_rhouz_l[2]+amdq_rhouz_r[2])+0.3535533905932737*amdq_rhouz_l[2]-0.223606797749979*flux_rho_r[1]*uz_surf_rl[1]+0.223606797749979*flux_rho_l[1]*uz_surf_lr[1]-0.223606797749979*flux_rho_r[1]*uz_surf_cr[1]+0.223606797749979*flux_rho_l[1]*uz_surf_cl[1])*dx1; 
  outrhou2[6] += ((-0.5000000000000001*(uz_surf_rl[1]+uz_surf_cr[1])*flux_rho_r[2])+0.5000000000000001*(uz_surf_lr[1]+uz_surf_cl[1])*flux_rho_l[2]-0.5590169943749476*flux_rho_r[0]*uz_surf_rl[1]+0.5590169943749476*flux_rho_l[0]*uz_surf_lr[1]-0.5590169943749476*flux_rho_r[0]*uz_surf_cr[1]+0.5590169943749476*flux_rho_l[0]*uz_surf_cl[1]-0.5590169943749476*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[1]+0.5590169943749476*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[1]-1.58113883008419*avg_p_ij_z_r[1]+1.58113883008419*avg_p_ij_z_l[1]+0.7905694150420948*apdq_rhouz_r[1]-0.7905694150420948*(apdq_rhouz_l[1]+amdq_rhouz_r[1])+0.7905694150420948*amdq_rhouz_l[1])*dx1; 
  outrhou2[7] += ((-0.4330127018922194*((uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[2]+(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[2]))-1.224744871391589*(avg_p_ij_z_r[2]+avg_p_ij_z_l[2])+0.6123724356957944*(apdq_rhouz_r[2]+apdq_rhouz_l[2])-0.6123724356957944*(amdq_rhouz_r[2]+amdq_rhouz_l[2])-0.3872983346207417*(flux_rho_r[1]*uz_surf_rl[1]+flux_rho_l[1]*uz_surf_lr[1]+flux_rho_r[1]*uz_surf_cr[1]+flux_rho_l[1]*uz_surf_cl[1]))*dx1; 
  outrhou2[8] += ((-0.5590169943749475*(uz_surf_rl[0]+uz_surf_cr[0])*flux_rho_r[2])+0.5590169943749475*(uz_surf_lr[0]+uz_surf_cl[0])*flux_rho_l[2]-1.58113883008419*avg_p_ij_z_r[2]+1.58113883008419*avg_p_ij_z_l[2]+0.7905694150420947*apdq_rhouz_r[2]-0.7905694150420947*(apdq_rhouz_l[2]+amdq_rhouz_r[2])+0.7905694150420947*amdq_rhouz_l[2]-0.5*flux_rho_r[1]*uz_surf_rl[1]+0.5*flux_rho_l[1]*uz_surf_lr[1]-0.5*flux_rho_r[1]*uz_surf_cr[1]+0.5*flux_rho_l[1]*uz_surf_cl[1])*dx1; 

  return 0.;

} 
