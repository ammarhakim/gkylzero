#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH int fpo_sgn_drag_coeff_1x3v_vx_ser_p2_lovx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface expansion of drag coefficient at LOWER cell boundary. 
  // sgn_drag_coeff_surf: Sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // returns const_sgn_drag_coeff: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *drag_coeff_surf_vx = &drag_coeff_surf[0]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[0]; 

  sgn_alpha_surf[0] = 0.0;
  sgn_alpha_surf[1] = 0.0;
  sgn_alpha_surf[2] = 0.0;
  sgn_alpha_surf[3] = 0.0;
  sgn_alpha_surf[4] = 0.0;
  sgn_alpha_surf[5] = 0.0;
  sgn_alpha_surf[6] = 0.0;
  sgn_alpha_surf[7] = 0.0;
  sgn_alpha_surf[8] = 0.0;
  sgn_alpha_surf[9] = 0.0;
  sgn_alpha_surf[10] = 0.0;
  sgn_alpha_surf[11] = 0.0;
  sgn_alpha_surf[12] = 0.0;
  sgn_alpha_surf[13] = 0.0;
  sgn_alpha_surf[14] = 0.0;
  sgn_alpha_surf[15] = 0.0;
  sgn_alpha_surf[16] = 0.0;
  sgn_alpha_surf[17] = 0.0;
  sgn_alpha_surf[18] = 0.0;
  sgn_alpha_surf[19] = 0.0;
  sgn_alpha_surf[20] = 0.0;
  sgn_alpha_surf[21] = 0.0;
  sgn_alpha_surf[22] = 0.0;
  sgn_alpha_surf[23] = 0.0;
  sgn_alpha_surf[24] = 0.0;
  sgn_alpha_surf[25] = 0.0;
  sgn_alpha_surf[26] = 0.0;
  int const_sgn_alpha_surf = 0; 
  return const_sgn_alpha_surf; 
} 

