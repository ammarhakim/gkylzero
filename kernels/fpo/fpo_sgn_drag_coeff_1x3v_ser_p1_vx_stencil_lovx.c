#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH int fpo_sgn_drag_coeff_1x3v_vx_ser_p1_lovx(const double *drag_coeff_surf, double *sgn_drag_coeff_surf) {
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
  int const_sgn_alpha_surf = 0; 
  return const_sgn_alpha_surf; 
} 

