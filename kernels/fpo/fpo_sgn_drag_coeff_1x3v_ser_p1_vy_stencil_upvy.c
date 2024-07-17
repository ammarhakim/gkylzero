#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH int fpo_sgn_drag_coeff_1x3v_vy_ser_p1_upvy(const double *drag_coeff_surf, double *sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface expansion of drag coefficient at LOWER cell boundary. 
  // sgn_drag_coeff_surf: Sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // returns const_sgn_drag_coeff: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *drag_coeff_surf_vy = &drag_coeff_surf[8]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[8]; 

  int const_sgn_alpha_surf = 1;  
  
  if (-(0.3535533905932734*drag_coeff_surf_vy[7])+0.3535533905932734*(drag_coeff_surf_vy[6]+drag_coeff_surf_vy[5]+drag_coeff_surf_vy[4])-0.3535533905932734*(drag_coeff_surf_vy[3]+drag_coeff_surf_vy[2]+drag_coeff_surf_vy[1])+0.3535533905932734*drag_coeff_surf_vy[0] > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (0.3535533905932734*drag_coeff_surf_vy[7]-0.3535533905932734*(drag_coeff_surf_vy[6]+drag_coeff_surf_vy[5])+0.3535533905932734*(drag_coeff_surf_vy[4]+drag_coeff_surf_vy[3])-0.3535533905932734*(drag_coeff_surf_vy[2]+drag_coeff_surf_vy[1])+0.3535533905932734*drag_coeff_surf_vy[0] > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*drag_coeff_surf_vy[7]-0.3535533905932734*drag_coeff_surf_vy[6]+0.3535533905932734*drag_coeff_surf_vy[5]-0.3535533905932734*(drag_coeff_surf_vy[4]+drag_coeff_surf_vy[3])+0.3535533905932734*drag_coeff_surf_vy[2]-0.3535533905932734*drag_coeff_surf_vy[1]+0.3535533905932734*drag_coeff_surf_vy[0] > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3535533905932734*drag_coeff_surf_vy[7])+0.3535533905932734*drag_coeff_surf_vy[6]-0.3535533905932734*(drag_coeff_surf_vy[5]+drag_coeff_surf_vy[4])+0.3535533905932734*(drag_coeff_surf_vy[3]+drag_coeff_surf_vy[2])-0.3535533905932734*drag_coeff_surf_vy[1]+0.3535533905932734*drag_coeff_surf_vy[0] > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(drag_coeff_surf_vy[7]+drag_coeff_surf_vy[6])-0.3535533905932734*(drag_coeff_surf_vy[5]+drag_coeff_surf_vy[4]+drag_coeff_surf_vy[3]+drag_coeff_surf_vy[2])+0.3535533905932734*(drag_coeff_surf_vy[1]+drag_coeff_surf_vy[0]) > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3535533905932734*(drag_coeff_surf_vy[7]+drag_coeff_surf_vy[6]))+0.3535533905932734*drag_coeff_surf_vy[5]-0.3535533905932734*drag_coeff_surf_vy[4]+0.3535533905932734*drag_coeff_surf_vy[3]-0.3535533905932734*drag_coeff_surf_vy[2]+0.3535533905932734*(drag_coeff_surf_vy[1]+drag_coeff_surf_vy[0]) > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3535533905932734*(drag_coeff_surf_vy[7]+drag_coeff_surf_vy[6]+drag_coeff_surf_vy[5]))+0.3535533905932734*drag_coeff_surf_vy[4]-0.3535533905932734*drag_coeff_surf_vy[3]+0.3535533905932734*(drag_coeff_surf_vy[2]+drag_coeff_surf_vy[1]+drag_coeff_surf_vy[0]) > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(drag_coeff_surf_vy[7]+drag_coeff_surf_vy[6]+drag_coeff_surf_vy[5]+drag_coeff_surf_vy[4]+drag_coeff_surf_vy[3]+drag_coeff_surf_vy[2]+drag_coeff_surf_vy[1]+drag_coeff_surf_vy[0]) > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 
} 

