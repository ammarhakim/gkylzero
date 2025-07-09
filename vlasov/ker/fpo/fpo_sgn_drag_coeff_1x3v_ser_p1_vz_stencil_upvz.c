#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_sgn_drag_coeff_1x3v_vz_ser_p1_upvz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf, int *const_sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // const_sgn_drag_coeff_surf: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *alpha_surf = &drag_coeff_surf[32]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[36]; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.42426406871192823*alpha_surf[15]-0.42426406871192823*alpha_surf[14]-0.31622776601683783*alpha_surf[13]+0.31622776601683783*alpha_surf[12]+0.42426406871192823*alpha_surf[11]-0.42426406871192857*alpha_surf[10]-0.31622776601683783*alpha_surf[9]+0.31622776601683783*alpha_surf[8]-0.6363961030678927*alpha_surf[7]+0.6363961030678927*alpha_surf[6]+0.4743416490252568*(alpha_surf[5]+alpha_surf[4])-0.4743416490252568*(alpha_surf[3]+alpha_surf[2])-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (-(0.5303300858899102*alpha_surf[15])+0.5303300858899102*alpha_surf[14]+0.3952847075210473*alpha_surf[13]-0.3952847075210471*alpha_surf[12]-0.31622776601683783*alpha_surf[9]+0.31622776601683783*alpha_surf[8]+0.4743416490252568*alpha_surf[4]-0.4743416490252568*alpha_surf[2]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.42426406871192823*alpha_surf[15]-0.42426406871192823*alpha_surf[14]-0.31622776601683783*alpha_surf[13]+0.31622776601683783*alpha_surf[12]-0.42426406871192823*alpha_surf[11]+0.42426406871192857*alpha_surf[10]-0.31622776601683783*alpha_surf[9]+0.31622776601683783*alpha_surf[8]+0.6363961030678927*alpha_surf[7]-0.6363961030678927*alpha_surf[6]-0.4743416490252568*alpha_surf[5]+0.4743416490252568*(alpha_surf[4]+alpha_surf[3])-0.4743416490252568*alpha_surf[2]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alpha_surf[13])+0.31622776601683783*alpha_surf[12]-0.5303300858899102*alpha_surf[11]+0.5303300858899102*alpha_surf[10]+0.3952847075210473*alpha_surf[9]-0.3952847075210471*alpha_surf[8]+0.4743416490252568*alpha_surf[5]-0.4743416490252568*alpha_surf[3]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3952847075210473*alpha_surf[13]-0.3952847075210471*alpha_surf[12]+0.3952847075210473*alpha_surf[9]-0.3952847075210471*alpha_surf[8]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alpha_surf[13])+0.31622776601683783*alpha_surf[12]+0.5303300858899102*alpha_surf[11]-0.5303300858899102*alpha_surf[10]+0.3952847075210473*alpha_surf[9]-0.3952847075210471*alpha_surf[8]-0.4743416490252568*alpha_surf[5]+0.4743416490252568*alpha_surf[3]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.42426406871192823*alpha_surf[15])+0.42426406871192823*alpha_surf[14]-0.31622776601683783*alpha_surf[13]+0.31622776601683783*alpha_surf[12]+0.42426406871192823*alpha_surf[11]-0.42426406871192857*alpha_surf[10]-0.31622776601683783*alpha_surf[9]+0.31622776601683783*alpha_surf[8]+0.6363961030678927*alpha_surf[7]-0.6363961030678927*alpha_surf[6]+0.4743416490252568*alpha_surf[5]-0.4743416490252568*(alpha_surf[4]+alpha_surf[3])+0.4743416490252568*alpha_surf[2]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alpha_surf[15]-0.5303300858899102*alpha_surf[14]+0.3952847075210473*alpha_surf[13]-0.3952847075210471*alpha_surf[12]-0.31622776601683783*alpha_surf[9]+0.31622776601683783*alpha_surf[8]-0.4743416490252568*alpha_surf[4]+0.4743416490252568*alpha_surf[2]-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.42426406871192823*alpha_surf[15])+0.42426406871192823*alpha_surf[14]-0.31622776601683783*alpha_surf[13]+0.31622776601683783*alpha_surf[12]-0.42426406871192823*alpha_surf[11]+0.42426406871192857*alpha_surf[10]-0.31622776601683783*alpha_surf[9]+0.31622776601683783*alpha_surf[8]-0.6363961030678927*alpha_surf[7]+0.6363961030678927*alpha_surf[6]-0.4743416490252568*(alpha_surf[5]+alpha_surf[4])+0.4743416490252568*(alpha_surf[3]+alpha_surf[2])-0.3535533905932734*alpha_surf[1]+0.3535533905932734*alpha_surf[0] > 0.) 
    sgn_alpha_surf[8] = 1.0; 
  else  
    sgn_alpha_surf[8] = -1.0; 
  
  if (sgn_alpha_surf[8] == sgn_alpha_surf[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.42426406871192823*(alpha_surf[15]+alpha_surf[14]))+0.31622776601683783*(alpha_surf[13]+alpha_surf[12])-0.42426406871192823*alpha_surf[11]-0.42426406871192857*alpha_surf[10]+0.31622776601683783*(alpha_surf[9]+alpha_surf[8])+0.6363961030678927*(alpha_surf[7]+alpha_surf[6])-0.4743416490252568*(alpha_surf[5]+alpha_surf[4]+alpha_surf[3]+alpha_surf[2])+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[9] = 1.0; 
  else  
    sgn_alpha_surf[9] = -1.0; 
  
  if (sgn_alpha_surf[9] == sgn_alpha_surf[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(alpha_surf[15]+alpha_surf[14])-0.3952847075210473*alpha_surf[13]-0.3952847075210471*alpha_surf[12]+0.31622776601683783*(alpha_surf[9]+alpha_surf[8])-0.4743416490252568*(alpha_surf[4]+alpha_surf[2])+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[10] = 1.0; 
  else  
    sgn_alpha_surf[10] = -1.0; 
  
  if (sgn_alpha_surf[10] == sgn_alpha_surf[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.42426406871192823*(alpha_surf[15]+alpha_surf[14]))+0.31622776601683783*(alpha_surf[13]+alpha_surf[12])+0.42426406871192823*alpha_surf[11]+0.42426406871192857*alpha_surf[10]+0.31622776601683783*(alpha_surf[9]+alpha_surf[8])-0.6363961030678927*(alpha_surf[7]+alpha_surf[6])+0.4743416490252568*alpha_surf[5]-0.4743416490252568*alpha_surf[4]+0.4743416490252568*alpha_surf[3]-0.4743416490252568*alpha_surf[2]+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[11] = 1.0; 
  else  
    sgn_alpha_surf[11] = -1.0; 
  
  if (sgn_alpha_surf[11] == sgn_alpha_surf[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alpha_surf[13]+alpha_surf[12])+0.5303300858899102*(alpha_surf[11]+alpha_surf[10])-0.3952847075210473*alpha_surf[9]-0.3952847075210471*alpha_surf[8]-0.4743416490252568*(alpha_surf[5]+alpha_surf[3])+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[12] = 1.0; 
  else  
    sgn_alpha_surf[12] = -1.0; 
  
  if (sgn_alpha_surf[12] == sgn_alpha_surf[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3952847075210473*alpha_surf[13])-0.3952847075210471*alpha_surf[12]-0.3952847075210473*alpha_surf[9]-0.3952847075210471*alpha_surf[8]+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[13] = 1.0; 
  else  
    sgn_alpha_surf[13] = -1.0; 
  
  if (sgn_alpha_surf[13] == sgn_alpha_surf[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alpha_surf[13]+alpha_surf[12])-0.5303300858899102*(alpha_surf[11]+alpha_surf[10])-0.3952847075210473*alpha_surf[9]-0.3952847075210471*alpha_surf[8]+0.4743416490252568*(alpha_surf[5]+alpha_surf[3])+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[14] = 1.0; 
  else  
    sgn_alpha_surf[14] = -1.0; 
  
  if (sgn_alpha_surf[14] == sgn_alpha_surf[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.42426406871192823*(alpha_surf[15]+alpha_surf[14])+0.31622776601683783*(alpha_surf[13]+alpha_surf[12])-0.42426406871192823*alpha_surf[11]-0.42426406871192857*alpha_surf[10]+0.31622776601683783*(alpha_surf[9]+alpha_surf[8])-0.6363961030678927*(alpha_surf[7]+alpha_surf[6])-0.4743416490252568*alpha_surf[5]+0.4743416490252568*alpha_surf[4]-0.4743416490252568*alpha_surf[3]+0.4743416490252568*alpha_surf[2]+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[15] = 1.0; 
  else  
    sgn_alpha_surf[15] = -1.0; 
  
  if (sgn_alpha_surf[15] == sgn_alpha_surf[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(alpha_surf[15]+alpha_surf[14]))-0.3952847075210473*alpha_surf[13]-0.3952847075210471*alpha_surf[12]+0.31622776601683783*(alpha_surf[9]+alpha_surf[8])+0.4743416490252568*(alpha_surf[4]+alpha_surf[2])+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[16] = 1.0; 
  else  
    sgn_alpha_surf[16] = -1.0; 
  
  if (sgn_alpha_surf[16] == sgn_alpha_surf[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.42426406871192823*(alpha_surf[15]+alpha_surf[14])+0.31622776601683783*(alpha_surf[13]+alpha_surf[12])+0.42426406871192823*alpha_surf[11]+0.42426406871192857*alpha_surf[10]+0.31622776601683783*(alpha_surf[9]+alpha_surf[8])+0.6363961030678927*(alpha_surf[7]+alpha_surf[6])+0.4743416490252568*(alpha_surf[5]+alpha_surf[4]+alpha_surf[3]+alpha_surf[2])+0.3535533905932734*(alpha_surf[1]+alpha_surf[0]) > 0.) 
    sgn_alpha_surf[17] = 1.0; 
  else  
    sgn_alpha_surf[17] = -1.0; 
  
  if (sgn_alpha_surf[17] == sgn_alpha_surf[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  *const_sgn_drag_coeff_surf = const_sgn_alpha_surf; 
} 
