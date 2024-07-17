#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH int fpo_sgn_drag_coeff_1x3v_vz_ser_p2_upvz(const double *drag_coeff_surf, double *sgn_drag_coeff_surf) {
  // drag_coeff_surf: Surface expansion of drag coefficient at LOWER cell boundary. 
  // sgn_drag_coeff_surf: Sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // returns const_sgn_drag_coeff: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  const double *drag_coeff_surf_vz = &drag_coeff_surf[40]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[54]; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.5692099788303081*(drag_coeff_surf_vz[19]+drag_coeff_surf_vz[18]+drag_coeff_surf_vz[17])-0.42426406871192823*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[15])-0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])-0.42426406871192823*drag_coeff_surf_vz[12]-0.42426406871192857*drag_coeff_surf_vz[11]-0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.6363961030678927*(drag_coeff_surf_vz[6]+drag_coeff_surf_vz[5]+drag_coeff_surf_vz[4])-0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[2]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (-(0.7115124735378848*drag_coeff_surf_vz[19])+0.5303300858899102*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[15])-0.42426406871192823*drag_coeff_surf_vz[12]-0.42426406871192857*drag_coeff_surf_vz[11]-0.3952847075210471*drag_coeff_surf_vz[9]+0.31622776601683783*(drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.6363961030678927*drag_coeff_surf_vz[4]-0.4743416490252568*(drag_coeff_surf_vz[2]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303081*drag_coeff_surf_vz[19]-0.5692099788303081*(drag_coeff_surf_vz[18]+drag_coeff_surf_vz[17])-0.42426406871192823*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[15])+0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])-0.42426406871192823*drag_coeff_surf_vz[12]-0.42426406871192857*drag_coeff_surf_vz[11]+0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.6363961030678927*(drag_coeff_surf_vz[6]+drag_coeff_surf_vz[5])+0.6363961030678927*drag_coeff_surf_vz[4]+0.4743416490252568*drag_coeff_surf_vz[3]-0.4743416490252568*(drag_coeff_surf_vz[2]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*drag_coeff_surf_vz[18])-0.42426406871192823*drag_coeff_surf_vz[15]+0.5303300858899102*drag_coeff_surf_vz[14]-0.42426406871192857*drag_coeff_surf_vz[13]+0.5303300858899102*drag_coeff_surf_vz[12]+0.31622776601683783*drag_coeff_surf_vz[9]-0.3952847075210471*drag_coeff_surf_vz[8]+0.31622776601683783*drag_coeff_surf_vz[7]+0.6363961030678927*drag_coeff_surf_vz[5]-0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(drag_coeff_surf_vz[15]+drag_coeff_surf_vz[12])-0.3952847075210471*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8])+0.31622776601683783*drag_coeff_surf_vz[7]-0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*drag_coeff_surf_vz[18]-0.42426406871192823*drag_coeff_surf_vz[15]-0.5303300858899102*drag_coeff_surf_vz[14]+0.42426406871192857*drag_coeff_surf_vz[13]+0.5303300858899102*drag_coeff_surf_vz[12]+0.31622776601683783*drag_coeff_surf_vz[9]-0.3952847075210471*drag_coeff_surf_vz[8]+0.31622776601683783*drag_coeff_surf_vz[7]-0.6363961030678927*drag_coeff_surf_vz[5]+0.4743416490252568*drag_coeff_surf_vz[3]-0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*drag_coeff_surf_vz[19])+0.5692099788303081*drag_coeff_surf_vz[18]-0.5692099788303081*drag_coeff_surf_vz[17]+0.42426406871192823*drag_coeff_surf_vz[16]-0.42426406871192823*drag_coeff_surf_vz[15]-0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])-0.42426406871192823*drag_coeff_surf_vz[12]+0.42426406871192857*drag_coeff_surf_vz[11]+0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.6363961030678927*drag_coeff_surf_vz[6]+0.6363961030678927*drag_coeff_surf_vz[5]-0.6363961030678927*drag_coeff_surf_vz[4]-0.4743416490252568*drag_coeff_surf_vz[3]+0.4743416490252568*drag_coeff_surf_vz[2]-0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*drag_coeff_surf_vz[19]-0.5303300858899102*drag_coeff_surf_vz[16]+0.5303300858899102*drag_coeff_surf_vz[15]-0.42426406871192823*drag_coeff_surf_vz[12]+0.42426406871192857*drag_coeff_surf_vz[11]-0.3952847075210471*drag_coeff_surf_vz[9]+0.31622776601683783*(drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.6363961030678927*drag_coeff_surf_vz[4]+0.4743416490252568*drag_coeff_surf_vz[2]-0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*(drag_coeff_surf_vz[19]+drag_coeff_surf_vz[18]))+0.5692099788303081*drag_coeff_surf_vz[17]+0.42426406871192823*drag_coeff_surf_vz[16]-0.42426406871192823*drag_coeff_surf_vz[15]+0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])-0.42426406871192823*drag_coeff_surf_vz[12]+0.42426406871192857*drag_coeff_surf_vz[11]-0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.6363961030678927*drag_coeff_surf_vz[6]-0.6363961030678927*(drag_coeff_surf_vz[5]+drag_coeff_surf_vz[4])+0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[2])-0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[8] = 1.0; 
  else  
    sgn_alpha_surf[8] = -1.0; 
  
  if (sgn_alpha_surf[8] == sgn_alpha_surf[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*drag_coeff_surf_vz[17])-0.42426406871192823*drag_coeff_surf_vz[16]-0.42426406871192857*drag_coeff_surf_vz[14]+0.5303300858899102*(drag_coeff_surf_vz[13]+drag_coeff_surf_vz[11])+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8])-0.3952847075210471*drag_coeff_surf_vz[7]+0.6363961030678927*drag_coeff_surf_vz[6]-0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[2])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[9] = 1.0; 
  else  
    sgn_alpha_surf[9] = -1.0; 
  
  if (sgn_alpha_surf[9] == sgn_alpha_surf[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[11])-0.3952847075210471*drag_coeff_surf_vz[9]+0.31622776601683783*drag_coeff_surf_vz[8]-0.3952847075210471*drag_coeff_surf_vz[7]-0.4743416490252568*drag_coeff_surf_vz[2]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[10] = 1.0; 
  else  
    sgn_alpha_surf[10] = -1.0; 
  
  if (sgn_alpha_surf[10] == sgn_alpha_surf[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*drag_coeff_surf_vz[17]-0.42426406871192823*drag_coeff_surf_vz[16]+0.42426406871192857*drag_coeff_surf_vz[14]-0.5303300858899102*drag_coeff_surf_vz[13]+0.5303300858899102*drag_coeff_surf_vz[11]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8])-0.3952847075210471*drag_coeff_surf_vz[7]-0.6363961030678927*drag_coeff_surf_vz[6]+0.4743416490252568*drag_coeff_surf_vz[3]-0.4743416490252568*drag_coeff_surf_vz[2]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[11] = 1.0; 
  else  
    sgn_alpha_surf[11] = -1.0; 
  
  if (sgn_alpha_surf[11] == sgn_alpha_surf[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])+0.31622776601683783*drag_coeff_surf_vz[9]-0.3952847075210471*(drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.4743416490252568*drag_coeff_surf_vz[3]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[12] = 1.0; 
  else  
    sgn_alpha_surf[12] = -1.0; 
  
  if (sgn_alpha_surf[12] == sgn_alpha_surf[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*drag_coeff_surf_vz[0]-0.3952847075210471*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7]) > 0.) 
    sgn_alpha_surf[13] = 1.0; 
  else  
    sgn_alpha_surf[13] = -1.0; 
  
  if (sgn_alpha_surf[13] == sgn_alpha_surf[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13]))+0.31622776601683783*drag_coeff_surf_vz[9]-0.3952847075210471*(drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.4743416490252568*drag_coeff_surf_vz[3]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[14] = 1.0; 
  else  
    sgn_alpha_surf[14] = -1.0; 
  
  if (sgn_alpha_surf[14] == sgn_alpha_surf[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*drag_coeff_surf_vz[17]+0.42426406871192823*drag_coeff_surf_vz[16]-0.42426406871192857*drag_coeff_surf_vz[14]+0.5303300858899102*drag_coeff_surf_vz[13]-0.5303300858899102*drag_coeff_surf_vz[11]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8])-0.3952847075210471*drag_coeff_surf_vz[7]-0.6363961030678927*drag_coeff_surf_vz[6]-0.4743416490252568*drag_coeff_surf_vz[3]+0.4743416490252568*drag_coeff_surf_vz[2]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[15] = 1.0; 
  else  
    sgn_alpha_surf[15] = -1.0; 
  
  if (sgn_alpha_surf[15] == sgn_alpha_surf[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[11]))-0.3952847075210471*drag_coeff_surf_vz[9]+0.31622776601683783*drag_coeff_surf_vz[8]-0.3952847075210471*drag_coeff_surf_vz[7]+0.4743416490252568*drag_coeff_surf_vz[2]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[16] = 1.0; 
  else  
    sgn_alpha_surf[16] = -1.0; 
  
  if (sgn_alpha_surf[16] == sgn_alpha_surf[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*drag_coeff_surf_vz[17])+0.42426406871192823*drag_coeff_surf_vz[16]+0.42426406871192857*drag_coeff_surf_vz[14]-0.5303300858899102*(drag_coeff_surf_vz[13]+drag_coeff_surf_vz[11])+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8])-0.3952847075210471*drag_coeff_surf_vz[7]+0.6363961030678927*drag_coeff_surf_vz[6]+0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[2])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[17] = 1.0; 
  else  
    sgn_alpha_surf[17] = -1.0; 
  
  if (sgn_alpha_surf[17] == sgn_alpha_surf[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*(drag_coeff_surf_vz[19]+drag_coeff_surf_vz[18]))+0.5692099788303081*drag_coeff_surf_vz[17]-0.42426406871192823*drag_coeff_surf_vz[16]+0.42426406871192823*drag_coeff_surf_vz[15]-0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])+0.42426406871192823*drag_coeff_surf_vz[12]-0.42426406871192857*drag_coeff_surf_vz[11]+0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.6363961030678927*drag_coeff_surf_vz[6]-0.6363961030678927*(drag_coeff_surf_vz[5]+drag_coeff_surf_vz[4])-0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[2])+0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[18] = 1.0; 
  else  
    sgn_alpha_surf[18] = -1.0; 
  
  if (sgn_alpha_surf[18] == sgn_alpha_surf[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*drag_coeff_surf_vz[19]+0.5303300858899102*drag_coeff_surf_vz[16]-0.5303300858899102*drag_coeff_surf_vz[15]+0.42426406871192823*drag_coeff_surf_vz[12]-0.42426406871192857*drag_coeff_surf_vz[11]-0.3952847075210471*drag_coeff_surf_vz[9]+0.31622776601683783*(drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.6363961030678927*drag_coeff_surf_vz[4]-0.4743416490252568*drag_coeff_surf_vz[2]+0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[19] = 1.0; 
  else  
    sgn_alpha_surf[19] = -1.0; 
  
  if (sgn_alpha_surf[19] == sgn_alpha_surf[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*drag_coeff_surf_vz[19])+0.5692099788303081*drag_coeff_surf_vz[18]-0.5692099788303081*drag_coeff_surf_vz[17]-0.42426406871192823*drag_coeff_surf_vz[16]+0.42426406871192823*drag_coeff_surf_vz[15]+0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])+0.42426406871192823*drag_coeff_surf_vz[12]-0.42426406871192857*drag_coeff_surf_vz[11]-0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.6363961030678927*drag_coeff_surf_vz[6]+0.6363961030678927*drag_coeff_surf_vz[5]-0.6363961030678927*drag_coeff_surf_vz[4]+0.4743416490252568*drag_coeff_surf_vz[3]-0.4743416490252568*drag_coeff_surf_vz[2]+0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[20] = 1.0; 
  else  
    sgn_alpha_surf[20] = -1.0; 
  
  if (sgn_alpha_surf[20] == sgn_alpha_surf[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*drag_coeff_surf_vz[18]+0.42426406871192823*drag_coeff_surf_vz[15]+0.5303300858899102*drag_coeff_surf_vz[14]-0.42426406871192857*drag_coeff_surf_vz[13]-0.5303300858899102*drag_coeff_surf_vz[12]+0.31622776601683783*drag_coeff_surf_vz[9]-0.3952847075210471*drag_coeff_surf_vz[8]+0.31622776601683783*drag_coeff_surf_vz[7]-0.6363961030678927*drag_coeff_surf_vz[5]-0.4743416490252568*drag_coeff_surf_vz[3]+0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[21] = 1.0; 
  else  
    sgn_alpha_surf[21] = -1.0; 
  
  if (sgn_alpha_surf[21] == sgn_alpha_surf[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(drag_coeff_surf_vz[15]+drag_coeff_surf_vz[12]))-0.3952847075210471*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8])+0.31622776601683783*drag_coeff_surf_vz[7]+0.4743416490252568*drag_coeff_surf_vz[1]+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[22] = 1.0; 
  else  
    sgn_alpha_surf[22] = -1.0; 
  
  if (sgn_alpha_surf[22] == sgn_alpha_surf[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*drag_coeff_surf_vz[18])+0.42426406871192823*drag_coeff_surf_vz[15]-0.5303300858899102*drag_coeff_surf_vz[14]+0.42426406871192857*drag_coeff_surf_vz[13]-0.5303300858899102*drag_coeff_surf_vz[12]+0.31622776601683783*drag_coeff_surf_vz[9]-0.3952847075210471*drag_coeff_surf_vz[8]+0.31622776601683783*drag_coeff_surf_vz[7]+0.6363961030678927*drag_coeff_surf_vz[5]+0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[23] = 1.0; 
  else  
    sgn_alpha_surf[23] = -1.0; 
  
  if (sgn_alpha_surf[23] == sgn_alpha_surf[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303081*drag_coeff_surf_vz[19]-0.5692099788303081*(drag_coeff_surf_vz[18]+drag_coeff_surf_vz[17])+0.42426406871192823*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[15])-0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])+0.42426406871192823*drag_coeff_surf_vz[12]+0.42426406871192857*drag_coeff_surf_vz[11]-0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])-0.6363961030678927*(drag_coeff_surf_vz[6]+drag_coeff_surf_vz[5])+0.6363961030678927*drag_coeff_surf_vz[4]-0.4743416490252568*drag_coeff_surf_vz[3]+0.4743416490252568*(drag_coeff_surf_vz[2]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[24] = 1.0; 
  else  
    sgn_alpha_surf[24] = -1.0; 
  
  if (sgn_alpha_surf[24] == sgn_alpha_surf[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*drag_coeff_surf_vz[19])-0.5303300858899102*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[15])+0.42426406871192823*drag_coeff_surf_vz[12]+0.42426406871192857*drag_coeff_surf_vz[11]-0.3952847075210471*drag_coeff_surf_vz[9]+0.31622776601683783*(drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.6363961030678927*drag_coeff_surf_vz[4]+0.4743416490252568*(drag_coeff_surf_vz[2]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[25] = 1.0; 
  else  
    sgn_alpha_surf[25] = -1.0; 
  
  if (sgn_alpha_surf[25] == sgn_alpha_surf[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303081*(drag_coeff_surf_vz[19]+drag_coeff_surf_vz[18]+drag_coeff_surf_vz[17])+0.42426406871192823*(drag_coeff_surf_vz[16]+drag_coeff_surf_vz[15])+0.42426406871192857*(drag_coeff_surf_vz[14]+drag_coeff_surf_vz[13])+0.42426406871192823*drag_coeff_surf_vz[12]+0.42426406871192857*drag_coeff_surf_vz[11]+0.8538149682454614*drag_coeff_surf_vz[10]+0.31622776601683783*(drag_coeff_surf_vz[9]+drag_coeff_surf_vz[8]+drag_coeff_surf_vz[7])+0.6363961030678927*(drag_coeff_surf_vz[6]+drag_coeff_surf_vz[5]+drag_coeff_surf_vz[4])+0.4743416490252568*(drag_coeff_surf_vz[3]+drag_coeff_surf_vz[2]+drag_coeff_surf_vz[1])+0.3535533905932734*drag_coeff_surf_vz[0] > 0.) 
    sgn_alpha_surf[26] = 1.0; 
  else  
    sgn_alpha_surf[26] = -1.0; 
  
  if (sgn_alpha_surf[26] == sgn_alpha_surf[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 
} 

