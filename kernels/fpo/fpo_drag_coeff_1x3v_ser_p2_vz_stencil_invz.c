#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH int fpo_drag_coeff_1x3v_vz_ser_p2_invz(const double *dxv, const double *gamma, 
    const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, 
    double *drag_coeff_surf, double *sgn_drag_coeff_surf
    ) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: Sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // returns const_sgn_drag_coeff: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[3]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  const double* H_R = fpo_h_stencil[2]; 
  
  const double *dHdv_surf_C = &fpo_dhdv_surf[40]; 
  
  double *out = &drag_coeff[96]; 
  double *out_surf = &drag_coeff_surf[40]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[54]; 
  
  out[0] = (0.489139870078079*H_R[14]-0.489139870078079*H_L[14]-0.7036456405748563*H_R[4]-0.7036456405748563*H_L[4]+1.4072912811497127*H_C[4]+0.5*H_R[0]-0.5*H_L[0])*dv1*gamma_avg; 
  out[1] = (0.4891398700780789*H_R[28]-0.4891398700780789*H_L[28]-0.7036456405748562*H_R[8]-0.7036456405748562*H_L[8]+1.4072912811497125*H_C[8]+0.5*H_R[1]-0.5*H_L[1])*dv1*gamma_avg; 
  out[2] = (0.4891398700780789*H_R[29]-0.4891398700780789*H_L[29]-0.7036456405748562*H_R[9]-0.7036456405748562*H_L[9]+1.4072912811497125*H_C[9]+0.5*H_R[2]-0.5*H_L[2])*dv1*gamma_avg; 
  out[3] = (0.4891398700780789*H_R[30]-0.4891398700780789*H_L[30]-0.7036456405748562*H_R[10]-0.7036456405748562*H_L[10]+1.4072912811497125*H_C[10]+0.5*H_R[3]-0.5*H_L[3])*dv1*gamma_avg; 
  out[4] = (0.8472151069828725*H_R[14]+0.8472151069828725*H_L[14]+1.694430213965745*H_C[14]-1.21875*H_R[4]+1.21875*H_L[4]+0.8660254037844386*H_R[0]+0.8660254037844386*H_L[0]-1.7320508075688772*H_C[0])*dv1*gamma_avg; 
  out[5] = (0.489139870078079*H_R[41]-0.489139870078079*H_L[41]-0.7036456405748563*H_R[16]-0.7036456405748563*H_L[16]+1.4072912811497127*H_C[16]+0.5*H_R[5]-0.5*H_L[5])*dv1*gamma_avg; 
  out[6] = (0.489139870078079*H_R[42]-0.489139870078079*H_L[42]-0.7036456405748563*H_R[17]-0.7036456405748563*H_L[17]+1.4072912811497127*H_C[17]+0.5*H_R[6]-0.5*H_L[6])*dv1*gamma_avg; 
  out[7] = (0.489139870078079*H_R[43]-0.489139870078079*H_L[43]-0.7036456405748563*H_R[18]-0.7036456405748563*H_L[18]+1.4072912811497127*H_C[18]+0.5*H_R[7]-0.5*H_L[7])*dv1*gamma_avg; 
  out[8] = (0.8472151069828725*H_R[28]+0.8472151069828725*H_L[28]+1.694430213965745*H_C[28]-1.21875*H_R[8]+1.21875*H_L[8]+0.8660254037844386*H_R[1]+0.8660254037844386*H_L[1]-1.7320508075688772*H_C[1])*dv1*gamma_avg; 
  out[9] = (0.8472151069828725*H_R[29]+0.8472151069828725*H_L[29]+1.694430213965745*H_C[29]-1.21875*H_R[9]+1.21875*H_L[9]+0.8660254037844386*H_R[2]+0.8660254037844386*H_L[2]-1.7320508075688772*H_C[2])*dv1*gamma_avg; 
  out[10] = (0.8472151069828725*H_R[30]+0.8472151069828725*H_L[30]+1.694430213965745*H_C[30]-1.21875*H_R[10]+1.21875*H_L[10]+0.8660254037844386*H_R[3]+0.8660254037844386*H_L[3]-1.7320508075688772*H_C[3])*dv1*gamma_avg; 
  out[11] = (-(0.7036456405748565*H_R[25])-0.7036456405748565*H_L[25]+1.407291281149713*H_C[25]+0.5*H_R[11]-0.5*H_L[11])*dv1*gamma_avg; 
  out[12] = (-(0.7036456405748565*H_R[26])-0.7036456405748565*H_L[26]+1.407291281149713*H_C[26]+0.5*H_R[12]-0.5*H_L[12])*dv1*gamma_avg; 
  out[13] = (-(0.7036456405748565*H_R[27])-0.7036456405748565*H_L[27]+1.407291281149713*H_C[27]+0.5*H_R[13]-0.5*H_L[13])*dv1*gamma_avg; 
  out[14] = (1.09375*H_R[14]-1.09375*H_L[14]-1.5733994843967631*H_R[4]-1.5733994843967631*H_L[4]-4.599167723621307*H_C[4]+1.118033988749895*H_R[0]-1.118033988749895*H_L[0])*dv1*gamma_avg; 
  out[15] = (0.4891398700780789*H_R[47]-0.4891398700780789*H_L[47]-0.7036456405748562*H_R[31]-0.7036456405748562*H_L[31]+1.4072912811497125*H_C[31]+0.5*H_R[15]-0.5*H_L[15])*dv1*gamma_avg; 
  out[16] = (0.8472151069828725*H_R[41]+0.8472151069828725*H_L[41]+1.694430213965745*H_C[41]-1.21875*H_R[16]+1.21875*H_L[16]+0.8660254037844386*H_R[5]+0.8660254037844386*H_L[5]-1.7320508075688772*H_C[5])*dv1*gamma_avg; 
  out[17] = (0.8472151069828725*H_R[42]+0.8472151069828725*H_L[42]+1.694430213965745*H_C[42]-1.21875*H_R[17]+1.21875*H_L[17]+0.8660254037844386*H_R[6]+0.8660254037844386*H_L[6]-1.7320508075688772*H_C[6])*dv1*gamma_avg; 
  out[18] = (0.8472151069828725*H_R[43]+0.8472151069828725*H_L[43]+1.694430213965745*H_C[43]-1.21875*H_R[18]+1.21875*H_L[18]+0.8660254037844386*H_R[7]+0.8660254037844386*H_L[7]-1.7320508075688772*H_C[7])*dv1*gamma_avg; 
  out[19] = (-(0.7036456405748565*H_R[35])-0.7036456405748565*H_L[35]+1.407291281149713*H_C[35]+0.5*H_R[19]-0.5*H_L[19])*dv1*gamma_avg; 
  out[20] = (-(0.7036456405748565*H_R[36])-0.7036456405748565*H_L[36]+1.407291281149713*H_C[36]+0.5*H_R[20]-0.5*H_L[20])*dv1*gamma_avg; 
  out[21] = (-(0.7036456405748565*H_R[37])-0.7036456405748565*H_L[37]+1.407291281149713*H_C[37]+0.5*H_R[21]-0.5*H_L[21])*dv1*gamma_avg; 
  out[22] = (-(0.7036456405748565*H_R[38])-0.7036456405748565*H_L[38]+1.407291281149713*H_C[38]+0.5*H_R[22]-0.5*H_L[22])*dv1*gamma_avg; 
  out[23] = (-(0.7036456405748565*H_R[39])-0.7036456405748565*H_L[39]+1.407291281149713*H_C[39]+0.5*H_R[23]-0.5*H_L[23])*dv1*gamma_avg; 
  out[24] = (-(0.7036456405748565*H_R[40])-0.7036456405748565*H_L[40]+1.407291281149713*H_C[40]+0.5*H_R[24]-0.5*H_L[24])*dv1*gamma_avg; 
  out[25] = (-(1.21875*H_R[25])+1.21875*H_L[25]+0.8660254037844387*H_R[11]+0.8660254037844387*H_L[11]-1.7320508075688774*H_C[11])*dv1*gamma_avg; 
  out[26] = (-(1.21875*H_R[26])+1.21875*H_L[26]+0.8660254037844387*H_R[12]+0.8660254037844387*H_L[12]-1.7320508075688774*H_C[12])*dv1*gamma_avg; 
  out[27] = (-(1.21875*H_R[27])+1.21875*H_L[27]+0.8660254037844387*H_R[13]+0.8660254037844387*H_L[13]-1.7320508075688774*H_C[13])*dv1*gamma_avg; 
  out[28] = (1.09375*H_R[28]-1.09375*H_L[28]-1.5733994843967631*H_R[8]-1.5733994843967631*H_L[8]-4.599167723621307*H_C[8]+1.118033988749895*H_R[1]-1.118033988749895*H_L[1])*dv1*gamma_avg; 
  out[29] = (1.09375*H_R[29]-1.09375*H_L[29]-1.5733994843967631*H_R[9]-1.5733994843967631*H_L[9]-4.599167723621307*H_C[9]+1.118033988749895*H_R[2]-1.118033988749895*H_L[2])*dv1*gamma_avg; 
  out[30] = (1.09375*H_R[30]-1.09375*H_L[30]-1.5733994843967631*H_R[10]-1.5733994843967631*H_L[10]-4.599167723621307*H_C[10]+1.118033988749895*H_R[3]-1.118033988749895*H_L[3])*dv1*gamma_avg; 
  out[31] = (0.8472151069828725*H_R[47]+0.8472151069828725*H_L[47]+1.694430213965745*H_C[47]-1.21875*H_R[31]+1.21875*H_L[31]+0.8660254037844386*H_R[15]+0.8660254037844386*H_L[15]-1.7320508075688772*H_C[15])*dv1*gamma_avg; 
  out[32] = (-(0.7036456405748565*H_R[44])-0.7036456405748565*H_L[44]+1.407291281149713*H_C[44]+0.5*H_R[32]-0.5*H_L[32])*dv1*gamma_avg; 
  out[33] = (-(0.7036456405748565*H_R[45])-0.7036456405748565*H_L[45]+1.407291281149713*H_C[45]+0.5*H_R[33]-0.5*H_L[33])*dv1*gamma_avg; 
  out[34] = (-(0.7036456405748565*H_R[46])-0.7036456405748565*H_L[46]+1.407291281149713*H_C[46]+0.5*H_R[34]-0.5*H_L[34])*dv1*gamma_avg; 
  out[35] = (-(1.21875*H_R[35])+1.21875*H_L[35]+0.8660254037844387*H_R[19]+0.8660254037844387*H_L[19]-1.7320508075688774*H_C[19])*dv1*gamma_avg; 
  out[36] = (-(1.21875*H_R[36])+1.21875*H_L[36]+0.8660254037844387*H_R[20]+0.8660254037844387*H_L[20]-1.7320508075688774*H_C[20])*dv1*gamma_avg; 
  out[37] = (-(1.21875*H_R[37])+1.21875*H_L[37]+0.8660254037844387*H_R[21]+0.8660254037844387*H_L[21]-1.7320508075688774*H_C[21])*dv1*gamma_avg; 
  out[38] = (-(1.21875*H_R[38])+1.21875*H_L[38]+0.8660254037844387*H_R[22]+0.8660254037844387*H_L[22]-1.7320508075688774*H_C[22])*dv1*gamma_avg; 
  out[39] = (-(1.21875*H_R[39])+1.21875*H_L[39]+0.8660254037844387*H_R[23]+0.8660254037844387*H_L[23]-1.7320508075688774*H_C[23])*dv1*gamma_avg; 
  out[40] = (-(1.21875*H_R[40])+1.21875*H_L[40]+0.8660254037844387*H_R[24]+0.8660254037844387*H_L[24]-1.7320508075688774*H_C[24])*dv1*gamma_avg; 
  out[41] = (1.09375*H_R[41]-1.09375*H_L[41]-1.5733994843967631*H_R[16]-1.5733994843967631*H_L[16]-4.599167723621307*H_C[16]+1.118033988749895*H_R[5]-1.118033988749895*H_L[5])*dv1*gamma_avg; 
  out[42] = (1.09375*H_R[42]-1.09375*H_L[42]-1.5733994843967631*H_R[17]-1.5733994843967631*H_L[17]-4.599167723621307*H_C[17]+1.118033988749895*H_R[6]-1.118033988749895*H_L[6])*dv1*gamma_avg; 
  out[43] = (1.09375*H_R[43]-1.09375*H_L[43]-1.5733994843967631*H_R[18]-1.5733994843967631*H_L[18]-4.599167723621307*H_C[18]+1.118033988749895*H_R[7]-1.118033988749895*H_L[7])*dv1*gamma_avg; 
  out[44] = (-(1.21875*H_R[44])+1.21875*H_L[44]+0.8660254037844387*H_R[32]+0.8660254037844387*H_L[32]-1.7320508075688774*H_C[32])*dv1*gamma_avg; 
  out[45] = (-(1.21875*H_R[45])+1.21875*H_L[45]+0.8660254037844387*H_R[33]+0.8660254037844387*H_L[33]-1.7320508075688774*H_C[33])*dv1*gamma_avg; 
  out[46] = (-(1.21875*H_R[46])+1.21875*H_L[46]+0.8660254037844387*H_R[34]+0.8660254037844387*H_L[34]-1.7320508075688774*H_C[34])*dv1*gamma_avg; 
  out[47] = (1.09375*H_R[47]-1.09375*H_L[47]-1.5733994843967631*H_R[31]-1.5733994843967631*H_L[31]-4.599167723621307*H_C[31]+1.118033988749895*H_R[15]-1.118033988749895*H_L[15])*dv1*gamma_avg; 

  out_surf[0] = -(0.03535533905932736*(53.665631459994955*H_L[14]-53.665631459994955*H_C[14]+95.26279441628824*H_L[4]+95.26279441628824*H_C[4]+75.0*H_L[0]-75.0*H_C[0])*dv1*gamma_avg); 
  out_surf[1] = -(0.03535533905932736*(53.66563145999495*H_L[28]-53.66563145999495*H_C[28]+95.26279441628824*H_L[8]+95.26279441628824*H_C[8]+75.0*H_L[1]-75.0*H_C[1])*dv1*gamma_avg); 
  out_surf[2] = -(0.03535533905932736*(53.66563145999495*H_L[29]-53.66563145999495*H_C[29]+95.26279441628824*H_L[9]+95.26279441628824*H_C[9]+75.0*H_L[2]-75.0*H_C[2])*dv1*gamma_avg); 
  out_surf[3] = -(0.03535533905932736*(53.66563145999495*H_L[30]-53.66563145999495*H_C[30]+95.26279441628824*H_L[10]+95.26279441628824*H_C[10]+75.0*H_L[3]-75.0*H_C[3])*dv1*gamma_avg); 
  out_surf[4] = -(0.03535533905932736*(53.665631459994955*H_L[41]-53.665631459994955*H_C[41]+95.26279441628824*H_L[16]+95.26279441628824*H_C[16]+75.0*H_L[5]-75.0*H_C[5])*dv1*gamma_avg); 
  out_surf[5] = -(0.03535533905932736*(53.665631459994955*H_L[42]-53.665631459994955*H_C[42]+95.26279441628824*H_L[17]+95.26279441628824*H_C[17]+75.0*H_L[6]-75.0*H_C[6])*dv1*gamma_avg); 
  out_surf[6] = -(0.03535533905932736*(53.665631459994955*H_L[43]-53.665631459994955*H_C[43]+95.26279441628824*H_L[18]+95.26279441628824*H_C[18]+75.0*H_L[7]-75.0*H_C[7])*dv1*gamma_avg); 
  out_surf[7] = -(0.03535533905932736*(95.26279441628826*H_L[25]+95.26279441628826*H_C[25]+75.0*H_L[11]-75.0*H_C[11])*dv1*gamma_avg); 
  out_surf[8] = -(0.03535533905932736*(95.26279441628826*H_L[26]+95.26279441628826*H_C[26]+75.0*H_L[12]-75.0*H_C[12])*dv1*gamma_avg); 
  out_surf[9] = -(0.03535533905932736*(95.26279441628826*H_L[27]+95.26279441628826*H_C[27]+75.0*H_L[13]-75.0*H_C[13])*dv1*gamma_avg); 
  out_surf[10] = -(0.03535533905932736*(53.66563145999495*H_L[47]-53.66563145999495*H_C[47]+95.26279441628824*H_L[31]+95.26279441628824*H_C[31]+75.0*H_L[15]-75.0*H_C[15])*dv1*gamma_avg); 
  out_surf[11] = -(0.03535533905932736*(95.26279441628826*H_L[35]+95.26279441628826*H_C[35]+75.0*H_L[19]-75.0*H_C[19])*dv1*gamma_avg); 
  out_surf[12] = -(0.03535533905932736*(95.26279441628826*H_L[36]+95.26279441628826*H_C[36]+75.0*H_L[20]-75.0*H_C[20])*dv1*gamma_avg); 
  out_surf[13] = -(0.03535533905932736*(95.26279441628826*H_L[37]+95.26279441628826*H_C[37]+75.0*H_L[21]-75.0*H_C[21])*dv1*gamma_avg); 
  out_surf[14] = -(0.03535533905932736*(95.26279441628826*H_L[38]+95.26279441628826*H_C[38]+75.0*H_L[22]-75.0*H_C[22])*dv1*gamma_avg); 
  out_surf[15] = -(0.03535533905932736*(95.26279441628826*H_L[39]+95.26279441628826*H_C[39]+75.0*H_L[23]-75.0*H_C[23])*dv1*gamma_avg); 
  out_surf[16] = -(0.03535533905932736*(95.26279441628826*H_L[40]+95.26279441628826*H_C[40]+75.0*H_L[24]-75.0*H_C[24])*dv1*gamma_avg); 
  out_surf[17] = -(0.03535533905932736*(95.26279441628826*H_L[44]+95.26279441628826*H_C[44]+75.0*H_L[32]-75.0*H_C[32])*dv1*gamma_avg); 
  out_surf[18] = -(0.03535533905932736*(95.26279441628826*H_L[45]+95.26279441628826*H_C[45]+75.0*H_L[33]-75.0*H_C[33])*dv1*gamma_avg); 
  out_surf[19] = -(0.03535533905932736*(95.26279441628826*H_L[46]+95.26279441628826*H_C[46]+75.0*H_L[34]-75.0*H_C[34])*dv1*gamma_avg); 

  int const_sgn_alpha_surf = 1;  
  
  if (0.5692099788303081*(out_surf[19]+out_surf[18]+out_surf[17])-0.42426406871192823*(out_surf[16]+out_surf[15])-0.42426406871192857*(out_surf[14]+out_surf[13])-0.42426406871192823*out_surf[12]-0.42426406871192857*out_surf[11]-0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])+0.6363961030678927*(out_surf[6]+out_surf[5]+out_surf[4])-0.4743416490252568*(out_surf[3]+out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (-(0.7115124735378848*out_surf[19])+0.5303300858899102*(out_surf[16]+out_surf[15])-0.42426406871192823*out_surf[12]-0.42426406871192857*out_surf[11]-0.3952847075210471*out_surf[9]+0.31622776601683783*(out_surf[8]+out_surf[7])+0.6363961030678927*out_surf[4]-0.4743416490252568*(out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303081*out_surf[19]-0.5692099788303081*(out_surf[18]+out_surf[17])-0.42426406871192823*(out_surf[16]+out_surf[15])+0.42426406871192857*(out_surf[14]+out_surf[13])-0.42426406871192823*out_surf[12]-0.42426406871192857*out_surf[11]+0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])-0.6363961030678927*(out_surf[6]+out_surf[5])+0.6363961030678927*out_surf[4]+0.4743416490252568*out_surf[3]-0.4743416490252568*(out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*out_surf[18])-0.42426406871192823*out_surf[15]+0.5303300858899102*out_surf[14]-0.42426406871192857*out_surf[13]+0.5303300858899102*out_surf[12]+0.31622776601683783*out_surf[9]-0.3952847075210471*out_surf[8]+0.31622776601683783*out_surf[7]+0.6363961030678927*out_surf[5]-0.4743416490252568*(out_surf[3]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(out_surf[15]+out_surf[12])-0.3952847075210471*(out_surf[9]+out_surf[8])+0.31622776601683783*out_surf[7]-0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*out_surf[18]-0.42426406871192823*out_surf[15]-0.5303300858899102*out_surf[14]+0.42426406871192857*out_surf[13]+0.5303300858899102*out_surf[12]+0.31622776601683783*out_surf[9]-0.3952847075210471*out_surf[8]+0.31622776601683783*out_surf[7]-0.6363961030678927*out_surf[5]+0.4743416490252568*out_surf[3]-0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*out_surf[19])+0.5692099788303081*out_surf[18]-0.5692099788303081*out_surf[17]+0.42426406871192823*out_surf[16]-0.42426406871192823*out_surf[15]-0.42426406871192857*(out_surf[14]+out_surf[13])-0.42426406871192823*out_surf[12]+0.42426406871192857*out_surf[11]+0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])-0.6363961030678927*out_surf[6]+0.6363961030678927*out_surf[5]-0.6363961030678927*out_surf[4]-0.4743416490252568*out_surf[3]+0.4743416490252568*out_surf[2]-0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*out_surf[19]-0.5303300858899102*out_surf[16]+0.5303300858899102*out_surf[15]-0.42426406871192823*out_surf[12]+0.42426406871192857*out_surf[11]-0.3952847075210471*out_surf[9]+0.31622776601683783*(out_surf[8]+out_surf[7])-0.6363961030678927*out_surf[4]+0.4743416490252568*out_surf[2]-0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*(out_surf[19]+out_surf[18]))+0.5692099788303081*out_surf[17]+0.42426406871192823*out_surf[16]-0.42426406871192823*out_surf[15]+0.42426406871192857*(out_surf[14]+out_surf[13])-0.42426406871192823*out_surf[12]+0.42426406871192857*out_surf[11]-0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])+0.6363961030678927*out_surf[6]-0.6363961030678927*(out_surf[5]+out_surf[4])+0.4743416490252568*(out_surf[3]+out_surf[2])-0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[8] = 1.0; 
  else  
    sgn_alpha_surf[8] = -1.0; 
  
  if (sgn_alpha_surf[8] == sgn_alpha_surf[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*out_surf[17])-0.42426406871192823*out_surf[16]-0.42426406871192857*out_surf[14]+0.5303300858899102*(out_surf[13]+out_surf[11])+0.31622776601683783*(out_surf[9]+out_surf[8])-0.3952847075210471*out_surf[7]+0.6363961030678927*out_surf[6]-0.4743416490252568*(out_surf[3]+out_surf[2])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[9] = 1.0; 
  else  
    sgn_alpha_surf[9] = -1.0; 
  
  if (sgn_alpha_surf[9] == sgn_alpha_surf[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(out_surf[16]+out_surf[11])-0.3952847075210471*out_surf[9]+0.31622776601683783*out_surf[8]-0.3952847075210471*out_surf[7]-0.4743416490252568*out_surf[2]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[10] = 1.0; 
  else  
    sgn_alpha_surf[10] = -1.0; 
  
  if (sgn_alpha_surf[10] == sgn_alpha_surf[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*out_surf[17]-0.42426406871192823*out_surf[16]+0.42426406871192857*out_surf[14]-0.5303300858899102*out_surf[13]+0.5303300858899102*out_surf[11]+0.31622776601683783*(out_surf[9]+out_surf[8])-0.3952847075210471*out_surf[7]-0.6363961030678927*out_surf[6]+0.4743416490252568*out_surf[3]-0.4743416490252568*out_surf[2]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[11] = 1.0; 
  else  
    sgn_alpha_surf[11] = -1.0; 
  
  if (sgn_alpha_surf[11] == sgn_alpha_surf[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*(out_surf[14]+out_surf[13])+0.31622776601683783*out_surf[9]-0.3952847075210471*(out_surf[8]+out_surf[7])-0.4743416490252568*out_surf[3]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[12] = 1.0; 
  else  
    sgn_alpha_surf[12] = -1.0; 
  
  if (sgn_alpha_surf[12] == sgn_alpha_surf[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*out_surf[0]-0.3952847075210471*(out_surf[9]+out_surf[8]+out_surf[7]) > 0.) 
    sgn_alpha_surf[13] = 1.0; 
  else  
    sgn_alpha_surf[13] = -1.0; 
  
  if (sgn_alpha_surf[13] == sgn_alpha_surf[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(out_surf[14]+out_surf[13]))+0.31622776601683783*out_surf[9]-0.3952847075210471*(out_surf[8]+out_surf[7])+0.4743416490252568*out_surf[3]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[14] = 1.0; 
  else  
    sgn_alpha_surf[14] = -1.0; 
  
  if (sgn_alpha_surf[14] == sgn_alpha_surf[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*out_surf[17]+0.42426406871192823*out_surf[16]-0.42426406871192857*out_surf[14]+0.5303300858899102*out_surf[13]-0.5303300858899102*out_surf[11]+0.31622776601683783*(out_surf[9]+out_surf[8])-0.3952847075210471*out_surf[7]-0.6363961030678927*out_surf[6]-0.4743416490252568*out_surf[3]+0.4743416490252568*out_surf[2]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[15] = 1.0; 
  else  
    sgn_alpha_surf[15] = -1.0; 
  
  if (sgn_alpha_surf[15] == sgn_alpha_surf[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(out_surf[16]+out_surf[11]))-0.3952847075210471*out_surf[9]+0.31622776601683783*out_surf[8]-0.3952847075210471*out_surf[7]+0.4743416490252568*out_surf[2]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[16] = 1.0; 
  else  
    sgn_alpha_surf[16] = -1.0; 
  
  if (sgn_alpha_surf[16] == sgn_alpha_surf[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*out_surf[17])+0.42426406871192823*out_surf[16]+0.42426406871192857*out_surf[14]-0.5303300858899102*(out_surf[13]+out_surf[11])+0.31622776601683783*(out_surf[9]+out_surf[8])-0.3952847075210471*out_surf[7]+0.6363961030678927*out_surf[6]+0.4743416490252568*(out_surf[3]+out_surf[2])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[17] = 1.0; 
  else  
    sgn_alpha_surf[17] = -1.0; 
  
  if (sgn_alpha_surf[17] == sgn_alpha_surf[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*(out_surf[19]+out_surf[18]))+0.5692099788303081*out_surf[17]-0.42426406871192823*out_surf[16]+0.42426406871192823*out_surf[15]-0.42426406871192857*(out_surf[14]+out_surf[13])+0.42426406871192823*out_surf[12]-0.42426406871192857*out_surf[11]+0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])+0.6363961030678927*out_surf[6]-0.6363961030678927*(out_surf[5]+out_surf[4])-0.4743416490252568*(out_surf[3]+out_surf[2])+0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[18] = 1.0; 
  else  
    sgn_alpha_surf[18] = -1.0; 
  
  if (sgn_alpha_surf[18] == sgn_alpha_surf[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*out_surf[19]+0.5303300858899102*out_surf[16]-0.5303300858899102*out_surf[15]+0.42426406871192823*out_surf[12]-0.42426406871192857*out_surf[11]-0.3952847075210471*out_surf[9]+0.31622776601683783*(out_surf[8]+out_surf[7])-0.6363961030678927*out_surf[4]-0.4743416490252568*out_surf[2]+0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[19] = 1.0; 
  else  
    sgn_alpha_surf[19] = -1.0; 
  
  if (sgn_alpha_surf[19] == sgn_alpha_surf[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5692099788303081*out_surf[19])+0.5692099788303081*out_surf[18]-0.5692099788303081*out_surf[17]-0.42426406871192823*out_surf[16]+0.42426406871192823*out_surf[15]+0.42426406871192857*(out_surf[14]+out_surf[13])+0.42426406871192823*out_surf[12]-0.42426406871192857*out_surf[11]-0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])-0.6363961030678927*out_surf[6]+0.6363961030678927*out_surf[5]-0.6363961030678927*out_surf[4]+0.4743416490252568*out_surf[3]-0.4743416490252568*out_surf[2]+0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[20] = 1.0; 
  else  
    sgn_alpha_surf[20] = -1.0; 
  
  if (sgn_alpha_surf[20] == sgn_alpha_surf[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.7115124735378848*out_surf[18]+0.42426406871192823*out_surf[15]+0.5303300858899102*out_surf[14]-0.42426406871192857*out_surf[13]-0.5303300858899102*out_surf[12]+0.31622776601683783*out_surf[9]-0.3952847075210471*out_surf[8]+0.31622776601683783*out_surf[7]-0.6363961030678927*out_surf[5]-0.4743416490252568*out_surf[3]+0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[21] = 1.0; 
  else  
    sgn_alpha_surf[21] = -1.0; 
  
  if (sgn_alpha_surf[21] == sgn_alpha_surf[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.5303300858899102*(out_surf[15]+out_surf[12]))-0.3952847075210471*(out_surf[9]+out_surf[8])+0.31622776601683783*out_surf[7]+0.4743416490252568*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[22] = 1.0; 
  else  
    sgn_alpha_surf[22] = -1.0; 
  
  if (sgn_alpha_surf[22] == sgn_alpha_surf[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*out_surf[18])+0.42426406871192823*out_surf[15]-0.5303300858899102*out_surf[14]+0.42426406871192857*out_surf[13]-0.5303300858899102*out_surf[12]+0.31622776601683783*out_surf[9]-0.3952847075210471*out_surf[8]+0.31622776601683783*out_surf[7]+0.6363961030678927*out_surf[5]+0.4743416490252568*(out_surf[3]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[23] = 1.0; 
  else  
    sgn_alpha_surf[23] = -1.0; 
  
  if (sgn_alpha_surf[23] == sgn_alpha_surf[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303081*out_surf[19]-0.5692099788303081*(out_surf[18]+out_surf[17])+0.42426406871192823*(out_surf[16]+out_surf[15])-0.42426406871192857*(out_surf[14]+out_surf[13])+0.42426406871192823*out_surf[12]+0.42426406871192857*out_surf[11]-0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])-0.6363961030678927*(out_surf[6]+out_surf[5])+0.6363961030678927*out_surf[4]-0.4743416490252568*out_surf[3]+0.4743416490252568*(out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[24] = 1.0; 
  else  
    sgn_alpha_surf[24] = -1.0; 
  
  if (sgn_alpha_surf[24] == sgn_alpha_surf[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.7115124735378848*out_surf[19])-0.5303300858899102*(out_surf[16]+out_surf[15])+0.42426406871192823*out_surf[12]+0.42426406871192857*out_surf[11]-0.3952847075210471*out_surf[9]+0.31622776601683783*(out_surf[8]+out_surf[7])+0.6363961030678927*out_surf[4]+0.4743416490252568*(out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[25] = 1.0; 
  else  
    sgn_alpha_surf[25] = -1.0; 
  
  if (sgn_alpha_surf[25] == sgn_alpha_surf[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5692099788303081*(out_surf[19]+out_surf[18]+out_surf[17])+0.42426406871192823*(out_surf[16]+out_surf[15])+0.42426406871192857*(out_surf[14]+out_surf[13])+0.42426406871192823*out_surf[12]+0.42426406871192857*out_surf[11]+0.8538149682454614*out_surf[10]+0.31622776601683783*(out_surf[9]+out_surf[8]+out_surf[7])+0.6363961030678927*(out_surf[6]+out_surf[5]+out_surf[4])+0.4743416490252568*(out_surf[3]+out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[26] = 1.0; 
  else  
    sgn_alpha_surf[26] = -1.0; 
  
  if (sgn_alpha_surf[26] == sgn_alpha_surf[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 
} 

