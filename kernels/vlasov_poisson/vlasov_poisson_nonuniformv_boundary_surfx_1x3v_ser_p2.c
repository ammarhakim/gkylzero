#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_boundary_surfx_1x3v_ser_p2(const double *w, const double *dxv, const double *vcoord, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // vcoord:      Discrete (DG) velocity coordinate.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double rdx2 = 2/dxv[0]; 
  double Ghat[20]; 

  if (edge == -1) { 

  if (0.3535533905932737*vcoord[0]-0.3952847075210473*vcoord[7]>0) { 

  Ghat[0] = 0.4330127018922193*vcoord[7]*fSkin[20]+0.5590169943749475*vcoord[1]*fSkin[19]+0.25*vcoord[7]*fSkin[12]+0.5590169943749475*vcoord[0]*fSkin[11]+vcoord[1]*(0.4330127018922193*fSkin[5]+0.25*fSkin[2])+vcoord[0]*(0.4330127018922193*fSkin[1]+0.25*fSkin[0]); 
  Ghat[1] = 0.3872983346207416*vcoord[1]*fSkin[20]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fSkin[19]+vcoord[1]*(0.223606797749979*fSkin[12]+0.5590169943749475*fSkin[11])+(0.3872983346207416*fSkin[5]+0.223606797749979*fSkin[2])*vcoord[7]+vcoord[0]*(0.4330127018922193*fSkin[5]+0.25*fSkin[2])+(0.4330127018922193*fSkin[1]+0.25*fSkin[0])*vcoord[1]; 
  Ghat[2] = 0.4330127018922193*vcoord[7]*fSkin[33]+0.5590169943749475*vcoord[1]*fSkin[32]+0.25*vcoord[7]*fSkin[22]+0.5590169943749475*vcoord[0]*fSkin[21]+vcoord[1]*(0.4330127018922193*fSkin[15]+0.25*fSkin[7])+vcoord[0]*(0.4330127018922193*fSkin[6]+0.25*fSkin[3]); 
  Ghat[3] = 0.4330127018922193*vcoord[7]*fSkin[36]+0.5590169943749475*vcoord[1]*fSkin[35]+0.25*vcoord[7]*fSkin[26]+0.5590169943749475*vcoord[0]*fSkin[25]+vcoord[1]*(0.4330127018922193*fSkin[16]+0.25*fSkin[9])+vcoord[0]*(0.4330127018922193*fSkin[8]+0.25*fSkin[4]); 
  Ghat[4] = 0.3872983346207416*vcoord[1]*fSkin[33]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fSkin[32]+vcoord[1]*(0.223606797749979*fSkin[22]+0.5590169943749475*fSkin[21])+(0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[15]+fSkin[7]*(0.223606797749979*vcoord[7]+0.25*vcoord[0])+vcoord[1]*(0.4330127018922193*fSkin[6]+0.25*fSkin[3]); 
  Ghat[5] = 0.3872983346207416*vcoord[1]*fSkin[36]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fSkin[35]+vcoord[1]*(0.223606797749979*fSkin[26]+0.5590169943749475*fSkin[25])+(0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[16]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fSkin[9]+vcoord[1]*(0.4330127018922193*fSkin[8]+0.25*fSkin[4]); 
  Ghat[6] = 0.4330127018922193*vcoord[7]*fSkin[45]+0.5590169943749475*vcoord[1]*fSkin[44]+0.25*vcoord[7]*fSkin[38]+0.5590169943749475*vcoord[0]*fSkin[37]+vcoord[1]*(0.4330127018922193*fSkin[31]+0.25*fSkin[18])+vcoord[0]*(0.4330127018922193*fSkin[17]+0.25*fSkin[10]); 
  Ghat[7] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[20]+0.5*vcoord[1]*fSkin[19]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fSkin[12]+vcoord[7]*(0.5590169943749475*fSkin[11]+0.4330127018922193*fSkin[1]+0.25*fSkin[0])+vcoord[1]*(0.3872983346207416*fSkin[5]+0.223606797749979*fSkin[2]); 
  Ghat[8] = vcoord[1]*(0.4330127018922193*fSkin[34]+0.25*fSkin[24])+vcoord[0]*(0.4330127018922193*fSkin[23]+0.25*fSkin[13]); 
  Ghat[9] = vcoord[1]*(0.4330127018922193*fSkin[41]+0.25*fSkin[29])+vcoord[0]*(0.4330127018922193*fSkin[28]+0.25*fSkin[14]); 
  Ghat[10] = 0.3872983346207416*vcoord[1]*fSkin[45]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fSkin[44]+vcoord[1]*(0.223606797749979*fSkin[38]+0.5590169943749475*fSkin[37])+(0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[31]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fSkin[18]+vcoord[1]*(0.4330127018922193*fSkin[17]+0.25*fSkin[10]); 
  Ghat[11] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[33]+0.5*vcoord[1]*fSkin[32]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fSkin[22]+0.5590169943749475*vcoord[7]*fSkin[21]+0.3872983346207416*vcoord[1]*fSkin[15]+(0.4330127018922193*fSkin[6]+0.25*fSkin[3])*vcoord[7]+0.223606797749979*vcoord[1]*fSkin[7]; 
  Ghat[12] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[34]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fSkin[24]+vcoord[1]*(0.4330127018922193*fSkin[23]+0.25*fSkin[13]); 
  Ghat[13] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[36]+0.5*vcoord[1]*fSkin[35]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fSkin[26]+0.5590169943749475*vcoord[7]*fSkin[25]+vcoord[1]*(0.3872983346207416*fSkin[16]+0.223606797749979*fSkin[9])+vcoord[7]*(0.4330127018922193*fSkin[8]+0.25*fSkin[4]); 
  Ghat[14] = vcoord[1]*(0.4330127018922193*fSkin[46]+0.25*fSkin[40])+vcoord[0]*(0.4330127018922193*fSkin[39]+0.25*fSkin[27]); 
  Ghat[15] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[41]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fSkin[29]+vcoord[1]*(0.4330127018922193*fSkin[28]+0.25*fSkin[14]); 
  Ghat[16] = vcoord[1]*(0.4330127018922193*fSkin[47]+0.25*fSkin[43])+vcoord[0]*(0.4330127018922193*fSkin[42]+0.25*fSkin[30]); 
  Ghat[17] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[45]+0.5*vcoord[1]*fSkin[44]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fSkin[38]+0.5590169943749475*vcoord[7]*fSkin[37]+vcoord[1]*(0.3872983346207416*fSkin[31]+0.223606797749979*fSkin[18])+vcoord[7]*(0.4330127018922193*fSkin[17]+0.25*fSkin[10]); 
  Ghat[18] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[46]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fSkin[40]+vcoord[1]*(0.4330127018922193*fSkin[39]+0.25*fSkin[27]); 
  Ghat[19] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fSkin[47]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fSkin[43]+vcoord[1]*(0.4330127018922193*fSkin[42]+0.25*fSkin[30]); 

  } else { 

  Ghat[0] = (-0.4330127018922194*vcoord[7]*fEdge[20])+0.5590169943749476*vcoord[1]*fEdge[19]+0.25*vcoord[7]*fEdge[12]+0.5590169943749475*vcoord[0]*fEdge[11]-0.4330127018922193*vcoord[1]*fEdge[5]+0.25*vcoord[1]*fEdge[2]-0.4330127018922193*vcoord[0]*fEdge[1]+0.25*fEdge[0]*vcoord[0]; 
  Ghat[1] = (-0.3872983346207417*vcoord[1]*fEdge[20])+0.5000000000000001*vcoord[7]*fEdge[19]+0.5590169943749476*vcoord[0]*fEdge[19]+0.223606797749979*vcoord[1]*fEdge[12]+0.5590169943749475*vcoord[1]*fEdge[11]-0.3872983346207416*fEdge[5]*vcoord[7]+0.223606797749979*fEdge[2]*vcoord[7]-0.4330127018922193*vcoord[0]*fEdge[5]+0.25*vcoord[0]*fEdge[2]-0.4330127018922193*fEdge[1]*vcoord[1]+0.25*fEdge[0]*vcoord[1]; 
  Ghat[2] = (-0.4330127018922193*vcoord[7]*fEdge[33])+0.5590169943749475*vcoord[1]*fEdge[32]+0.2500000000000001*vcoord[7]*fEdge[22]+0.5590169943749476*vcoord[0]*fEdge[21]-0.4330127018922193*vcoord[1]*fEdge[15]+0.25*vcoord[1]*fEdge[7]-0.4330127018922193*vcoord[0]*fEdge[6]+0.25*vcoord[0]*fEdge[3]; 
  Ghat[3] = (-0.4330127018922193*vcoord[7]*fEdge[36])+0.5590169943749475*vcoord[1]*fEdge[35]+0.2500000000000001*vcoord[7]*fEdge[26]+0.5590169943749476*vcoord[0]*fEdge[25]-0.4330127018922193*vcoord[1]*fEdge[16]+0.25*vcoord[1]*fEdge[9]-0.4330127018922193*vcoord[0]*fEdge[8]+0.25*vcoord[0]*fEdge[4]; 
  Ghat[4] = (-0.3872983346207416*vcoord[1]*fEdge[33])+0.5*vcoord[7]*fEdge[32]+0.5590169943749475*vcoord[0]*fEdge[32]+0.223606797749979*vcoord[1]*fEdge[22]+0.5590169943749476*vcoord[1]*fEdge[21]-0.3872983346207416*vcoord[7]*fEdge[15]-0.4330127018922193*vcoord[0]*fEdge[15]+0.223606797749979*fEdge[7]*vcoord[7]+0.25*vcoord[0]*fEdge[7]-0.4330127018922193*vcoord[1]*fEdge[6]+0.25*vcoord[1]*fEdge[3]; 
  Ghat[5] = (-0.3872983346207416*vcoord[1]*fEdge[36])+0.5*vcoord[7]*fEdge[35]+0.5590169943749475*vcoord[0]*fEdge[35]+0.223606797749979*vcoord[1]*fEdge[26]+0.5590169943749476*vcoord[1]*fEdge[25]-0.3872983346207416*vcoord[7]*fEdge[16]-0.4330127018922193*vcoord[0]*fEdge[16]+0.223606797749979*vcoord[7]*fEdge[9]+0.25*vcoord[0]*fEdge[9]-0.4330127018922193*vcoord[1]*fEdge[8]+0.25*vcoord[1]*fEdge[4]; 
  Ghat[6] = (-0.4330127018922194*vcoord[7]*fEdge[45])+0.5590169943749476*vcoord[1]*fEdge[44]+0.25*vcoord[7]*fEdge[38]+0.5590169943749475*vcoord[0]*fEdge[37]-0.4330127018922193*vcoord[1]*fEdge[31]+0.25*vcoord[1]*fEdge[18]-0.4330127018922193*vcoord[0]*fEdge[17]+0.25*vcoord[0]*fEdge[10]; 
  Ghat[7] = (-0.276641667586244*vcoord[7]*fEdge[20])-0.4330127018922194*vcoord[0]*fEdge[20]+0.5000000000000001*vcoord[1]*fEdge[19]+0.159719141249985*vcoord[7]*fEdge[12]+0.25*vcoord[0]*fEdge[12]+0.5590169943749475*vcoord[7]*fEdge[11]-0.4330127018922193*fEdge[1]*vcoord[7]+0.25*fEdge[0]*vcoord[7]-0.3872983346207416*vcoord[1]*fEdge[5]+0.223606797749979*vcoord[1]*fEdge[2]; 
  Ghat[8] = (-0.4330127018922193*vcoord[1]*fEdge[34])+0.2500000000000001*vcoord[1]*fEdge[24]-0.4330127018922194*vcoord[0]*fEdge[23]+0.25*vcoord[0]*fEdge[13]; 
  Ghat[9] = (-0.4330127018922193*vcoord[1]*fEdge[41])+0.2500000000000001*vcoord[1]*fEdge[29]-0.4330127018922194*vcoord[0]*fEdge[28]+0.25*vcoord[0]*fEdge[14]; 
  Ghat[10] = (-0.3872983346207417*vcoord[1]*fEdge[45])+0.5000000000000001*vcoord[7]*fEdge[44]+0.5590169943749476*vcoord[0]*fEdge[44]+0.223606797749979*vcoord[1]*fEdge[38]+0.5590169943749475*vcoord[1]*fEdge[37]-0.3872983346207416*vcoord[7]*fEdge[31]-0.4330127018922193*vcoord[0]*fEdge[31]+0.223606797749979*vcoord[7]*fEdge[18]+0.25*vcoord[0]*fEdge[18]-0.4330127018922193*vcoord[1]*fEdge[17]+0.25*vcoord[1]*fEdge[10]; 
  Ghat[11] = (-0.276641667586244*vcoord[7]*fEdge[33])-0.4330127018922194*vcoord[0]*fEdge[33]+0.5000000000000001*vcoord[1]*fEdge[32]+0.159719141249985*vcoord[7]*fEdge[22]+0.25*vcoord[0]*fEdge[22]+0.5590169943749475*vcoord[7]*fEdge[21]-0.3872983346207417*vcoord[1]*fEdge[15]-0.4330127018922194*fEdge[6]*vcoord[7]+0.2500000000000001*fEdge[3]*vcoord[7]+0.223606797749979*vcoord[1]*fEdge[7]; 
  Ghat[12] = (-0.3872983346207417*vcoord[7]*fEdge[34])-0.4330127018922194*vcoord[0]*fEdge[34]+0.223606797749979*vcoord[7]*fEdge[24]+0.25*vcoord[0]*fEdge[24]-0.4330127018922193*vcoord[1]*fEdge[23]+0.2500000000000001*vcoord[1]*fEdge[13]; 
  Ghat[13] = (-0.276641667586244*vcoord[7]*fEdge[36])-0.4330127018922194*vcoord[0]*fEdge[36]+0.5000000000000001*vcoord[1]*fEdge[35]+0.159719141249985*vcoord[7]*fEdge[26]+0.25*vcoord[0]*fEdge[26]+0.5590169943749475*vcoord[7]*fEdge[25]-0.3872983346207417*vcoord[1]*fEdge[16]+0.223606797749979*vcoord[1]*fEdge[9]-0.4330127018922194*vcoord[7]*fEdge[8]+0.2500000000000001*fEdge[4]*vcoord[7]; 
  Ghat[14] = (-0.4330127018922193*vcoord[1]*fEdge[46])+0.2500000000000001*vcoord[1]*fEdge[40]-0.4330127018922194*vcoord[0]*fEdge[39]+0.25*vcoord[0]*fEdge[27]; 
  Ghat[15] = (-0.3872983346207417*vcoord[7]*fEdge[41])-0.4330127018922194*vcoord[0]*fEdge[41]+0.223606797749979*vcoord[7]*fEdge[29]+0.25*vcoord[0]*fEdge[29]-0.4330127018922193*vcoord[1]*fEdge[28]+0.2500000000000001*vcoord[1]*fEdge[14]; 
  Ghat[16] = (-0.4330127018922193*vcoord[1]*fEdge[47])+0.2500000000000001*vcoord[1]*fEdge[43]-0.4330127018922194*vcoord[0]*fEdge[42]+0.25*vcoord[0]*fEdge[30]; 
  Ghat[17] = (-0.276641667586244*vcoord[7]*fEdge[45])-0.4330127018922194*vcoord[0]*fEdge[45]+0.5000000000000001*vcoord[1]*fEdge[44]+0.159719141249985*vcoord[7]*fEdge[38]+0.25*vcoord[0]*fEdge[38]+0.5590169943749475*vcoord[7]*fEdge[37]-0.3872983346207416*vcoord[1]*fEdge[31]+0.223606797749979*vcoord[1]*fEdge[18]-0.4330127018922193*vcoord[7]*fEdge[17]+0.25*vcoord[7]*fEdge[10]; 
  Ghat[18] = (-0.3872983346207417*vcoord[7]*fEdge[46])-0.4330127018922194*vcoord[0]*fEdge[46]+0.223606797749979*vcoord[7]*fEdge[40]+0.25*vcoord[0]*fEdge[40]-0.4330127018922193*vcoord[1]*fEdge[39]+0.2500000000000001*vcoord[1]*fEdge[27]; 
  Ghat[19] = (-0.3872983346207417*vcoord[7]*fEdge[47])-0.4330127018922194*vcoord[0]*fEdge[47]+0.223606797749979*vcoord[7]*fEdge[43]+0.25*vcoord[0]*fEdge[43]-0.4330127018922193*vcoord[1]*fEdge[42]+0.2500000000000001*vcoord[1]*fEdge[30]; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*rdx2; 
  out[1] += -1.224744871391589*Ghat[0]*rdx2; 
  out[2] += -0.7071067811865475*Ghat[1]*rdx2; 
  out[3] += -0.7071067811865475*Ghat[2]*rdx2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdx2; 
  out[5] += -1.224744871391589*Ghat[1]*rdx2; 
  out[6] += -1.224744871391589*Ghat[2]*rdx2; 
  out[7] += -0.7071067811865475*Ghat[4]*rdx2; 
  out[8] += -1.224744871391589*Ghat[3]*rdx2; 
  out[9] += -0.7071067811865475*Ghat[5]*rdx2; 
  out[10] += -0.7071067811865475*Ghat[6]*rdx2; 
  out[11] += -1.58113883008419*Ghat[0]*rdx2; 
  out[12] += -0.7071067811865475*Ghat[7]*rdx2; 
  out[13] += -0.7071067811865475*Ghat[8]*rdx2; 
  out[14] += -0.7071067811865475*Ghat[9]*rdx2; 
  out[15] += -1.224744871391589*Ghat[4]*rdx2; 
  out[16] += -1.224744871391589*Ghat[5]*rdx2; 
  out[17] += -1.224744871391589*Ghat[6]*rdx2; 
  out[18] += -0.7071067811865475*Ghat[10]*rdx2; 
  out[19] += -1.58113883008419*Ghat[1]*rdx2; 
  out[20] += -1.224744871391589*Ghat[7]*rdx2; 
  out[21] += -1.58113883008419*Ghat[2]*rdx2; 
  out[22] += -0.7071067811865475*Ghat[11]*rdx2; 
  out[23] += -1.224744871391589*Ghat[8]*rdx2; 
  out[24] += -0.7071067811865475*Ghat[12]*rdx2; 
  out[25] += -1.58113883008419*Ghat[3]*rdx2; 
  out[26] += -0.7071067811865475*Ghat[13]*rdx2; 
  out[27] += -0.7071067811865475*Ghat[14]*rdx2; 
  out[28] += -1.224744871391589*Ghat[9]*rdx2; 
  out[29] += -0.7071067811865475*Ghat[15]*rdx2; 
  out[30] += -0.7071067811865475*Ghat[16]*rdx2; 
  out[31] += -1.224744871391589*Ghat[10]*rdx2; 
  out[32] += -1.58113883008419*Ghat[4]*rdx2; 
  out[33] += -1.224744871391589*Ghat[11]*rdx2; 
  out[34] += -1.224744871391589*Ghat[12]*rdx2; 
  out[35] += -1.58113883008419*Ghat[5]*rdx2; 
  out[36] += -1.224744871391589*Ghat[13]*rdx2; 
  out[37] += -1.58113883008419*Ghat[6]*rdx2; 
  out[38] += -0.7071067811865475*Ghat[17]*rdx2; 
  out[39] += -1.224744871391589*Ghat[14]*rdx2; 
  out[40] += -0.7071067811865475*Ghat[18]*rdx2; 
  out[41] += -1.224744871391589*Ghat[15]*rdx2; 
  out[42] += -1.224744871391589*Ghat[16]*rdx2; 
  out[43] += -0.7071067811865475*Ghat[19]*rdx2; 
  out[44] += -1.58113883008419*Ghat[10]*rdx2; 
  out[45] += -1.224744871391589*Ghat[17]*rdx2; 
  out[46] += -1.224744871391589*Ghat[18]*rdx2; 
  out[47] += -1.224744871391589*Ghat[19]*rdx2; 

  } else { 

  if (0.3535533905932737*vcoord[0]-0.3952847075210473*vcoord[7]>0) { 

  Ghat[0] = 0.4330127018922193*vcoord[7]*fEdge[20]+0.5590169943749475*vcoord[1]*fEdge[19]+0.25*vcoord[7]*fEdge[12]+0.5590169943749475*vcoord[0]*fEdge[11]+vcoord[1]*(0.4330127018922193*fEdge[5]+0.25*fEdge[2])+vcoord[0]*(0.4330127018922193*fEdge[1]+0.25*fEdge[0]); 
  Ghat[1] = 0.3872983346207416*vcoord[1]*fEdge[20]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fEdge[19]+vcoord[1]*(0.223606797749979*fEdge[12]+0.5590169943749475*fEdge[11])+(0.3872983346207416*fEdge[5]+0.223606797749979*fEdge[2])*vcoord[7]+vcoord[0]*(0.4330127018922193*fEdge[5]+0.25*fEdge[2])+(0.4330127018922193*fEdge[1]+0.25*fEdge[0])*vcoord[1]; 
  Ghat[2] = 0.4330127018922193*vcoord[7]*fEdge[33]+0.5590169943749475*vcoord[1]*fEdge[32]+0.25*vcoord[7]*fEdge[22]+0.5590169943749475*vcoord[0]*fEdge[21]+vcoord[1]*(0.4330127018922193*fEdge[15]+0.25*fEdge[7])+vcoord[0]*(0.4330127018922193*fEdge[6]+0.25*fEdge[3]); 
  Ghat[3] = 0.4330127018922193*vcoord[7]*fEdge[36]+0.5590169943749475*vcoord[1]*fEdge[35]+0.25*vcoord[7]*fEdge[26]+0.5590169943749475*vcoord[0]*fEdge[25]+vcoord[1]*(0.4330127018922193*fEdge[16]+0.25*fEdge[9])+vcoord[0]*(0.4330127018922193*fEdge[8]+0.25*fEdge[4]); 
  Ghat[4] = 0.3872983346207416*vcoord[1]*fEdge[33]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fEdge[32]+vcoord[1]*(0.223606797749979*fEdge[22]+0.5590169943749475*fEdge[21])+(0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[15]+fEdge[7]*(0.223606797749979*vcoord[7]+0.25*vcoord[0])+vcoord[1]*(0.4330127018922193*fEdge[6]+0.25*fEdge[3]); 
  Ghat[5] = 0.3872983346207416*vcoord[1]*fEdge[36]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fEdge[35]+vcoord[1]*(0.223606797749979*fEdge[26]+0.5590169943749475*fEdge[25])+(0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[16]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fEdge[9]+vcoord[1]*(0.4330127018922193*fEdge[8]+0.25*fEdge[4]); 
  Ghat[6] = 0.4330127018922193*vcoord[7]*fEdge[45]+0.5590169943749475*vcoord[1]*fEdge[44]+0.25*vcoord[7]*fEdge[38]+0.5590169943749475*vcoord[0]*fEdge[37]+vcoord[1]*(0.4330127018922193*fEdge[31]+0.25*fEdge[18])+vcoord[0]*(0.4330127018922193*fEdge[17]+0.25*fEdge[10]); 
  Ghat[7] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[20]+0.5*vcoord[1]*fEdge[19]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fEdge[12]+vcoord[7]*(0.5590169943749475*fEdge[11]+0.4330127018922193*fEdge[1]+0.25*fEdge[0])+vcoord[1]*(0.3872983346207416*fEdge[5]+0.223606797749979*fEdge[2]); 
  Ghat[8] = vcoord[1]*(0.4330127018922193*fEdge[34]+0.25*fEdge[24])+vcoord[0]*(0.4330127018922193*fEdge[23]+0.25*fEdge[13]); 
  Ghat[9] = vcoord[1]*(0.4330127018922193*fEdge[41]+0.25*fEdge[29])+vcoord[0]*(0.4330127018922193*fEdge[28]+0.25*fEdge[14]); 
  Ghat[10] = 0.3872983346207416*vcoord[1]*fEdge[45]+(0.5*vcoord[7]+0.5590169943749475*vcoord[0])*fEdge[44]+vcoord[1]*(0.223606797749979*fEdge[38]+0.5590169943749475*fEdge[37])+(0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[31]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fEdge[18]+vcoord[1]*(0.4330127018922193*fEdge[17]+0.25*fEdge[10]); 
  Ghat[11] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[33]+0.5*vcoord[1]*fEdge[32]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fEdge[22]+0.5590169943749475*vcoord[7]*fEdge[21]+0.3872983346207416*vcoord[1]*fEdge[15]+(0.4330127018922193*fEdge[6]+0.25*fEdge[3])*vcoord[7]+0.223606797749979*vcoord[1]*fEdge[7]; 
  Ghat[12] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[34]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fEdge[24]+vcoord[1]*(0.4330127018922193*fEdge[23]+0.25*fEdge[13]); 
  Ghat[13] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[36]+0.5*vcoord[1]*fEdge[35]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fEdge[26]+0.5590169943749475*vcoord[7]*fEdge[25]+vcoord[1]*(0.3872983346207416*fEdge[16]+0.223606797749979*fEdge[9])+vcoord[7]*(0.4330127018922193*fEdge[8]+0.25*fEdge[4]); 
  Ghat[14] = vcoord[1]*(0.4330127018922193*fEdge[46]+0.25*fEdge[40])+vcoord[0]*(0.4330127018922193*fEdge[39]+0.25*fEdge[27]); 
  Ghat[15] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[41]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fEdge[29]+vcoord[1]*(0.4330127018922193*fEdge[28]+0.25*fEdge[14]); 
  Ghat[16] = vcoord[1]*(0.4330127018922193*fEdge[47]+0.25*fEdge[43])+vcoord[0]*(0.4330127018922193*fEdge[42]+0.25*fEdge[30]); 
  Ghat[17] = (0.276641667586244*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[45]+0.5*vcoord[1]*fEdge[44]+(0.159719141249985*vcoord[7]+0.25*vcoord[0])*fEdge[38]+0.5590169943749475*vcoord[7]*fEdge[37]+vcoord[1]*(0.3872983346207416*fEdge[31]+0.223606797749979*fEdge[18])+vcoord[7]*(0.4330127018922193*fEdge[17]+0.25*fEdge[10]); 
  Ghat[18] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[46]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fEdge[40]+vcoord[1]*(0.4330127018922193*fEdge[39]+0.25*fEdge[27]); 
  Ghat[19] = (0.3872983346207416*vcoord[7]+0.4330127018922193*vcoord[0])*fEdge[47]+(0.223606797749979*vcoord[7]+0.25*vcoord[0])*fEdge[43]+vcoord[1]*(0.4330127018922193*fEdge[42]+0.25*fEdge[30]); 

  } else { 

  Ghat[0] = (-0.4330127018922194*vcoord[7]*fSkin[20])+0.5590169943749476*vcoord[1]*fSkin[19]+0.25*vcoord[7]*fSkin[12]+0.5590169943749475*vcoord[0]*fSkin[11]-0.4330127018922193*vcoord[1]*fSkin[5]+0.25*vcoord[1]*fSkin[2]-0.4330127018922193*vcoord[0]*fSkin[1]+0.25*fSkin[0]*vcoord[0]; 
  Ghat[1] = (-0.3872983346207417*vcoord[1]*fSkin[20])+0.5000000000000001*vcoord[7]*fSkin[19]+0.5590169943749476*vcoord[0]*fSkin[19]+0.223606797749979*vcoord[1]*fSkin[12]+0.5590169943749475*vcoord[1]*fSkin[11]-0.3872983346207416*fSkin[5]*vcoord[7]+0.223606797749979*fSkin[2]*vcoord[7]-0.4330127018922193*vcoord[0]*fSkin[5]+0.25*vcoord[0]*fSkin[2]-0.4330127018922193*fSkin[1]*vcoord[1]+0.25*fSkin[0]*vcoord[1]; 
  Ghat[2] = (-0.4330127018922193*vcoord[7]*fSkin[33])+0.5590169943749475*vcoord[1]*fSkin[32]+0.2500000000000001*vcoord[7]*fSkin[22]+0.5590169943749476*vcoord[0]*fSkin[21]-0.4330127018922193*vcoord[1]*fSkin[15]+0.25*vcoord[1]*fSkin[7]-0.4330127018922193*vcoord[0]*fSkin[6]+0.25*vcoord[0]*fSkin[3]; 
  Ghat[3] = (-0.4330127018922193*vcoord[7]*fSkin[36])+0.5590169943749475*vcoord[1]*fSkin[35]+0.2500000000000001*vcoord[7]*fSkin[26]+0.5590169943749476*vcoord[0]*fSkin[25]-0.4330127018922193*vcoord[1]*fSkin[16]+0.25*vcoord[1]*fSkin[9]-0.4330127018922193*vcoord[0]*fSkin[8]+0.25*vcoord[0]*fSkin[4]; 
  Ghat[4] = (-0.3872983346207416*vcoord[1]*fSkin[33])+0.5*vcoord[7]*fSkin[32]+0.5590169943749475*vcoord[0]*fSkin[32]+0.223606797749979*vcoord[1]*fSkin[22]+0.5590169943749476*vcoord[1]*fSkin[21]-0.3872983346207416*vcoord[7]*fSkin[15]-0.4330127018922193*vcoord[0]*fSkin[15]+0.223606797749979*fSkin[7]*vcoord[7]+0.25*vcoord[0]*fSkin[7]-0.4330127018922193*vcoord[1]*fSkin[6]+0.25*vcoord[1]*fSkin[3]; 
  Ghat[5] = (-0.3872983346207416*vcoord[1]*fSkin[36])+0.5*vcoord[7]*fSkin[35]+0.5590169943749475*vcoord[0]*fSkin[35]+0.223606797749979*vcoord[1]*fSkin[26]+0.5590169943749476*vcoord[1]*fSkin[25]-0.3872983346207416*vcoord[7]*fSkin[16]-0.4330127018922193*vcoord[0]*fSkin[16]+0.223606797749979*vcoord[7]*fSkin[9]+0.25*vcoord[0]*fSkin[9]-0.4330127018922193*vcoord[1]*fSkin[8]+0.25*vcoord[1]*fSkin[4]; 
  Ghat[6] = (-0.4330127018922194*vcoord[7]*fSkin[45])+0.5590169943749476*vcoord[1]*fSkin[44]+0.25*vcoord[7]*fSkin[38]+0.5590169943749475*vcoord[0]*fSkin[37]-0.4330127018922193*vcoord[1]*fSkin[31]+0.25*vcoord[1]*fSkin[18]-0.4330127018922193*vcoord[0]*fSkin[17]+0.25*vcoord[0]*fSkin[10]; 
  Ghat[7] = (-0.276641667586244*vcoord[7]*fSkin[20])-0.4330127018922194*vcoord[0]*fSkin[20]+0.5000000000000001*vcoord[1]*fSkin[19]+0.159719141249985*vcoord[7]*fSkin[12]+0.25*vcoord[0]*fSkin[12]+0.5590169943749475*vcoord[7]*fSkin[11]-0.4330127018922193*fSkin[1]*vcoord[7]+0.25*fSkin[0]*vcoord[7]-0.3872983346207416*vcoord[1]*fSkin[5]+0.223606797749979*vcoord[1]*fSkin[2]; 
  Ghat[8] = (-0.4330127018922193*vcoord[1]*fSkin[34])+0.2500000000000001*vcoord[1]*fSkin[24]-0.4330127018922194*vcoord[0]*fSkin[23]+0.25*vcoord[0]*fSkin[13]; 
  Ghat[9] = (-0.4330127018922193*vcoord[1]*fSkin[41])+0.2500000000000001*vcoord[1]*fSkin[29]-0.4330127018922194*vcoord[0]*fSkin[28]+0.25*vcoord[0]*fSkin[14]; 
  Ghat[10] = (-0.3872983346207417*vcoord[1]*fSkin[45])+0.5000000000000001*vcoord[7]*fSkin[44]+0.5590169943749476*vcoord[0]*fSkin[44]+0.223606797749979*vcoord[1]*fSkin[38]+0.5590169943749475*vcoord[1]*fSkin[37]-0.3872983346207416*vcoord[7]*fSkin[31]-0.4330127018922193*vcoord[0]*fSkin[31]+0.223606797749979*vcoord[7]*fSkin[18]+0.25*vcoord[0]*fSkin[18]-0.4330127018922193*vcoord[1]*fSkin[17]+0.25*vcoord[1]*fSkin[10]; 
  Ghat[11] = (-0.276641667586244*vcoord[7]*fSkin[33])-0.4330127018922194*vcoord[0]*fSkin[33]+0.5000000000000001*vcoord[1]*fSkin[32]+0.159719141249985*vcoord[7]*fSkin[22]+0.25*vcoord[0]*fSkin[22]+0.5590169943749475*vcoord[7]*fSkin[21]-0.3872983346207417*vcoord[1]*fSkin[15]-0.4330127018922194*fSkin[6]*vcoord[7]+0.2500000000000001*fSkin[3]*vcoord[7]+0.223606797749979*vcoord[1]*fSkin[7]; 
  Ghat[12] = (-0.3872983346207417*vcoord[7]*fSkin[34])-0.4330127018922194*vcoord[0]*fSkin[34]+0.223606797749979*vcoord[7]*fSkin[24]+0.25*vcoord[0]*fSkin[24]-0.4330127018922193*vcoord[1]*fSkin[23]+0.2500000000000001*vcoord[1]*fSkin[13]; 
  Ghat[13] = (-0.276641667586244*vcoord[7]*fSkin[36])-0.4330127018922194*vcoord[0]*fSkin[36]+0.5000000000000001*vcoord[1]*fSkin[35]+0.159719141249985*vcoord[7]*fSkin[26]+0.25*vcoord[0]*fSkin[26]+0.5590169943749475*vcoord[7]*fSkin[25]-0.3872983346207417*vcoord[1]*fSkin[16]+0.223606797749979*vcoord[1]*fSkin[9]-0.4330127018922194*vcoord[7]*fSkin[8]+0.2500000000000001*fSkin[4]*vcoord[7]; 
  Ghat[14] = (-0.4330127018922193*vcoord[1]*fSkin[46])+0.2500000000000001*vcoord[1]*fSkin[40]-0.4330127018922194*vcoord[0]*fSkin[39]+0.25*vcoord[0]*fSkin[27]; 
  Ghat[15] = (-0.3872983346207417*vcoord[7]*fSkin[41])-0.4330127018922194*vcoord[0]*fSkin[41]+0.223606797749979*vcoord[7]*fSkin[29]+0.25*vcoord[0]*fSkin[29]-0.4330127018922193*vcoord[1]*fSkin[28]+0.2500000000000001*vcoord[1]*fSkin[14]; 
  Ghat[16] = (-0.4330127018922193*vcoord[1]*fSkin[47])+0.2500000000000001*vcoord[1]*fSkin[43]-0.4330127018922194*vcoord[0]*fSkin[42]+0.25*vcoord[0]*fSkin[30]; 
  Ghat[17] = (-0.276641667586244*vcoord[7]*fSkin[45])-0.4330127018922194*vcoord[0]*fSkin[45]+0.5000000000000001*vcoord[1]*fSkin[44]+0.159719141249985*vcoord[7]*fSkin[38]+0.25*vcoord[0]*fSkin[38]+0.5590169943749475*vcoord[7]*fSkin[37]-0.3872983346207416*vcoord[1]*fSkin[31]+0.223606797749979*vcoord[1]*fSkin[18]-0.4330127018922193*vcoord[7]*fSkin[17]+0.25*vcoord[7]*fSkin[10]; 
  Ghat[18] = (-0.3872983346207417*vcoord[7]*fSkin[46])-0.4330127018922194*vcoord[0]*fSkin[46]+0.223606797749979*vcoord[7]*fSkin[40]+0.25*vcoord[0]*fSkin[40]-0.4330127018922193*vcoord[1]*fSkin[39]+0.2500000000000001*vcoord[1]*fSkin[27]; 
  Ghat[19] = (-0.3872983346207417*vcoord[7]*fSkin[47])-0.4330127018922194*vcoord[0]*fSkin[47]+0.223606797749979*vcoord[7]*fSkin[43]+0.25*vcoord[0]*fSkin[43]-0.4330127018922193*vcoord[1]*fSkin[42]+0.2500000000000001*vcoord[1]*fSkin[30]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*rdx2; 
  out[1] += -1.224744871391589*Ghat[0]*rdx2; 
  out[2] += 0.7071067811865475*Ghat[1]*rdx2; 
  out[3] += 0.7071067811865475*Ghat[2]*rdx2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdx2; 
  out[5] += -1.224744871391589*Ghat[1]*rdx2; 
  out[6] += -1.224744871391589*Ghat[2]*rdx2; 
  out[7] += 0.7071067811865475*Ghat[4]*rdx2; 
  out[8] += -1.224744871391589*Ghat[3]*rdx2; 
  out[9] += 0.7071067811865475*Ghat[5]*rdx2; 
  out[10] += 0.7071067811865475*Ghat[6]*rdx2; 
  out[11] += 1.58113883008419*Ghat[0]*rdx2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdx2; 
  out[13] += 0.7071067811865475*Ghat[8]*rdx2; 
  out[14] += 0.7071067811865475*Ghat[9]*rdx2; 
  out[15] += -1.224744871391589*Ghat[4]*rdx2; 
  out[16] += -1.224744871391589*Ghat[5]*rdx2; 
  out[17] += -1.224744871391589*Ghat[6]*rdx2; 
  out[18] += 0.7071067811865475*Ghat[10]*rdx2; 
  out[19] += 1.58113883008419*Ghat[1]*rdx2; 
  out[20] += -1.224744871391589*Ghat[7]*rdx2; 
  out[21] += 1.58113883008419*Ghat[2]*rdx2; 
  out[22] += 0.7071067811865475*Ghat[11]*rdx2; 
  out[23] += -1.224744871391589*Ghat[8]*rdx2; 
  out[24] += 0.7071067811865475*Ghat[12]*rdx2; 
  out[25] += 1.58113883008419*Ghat[3]*rdx2; 
  out[26] += 0.7071067811865475*Ghat[13]*rdx2; 
  out[27] += 0.7071067811865475*Ghat[14]*rdx2; 
  out[28] += -1.224744871391589*Ghat[9]*rdx2; 
  out[29] += 0.7071067811865475*Ghat[15]*rdx2; 
  out[30] += 0.7071067811865475*Ghat[16]*rdx2; 
  out[31] += -1.224744871391589*Ghat[10]*rdx2; 
  out[32] += 1.58113883008419*Ghat[4]*rdx2; 
  out[33] += -1.224744871391589*Ghat[11]*rdx2; 
  out[34] += -1.224744871391589*Ghat[12]*rdx2; 
  out[35] += 1.58113883008419*Ghat[5]*rdx2; 
  out[36] += -1.224744871391589*Ghat[13]*rdx2; 
  out[37] += 1.58113883008419*Ghat[6]*rdx2; 
  out[38] += 0.7071067811865475*Ghat[17]*rdx2; 
  out[39] += -1.224744871391589*Ghat[14]*rdx2; 
  out[40] += 0.7071067811865475*Ghat[18]*rdx2; 
  out[41] += -1.224744871391589*Ghat[15]*rdx2; 
  out[42] += -1.224744871391589*Ghat[16]*rdx2; 
  out[43] += 0.7071067811865475*Ghat[19]*rdx2; 
  out[44] += 1.58113883008419*Ghat[10]*rdx2; 
  out[45] += -1.224744871391589*Ghat[17]*rdx2; 
  out[46] += -1.224744871391589*Ghat[18]*rdx2; 
  out[47] += -1.224744871391589*Ghat[19]*rdx2; 

  } 
  return 0.;

} 
