#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // amax:        amax in global lax flux.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[12]; 
  const double *B1 = &qmem[16]; 
  const double *B2 = &qmem[20]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[3]*dv3; 

  double amid = 0.25*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[3]-1.224744871391589*fEdge[3]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[7]-1.224744871391589*fEdge[7]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = 1.224744871391589*fSkin[8]-1.224744871391589*fEdge[8]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = 1.224744871391589*fSkin[11]-1.224744871391589*fEdge[11]+0.7071067811865475*(fSkin[4]+fEdge[4]); 
  favg[4] = 1.224744871391589*fSkin[14]-1.224744871391589*fEdge[14]+0.7071067811865475*(fSkin[5]+fEdge[5]); 
  favg[5] = 1.224744871391589*fSkin[16]-1.224744871391589*fEdge[16]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = 1.224744871391589*fSkin[18]-1.224744871391589*fEdge[18]+0.7071067811865475*(fSkin[9]+fEdge[9]); 
  favg[7] = 1.224744871391589*fSkin[19]-1.224744871391589*fEdge[19]+0.7071067811865475*(fSkin[10]+fEdge[10]); 
  favg[8] = 1.224744871391589*fSkin[21]-1.224744871391589*fEdge[21]+0.7071067811865475*(fSkin[12]+fEdge[12]); 
  favg[9] = 1.224744871391589*fSkin[22]-1.224744871391589*fEdge[22]+0.7071067811865475*(fSkin[13]+fEdge[13]); 
  favg[10] = 1.224744871391589*fSkin[25]-1.224744871391589*fEdge[25]+0.7071067811865475*(fSkin[15]+fEdge[15]); 
  favg[11] = 1.224744871391589*fSkin[26]-1.224744871391589*fEdge[26]+0.7071067811865475*(fSkin[17]+fEdge[17]); 
  favg[12] = 1.224744871391589*fSkin[27]-1.224744871391589*fEdge[27]+0.7071067811865475*(fSkin[20]+fEdge[20]); 
  favg[13] = 1.224744871391589*fSkin[29]-1.224744871391589*fEdge[29]+0.7071067811865475*(fSkin[23]+fEdge[23]); 
  favg[14] = 1.224744871391589*fSkin[30]-1.224744871391589*fEdge[30]+0.7071067811865475*(fSkin[24]+fEdge[24]); 
  favg[15] = 1.224744871391589*fSkin[31]-1.224744871391589*fEdge[31]+0.7071067811865475*(fSkin[28]+fEdge[28]); 

  Ghat[0] = (0.6123724356957944*(fSkin[3]+fEdge[3])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[7]+fEdge[7])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[8]+fEdge[8])+0.3535533905932737*fSkin[2]-0.3535533905932737*fEdge[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[11]+fEdge[11])+0.3535533905932737*fSkin[4]-0.3535533905932737*fEdge[4])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[14]+fEdge[14])+0.3535533905932737*fSkin[5]-0.3535533905932737*fEdge[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = (0.6123724356957944*(fSkin[16]+fEdge[16])+0.3535533905932737*fSkin[6]-0.3535533905932737*fEdge[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = (0.6123724356957944*(fSkin[18]+fEdge[18])+0.3535533905932737*fSkin[9]-0.3535533905932737*fEdge[9])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[19]+fEdge[19])+0.3535533905932737*fSkin[10]-0.3535533905932737*fEdge[10])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = (0.6123724356957944*(fSkin[21]+fEdge[21])+0.3535533905932737*fSkin[12]-0.3535533905932737*fEdge[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = (0.6123724356957944*(fSkin[22]+fEdge[22])+0.3535533905932737*fSkin[13]-0.3535533905932737*fEdge[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = (0.6123724356957944*(fSkin[25]+fEdge[25])+0.3535533905932737*fSkin[15]-0.3535533905932737*fEdge[15])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = (0.6123724356957944*(fSkin[26]+fEdge[26])+0.3535533905932737*fSkin[17]-0.3535533905932737*fEdge[17])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = (0.6123724356957944*(fSkin[27]+fEdge[27])+0.3535533905932737*fSkin[20]-0.3535533905932737*fEdge[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = (0.6123724356957944*(fSkin[29]+fEdge[29])+0.3535533905932737*fSkin[23]-0.3535533905932737*fEdge[23])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = (0.6123724356957944*(fSkin[30]+fEdge[30])+0.3535533905932737*fSkin[24]-0.3535533905932737*fEdge[24])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = (0.6123724356957944*(fSkin[31]+fEdge[31])+0.3535533905932737*fSkin[28]-0.3535533905932737*fEdge[28])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += -0.7071067811865475*Ghat[8]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -0.7071067811865475*Ghat[11]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 
  out[20] += -0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[8]*dv10; 
  out[22] += -1.224744871391589*Ghat[9]*dv10; 
  out[23] += -0.7071067811865475*Ghat[13]*dv10; 
  out[24] += -0.7071067811865475*Ghat[14]*dv10; 
  out[25] += -1.224744871391589*Ghat[10]*dv10; 
  out[26] += -1.224744871391589*Ghat[11]*dv10; 
  out[27] += -1.224744871391589*Ghat[12]*dv10; 
  out[28] += -0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -1.224744871391589*Ghat[13]*dv10; 
  out[30] += -1.224744871391589*Ghat[14]*dv10; 
  out[31] += -1.224744871391589*Ghat[15]*dv10; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[3])+1.224744871391589*fEdge[3]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[7])+1.224744871391589*fEdge[7]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = (-1.224744871391589*fSkin[8])+1.224744871391589*fEdge[8]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = (-1.224744871391589*fSkin[11])+1.224744871391589*fEdge[11]+0.7071067811865475*(fSkin[4]+fEdge[4]); 
  favg[4] = (-1.224744871391589*fSkin[14])+1.224744871391589*fEdge[14]+0.7071067811865475*(fSkin[5]+fEdge[5]); 
  favg[5] = (-1.224744871391589*fSkin[16])+1.224744871391589*fEdge[16]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = (-1.224744871391589*fSkin[18])+1.224744871391589*fEdge[18]+0.7071067811865475*(fSkin[9]+fEdge[9]); 
  favg[7] = (-1.224744871391589*fSkin[19])+1.224744871391589*fEdge[19]+0.7071067811865475*(fSkin[10]+fEdge[10]); 
  favg[8] = (-1.224744871391589*fSkin[21])+1.224744871391589*fEdge[21]+0.7071067811865475*(fSkin[12]+fEdge[12]); 
  favg[9] = (-1.224744871391589*fSkin[22])+1.224744871391589*fEdge[22]+0.7071067811865475*(fSkin[13]+fEdge[13]); 
  favg[10] = (-1.224744871391589*fSkin[25])+1.224744871391589*fEdge[25]+0.7071067811865475*(fSkin[15]+fEdge[15]); 
  favg[11] = (-1.224744871391589*fSkin[26])+1.224744871391589*fEdge[26]+0.7071067811865475*(fSkin[17]+fEdge[17]); 
  favg[12] = (-1.224744871391589*fSkin[27])+1.224744871391589*fEdge[27]+0.7071067811865475*(fSkin[20]+fEdge[20]); 
  favg[13] = (-1.224744871391589*fSkin[29])+1.224744871391589*fEdge[29]+0.7071067811865475*(fSkin[23]+fEdge[23]); 
  favg[14] = (-1.224744871391589*fSkin[30])+1.224744871391589*fEdge[30]+0.7071067811865475*(fSkin[24]+fEdge[24]); 
  favg[15] = (-1.224744871391589*fSkin[31])+1.224744871391589*fEdge[31]+0.7071067811865475*(fSkin[28]+fEdge[28]); 

  Ghat[0] = (0.6123724356957944*(fSkin[3]+fEdge[3])-0.3535533905932737*fSkin[0]+0.3535533905932737*fEdge[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[7]+fEdge[7])-0.3535533905932737*fSkin[1]+0.3535533905932737*fEdge[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[8]+fEdge[8])-0.3535533905932737*fSkin[2]+0.3535533905932737*fEdge[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[11]+fEdge[11])-0.3535533905932737*fSkin[4]+0.3535533905932737*fEdge[4])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[14]+fEdge[14])-0.3535533905932737*fSkin[5]+0.3535533905932737*fEdge[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = (0.6123724356957944*(fSkin[16]+fEdge[16])-0.3535533905932737*fSkin[6]+0.3535533905932737*fEdge[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = (0.6123724356957944*(fSkin[18]+fEdge[18])-0.3535533905932737*fSkin[9]+0.3535533905932737*fEdge[9])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[19]+fEdge[19])-0.3535533905932737*fSkin[10]+0.3535533905932737*fEdge[10])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = (0.6123724356957944*(fSkin[21]+fEdge[21])-0.3535533905932737*fSkin[12]+0.3535533905932737*fEdge[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = (0.6123724356957944*(fSkin[22]+fEdge[22])-0.3535533905932737*fSkin[13]+0.3535533905932737*fEdge[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = (0.6123724356957944*(fSkin[25]+fEdge[25])-0.3535533905932737*fSkin[15]+0.3535533905932737*fEdge[15])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = (0.6123724356957944*(fSkin[26]+fEdge[26])-0.3535533905932737*fSkin[17]+0.3535533905932737*fEdge[17])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = (0.6123724356957944*(fSkin[27]+fEdge[27])-0.3535533905932737*fSkin[20]+0.3535533905932737*fEdge[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = (0.6123724356957944*(fSkin[29]+fEdge[29])-0.3535533905932737*fSkin[23]+0.3535533905932737*fEdge[23])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = (0.6123724356957944*(fSkin[30]+fEdge[30])-0.3535533905932737*fSkin[24]+0.3535533905932737*fEdge[24])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = (0.6123724356957944*(fSkin[31]+fEdge[31])-0.3535533905932737*fSkin[28]+0.3535533905932737*fEdge[28])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += 0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += 0.7071067811865475*Ghat[11]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 
  out[20] += 0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[8]*dv10; 
  out[22] += -1.224744871391589*Ghat[9]*dv10; 
  out[23] += 0.7071067811865475*Ghat[13]*dv10; 
  out[24] += 0.7071067811865475*Ghat[14]*dv10; 
  out[25] += -1.224744871391589*Ghat[10]*dv10; 
  out[26] += -1.224744871391589*Ghat[11]*dv10; 
  out[27] += -1.224744871391589*Ghat[12]*dv10; 
  out[28] += 0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -1.224744871391589*Ghat[13]*dv10; 
  out[30] += -1.224744871391589*Ghat[14]*dv10; 
  out[31] += -1.224744871391589*Ghat[15]*dv10; 

  } 
  return fabs(amid); 
} 
GKYL_CU_DH double vlasov_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // amax:        amax in global lax flux.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E1 = &qmem[4]; 
  const double *B0 = &qmem[12]; 
  const double *B1 = &qmem[16]; 
  const double *B2 = &qmem[20]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = 0.5773502691896258*B0[0]*dv3; 
  alpha[5] = 2.0*B0[3]*wv3-2.0*B2[3]*wv1+2.0*E1[3]; 
  alpha[6] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = -0.5773502691896258*B2[2]*dv1; 
  alpha[8] = 0.5773502691896258*B0[1]*dv3; 
  alpha[9] = 0.5773502691896258*B0[2]*dv3; 
  alpha[11] = -0.5773502691896258*B2[3]*dv1; 
  alpha[12] = 0.5773502691896258*B0[3]*dv3; 

  double amid = 0.25*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[4]-1.224744871391589*fEdge[4]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[9]-1.224744871391589*fEdge[9]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = 1.224744871391589*fSkin[10]-1.224744871391589*fEdge[10]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = 1.224744871391589*fSkin[11]-1.224744871391589*fEdge[11]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[4] = 1.224744871391589*fSkin[15]-1.224744871391589*fEdge[15]+0.7071067811865475*(fSkin[5]+fEdge[5]); 
  favg[5] = 1.224744871391589*fSkin[17]-1.224744871391589*fEdge[17]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = 1.224744871391589*fSkin[18]-1.224744871391589*fEdge[18]+0.7071067811865475*(fSkin[7]+fEdge[7]); 
  favg[7] = 1.224744871391589*fSkin[19]-1.224744871391589*fEdge[19]+0.7071067811865475*(fSkin[8]+fEdge[8]); 
  favg[8] = 1.224744871391589*fSkin[23]-1.224744871391589*fEdge[23]+0.7071067811865475*(fSkin[12]+fEdge[12]); 
  favg[9] = 1.224744871391589*fSkin[24]-1.224744871391589*fEdge[24]+0.7071067811865475*(fSkin[13]+fEdge[13]); 
  favg[10] = 1.224744871391589*fSkin[25]-1.224744871391589*fEdge[25]+0.7071067811865475*(fSkin[14]+fEdge[14]); 
  favg[11] = 1.224744871391589*fSkin[26]-1.224744871391589*fEdge[26]+0.7071067811865475*(fSkin[16]+fEdge[16]); 
  favg[12] = 1.224744871391589*fSkin[28]-1.224744871391589*fEdge[28]+0.7071067811865475*(fSkin[20]+fEdge[20]); 
  favg[13] = 1.224744871391589*fSkin[29]-1.224744871391589*fEdge[29]+0.7071067811865475*(fSkin[21]+fEdge[21]); 
  favg[14] = 1.224744871391589*fSkin[30]-1.224744871391589*fEdge[30]+0.7071067811865475*(fSkin[22]+fEdge[22]); 
  favg[15] = 1.224744871391589*fSkin[31]-1.224744871391589*fEdge[31]+0.7071067811865475*(fSkin[27]+fEdge[27]); 

  Ghat[0] = (0.6123724356957944*(fSkin[4]+fEdge[4])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[9]+fEdge[9])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[10]+fEdge[10])+0.3535533905932737*fSkin[2]-0.3535533905932737*fEdge[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[11]+fEdge[11])+0.3535533905932737*fSkin[3]-0.3535533905932737*fEdge[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[15]+fEdge[15])+0.3535533905932737*fSkin[5]-0.3535533905932737*fEdge[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = (0.6123724356957944*(fSkin[17]+fEdge[17])+0.3535533905932737*fSkin[6]-0.3535533905932737*fEdge[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = (0.6123724356957944*(fSkin[18]+fEdge[18])+0.3535533905932737*fSkin[7]-0.3535533905932737*fEdge[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[19]+fEdge[19])+0.3535533905932737*fSkin[8]-0.3535533905932737*fEdge[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = (0.6123724356957944*(fSkin[23]+fEdge[23])+0.3535533905932737*fSkin[12]-0.3535533905932737*fEdge[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = (0.6123724356957944*(fSkin[24]+fEdge[24])+0.3535533905932737*fSkin[13]-0.3535533905932737*fEdge[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = (0.6123724356957944*(fSkin[25]+fEdge[25])+0.3535533905932737*fSkin[14]-0.3535533905932737*fEdge[14])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = (0.6123724356957944*(fSkin[26]+fEdge[26])+0.3535533905932737*fSkin[16]-0.3535533905932737*fEdge[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = (0.6123724356957944*(fSkin[28]+fEdge[28])+0.3535533905932737*fSkin[20]-0.3535533905932737*fEdge[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = (0.6123724356957944*(fSkin[29]+fEdge[29])+0.3535533905932737*fSkin[21]-0.3535533905932737*fEdge[21])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = (0.6123724356957944*(fSkin[30]+fEdge[30])+0.3535533905932737*fSkin[22]-0.3535533905932737*fEdge[22])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = (0.6123724356957944*(fSkin[31]+fEdge[31])+0.3535533905932737*fSkin[27]-0.3535533905932737*fEdge[27])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += -0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -1.224744871391589*Ghat[1]*dv11; 
  out[10] += -1.224744871391589*Ghat[2]*dv11; 
  out[11] += -1.224744871391589*Ghat[3]*dv11; 
  out[12] += -0.7071067811865475*Ghat[8]*dv11; 
  out[13] += -0.7071067811865475*Ghat[9]*dv11; 
  out[14] += -0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -1.224744871391589*Ghat[4]*dv11; 
  out[16] += -0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -1.224744871391589*Ghat[7]*dv11; 
  out[20] += -0.7071067811865475*Ghat[12]*dv11; 
  out[21] += -0.7071067811865475*Ghat[13]*dv11; 
  out[22] += -0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.224744871391589*Ghat[9]*dv11; 
  out[25] += -1.224744871391589*Ghat[10]*dv11; 
  out[26] += -1.224744871391589*Ghat[11]*dv11; 
  out[27] += -0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -1.224744871391589*Ghat[12]*dv11; 
  out[29] += -1.224744871391589*Ghat[13]*dv11; 
  out[30] += -1.224744871391589*Ghat[14]*dv11; 
  out[31] += -1.224744871391589*Ghat[15]*dv11; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[4])+1.224744871391589*fEdge[4]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[9])+1.224744871391589*fEdge[9]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = (-1.224744871391589*fSkin[10])+1.224744871391589*fEdge[10]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = (-1.224744871391589*fSkin[11])+1.224744871391589*fEdge[11]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[4] = (-1.224744871391589*fSkin[15])+1.224744871391589*fEdge[15]+0.7071067811865475*(fSkin[5]+fEdge[5]); 
  favg[5] = (-1.224744871391589*fSkin[17])+1.224744871391589*fEdge[17]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = (-1.224744871391589*fSkin[18])+1.224744871391589*fEdge[18]+0.7071067811865475*(fSkin[7]+fEdge[7]); 
  favg[7] = (-1.224744871391589*fSkin[19])+1.224744871391589*fEdge[19]+0.7071067811865475*(fSkin[8]+fEdge[8]); 
  favg[8] = (-1.224744871391589*fSkin[23])+1.224744871391589*fEdge[23]+0.7071067811865475*(fSkin[12]+fEdge[12]); 
  favg[9] = (-1.224744871391589*fSkin[24])+1.224744871391589*fEdge[24]+0.7071067811865475*(fSkin[13]+fEdge[13]); 
  favg[10] = (-1.224744871391589*fSkin[25])+1.224744871391589*fEdge[25]+0.7071067811865475*(fSkin[14]+fEdge[14]); 
  favg[11] = (-1.224744871391589*fSkin[26])+1.224744871391589*fEdge[26]+0.7071067811865475*(fSkin[16]+fEdge[16]); 
  favg[12] = (-1.224744871391589*fSkin[28])+1.224744871391589*fEdge[28]+0.7071067811865475*(fSkin[20]+fEdge[20]); 
  favg[13] = (-1.224744871391589*fSkin[29])+1.224744871391589*fEdge[29]+0.7071067811865475*(fSkin[21]+fEdge[21]); 
  favg[14] = (-1.224744871391589*fSkin[30])+1.224744871391589*fEdge[30]+0.7071067811865475*(fSkin[22]+fEdge[22]); 
  favg[15] = (-1.224744871391589*fSkin[31])+1.224744871391589*fEdge[31]+0.7071067811865475*(fSkin[27]+fEdge[27]); 

  Ghat[0] = (0.6123724356957944*(fSkin[4]+fEdge[4])-0.3535533905932737*fSkin[0]+0.3535533905932737*fEdge[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[9]+fEdge[9])-0.3535533905932737*fSkin[1]+0.3535533905932737*fEdge[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[10]+fEdge[10])-0.3535533905932737*fSkin[2]+0.3535533905932737*fEdge[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[11]+fEdge[11])-0.3535533905932737*fSkin[3]+0.3535533905932737*fEdge[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[15]+fEdge[15])-0.3535533905932737*fSkin[5]+0.3535533905932737*fEdge[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = (0.6123724356957944*(fSkin[17]+fEdge[17])-0.3535533905932737*fSkin[6]+0.3535533905932737*fEdge[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = (0.6123724356957944*(fSkin[18]+fEdge[18])-0.3535533905932737*fSkin[7]+0.3535533905932737*fEdge[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[19]+fEdge[19])-0.3535533905932737*fSkin[8]+0.3535533905932737*fEdge[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = (0.6123724356957944*(fSkin[23]+fEdge[23])-0.3535533905932737*fSkin[12]+0.3535533905932737*fEdge[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = (0.6123724356957944*(fSkin[24]+fEdge[24])-0.3535533905932737*fSkin[13]+0.3535533905932737*fEdge[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = (0.6123724356957944*(fSkin[25]+fEdge[25])-0.3535533905932737*fSkin[14]+0.3535533905932737*fEdge[14])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = (0.6123724356957944*(fSkin[26]+fEdge[26])-0.3535533905932737*fSkin[16]+0.3535533905932737*fEdge[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = (0.6123724356957944*(fSkin[28]+fEdge[28])-0.3535533905932737*fSkin[20]+0.3535533905932737*fEdge[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = (0.6123724356957944*(fSkin[29]+fEdge[29])-0.3535533905932737*fSkin[21]+0.3535533905932737*fEdge[21])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = (0.6123724356957944*(fSkin[30]+fEdge[30])-0.3535533905932737*fSkin[22]+0.3535533905932737*fEdge[22])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = (0.6123724356957944*(fSkin[31]+fEdge[31])-0.3535533905932737*fSkin[27]+0.3535533905932737*fEdge[27])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += 0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -1.224744871391589*Ghat[1]*dv11; 
  out[10] += -1.224744871391589*Ghat[2]*dv11; 
  out[11] += -1.224744871391589*Ghat[3]*dv11; 
  out[12] += 0.7071067811865475*Ghat[8]*dv11; 
  out[13] += 0.7071067811865475*Ghat[9]*dv11; 
  out[14] += 0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -1.224744871391589*Ghat[4]*dv11; 
  out[16] += 0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -1.224744871391589*Ghat[7]*dv11; 
  out[20] += 0.7071067811865475*Ghat[12]*dv11; 
  out[21] += 0.7071067811865475*Ghat[13]*dv11; 
  out[22] += 0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.224744871391589*Ghat[9]*dv11; 
  out[25] += -1.224744871391589*Ghat[10]*dv11; 
  out[26] += -1.224744871391589*Ghat[11]*dv11; 
  out[27] += 0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -1.224744871391589*Ghat[12]*dv11; 
  out[29] += -1.224744871391589*Ghat[13]*dv11; 
  out[30] += -1.224744871391589*Ghat[14]*dv11; 
  out[31] += -1.224744871391589*Ghat[15]*dv11; 

  } 
  return fabs(amid); 
} 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // amax:        amax in global lax flux.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv12 = 2/dxv[4]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E2 = &qmem[8]; 
  const double *B0 = &qmem[12]; 
  const double *B1 = &qmem[16]; 
  const double *B2 = &qmem[20]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 2.0*(B1[2]*wv1+E2[2])-2.0*B0[2]*wv2; 
  alpha[3] = 0.5773502691896258*B1[0]*dv1; 
  alpha[4] = -0.5773502691896258*B0[0]*dv2; 
  alpha[5] = 2.0*(B1[3]*wv1+E2[3])-2.0*B0[3]*wv2; 
  alpha[6] = 0.5773502691896258*B1[1]*dv1; 
  alpha[7] = 0.5773502691896258*B1[2]*dv1; 
  alpha[8] = -0.5773502691896258*B0[1]*dv2; 
  alpha[9] = -0.5773502691896258*B0[2]*dv2; 
  alpha[11] = 0.5773502691896258*B1[3]*dv1; 
  alpha[12] = -0.5773502691896258*B0[3]*dv2; 

  double amid = 0.25*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[5]-1.224744871391589*fEdge[5]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[12]-1.224744871391589*fEdge[12]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = 1.224744871391589*fSkin[13]-1.224744871391589*fEdge[13]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = 1.224744871391589*fSkin[14]-1.224744871391589*fEdge[14]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[4] = 1.224744871391589*fSkin[15]-1.224744871391589*fEdge[15]+0.7071067811865475*(fSkin[4]+fEdge[4]); 
  favg[5] = 1.224744871391589*fSkin[20]-1.224744871391589*fEdge[20]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = 1.224744871391589*fSkin[21]-1.224744871391589*fEdge[21]+0.7071067811865475*(fSkin[7]+fEdge[7]); 
  favg[7] = 1.224744871391589*fSkin[22]-1.224744871391589*fEdge[22]+0.7071067811865475*(fSkin[8]+fEdge[8]); 
  favg[8] = 1.224744871391589*fSkin[23]-1.224744871391589*fEdge[23]+0.7071067811865475*(fSkin[9]+fEdge[9]); 
  favg[9] = 1.224744871391589*fSkin[24]-1.224744871391589*fEdge[24]+0.7071067811865475*(fSkin[10]+fEdge[10]); 
  favg[10] = 1.224744871391589*fSkin[25]-1.224744871391589*fEdge[25]+0.7071067811865475*(fSkin[11]+fEdge[11]); 
  favg[11] = 1.224744871391589*fSkin[27]-1.224744871391589*fEdge[27]+0.7071067811865475*(fSkin[16]+fEdge[16]); 
  favg[12] = 1.224744871391589*fSkin[28]-1.224744871391589*fEdge[28]+0.7071067811865475*(fSkin[17]+fEdge[17]); 
  favg[13] = 1.224744871391589*fSkin[29]-1.224744871391589*fEdge[29]+0.7071067811865475*(fSkin[18]+fEdge[18]); 
  favg[14] = 1.224744871391589*fSkin[30]-1.224744871391589*fEdge[30]+0.7071067811865475*(fSkin[19]+fEdge[19]); 
  favg[15] = 1.224744871391589*fSkin[31]-1.224744871391589*fEdge[31]+0.7071067811865475*(fSkin[26]+fEdge[26]); 

  Ghat[0] = (0.6123724356957944*(fSkin[5]+fEdge[5])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[12]+fEdge[12])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[13]+fEdge[13])+0.3535533905932737*fSkin[2]-0.3535533905932737*fEdge[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[14]+fEdge[14])+0.3535533905932737*fSkin[3]-0.3535533905932737*fEdge[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[15]+fEdge[15])+0.3535533905932737*fSkin[4]-0.3535533905932737*fEdge[4])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = (0.6123724356957944*(fSkin[20]+fEdge[20])+0.3535533905932737*fSkin[6]-0.3535533905932737*fEdge[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = (0.6123724356957944*(fSkin[21]+fEdge[21])+0.3535533905932737*fSkin[7]-0.3535533905932737*fEdge[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[22]+fEdge[22])+0.3535533905932737*fSkin[8]-0.3535533905932737*fEdge[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = (0.6123724356957944*(fSkin[23]+fEdge[23])+0.3535533905932737*fSkin[9]-0.3535533905932737*fEdge[9])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = (0.6123724356957944*(fSkin[24]+fEdge[24])+0.3535533905932737*fSkin[10]-0.3535533905932737*fEdge[10])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = (0.6123724356957944*(fSkin[25]+fEdge[25])+0.3535533905932737*fSkin[11]-0.3535533905932737*fEdge[11])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = (0.6123724356957944*(fSkin[27]+fEdge[27])+0.3535533905932737*fSkin[16]-0.3535533905932737*fEdge[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = (0.6123724356957944*(fSkin[28]+fEdge[28])+0.3535533905932737*fSkin[17]-0.3535533905932737*fEdge[17])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = (0.6123724356957944*(fSkin[29]+fEdge[29])+0.3535533905932737*fSkin[18]-0.3535533905932737*fEdge[18])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = (0.6123724356957944*(fSkin[30]+fEdge[30])+0.3535533905932737*fSkin[19]-0.3535533905932737*fEdge[19])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = (0.6123724356957944*(fSkin[31]+fEdge[31])+0.3535533905932737*fSkin[26]-0.3535533905932737*fEdge[26])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -0.7071067811865475*Ghat[4]*dv12; 
  out[5] += -1.224744871391589*Ghat[0]*dv12; 
  out[6] += -0.7071067811865475*Ghat[5]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -0.7071067811865475*Ghat[7]*dv12; 
  out[9] += -0.7071067811865475*Ghat[8]*dv12; 
  out[10] += -0.7071067811865475*Ghat[9]*dv12; 
  out[11] += -0.7071067811865475*Ghat[10]*dv12; 
  out[12] += -1.224744871391589*Ghat[1]*dv12; 
  out[13] += -1.224744871391589*Ghat[2]*dv12; 
  out[14] += -1.224744871391589*Ghat[3]*dv12; 
  out[15] += -1.224744871391589*Ghat[4]*dv12; 
  out[16] += -0.7071067811865475*Ghat[11]*dv12; 
  out[17] += -0.7071067811865475*Ghat[12]*dv12; 
  out[18] += -0.7071067811865475*Ghat[13]*dv12; 
  out[19] += -0.7071067811865475*Ghat[14]*dv12; 
  out[20] += -1.224744871391589*Ghat[5]*dv12; 
  out[21] += -1.224744871391589*Ghat[6]*dv12; 
  out[22] += -1.224744871391589*Ghat[7]*dv12; 
  out[23] += -1.224744871391589*Ghat[8]*dv12; 
  out[24] += -1.224744871391589*Ghat[9]*dv12; 
  out[25] += -1.224744871391589*Ghat[10]*dv12; 
  out[26] += -0.7071067811865475*Ghat[15]*dv12; 
  out[27] += -1.224744871391589*Ghat[11]*dv12; 
  out[28] += -1.224744871391589*Ghat[12]*dv12; 
  out[29] += -1.224744871391589*Ghat[13]*dv12; 
  out[30] += -1.224744871391589*Ghat[14]*dv12; 
  out[31] += -1.224744871391589*Ghat[15]*dv12; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[5])+1.224744871391589*fEdge[5]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[12])+1.224744871391589*fEdge[12]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = (-1.224744871391589*fSkin[13])+1.224744871391589*fEdge[13]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = (-1.224744871391589*fSkin[14])+1.224744871391589*fEdge[14]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[4] = (-1.224744871391589*fSkin[15])+1.224744871391589*fEdge[15]+0.7071067811865475*(fSkin[4]+fEdge[4]); 
  favg[5] = (-1.224744871391589*fSkin[20])+1.224744871391589*fEdge[20]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = (-1.224744871391589*fSkin[21])+1.224744871391589*fEdge[21]+0.7071067811865475*(fSkin[7]+fEdge[7]); 
  favg[7] = (-1.224744871391589*fSkin[22])+1.224744871391589*fEdge[22]+0.7071067811865475*(fSkin[8]+fEdge[8]); 
  favg[8] = (-1.224744871391589*fSkin[23])+1.224744871391589*fEdge[23]+0.7071067811865475*(fSkin[9]+fEdge[9]); 
  favg[9] = (-1.224744871391589*fSkin[24])+1.224744871391589*fEdge[24]+0.7071067811865475*(fSkin[10]+fEdge[10]); 
  favg[10] = (-1.224744871391589*fSkin[25])+1.224744871391589*fEdge[25]+0.7071067811865475*(fSkin[11]+fEdge[11]); 
  favg[11] = (-1.224744871391589*fSkin[27])+1.224744871391589*fEdge[27]+0.7071067811865475*(fSkin[16]+fEdge[16]); 
  favg[12] = (-1.224744871391589*fSkin[28])+1.224744871391589*fEdge[28]+0.7071067811865475*(fSkin[17]+fEdge[17]); 
  favg[13] = (-1.224744871391589*fSkin[29])+1.224744871391589*fEdge[29]+0.7071067811865475*(fSkin[18]+fEdge[18]); 
  favg[14] = (-1.224744871391589*fSkin[30])+1.224744871391589*fEdge[30]+0.7071067811865475*(fSkin[19]+fEdge[19]); 
  favg[15] = (-1.224744871391589*fSkin[31])+1.224744871391589*fEdge[31]+0.7071067811865475*(fSkin[26]+fEdge[26]); 

  Ghat[0] = (0.6123724356957944*(fSkin[5]+fEdge[5])-0.3535533905932737*fSkin[0]+0.3535533905932737*fEdge[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[12]+fEdge[12])-0.3535533905932737*fSkin[1]+0.3535533905932737*fEdge[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[13]+fEdge[13])-0.3535533905932737*fSkin[2]+0.3535533905932737*fEdge[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[14]+fEdge[14])-0.3535533905932737*fSkin[3]+0.3535533905932737*fEdge[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[15]+fEdge[15])-0.3535533905932737*fSkin[4]+0.3535533905932737*fEdge[4])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = (0.6123724356957944*(fSkin[20]+fEdge[20])-0.3535533905932737*fSkin[6]+0.3535533905932737*fEdge[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = (0.6123724356957944*(fSkin[21]+fEdge[21])-0.3535533905932737*fSkin[7]+0.3535533905932737*fEdge[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[22]+fEdge[22])-0.3535533905932737*fSkin[8]+0.3535533905932737*fEdge[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = (0.6123724356957944*(fSkin[23]+fEdge[23])-0.3535533905932737*fSkin[9]+0.3535533905932737*fEdge[9])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = (0.6123724356957944*(fSkin[24]+fEdge[24])-0.3535533905932737*fSkin[10]+0.3535533905932737*fEdge[10])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = (0.6123724356957944*(fSkin[25]+fEdge[25])-0.3535533905932737*fSkin[11]+0.3535533905932737*fEdge[11])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = (0.6123724356957944*(fSkin[27]+fEdge[27])-0.3535533905932737*fSkin[16]+0.3535533905932737*fEdge[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = (0.6123724356957944*(fSkin[28]+fEdge[28])-0.3535533905932737*fSkin[17]+0.3535533905932737*fEdge[17])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = (0.6123724356957944*(fSkin[29]+fEdge[29])-0.3535533905932737*fSkin[18]+0.3535533905932737*fEdge[18])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = (0.6123724356957944*(fSkin[30]+fEdge[30])-0.3535533905932737*fSkin[19]+0.3535533905932737*fEdge[19])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = (0.6123724356957944*(fSkin[31]+fEdge[31])-0.3535533905932737*fSkin[26]+0.3535533905932737*fEdge[26])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += 0.7071067811865475*Ghat[4]*dv12; 
  out[5] += -1.224744871391589*Ghat[0]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += 0.7071067811865475*Ghat[7]*dv12; 
  out[9] += 0.7071067811865475*Ghat[8]*dv12; 
  out[10] += 0.7071067811865475*Ghat[9]*dv12; 
  out[11] += 0.7071067811865475*Ghat[10]*dv12; 
  out[12] += -1.224744871391589*Ghat[1]*dv12; 
  out[13] += -1.224744871391589*Ghat[2]*dv12; 
  out[14] += -1.224744871391589*Ghat[3]*dv12; 
  out[15] += -1.224744871391589*Ghat[4]*dv12; 
  out[16] += 0.7071067811865475*Ghat[11]*dv12; 
  out[17] += 0.7071067811865475*Ghat[12]*dv12; 
  out[18] += 0.7071067811865475*Ghat[13]*dv12; 
  out[19] += 0.7071067811865475*Ghat[14]*dv12; 
  out[20] += -1.224744871391589*Ghat[5]*dv12; 
  out[21] += -1.224744871391589*Ghat[6]*dv12; 
  out[22] += -1.224744871391589*Ghat[7]*dv12; 
  out[23] += -1.224744871391589*Ghat[8]*dv12; 
  out[24] += -1.224744871391589*Ghat[9]*dv12; 
  out[25] += -1.224744871391589*Ghat[10]*dv12; 
  out[26] += 0.7071067811865475*Ghat[15]*dv12; 
  out[27] += -1.224744871391589*Ghat[11]*dv12; 
  out[28] += -1.224744871391589*Ghat[12]*dv12; 
  out[29] += -1.224744871391589*Ghat[13]*dv12; 
  out[30] += -1.224744871391589*Ghat[14]*dv12; 
  out[31] += -1.224744871391589*Ghat[15]*dv12; 

  } 
  return fabs(amid); 
} 
