#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* restrict out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // amax:        amax in global lax flux.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv12 = 2/dxv[3]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E2 = &qmem[4]; 
  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double Ghat[8]; 
  double favg[8]; 
  double alpha[8]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 0.5773502691896258*B1[0]*dv1; 
  alpha[3] = -0.5773502691896258*B0[0]*dv2; 
  alpha[4] = 0.5773502691896258*B1[1]*dv1; 
  alpha[5] = -0.5773502691896258*B0[1]*dv2; 

  double amid = 0.3535533905932737*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[4]-1.224744871391589*fEdge[4]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[8]-1.224744871391589*fEdge[8]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = 1.224744871391589*fSkin[9]-1.224744871391589*fEdge[9]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = 1.224744871391589*fSkin[10]-1.224744871391589*fEdge[10]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[4] = 1.224744871391589*fSkin[12]-1.224744871391589*fEdge[12]+0.7071067811865475*(fSkin[5]+fEdge[5]); 
  favg[5] = 1.224744871391589*fSkin[13]-1.224744871391589*fEdge[13]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = 1.224744871391589*fSkin[14]-1.224744871391589*fEdge[14]+0.7071067811865475*(fSkin[7]+fEdge[7]); 
  favg[7] = 1.224744871391589*fSkin[15]-1.224744871391589*fEdge[15]+0.7071067811865475*(fSkin[11]+fEdge[11]); 

  Ghat[0] = (0.6123724356957944*(fSkin[4]+fEdge[4])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.1767766952966368*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[8]+fEdge[8])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[9]+fEdge[9])+0.3535533905932737*fSkin[2]-0.3535533905932737*fEdge[2])*amax+0.1767766952966368*(alpha[5]*favg[7]+alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[10]+fEdge[10])+0.3535533905932737*fSkin[3]-0.3535533905932737*fEdge[3])*amax+0.1767766952966368*(alpha[4]*favg[7]+alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[12]+fEdge[12])+0.3535533905932737*fSkin[5]-0.3535533905932737*fEdge[5])*amax+0.1767766952966368*(alpha[3]*favg[7]+alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = (0.6123724356957944*(fSkin[13]+fEdge[13])+0.3535533905932737*fSkin[6]-0.3535533905932737*fEdge[6])*amax+0.1767766952966368*(alpha[2]*favg[7]+alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = (0.6123724356957944*(fSkin[14]+fEdge[14])+0.3535533905932737*fSkin[7]-0.3535533905932737*fEdge[7])*amax+0.1767766952966368*(alpha[1]*favg[7]+alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[15]+fEdge[15])+0.3535533905932737*fSkin[11]-0.3535533905932737*fEdge[11])*amax+0.1767766952966368*(alpha[0]*favg[7]+alpha[1]*favg[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += -0.7071067811865475*Ghat[4]*dv12; 
  out[6] += -0.7071067811865475*Ghat[5]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += -0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -1.224744871391589*Ghat[4]*dv12; 
  out[13] += -1.224744871391589*Ghat[5]*dv12; 
  out[14] += -1.224744871391589*Ghat[6]*dv12; 
  out[15] += -1.224744871391589*Ghat[7]*dv12; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[4])+1.224744871391589*fEdge[4]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[8])+1.224744871391589*fEdge[8]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = (-1.224744871391589*fSkin[9])+1.224744871391589*fEdge[9]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = (-1.224744871391589*fSkin[10])+1.224744871391589*fEdge[10]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[4] = (-1.224744871391589*fSkin[12])+1.224744871391589*fEdge[12]+0.7071067811865475*(fSkin[5]+fEdge[5]); 
  favg[5] = (-1.224744871391589*fSkin[13])+1.224744871391589*fEdge[13]+0.7071067811865475*(fSkin[6]+fEdge[6]); 
  favg[6] = (-1.224744871391589*fSkin[14])+1.224744871391589*fEdge[14]+0.7071067811865475*(fSkin[7]+fEdge[7]); 
  favg[7] = (-1.224744871391589*fSkin[15])+1.224744871391589*fEdge[15]+0.7071067811865475*(fSkin[11]+fEdge[11]); 

  Ghat[0] = (0.6123724356957944*(fSkin[4]+fEdge[4])-0.3535533905932737*fSkin[0]+0.3535533905932737*fEdge[0])*amax+0.1767766952966368*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[8]+fEdge[8])-0.3535533905932737*fSkin[1]+0.3535533905932737*fEdge[1])*amax+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[9]+fEdge[9])-0.3535533905932737*fSkin[2]+0.3535533905932737*fEdge[2])*amax+0.1767766952966368*(alpha[5]*favg[7]+alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[10]+fEdge[10])-0.3535533905932737*fSkin[3]+0.3535533905932737*fEdge[3])*amax+0.1767766952966368*(alpha[4]*favg[7]+alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = (0.6123724356957944*(fSkin[12]+fEdge[12])-0.3535533905932737*fSkin[5]+0.3535533905932737*fEdge[5])*amax+0.1767766952966368*(alpha[3]*favg[7]+alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = (0.6123724356957944*(fSkin[13]+fEdge[13])-0.3535533905932737*fSkin[6]+0.3535533905932737*fEdge[6])*amax+0.1767766952966368*(alpha[2]*favg[7]+alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = (0.6123724356957944*(fSkin[14]+fEdge[14])-0.3535533905932737*fSkin[7]+0.3535533905932737*fEdge[7])*amax+0.1767766952966368*(alpha[1]*favg[7]+alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = (0.6123724356957944*(fSkin[15]+fEdge[15])-0.3535533905932737*fSkin[11]+0.3535533905932737*fEdge[11])*amax+0.1767766952966368*(alpha[0]*favg[7]+alpha[1]*favg[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -1.224744871391589*Ghat[0]*dv12; 
  out[5] += 0.7071067811865475*Ghat[4]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -1.224744871391589*Ghat[1]*dv12; 
  out[9] += -1.224744871391589*Ghat[2]*dv12; 
  out[10] += -1.224744871391589*Ghat[3]*dv12; 
  out[11] += 0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -1.224744871391589*Ghat[4]*dv12; 
  out[13] += -1.224744871391589*Ghat[5]*dv12; 
  out[14] += -1.224744871391589*Ghat[6]*dv12; 
  out[15] += -1.224744871391589*Ghat[7]*dv12; 

  } 
  return fabs(amid); 
} 
