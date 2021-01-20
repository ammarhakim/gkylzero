#include <gkyl_vlasov_kernels.h> 
double vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* restrict out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // amax:        amax in global lax flux.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[10]; 

  double Ghat[4]; 
  double favg[4]; 
  double alpha[4]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double amid = 0.5*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[2]-1.224744871391589*fEdge[2]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[4]-1.224744871391589*fEdge[4]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = 1.224744871391589*fSkin[6]-1.224744871391589*fEdge[6]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[3] = 1.224744871391589*fSkin[7]-1.224744871391589*fEdge[7]+0.7071067811865475*(fSkin[5]+fEdge[5]); 

  Ghat[0] = (0.6123724356957944*(fSkin[2]+fEdge[2])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[4]+fEdge[4])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.25*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[6]+fEdge[6])+0.3535533905932737*fSkin[3]-0.3535533905932737*fEdge[3])*amax+0.25*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[7]+fEdge[7])+0.3535533905932737*fSkin[5]-0.3535533905932737*fEdge[5])*amax+0.25*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[2])+1.224744871391589*fEdge[2]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[4])+1.224744871391589*fEdge[4]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = (-1.224744871391589*fSkin[6])+1.224744871391589*fEdge[6]+0.7071067811865475*(fSkin[3]+fEdge[3]); 
  favg[3] = (-1.224744871391589*fSkin[7])+1.224744871391589*fEdge[7]+0.7071067811865475*(fSkin[5]+fEdge[5]); 

  Ghat[0] = (0.6123724356957944*(fSkin[2]+fEdge[2])-0.3535533905932737*fSkin[0]+0.3535533905932737*fEdge[0])*amax+0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[4]+fEdge[4])-0.3535533905932737*fSkin[1]+0.3535533905932737*fEdge[1])*amax+0.25*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[6]+fEdge[6])-0.3535533905932737*fSkin[3]+0.3535533905932737*fEdge[3])*amax+0.25*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[7]+fEdge[7])-0.3535533905932737*fSkin[5]+0.3535533905932737*fEdge[5])*amax+0.25*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 

  } 
  return fabs(amid); 
} 
double vlasov_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* restrict out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // amax:        amax in global lax flux.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E1 = &qmem[2]; 
  const double *B2 = &qmem[10]; 

  double Ghat[4]; 
  double favg[4]; 
  double alpha[4]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 

  double amid = 0.5*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[3]-1.224744871391589*fEdge[3]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[5]-1.224744871391589*fEdge[5]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = 1.224744871391589*fSkin[6]-1.224744871391589*fEdge[6]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = 1.224744871391589*fSkin[7]-1.224744871391589*fEdge[7]+0.7071067811865475*(fSkin[4]+fEdge[4]); 

  Ghat[0] = (0.6123724356957944*(fSkin[3]+fEdge[3])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[5]+fEdge[5])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.25*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[6]+fEdge[6])+0.3535533905932737*fSkin[2]-0.3535533905932737*fEdge[2])*amax+0.25*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[7]+fEdge[7])+0.3535533905932737*fSkin[4]-0.3535533905932737*fEdge[4])*amax+0.25*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += -0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -1.224744871391589*Ghat[3]*dv11; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[3])+1.224744871391589*fEdge[3]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[5])+1.224744871391589*fEdge[5]+0.7071067811865475*(fSkin[1]+fEdge[1]); 
  favg[2] = (-1.224744871391589*fSkin[6])+1.224744871391589*fEdge[6]+0.7071067811865475*(fSkin[2]+fEdge[2]); 
  favg[3] = (-1.224744871391589*fSkin[7])+1.224744871391589*fEdge[7]+0.7071067811865475*(fSkin[4]+fEdge[4]); 

  Ghat[0] = (0.6123724356957944*(fSkin[3]+fEdge[3])-0.3535533905932737*fSkin[0]+0.3535533905932737*fEdge[0])*amax+0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[5]+fEdge[5])-0.3535533905932737*fSkin[1]+0.3535533905932737*fEdge[1])*amax+0.25*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = (0.6123724356957944*(fSkin[6]+fEdge[6])-0.3535533905932737*fSkin[2]+0.3535533905932737*fEdge[2])*amax+0.25*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = (0.6123724356957944*(fSkin[7]+fEdge[7])-0.3535533905932737*fSkin[4]+0.3535533905932737*fEdge[4])*amax+0.25*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -1.224744871391589*Ghat[3]*dv11; 

  } 
  return fabs(amid); 
} 
