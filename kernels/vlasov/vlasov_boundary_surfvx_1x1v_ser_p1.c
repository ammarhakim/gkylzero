#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* restrict out) 
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
  const double *E0 = &qmem[0]; 

  double Ghat[2]; 
  double favg[2]; 
  double alpha[2]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  double amid = 0.7071067811865475*alpha[0]; 

  if (edge == -1) { 

  favg[0] = 1.224744871391589*fSkin[2]-1.224744871391589*fEdge[2]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = 1.224744871391589*fSkin[3]-1.224744871391589*fEdge[3]+0.7071067811865475*(fSkin[1]+fEdge[1]); 

  Ghat[0] = (0.6123724356957944*(fSkin[2]+fEdge[2])+0.3535533905932737*fSkin[0]-0.3535533905932737*fEdge[0])*amax+0.3535533905932737*(alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[3]+fEdge[3])+0.3535533905932737*fSkin[1]-0.3535533905932737*fEdge[1])*amax+0.3535533905932737*(alpha[0]*favg[1]+favg[0]*alpha[1]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 

  } else { 

  favg[0] = (-1.224744871391589*fSkin[2])+1.224744871391589*fEdge[2]+0.7071067811865475*(fSkin[0]+fEdge[0]); 
  favg[1] = (-1.224744871391589*fSkin[3])+1.224744871391589*fEdge[3]+0.7071067811865475*(fSkin[1]+fEdge[1]); 

  Ghat[0] = (0.6123724356957944*(fSkin[2]+fEdge[2])-0.3535533905932737*fSkin[0])*amax+0.3535533905932737*(fEdge[0]*amax+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = (0.6123724356957944*(fSkin[3]+fEdge[3])-0.3535533905932737*fSkin[1])*amax+0.3535533905932737*(fEdge[1]*amax+alpha[0]*favg[1]+favg[0]*alpha[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 

  } 
  return fabs(amid); 
} 
