#include <gkyl_vlasov_bflux_kernels.h> 
GKYL_CU_DH void vlasov_bflux_1x1v_ser_p2(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fIn:  Input Distribution function in cell.
  // out:       Incremented distribution function in center cell.
 if (edge == GKYL_VX_UPPER) {

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[3]; 
  Ghat[0] = (1.581138830084189*fIn[4]+1.224744871391589*fIn[1]+0.7071067811865475*fIn[0])*wv+(0.4564354645876384*fIn[6]+0.3535533905932737*fIn[3]+0.2041241452319314*fIn[2])*dv; 
  Ghat[1] = (1.581138830084189*fIn[6]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[2])*wv+(0.3162277660168379*fIn[7]+0.1825741858350553*fIn[5]+0.4564354645876384*fIn[4]+0.3535533905932737*fIn[1]+0.2041241452319314*fIn[0])*dv; 
  Ghat[2] = (1.224744871391589*fIn[7]+0.7071067811865475*fIn[5])*wv+(0.408248290463863*fIn[6]+0.3162277660168379*fIn[3]+0.1825741858350553*fIn[2])*dv; 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += -1.581138830084189*Ghat[0]*dx10; 
  out[5] += -0.7071067811865475*Ghat[2]*dx10; 
  out[6] += -1.581138830084189*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 

 } else if (edge == GKYL_VX_LOWER) {

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[3]; 
  Ghat[0] = 1.581138830084189*fIn[4]*wv-1.224744871391589*fIn[1]*wv+0.7071067811865475*fIn[0]*wv+0.4564354645876383*fIn[6]*dv-0.3535533905932737*fIn[3]*dv+0.2041241452319314*fIn[2]*dv; 
  Ghat[1] = 1.581138830084189*fIn[6]*wv-1.224744871391589*fIn[3]*wv+0.7071067811865475*fIn[2]*wv-0.3162277660168379*fIn[7]*dv+0.1825741858350553*fIn[5]*dv+0.4564354645876384*fIn[4]*dv-0.3535533905932737*fIn[1]*dv+0.2041241452319314*fIn[0]*dv; 
  Ghat[2] = (-1.224744871391589*fIn[7]*wv)+0.7071067811865475*fIn[5]*wv+0.4082482904638629*fIn[6]*dv-0.3162277660168379*fIn[3]*dv+0.1825741858350553*fIn[2]*dv; 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += 1.581138830084189*Ghat[0]*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dx10; 
  out[6] += 1.581138830084189*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 
  } 

} 

