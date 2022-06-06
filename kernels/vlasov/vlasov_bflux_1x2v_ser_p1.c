#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_bflux_1x2v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fIn:  Input Distribution function in cell.
  // out:       Incremented distribution function in center cell.
 if (edge == GKYL_VX_UPPER) {

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[4]; 
  Ghat[0] = (1.224744871391589*fIn[1]+0.7071067811865475*fIn[0])*wv+(0.3535533905932737*fIn[4]+0.2041241452319314*fIn[2])*dv; 
  Ghat[1] = (1.224744871391589*fIn[4]+0.7071067811865475*fIn[2])*wv+(0.3535533905932737*fIn[1]+0.2041241452319314*fIn[0])*dv; 
  Ghat[2] = (1.224744871391589*fIn[5]+0.7071067811865475*fIn[3])*wv+(0.3535533905932737*fIn[7]+0.2041241452319314*fIn[6])*dv; 
  Ghat[3] = (1.224744871391589*fIn[7]+0.7071067811865475*fIn[6])*wv+(0.3535533905932737*fIn[5]+0.2041241452319314*fIn[3])*dv; 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += -0.7071067811865475*Ghat[3]*dx10; 
  out[7] += -1.224744871391589*Ghat[3]*dx10; 

 } else if (edge == GKYL_VX_LOWER) {

  Ghat[0] = -0.08333333333333333*((14.69693845669906*fIn[1]-8.485281374238572*fIn[0])*wv+(4.242640687119286*fIn[4]-2.449489742783178*fIn[2])*dv); 
  Ghat[1] = -0.08333333333333333*((14.69693845669906*fIn[4]-8.485281374238572*fIn[2])*wv+(4.242640687119286*fIn[1]-2.449489742783178*fIn[0])*dv); 
  Ghat[2] = -0.08333333333333333*((14.69693845669906*fIn[5]-8.485281374238572*fIn[3])*wv+(4.242640687119286*fIn[7]-2.449489742783178*fIn[6])*dv); 
  Ghat[3] = -0.08333333333333333*((14.69693845669906*fIn[7]-8.485281374238572*fIn[6])*wv+(4.242640687119286*fIn[5]-2.449489742783178*fIn[3])*dv); 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += 0.7071067811865475*Ghat[3]*dx10; 
  out[7] += -1.224744871391589*Ghat[3]*dx10; 
  } 

} 

