#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_bflux_1x1v_ser_p1(const int *idx, const double *w, enum gkyl_vel_edge edge, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fIn:  Input Distribution function in cell.
  // out:       Incremented distribution function in center cell.
 if (edge == GKYL_VX_UPPER) {

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[2]; 
  Ghat[0] = (1.224744871391589*fIn[1]+0.7071067811865475*fIn[0])*wv+(0.3535533905932737*fIn[3]+0.2041241452319314*fIn[2])*dv; 
  Ghat[1] = (1.224744871391589*fIn[3]+0.7071067811865475*fIn[2])*wv+(0.3535533905932737*fIn[1]+0.2041241452319314*fIn[0])*dv; 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 

 } else if (edge == GKYL_VX_LOWER) {

  Ghat[0] = -0.1178511301977578*((10.39230484541326*fIn[1]-6.0*fIn[0])*wv+(3.0*fIn[3]-1.732050807568877*fIn[2])*dv); 
  Ghat[1] = -0.1178511301977578*((10.39230484541326*fIn[3]-6.0*fIn[2])*wv+(3.0*fIn[1]-1.732050807568877*fIn[0])*dv); 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  } 

} 

