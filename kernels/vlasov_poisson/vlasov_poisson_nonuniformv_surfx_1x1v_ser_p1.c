#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *vcoord, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vcoord:    Discrete (DG) velocity coordinate.
  // alpha_geo: Fields used only for general geometry.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double rdx2 = 2/dxv[0]; 
  double Ghat_r[3]; 
  double Ghat_l[3]; 
  if (0.7071067811865475*vcoord[0]-0.7905694150420947*vcoord[2]>0) { 

  Ghat_r[0] = vcoord[2]*(0.8660254037844386*fc[5]+0.5*fc[4])+vcoord[1]*(0.8660254037844386*fc[3]+0.5*fc[2])+vcoord[0]*(0.8660254037844386*fc[1]+0.5*fc[0]); 
  Ghat_r[1] = vcoord[1]*(0.7745966692414833*fc[5]+0.4472135954999579*fc[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fc[3]+fc[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fc[1]+0.5*fc[0])*vcoord[1]; 
  Ghat_r[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fc[5]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fc[4]+0.7745966692414833*vcoord[1]*fc[3]+(0.8660254037844386*fc[1]+0.5*fc[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fc[2]; 

  Ghat_l[0] = vcoord[2]*(0.8660254037844386*fl[5]+0.5*fl[4])+vcoord[1]*(0.8660254037844386*fl[3]+0.5*fl[2])+vcoord[0]*(0.8660254037844386*fl[1]+0.5*fl[0]); 
  Ghat_l[1] = vcoord[1]*(0.7745966692414833*fl[5]+0.4472135954999579*fl[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fl[3]+fl[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fl[1]+0.5*fl[0])*vcoord[1]; 
  Ghat_l[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fl[5]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fl[4]+0.7745966692414833*vcoord[1]*fl[3]+(0.8660254037844386*fl[1]+0.5*fl[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fl[2]; 

  } else { 

  Ghat_r[0] = -0.1*(vcoord[2]*(8.660254037844387*fr[5]-5.0*fr[4])+vcoord[1]*(8.660254037844386*fr[3]-5.0*fr[2])+vcoord[0]*(8.660254037844386*fr[1]-5.0*fr[0])); 
  Ghat_r[1] = -0.1*(vcoord[1]*(7.745966692414834*fr[5]-4.47213595499958*fr[4])+(7.745966692414834*vcoord[2]+8.660254037844386*vcoord[0])*fr[3]+fr[2]*((-4.47213595499958*vcoord[2])-5.0*vcoord[0])+(8.660254037844386*fr[1]-5.0*fr[0])*vcoord[1]); 
  Ghat_r[2] = -0.01428571428571429*((38.72983346207417*vcoord[2]+60.62177826491071*vcoord[0])*fr[5]+((-22.3606797749979*vcoord[2])-35.0*vcoord[0])*fr[4]+54.22176684690384*vcoord[1]*fr[3]+(60.6217782649107*fr[1]-35.0*fr[0])*vcoord[2]-31.30495168499706*vcoord[1]*fr[2]); 

  Ghat_l[0] = -0.1*(vcoord[2]*(8.660254037844387*fc[5]-5.0*fc[4])+vcoord[1]*(8.660254037844386*fc[3]-5.0*fc[2])+vcoord[0]*(8.660254037844386*fc[1]-5.0*fc[0])); 
  Ghat_l[1] = -0.1*(vcoord[1]*(7.745966692414834*fc[5]-4.47213595499958*fc[4])+(7.745966692414834*vcoord[2]+8.660254037844386*vcoord[0])*fc[3]+fc[2]*((-4.47213595499958*vcoord[2])-5.0*vcoord[0])+(8.660254037844386*fc[1]-5.0*fc[0])*vcoord[1]); 
  Ghat_l[2] = -0.01428571428571429*((38.72983346207417*vcoord[2]+60.62177826491071*vcoord[0])*fc[5]+((-22.3606797749979*vcoord[2])-35.0*vcoord[0])*fc[4]+54.22176684690384*vcoord[1]*fc[3]+(60.6217782649107*fc[1]-35.0*fc[0])*vcoord[2]-31.30495168499706*vcoord[1]*fc[2]); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*rdx2; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdx2; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*rdx2; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdx2; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*rdx2; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdx2; 

  return 0.;

} 
