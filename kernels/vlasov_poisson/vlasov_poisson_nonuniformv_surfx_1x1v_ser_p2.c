#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *vcoord, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  Ghat_r[0] = 0.8660254037844386*vcoord[2]*fc[7]+1.118033988749895*vcoord[1]*fc[6]+0.5*vcoord[2]*fc[5]+1.118033988749895*vcoord[0]*fc[4]+vcoord[1]*(0.8660254037844386*fc[3]+0.5*fc[2])+vcoord[0]*(0.8660254037844386*fc[1]+0.5*fc[0]); 
  Ghat_r[1] = 0.7745966692414833*vcoord[1]*fc[7]+(vcoord[2]+1.118033988749895*vcoord[0])*fc[6]+vcoord[1]*(0.4472135954999579*fc[5]+1.118033988749895*fc[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fc[3]+fc[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fc[1]+0.5*fc[0])*vcoord[1]; 
  Ghat_r[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fc[7]+vcoord[1]*fc[6]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fc[5]+1.118033988749895*vcoord[2]*fc[4]+0.7745966692414833*vcoord[1]*fc[3]+(0.8660254037844386*fc[1]+0.5*fc[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fc[2]; 

  Ghat_l[0] = 0.8660254037844386*vcoord[2]*fl[7]+1.118033988749895*vcoord[1]*fl[6]+0.5*vcoord[2]*fl[5]+1.118033988749895*vcoord[0]*fl[4]+vcoord[1]*(0.8660254037844386*fl[3]+0.5*fl[2])+vcoord[0]*(0.8660254037844386*fl[1]+0.5*fl[0]); 
  Ghat_l[1] = 0.7745966692414833*vcoord[1]*fl[7]+(vcoord[2]+1.118033988749895*vcoord[0])*fl[6]+vcoord[1]*(0.4472135954999579*fl[5]+1.118033988749895*fl[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fl[3]+fl[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fl[1]+0.5*fl[0])*vcoord[1]; 
  Ghat_l[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fl[7]+vcoord[1]*fl[6]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fl[5]+1.118033988749895*vcoord[2]*fl[4]+0.7745966692414833*vcoord[1]*fl[3]+(0.8660254037844386*fl[1]+0.5*fl[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fl[2]; 

  } else { 

  Ghat_r[0] = (-0.8660254037844387*vcoord[2]*fr[7])+1.118033988749895*vcoord[1]*fr[6]+0.5*vcoord[2]*fr[5]+1.118033988749895*vcoord[0]*fr[4]-0.8660254037844386*vcoord[1]*fr[3]+0.5*vcoord[1]*fr[2]-0.8660254037844386*vcoord[0]*fr[1]+0.5*fr[0]*vcoord[0]; 
  Ghat_r[1] = (-0.7745966692414834*vcoord[1]*fr[7])+1.0*vcoord[2]*fr[6]+1.118033988749895*vcoord[0]*fr[6]+0.4472135954999579*vcoord[1]*fr[5]+1.118033988749895*vcoord[1]*fr[4]-0.7745966692414833*vcoord[2]*fr[3]-0.8660254037844386*vcoord[0]*fr[3]+0.4472135954999579*fr[2]*vcoord[2]+0.5*vcoord[0]*fr[2]-0.8660254037844386*fr[1]*vcoord[1]+0.5*fr[0]*vcoord[1]; 
  Ghat_r[2] = (-0.5532833351724881*vcoord[2]*fr[7])-0.8660254037844387*vcoord[0]*fr[7]+1.0*vcoord[1]*fr[6]+0.31943828249997*vcoord[2]*fr[5]+0.5*vcoord[0]*fr[5]+1.118033988749895*vcoord[2]*fr[4]-0.7745966692414833*vcoord[1]*fr[3]-0.8660254037844386*fr[1]*vcoord[2]+0.5*fr[0]*vcoord[2]+0.4472135954999579*vcoord[1]*fr[2]; 

  Ghat_l[0] = (-0.8660254037844387*vcoord[2]*fc[7])+1.118033988749895*vcoord[1]*fc[6]+0.5*vcoord[2]*fc[5]+1.118033988749895*vcoord[0]*fc[4]-0.8660254037844386*vcoord[1]*fc[3]+0.5*vcoord[1]*fc[2]-0.8660254037844386*vcoord[0]*fc[1]+0.5*fc[0]*vcoord[0]; 
  Ghat_l[1] = (-0.7745966692414834*vcoord[1]*fc[7])+1.0*vcoord[2]*fc[6]+1.118033988749895*vcoord[0]*fc[6]+0.4472135954999579*vcoord[1]*fc[5]+1.118033988749895*vcoord[1]*fc[4]-0.7745966692414833*vcoord[2]*fc[3]-0.8660254037844386*vcoord[0]*fc[3]+0.4472135954999579*fc[2]*vcoord[2]+0.5*vcoord[0]*fc[2]-0.8660254037844386*fc[1]*vcoord[1]+0.5*fc[0]*vcoord[1]; 
  Ghat_l[2] = (-0.5532833351724881*vcoord[2]*fc[7])-0.8660254037844387*vcoord[0]*fc[7]+1.0*vcoord[1]*fc[6]+0.31943828249997*vcoord[2]*fc[5]+0.5*vcoord[0]*fc[5]+1.118033988749895*vcoord[2]*fc[4]-0.7745966692414833*vcoord[1]*fc[3]-0.8660254037844386*fc[1]*vcoord[2]+0.5*fc[0]*vcoord[2]+0.4472135954999579*vcoord[1]*fc[2]; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*rdx2; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdx2; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*rdx2; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdx2; 
  out[4] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*rdx2; 
  out[5] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*rdx2; 
  out[6] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*rdx2; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdx2; 

  return 0.;

} 
