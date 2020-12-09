#include <gkyl_vlasov_kernels.h> 
void vlasov_surf_1x1v_x_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[8]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[4]+1.732050807568877*fl[1]+fl[0])*wxl+(0.6454972243679028*fl[6]+0.5*fl[3]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[4])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-1.118033988749895*fl[6])-0.8660254037844386*fl[3]-0.5*fl[2])*dvxl; 
  incr[2] = (2.23606797749979*fl[6]+1.732050807568877*fl[3]+fl[2])*wxl+(0.447213595499958*fl[7]+0.2581988897471612*fl[5]+0.6454972243679029*fl[4]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = ((-3.872983346207417*fl[6])-3.0*fl[3]-1.732050807568877*fl[2])*wxl+((-0.7745966692414834*fl[7])-0.4472135954999579*fl[5]-1.118033988749895*fl[4]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[4] = (5.0*fl[4]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.443375672974065*fl[6]+1.118033988749895*fl[3]+0.6454972243679029*fl[2])*dvxl; 
  incr[5] = (1.732050807568877*fl[7]+fl[5])*wxl+(0.5773502691896257*fl[6]+0.4472135954999579*fl[3]+0.2581988897471612*fl[2])*dvxl; 
  incr[6] = (5.0*fl[6]+3.872983346207417*fl[3]+2.23606797749979*fl[2])*wxl+(fl[7]+0.5773502691896257*fl[5]+1.443375672974065*fl[4]+1.118033988749895*fl[1]+0.6454972243679028*fl[0])*dvxl; 
  incr[7] = ((-3.0*fl[7])-1.732050807568877*fl[5])*wxl+((-1.0*fl[6])-0.7745966692414834*fl[3]-0.447213595499958*fl[2])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[4]-1.732050807568877*fr[1]+fr[0])*wxr+(0.6454972243679028*fr[6]-0.5*fr[3]+0.2886751345948129*fr[2])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[4])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+((-1.118033988749895*fr[6])+0.8660254037844386*fr[3]-0.5*fr[2])*dvxr; 
  incr[2] = (2.23606797749979*fr[6]-1.732050807568877*fr[3]+fr[2])*wxr+((-0.447213595499958*fr[7])+0.2581988897471612*fr[5]+0.6454972243679029*fr[4]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = ((-3.872983346207417*fr[6])+3.0*fr[3]-1.732050807568877*fr[2])*wxr+(0.7745966692414834*fr[7]-0.4472135954999579*fr[5]-1.118033988749895*fr[4]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[4] = (5.0*fr[4]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(1.443375672974065*fr[6]-1.118033988749895*fr[3]+0.6454972243679029*fr[2])*dvxr; 
  incr[5] = (fr[5]-1.732050807568877*fr[7])*wxr+(0.5773502691896257*fr[6]-0.4472135954999579*fr[3]+0.2581988897471612*fr[2])*dvxr; 
  incr[6] = (5.0*fr[6]-3.872983346207417*fr[3]+2.23606797749979*fr[2])*wxr+((-1.0*fr[7])+0.5773502691896257*fr[5]+1.443375672974065*fr[4]-1.118033988749895*fr[1]+0.6454972243679028*fr[0])*dvxr; 
  incr[7] = (3.0*fr[7]-1.732050807568877*fr[5])*wxr+((-1.0*fr[6])+0.7745966692414834*fr[3]-0.447213595499958*fr[2])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  } 
} 
double vlasov_surf_1x1v_vx_p2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 

  double Ghat[3]; 
  double favg[3]; 
  double alpha[3]; 

  favg[0] = 1.58113883008419*fr[5]+1.58113883008419*fl[5]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 1.58113883008419*fr[7]+1.58113883008419*fl[7]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 

  const double amid = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 

  Ghat[0] = 0.3535533905932737*(alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[5]-1.0*(2.23606797749979*fl[5]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.07071067811865475*(4.47213595499958*(alpha[1]*favg[2]+favg[1]*alpha[2])+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[7]-1.0*(3.872983346207417*fl[7]+3.0*(fr[3]+fl[3])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.07071067811865475*(8.660254037844387*(fr[6]+fl[6])+5.0*(fl[4]-1.0*fr[4]))*amax+0.01010152544552211*((22.3606797749979*alpha[2]+35.0*alpha[0])*favg[2]+35.0*favg[0]*alpha[2]+31.30495168499706*alpha[1]*favg[1]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[5] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += 1.58113883008419*Ghat[1]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[5] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -1.58113883008419*Ghat[1]*dv10l; 

  return fabs(amid); 
} 
