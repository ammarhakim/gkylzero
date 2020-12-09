#include <gkyl_vlasov_kernels.h> 
void vlasov_surf_1x3v_x_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[16]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[1]+fl[0])*wxl+(0.5*fl[5]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.0*fl[1])-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[5])-0.5*fl[2])*dvxl; 
  incr[2] = (1.732050807568877*fl[5]+fl[2])*wxl+(0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = (1.732050807568877*fl[6]+fl[3])*wxl+(0.5*fl[11]+0.2886751345948129*fl[7])*dvxl; 
  incr[4] = (1.732050807568877*fl[8]+fl[4])*wxl+(0.5*fl[12]+0.2886751345948129*fl[9])*dvxl; 
  incr[5] = ((-3.0*fl[5])-1.732050807568877*fl[2])*wxl+((-0.8660254037844386*fl[1])-0.5*fl[0])*dvxl; 
  incr[6] = ((-3.0*fl[6])-1.732050807568877*fl[3])*wxl+((-0.8660254037844386*fl[11])-0.5*fl[7])*dvxl; 
  incr[7] = (1.732050807568877*fl[11]+fl[7])*wxl+(0.5*fl[6]+0.2886751345948129*fl[3])*dvxl; 
  incr[8] = ((-3.0*fl[8])-1.732050807568877*fl[4])*wxl+((-0.8660254037844386*fl[12])-0.5*fl[9])*dvxl; 
  incr[9] = (1.732050807568877*fl[12]+fl[9])*wxl+(0.5*fl[8]+0.2886751345948129*fl[4])*dvxl; 
  incr[10] = (1.732050807568877*fl[13]+fl[10])*wxl+(0.5*fl[15]+0.2886751345948129*fl[14])*dvxl; 
  incr[11] = ((-3.0*fl[11])-1.732050807568877*fl[7])*wxl+((-0.8660254037844386*fl[6])-0.5*fl[3])*dvxl; 
  incr[12] = ((-3.0*fl[12])-1.732050807568877*fl[9])*wxl+((-0.8660254037844386*fl[8])-0.5*fl[4])*dvxl; 
  incr[13] = ((-3.0*fl[13])-1.732050807568877*fl[10])*wxl+((-0.8660254037844386*fl[15])-0.5*fl[14])*dvxl; 
  incr[14] = (1.732050807568877*fl[15]+fl[14])*wxl+(0.5*fl[13]+0.2886751345948129*fl[10])*dvxl; 
  incr[15] = ((-3.0*fl[15])-1.732050807568877*fl[14])*wxl+((-0.8660254037844386*fl[13])-0.5*fl[10])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[1])*wxr+(0.2886751345948129*fr[2]-0.5*fr[5])*dvxr; 
  incr[1] = (3.0*fr[1]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[5]-0.5*fr[2])*dvxr; 
  incr[2] = (fr[2]-1.732050807568877*fr[5])*wxr+(0.2886751345948129*fr[0]-0.5*fr[1])*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[6])*wxr+(0.2886751345948129*fr[7]-0.5*fr[11])*dvxr; 
  incr[4] = (fr[4]-1.732050807568877*fr[8])*wxr+(0.2886751345948129*fr[9]-0.5*fr[12])*dvxr; 
  incr[5] = (3.0*fr[5]-1.732050807568877*fr[2])*wxr+(0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[6] = (3.0*fr[6]-1.732050807568877*fr[3])*wxr+(0.8660254037844386*fr[11]-0.5*fr[7])*dvxr; 
  incr[7] = (fr[7]-1.732050807568877*fr[11])*wxr+(0.2886751345948129*fr[3]-0.5*fr[6])*dvxr; 
  incr[8] = (3.0*fr[8]-1.732050807568877*fr[4])*wxr+(0.8660254037844386*fr[12]-0.5*fr[9])*dvxr; 
  incr[9] = (fr[9]-1.732050807568877*fr[12])*wxr+(0.2886751345948129*fr[4]-0.5*fr[8])*dvxr; 
  incr[10] = (fr[10]-1.732050807568877*fr[13])*wxr+(0.2886751345948129*fr[14]-0.5*fr[15])*dvxr; 
  incr[11] = (3.0*fr[11]-1.732050807568877*fr[7])*wxr+(0.8660254037844386*fr[6]-0.5*fr[3])*dvxr; 
  incr[12] = (3.0*fr[12]-1.732050807568877*fr[9])*wxr+(0.8660254037844386*fr[8]-0.5*fr[4])*dvxr; 
  incr[13] = (3.0*fr[13]-1.732050807568877*fr[10])*wxr+(0.8660254037844386*fr[15]-0.5*fr[14])*dvxr; 
  incr[14] = (fr[14]-1.732050807568877*fr[15])*wxr+(0.2886751345948129*fr[10]-0.5*fr[13])*dvxr; 
  incr[15] = (3.0*fr[15]-1.732050807568877*fr[14])*wxr+(0.8660254037844386*fr[13]-0.5*fr[10])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  } 
} 
double vlasov_surf_1x3v_vx_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double dv3 = dxvr[3], wv3 = wr[3]; 
  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double Ghat[8]; 
  double favg[8]; 
  double alpha[8]; 

  favg[0] = (-1.224744871391589*fr[2])+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[3] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[5] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[7] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 0.5773502691896258*B2[0]*dv2; 
  alpha[3] = -0.5773502691896258*B1[0]*dv3; 
  alpha[4] = 0.5773502691896258*B2[1]*dv2; 
  alpha[5] = -0.5773502691896258*B1[1]*dv3; 

  const double amid = 0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[2]+fl[2])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[1]+fl[1])*amax+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[3]+fl[3])*amax+0.1767766952966368*(alpha[5]*favg[7]+alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[4]+fl[4])*amax+0.1767766952966368*(alpha[4]*favg[7]+alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[6]+fl[6])*amax+0.1767766952966368*(alpha[3]*favg[7]+alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[8]+fl[8])*amax+0.1767766952966368*(alpha[2]*favg[7]+alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[10]+fl[10])*amax+0.1767766952966368*(alpha[1]*favg[7]+alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[13]+fl[13])*amax+0.1767766952966368*(alpha[0]*favg[7]+alpha[1]*favg[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[6] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[9] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[10] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[11] += -1.224744871391589*Ghat[4]*dv10r; 
  outr[12] += -1.224744871391589*Ghat[5]*dv10r; 
  outr[13] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[14] += -1.224744871391589*Ghat[6]*dv10r; 
  outr[15] += -1.224744871391589*Ghat[7]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[6] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[9] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[10] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[11] += -1.224744871391589*Ghat[4]*dv10l; 
  outl[12] += -1.224744871391589*Ghat[5]*dv10l; 
  outl[13] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[14] += -1.224744871391589*Ghat[6]*dv10l; 
  outl[15] += -1.224744871391589*Ghat[7]*dv10l; 

  return fabs(amid); 
} 
double vlasov_surf_1x3v_vy_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[2]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double dv3 = dxvr[3], wv3 = wr[3]; 
  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double Ghat[8]; 
  double favg[8]; 
  double alpha[8]; 

  favg[0] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[6])+1.224744871391589*fl[6]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[3] = 0.5773502691896258*B0[0]*dv3; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  alpha[5] = 0.5773502691896258*B0[1]*dv3; 

  const double amid = 0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[6]+fl[6])-1.0*fr[1]+fl[1])*amax+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[2]+fl[2])*amax+0.1767766952966368*(alpha[5]*favg[7]+alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[4]+fl[4])*amax+0.1767766952966368*(alpha[4]*favg[7]+alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[5]+fl[5])*amax+0.1767766952966368*(alpha[3]*favg[7]+alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[8]+fl[8])*amax+0.1767766952966368*(alpha[2]*favg[7]+alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[9]+fl[9])*amax+0.1767766952966368*(alpha[1]*favg[7]+alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[12]+fl[12])*amax+0.1767766952966368*(alpha[0]*favg[7]+alpha[1]*favg[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[6] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[7] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[9] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[11] += -1.224744871391589*Ghat[4]*dv11r; 
  outr[12] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[13] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[14] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[15] += -1.224744871391589*Ghat[7]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[6] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[7] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[9] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[11] += -1.224744871391589*Ghat[4]*dv11l; 
  outl[12] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[13] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[14] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[15] += -1.224744871391589*Ghat[7]*dv11l; 

  return fabs(amid); 
} 
double vlasov_surf_1x3v_vz_p1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double* restrict outl, double* restrict outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[3]; 
  double dv12r = 2/dxvr[3]; 
  const double *E2 = &EM[4]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double dv3 = dxvr[3], wv3 = wr[3]; 
  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double Ghat[8]; 
  double favg[8]; 
  double alpha[8]; 

  favg[0] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[8])+1.224744871391589*fl[8]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 0.5773502691896258*B1[0]*dv1; 
  alpha[3] = -0.5773502691896258*B0[0]*dv2; 
  alpha[4] = 0.5773502691896258*B1[1]*dv1; 
  alpha[5] = -0.5773502691896258*B0[1]*dv2; 

  const double amid = 0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[1]+fl[1])*amax+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[2]+fl[2])*amax+0.1767766952966368*(alpha[5]*favg[7]+alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[3]+fl[3])*amax+0.1767766952966368*(alpha[4]*favg[7]+alpha[2]*favg[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[5]+fl[5])*amax+0.1767766952966368*(alpha[3]*favg[7]+alpha[5]*favg[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[6]+fl[6])*amax+0.1767766952966368*(alpha[2]*favg[7]+alpha[4]*favg[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[7]+fl[7])*amax+0.1767766952966368*(alpha[1]*favg[7]+alpha[0]*favg[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[11]+fl[11])*amax+0.1767766952966368*(alpha[0]*favg[7]+alpha[1]*favg[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv12r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv12r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv12r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv12r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv12r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv12r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv12r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv12r; 
  outr[8] += -1.224744871391589*Ghat[1]*dv12r; 
  outr[9] += -1.224744871391589*Ghat[2]*dv12r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv12r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv12r; 
  outr[12] += -1.224744871391589*Ghat[4]*dv12r; 
  outr[13] += -1.224744871391589*Ghat[5]*dv12r; 
  outr[14] += -1.224744871391589*Ghat[6]*dv12r; 
  outr[15] += -1.224744871391589*Ghat[7]*dv12r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv12l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv12l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv12l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv12l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv12l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv12l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv12l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv12l; 
  outl[8] += -1.224744871391589*Ghat[1]*dv12l; 
  outl[9] += -1.224744871391589*Ghat[2]*dv12l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv12l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv12l; 
  outl[12] += -1.224744871391589*Ghat[4]*dv12l; 
  outl[13] += -1.224744871391589*Ghat[5]*dv12l; 
  outl[14] += -1.224744871391589*Ghat[6]*dv12l; 
  outl[15] += -1.224744871391589*Ghat[7]*dv12l; 

  return fabs(amid); 
} 
