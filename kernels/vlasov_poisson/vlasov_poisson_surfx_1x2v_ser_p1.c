#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_surfx_1x2v_ser_p1(const double *w, const double *dxv, 
    const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[8]; 
  double Ghat_l[8]; 

  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[4]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[4]+0.7071067811865475*fc[2])*wv+(0.3162277660168379*fc[9]+0.1825741858350554*fc[8]+0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[5]+0.7071067811865475*fc[3])*wv+(0.3535533905932737*fc[7]+0.2041241452319315*fc[6])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[7]+0.7071067811865475*fc[6])*wv+(0.3162277660168379*fc[11]+0.1825741858350554*fc[10]+0.3535533905932737*fc[5]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[9]+0.7071067811865475*fc[8])*wv+(0.3162277660168379*fc[4]+0.1825741858350554*fc[2])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[13]+0.7071067811865475*fc[12])*wv+(0.3535533905932737*fc[15]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[11]+0.7071067811865475*fc[10])*wv+(0.3162277660168379*fc[7]+0.1825741858350554*fc[6])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[15]+0.7071067811865475*fc[14])*wv+(0.3535533905932737*fc[13]+0.2041241452319315*fc[12])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[4]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[4]+0.7071067811865475*fl[2])*wv+(0.3162277660168379*fl[9]+0.1825741858350554*fl[8]+0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[5]+0.7071067811865475*fl[3])*wv+(0.3535533905932737*fl[7]+0.2041241452319315*fl[6])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[7]+0.7071067811865475*fl[6])*wv+(0.3162277660168379*fl[11]+0.1825741858350554*fl[10]+0.3535533905932737*fl[5]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[9]+0.7071067811865475*fl[8])*wv+(0.3162277660168379*fl[4]+0.1825741858350554*fl[2])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[13]+0.7071067811865475*fl[12])*wv+(0.3535533905932737*fl[15]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[11]+0.7071067811865475*fl[10])*wv+(0.3162277660168379*fl[7]+0.1825741858350554*fl[6])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[15]+0.7071067811865475*fl[14])*wv+(0.3535533905932737*fl[13]+0.2041241452319315*fl[12])*dv; 

  } else { 

  Ghat_r[0] = -0.08333333333333333*((14.69693845669907*fr[1]-8.485281374238571*fr[0])*wv+(4.242640687119286*fr[4]-2.449489742783178*fr[2])*dv); 
  Ghat_r[1] = -0.01666666666666667*((73.48469228349535*fr[4]-42.42640687119286*fr[2])*wv+(18.97366596101028*fr[9]-10.95445115010332*fr[8]+21.21320343559643*fr[1]-12.24744871391589*fr[0])*dv); 
  Ghat_r[2] = -0.08333333333333333*((14.69693845669907*fr[5]-8.485281374238571*fr[3])*wv+(4.242640687119286*fr[7]-2.449489742783178*fr[6])*dv); 
  Ghat_r[3] = -0.01666666666666667*((73.48469228349535*fr[7]-42.42640687119286*fr[6])*wv+(18.97366596101028*fr[11]-10.95445115010333*fr[10]+21.21320343559643*fr[5]-12.24744871391589*fr[3])*dv); 
  Ghat_r[4] = -0.03333333333333333*((36.74234614174768*fr[9]-21.21320343559643*fr[8])*wv+(9.48683298050514*fr[4]-5.477225575051662*fr[2])*dv); 
  Ghat_r[5] = -0.01666666666666667*((73.48469228349536*fr[13]-42.42640687119286*fr[12])*wv+(21.21320343559643*fr[15]-12.24744871391589*fr[14])*dv); 
  Ghat_r[6] = -0.03333333333333333*((36.74234614174768*fr[11]-21.21320343559643*fr[10])*wv+(9.48683298050514*fr[7]-5.477225575051662*fr[6])*dv); 
  Ghat_r[7] = -0.01666666666666667*((73.48469228349536*fr[15]-42.42640687119286*fr[14])*wv+(21.21320343559643*fr[13]-12.24744871391589*fr[12])*dv); 

  Ghat_l[0] = -0.08333333333333333*((14.69693845669907*fc[1]-8.485281374238571*fc[0])*wv+(4.242640687119286*fc[4]-2.449489742783178*fc[2])*dv); 
  Ghat_l[1] = -0.01666666666666667*((73.48469228349535*fc[4]-42.42640687119286*fc[2])*wv+(18.97366596101028*fc[9]-10.95445115010332*fc[8]+21.21320343559643*fc[1]-12.24744871391589*fc[0])*dv); 
  Ghat_l[2] = -0.08333333333333333*((14.69693845669907*fc[5]-8.485281374238571*fc[3])*wv+(4.242640687119286*fc[7]-2.449489742783178*fc[6])*dv); 
  Ghat_l[3] = -0.01666666666666667*((73.48469228349535*fc[7]-42.42640687119286*fc[6])*wv+(18.97366596101028*fc[11]-10.95445115010333*fc[10]+21.21320343559643*fc[5]-12.24744871391589*fc[3])*dv); 
  Ghat_l[4] = -0.03333333333333333*((36.74234614174768*fc[9]-21.21320343559643*fc[8])*wv+(9.48683298050514*fc[4]-5.477225575051662*fc[2])*dv); 
  Ghat_l[5] = -0.01666666666666667*((73.48469228349536*fc[13]-42.42640687119286*fc[12])*wv+(21.21320343559643*fc[15]-12.24744871391589*fc[14])*dv); 
  Ghat_l[6] = -0.03333333333333333*((36.74234614174768*fc[11]-21.21320343559643*fc[10])*wv+(9.48683298050514*fc[7]-5.477225575051662*fc[6])*dv); 
  Ghat_l[7] = -0.01666666666666667*((73.48469228349536*fc[15]-42.42640687119286*fc[14])*wv+(21.21320343559643*fc[13]-12.24744871391589*fc[12])*dv); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[9] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[12] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 

  return 0.;

} 