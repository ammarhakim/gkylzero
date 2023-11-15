#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:    cell-center coordinates. 
  // dxv[3]:  cell spacing. 
  // nI:       ion density 
  // nuField:    sqrt(mu*me/2B)*v^2*nu(v) field dg represenation (v''(v||,mu) in notes
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double rdv2 = 2.0/dxv[2]; 

  double Ghat_r[6] = {0.0}; 
  double Ghat_l[6] = {0.0}; 
  double alphaDrSurf_l[6] = {0.0}; 
  alphaDrSurf_l[0] = -(0.8660254037844386*nI[1]*nuField[5])-0.8660254037844386*nI[0]*nuField[3]+0.5*nI[1]*nuField[1]+0.5*nI[0]*nuField[0]; 
  alphaDrSurf_l[1] = -(0.8660254037844386*nI[0]*nuField[5])-0.8660254037844386*nI[1]*nuField[3]+0.5*nI[0]*nuField[1]+0.5*nuField[0]*nI[1]; 
  alphaDrSurf_l[2] = -(0.8660254037844386*nI[1]*nuField[7])-0.8660254037844386*nI[0]*nuField[6]+0.5*nI[1]*nuField[4]+0.5*nI[0]*nuField[2]; 
  alphaDrSurf_l[3] = -(0.8660254037844386*nI[0]*nuField[7])-0.8660254037844386*nI[1]*nuField[6]+0.5*nI[0]*nuField[4]+0.5*nI[1]*nuField[2]; 
  alphaDrSurf_l[4] = -(0.8660254037844386*nI[1]*nuField[11])-0.8660254037844387*nI[0]*nuField[10]+0.5000000000000001*nI[1]*nuField[9]+0.5*nI[0]*nuField[8]; 
  alphaDrSurf_l[5] = -(0.8660254037844387*nI[0]*nuField[11])-0.8660254037844386*nI[1]*nuField[10]+0.5*nI[0]*nuField[9]+0.5000000000000001*nI[1]*nuField[8]; 

  double alphaDrSurf_r[6] = {0.0}; 
  alphaDrSurf_r[0] = 0.8660254037844386*nI[1]*nuField[5]+0.8660254037844386*nI[0]*nuField[3]+0.5*nI[1]*nuField[1]+0.5*nI[0]*nuField[0]; 
  alphaDrSurf_r[1] = 0.8660254037844386*nI[0]*nuField[5]+0.8660254037844386*nI[1]*nuField[3]+0.5*nI[0]*nuField[1]+0.5*nuField[0]*nI[1]; 
  alphaDrSurf_r[2] = 0.8660254037844386*nI[1]*nuField[7]+0.8660254037844386*nI[0]*nuField[6]+0.5*nI[1]*nuField[4]+0.5*nI[0]*nuField[2]; 
  alphaDrSurf_r[3] = 0.8660254037844386*nI[0]*nuField[7]+0.8660254037844386*nI[1]*nuField[6]+0.5*nI[0]*nuField[4]+0.5*nI[1]*nuField[2]; 
  alphaDrSurf_r[4] = 0.8660254037844386*nI[1]*nuField[11]+0.8660254037844387*nI[0]*nuField[10]+0.5000000000000001*nI[1]*nuField[9]+0.5*nI[0]*nuField[8]; 
  alphaDrSurf_r[5] = 0.8660254037844387*nI[0]*nuField[11]+0.8660254037844386*nI[1]*nuField[10]+0.5*nI[0]*nuField[9]+0.5000000000000001*nI[1]*nuField[8]; 

  Ghat_l[0] = -(0.05*(12.247448713915892*(alphaDrSurf_l[5]*fc[11]+alphaDrSurf_l[4]*fc[10])-7.0710678118654755*(alphaDrSurf_l[5]*fc[9]+alphaDrSurf_l[4]*fc[8])+12.24744871391589*(alphaDrSurf_l[3]*fc[7]+alphaDrSurf_l[2]*fc[6]+alphaDrSurf_l[1]*fc[5])-7.0710678118654755*alphaDrSurf_l[3]*fc[4]+12.24744871391589*alphaDrSurf_l[0]*fc[3]-7.0710678118654755*(alphaDrSurf_l[2]*fc[2]+alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0]))); 
  Ghat_l[1] = -(0.016666666666666666*(36.74234614174767*(alphaDrSurf_l[4]*fc[11]+alphaDrSurf_l[5]*fc[10])-21.21320343559643*(alphaDrSurf_l[4]*fc[9]+alphaDrSurf_l[5]*fc[8])+36.74234614174767*(alphaDrSurf_l[2]*fc[7]+alphaDrSurf_l[3]*fc[6]+alphaDrSurf_l[0]*fc[5])-21.213203435596427*alphaDrSurf_l[2]*fc[4]+36.74234614174767*alphaDrSurf_l[1]*fc[3]-21.213203435596427*(fc[2]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1]))); 
  Ghat_l[2] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf_l[3]*fc[11]+alphaDrSurf_l[2]*fc[10])-18.97366596101028*(alphaDrSurf_l[3]*fc[9]+alphaDrSurf_l[2]*fc[8])+(32.86335345030997*alphaDrSurf_l[5]+36.74234614174767*alphaDrSurf_l[1])*fc[7]+(32.86335345030997*alphaDrSurf_l[4]+36.74234614174767*alphaDrSurf_l[0])*fc[6]+36.74234614174767*alphaDrSurf_l[3]*fc[5]+fc[4]*(-(18.97366596101028*alphaDrSurf_l[5])-21.213203435596427*alphaDrSurf_l[1])-18.97366596101028*fc[2]*alphaDrSurf_l[4]+36.74234614174767*alphaDrSurf_l[2]*fc[3]-21.213203435596427*(fc[1]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*fc[2]+fc[0]*alphaDrSurf_l[2]))); 
  Ghat_l[3] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf_l[2]*fc[11]+alphaDrSurf_l[3]*fc[10])-18.97366596101028*(alphaDrSurf_l[2]*fc[9]+alphaDrSurf_l[3]*fc[8])+(32.86335345030997*alphaDrSurf_l[4]+36.74234614174767*alphaDrSurf_l[0])*fc[7]+(32.86335345030997*alphaDrSurf_l[5]+36.74234614174767*alphaDrSurf_l[1])*fc[6]+36.74234614174767*alphaDrSurf_l[2]*fc[5]-18.97366596101028*fc[2]*alphaDrSurf_l[5]+(-(18.97366596101028*alphaDrSurf_l[4])-21.213203435596427*alphaDrSurf_l[0])*fc[4]+36.74234614174767*alphaDrSurf_l[3]*fc[3]-21.213203435596427*(fc[0]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*fc[2]+fc[1]*alphaDrSurf_l[2]))); 
  Ghat_l[4] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf_l[5]+257.1964229922337*alphaDrSurf_l[1])*fc[11]+(164.3167672515499*alphaDrSurf_l[4]+257.19642299223375*alphaDrSurf_l[0])*fc[10]+(-(94.86832980505142*alphaDrSurf_l[5])-148.49242404917499*alphaDrSurf_l[1])*fc[9]+(-(94.86832980505142*alphaDrSurf_l[4])-148.49242404917499*alphaDrSurf_l[0])*fc[8]+230.04347415216978*(alphaDrSurf_l[3]*fc[7]+alphaDrSurf_l[2]*fc[6])+alphaDrSurf_l[5]*(257.19642299223375*fc[5]-148.49242404917499*fc[1])-132.81566172707196*alphaDrSurf_l[3]*fc[4]+(257.1964229922337*fc[3]-148.49242404917499*fc[0])*alphaDrSurf_l[4]-132.81566172707196*alphaDrSurf_l[2]*fc[2])); 
  Ghat_l[5] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf_l[4]+257.19642299223375*alphaDrSurf_l[0])*fc[11]+(164.3167672515499*alphaDrSurf_l[5]+257.1964229922337*alphaDrSurf_l[1])*fc[10]+(-(94.86832980505142*alphaDrSurf_l[4])-148.49242404917499*alphaDrSurf_l[0])*fc[9]+(-(94.86832980505142*alphaDrSurf_l[5])-148.49242404917499*alphaDrSurf_l[1])*fc[8]+230.0434741521698*(alphaDrSurf_l[2]*fc[7]+alphaDrSurf_l[3]*fc[6])+257.19642299223375*alphaDrSurf_l[4]*fc[5]+(257.1964229922337*fc[3]-148.49242404917499*fc[0])*alphaDrSurf_l[5]-132.81566172707196*alphaDrSurf_l[2]*fc[4]-148.49242404917499*fc[1]*alphaDrSurf_l[4]-132.81566172707196*fc[2]*alphaDrSurf_l[3])); 

  Ghat_r[0] = -(0.05*(12.247448713915892*(alphaDrSurf_r[5]*fr[11]+alphaDrSurf_r[4]*fr[10])-7.0710678118654755*(alphaDrSurf_r[5]*fr[9]+alphaDrSurf_r[4]*fr[8])+12.24744871391589*(alphaDrSurf_r[3]*fr[7]+alphaDrSurf_r[2]*fr[6]+alphaDrSurf_r[1]*fr[5])-7.0710678118654755*alphaDrSurf_r[3]*fr[4]+12.24744871391589*alphaDrSurf_r[0]*fr[3]-7.0710678118654755*(alphaDrSurf_r[2]*fr[2]+alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0]))); 
  Ghat_r[1] = -(0.016666666666666666*(36.74234614174767*(alphaDrSurf_r[4]*fr[11]+alphaDrSurf_r[5]*fr[10])-21.21320343559643*(alphaDrSurf_r[4]*fr[9]+alphaDrSurf_r[5]*fr[8])+36.74234614174767*(alphaDrSurf_r[2]*fr[7]+alphaDrSurf_r[3]*fr[6]+alphaDrSurf_r[0]*fr[5])-21.213203435596427*alphaDrSurf_r[2]*fr[4]+36.74234614174767*alphaDrSurf_r[1]*fr[3]-21.213203435596427*(fr[2]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1]))); 
  Ghat_r[2] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf_r[3]*fr[11]+alphaDrSurf_r[2]*fr[10])-18.97366596101028*(alphaDrSurf_r[3]*fr[9]+alphaDrSurf_r[2]*fr[8])+(32.86335345030997*alphaDrSurf_r[5]+36.74234614174767*alphaDrSurf_r[1])*fr[7]+(32.86335345030997*alphaDrSurf_r[4]+36.74234614174767*alphaDrSurf_r[0])*fr[6]+36.74234614174767*alphaDrSurf_r[3]*fr[5]+fr[4]*(-(18.97366596101028*alphaDrSurf_r[5])-21.213203435596427*alphaDrSurf_r[1])-18.97366596101028*fr[2]*alphaDrSurf_r[4]+36.74234614174767*alphaDrSurf_r[2]*fr[3]-21.213203435596427*(fr[1]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*fr[2]+fr[0]*alphaDrSurf_r[2]))); 
  Ghat_r[3] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf_r[2]*fr[11]+alphaDrSurf_r[3]*fr[10])-18.97366596101028*(alphaDrSurf_r[2]*fr[9]+alphaDrSurf_r[3]*fr[8])+(32.86335345030997*alphaDrSurf_r[4]+36.74234614174767*alphaDrSurf_r[0])*fr[7]+(32.86335345030997*alphaDrSurf_r[5]+36.74234614174767*alphaDrSurf_r[1])*fr[6]+36.74234614174767*alphaDrSurf_r[2]*fr[5]-18.97366596101028*fr[2]*alphaDrSurf_r[5]+(-(18.97366596101028*alphaDrSurf_r[4])-21.213203435596427*alphaDrSurf_r[0])*fr[4]+36.74234614174767*alphaDrSurf_r[3]*fr[3]-21.213203435596427*(fr[0]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*fr[2]+fr[1]*alphaDrSurf_r[2]))); 
  Ghat_r[4] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf_r[5]+257.1964229922337*alphaDrSurf_r[1])*fr[11]+(164.3167672515499*alphaDrSurf_r[4]+257.19642299223375*alphaDrSurf_r[0])*fr[10]+(-(94.86832980505142*alphaDrSurf_r[5])-148.49242404917499*alphaDrSurf_r[1])*fr[9]+(-(94.86832980505142*alphaDrSurf_r[4])-148.49242404917499*alphaDrSurf_r[0])*fr[8]+230.04347415216978*(alphaDrSurf_r[3]*fr[7]+alphaDrSurf_r[2]*fr[6])+alphaDrSurf_r[5]*(257.19642299223375*fr[5]-148.49242404917499*fr[1])-132.81566172707196*alphaDrSurf_r[3]*fr[4]+(257.1964229922337*fr[3]-148.49242404917499*fr[0])*alphaDrSurf_r[4]-132.81566172707196*alphaDrSurf_r[2]*fr[2])); 
  Ghat_r[5] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf_r[4]+257.19642299223375*alphaDrSurf_r[0])*fr[11]+(164.3167672515499*alphaDrSurf_r[5]+257.1964229922337*alphaDrSurf_r[1])*fr[10]+(-(94.86832980505142*alphaDrSurf_r[4])-148.49242404917499*alphaDrSurf_r[0])*fr[9]+(-(94.86832980505142*alphaDrSurf_r[5])-148.49242404917499*alphaDrSurf_r[1])*fr[8]+230.0434741521698*(alphaDrSurf_r[2]*fr[7]+alphaDrSurf_r[3]*fr[6])+257.19642299223375*alphaDrSurf_r[4]*fr[5]+(257.1964229922337*fr[3]-148.49242404917499*fr[0])*alphaDrSurf_r[5]-132.81566172707196*alphaDrSurf_r[2]*fr[4]-148.49242404917499*fr[1]*alphaDrSurf_r[4]-132.81566172707196*fr[2]*alphaDrSurf_r[3])); 

  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*rdv2; 
  out[3] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[4] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*rdv2; 
  out[5] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[7] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[8] += (0.7071067811865475*Ghat_r[4]-0.7071067811865475*Ghat_l[4])*rdv2; 
  out[9] += (0.7071067811865475*Ghat_r[5]-0.7071067811865475*Ghat_l[5])*rdv2; 
  out[10] += 1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[11] += 1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 

  return 0.;

} 
