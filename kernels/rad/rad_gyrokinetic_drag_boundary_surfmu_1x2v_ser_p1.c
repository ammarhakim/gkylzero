#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // nI:        ion density. 
  // nuField:       2/pi*v*nu(v) field dg representation (v''(v||,mu) in notes
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  double rdv2 = 2.0/dxv[2]; 
  double alphaDrSurf[6] = {0.0}; 
  double Ghat[6] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(1.7320508075688772*(nI[1]*nuField[5]+nI[0]*nuField[3])+nI[1]*nuField[1]+nI[0]*nuField[0]); 
  alphaDrSurf[1] = 0.5*(1.7320508075688772*(nI[0]*nuField[5]+nI[1]*nuField[3])+nI[0]*nuField[1]+nuField[0]*nI[1]); 
  alphaDrSurf[2] = 0.5*(1.7320508075688772*(nI[1]*nuField[7]+nI[0]*nuField[6])+nI[1]*nuField[4]+nI[0]*nuField[2]); 
  alphaDrSurf[3] = 0.5*(1.7320508075688772*(nI[0]*nuField[7]+nI[1]*nuField[6])+nI[0]*nuField[4]+nI[1]*nuField[2]); 
  alphaDrSurf[4] = 0.03333333333333333*(25.980762113533157*nI[1]*nuField[11]+25.98076211353316*nI[0]*nuField[10]+15.000000000000002*nI[1]*nuField[9]+15.0*nI[0]*nuField[8]); 
  alphaDrSurf[5] = 0.03333333333333333*(25.98076211353316*nI[0]*nuField[11]+25.980762113533157*nI[1]*nuField[10]+15.0*nI[0]*nuField[9]+15.000000000000002*nI[1]*nuField[8]); 

  Ghat[0] = -(0.05*(12.247448713915892*(alphaDrSurf[5]*fEdge[11]+alphaDrSurf[4]*fEdge[10])-7.0710678118654755*(alphaDrSurf[5]*fEdge[9]+alphaDrSurf[4]*fEdge[8])+12.24744871391589*(alphaDrSurf[3]*fEdge[7]+alphaDrSurf[2]*fEdge[6]+alphaDrSurf[1]*fEdge[5])-7.0710678118654755*alphaDrSurf[3]*fEdge[4]+12.24744871391589*alphaDrSurf[0]*fEdge[3]-7.0710678118654755*(alphaDrSurf[2]*fEdge[2]+alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0]))); 
  Ghat[1] = -(0.016666666666666666*(36.74234614174767*(alphaDrSurf[4]*fEdge[11]+alphaDrSurf[5]*fEdge[10])-21.21320343559643*(alphaDrSurf[4]*fEdge[9]+alphaDrSurf[5]*fEdge[8])+36.74234614174767*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[3]*fEdge[6]+alphaDrSurf[0]*fEdge[5])-21.213203435596427*alphaDrSurf[2]*fEdge[4]+36.74234614174767*alphaDrSurf[1]*fEdge[3]-21.213203435596427*(fEdge[2]*alphaDrSurf[3]+alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1]))); 
  Ghat[2] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf[3]*fEdge[11]+alphaDrSurf[2]*fEdge[10])-18.97366596101028*(alphaDrSurf[3]*fEdge[9]+alphaDrSurf[2]*fEdge[8])+(32.86335345030997*alphaDrSurf[5]+36.74234614174767*alphaDrSurf[1])*fEdge[7]+(32.86335345030997*alphaDrSurf[4]+36.74234614174767*alphaDrSurf[0])*fEdge[6]+36.74234614174767*alphaDrSurf[3]*fEdge[5]+fEdge[4]*(-(18.97366596101028*alphaDrSurf[5])-21.213203435596427*alphaDrSurf[1])-18.97366596101028*fEdge[2]*alphaDrSurf[4]+36.74234614174767*alphaDrSurf[2]*fEdge[3]-21.213203435596427*(fEdge[1]*alphaDrSurf[3]+alphaDrSurf[0]*fEdge[2]+fEdge[0]*alphaDrSurf[2]))); 
  Ghat[3] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf[2]*fEdge[11]+alphaDrSurf[3]*fEdge[10])-18.97366596101028*(alphaDrSurf[2]*fEdge[9]+alphaDrSurf[3]*fEdge[8])+(32.86335345030997*alphaDrSurf[4]+36.74234614174767*alphaDrSurf[0])*fEdge[7]+(32.86335345030997*alphaDrSurf[5]+36.74234614174767*alphaDrSurf[1])*fEdge[6]+36.74234614174767*alphaDrSurf[2]*fEdge[5]-18.97366596101028*fEdge[2]*alphaDrSurf[5]+(-(18.97366596101028*alphaDrSurf[4])-21.213203435596427*alphaDrSurf[0])*fEdge[4]+36.74234614174767*alphaDrSurf[3]*fEdge[3]-21.213203435596427*(fEdge[0]*alphaDrSurf[3]+alphaDrSurf[1]*fEdge[2]+fEdge[1]*alphaDrSurf[2]))); 
  Ghat[4] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf[5]+257.1964229922337*alphaDrSurf[1])*fEdge[11]+(164.3167672515499*alphaDrSurf[4]+257.19642299223375*alphaDrSurf[0])*fEdge[10]+(-(94.86832980505142*alphaDrSurf[5])-148.49242404917499*alphaDrSurf[1])*fEdge[9]+(-(94.86832980505142*alphaDrSurf[4])-148.49242404917499*alphaDrSurf[0])*fEdge[8]+230.04347415216978*(alphaDrSurf[3]*fEdge[7]+alphaDrSurf[2]*fEdge[6])+alphaDrSurf[5]*(257.19642299223375*fEdge[5]-148.49242404917499*fEdge[1])-132.81566172707196*alphaDrSurf[3]*fEdge[4]+(257.1964229922337*fEdge[3]-148.49242404917499*fEdge[0])*alphaDrSurf[4]-132.81566172707196*alphaDrSurf[2]*fEdge[2])); 
  Ghat[5] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf[4]+257.19642299223375*alphaDrSurf[0])*fEdge[11]+(164.3167672515499*alphaDrSurf[5]+257.1964229922337*alphaDrSurf[1])*fEdge[10]+(-(94.86832980505142*alphaDrSurf[4])-148.49242404917499*alphaDrSurf[0])*fEdge[9]+(-(94.86832980505142*alphaDrSurf[5])-148.49242404917499*alphaDrSurf[1])*fEdge[8]+230.0434741521698*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[3]*fEdge[6])+257.19642299223375*alphaDrSurf[4]*fEdge[5]+(257.1964229922337*fEdge[3]-148.49242404917499*fEdge[0])*alphaDrSurf[5]-132.81566172707196*alphaDrSurf[2]*fEdge[4]-148.49242404917499*fEdge[1]*alphaDrSurf[4]-132.81566172707196*fEdge[2]*alphaDrSurf[3])); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[10] += 1.224744871391589*Ghat[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat[5]*rdv2; 

  } else { 

  alphaDrSurf[0] = -(0.5*(1.7320508075688772*(nI[1]*nuField[5]+nI[0]*nuField[3])-1.0*(nI[1]*nuField[1]+nI[0]*nuField[0]))); 
  alphaDrSurf[1] = -(0.5*(1.7320508075688772*(nI[0]*nuField[5]+nI[1]*nuField[3])-1.0*(nI[0]*nuField[1]+nuField[0]*nI[1]))); 
  alphaDrSurf[2] = -(0.5*(1.7320508075688772*(nI[1]*nuField[7]+nI[0]*nuField[6])-1.0*(nI[1]*nuField[4]+nI[0]*nuField[2]))); 
  alphaDrSurf[3] = -(0.5*(1.7320508075688772*(nI[0]*nuField[7]+nI[1]*nuField[6])-1.0*(nI[0]*nuField[4]+nI[1]*nuField[2]))); 
  alphaDrSurf[4] = -(0.03333333333333333*(25.980762113533157*nI[1]*nuField[11]+25.98076211353316*nI[0]*nuField[10]-15.000000000000002*nI[1]*nuField[9]-15.0*nI[0]*nuField[8])); 
  alphaDrSurf[5] = -(0.03333333333333333*(25.98076211353316*nI[0]*nuField[11]+25.980762113533157*nI[1]*nuField[10]-15.0*nI[0]*nuField[9]-15.000000000000002*nI[1]*nuField[8])); 

  Ghat[0] = -(0.05*(12.247448713915892*(alphaDrSurf[5]*fSkin[11]+alphaDrSurf[4]*fSkin[10])-7.0710678118654755*(alphaDrSurf[5]*fSkin[9]+alphaDrSurf[4]*fSkin[8])+12.24744871391589*(alphaDrSurf[3]*fSkin[7]+alphaDrSurf[2]*fSkin[6]+alphaDrSurf[1]*fSkin[5])-7.0710678118654755*alphaDrSurf[3]*fSkin[4]+12.24744871391589*alphaDrSurf[0]*fSkin[3]-7.0710678118654755*(alphaDrSurf[2]*fSkin[2]+alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0]))); 
  Ghat[1] = -(0.016666666666666666*(36.74234614174767*(alphaDrSurf[4]*fSkin[11]+alphaDrSurf[5]*fSkin[10])-21.21320343559643*(alphaDrSurf[4]*fSkin[9]+alphaDrSurf[5]*fSkin[8])+36.74234614174767*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[3]*fSkin[6]+alphaDrSurf[0]*fSkin[5])-21.213203435596427*alphaDrSurf[2]*fSkin[4]+36.74234614174767*alphaDrSurf[1]*fSkin[3]-21.213203435596427*(fSkin[2]*alphaDrSurf[3]+alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1]))); 
  Ghat[2] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf[3]*fSkin[11]+alphaDrSurf[2]*fSkin[10])-18.97366596101028*(alphaDrSurf[3]*fSkin[9]+alphaDrSurf[2]*fSkin[8])+(32.86335345030997*alphaDrSurf[5]+36.74234614174767*alphaDrSurf[1])*fSkin[7]+(32.86335345030997*alphaDrSurf[4]+36.74234614174767*alphaDrSurf[0])*fSkin[6]+36.74234614174767*alphaDrSurf[3]*fSkin[5]+fSkin[4]*(-(18.97366596101028*alphaDrSurf[5])-21.213203435596427*alphaDrSurf[1])-18.97366596101028*fSkin[2]*alphaDrSurf[4]+36.74234614174767*alphaDrSurf[2]*fSkin[3]-21.213203435596427*(fSkin[1]*alphaDrSurf[3]+alphaDrSurf[0]*fSkin[2]+fSkin[0]*alphaDrSurf[2]))); 
  Ghat[3] = -(0.016666666666666666*(32.86335345030997*(alphaDrSurf[2]*fSkin[11]+alphaDrSurf[3]*fSkin[10])-18.97366596101028*(alphaDrSurf[2]*fSkin[9]+alphaDrSurf[3]*fSkin[8])+(32.86335345030997*alphaDrSurf[4]+36.74234614174767*alphaDrSurf[0])*fSkin[7]+(32.86335345030997*alphaDrSurf[5]+36.74234614174767*alphaDrSurf[1])*fSkin[6]+36.74234614174767*alphaDrSurf[2]*fSkin[5]-18.97366596101028*fSkin[2]*alphaDrSurf[5]+(-(18.97366596101028*alphaDrSurf[4])-21.213203435596427*alphaDrSurf[0])*fSkin[4]+36.74234614174767*alphaDrSurf[3]*fSkin[3]-21.213203435596427*(fSkin[0]*alphaDrSurf[3]+alphaDrSurf[1]*fSkin[2]+fSkin[1]*alphaDrSurf[2]))); 
  Ghat[4] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf[5]+257.1964229922337*alphaDrSurf[1])*fSkin[11]+(164.3167672515499*alphaDrSurf[4]+257.19642299223375*alphaDrSurf[0])*fSkin[10]+(-(94.86832980505142*alphaDrSurf[5])-148.49242404917499*alphaDrSurf[1])*fSkin[9]+(-(94.86832980505142*alphaDrSurf[4])-148.49242404917499*alphaDrSurf[0])*fSkin[8]+230.04347415216978*(alphaDrSurf[3]*fSkin[7]+alphaDrSurf[2]*fSkin[6])+alphaDrSurf[5]*(257.19642299223375*fSkin[5]-148.49242404917499*fSkin[1])-132.81566172707196*alphaDrSurf[3]*fSkin[4]+(257.1964229922337*fSkin[3]-148.49242404917499*fSkin[0])*alphaDrSurf[4]-132.81566172707196*alphaDrSurf[2]*fSkin[2])); 
  Ghat[5] = -(0.002380952380952381*((164.3167672515499*alphaDrSurf[4]+257.19642299223375*alphaDrSurf[0])*fSkin[11]+(164.3167672515499*alphaDrSurf[5]+257.1964229922337*alphaDrSurf[1])*fSkin[10]+(-(94.86832980505142*alphaDrSurf[4])-148.49242404917499*alphaDrSurf[0])*fSkin[9]+(-(94.86832980505142*alphaDrSurf[5])-148.49242404917499*alphaDrSurf[1])*fSkin[8]+230.0434741521698*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[3]*fSkin[6])+257.19642299223375*alphaDrSurf[4]*fSkin[5]+(257.1964229922337*fSkin[3]-148.49242404917499*fSkin[0])*alphaDrSurf[5]-132.81566172707196*alphaDrSurf[2]*fSkin[4]-148.49242404917499*fSkin[1]*alphaDrSurf[4]-132.81566172707196*fSkin[2]*alphaDrSurf[3])); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += -(0.7071067811865475*Ghat[2]*rdv2); 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -(0.7071067811865475*Ghat[3]*rdv2); 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += -(0.7071067811865475*Ghat[4]*rdv2); 
  out[9] += -(0.7071067811865475*Ghat[5]*rdv2); 
  out[10] += 1.224744871391589*Ghat[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat[5]*rdv2; 

  } 

  return 0.;

} 
