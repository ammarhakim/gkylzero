#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
   double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvx = w[1];
  double rdvx2 = 2.0/dxv[1];
  double wvy = w[2];
  double rdvy2 = 2.0/dxv[2];
  double wvz = w[3];
  double rdvz2 = 2.0/dxv[3];

  double *alphaR = &alpha_surf[0];
  double *sgn_alpha_surfR = &sgn_alpha_surf[0];
  alphaR[0] = 2.121320343559642*hamil[5]*rdvx2+1.224744871391589*hamil[2]*rdvx2; 
  alphaR[1] = 4.743416490252569*hamil[17]*rdvx2+2.738612787525831*hamil[16]*rdvx2; 
  alphaR[2] = 2.121320343559642*hamil[11]*rdvx2+1.224744871391589*hamil[7]*rdvx2; 
  alphaR[3] = 2.121320343559642*hamil[12]*rdvx2+1.224744871391589*hamil[9]*rdvx2; 
  alphaR[4] = 4.743416490252569*hamil[20]*rdvx2+2.738612787525831*hamil[18]*rdvx2; 
  alphaR[5] = 4.743416490252569*hamil[21]*rdvx2+2.738612787525831*hamil[19]*rdvx2; 
  alphaR[6] = 2.121320343559642*hamil[15]*rdvx2+1.224744871391589*hamil[14]*rdvx2; 
  alphaR[8] = 2.121320343559642*hamil[28]*rdvx2+1.224744871391589*hamil[26]*rdvx2; 
  alphaR[9] = 2.121320343559642*hamil[36]*rdvx2+1.224744871391589*hamil[34]*rdvx2; 
  alphaR[10] = 4.743416490252569*hamil[23]*rdvx2+2.738612787525831*hamil[22]*rdvx2; 
  alphaR[14] = 2.121320343559642*hamil[31]*rdvx2+1.224744871391589*hamil[30]*rdvx2; 
  alphaR[16] = 2.121320343559642*hamil[39]*rdvx2+1.224744871391589*hamil[38]*rdvx2; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.4242640687119282*alphaR[16])-0.4242640687119286*alphaR[14]-0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])+0.6363961030678927*(alphaR[6]+alphaR[5]+alphaR[4])-0.4743416490252568*(alphaR[3]+alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (0.5303300858899102*alphaR[16]-0.3952847075210471*alphaR[9]+0.3162277660168378*alphaR[8]+0.6363961030678927*alphaR[4]-0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaR[16])+0.4242640687119286*alphaR[14]+0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])-0.6363961030678927*(alphaR[6]+alphaR[5])+0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[3]-0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[14]+0.3162277660168378*alphaR[9]-0.3952847075210471*alphaR[8]+0.6363961030678927*alphaR[5]-0.4743416490252568*(alphaR[3]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3952847075210471*(alphaR[9]+alphaR[8]))-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[14])+0.3162277660168378*alphaR[9]-0.3952847075210471*alphaR[8]-0.6363961030678927*alphaR[5]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaR[16]-0.4242640687119286*alphaR[14]+0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])-0.6363961030678927*alphaR[6]+0.6363961030678927*alphaR[5]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[16])-0.3952847075210471*alphaR[9]+0.3162277660168378*alphaR[8]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaR[16]+0.4242640687119286*alphaR[14]-0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])+0.6363961030678927*alphaR[6]-0.6363961030678927*(alphaR[5]+alphaR[4])+0.4743416490252568*(alphaR[3]+alphaR[2])-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaR[16])-0.4242640687119286*alphaR[14]+0.3162277660168378*(alphaR[9]+alphaR[8])+0.6363961030678927*alphaR[6]-0.4743416490252568*(alphaR[3]+alphaR[2])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[16]-0.3952847075210471*alphaR[9]+0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaR[16])+0.4242640687119286*alphaR[14]+0.3162277660168378*(alphaR[9]+alphaR[8])-0.6363961030678927*alphaR[6]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[14]+0.3162277660168378*alphaR[9]-0.3952847075210471*alphaR[8]-0.4743416490252568*alphaR[3]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR[0]-0.3952847075210471*(alphaR[9]+alphaR[8]) > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[14])+0.3162277660168378*alphaR[9]-0.3952847075210471*alphaR[8]+0.4743416490252568*alphaR[3]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaR[16]-0.4242640687119286*alphaR[14]+0.3162277660168378*(alphaR[9]+alphaR[8])-0.6363961030678927*alphaR[6]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[16])-0.3952847075210471*alphaR[9]+0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaR[16]+0.4242640687119286*alphaR[14]+0.3162277660168378*(alphaR[9]+alphaR[8])+0.6363961030678927*alphaR[6]+0.4743416490252568*(alphaR[3]+alphaR[2])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaR[16])-0.4242640687119286*alphaR[14]+0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])+0.6363961030678927*alphaR[6]-0.6363961030678927*(alphaR[5]+alphaR[4])-0.4743416490252568*(alphaR[3]+alphaR[2])+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[16]-0.3952847075210471*alphaR[9]+0.3162277660168378*alphaR[8]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119282*alphaR[16])+0.4242640687119286*alphaR[14]-0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])-0.6363961030678927*alphaR[6]+0.6363961030678927*alphaR[5]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[14]+0.3162277660168378*alphaR[9]-0.3952847075210471*alphaR[8]-0.6363961030678927*alphaR[5]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3952847075210471*(alphaR[9]+alphaR[8]))+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[14])+0.3162277660168378*alphaR[9]-0.3952847075210471*alphaR[8]+0.6363961030678927*alphaR[5]+0.4743416490252568*(alphaR[3]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaR[16]-0.4242640687119286*alphaR[14]-0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])-0.6363961030678927*(alphaR[6]+alphaR[5])+0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[3]+0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[16])-0.3952847075210471*alphaR[9]+0.3162277660168378*alphaR[8]+0.6363961030678927*alphaR[4]+0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119282*alphaR[16]+0.4242640687119286*alphaR[14]+0.8538149682454614*alphaR[10]+0.3162277660168378*(alphaR[9]+alphaR[8])+0.6363961030678927*(alphaR[6]+alphaR[5]+alphaR[4])+0.4743416490252568*(alphaR[3]+alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
