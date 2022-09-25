#include <gkyl_mom_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void mom_vlasov_pkpm_surfx_3x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *u_i, const double *bvar, const double *fl, const double *fr, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]/2.0; 
  const double wvpar = w[0], dvpar = dxv[0]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *ul = &u_i[0]; 
  const double *bl = &bvar[0]; 
  double *rho_flux = &out[0]; 
  double *heat_flux = &out[24]; 
  double u[20] = {0.0}; 
  double q[20] = {0.0}; 

  u[0] = 2.23606797749979*ul[7]+1.732050807568877*ul[1]+ul[0]; 
  u[1] = 2.23606797749979*ul[11]+1.732050807568877*ul[4]+ul[2]; 
  u[2] = 2.23606797749979*ul[13]+1.732050807568877*ul[5]+ul[3]; 
  u[4] = 2.23606797749979*ul[17]+1.732050807568877*ul[10]+ul[6]; 
  u[7] = 1.732050807568877*ul[12]+ul[8]; 
  u[8] = 1.732050807568877*ul[15]+ul[9]; 
  u[11] = 1.732050807568877*ul[18]+ul[14]; 
  u[12] = 1.732050807568877*ul[19]+ul[16]; 

  q[0] = 1.118033988749895*bl[7]*wvpar_cu+0.8660254037844386*bl[1]*wvpar_cu+0.5*bl[0]*wvpar_cu+0.2795084971874737*bl[7]*dvpar_sq*wvpar+0.2165063509461096*bl[1]*dvpar_sq*wvpar+0.125*bl[0]*dvpar_sq*wvpar; 
  q[1] = 1.118033988749895*bl[11]*wvpar_cu+0.8660254037844386*bl[4]*wvpar_cu+0.5*bl[2]*wvpar_cu+0.2795084971874738*bl[11]*dvpar_sq*wvpar+0.2165063509461096*bl[4]*dvpar_sq*wvpar+0.125*bl[2]*dvpar_sq*wvpar; 
  q[2] = 1.118033988749895*bl[13]*wvpar_cu+0.8660254037844386*bl[5]*wvpar_cu+0.5*bl[3]*wvpar_cu+0.2795084971874738*bl[13]*dvpar_sq*wvpar+0.2165063509461096*bl[5]*dvpar_sq*wvpar+0.125*bl[3]*dvpar_sq*wvpar; 
  q[3] = 0.9682458365518543*bl[7]*dvpar*wvpar_sq+0.75*bl[1]*dvpar*wvpar_sq+0.4330127018922193*bl[0]*dvpar*wvpar_sq+0.04841229182759271*bl[7]*dvpar_cu+0.0375*bl[1]*dvpar_cu+0.02165063509461097*bl[0]*dvpar_cu; 
  q[4] = 1.118033988749895*bl[17]*wvpar_cu+0.8660254037844386*bl[10]*wvpar_cu+0.5*bl[6]*wvpar_cu+0.2795084971874737*bl[17]*dvpar_sq*wvpar+0.2165063509461096*bl[10]*dvpar_sq*wvpar+0.125*bl[6]*dvpar_sq*wvpar; 
  q[5] = 0.9682458365518543*bl[11]*dvpar*wvpar_sq+0.75*bl[4]*dvpar*wvpar_sq+0.4330127018922193*bl[2]*dvpar*wvpar_sq+0.04841229182759271*bl[11]*dvpar_cu+0.0375*bl[4]*dvpar_cu+0.02165063509461097*bl[2]*dvpar_cu; 
  q[6] = 0.9682458365518543*bl[13]*dvpar*wvpar_sq+0.75*bl[5]*dvpar*wvpar_sq+0.4330127018922193*bl[3]*dvpar*wvpar_sq+0.04841229182759271*bl[13]*dvpar_cu+0.0375*bl[5]*dvpar_cu+0.02165063509461097*bl[3]*dvpar_cu; 
  q[7] = 0.8660254037844387*bl[12]*wvpar_cu+0.5*bl[8]*wvpar_cu+0.2165063509461097*bl[12]*dvpar_sq*wvpar+0.125*bl[8]*dvpar_sq*wvpar; 
  q[8] = 0.8660254037844387*bl[15]*wvpar_cu+0.5*bl[9]*wvpar_cu+0.2165063509461097*bl[15]*dvpar_sq*wvpar+0.125*bl[9]*dvpar_sq*wvpar; 
  q[9] = 0.25*bl[7]*dvpar_sq*wvpar+0.1936491673103708*bl[1]*dvpar_sq*wvpar+0.1118033988749895*bl[0]*dvpar_sq*wvpar; 
  q[10] = 0.9682458365518543*bl[17]*dvpar*wvpar_sq+0.75*bl[10]*dvpar*wvpar_sq+0.4330127018922193*bl[6]*dvpar*wvpar_sq+0.04841229182759271*bl[17]*dvpar_cu+0.0375*bl[10]*dvpar_cu+0.02165063509461097*bl[6]*dvpar_cu; 
  q[11] = 0.8660254037844387*bl[18]*wvpar_cu+0.5*bl[14]*wvpar_cu+0.2165063509461097*bl[18]*dvpar_sq*wvpar+0.125*bl[14]*dvpar_sq*wvpar; 
  q[12] = 0.8660254037844387*bl[19]*wvpar_cu+0.5*bl[16]*wvpar_cu+0.2165063509461097*bl[19]*dvpar_sq*wvpar+0.125*bl[16]*dvpar_sq*wvpar; 
  q[13] = 0.75*bl[12]*dvpar*wvpar_sq+0.4330127018922194*bl[8]*dvpar*wvpar_sq+0.0375*bl[12]*dvpar_cu+0.02165063509461096*bl[8]*dvpar_cu; 
  q[14] = 0.75*bl[15]*dvpar*wvpar_sq+0.4330127018922194*bl[9]*dvpar*wvpar_sq+0.0375*bl[15]*dvpar_cu+0.02165063509461096*bl[9]*dvpar_cu; 
  q[15] = 0.25*bl[11]*dvpar_sq*wvpar+0.1936491673103709*bl[4]*dvpar_sq*wvpar+0.1118033988749895*bl[2]*dvpar_sq*wvpar; 
  q[16] = 0.25*bl[13]*dvpar_sq*wvpar+0.1936491673103709*bl[5]*dvpar_sq*wvpar+0.1118033988749895*bl[3]*dvpar_sq*wvpar; 
  q[17] = 0.75*bl[18]*dvpar*wvpar_sq+0.4330127018922194*bl[14]*dvpar*wvpar_sq+0.0375*bl[18]*dvpar_cu+0.02165063509461096*bl[14]*dvpar_cu; 
  q[18] = 0.75*bl[19]*dvpar*wvpar_sq+0.4330127018922194*bl[16]*dvpar*wvpar_sq+0.0375*bl[19]*dvpar_cu+0.02165063509461096*bl[16]*dvpar_cu; 
  q[19] = 0.25*bl[17]*dvpar_sq*wvpar+0.1936491673103708*bl[10]*dvpar_sq*wvpar+0.1118033988749895*bl[6]*dvpar_sq*wvpar; 

  double fUpwindQuad_u[27] = {0.0};
  double fUpwindQuad_q[27] = {0.0};
  double fUpwind_u[20] = {0.0};;
  double fUpwind_q[20] = {0.0};
  double Ghat_u[20] = {0.0}; 
  double Ghat_q[20] = {0.0}; 

  if ((-0.4242640687119281*u[12])-0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])+0.6363961030678926*u[4]-0.4743416490252568*(u[2]+u[1])+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_u[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.5692099788303082*(q[19]+q[18]+q[17])-0.4242640687119281*(q[16]+q[15])-0.4242640687119285*(q[14]+q[13])-0.4242640687119281*q[12]-0.4242640687119285*q[11]-0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])+0.6363961030678926*(q[6]+q[5]+q[4])-0.4743416490252568*(q[3]+q[2]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_q[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if ((-0.4242640687119281*u[12])-0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])+0.6363961030678926*u[4]-0.4743416490252568*(u[2]+u[1])+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_u[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if ((-0.711512473537885*q[19])+0.5303300858899104*(q[16]+q[15])-0.4242640687119281*q[12]-0.4242640687119285*q[11]-0.3952847075210473*q[9]+0.3162277660168379*(q[8]+q[7])+0.6363961030678926*q[4]-0.4743416490252568*(q[2]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_q[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if ((-0.4242640687119281*u[12])-0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])+0.6363961030678926*u[4]-0.4743416490252568*(u[2]+u[1])+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_u[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.5692099788303082*q[19]-0.5692099788303082*(q[18]+q[17])-0.4242640687119281*(q[16]+q[15])+0.4242640687119285*(q[14]+q[13])-0.4242640687119281*q[12]-0.4242640687119285*q[11]+0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])-0.6363961030678926*(q[6]+q[5])+0.6363961030678926*q[4]+0.4743416490252568*q[3]-0.4743416490252568*(q[2]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_q[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.5303300858899104*u[12]-0.3952847075210473*u[8]+0.3162277660168379*u[7]-0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_u[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  if ((-0.711512473537885*q[18])-0.4242640687119281*q[15]+0.5303300858899104*q[14]-0.4242640687119285*q[13]+0.5303300858899104*q[12]+0.3162277660168379*q[9]-0.3952847075210473*q[8]+0.3162277660168379*q[7]+0.6363961030678926*q[5]-0.4743416490252568*(q[3]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_q[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.5303300858899104*u[12]-0.3952847075210473*u[8]+0.3162277660168379*u[7]-0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_u[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.5303300858899104*(q[15]+q[12])-0.3952847075210473*(q[9]+q[8])+0.3162277660168379*q[7]-0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_q[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.5303300858899104*u[12]-0.3952847075210473*u[8]+0.3162277660168379*u[7]-0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_u[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  if (0.711512473537885*q[18]-0.4242640687119281*q[15]-0.5303300858899104*q[14]+0.4242640687119285*q[13]+0.5303300858899104*q[12]+0.3162277660168379*q[9]-0.3952847075210473*q[8]+0.3162277660168379*q[7]-0.6363961030678926*q[5]+0.4743416490252568*q[3]-0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_q[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  if ((-0.4242640687119281*u[12])+0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])-0.6363961030678926*u[4]+0.4743416490252568*u[2]-0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_u[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  if ((-0.5692099788303082*q[19])+0.5692099788303082*q[18]-0.5692099788303082*q[17]+0.4242640687119281*q[16]-0.4242640687119281*q[15]-0.4242640687119285*(q[14]+q[13])-0.4242640687119281*q[12]+0.4242640687119285*q[11]+0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])-0.6363961030678926*q[6]+0.6363961030678926*q[5]-0.6363961030678926*q[4]-0.4743416490252568*q[3]+0.4743416490252568*q[2]-0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_q[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  if ((-0.4242640687119281*u[12])+0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])-0.6363961030678926*u[4]+0.4743416490252568*u[2]-0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_u[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  if (0.711512473537885*q[19]-0.5303300858899104*q[16]+0.5303300858899104*q[15]-0.4242640687119281*q[12]+0.4242640687119285*q[11]-0.3952847075210473*q[9]+0.3162277660168379*(q[8]+q[7])-0.6363961030678926*q[4]+0.4743416490252568*q[2]-0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_q[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  if ((-0.4242640687119281*u[12])+0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])-0.6363961030678926*u[4]+0.4743416490252568*u[2]-0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_u[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 
  if ((-0.5692099788303082*(q[19]+q[18]))+0.5692099788303082*q[17]+0.4242640687119281*q[16]-0.4242640687119281*q[15]+0.4242640687119285*(q[14]+q[13])-0.4242640687119281*q[12]+0.4242640687119285*q[11]-0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])+0.6363961030678926*q[6]-0.6363961030678926*(q[5]+q[4])+0.4743416490252568*(q[3]+q[2])-0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_q[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 
  if (0.5303300858899104*u[11]+0.3162277660168379*u[8]-0.3952847075210473*u[7]-0.4743416490252568*u[2]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_u[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fr); 
  } 
  if ((-0.711512473537885*q[17])-0.4242640687119281*q[16]-0.4242640687119285*q[14]+0.5303300858899104*(q[13]+q[11])+0.3162277660168379*(q[9]+q[8])-0.3952847075210473*q[7]+0.6363961030678926*q[6]-0.4743416490252568*(q[3]+q[2])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_q[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*u[11]+0.3162277660168379*u[8]-0.3952847075210473*u[7]-0.4743416490252568*u[2]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_u[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fr); 
  } 
  if (0.5303300858899104*(q[16]+q[11])-0.3952847075210473*q[9]+0.3162277660168379*q[8]-0.3952847075210473*q[7]-0.4743416490252568*q[2]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_q[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fr); 
  } 
  if (0.5303300858899104*u[11]+0.3162277660168379*u[8]-0.3952847075210473*u[7]-0.4743416490252568*u[2]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_u[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fr); 
  } 
  if (0.711512473537885*q[17]-0.4242640687119281*q[16]+0.4242640687119285*q[14]-0.5303300858899104*q[13]+0.5303300858899104*q[11]+0.3162277660168379*(q[9]+q[8])-0.3952847075210473*q[7]-0.6363961030678926*q[6]+0.4743416490252568*q[3]-0.4743416490252568*q[2]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_q[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fr); 
  } 
  if (0.3535533905932737*u[0]-0.3952847075210473*(u[8]+u[7]) > 0) { 
    fUpwindQuad_u[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_u[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fr); 
  } 
  if (0.5303300858899104*(q[14]+q[13])+0.3162277660168379*q[9]-0.3952847075210473*(q[8]+q[7])-0.4743416490252568*q[3]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_q[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*u[0]-0.3952847075210473*(u[8]+u[7]) > 0) { 
    fUpwindQuad_u[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_u[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fr); 
  } 
  if (0.3535533905932737*q[0]-0.3952847075210473*(q[9]+q[8]+q[7]) > 0) { 
    fUpwindQuad_q[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_q[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fr); 
  } 
  if (0.3535533905932737*u[0]-0.3952847075210473*(u[8]+u[7]) > 0) { 
    fUpwindQuad_u[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_u[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fr); 
  } 
  if ((-0.5303300858899104*(q[14]+q[13]))+0.3162277660168379*q[9]-0.3952847075210473*(q[8]+q[7])+0.4743416490252568*q[3]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_q[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fr); 
  } 
  if ((-0.5303300858899104*u[11])+0.3162277660168379*u[8]-0.3952847075210473*u[7]+0.4743416490252568*u[2]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_u[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fr); 
  } 
  if (0.711512473537885*q[17]+0.4242640687119281*q[16]-0.4242640687119285*q[14]+0.5303300858899104*q[13]-0.5303300858899104*q[11]+0.3162277660168379*(q[9]+q[8])-0.3952847075210473*q[7]-0.6363961030678926*q[6]-0.4743416490252568*q[3]+0.4743416490252568*q[2]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_q[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fr); 
  } 
  if ((-0.5303300858899104*u[11])+0.3162277660168379*u[8]-0.3952847075210473*u[7]+0.4743416490252568*u[2]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_u[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fr); 
  } 
  if ((-0.5303300858899104*(q[16]+q[11]))-0.3952847075210473*q[9]+0.3162277660168379*q[8]-0.3952847075210473*q[7]+0.4743416490252568*q[2]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_q[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fr); 
  } 
  if ((-0.5303300858899104*u[11])+0.3162277660168379*u[8]-0.3952847075210473*u[7]+0.4743416490252568*u[2]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_u[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fr); 
  } 
  if ((-0.711512473537885*q[17])+0.4242640687119281*q[16]+0.4242640687119285*q[14]-0.5303300858899104*(q[13]+q[11])+0.3162277660168379*(q[9]+q[8])-0.3952847075210473*q[7]+0.6363961030678926*q[6]+0.4743416490252568*(q[3]+q[2])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_q[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fr); 
  } 
  if (0.4242640687119281*u[12]-0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])-0.6363961030678926*u[4]-0.4743416490252568*u[2]+0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fl); 
  } else { 
    fUpwindQuad_u[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fr); 
  } 
  if ((-0.5692099788303082*(q[19]+q[18]))+0.5692099788303082*q[17]-0.4242640687119281*q[16]+0.4242640687119281*q[15]-0.4242640687119285*(q[14]+q[13])+0.4242640687119281*q[12]-0.4242640687119285*q[11]+0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])+0.6363961030678926*q[6]-0.6363961030678926*(q[5]+q[4])-0.4743416490252568*(q[3]+q[2])+0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fl); 
  } else { 
    fUpwindQuad_q[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fr); 
  } 
  if (0.4242640687119281*u[12]-0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])-0.6363961030678926*u[4]-0.4743416490252568*u[2]+0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_u[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fr); 
  } 
  if (0.711512473537885*q[19]+0.5303300858899104*q[16]-0.5303300858899104*q[15]+0.4242640687119281*q[12]-0.4242640687119285*q[11]-0.3952847075210473*q[9]+0.3162277660168379*(q[8]+q[7])-0.6363961030678926*q[4]-0.4743416490252568*q[2]+0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_q[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fr); 
  } 
  if (0.4242640687119281*u[12]-0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])-0.6363961030678926*u[4]-0.4743416490252568*u[2]+0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fl); 
  } else { 
    fUpwindQuad_u[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fr); 
  } 
  if ((-0.5692099788303082*q[19])+0.5692099788303082*q[18]-0.5692099788303082*q[17]-0.4242640687119281*q[16]+0.4242640687119281*q[15]+0.4242640687119285*(q[14]+q[13])+0.4242640687119281*q[12]-0.4242640687119285*q[11]-0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])-0.6363961030678926*q[6]+0.6363961030678926*q[5]-0.6363961030678926*q[4]+0.4743416490252568*q[3]-0.4743416490252568*q[2]+0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fl); 
  } else { 
    fUpwindQuad_q[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fr); 
  } 
  if ((-0.5303300858899104*u[12])-0.3952847075210473*u[8]+0.3162277660168379*u[7]+0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fl); 
  } else { 
    fUpwindQuad_u[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fr); 
  } 
  if (0.711512473537885*q[18]+0.4242640687119281*q[15]+0.5303300858899104*q[14]-0.4242640687119285*q[13]-0.5303300858899104*q[12]+0.3162277660168379*q[9]-0.3952847075210473*q[8]+0.3162277660168379*q[7]-0.6363961030678926*q[5]-0.4743416490252568*q[3]+0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fl); 
  } else { 
    fUpwindQuad_q[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fr); 
  } 
  if ((-0.5303300858899104*u[12])-0.3952847075210473*u[8]+0.3162277660168379*u[7]+0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fl); 
  } else { 
    fUpwindQuad_u[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fr); 
  } 
  if ((-0.5303300858899104*(q[15]+q[12]))-0.3952847075210473*(q[9]+q[8])+0.3162277660168379*q[7]+0.4743416490252568*q[1]+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fl); 
  } else { 
    fUpwindQuad_q[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fr); 
  } 
  if ((-0.5303300858899104*u[12])-0.3952847075210473*u[8]+0.3162277660168379*u[7]+0.4743416490252568*u[1]+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_u[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fr); 
  } 
  if ((-0.711512473537885*q[18])+0.4242640687119281*q[15]-0.5303300858899104*q[14]+0.4242640687119285*q[13]-0.5303300858899104*q[12]+0.3162277660168379*q[9]-0.3952847075210473*q[8]+0.3162277660168379*q[7]+0.6363961030678926*q[5]+0.4743416490252568*(q[3]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_q[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fr); 
  } 
  if (0.4242640687119281*u[12]+0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])+0.6363961030678926*u[4]+0.4743416490252568*(u[2]+u[1])+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fl); 
  } else { 
    fUpwindQuad_u[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fr); 
  } 
  if (0.5692099788303082*q[19]-0.5692099788303082*(q[18]+q[17])+0.4242640687119281*(q[16]+q[15])-0.4242640687119285*(q[14]+q[13])+0.4242640687119281*q[12]+0.4242640687119285*q[11]-0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])-0.6363961030678926*(q[6]+q[5])+0.6363961030678926*q[4]-0.4743416490252568*q[3]+0.4743416490252568*(q[2]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fl); 
  } else { 
    fUpwindQuad_q[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fr); 
  } 
  if (0.4242640687119281*u[12]+0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])+0.6363961030678926*u[4]+0.4743416490252568*(u[2]+u[1])+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fl); 
  } else { 
    fUpwindQuad_u[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fr); 
  } 
  if ((-0.711512473537885*q[19])-0.5303300858899104*(q[16]+q[15])+0.4242640687119281*q[12]+0.4242640687119285*q[11]-0.3952847075210473*q[9]+0.3162277660168379*(q[8]+q[7])+0.6363961030678926*q[4]+0.4743416490252568*(q[2]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fl); 
  } else { 
    fUpwindQuad_q[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fr); 
  } 
  if (0.4242640687119281*u[12]+0.4242640687119285*u[11]+0.3162277660168379*(u[8]+u[7])+0.6363961030678926*u[4]+0.4743416490252568*(u[2]+u[1])+0.3535533905932737*u[0] > 0) { 
    fUpwindQuad_u[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_u[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fr); 
  } 
  if (0.5692099788303082*(q[19]+q[18]+q[17])+0.4242640687119281*(q[16]+q[15])+0.4242640687119285*(q[14]+q[13])+0.4242640687119281*q[12]+0.4242640687119285*q[11]+0.853814968245462*q[10]+0.3162277660168379*(q[9]+q[8]+q[7])+0.6363961030678926*(q[6]+q[5]+q[4])+0.4743416490252568*(q[3]+q[2]+q[1])+0.3535533905932737*q[0] > 0) { 
    fUpwindQuad_q[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_q[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_u, fUpwind_u); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_q, fUpwind_q); 

  rho_flux[0] += (0.5*fUpwind_u[12]*u[12]+0.5*fUpwind_u[11]*u[11]+0.5*fUpwind_u[8]*u[8]+0.5*fUpwind_u[7]*u[7]+0.5*fUpwind_u[4]*u[4]+0.5*fUpwind_u[2]*u[2]+0.5*fUpwind_u[1]*u[1]+0.5*fUpwind_u[0]*u[0])*mass*volFact; 
  rho_flux[1] += (0.5000000000000001*fUpwind_u[8]*u[12]+0.5000000000000001*u[8]*fUpwind_u[12]+0.447213595499958*fUpwind_u[4]*u[11]+0.447213595499958*u[4]*fUpwind_u[11]+0.4472135954999579*fUpwind_u[1]*u[7]+0.4472135954999579*u[1]*fUpwind_u[7]+0.5*fUpwind_u[2]*u[4]+0.5*u[2]*fUpwind_u[4]+0.5*fUpwind_u[0]*u[1]+0.5*u[0]*fUpwind_u[1])*mass*volFact; 
  rho_flux[2] += (0.447213595499958*fUpwind_u[4]*u[12]+0.447213595499958*u[4]*fUpwind_u[12]+0.5000000000000001*fUpwind_u[7]*u[11]+0.5000000000000001*u[7]*fUpwind_u[11]+0.4472135954999579*fUpwind_u[2]*u[8]+0.4472135954999579*u[2]*fUpwind_u[8]+0.5*fUpwind_u[1]*u[4]+0.5*u[1]*fUpwind_u[4]+0.5*fUpwind_u[0]*u[2]+0.5*u[0]*fUpwind_u[2])*mass*volFact; 
  rho_flux[3] += (0.4*fUpwind_u[11]*u[12]+0.447213595499958*fUpwind_u[2]*u[12]+0.4*u[11]*fUpwind_u[12]+0.447213595499958*u[2]*fUpwind_u[12]+0.447213595499958*fUpwind_u[1]*u[11]+0.447213595499958*u[1]*fUpwind_u[11]+0.4472135954999579*fUpwind_u[4]*u[8]+0.4472135954999579*u[4]*fUpwind_u[8]+0.4472135954999579*fUpwind_u[4]*u[7]+0.4472135954999579*u[4]*fUpwind_u[7]+0.5*fUpwind_u[0]*u[4]+0.5*u[0]*fUpwind_u[4]+0.5*fUpwind_u[1]*u[2]+0.5*u[1]*fUpwind_u[2])*mass*volFact; 
  rho_flux[4] += (0.4472135954999579*fUpwind_u[12]*u[12]+0.31943828249997*fUpwind_u[11]*u[11]+0.5000000000000001*fUpwind_u[2]*u[11]+0.5000000000000001*u[2]*fUpwind_u[11]+0.31943828249997*fUpwind_u[7]*u[7]+0.5*fUpwind_u[0]*u[7]+0.5*u[0]*fUpwind_u[7]+0.4472135954999579*fUpwind_u[4]*u[4]+0.4472135954999579*fUpwind_u[1]*u[1])*mass*volFact; 
  rho_flux[5] += (0.31943828249997*fUpwind_u[12]*u[12]+0.5000000000000001*fUpwind_u[1]*u[12]+0.5000000000000001*u[1]*fUpwind_u[12]+0.4472135954999579*fUpwind_u[11]*u[11]+0.31943828249997*fUpwind_u[8]*u[8]+0.5*fUpwind_u[0]*u[8]+0.5*u[0]*fUpwind_u[8]+0.4472135954999579*fUpwind_u[4]*u[4]+0.4472135954999579*fUpwind_u[2]*u[2])*mass*volFact; 
  rho_flux[6] += (0.4*fUpwind_u[4]*u[12]+0.4*u[4]*fUpwind_u[12]+0.4472135954999579*fUpwind_u[8]*u[11]+0.31943828249997*fUpwind_u[7]*u[11]+0.5*fUpwind_u[0]*u[11]+0.4472135954999579*u[8]*fUpwind_u[11]+0.31943828249997*u[7]*fUpwind_u[11]+0.5*u[0]*fUpwind_u[11]+0.5000000000000001*fUpwind_u[2]*u[7]+0.5000000000000001*u[2]*fUpwind_u[7]+0.447213595499958*fUpwind_u[1]*u[4]+0.447213595499958*u[1]*fUpwind_u[4])*mass*volFact; 
  rho_flux[7] += (0.31943828249997*fUpwind_u[8]*u[12]+0.4472135954999579*fUpwind_u[7]*u[12]+0.5*fUpwind_u[0]*u[12]+0.31943828249997*u[8]*fUpwind_u[12]+0.4472135954999579*u[7]*fUpwind_u[12]+0.5*u[0]*fUpwind_u[12]+0.4*fUpwind_u[4]*u[11]+0.4*u[4]*fUpwind_u[11]+0.5000000000000001*fUpwind_u[1]*u[8]+0.5000000000000001*u[1]*fUpwind_u[8]+0.447213595499958*fUpwind_u[2]*u[4]+0.447213595499958*u[2]*fUpwind_u[4])*mass*volFact; 
  heat_flux[0] += (0.5*fUpwind_q[19]*q[19]+0.5*fUpwind_q[18]*q[18]+0.5*fUpwind_q[17]*q[17]+0.5*fUpwind_q[16]*q[16]+0.5*fUpwind_q[15]*q[15]+0.5*fUpwind_q[14]*q[14]+0.5*fUpwind_q[13]*q[13]+0.5*fUpwind_q[12]*q[12]+0.5*fUpwind_q[11]*q[11]+0.5*fUpwind_q[10]*q[10]+0.5*fUpwind_q[9]*q[9]+0.5*fUpwind_q[8]*q[8]+0.5*fUpwind_q[7]*q[7]+0.5*fUpwind_q[6]*q[6]+0.5*fUpwind_q[5]*q[5]+0.5*fUpwind_q[4]*q[4]+0.5*fUpwind_q[3]*q[3]+0.5*fUpwind_q[2]*q[2]+0.5*fUpwind_q[1]*q[1]+0.5*fUpwind_q[0]*q[0])*mass*volFact; 
  heat_flux[1] += (0.5000000000000001*fUpwind_q[16]*q[19]+0.5000000000000001*q[16]*fUpwind_q[19]+0.5000000000000001*fUpwind_q[14]*q[18]+0.5000000000000001*q[14]*fUpwind_q[18]+0.4472135954999579*fUpwind_q[10]*q[17]+0.4472135954999579*q[10]*fUpwind_q[17]+0.5000000000000001*fUpwind_q[9]*q[15]+0.5000000000000001*q[9]*fUpwind_q[15]+0.447213595499958*fUpwind_q[5]*q[13]+0.447213595499958*q[5]*fUpwind_q[13]+0.5000000000000001*fUpwind_q[8]*q[12]+0.5000000000000001*q[8]*fUpwind_q[12]+0.447213595499958*fUpwind_q[4]*q[11]+0.447213595499958*q[4]*fUpwind_q[11]+0.5*fUpwind_q[6]*q[10]+0.5*q[6]*fUpwind_q[10]+0.4472135954999579*fUpwind_q[1]*q[7]+0.4472135954999579*q[1]*fUpwind_q[7]+0.5*fUpwind_q[3]*q[5]+0.5*q[3]*fUpwind_q[5]+0.5*fUpwind_q[2]*q[4]+0.5*q[2]*fUpwind_q[4]+0.5*fUpwind_q[0]*q[1]+0.5*q[0]*fUpwind_q[1])*mass*volFact; 
  heat_flux[2] += (0.5000000000000001*fUpwind_q[15]*q[19]+0.5000000000000001*q[15]*fUpwind_q[19]+0.4472135954999579*fUpwind_q[10]*q[18]+0.4472135954999579*q[10]*fUpwind_q[18]+0.5000000000000001*fUpwind_q[13]*q[17]+0.5000000000000001*q[13]*fUpwind_q[17]+0.5000000000000001*fUpwind_q[9]*q[16]+0.5000000000000001*q[9]*fUpwind_q[16]+0.447213595499958*fUpwind_q[6]*q[14]+0.447213595499958*q[6]*fUpwind_q[14]+0.447213595499958*fUpwind_q[4]*q[12]+0.447213595499958*q[4]*fUpwind_q[12]+0.5000000000000001*fUpwind_q[7]*q[11]+0.5000000000000001*q[7]*fUpwind_q[11]+0.5*fUpwind_q[5]*q[10]+0.5*q[5]*fUpwind_q[10]+0.4472135954999579*fUpwind_q[2]*q[8]+0.4472135954999579*q[2]*fUpwind_q[8]+0.5*fUpwind_q[3]*q[6]+0.5*q[3]*fUpwind_q[6]+0.5*fUpwind_q[1]*q[4]+0.5*q[1]*fUpwind_q[4]+0.5*fUpwind_q[0]*q[2]+0.5*q[0]*fUpwind_q[2])*mass*volFact; 
  heat_flux[3] += (0.5*fUpwind_q[9]*q[19]+0.5*q[9]*fUpwind_q[19]+0.4*fUpwind_q[17]*q[18]+0.4472135954999579*fUpwind_q[6]*q[18]+0.4*q[17]*fUpwind_q[18]+0.4472135954999579*q[6]*fUpwind_q[18]+0.4472135954999579*fUpwind_q[5]*q[17]+0.4472135954999579*q[5]*fUpwind_q[17]+0.5*fUpwind_q[15]*q[16]+0.5*q[15]*fUpwind_q[16]+0.447213595499958*fUpwind_q[10]*q[14]+0.447213595499958*q[10]*fUpwind_q[14]+0.447213595499958*fUpwind_q[10]*q[13]+0.447213595499958*q[10]*fUpwind_q[13]+0.4*fUpwind_q[11]*q[12]+0.447213595499958*fUpwind_q[2]*q[12]+0.4*q[11]*fUpwind_q[12]+0.447213595499958*q[2]*fUpwind_q[12]+0.447213595499958*fUpwind_q[1]*q[11]+0.447213595499958*q[1]*fUpwind_q[11]+0.5*fUpwind_q[3]*q[10]+0.5*q[3]*fUpwind_q[10]+0.4472135954999579*fUpwind_q[4]*q[8]+0.4472135954999579*q[4]*fUpwind_q[8]+0.4472135954999579*fUpwind_q[4]*q[7]+0.4472135954999579*q[4]*fUpwind_q[7]+0.5*fUpwind_q[5]*q[6]+0.5*q[5]*fUpwind_q[6]+0.5*fUpwind_q[0]*q[4]+0.5*q[0]*fUpwind_q[4]+0.5*fUpwind_q[1]*q[2]+0.5*q[1]*fUpwind_q[2])*mass*volFact; 
  heat_flux[4] += (0.4472135954999579*fUpwind_q[19]*q[19]+0.4472135954999579*fUpwind_q[18]*q[18]+0.31943828249997*fUpwind_q[17]*q[17]+0.5*fUpwind_q[6]*q[17]+0.5*q[6]*fUpwind_q[17]+0.4472135954999579*fUpwind_q[15]*q[15]+0.31943828249997*fUpwind_q[13]*q[13]+0.5000000000000001*fUpwind_q[3]*q[13]+0.5000000000000001*q[3]*fUpwind_q[13]+0.4472135954999579*fUpwind_q[12]*q[12]+0.31943828249997*fUpwind_q[11]*q[11]+0.5000000000000001*fUpwind_q[2]*q[11]+0.5000000000000001*q[2]*fUpwind_q[11]+0.4472135954999579*fUpwind_q[10]*q[10]+0.31943828249997*fUpwind_q[7]*q[7]+0.5*fUpwind_q[0]*q[7]+0.5*q[0]*fUpwind_q[7]+0.4472135954999579*fUpwind_q[5]*q[5]+0.4472135954999579*fUpwind_q[4]*q[4]+0.4472135954999579*fUpwind_q[1]*q[1])*mass*volFact; 
  heat_flux[5] += (0.4472135954999579*fUpwind_q[19]*q[19]+0.31943828249997*fUpwind_q[18]*q[18]+0.5*fUpwind_q[5]*q[18]+0.5*q[5]*fUpwind_q[18]+0.4472135954999579*fUpwind_q[17]*q[17]+0.4472135954999579*fUpwind_q[16]*q[16]+0.31943828249997*fUpwind_q[14]*q[14]+0.5000000000000001*fUpwind_q[3]*q[14]+0.5000000000000001*q[3]*fUpwind_q[14]+0.31943828249997*fUpwind_q[12]*q[12]+0.5000000000000001*fUpwind_q[1]*q[12]+0.5000000000000001*q[1]*fUpwind_q[12]+0.4472135954999579*fUpwind_q[11]*q[11]+0.4472135954999579*fUpwind_q[10]*q[10]+0.31943828249997*fUpwind_q[8]*q[8]+0.5*fUpwind_q[0]*q[8]+0.5*q[0]*fUpwind_q[8]+0.4472135954999579*fUpwind_q[6]*q[6]+0.4472135954999579*fUpwind_q[4]*q[4]+0.4472135954999579*fUpwind_q[2]*q[2])*mass*volFact; 
  heat_flux[6] += (0.4472135954999579*fUpwind_q[15]*q[19]+0.4472135954999579*q[15]*fUpwind_q[19]+0.4*fUpwind_q[10]*q[18]+0.4*q[10]*fUpwind_q[18]+0.4472135954999579*fUpwind_q[14]*q[17]+0.31943828249997*fUpwind_q[13]*q[17]+0.5000000000000001*fUpwind_q[3]*q[17]+0.4472135954999579*q[14]*fUpwind_q[17]+0.31943828249997*q[13]*fUpwind_q[17]+0.5000000000000001*q[3]*fUpwind_q[17]+0.5*fUpwind_q[6]*q[13]+0.5*q[6]*fUpwind_q[13]+0.4*fUpwind_q[4]*q[12]+0.4*q[4]*fUpwind_q[12]+0.4472135954999579*fUpwind_q[8]*q[11]+0.31943828249997*fUpwind_q[7]*q[11]+0.5*fUpwind_q[0]*q[11]+0.4472135954999579*q[8]*fUpwind_q[11]+0.31943828249997*q[7]*fUpwind_q[11]+0.5*q[0]*fUpwind_q[11]+0.447213595499958*fUpwind_q[5]*q[10]+0.447213595499958*q[5]*fUpwind_q[10]+0.5000000000000001*fUpwind_q[2]*q[7]+0.5000000000000001*q[2]*fUpwind_q[7]+0.447213595499958*fUpwind_q[1]*q[4]+0.447213595499958*q[1]*fUpwind_q[4])*mass*volFact; 
  heat_flux[7] += (0.4472135954999579*fUpwind_q[16]*q[19]+0.4472135954999579*q[16]*fUpwind_q[19]+0.31943828249997*fUpwind_q[14]*q[18]+0.4472135954999579*fUpwind_q[13]*q[18]+0.5000000000000001*fUpwind_q[3]*q[18]+0.31943828249997*q[14]*fUpwind_q[18]+0.4472135954999579*q[13]*fUpwind_q[18]+0.5000000000000001*q[3]*fUpwind_q[18]+0.4*fUpwind_q[10]*q[17]+0.4*q[10]*fUpwind_q[17]+0.5*fUpwind_q[5]*q[14]+0.5*q[5]*fUpwind_q[14]+0.31943828249997*fUpwind_q[8]*q[12]+0.4472135954999579*fUpwind_q[7]*q[12]+0.5*fUpwind_q[0]*q[12]+0.31943828249997*q[8]*fUpwind_q[12]+0.4472135954999579*q[7]*fUpwind_q[12]+0.5*q[0]*fUpwind_q[12]+0.4*fUpwind_q[4]*q[11]+0.4*q[4]*fUpwind_q[11]+0.447213595499958*fUpwind_q[6]*q[10]+0.447213595499958*q[6]*fUpwind_q[10]+0.5000000000000001*fUpwind_q[1]*q[8]+0.5000000000000001*q[1]*fUpwind_q[8]+0.447213595499958*fUpwind_q[2]*q[4]+0.447213595499958*q[2]*fUpwind_q[4])*mass*volFact; 
} 
