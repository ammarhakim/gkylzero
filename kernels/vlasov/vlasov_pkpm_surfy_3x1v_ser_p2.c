#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfy_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *u_i, const double *bvar, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i: Input bulk velocity (ux,uy,uz) in cell being updated (ASSUMED TO BE CONTINUOUS).
  // bvar: Input magnetic field unit vector in cell being updated (ASSUMED TO BE CONTINUOUS).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[1]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *uc = &u_i[20]; 
  const double *bc = &bvar[20]; 
  double alphaSurf_l[20] = {0.0}; 
  alphaSurf_l[0] = 2.23606797749979*bc[8]*wvpar-1.732050807568877*bc[2]*wvpar+bc[0]*wvpar+2.23606797749979*uc[8]-1.732050807568877*uc[2]+uc[0]; 
  alphaSurf_l[1] = 2.23606797749979*bc[12]*wvpar-1.732050807568877*bc[4]*wvpar+bc[1]*wvpar+2.23606797749979*uc[12]-1.732050807568877*uc[4]+uc[1]; 
  alphaSurf_l[2] = 2.23606797749979*bc[14]*wvpar-1.732050807568877*bc[6]*wvpar+bc[3]*wvpar+2.23606797749979*uc[14]-1.732050807568877*uc[6]+uc[3]; 
  alphaSurf_l[3] = 0.6454972243679029*bc[8]*dvpar-0.5*bc[2]*dvpar+0.2886751345948129*bc[0]*dvpar; 
  alphaSurf_l[4] = 2.23606797749979*bc[18]*wvpar-1.732050807568877*bc[10]*wvpar+bc[5]*wvpar+2.23606797749979*uc[18]-1.732050807568877*uc[10]+uc[5]; 
  alphaSurf_l[5] = 0.6454972243679028*bc[12]*dvpar-0.5*bc[4]*dvpar+0.2886751345948129*bc[1]*dvpar; 
  alphaSurf_l[6] = 0.6454972243679028*bc[14]*dvpar-0.5*bc[6]*dvpar+0.2886751345948129*bc[3]*dvpar; 
  alphaSurf_l[7] = (-1.732050807568877*bc[11]*wvpar)+bc[7]*wvpar-1.732050807568877*uc[11]+uc[7]; 
  alphaSurf_l[8] = (-1.732050807568877*bc[16]*wvpar)+bc[9]*wvpar-1.732050807568877*uc[16]+uc[9]; 
  alphaSurf_l[10] = 0.6454972243679029*bc[18]*dvpar-0.5*bc[10]*dvpar+0.2886751345948129*bc[5]*dvpar; 
  alphaSurf_l[11] = (-1.732050807568877*bc[17]*wvpar)+bc[13]*wvpar-1.732050807568877*uc[17]+uc[13]; 
  alphaSurf_l[12] = (-1.732050807568877*bc[19]*wvpar)+bc[15]*wvpar-1.732050807568877*uc[19]+uc[15]; 
  alphaSurf_l[13] = 0.2886751345948129*bc[7]*dvpar-0.5*bc[11]*dvpar; 
  alphaSurf_l[14] = 0.2886751345948129*bc[9]*dvpar-0.5*bc[16]*dvpar; 
  alphaSurf_l[17] = 0.2886751345948129*bc[13]*dvpar-0.5*bc[17]*dvpar; 
  alphaSurf_l[18] = 0.2886751345948129*bc[15]*dvpar-0.5*bc[19]*dvpar; 

  double alphaSurf_r[20] = {0.0}; 
  alphaSurf_r[0] = 2.23606797749979*bc[8]*wvpar+1.732050807568877*bc[2]*wvpar+bc[0]*wvpar+2.23606797749979*uc[8]+1.732050807568877*uc[2]+uc[0]; 
  alphaSurf_r[1] = 2.23606797749979*bc[12]*wvpar+1.732050807568877*bc[4]*wvpar+bc[1]*wvpar+2.23606797749979*uc[12]+1.732050807568877*uc[4]+uc[1]; 
  alphaSurf_r[2] = 2.23606797749979*bc[14]*wvpar+1.732050807568877*bc[6]*wvpar+bc[3]*wvpar+2.23606797749979*uc[14]+1.732050807568877*uc[6]+uc[3]; 
  alphaSurf_r[3] = 0.6454972243679029*bc[8]*dvpar+0.5*bc[2]*dvpar+0.2886751345948129*bc[0]*dvpar; 
  alphaSurf_r[4] = 2.23606797749979*bc[18]*wvpar+1.732050807568877*bc[10]*wvpar+bc[5]*wvpar+2.23606797749979*uc[18]+1.732050807568877*uc[10]+uc[5]; 
  alphaSurf_r[5] = 0.6454972243679028*bc[12]*dvpar+0.5*bc[4]*dvpar+0.2886751345948129*bc[1]*dvpar; 
  alphaSurf_r[6] = 0.6454972243679028*bc[14]*dvpar+0.5*bc[6]*dvpar+0.2886751345948129*bc[3]*dvpar; 
  alphaSurf_r[7] = 1.732050807568877*bc[11]*wvpar+bc[7]*wvpar+1.732050807568877*uc[11]+uc[7]; 
  alphaSurf_r[8] = 1.732050807568877*bc[16]*wvpar+bc[9]*wvpar+1.732050807568877*uc[16]+uc[9]; 
  alphaSurf_r[10] = 0.6454972243679029*bc[18]*dvpar+0.5*bc[10]*dvpar+0.2886751345948129*bc[5]*dvpar; 
  alphaSurf_r[11] = 1.732050807568877*bc[17]*wvpar+bc[13]*wvpar+1.732050807568877*uc[17]+uc[13]; 
  alphaSurf_r[12] = 1.732050807568877*bc[19]*wvpar+bc[15]*wvpar+1.732050807568877*uc[19]+uc[15]; 
  alphaSurf_r[13] = 0.5*bc[11]*dvpar+0.2886751345948129*bc[7]*dvpar; 
  alphaSurf_r[14] = 0.5*bc[16]*dvpar+0.2886751345948129*bc[9]*dvpar; 
  alphaSurf_r[17] = 0.5*bc[17]*dvpar+0.2886751345948129*bc[13]*dvpar; 
  alphaSurf_r[18] = 0.5*bc[19]*dvpar+0.2886751345948129*bc[15]*dvpar; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17])-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5]+alphaSurf_l[4])-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17])-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5]+alphaSurf_r[4])-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  if ((-0.4242640687119281*alphaSurf_l[12])-0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fc); 
  } 
  if ((-0.4242640687119281*alphaSurf_r[12])-0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  if ((-0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17]))+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5])+0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fc); 
  } 
  if ((-0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17]))+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5])+0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[18])+0.5303300858899104*alphaSurf_l[14]-0.4242640687119285*alphaSurf_l[13]+0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[5]-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[18])+0.5303300858899104*alphaSurf_r[14]-0.4242640687119285*alphaSurf_r[13]+0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[5]-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fr); 
  } 
  if (0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fc); 
  } 
  if (0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[18]-0.5303300858899104*alphaSurf_l[14]+0.4242640687119285*alphaSurf_l[13]+0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[5]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[18]-0.5303300858899104*alphaSurf_r[14]+0.4242640687119285*alphaSurf_r[13]+0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[5]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fr); 
  } 
  if (0.5692099788303082*alphaSurf_l[18]-0.5692099788303082*alphaSurf_l[17]-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[6]+0.6363961030678926*alphaSurf_l[5]-0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*alphaSurf_l[2]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fc); 
  } 
  if (0.5692099788303082*alphaSurf_r[18]-0.5692099788303082*alphaSurf_r[17]-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[6]+0.6363961030678926*alphaSurf_r[5]-0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*alphaSurf_r[2]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fr); 
  } 
  if ((-0.4242640687119281*alphaSurf_l[12])+0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[2]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fc); 
  } 
  if ((-0.4242640687119281*alphaSurf_r[12])+0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[2]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fr); 
  } 
  if ((-0.5692099788303082*alphaSurf_l[18])+0.5692099788303082*alphaSurf_l[17]+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[6]-0.6363961030678926*(alphaSurf_l[5]+alphaSurf_l[4])+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fc); 
  } 
  if ((-0.5692099788303082*alphaSurf_r[18])+0.5692099788303082*alphaSurf_r[17]+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[6]-0.6363961030678926*(alphaSurf_r[5]+alphaSurf_r[4])+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[17])-0.4242640687119285*alphaSurf_l[14]+0.5303300858899104*(alphaSurf_l[13]+alphaSurf_l[11])+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[6]-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[17])-0.4242640687119285*alphaSurf_r[14]+0.5303300858899104*(alphaSurf_r[13]+alphaSurf_r[11])+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[6]-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fc); 
  } 
  if (0.5303300858899104*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[17]+0.4242640687119285*alphaSurf_l[14]-0.5303300858899104*alphaSurf_l[13]+0.5303300858899104*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[6]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[17]+0.4242640687119285*alphaSurf_r[14]-0.5303300858899104*alphaSurf_r[13]+0.5303300858899104*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[6]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*(alphaSurf_l[14]+alphaSurf_l[13])-0.3952847075210473*(alphaSurf_l[8]+alphaSurf_l[7])-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fc); 
  } 
  if (0.5303300858899104*(alphaSurf_r[14]+alphaSurf_r[13])-0.3952847075210473*(alphaSurf_r[8]+alphaSurf_r[7])-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alphaSurf_l[0]-0.3952847075210473*(alphaSurf_l[8]+alphaSurf_l[7]) > 0) { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fc); 
  } 
  if (0.3535533905932737*alphaSurf_r[0]-0.3952847075210473*(alphaSurf_r[8]+alphaSurf_r[7]) > 0) { 
    fUpwindQuad_r[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fr); 
  } 
  if ((-0.5303300858899104*(alphaSurf_l[14]+alphaSurf_l[13]))-0.3952847075210473*(alphaSurf_l[8]+alphaSurf_l[7])+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fc); 
  } 
  if ((-0.5303300858899104*(alphaSurf_r[14]+alphaSurf_r[13]))-0.3952847075210473*(alphaSurf_r[8]+alphaSurf_r[7])+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[17]-0.4242640687119285*alphaSurf_l[14]+0.5303300858899104*alphaSurf_l[13]-0.5303300858899104*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[6]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[17]-0.4242640687119285*alphaSurf_r[14]+0.5303300858899104*alphaSurf_r[13]-0.5303300858899104*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[6]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fr); 
  } 
  if ((-0.5303300858899104*alphaSurf_l[11])+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]+0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fc); 
  } 
  if ((-0.5303300858899104*alphaSurf_r[11])+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]+0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_r[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[17])+0.4242640687119285*alphaSurf_l[14]-0.5303300858899104*(alphaSurf_l[13]+alphaSurf_l[11])+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[6]+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[17])+0.4242640687119285*alphaSurf_r[14]-0.5303300858899104*(alphaSurf_r[13]+alphaSurf_r[11])+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[6]+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fr); 
  } 
  if ((-0.5692099788303082*alphaSurf_l[18])+0.5692099788303082*alphaSurf_l[17]-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[6]-0.6363961030678926*(alphaSurf_l[5]+alphaSurf_l[4])-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fl); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fc); 
  } 
  if ((-0.5692099788303082*alphaSurf_r[18])+0.5692099788303082*alphaSurf_r[17]-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[6]-0.6363961030678926*(alphaSurf_r[5]+alphaSurf_r[4])-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_r[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fr); 
  } 
  if (0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[2]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fc); 
  } 
  if (0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[2]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_r[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fr); 
  } 
  if (0.5692099788303082*alphaSurf_l[18]-0.5692099788303082*alphaSurf_l[17]+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[6]+0.6363961030678926*alphaSurf_l[5]-0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*alphaSurf_l[2]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fl); 
  } else { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fc); 
  } 
  if (0.5692099788303082*alphaSurf_r[18]-0.5692099788303082*alphaSurf_r[17]+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[6]+0.6363961030678926*alphaSurf_r[5]-0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*alphaSurf_r[2]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_r[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[18]+0.5303300858899104*alphaSurf_l[14]-0.4242640687119285*alphaSurf_l[13]-0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[5]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fl); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[18]+0.5303300858899104*alphaSurf_r[14]-0.4242640687119285*alphaSurf_r[13]-0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[5]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_r[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fr); 
  } 
  if ((-0.5303300858899104*alphaSurf_l[12])-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fl); 
  } else { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fc); 
  } 
  if ((-0.5303300858899104*alphaSurf_r[12])-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_r[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[18])-0.5303300858899104*alphaSurf_l[14]+0.4242640687119285*alphaSurf_l[13]-0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[5]+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[18])-0.5303300858899104*alphaSurf_r[14]+0.4242640687119285*alphaSurf_r[13]-0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[5]+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_r[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fr); 
  } 
  if ((-0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17]))-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5])+0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fl); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fc); 
  } 
  if ((-0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17]))-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5])+0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_r[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fr); 
  } 
  if (0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fl); 
  } else { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fc); 
  } 
  if (0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_r[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fr); 
  } 
  if (0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17])+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5]+alphaSurf_l[4])+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fc); 
  } 
  if (0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17])+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5]+alphaSurf_r[4])+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_r[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alphaSurf_l[18]*fUpwind_l[18]+0.3535533905932737*alphaSurf_l[17]*fUpwind_l[17]+0.3535533905932737*alphaSurf_l[14]*fUpwind_l[14]+0.3535533905932737*alphaSurf_l[13]*fUpwind_l[13]+0.3535533905932737*alphaSurf_l[12]*fUpwind_l[12]+0.3535533905932737*alphaSurf_l[11]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[10]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[8]*fUpwind_l[8]+0.3535533905932737*alphaSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[6]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[5]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[4]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[3]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[2]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alphaSurf_l[14]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[14]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[10]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[13]+0.3535533905932737*alphaSurf_l[8]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[8]*alphaSurf_l[12]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[11]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alphaSurf_l[10]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.3162277660168379*alphaSurf_l[10]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[10]*alphaSurf_l[18]+0.3535533905932737*alphaSurf_l[13]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[13]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[14]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[12]+0.3535533905932737*alphaSurf_l[7]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[7]*alphaSurf_l[11]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[10]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[8]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[2]; 
  Ghat_l[3] = 0.3162277660168379*alphaSurf_l[10]*fUpwind_l[19]+0.3535533905932737*alphaSurf_l[12]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[12]*alphaSurf_l[18]+0.3535533905932737*alphaSurf_l[11]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[11]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[15]+0.3535533905932737*alphaSurf_l[8]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[8]*alphaSurf_l[14]+0.3535533905932737*alphaSurf_l[7]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[7]*alphaSurf_l[13]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[4]*alphaSurf_l[10]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[9]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[3]; 
  Ghat_l[4] = 0.2828427124746191*alphaSurf_l[17]*fUpwind_l[18]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[18]+0.2828427124746191*fUpwind_l[17]*alphaSurf_l[18]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[10]*alphaSurf_l[14]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[10]*alphaSurf_l[13]+0.2828427124746191*alphaSurf_l[11]*fUpwind_l[12]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[12]+0.2828427124746191*fUpwind_l[11]*alphaSurf_l[12]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[12]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[11]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[10]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[2]; 
  Ghat_l[5] = 0.2828427124746191*alphaSurf_l[17]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[19]+0.3535533905932737*alphaSurf_l[8]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[8]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[16]+0.2828427124746191*alphaSurf_l[13]*fUpwind_l[15]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[15]+0.3535533905932737*alphaSurf_l[12]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[12]*alphaSurf_l[14]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[13]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[10]*alphaSurf_l[11]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[10]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[3]; 
  Ghat_l[6] = 0.2828427124746191*alphaSurf_l[18]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[18]+0.3535533905932737*alphaSurf_l[7]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[7]*alphaSurf_l[17]+0.2828427124746191*alphaSurf_l[14]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[15]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[14]+0.3535533905932737*alphaSurf_l[11]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[11]*alphaSurf_l[13]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[10]*alphaSurf_l[12]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[10]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[8]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[3]; 
  Ghat_l[7] = 0.3162277660168379*alphaSurf_l[18]*fUpwind_l[18]+0.2258769757263128*alphaSurf_l[17]*fUpwind_l[17]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[6]*alphaSurf_l[17]+0.2258769757263128*alphaSurf_l[13]*fUpwind_l[13]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[13]+0.3162277660168379*alphaSurf_l[12]*fUpwind_l[12]+0.2258769757263128*alphaSurf_l[11]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[10]+0.2258769757263128*alphaSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[5]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[4]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.2258769757263128*alphaSurf_l[18]*fUpwind_l[18]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[17]*fUpwind_l[17]+0.2258769757263128*alphaSurf_l[14]*fUpwind_l[14]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[14]+0.2258769757263128*alphaSurf_l[12]*fUpwind_l[12]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[12]+0.3162277660168379*alphaSurf_l[11]*fUpwind_l[11]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[10]+0.2258769757263128*alphaSurf_l[8]*fUpwind_l[8]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[6]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[4]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alphaSurf_l[4]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[18]*fUpwind_l[18]+0.3162277660168379*alphaSurf_l[17]*fUpwind_l[17]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[16]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[15]+0.3162277660168379*alphaSurf_l[14]*fUpwind_l[14]+0.3162277660168379*alphaSurf_l[13]*fUpwind_l[13]+0.3162277660168379*alphaSurf_l[10]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[6]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[5]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[3]; 
  Ghat_l[10] = 0.282842712474619*alphaSurf_l[14]*fUpwind_l[19]+0.282842712474619*alphaSurf_l[13]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[19]+0.282842712474619*alphaSurf_l[11]*fUpwind_l[18]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[18]+0.282842712474619*fUpwind_l[16]*alphaSurf_l[18]+0.282842712474619*fUpwind_l[11]*alphaSurf_l[18]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[18]+0.282842712474619*alphaSurf_l[12]*fUpwind_l[17]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[17]+0.282842712474619*fUpwind_l[15]*alphaSurf_l[17]+0.282842712474619*fUpwind_l[12]*alphaSurf_l[17]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[15]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[14]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[13]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[12]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]*fUpwind_l[10]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[10]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[9]*alphaSurf_l[10]+0.3162277660168379*fUpwind_l[8]*alphaSurf_l[10]+0.3162277660168379*fUpwind_l[7]*alphaSurf_l[10]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[4]; 
  Ghat_l[11] = 0.282842712474619*alphaSurf_l[10]*fUpwind_l[18]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[14]*fUpwind_l[17]+0.2258769757263128*alphaSurf_l[13]*fUpwind_l[17]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[14]*alphaSurf_l[17]+0.2258769757263128*fUpwind_l[13]*alphaSurf_l[17]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[17]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[6]*alphaSurf_l[13]+0.2828427124746191*alphaSurf_l[4]*fUpwind_l[12]+0.2828427124746191*fUpwind_l[4]*alphaSurf_l[12]+0.3162277660168379*alphaSurf_l[8]*fUpwind_l[11]+0.2258769757263128*alphaSurf_l[7]*fUpwind_l[11]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[8]*alphaSurf_l[11]+0.2258769757263128*fUpwind_l[7]*alphaSurf_l[11]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[4]; 
  Ghat_l[12] = 0.2258769757263128*alphaSurf_l[14]*fUpwind_l[18]+0.3162277660168379*alphaSurf_l[13]*fUpwind_l[18]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[18]+0.2258769757263128*fUpwind_l[14]*alphaSurf_l[18]+0.3162277660168379*fUpwind_l[13]*alphaSurf_l[18]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[18]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[17]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[17]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[14]+0.2258769757263128*alphaSurf_l[8]*fUpwind_l[12]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[12]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[12]+0.2258769757263128*fUpwind_l[8]*alphaSurf_l[12]+0.3162277660168379*fUpwind_l[7]*alphaSurf_l[12]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[12]+0.2828427124746191*alphaSurf_l[4]*fUpwind_l[11]+0.2828427124746191*fUpwind_l[4]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[4]; 
  Ghat_l[13] = 0.282842712474619*alphaSurf_l[10]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[12]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[12]*alphaSurf_l[18]+0.2258769757263128*alphaSurf_l[11]*fUpwind_l[17]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[16]*alphaSurf_l[17]+0.2258769757263128*fUpwind_l[11]*alphaSurf_l[17]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[17]+0.2828427124746191*alphaSurf_l[5]*fUpwind_l[15]+0.2258769757263128*alphaSurf_l[7]*fUpwind_l[13]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[9]*alphaSurf_l[13]+0.2258769757263128*fUpwind_l[7]*alphaSurf_l[13]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[13]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[6]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[5]; 
  Ghat_l[14] = 0.282842712474619*alphaSurf_l[10]*fUpwind_l[19]+0.2258769757263128*alphaSurf_l[12]*fUpwind_l[18]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[15]*alphaSurf_l[18]+0.2258769757263128*fUpwind_l[12]*alphaSurf_l[18]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[11]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[11]*alphaSurf_l[17]+0.2828427124746191*alphaSurf_l[6]*fUpwind_l[16]+0.2258769757263128*alphaSurf_l[8]*fUpwind_l[14]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[9]*alphaSurf_l[14]+0.2258769757263128*fUpwind_l[8]*alphaSurf_l[14]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[14]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[12]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[6]; 
  Ghat_l[15] = 0.3162277660168379*alphaSurf_l[11]*fUpwind_l[19]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[14]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[14]*alphaSurf_l[18]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[17]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[17]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[15]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[15]+0.2828427124746191*alphaSurf_l[5]*fUpwind_l[13]+0.2828427124746191*fUpwind_l[5]*alphaSurf_l[13]+0.3162277660168379*alphaSurf_l[6]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[6]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[3]*alphaSurf_l[5]; 
  Ghat_l[16] = 0.3162277660168379*alphaSurf_l[12]*fUpwind_l[19]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[19]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[18]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[13]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[13]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[8]*fUpwind_l[16]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[16]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[15]+0.2828427124746191*alphaSurf_l[6]*fUpwind_l[14]+0.2828427124746191*fUpwind_l[6]*alphaSurf_l[14]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[3]*alphaSurf_l[6]; 
  Ghat_l[17] = 0.2529822128134704*alphaSurf_l[18]*fUpwind_l[19]+0.2828427124746191*alphaSurf_l[5]*fUpwind_l[19]+0.2828427124746191*alphaSurf_l[4]*fUpwind_l[18]+0.2828427124746191*fUpwind_l[4]*alphaSurf_l[18]+0.3162277660168379*alphaSurf_l[8]*fUpwind_l[17]+0.2258769757263128*alphaSurf_l[7]*fUpwind_l[17]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[9]*alphaSurf_l[17]+0.3162277660168379*fUpwind_l[8]*alphaSurf_l[17]+0.2258769757263128*fUpwind_l[7]*alphaSurf_l[17]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[13]*fUpwind_l[16]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[15]+0.3162277660168379*alphaSurf_l[11]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[11]*alphaSurf_l[14]+0.2258769757263128*alphaSurf_l[11]*fUpwind_l[13]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[13]+0.2258769757263128*fUpwind_l[11]*alphaSurf_l[13]+0.3535533905932737*fUpwind_l[2]*alphaSurf_l[13]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[12]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[12]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[1]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[1]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[5]; 
  Ghat_l[18] = 0.2529822128134704*alphaSurf_l[17]*fUpwind_l[19]+0.2828427124746191*alphaSurf_l[6]*fUpwind_l[19]+0.2258769757263128*alphaSurf_l[8]*fUpwind_l[18]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[18]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[9]*alphaSurf_l[18]+0.2258769757263128*fUpwind_l[8]*alphaSurf_l[18]+0.3162277660168379*fUpwind_l[7]*alphaSurf_l[18]+0.3535533905932737*fUpwind_l[0]*alphaSurf_l[18]+0.2828427124746191*alphaSurf_l[4]*fUpwind_l[17]+0.2828427124746191*fUpwind_l[4]*alphaSurf_l[17]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[14]*fUpwind_l[15]+0.2258769757263128*alphaSurf_l[12]*fUpwind_l[14]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[14]+0.2258769757263128*fUpwind_l[12]*alphaSurf_l[14]+0.3535533905932737*fUpwind_l[1]*alphaSurf_l[14]+0.3162277660168379*alphaSurf_l[12]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[12]*alphaSurf_l[13]+0.3535533905932737*alphaSurf_l[3]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[3]*alphaSurf_l[12]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[11]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[2]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[2]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[5]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[5]*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[4]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[4]*alphaSurf_l[6]; 
  Ghat_l[19] = 0.3162277660168379*alphaSurf_l[8]*fUpwind_l[19]+0.3162277660168379*alphaSurf_l[7]*fUpwind_l[19]+0.3535533905932737*alphaSurf_l[0]*fUpwind_l[19]+0.2529822128134704*alphaSurf_l[17]*fUpwind_l[18]+0.2828427124746191*alphaSurf_l[6]*fUpwind_l[18]+0.2529822128134704*fUpwind_l[17]*alphaSurf_l[18]+0.2828427124746191*fUpwind_l[6]*alphaSurf_l[18]+0.2828427124746191*alphaSurf_l[5]*fUpwind_l[17]+0.2828427124746191*fUpwind_l[5]*alphaSurf_l[17]+0.3162277660168379*alphaSurf_l[12]*fUpwind_l[16]+0.3535533905932737*alphaSurf_l[1]*fUpwind_l[16]+0.3162277660168379*alphaSurf_l[11]*fUpwind_l[15]+0.3535533905932737*alphaSurf_l[2]*fUpwind_l[15]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[14]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[14]+0.282842712474619*alphaSurf_l[10]*fUpwind_l[13]+0.282842712474619*fUpwind_l[10]*alphaSurf_l[13]+0.3162277660168379*alphaSurf_l[3]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[3]*alphaSurf_l[10]+0.3535533905932737*alphaSurf_l[4]*fUpwind_l[9]+0.3162277660168379*alphaSurf_l[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alphaSurf_l[6]; 

  Ghat_r[0] = 0.3535533905932737*alphaSurf_r[18]*fUpwind_r[18]+0.3535533905932737*alphaSurf_r[17]*fUpwind_r[17]+0.3535533905932737*alphaSurf_r[14]*fUpwind_r[14]+0.3535533905932737*alphaSurf_r[13]*fUpwind_r[13]+0.3535533905932737*alphaSurf_r[12]*fUpwind_r[12]+0.3535533905932737*alphaSurf_r[11]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[10]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[8]*fUpwind_r[8]+0.3535533905932737*alphaSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[6]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[5]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[4]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[3]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[2]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alphaSurf_r[14]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[14]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[10]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[13]+0.3535533905932737*alphaSurf_r[8]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[8]*alphaSurf_r[12]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[11]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alphaSurf_r[10]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.3162277660168379*alphaSurf_r[10]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[10]*alphaSurf_r[18]+0.3535533905932737*alphaSurf_r[13]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[13]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[14]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[12]+0.3535533905932737*alphaSurf_r[7]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[7]*alphaSurf_r[11]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[10]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[8]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[2]; 
  Ghat_r[3] = 0.3162277660168379*alphaSurf_r[10]*fUpwind_r[19]+0.3535533905932737*alphaSurf_r[12]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[12]*alphaSurf_r[18]+0.3535533905932737*alphaSurf_r[11]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[11]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[15]+0.3535533905932737*alphaSurf_r[8]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[8]*alphaSurf_r[14]+0.3535533905932737*alphaSurf_r[7]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[7]*alphaSurf_r[13]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[4]*alphaSurf_r[10]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[9]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[3]; 
  Ghat_r[4] = 0.2828427124746191*alphaSurf_r[17]*fUpwind_r[18]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[18]+0.2828427124746191*fUpwind_r[17]*alphaSurf_r[18]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[10]*alphaSurf_r[14]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[10]*alphaSurf_r[13]+0.2828427124746191*alphaSurf_r[11]*fUpwind_r[12]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[12]+0.2828427124746191*fUpwind_r[11]*alphaSurf_r[12]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[12]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[11]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[10]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[2]; 
  Ghat_r[5] = 0.2828427124746191*alphaSurf_r[17]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[19]+0.3535533905932737*alphaSurf_r[8]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[8]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[16]+0.2828427124746191*alphaSurf_r[13]*fUpwind_r[15]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[15]+0.3535533905932737*alphaSurf_r[12]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[12]*alphaSurf_r[14]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[13]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[10]*alphaSurf_r[11]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[10]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[4]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[3]; 
  Ghat_r[6] = 0.2828427124746191*alphaSurf_r[18]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[18]+0.3535533905932737*alphaSurf_r[7]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[7]*alphaSurf_r[17]+0.2828427124746191*alphaSurf_r[14]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[15]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[14]+0.3535533905932737*alphaSurf_r[11]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[11]*alphaSurf_r[13]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[10]*alphaSurf_r[12]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[10]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[8]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[3]; 
  Ghat_r[7] = 0.3162277660168379*alphaSurf_r[18]*fUpwind_r[18]+0.2258769757263128*alphaSurf_r[17]*fUpwind_r[17]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[6]*alphaSurf_r[17]+0.2258769757263128*alphaSurf_r[13]*fUpwind_r[13]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[13]+0.3162277660168379*alphaSurf_r[12]*fUpwind_r[12]+0.2258769757263128*alphaSurf_r[11]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[10]+0.2258769757263128*alphaSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[5]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[4]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.2258769757263128*alphaSurf_r[18]*fUpwind_r[18]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[17]*fUpwind_r[17]+0.2258769757263128*alphaSurf_r[14]*fUpwind_r[14]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[14]+0.2258769757263128*alphaSurf_r[12]*fUpwind_r[12]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[12]+0.3162277660168379*alphaSurf_r[11]*fUpwind_r[11]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[10]+0.2258769757263128*alphaSurf_r[8]*fUpwind_r[8]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[6]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[4]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alphaSurf_r[4]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[18]*fUpwind_r[18]+0.3162277660168379*alphaSurf_r[17]*fUpwind_r[17]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[16]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[15]+0.3162277660168379*alphaSurf_r[14]*fUpwind_r[14]+0.3162277660168379*alphaSurf_r[13]*fUpwind_r[13]+0.3162277660168379*alphaSurf_r[10]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[6]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[5]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[3]; 
  Ghat_r[10] = 0.282842712474619*alphaSurf_r[14]*fUpwind_r[19]+0.282842712474619*alphaSurf_r[13]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[19]+0.282842712474619*alphaSurf_r[11]*fUpwind_r[18]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[18]+0.282842712474619*fUpwind_r[16]*alphaSurf_r[18]+0.282842712474619*fUpwind_r[11]*alphaSurf_r[18]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[18]+0.282842712474619*alphaSurf_r[12]*fUpwind_r[17]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[17]+0.282842712474619*fUpwind_r[15]*alphaSurf_r[17]+0.282842712474619*fUpwind_r[12]*alphaSurf_r[17]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[15]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[14]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[13]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[12]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]*fUpwind_r[10]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[10]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[9]*alphaSurf_r[10]+0.3162277660168379*fUpwind_r[8]*alphaSurf_r[10]+0.3162277660168379*fUpwind_r[7]*alphaSurf_r[10]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[4]; 
  Ghat_r[11] = 0.282842712474619*alphaSurf_r[10]*fUpwind_r[18]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[14]*fUpwind_r[17]+0.2258769757263128*alphaSurf_r[13]*fUpwind_r[17]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[14]*alphaSurf_r[17]+0.2258769757263128*fUpwind_r[13]*alphaSurf_r[17]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[17]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[6]*alphaSurf_r[13]+0.2828427124746191*alphaSurf_r[4]*fUpwind_r[12]+0.2828427124746191*fUpwind_r[4]*alphaSurf_r[12]+0.3162277660168379*alphaSurf_r[8]*fUpwind_r[11]+0.2258769757263128*alphaSurf_r[7]*fUpwind_r[11]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[8]*alphaSurf_r[11]+0.2258769757263128*fUpwind_r[7]*alphaSurf_r[11]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[4]; 
  Ghat_r[12] = 0.2258769757263128*alphaSurf_r[14]*fUpwind_r[18]+0.3162277660168379*alphaSurf_r[13]*fUpwind_r[18]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[18]+0.2258769757263128*fUpwind_r[14]*alphaSurf_r[18]+0.3162277660168379*fUpwind_r[13]*alphaSurf_r[18]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[18]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[17]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[17]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[14]+0.2258769757263128*alphaSurf_r[8]*fUpwind_r[12]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[12]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[12]+0.2258769757263128*fUpwind_r[8]*alphaSurf_r[12]+0.3162277660168379*fUpwind_r[7]*alphaSurf_r[12]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[12]+0.2828427124746191*alphaSurf_r[4]*fUpwind_r[11]+0.2828427124746191*fUpwind_r[4]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[4]; 
  Ghat_r[13] = 0.282842712474619*alphaSurf_r[10]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[12]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[12]*alphaSurf_r[18]+0.2258769757263128*alphaSurf_r[11]*fUpwind_r[17]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[16]*alphaSurf_r[17]+0.2258769757263128*fUpwind_r[11]*alphaSurf_r[17]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[17]+0.2828427124746191*alphaSurf_r[5]*fUpwind_r[15]+0.2258769757263128*alphaSurf_r[7]*fUpwind_r[13]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[9]*alphaSurf_r[13]+0.2258769757263128*fUpwind_r[7]*alphaSurf_r[13]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[13]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[6]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[5]; 
  Ghat_r[14] = 0.282842712474619*alphaSurf_r[10]*fUpwind_r[19]+0.2258769757263128*alphaSurf_r[12]*fUpwind_r[18]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[15]*alphaSurf_r[18]+0.2258769757263128*fUpwind_r[12]*alphaSurf_r[18]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[11]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[11]*alphaSurf_r[17]+0.2828427124746191*alphaSurf_r[6]*fUpwind_r[16]+0.2258769757263128*alphaSurf_r[8]*fUpwind_r[14]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[9]*alphaSurf_r[14]+0.2258769757263128*fUpwind_r[8]*alphaSurf_r[14]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[14]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[12]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[6]; 
  Ghat_r[15] = 0.3162277660168379*alphaSurf_r[11]*fUpwind_r[19]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[14]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[14]*alphaSurf_r[18]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[17]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[17]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[15]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[15]+0.2828427124746191*alphaSurf_r[5]*fUpwind_r[13]+0.2828427124746191*fUpwind_r[5]*alphaSurf_r[13]+0.3162277660168379*alphaSurf_r[6]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[6]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[3]*alphaSurf_r[5]; 
  Ghat_r[16] = 0.3162277660168379*alphaSurf_r[12]*fUpwind_r[19]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[19]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[18]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[13]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[13]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[8]*fUpwind_r[16]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[16]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[15]+0.2828427124746191*alphaSurf_r[6]*fUpwind_r[14]+0.2828427124746191*fUpwind_r[6]*alphaSurf_r[14]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[3]*alphaSurf_r[6]; 
  Ghat_r[17] = 0.2529822128134704*alphaSurf_r[18]*fUpwind_r[19]+0.2828427124746191*alphaSurf_r[5]*fUpwind_r[19]+0.2828427124746191*alphaSurf_r[4]*fUpwind_r[18]+0.2828427124746191*fUpwind_r[4]*alphaSurf_r[18]+0.3162277660168379*alphaSurf_r[8]*fUpwind_r[17]+0.2258769757263128*alphaSurf_r[7]*fUpwind_r[17]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[9]*alphaSurf_r[17]+0.3162277660168379*fUpwind_r[8]*alphaSurf_r[17]+0.2258769757263128*fUpwind_r[7]*alphaSurf_r[17]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[13]*fUpwind_r[16]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[15]+0.3162277660168379*alphaSurf_r[11]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[11]*alphaSurf_r[14]+0.2258769757263128*alphaSurf_r[11]*fUpwind_r[13]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[13]+0.2258769757263128*fUpwind_r[11]*alphaSurf_r[13]+0.3535533905932737*fUpwind_r[2]*alphaSurf_r[13]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[12]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[12]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[1]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[1]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[6]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[6]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[5]; 
  Ghat_r[18] = 0.2529822128134704*alphaSurf_r[17]*fUpwind_r[19]+0.2828427124746191*alphaSurf_r[6]*fUpwind_r[19]+0.2258769757263128*alphaSurf_r[8]*fUpwind_r[18]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[18]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[9]*alphaSurf_r[18]+0.2258769757263128*fUpwind_r[8]*alphaSurf_r[18]+0.3162277660168379*fUpwind_r[7]*alphaSurf_r[18]+0.3535533905932737*fUpwind_r[0]*alphaSurf_r[18]+0.2828427124746191*alphaSurf_r[4]*fUpwind_r[17]+0.2828427124746191*fUpwind_r[4]*alphaSurf_r[17]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[14]*fUpwind_r[15]+0.2258769757263128*alphaSurf_r[12]*fUpwind_r[14]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[14]+0.2258769757263128*fUpwind_r[12]*alphaSurf_r[14]+0.3535533905932737*fUpwind_r[1]*alphaSurf_r[14]+0.3162277660168379*alphaSurf_r[12]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[12]*alphaSurf_r[13]+0.3535533905932737*alphaSurf_r[3]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[3]*alphaSurf_r[12]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[11]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[2]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[2]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[5]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[5]*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[4]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[4]*alphaSurf_r[6]; 
  Ghat_r[19] = 0.3162277660168379*alphaSurf_r[8]*fUpwind_r[19]+0.3162277660168379*alphaSurf_r[7]*fUpwind_r[19]+0.3535533905932737*alphaSurf_r[0]*fUpwind_r[19]+0.2529822128134704*alphaSurf_r[17]*fUpwind_r[18]+0.2828427124746191*alphaSurf_r[6]*fUpwind_r[18]+0.2529822128134704*fUpwind_r[17]*alphaSurf_r[18]+0.2828427124746191*fUpwind_r[6]*alphaSurf_r[18]+0.2828427124746191*alphaSurf_r[5]*fUpwind_r[17]+0.2828427124746191*fUpwind_r[5]*alphaSurf_r[17]+0.3162277660168379*alphaSurf_r[12]*fUpwind_r[16]+0.3535533905932737*alphaSurf_r[1]*fUpwind_r[16]+0.3162277660168379*alphaSurf_r[11]*fUpwind_r[15]+0.3535533905932737*alphaSurf_r[2]*fUpwind_r[15]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[14]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[14]+0.282842712474619*alphaSurf_r[10]*fUpwind_r[13]+0.282842712474619*fUpwind_r[10]*alphaSurf_r[13]+0.3162277660168379*alphaSurf_r[3]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[3]*alphaSurf_r[10]+0.3535533905932737*alphaSurf_r[4]*fUpwind_r[9]+0.3162277660168379*alphaSurf_r[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alphaSurf_r[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx1; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx1; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx1; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx1; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx1; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx1; 
  out[12] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx1; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx1; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx1; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx1; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx1; 
  out[17] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx1; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx1; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx1; 
  out[20] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx1; 
  out[21] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx1; 
  out[22] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dx1; 
  out[23] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx1; 
  out[24] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx1; 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx1; 
  out[26] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dx1; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx1; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx1; 
  out[29] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx1; 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dx1; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx1; 
  out[32] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx1; 
  out[33] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dx1; 
  out[34] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx1; 
  out[35] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx1; 
  out[36] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dx1; 
  out[37] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dx1; 
  out[38] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dx1; 
  out[39] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dx1; 
  out[40] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx1; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx1; 
  out[42] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dx1; 
  out[43] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dx1; 
  out[44] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dx1; 
  out[45] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dx1; 
  out[46] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dx1; 
  out[47] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dx1; 

} 
