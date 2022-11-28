#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *bb_grad_u, const double *p_force, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:       Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dx1 = 2.0/dxv[1]; 
  const double dx2 = 2.0/dxv[2]; 
  const double dv1par = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  double alphaSurf[20] = {0.0}; 
  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};;
  double Ghat[20] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force[7]; 
  alphaSurf[8] = (-1.0*bb_grad_u[8]*wvpar)-0.5*bb_grad_u[8]*dvpar+p_force[8]; 
  alphaSurf[9] = (-1.0*bb_grad_u[9]*wvpar)-0.5*bb_grad_u[9]*dvpar+p_force[9]; 
  alphaSurf[10] = (-1.0*bb_grad_u[10]*wvpar)-0.5*bb_grad_u[10]*dvpar+p_force[10]; 
  alphaSurf[11] = (-1.0*bb_grad_u[11]*wvpar)-0.5*bb_grad_u[11]*dvpar+p_force[11]; 
  alphaSurf[12] = (-1.0*bb_grad_u[12]*wvpar)-0.5*bb_grad_u[12]*dvpar+p_force[12]; 
  alphaSurf[13] = (-1.0*bb_grad_u[13]*wvpar)-0.5*bb_grad_u[13]*dvpar+p_force[13]; 
  alphaSurf[14] = (-1.0*bb_grad_u[14]*wvpar)-0.5*bb_grad_u[14]*dvpar+p_force[14]; 
  alphaSurf[15] = (-1.0*bb_grad_u[15]*wvpar)-0.5*bb_grad_u[15]*dvpar+p_force[15]; 
  alphaSurf[16] = (-1.0*bb_grad_u[16]*wvpar)-0.5*bb_grad_u[16]*dvpar+p_force[16]; 
  alphaSurf[17] = (-1.0*bb_grad_u[17]*wvpar)-0.5*bb_grad_u[17]*dvpar+p_force[17]; 
  alphaSurf[18] = (-1.0*bb_grad_u[18]*wvpar)-0.5*bb_grad_u[18]*dvpar+p_force[18]; 
  alphaSurf[19] = (-1.0*bb_grad_u[19]*wvpar)-0.5*bb_grad_u[19]*dvpar+p_force[19]; 

  if (0.5692099788303082*(alphaSurf[19]+alphaSurf[18]+alphaSurf[17])-0.4242640687119281*(alphaSurf[16]+alphaSurf[15])-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*(alphaSurf[6]+alphaSurf[5]+alphaSurf[4])-0.4743416490252568*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if ((-0.711512473537885*alphaSurf[19])+0.5303300858899104*(alphaSurf[16]+alphaSurf[15])-0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[4]-0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (0.5692099788303082*alphaSurf[19]-0.5692099788303082*(alphaSurf[18]+alphaSurf[17])-0.4242640687119281*(alphaSurf[16]+alphaSurf[15])+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*(alphaSurf[6]+alphaSurf[5])+0.6363961030678926*alphaSurf[4]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if ((-0.711512473537885*alphaSurf[18])-0.4242640687119281*alphaSurf[15]+0.5303300858899104*alphaSurf[14]-0.4242640687119285*alphaSurf[13]+0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]+0.6363961030678926*alphaSurf[5]-0.4743416490252568*(alphaSurf[3]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (0.5303300858899104*(alphaSurf[15]+alphaSurf[12])-0.3952847075210473*(alphaSurf[9]+alphaSurf[8])+0.3162277660168379*alphaSurf[7]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if (0.711512473537885*alphaSurf[18]-0.4242640687119281*alphaSurf[15]-0.5303300858899104*alphaSurf[14]+0.4242640687119285*alphaSurf[13]+0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]-0.6363961030678926*alphaSurf[5]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.5692099788303082*alphaSurf[19])+0.5692099788303082*alphaSurf[18]-0.5692099788303082*alphaSurf[17]+0.4242640687119281*alphaSurf[16]-0.4242640687119281*alphaSurf[15]-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[6]+0.6363961030678926*alphaSurf[5]-0.6363961030678926*alphaSurf[4]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*alphaSurf[2]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if (0.711512473537885*alphaSurf[19]-0.5303300858899104*alphaSurf[16]+0.5303300858899104*alphaSurf[15]-0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[4]+0.4743416490252568*alphaSurf[2]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if ((-0.5692099788303082*(alphaSurf[19]+alphaSurf[18]))+0.5692099788303082*alphaSurf[17]+0.4242640687119281*alphaSurf[16]-0.4242640687119281*alphaSurf[15]+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[6]-0.6363961030678926*(alphaSurf[5]+alphaSurf[4])+0.4743416490252568*(alphaSurf[3]+alphaSurf[2])-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if ((-0.711512473537885*alphaSurf[17])-0.4242640687119281*alphaSurf[16]-0.4242640687119285*alphaSurf[14]+0.5303300858899104*(alphaSurf[13]+alphaSurf[11])+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]+0.6363961030678926*alphaSurf[6]-0.4743416490252568*(alphaSurf[3]+alphaSurf[2])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fEdge); 
  } 
  if (0.5303300858899104*(alphaSurf[16]+alphaSurf[11])-0.3952847075210473*alphaSurf[9]+0.3162277660168379*alphaSurf[8]-0.3952847075210473*alphaSurf[7]-0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fEdge); 
  } 
  if (0.711512473537885*alphaSurf[17]-0.4242640687119281*alphaSurf[16]+0.4242640687119285*alphaSurf[14]-0.5303300858899104*alphaSurf[13]+0.5303300858899104*alphaSurf[11]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]-0.6363961030678926*alphaSurf[6]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (0.5303300858899104*(alphaSurf[14]+alphaSurf[13])+0.3162277660168379*alphaSurf[9]-0.3952847075210473*(alphaSurf[8]+alphaSurf[7])-0.4743416490252568*alphaSurf[3]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fEdge); 
  } 
  if (0.3535533905932737*alphaSurf[0]-0.3952847075210473*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7]) > 0) { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fEdge); 
  } 
  if ((-0.5303300858899104*(alphaSurf[14]+alphaSurf[13]))+0.3162277660168379*alphaSurf[9]-0.3952847075210473*(alphaSurf[8]+alphaSurf[7])+0.4743416490252568*alphaSurf[3]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (0.711512473537885*alphaSurf[17]+0.4242640687119281*alphaSurf[16]-0.4242640687119285*alphaSurf[14]+0.5303300858899104*alphaSurf[13]-0.5303300858899104*alphaSurf[11]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]-0.6363961030678926*alphaSurf[6]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fEdge); 
  } 
  if ((-0.5303300858899104*(alphaSurf[16]+alphaSurf[11]))-0.3952847075210473*alphaSurf[9]+0.3162277660168379*alphaSurf[8]-0.3952847075210473*alphaSurf[7]+0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fEdge); 
  } 
  if ((-0.711512473537885*alphaSurf[17])+0.4242640687119281*alphaSurf[16]+0.4242640687119285*alphaSurf[14]-0.5303300858899104*(alphaSurf[13]+alphaSurf[11])+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]+0.6363961030678926*alphaSurf[6]+0.4743416490252568*(alphaSurf[3]+alphaSurf[2])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fEdge); 
  } 
  if ((-0.5692099788303082*(alphaSurf[19]+alphaSurf[18]))+0.5692099788303082*alphaSurf[17]-0.4242640687119281*alphaSurf[16]+0.4242640687119281*alphaSurf[15]-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[6]-0.6363961030678926*(alphaSurf[5]+alphaSurf[4])-0.4743416490252568*(alphaSurf[3]+alphaSurf[2])+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fEdge); 
  } 
  if (0.711512473537885*alphaSurf[19]+0.5303300858899104*alphaSurf[16]-0.5303300858899104*alphaSurf[15]+0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[4]-0.4743416490252568*alphaSurf[2]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fEdge); 
  } 
  if ((-0.5692099788303082*alphaSurf[19])+0.5692099788303082*alphaSurf[18]-0.5692099788303082*alphaSurf[17]-0.4242640687119281*alphaSurf[16]+0.4242640687119281*alphaSurf[15]+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[6]+0.6363961030678926*alphaSurf[5]-0.6363961030678926*alphaSurf[4]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*alphaSurf[2]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fEdge); 
  } 
  if (0.711512473537885*alphaSurf[18]+0.4242640687119281*alphaSurf[15]+0.5303300858899104*alphaSurf[14]-0.4242640687119285*alphaSurf[13]-0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]-0.6363961030678926*alphaSurf[5]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fEdge); 
  } 
  if ((-0.5303300858899104*(alphaSurf[15]+alphaSurf[12]))-0.3952847075210473*(alphaSurf[9]+alphaSurf[8])+0.3162277660168379*alphaSurf[7]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fSkin); 
  } else { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fEdge); 
  } 
  if ((-0.711512473537885*alphaSurf[18])+0.4242640687119281*alphaSurf[15]-0.5303300858899104*alphaSurf[14]+0.4242640687119285*alphaSurf[13]-0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]+0.6363961030678926*alphaSurf[5]+0.4743416490252568*(alphaSurf[3]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if (0.5692099788303082*alphaSurf[19]-0.5692099788303082*(alphaSurf[18]+alphaSurf[17])+0.4242640687119281*(alphaSurf[16]+alphaSurf[15])-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*(alphaSurf[6]+alphaSurf[5])+0.6363961030678926*alphaSurf[4]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fEdge); 
  } 
  if ((-0.711512473537885*alphaSurf[19])-0.5303300858899104*(alphaSurf[16]+alphaSurf[15])+0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[4]+0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fSkin); 
  } else { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fEdge); 
  } 
  if (0.5692099788303082*(alphaSurf[19]+alphaSurf[18]+alphaSurf[17])+0.4242640687119281*(alphaSurf[16]+alphaSurf[15])+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*(alphaSurf[6]+alphaSurf[5]+alphaSurf[4])+0.4743416490252568*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaSurf[19]*fUpwind[19]+0.3535533905932737*alphaSurf[18]*fUpwind[18]+0.3535533905932737*alphaSurf[17]*fUpwind[17]+0.3535533905932737*alphaSurf[16]*fUpwind[16]+0.3535533905932737*alphaSurf[15]*fUpwind[15]+0.3535533905932737*alphaSurf[14]*fUpwind[14]+0.3535533905932737*alphaSurf[13]*fUpwind[13]+0.3535533905932737*alphaSurf[12]*fUpwind[12]+0.3535533905932737*alphaSurf[11]*fUpwind[11]+0.3535533905932737*alphaSurf[10]*fUpwind[10]+0.3535533905932737*alphaSurf[9]*fUpwind[9]+0.3535533905932737*alphaSurf[8]*fUpwind[8]+0.3535533905932737*alphaSurf[7]*fUpwind[7]+0.3535533905932737*alphaSurf[6]*fUpwind[6]+0.3535533905932737*alphaSurf[5]*fUpwind[5]+0.3535533905932737*alphaSurf[4]*fUpwind[4]+0.3535533905932737*alphaSurf[3]*fUpwind[3]+0.3535533905932737*alphaSurf[2]*fUpwind[2]+0.3535533905932737*alphaSurf[1]*fUpwind[1]+0.3535533905932737*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaSurf[16]*fUpwind[19]+0.3535533905932737*fUpwind[16]*alphaSurf[19]+0.3535533905932737*alphaSurf[14]*fUpwind[18]+0.3535533905932737*fUpwind[14]*alphaSurf[18]+0.3162277660168379*alphaSurf[10]*fUpwind[17]+0.3162277660168379*fUpwind[10]*alphaSurf[17]+0.3535533905932737*alphaSurf[9]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaSurf[15]+0.3162277660168379*alphaSurf[5]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alphaSurf[13]+0.3535533905932737*alphaSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaSurf[12]+0.3162277660168379*alphaSurf[4]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaSurf[11]+0.3535533905932737*alphaSurf[6]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaSurf[10]+0.3162277660168379*alphaSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaSurf[7]+0.3535533905932737*alphaSurf[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.3535533905932737*alphaSurf[15]*fUpwind[19]+0.3535533905932737*fUpwind[15]*alphaSurf[19]+0.3162277660168379*alphaSurf[10]*fUpwind[18]+0.3162277660168379*fUpwind[10]*alphaSurf[18]+0.3535533905932737*alphaSurf[13]*fUpwind[17]+0.3535533905932737*fUpwind[13]*alphaSurf[17]+0.3535533905932737*alphaSurf[9]*fUpwind[16]+0.3535533905932737*fUpwind[9]*alphaSurf[16]+0.3162277660168379*alphaSurf[6]*fUpwind[14]+0.3162277660168379*fUpwind[6]*alphaSurf[14]+0.3162277660168379*alphaSurf[4]*fUpwind[12]+0.3162277660168379*fUpwind[4]*alphaSurf[12]+0.3535533905932737*alphaSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaSurf[11]+0.3535533905932737*alphaSurf[5]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaSurf[10]+0.3162277660168379*alphaSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaSurf[8]+0.3535533905932737*alphaSurf[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaSurf[2]; 
  Ghat[3] = 0.3162277660168379*alphaSurf[10]*fUpwind[19]+0.3162277660168379*fUpwind[10]*alphaSurf[19]+0.3535533905932737*alphaSurf[12]*fUpwind[18]+0.3535533905932737*fUpwind[12]*alphaSurf[18]+0.3535533905932737*alphaSurf[11]*fUpwind[17]+0.3535533905932737*fUpwind[11]*alphaSurf[17]+0.3162277660168379*alphaSurf[6]*fUpwind[16]+0.3162277660168379*fUpwind[6]*alphaSurf[16]+0.3162277660168379*alphaSurf[5]*fUpwind[15]+0.3162277660168379*fUpwind[5]*alphaSurf[15]+0.3535533905932737*alphaSurf[8]*fUpwind[14]+0.3535533905932737*fUpwind[8]*alphaSurf[14]+0.3535533905932737*alphaSurf[7]*fUpwind[13]+0.3535533905932737*fUpwind[7]*alphaSurf[13]+0.3535533905932737*alphaSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[4]*alphaSurf[10]+0.3162277660168379*alphaSurf[3]*fUpwind[9]+0.3162277660168379*fUpwind[3]*alphaSurf[9]+0.3535533905932737*alphaSurf[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alphaSurf[5]+0.3535533905932737*alphaSurf[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alphaSurf[3]; 
  Ghat[4] = 0.3535533905932737*alphaSurf[9]*fUpwind[19]+0.3535533905932737*fUpwind[9]*alphaSurf[19]+0.2828427124746191*alphaSurf[17]*fUpwind[18]+0.3162277660168379*alphaSurf[6]*fUpwind[18]+0.2828427124746191*fUpwind[17]*alphaSurf[18]+0.3162277660168379*fUpwind[6]*alphaSurf[18]+0.3162277660168379*alphaSurf[5]*fUpwind[17]+0.3162277660168379*fUpwind[5]*alphaSurf[17]+0.3535533905932737*alphaSurf[15]*fUpwind[16]+0.3535533905932737*fUpwind[15]*alphaSurf[16]+0.3162277660168379*alphaSurf[10]*fUpwind[14]+0.3162277660168379*fUpwind[10]*alphaSurf[14]+0.3162277660168379*alphaSurf[10]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaSurf[13]+0.2828427124746191*alphaSurf[11]*fUpwind[12]+0.3162277660168379*alphaSurf[2]*fUpwind[12]+0.2828427124746191*fUpwind[11]*alphaSurf[12]+0.3162277660168379*fUpwind[2]*alphaSurf[12]+0.3162277660168379*alphaSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[1]*alphaSurf[11]+0.3535533905932737*alphaSurf[3]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaSurf[10]+0.3162277660168379*alphaSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaSurf[8]+0.3162277660168379*alphaSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaSurf[7]+0.3535533905932737*alphaSurf[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaSurf[4]+0.3535533905932737*alphaSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaSurf[2]; 
  Ghat[5] = 0.2828427124746191*alphaSurf[17]*fUpwind[19]+0.3162277660168379*alphaSurf[6]*fUpwind[19]+0.2828427124746191*fUpwind[17]*alphaSurf[19]+0.3162277660168379*fUpwind[6]*alphaSurf[19]+0.3535533905932737*alphaSurf[8]*fUpwind[18]+0.3535533905932737*fUpwind[8]*alphaSurf[18]+0.3162277660168379*alphaSurf[4]*fUpwind[17]+0.3162277660168379*fUpwind[4]*alphaSurf[17]+0.3162277660168379*alphaSurf[10]*fUpwind[16]+0.3162277660168379*fUpwind[10]*alphaSurf[16]+0.2828427124746191*alphaSurf[13]*fUpwind[15]+0.3162277660168379*alphaSurf[3]*fUpwind[15]+0.2828427124746191*fUpwind[13]*alphaSurf[15]+0.3162277660168379*fUpwind[3]*alphaSurf[15]+0.3535533905932737*alphaSurf[12]*fUpwind[14]+0.3535533905932737*fUpwind[12]*alphaSurf[14]+0.3162277660168379*alphaSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[1]*alphaSurf[13]+0.3162277660168379*alphaSurf[10]*fUpwind[11]+0.3162277660168379*fUpwind[10]*alphaSurf[11]+0.3535533905932737*alphaSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[2]*alphaSurf[10]+0.3162277660168379*alphaSurf[5]*fUpwind[9]+0.3162277660168379*fUpwind[5]*alphaSurf[9]+0.3162277660168379*alphaSurf[5]*fUpwind[7]+0.3162277660168379*fUpwind[5]*alphaSurf[7]+0.3535533905932737*alphaSurf[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alphaSurf[5]+0.3535533905932737*alphaSurf[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alphaSurf[3]; 
  Ghat[6] = 0.2828427124746191*alphaSurf[18]*fUpwind[19]+0.3162277660168379*alphaSurf[5]*fUpwind[19]+0.2828427124746191*fUpwind[18]*alphaSurf[19]+0.3162277660168379*fUpwind[5]*alphaSurf[19]+0.3162277660168379*alphaSurf[4]*fUpwind[18]+0.3162277660168379*fUpwind[4]*alphaSurf[18]+0.3535533905932737*alphaSurf[7]*fUpwind[17]+0.3535533905932737*fUpwind[7]*alphaSurf[17]+0.2828427124746191*alphaSurf[14]*fUpwind[16]+0.3162277660168379*alphaSurf[3]*fUpwind[16]+0.2828427124746191*fUpwind[14]*alphaSurf[16]+0.3162277660168379*fUpwind[3]*alphaSurf[16]+0.3162277660168379*alphaSurf[10]*fUpwind[15]+0.3162277660168379*fUpwind[10]*alphaSurf[15]+0.3162277660168379*alphaSurf[2]*fUpwind[14]+0.3162277660168379*fUpwind[2]*alphaSurf[14]+0.3535533905932737*alphaSurf[11]*fUpwind[13]+0.3535533905932737*fUpwind[11]*alphaSurf[13]+0.3162277660168379*alphaSurf[10]*fUpwind[12]+0.3162277660168379*fUpwind[10]*alphaSurf[12]+0.3535533905932737*alphaSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[1]*alphaSurf[10]+0.3162277660168379*alphaSurf[6]*fUpwind[9]+0.3162277660168379*fUpwind[6]*alphaSurf[9]+0.3162277660168379*alphaSurf[6]*fUpwind[8]+0.3162277660168379*fUpwind[6]*alphaSurf[8]+0.3535533905932737*alphaSurf[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alphaSurf[6]+0.3535533905932737*alphaSurf[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alphaSurf[3]; 
  Ghat[7] = 0.3162277660168379*alphaSurf[19]*fUpwind[19]+0.3162277660168379*alphaSurf[18]*fUpwind[18]+0.2258769757263128*alphaSurf[17]*fUpwind[17]+0.3535533905932737*alphaSurf[6]*fUpwind[17]+0.3535533905932737*fUpwind[6]*alphaSurf[17]+0.3162277660168379*alphaSurf[15]*fUpwind[15]+0.2258769757263128*alphaSurf[13]*fUpwind[13]+0.3535533905932737*alphaSurf[3]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaSurf[13]+0.3162277660168379*alphaSurf[12]*fUpwind[12]+0.2258769757263128*alphaSurf[11]*fUpwind[11]+0.3535533905932737*alphaSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaSurf[11]+0.3162277660168379*alphaSurf[10]*fUpwind[10]+0.2258769757263128*alphaSurf[7]*fUpwind[7]+0.3535533905932737*alphaSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaSurf[7]+0.3162277660168379*alphaSurf[5]*fUpwind[5]+0.3162277660168379*alphaSurf[4]*fUpwind[4]+0.3162277660168379*alphaSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.3162277660168379*alphaSurf[19]*fUpwind[19]+0.2258769757263128*alphaSurf[18]*fUpwind[18]+0.3535533905932737*alphaSurf[5]*fUpwind[18]+0.3535533905932737*fUpwind[5]*alphaSurf[18]+0.3162277660168379*alphaSurf[17]*fUpwind[17]+0.3162277660168379*alphaSurf[16]*fUpwind[16]+0.2258769757263128*alphaSurf[14]*fUpwind[14]+0.3535533905932737*alphaSurf[3]*fUpwind[14]+0.3535533905932737*fUpwind[3]*alphaSurf[14]+0.2258769757263128*alphaSurf[12]*fUpwind[12]+0.3535533905932737*alphaSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaSurf[12]+0.3162277660168379*alphaSurf[11]*fUpwind[11]+0.3162277660168379*alphaSurf[10]*fUpwind[10]+0.2258769757263128*alphaSurf[8]*fUpwind[8]+0.3535533905932737*alphaSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaSurf[8]+0.3162277660168379*alphaSurf[6]*fUpwind[6]+0.3162277660168379*alphaSurf[4]*fUpwind[4]+0.3162277660168379*alphaSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.2258769757263128*alphaSurf[19]*fUpwind[19]+0.3535533905932737*alphaSurf[4]*fUpwind[19]+0.3535533905932737*fUpwind[4]*alphaSurf[19]+0.3162277660168379*alphaSurf[18]*fUpwind[18]+0.3162277660168379*alphaSurf[17]*fUpwind[17]+0.2258769757263128*alphaSurf[16]*fUpwind[16]+0.3535533905932737*alphaSurf[2]*fUpwind[16]+0.3535533905932737*fUpwind[2]*alphaSurf[16]+0.2258769757263128*alphaSurf[15]*fUpwind[15]+0.3535533905932737*alphaSurf[1]*fUpwind[15]+0.3535533905932737*fUpwind[1]*alphaSurf[15]+0.3162277660168379*alphaSurf[14]*fUpwind[14]+0.3162277660168379*alphaSurf[13]*fUpwind[13]+0.3162277660168379*alphaSurf[10]*fUpwind[10]+0.2258769757263128*alphaSurf[9]*fUpwind[9]+0.3535533905932737*alphaSurf[0]*fUpwind[9]+0.3535533905932737*fUpwind[0]*alphaSurf[9]+0.3162277660168379*alphaSurf[6]*fUpwind[6]+0.3162277660168379*alphaSurf[5]*fUpwind[5]+0.3162277660168379*alphaSurf[3]*fUpwind[3]; 
  Ghat[10] = 0.282842712474619*alphaSurf[14]*fUpwind[19]+0.282842712474619*alphaSurf[13]*fUpwind[19]+0.3162277660168379*alphaSurf[3]*fUpwind[19]+0.282842712474619*fUpwind[14]*alphaSurf[19]+0.282842712474619*fUpwind[13]*alphaSurf[19]+0.3162277660168379*fUpwind[3]*alphaSurf[19]+0.282842712474619*alphaSurf[16]*fUpwind[18]+0.282842712474619*alphaSurf[11]*fUpwind[18]+0.3162277660168379*alphaSurf[2]*fUpwind[18]+0.282842712474619*fUpwind[16]*alphaSurf[18]+0.282842712474619*fUpwind[11]*alphaSurf[18]+0.3162277660168379*fUpwind[2]*alphaSurf[18]+0.282842712474619*alphaSurf[15]*fUpwind[17]+0.282842712474619*alphaSurf[12]*fUpwind[17]+0.3162277660168379*alphaSurf[1]*fUpwind[17]+0.282842712474619*fUpwind[15]*alphaSurf[17]+0.282842712474619*fUpwind[12]*alphaSurf[17]+0.3162277660168379*fUpwind[1]*alphaSurf[17]+0.3162277660168379*alphaSurf[5]*fUpwind[16]+0.3162277660168379*fUpwind[5]*alphaSurf[16]+0.3162277660168379*alphaSurf[6]*fUpwind[15]+0.3162277660168379*fUpwind[6]*alphaSurf[15]+0.3162277660168379*alphaSurf[4]*fUpwind[14]+0.3162277660168379*fUpwind[4]*alphaSurf[14]+0.3162277660168379*alphaSurf[4]*fUpwind[13]+0.3162277660168379*fUpwind[4]*alphaSurf[13]+0.3162277660168379*alphaSurf[6]*fUpwind[12]+0.3162277660168379*fUpwind[6]*alphaSurf[12]+0.3162277660168379*alphaSurf[5]*fUpwind[11]+0.3162277660168379*fUpwind[5]*alphaSurf[11]+0.3162277660168379*alphaSurf[9]*fUpwind[10]+0.3162277660168379*alphaSurf[8]*fUpwind[10]+0.3162277660168379*alphaSurf[7]*fUpwind[10]+0.3535533905932737*alphaSurf[0]*fUpwind[10]+0.3162277660168379*fUpwind[9]*alphaSurf[10]+0.3162277660168379*fUpwind[8]*alphaSurf[10]+0.3162277660168379*fUpwind[7]*alphaSurf[10]+0.3535533905932737*fUpwind[0]*alphaSurf[10]+0.3535533905932737*alphaSurf[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alphaSurf[6]+0.3535533905932737*alphaSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alphaSurf[5]+0.3535533905932737*alphaSurf[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alphaSurf[4]; 
  Ghat[11] = 0.3162277660168379*alphaSurf[15]*fUpwind[19]+0.3162277660168379*fUpwind[15]*alphaSurf[19]+0.282842712474619*alphaSurf[10]*fUpwind[18]+0.282842712474619*fUpwind[10]*alphaSurf[18]+0.3162277660168379*alphaSurf[14]*fUpwind[17]+0.2258769757263128*alphaSurf[13]*fUpwind[17]+0.3535533905932737*alphaSurf[3]*fUpwind[17]+0.3162277660168379*fUpwind[14]*alphaSurf[17]+0.2258769757263128*fUpwind[13]*alphaSurf[17]+0.3535533905932737*fUpwind[3]*alphaSurf[17]+0.3535533905932737*alphaSurf[6]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaSurf[13]+0.2828427124746191*alphaSurf[4]*fUpwind[12]+0.2828427124746191*fUpwind[4]*alphaSurf[12]+0.3162277660168379*alphaSurf[8]*fUpwind[11]+0.2258769757263128*alphaSurf[7]*fUpwind[11]+0.3535533905932737*alphaSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaSurf[11]+0.2258769757263128*fUpwind[7]*alphaSurf[11]+0.3535533905932737*fUpwind[0]*alphaSurf[11]+0.3162277660168379*alphaSurf[5]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaSurf[10]+0.3535533905932737*alphaSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaSurf[7]+0.3162277660168379*alphaSurf[1]*fUpwind[4]+0.3162277660168379*fUpwind[1]*alphaSurf[4]; 
  Ghat[12] = 0.3162277660168379*alphaSurf[16]*fUpwind[19]+0.3162277660168379*fUpwind[16]*alphaSurf[19]+0.2258769757263128*alphaSurf[14]*fUpwind[18]+0.3162277660168379*alphaSurf[13]*fUpwind[18]+0.3535533905932737*alphaSurf[3]*fUpwind[18]+0.2258769757263128*fUpwind[14]*alphaSurf[18]+0.3162277660168379*fUpwind[13]*alphaSurf[18]+0.3535533905932737*fUpwind[3]*alphaSurf[18]+0.282842712474619*alphaSurf[10]*fUpwind[17]+0.282842712474619*fUpwind[10]*alphaSurf[17]+0.3535533905932737*alphaSurf[5]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaSurf[14]+0.2258769757263128*alphaSurf[8]*fUpwind[12]+0.3162277660168379*alphaSurf[7]*fUpwind[12]+0.3535533905932737*alphaSurf[0]*fUpwind[12]+0.2258769757263128*fUpwind[8]*alphaSurf[12]+0.3162277660168379*fUpwind[7]*alphaSurf[12]+0.3535533905932737*fUpwind[0]*alphaSurf[12]+0.2828427124746191*alphaSurf[4]*fUpwind[11]+0.2828427124746191*fUpwind[4]*alphaSurf[11]+0.3162277660168379*alphaSurf[6]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaSurf[10]+0.3535533905932737*alphaSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaSurf[8]+0.3162277660168379*alphaSurf[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alphaSurf[4]; 
  Ghat[13] = 0.282842712474619*alphaSurf[10]*fUpwind[19]+0.282842712474619*fUpwind[10]*alphaSurf[19]+0.3162277660168379*alphaSurf[12]*fUpwind[18]+0.3162277660168379*fUpwind[12]*alphaSurf[18]+0.3162277660168379*alphaSurf[16]*fUpwind[17]+0.2258769757263128*alphaSurf[11]*fUpwind[17]+0.3535533905932737*alphaSurf[2]*fUpwind[17]+0.3162277660168379*fUpwind[16]*alphaSurf[17]+0.2258769757263128*fUpwind[11]*alphaSurf[17]+0.3535533905932737*fUpwind[2]*alphaSurf[17]+0.2828427124746191*alphaSurf[5]*fUpwind[15]+0.2828427124746191*fUpwind[5]*alphaSurf[15]+0.3162277660168379*alphaSurf[9]*fUpwind[13]+0.2258769757263128*alphaSurf[7]*fUpwind[13]+0.3535533905932737*alphaSurf[0]*fUpwind[13]+0.3162277660168379*fUpwind[9]*alphaSurf[13]+0.2258769757263128*fUpwind[7]*alphaSurf[13]+0.3535533905932737*fUpwind[0]*alphaSurf[13]+0.3535533905932737*alphaSurf[6]*fUpwind[11]+0.3535533905932737*fUpwind[6]*alphaSurf[11]+0.3162277660168379*alphaSurf[4]*fUpwind[10]+0.3162277660168379*fUpwind[4]*alphaSurf[10]+0.3535533905932737*alphaSurf[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alphaSurf[7]+0.3162277660168379*alphaSurf[1]*fUpwind[5]+0.3162277660168379*fUpwind[1]*alphaSurf[5]; 
  Ghat[14] = 0.282842712474619*alphaSurf[10]*fUpwind[19]+0.282842712474619*fUpwind[10]*alphaSurf[19]+0.3162277660168379*alphaSurf[15]*fUpwind[18]+0.2258769757263128*alphaSurf[12]*fUpwind[18]+0.3535533905932737*alphaSurf[1]*fUpwind[18]+0.3162277660168379*fUpwind[15]*alphaSurf[18]+0.2258769757263128*fUpwind[12]*alphaSurf[18]+0.3535533905932737*fUpwind[1]*alphaSurf[18]+0.3162277660168379*alphaSurf[11]*fUpwind[17]+0.3162277660168379*fUpwind[11]*alphaSurf[17]+0.2828427124746191*alphaSurf[6]*fUpwind[16]+0.2828427124746191*fUpwind[6]*alphaSurf[16]+0.3162277660168379*alphaSurf[9]*fUpwind[14]+0.2258769757263128*alphaSurf[8]*fUpwind[14]+0.3535533905932737*alphaSurf[0]*fUpwind[14]+0.3162277660168379*fUpwind[9]*alphaSurf[14]+0.2258769757263128*fUpwind[8]*alphaSurf[14]+0.3535533905932737*fUpwind[0]*alphaSurf[14]+0.3535533905932737*alphaSurf[5]*fUpwind[12]+0.3535533905932737*fUpwind[5]*alphaSurf[12]+0.3162277660168379*alphaSurf[4]*fUpwind[10]+0.3162277660168379*fUpwind[4]*alphaSurf[10]+0.3535533905932737*alphaSurf[3]*fUpwind[8]+0.3535533905932737*fUpwind[3]*alphaSurf[8]+0.3162277660168379*alphaSurf[2]*fUpwind[6]+0.3162277660168379*fUpwind[2]*alphaSurf[6]; 
  Ghat[15] = 0.2258769757263128*alphaSurf[16]*fUpwind[19]+0.3162277660168379*alphaSurf[11]*fUpwind[19]+0.3535533905932737*alphaSurf[2]*fUpwind[19]+0.2258769757263128*fUpwind[16]*alphaSurf[19]+0.3162277660168379*fUpwind[11]*alphaSurf[19]+0.3535533905932737*fUpwind[2]*alphaSurf[19]+0.3162277660168379*alphaSurf[14]*fUpwind[18]+0.3162277660168379*fUpwind[14]*alphaSurf[18]+0.282842712474619*alphaSurf[10]*fUpwind[17]+0.282842712474619*fUpwind[10]*alphaSurf[17]+0.3535533905932737*alphaSurf[4]*fUpwind[16]+0.3535533905932737*fUpwind[4]*alphaSurf[16]+0.2258769757263128*alphaSurf[9]*fUpwind[15]+0.3162277660168379*alphaSurf[7]*fUpwind[15]+0.3535533905932737*alphaSurf[0]*fUpwind[15]+0.2258769757263128*fUpwind[9]*alphaSurf[15]+0.3162277660168379*fUpwind[7]*alphaSurf[15]+0.3535533905932737*fUpwind[0]*alphaSurf[15]+0.2828427124746191*alphaSurf[5]*fUpwind[13]+0.2828427124746191*fUpwind[5]*alphaSurf[13]+0.3162277660168379*alphaSurf[6]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaSurf[10]+0.3535533905932737*alphaSurf[1]*fUpwind[9]+0.3535533905932737*fUpwind[1]*alphaSurf[9]+0.3162277660168379*alphaSurf[3]*fUpwind[5]+0.3162277660168379*fUpwind[3]*alphaSurf[5]; 
  Ghat[16] = 0.2258769757263128*alphaSurf[15]*fUpwind[19]+0.3162277660168379*alphaSurf[12]*fUpwind[19]+0.3535533905932737*alphaSurf[1]*fUpwind[19]+0.2258769757263128*fUpwind[15]*alphaSurf[19]+0.3162277660168379*fUpwind[12]*alphaSurf[19]+0.3535533905932737*fUpwind[1]*alphaSurf[19]+0.282842712474619*alphaSurf[10]*fUpwind[18]+0.282842712474619*fUpwind[10]*alphaSurf[18]+0.3162277660168379*alphaSurf[13]*fUpwind[17]+0.3162277660168379*fUpwind[13]*alphaSurf[17]+0.2258769757263128*alphaSurf[9]*fUpwind[16]+0.3162277660168379*alphaSurf[8]*fUpwind[16]+0.3535533905932737*alphaSurf[0]*fUpwind[16]+0.2258769757263128*fUpwind[9]*alphaSurf[16]+0.3162277660168379*fUpwind[8]*alphaSurf[16]+0.3535533905932737*fUpwind[0]*alphaSurf[16]+0.3535533905932737*alphaSurf[4]*fUpwind[15]+0.3535533905932737*fUpwind[4]*alphaSurf[15]+0.2828427124746191*alphaSurf[6]*fUpwind[14]+0.2828427124746191*fUpwind[6]*alphaSurf[14]+0.3162277660168379*alphaSurf[5]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaSurf[10]+0.3535533905932737*alphaSurf[2]*fUpwind[9]+0.3535533905932737*fUpwind[2]*alphaSurf[9]+0.3162277660168379*alphaSurf[3]*fUpwind[6]+0.3162277660168379*fUpwind[3]*alphaSurf[6]; 
  Ghat[17] = 0.2529822128134704*alphaSurf[18]*fUpwind[19]+0.2828427124746191*alphaSurf[5]*fUpwind[19]+0.2529822128134704*fUpwind[18]*alphaSurf[19]+0.2828427124746191*fUpwind[5]*alphaSurf[19]+0.2828427124746191*alphaSurf[4]*fUpwind[18]+0.2828427124746191*fUpwind[4]*alphaSurf[18]+0.3162277660168379*alphaSurf[9]*fUpwind[17]+0.3162277660168379*alphaSurf[8]*fUpwind[17]+0.2258769757263128*alphaSurf[7]*fUpwind[17]+0.3535533905932737*alphaSurf[0]*fUpwind[17]+0.3162277660168379*fUpwind[9]*alphaSurf[17]+0.3162277660168379*fUpwind[8]*alphaSurf[17]+0.2258769757263128*fUpwind[7]*alphaSurf[17]+0.3535533905932737*fUpwind[0]*alphaSurf[17]+0.3162277660168379*alphaSurf[13]*fUpwind[16]+0.3162277660168379*fUpwind[13]*alphaSurf[16]+0.282842712474619*alphaSurf[10]*fUpwind[15]+0.282842712474619*fUpwind[10]*alphaSurf[15]+0.3162277660168379*alphaSurf[11]*fUpwind[14]+0.3162277660168379*fUpwind[11]*alphaSurf[14]+0.2258769757263128*alphaSurf[11]*fUpwind[13]+0.3535533905932737*alphaSurf[2]*fUpwind[13]+0.2258769757263128*fUpwind[11]*alphaSurf[13]+0.3535533905932737*fUpwind[2]*alphaSurf[13]+0.282842712474619*alphaSurf[10]*fUpwind[12]+0.282842712474619*fUpwind[10]*alphaSurf[12]+0.3535533905932737*alphaSurf[3]*fUpwind[11]+0.3535533905932737*fUpwind[3]*alphaSurf[11]+0.3162277660168379*alphaSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[1]*alphaSurf[10]+0.3535533905932737*alphaSurf[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alphaSurf[7]+0.3162277660168379*alphaSurf[4]*fUpwind[5]+0.3162277660168379*fUpwind[4]*alphaSurf[5]; 
  Ghat[18] = 0.2529822128134704*alphaSurf[17]*fUpwind[19]+0.2828427124746191*alphaSurf[6]*fUpwind[19]+0.2529822128134704*fUpwind[17]*alphaSurf[19]+0.2828427124746191*fUpwind[6]*alphaSurf[19]+0.3162277660168379*alphaSurf[9]*fUpwind[18]+0.2258769757263128*alphaSurf[8]*fUpwind[18]+0.3162277660168379*alphaSurf[7]*fUpwind[18]+0.3535533905932737*alphaSurf[0]*fUpwind[18]+0.3162277660168379*fUpwind[9]*alphaSurf[18]+0.2258769757263128*fUpwind[8]*alphaSurf[18]+0.3162277660168379*fUpwind[7]*alphaSurf[18]+0.3535533905932737*fUpwind[0]*alphaSurf[18]+0.2828427124746191*alphaSurf[4]*fUpwind[17]+0.2828427124746191*fUpwind[4]*alphaSurf[17]+0.282842712474619*alphaSurf[10]*fUpwind[16]+0.282842712474619*fUpwind[10]*alphaSurf[16]+0.3162277660168379*alphaSurf[14]*fUpwind[15]+0.3162277660168379*fUpwind[14]*alphaSurf[15]+0.2258769757263128*alphaSurf[12]*fUpwind[14]+0.3535533905932737*alphaSurf[1]*fUpwind[14]+0.2258769757263128*fUpwind[12]*alphaSurf[14]+0.3535533905932737*fUpwind[1]*alphaSurf[14]+0.3162277660168379*alphaSurf[12]*fUpwind[13]+0.3162277660168379*fUpwind[12]*alphaSurf[13]+0.3535533905932737*alphaSurf[3]*fUpwind[12]+0.3535533905932737*fUpwind[3]*alphaSurf[12]+0.282842712474619*alphaSurf[10]*fUpwind[11]+0.282842712474619*fUpwind[10]*alphaSurf[11]+0.3162277660168379*alphaSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[2]*alphaSurf[10]+0.3535533905932737*alphaSurf[5]*fUpwind[8]+0.3535533905932737*fUpwind[5]*alphaSurf[8]+0.3162277660168379*alphaSurf[4]*fUpwind[6]+0.3162277660168379*fUpwind[4]*alphaSurf[6]; 
  Ghat[19] = 0.2258769757263128*alphaSurf[9]*fUpwind[19]+0.3162277660168379*alphaSurf[8]*fUpwind[19]+0.3162277660168379*alphaSurf[7]*fUpwind[19]+0.3535533905932737*alphaSurf[0]*fUpwind[19]+0.2258769757263128*fUpwind[9]*alphaSurf[19]+0.3162277660168379*fUpwind[8]*alphaSurf[19]+0.3162277660168379*fUpwind[7]*alphaSurf[19]+0.3535533905932737*fUpwind[0]*alphaSurf[19]+0.2529822128134704*alphaSurf[17]*fUpwind[18]+0.2828427124746191*alphaSurf[6]*fUpwind[18]+0.2529822128134704*fUpwind[17]*alphaSurf[18]+0.2828427124746191*fUpwind[6]*alphaSurf[18]+0.2828427124746191*alphaSurf[5]*fUpwind[17]+0.2828427124746191*fUpwind[5]*alphaSurf[17]+0.2258769757263128*alphaSurf[15]*fUpwind[16]+0.3162277660168379*alphaSurf[12]*fUpwind[16]+0.3535533905932737*alphaSurf[1]*fUpwind[16]+0.2258769757263128*fUpwind[15]*alphaSurf[16]+0.3162277660168379*fUpwind[12]*alphaSurf[16]+0.3535533905932737*fUpwind[1]*alphaSurf[16]+0.3162277660168379*alphaSurf[11]*fUpwind[15]+0.3535533905932737*alphaSurf[2]*fUpwind[15]+0.3162277660168379*fUpwind[11]*alphaSurf[15]+0.3535533905932737*fUpwind[2]*alphaSurf[15]+0.282842712474619*alphaSurf[10]*fUpwind[14]+0.282842712474619*fUpwind[10]*alphaSurf[14]+0.282842712474619*alphaSurf[10]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaSurf[13]+0.3162277660168379*alphaSurf[3]*fUpwind[10]+0.3162277660168379*fUpwind[3]*alphaSurf[10]+0.3535533905932737*alphaSurf[4]*fUpwind[9]+0.3535533905932737*fUpwind[4]*alphaSurf[9]+0.3162277660168379*alphaSurf[5]*fUpwind[6]+0.3162277660168379*fUpwind[5]*alphaSurf[6]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += -0.7071067811865475*Ghat[3]*dv1par; 
  out[4] += -1.224744871391589*Ghat[0]*dv1par; 
  out[5] += -0.7071067811865475*Ghat[4]*dv1par; 
  out[6] += -0.7071067811865475*Ghat[5]*dv1par; 
  out[7] += -0.7071067811865475*Ghat[6]*dv1par; 
  out[8] += -1.224744871391589*Ghat[1]*dv1par; 
  out[9] += -1.224744871391589*Ghat[2]*dv1par; 
  out[10] += -1.224744871391589*Ghat[3]*dv1par; 
  out[11] += -0.7071067811865475*Ghat[7]*dv1par; 
  out[12] += -0.7071067811865475*Ghat[8]*dv1par; 
  out[13] += -0.7071067811865475*Ghat[9]*dv1par; 
  out[14] += -1.58113883008419*Ghat[0]*dv1par; 
  out[15] += -0.7071067811865475*Ghat[10]*dv1par; 
  out[16] += -1.224744871391589*Ghat[4]*dv1par; 
  out[17] += -1.224744871391589*Ghat[5]*dv1par; 
  out[18] += -1.224744871391589*Ghat[6]*dv1par; 
  out[19] += -0.7071067811865475*Ghat[11]*dv1par; 
  out[20] += -0.7071067811865475*Ghat[12]*dv1par; 
  out[21] += -0.7071067811865475*Ghat[13]*dv1par; 
  out[22] += -0.7071067811865475*Ghat[14]*dv1par; 
  out[23] += -0.7071067811865475*Ghat[15]*dv1par; 
  out[24] += -0.7071067811865475*Ghat[16]*dv1par; 
  out[25] += -1.224744871391589*Ghat[7]*dv1par; 
  out[26] += -1.224744871391589*Ghat[8]*dv1par; 
  out[27] += -1.224744871391589*Ghat[9]*dv1par; 
  out[28] += -1.58113883008419*Ghat[1]*dv1par; 
  out[29] += -1.58113883008419*Ghat[2]*dv1par; 
  out[30] += -1.58113883008419*Ghat[3]*dv1par; 
  out[31] += -1.224744871391589*Ghat[10]*dv1par; 
  out[32] += -0.7071067811865475*Ghat[17]*dv1par; 
  out[33] += -0.7071067811865475*Ghat[18]*dv1par; 
  out[34] += -0.7071067811865475*Ghat[19]*dv1par; 
  out[35] += -1.224744871391589*Ghat[11]*dv1par; 
  out[36] += -1.224744871391589*Ghat[12]*dv1par; 
  out[37] += -1.224744871391589*Ghat[13]*dv1par; 
  out[38] += -1.224744871391589*Ghat[14]*dv1par; 
  out[39] += -1.224744871391589*Ghat[15]*dv1par; 
  out[40] += -1.224744871391589*Ghat[16]*dv1par; 
  out[41] += -1.58113883008419*Ghat[4]*dv1par; 
  out[42] += -1.58113883008419*Ghat[5]*dv1par; 
  out[43] += -1.58113883008419*Ghat[6]*dv1par; 
  out[44] += -1.224744871391589*Ghat[17]*dv1par; 
  out[45] += -1.224744871391589*Ghat[18]*dv1par; 
  out[46] += -1.224744871391589*Ghat[19]*dv1par; 
  out[47] += -1.58113883008419*Ghat[10]*dv1par; 

  } else { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force[7]; 
  alphaSurf[8] = (-1.0*bb_grad_u[8]*wvpar)+0.5*bb_grad_u[8]*dvpar+p_force[8]; 
  alphaSurf[9] = (-1.0*bb_grad_u[9]*wvpar)+0.5*bb_grad_u[9]*dvpar+p_force[9]; 
  alphaSurf[10] = (-1.0*bb_grad_u[10]*wvpar)+0.5*bb_grad_u[10]*dvpar+p_force[10]; 
  alphaSurf[11] = (-1.0*bb_grad_u[11]*wvpar)+0.5*bb_grad_u[11]*dvpar+p_force[11]; 
  alphaSurf[12] = (-1.0*bb_grad_u[12]*wvpar)+0.5*bb_grad_u[12]*dvpar+p_force[12]; 
  alphaSurf[13] = (-1.0*bb_grad_u[13]*wvpar)+0.5*bb_grad_u[13]*dvpar+p_force[13]; 
  alphaSurf[14] = (-1.0*bb_grad_u[14]*wvpar)+0.5*bb_grad_u[14]*dvpar+p_force[14]; 
  alphaSurf[15] = (-1.0*bb_grad_u[15]*wvpar)+0.5*bb_grad_u[15]*dvpar+p_force[15]; 
  alphaSurf[16] = (-1.0*bb_grad_u[16]*wvpar)+0.5*bb_grad_u[16]*dvpar+p_force[16]; 
  alphaSurf[17] = (-1.0*bb_grad_u[17]*wvpar)+0.5*bb_grad_u[17]*dvpar+p_force[17]; 
  alphaSurf[18] = (-1.0*bb_grad_u[18]*wvpar)+0.5*bb_grad_u[18]*dvpar+p_force[18]; 
  alphaSurf[19] = (-1.0*bb_grad_u[19]*wvpar)+0.5*bb_grad_u[19]*dvpar+p_force[19]; 

  if (0.5692099788303082*(alphaSurf[19]+alphaSurf[18]+alphaSurf[17])-0.4242640687119281*(alphaSurf[16]+alphaSurf[15])-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*(alphaSurf[6]+alphaSurf[5]+alphaSurf[4])-0.4743416490252568*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if ((-0.711512473537885*alphaSurf[19])+0.5303300858899104*(alphaSurf[16]+alphaSurf[15])-0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[4]-0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (0.5692099788303082*alphaSurf[19]-0.5692099788303082*(alphaSurf[18]+alphaSurf[17])-0.4242640687119281*(alphaSurf[16]+alphaSurf[15])+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*(alphaSurf[6]+alphaSurf[5])+0.6363961030678926*alphaSurf[4]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if ((-0.711512473537885*alphaSurf[18])-0.4242640687119281*alphaSurf[15]+0.5303300858899104*alphaSurf[14]-0.4242640687119285*alphaSurf[13]+0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]+0.6363961030678926*alphaSurf[5]-0.4743416490252568*(alphaSurf[3]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (0.5303300858899104*(alphaSurf[15]+alphaSurf[12])-0.3952847075210473*(alphaSurf[9]+alphaSurf[8])+0.3162277660168379*alphaSurf[7]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if (0.711512473537885*alphaSurf[18]-0.4242640687119281*alphaSurf[15]-0.5303300858899104*alphaSurf[14]+0.4242640687119285*alphaSurf[13]+0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]-0.6363961030678926*alphaSurf[5]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.5692099788303082*alphaSurf[19])+0.5692099788303082*alphaSurf[18]-0.5692099788303082*alphaSurf[17]+0.4242640687119281*alphaSurf[16]-0.4242640687119281*alphaSurf[15]-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[6]+0.6363961030678926*alphaSurf[5]-0.6363961030678926*alphaSurf[4]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*alphaSurf[2]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if (0.711512473537885*alphaSurf[19]-0.5303300858899104*alphaSurf[16]+0.5303300858899104*alphaSurf[15]-0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[4]+0.4743416490252568*alphaSurf[2]-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if ((-0.5692099788303082*(alphaSurf[19]+alphaSurf[18]))+0.5692099788303082*alphaSurf[17]+0.4242640687119281*alphaSurf[16]-0.4242640687119281*alphaSurf[15]+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])-0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[6]-0.6363961030678926*(alphaSurf[5]+alphaSurf[4])+0.4743416490252568*(alphaSurf[3]+alphaSurf[2])-0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if ((-0.711512473537885*alphaSurf[17])-0.4242640687119281*alphaSurf[16]-0.4242640687119285*alphaSurf[14]+0.5303300858899104*(alphaSurf[13]+alphaSurf[11])+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]+0.6363961030678926*alphaSurf[6]-0.4743416490252568*(alphaSurf[3]+alphaSurf[2])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fSkin); 
  } 
  if (0.5303300858899104*(alphaSurf[16]+alphaSurf[11])-0.3952847075210473*alphaSurf[9]+0.3162277660168379*alphaSurf[8]-0.3952847075210473*alphaSurf[7]-0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fSkin); 
  } 
  if (0.711512473537885*alphaSurf[17]-0.4242640687119281*alphaSurf[16]+0.4242640687119285*alphaSurf[14]-0.5303300858899104*alphaSurf[13]+0.5303300858899104*alphaSurf[11]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]-0.6363961030678926*alphaSurf[6]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (0.5303300858899104*(alphaSurf[14]+alphaSurf[13])+0.3162277660168379*alphaSurf[9]-0.3952847075210473*(alphaSurf[8]+alphaSurf[7])-0.4743416490252568*alphaSurf[3]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fSkin); 
  } 
  if (0.3535533905932737*alphaSurf[0]-0.3952847075210473*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7]) > 0) { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fSkin); 
  } 
  if ((-0.5303300858899104*(alphaSurf[14]+alphaSurf[13]))+0.3162277660168379*alphaSurf[9]-0.3952847075210473*(alphaSurf[8]+alphaSurf[7])+0.4743416490252568*alphaSurf[3]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (0.711512473537885*alphaSurf[17]+0.4242640687119281*alphaSurf[16]-0.4242640687119285*alphaSurf[14]+0.5303300858899104*alphaSurf[13]-0.5303300858899104*alphaSurf[11]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]-0.6363961030678926*alphaSurf[6]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fSkin); 
  } 
  if ((-0.5303300858899104*(alphaSurf[16]+alphaSurf[11]))-0.3952847075210473*alphaSurf[9]+0.3162277660168379*alphaSurf[8]-0.3952847075210473*alphaSurf[7]+0.4743416490252568*alphaSurf[2]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fSkin); 
  } 
  if ((-0.711512473537885*alphaSurf[17])+0.4242640687119281*alphaSurf[16]+0.4242640687119285*alphaSurf[14]-0.5303300858899104*(alphaSurf[13]+alphaSurf[11])+0.3162277660168379*(alphaSurf[9]+alphaSurf[8])-0.3952847075210473*alphaSurf[7]+0.6363961030678926*alphaSurf[6]+0.4743416490252568*(alphaSurf[3]+alphaSurf[2])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fSkin); 
  } 
  if ((-0.5692099788303082*(alphaSurf[19]+alphaSurf[18]))+0.5692099788303082*alphaSurf[17]-0.4242640687119281*alphaSurf[16]+0.4242640687119281*alphaSurf[15]-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[6]-0.6363961030678926*(alphaSurf[5]+alphaSurf[4])-0.4743416490252568*(alphaSurf[3]+alphaSurf[2])+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fSkin); 
  } 
  if (0.711512473537885*alphaSurf[19]+0.5303300858899104*alphaSurf[16]-0.5303300858899104*alphaSurf[15]+0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[4]-0.4743416490252568*alphaSurf[2]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fSkin); 
  } 
  if ((-0.5692099788303082*alphaSurf[19])+0.5692099788303082*alphaSurf[18]-0.5692099788303082*alphaSurf[17]-0.4242640687119281*alphaSurf[16]+0.4242640687119281*alphaSurf[15]+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]-0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*alphaSurf[6]+0.6363961030678926*alphaSurf[5]-0.6363961030678926*alphaSurf[4]+0.4743416490252568*alphaSurf[3]-0.4743416490252568*alphaSurf[2]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fSkin); 
  } 
  if (0.711512473537885*alphaSurf[18]+0.4242640687119281*alphaSurf[15]+0.5303300858899104*alphaSurf[14]-0.4242640687119285*alphaSurf[13]-0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]-0.6363961030678926*alphaSurf[5]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fSkin); 
  } 
  if ((-0.5303300858899104*(alphaSurf[15]+alphaSurf[12]))-0.3952847075210473*(alphaSurf[9]+alphaSurf[8])+0.3162277660168379*alphaSurf[7]+0.4743416490252568*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fEdge); 
  } else { 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fSkin); 
  } 
  if ((-0.711512473537885*alphaSurf[18])+0.4242640687119281*alphaSurf[15]-0.5303300858899104*alphaSurf[14]+0.4242640687119285*alphaSurf[13]-0.5303300858899104*alphaSurf[12]+0.3162277660168379*alphaSurf[9]-0.3952847075210473*alphaSurf[8]+0.3162277660168379*alphaSurf[7]+0.6363961030678926*alphaSurf[5]+0.4743416490252568*(alphaSurf[3]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if (0.5692099788303082*alphaSurf[19]-0.5692099788303082*(alphaSurf[18]+alphaSurf[17])+0.4242640687119281*(alphaSurf[16]+alphaSurf[15])-0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])-0.6363961030678926*(alphaSurf[6]+alphaSurf[5])+0.6363961030678926*alphaSurf[4]-0.4743416490252568*alphaSurf[3]+0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fSkin); 
  } 
  if ((-0.711512473537885*alphaSurf[19])-0.5303300858899104*(alphaSurf[16]+alphaSurf[15])+0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]-0.3952847075210473*alphaSurf[9]+0.3162277660168379*(alphaSurf[8]+alphaSurf[7])+0.6363961030678926*alphaSurf[4]+0.4743416490252568*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fEdge); 
  } else { 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fSkin); 
  } 
  if (0.5692099788303082*(alphaSurf[19]+alphaSurf[18]+alphaSurf[17])+0.4242640687119281*(alphaSurf[16]+alphaSurf[15])+0.4242640687119285*(alphaSurf[14]+alphaSurf[13])+0.4242640687119281*alphaSurf[12]+0.4242640687119285*alphaSurf[11]+0.853814968245462*alphaSurf[10]+0.3162277660168379*(alphaSurf[9]+alphaSurf[8]+alphaSurf[7])+0.6363961030678926*(alphaSurf[6]+alphaSurf[5]+alphaSurf[4])+0.4743416490252568*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaSurf[19]*fUpwind[19]+0.3535533905932737*alphaSurf[18]*fUpwind[18]+0.3535533905932737*alphaSurf[17]*fUpwind[17]+0.3535533905932737*alphaSurf[16]*fUpwind[16]+0.3535533905932737*alphaSurf[15]*fUpwind[15]+0.3535533905932737*alphaSurf[14]*fUpwind[14]+0.3535533905932737*alphaSurf[13]*fUpwind[13]+0.3535533905932737*alphaSurf[12]*fUpwind[12]+0.3535533905932737*alphaSurf[11]*fUpwind[11]+0.3535533905932737*alphaSurf[10]*fUpwind[10]+0.3535533905932737*alphaSurf[9]*fUpwind[9]+0.3535533905932737*alphaSurf[8]*fUpwind[8]+0.3535533905932737*alphaSurf[7]*fUpwind[7]+0.3535533905932737*alphaSurf[6]*fUpwind[6]+0.3535533905932737*alphaSurf[5]*fUpwind[5]+0.3535533905932737*alphaSurf[4]*fUpwind[4]+0.3535533905932737*alphaSurf[3]*fUpwind[3]+0.3535533905932737*alphaSurf[2]*fUpwind[2]+0.3535533905932737*alphaSurf[1]*fUpwind[1]+0.3535533905932737*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaSurf[16]*fUpwind[19]+0.3535533905932737*fUpwind[16]*alphaSurf[19]+0.3535533905932737*alphaSurf[14]*fUpwind[18]+0.3535533905932737*fUpwind[14]*alphaSurf[18]+0.3162277660168379*alphaSurf[10]*fUpwind[17]+0.3162277660168379*fUpwind[10]*alphaSurf[17]+0.3535533905932737*alphaSurf[9]*fUpwind[15]+0.3535533905932737*fUpwind[9]*alphaSurf[15]+0.3162277660168379*alphaSurf[5]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alphaSurf[13]+0.3535533905932737*alphaSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaSurf[12]+0.3162277660168379*alphaSurf[4]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaSurf[11]+0.3535533905932737*alphaSurf[6]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaSurf[10]+0.3162277660168379*alphaSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaSurf[7]+0.3535533905932737*alphaSurf[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.3535533905932737*alphaSurf[15]*fUpwind[19]+0.3535533905932737*fUpwind[15]*alphaSurf[19]+0.3162277660168379*alphaSurf[10]*fUpwind[18]+0.3162277660168379*fUpwind[10]*alphaSurf[18]+0.3535533905932737*alphaSurf[13]*fUpwind[17]+0.3535533905932737*fUpwind[13]*alphaSurf[17]+0.3535533905932737*alphaSurf[9]*fUpwind[16]+0.3535533905932737*fUpwind[9]*alphaSurf[16]+0.3162277660168379*alphaSurf[6]*fUpwind[14]+0.3162277660168379*fUpwind[6]*alphaSurf[14]+0.3162277660168379*alphaSurf[4]*fUpwind[12]+0.3162277660168379*fUpwind[4]*alphaSurf[12]+0.3535533905932737*alphaSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaSurf[11]+0.3535533905932737*alphaSurf[5]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaSurf[10]+0.3162277660168379*alphaSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaSurf[8]+0.3535533905932737*alphaSurf[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaSurf[2]; 
  Ghat[3] = 0.3162277660168379*alphaSurf[10]*fUpwind[19]+0.3162277660168379*fUpwind[10]*alphaSurf[19]+0.3535533905932737*alphaSurf[12]*fUpwind[18]+0.3535533905932737*fUpwind[12]*alphaSurf[18]+0.3535533905932737*alphaSurf[11]*fUpwind[17]+0.3535533905932737*fUpwind[11]*alphaSurf[17]+0.3162277660168379*alphaSurf[6]*fUpwind[16]+0.3162277660168379*fUpwind[6]*alphaSurf[16]+0.3162277660168379*alphaSurf[5]*fUpwind[15]+0.3162277660168379*fUpwind[5]*alphaSurf[15]+0.3535533905932737*alphaSurf[8]*fUpwind[14]+0.3535533905932737*fUpwind[8]*alphaSurf[14]+0.3535533905932737*alphaSurf[7]*fUpwind[13]+0.3535533905932737*fUpwind[7]*alphaSurf[13]+0.3535533905932737*alphaSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[4]*alphaSurf[10]+0.3162277660168379*alphaSurf[3]*fUpwind[9]+0.3162277660168379*fUpwind[3]*alphaSurf[9]+0.3535533905932737*alphaSurf[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alphaSurf[5]+0.3535533905932737*alphaSurf[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alphaSurf[3]; 
  Ghat[4] = 0.3535533905932737*alphaSurf[9]*fUpwind[19]+0.3535533905932737*fUpwind[9]*alphaSurf[19]+0.2828427124746191*alphaSurf[17]*fUpwind[18]+0.3162277660168379*alphaSurf[6]*fUpwind[18]+0.2828427124746191*fUpwind[17]*alphaSurf[18]+0.3162277660168379*fUpwind[6]*alphaSurf[18]+0.3162277660168379*alphaSurf[5]*fUpwind[17]+0.3162277660168379*fUpwind[5]*alphaSurf[17]+0.3535533905932737*alphaSurf[15]*fUpwind[16]+0.3535533905932737*fUpwind[15]*alphaSurf[16]+0.3162277660168379*alphaSurf[10]*fUpwind[14]+0.3162277660168379*fUpwind[10]*alphaSurf[14]+0.3162277660168379*alphaSurf[10]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaSurf[13]+0.2828427124746191*alphaSurf[11]*fUpwind[12]+0.3162277660168379*alphaSurf[2]*fUpwind[12]+0.2828427124746191*fUpwind[11]*alphaSurf[12]+0.3162277660168379*fUpwind[2]*alphaSurf[12]+0.3162277660168379*alphaSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[1]*alphaSurf[11]+0.3535533905932737*alphaSurf[3]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaSurf[10]+0.3162277660168379*alphaSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaSurf[8]+0.3162277660168379*alphaSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaSurf[7]+0.3535533905932737*alphaSurf[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaSurf[4]+0.3535533905932737*alphaSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaSurf[2]; 
  Ghat[5] = 0.2828427124746191*alphaSurf[17]*fUpwind[19]+0.3162277660168379*alphaSurf[6]*fUpwind[19]+0.2828427124746191*fUpwind[17]*alphaSurf[19]+0.3162277660168379*fUpwind[6]*alphaSurf[19]+0.3535533905932737*alphaSurf[8]*fUpwind[18]+0.3535533905932737*fUpwind[8]*alphaSurf[18]+0.3162277660168379*alphaSurf[4]*fUpwind[17]+0.3162277660168379*fUpwind[4]*alphaSurf[17]+0.3162277660168379*alphaSurf[10]*fUpwind[16]+0.3162277660168379*fUpwind[10]*alphaSurf[16]+0.2828427124746191*alphaSurf[13]*fUpwind[15]+0.3162277660168379*alphaSurf[3]*fUpwind[15]+0.2828427124746191*fUpwind[13]*alphaSurf[15]+0.3162277660168379*fUpwind[3]*alphaSurf[15]+0.3535533905932737*alphaSurf[12]*fUpwind[14]+0.3535533905932737*fUpwind[12]*alphaSurf[14]+0.3162277660168379*alphaSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[1]*alphaSurf[13]+0.3162277660168379*alphaSurf[10]*fUpwind[11]+0.3162277660168379*fUpwind[10]*alphaSurf[11]+0.3535533905932737*alphaSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[2]*alphaSurf[10]+0.3162277660168379*alphaSurf[5]*fUpwind[9]+0.3162277660168379*fUpwind[5]*alphaSurf[9]+0.3162277660168379*alphaSurf[5]*fUpwind[7]+0.3162277660168379*fUpwind[5]*alphaSurf[7]+0.3535533905932737*alphaSurf[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alphaSurf[5]+0.3535533905932737*alphaSurf[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alphaSurf[3]; 
  Ghat[6] = 0.2828427124746191*alphaSurf[18]*fUpwind[19]+0.3162277660168379*alphaSurf[5]*fUpwind[19]+0.2828427124746191*fUpwind[18]*alphaSurf[19]+0.3162277660168379*fUpwind[5]*alphaSurf[19]+0.3162277660168379*alphaSurf[4]*fUpwind[18]+0.3162277660168379*fUpwind[4]*alphaSurf[18]+0.3535533905932737*alphaSurf[7]*fUpwind[17]+0.3535533905932737*fUpwind[7]*alphaSurf[17]+0.2828427124746191*alphaSurf[14]*fUpwind[16]+0.3162277660168379*alphaSurf[3]*fUpwind[16]+0.2828427124746191*fUpwind[14]*alphaSurf[16]+0.3162277660168379*fUpwind[3]*alphaSurf[16]+0.3162277660168379*alphaSurf[10]*fUpwind[15]+0.3162277660168379*fUpwind[10]*alphaSurf[15]+0.3162277660168379*alphaSurf[2]*fUpwind[14]+0.3162277660168379*fUpwind[2]*alphaSurf[14]+0.3535533905932737*alphaSurf[11]*fUpwind[13]+0.3535533905932737*fUpwind[11]*alphaSurf[13]+0.3162277660168379*alphaSurf[10]*fUpwind[12]+0.3162277660168379*fUpwind[10]*alphaSurf[12]+0.3535533905932737*alphaSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[1]*alphaSurf[10]+0.3162277660168379*alphaSurf[6]*fUpwind[9]+0.3162277660168379*fUpwind[6]*alphaSurf[9]+0.3162277660168379*alphaSurf[6]*fUpwind[8]+0.3162277660168379*fUpwind[6]*alphaSurf[8]+0.3535533905932737*alphaSurf[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alphaSurf[6]+0.3535533905932737*alphaSurf[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alphaSurf[3]; 
  Ghat[7] = 0.3162277660168379*alphaSurf[19]*fUpwind[19]+0.3162277660168379*alphaSurf[18]*fUpwind[18]+0.2258769757263128*alphaSurf[17]*fUpwind[17]+0.3535533905932737*alphaSurf[6]*fUpwind[17]+0.3535533905932737*fUpwind[6]*alphaSurf[17]+0.3162277660168379*alphaSurf[15]*fUpwind[15]+0.2258769757263128*alphaSurf[13]*fUpwind[13]+0.3535533905932737*alphaSurf[3]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaSurf[13]+0.3162277660168379*alphaSurf[12]*fUpwind[12]+0.2258769757263128*alphaSurf[11]*fUpwind[11]+0.3535533905932737*alphaSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaSurf[11]+0.3162277660168379*alphaSurf[10]*fUpwind[10]+0.2258769757263128*alphaSurf[7]*fUpwind[7]+0.3535533905932737*alphaSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaSurf[7]+0.3162277660168379*alphaSurf[5]*fUpwind[5]+0.3162277660168379*alphaSurf[4]*fUpwind[4]+0.3162277660168379*alphaSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.3162277660168379*alphaSurf[19]*fUpwind[19]+0.2258769757263128*alphaSurf[18]*fUpwind[18]+0.3535533905932737*alphaSurf[5]*fUpwind[18]+0.3535533905932737*fUpwind[5]*alphaSurf[18]+0.3162277660168379*alphaSurf[17]*fUpwind[17]+0.3162277660168379*alphaSurf[16]*fUpwind[16]+0.2258769757263128*alphaSurf[14]*fUpwind[14]+0.3535533905932737*alphaSurf[3]*fUpwind[14]+0.3535533905932737*fUpwind[3]*alphaSurf[14]+0.2258769757263128*alphaSurf[12]*fUpwind[12]+0.3535533905932737*alphaSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaSurf[12]+0.3162277660168379*alphaSurf[11]*fUpwind[11]+0.3162277660168379*alphaSurf[10]*fUpwind[10]+0.2258769757263128*alphaSurf[8]*fUpwind[8]+0.3535533905932737*alphaSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaSurf[8]+0.3162277660168379*alphaSurf[6]*fUpwind[6]+0.3162277660168379*alphaSurf[4]*fUpwind[4]+0.3162277660168379*alphaSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.2258769757263128*alphaSurf[19]*fUpwind[19]+0.3535533905932737*alphaSurf[4]*fUpwind[19]+0.3535533905932737*fUpwind[4]*alphaSurf[19]+0.3162277660168379*alphaSurf[18]*fUpwind[18]+0.3162277660168379*alphaSurf[17]*fUpwind[17]+0.2258769757263128*alphaSurf[16]*fUpwind[16]+0.3535533905932737*alphaSurf[2]*fUpwind[16]+0.3535533905932737*fUpwind[2]*alphaSurf[16]+0.2258769757263128*alphaSurf[15]*fUpwind[15]+0.3535533905932737*alphaSurf[1]*fUpwind[15]+0.3535533905932737*fUpwind[1]*alphaSurf[15]+0.3162277660168379*alphaSurf[14]*fUpwind[14]+0.3162277660168379*alphaSurf[13]*fUpwind[13]+0.3162277660168379*alphaSurf[10]*fUpwind[10]+0.2258769757263128*alphaSurf[9]*fUpwind[9]+0.3535533905932737*alphaSurf[0]*fUpwind[9]+0.3535533905932737*fUpwind[0]*alphaSurf[9]+0.3162277660168379*alphaSurf[6]*fUpwind[6]+0.3162277660168379*alphaSurf[5]*fUpwind[5]+0.3162277660168379*alphaSurf[3]*fUpwind[3]; 
  Ghat[10] = 0.282842712474619*alphaSurf[14]*fUpwind[19]+0.282842712474619*alphaSurf[13]*fUpwind[19]+0.3162277660168379*alphaSurf[3]*fUpwind[19]+0.282842712474619*fUpwind[14]*alphaSurf[19]+0.282842712474619*fUpwind[13]*alphaSurf[19]+0.3162277660168379*fUpwind[3]*alphaSurf[19]+0.282842712474619*alphaSurf[16]*fUpwind[18]+0.282842712474619*alphaSurf[11]*fUpwind[18]+0.3162277660168379*alphaSurf[2]*fUpwind[18]+0.282842712474619*fUpwind[16]*alphaSurf[18]+0.282842712474619*fUpwind[11]*alphaSurf[18]+0.3162277660168379*fUpwind[2]*alphaSurf[18]+0.282842712474619*alphaSurf[15]*fUpwind[17]+0.282842712474619*alphaSurf[12]*fUpwind[17]+0.3162277660168379*alphaSurf[1]*fUpwind[17]+0.282842712474619*fUpwind[15]*alphaSurf[17]+0.282842712474619*fUpwind[12]*alphaSurf[17]+0.3162277660168379*fUpwind[1]*alphaSurf[17]+0.3162277660168379*alphaSurf[5]*fUpwind[16]+0.3162277660168379*fUpwind[5]*alphaSurf[16]+0.3162277660168379*alphaSurf[6]*fUpwind[15]+0.3162277660168379*fUpwind[6]*alphaSurf[15]+0.3162277660168379*alphaSurf[4]*fUpwind[14]+0.3162277660168379*fUpwind[4]*alphaSurf[14]+0.3162277660168379*alphaSurf[4]*fUpwind[13]+0.3162277660168379*fUpwind[4]*alphaSurf[13]+0.3162277660168379*alphaSurf[6]*fUpwind[12]+0.3162277660168379*fUpwind[6]*alphaSurf[12]+0.3162277660168379*alphaSurf[5]*fUpwind[11]+0.3162277660168379*fUpwind[5]*alphaSurf[11]+0.3162277660168379*alphaSurf[9]*fUpwind[10]+0.3162277660168379*alphaSurf[8]*fUpwind[10]+0.3162277660168379*alphaSurf[7]*fUpwind[10]+0.3535533905932737*alphaSurf[0]*fUpwind[10]+0.3162277660168379*fUpwind[9]*alphaSurf[10]+0.3162277660168379*fUpwind[8]*alphaSurf[10]+0.3162277660168379*fUpwind[7]*alphaSurf[10]+0.3535533905932737*fUpwind[0]*alphaSurf[10]+0.3535533905932737*alphaSurf[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alphaSurf[6]+0.3535533905932737*alphaSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alphaSurf[5]+0.3535533905932737*alphaSurf[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alphaSurf[4]; 
  Ghat[11] = 0.3162277660168379*alphaSurf[15]*fUpwind[19]+0.3162277660168379*fUpwind[15]*alphaSurf[19]+0.282842712474619*alphaSurf[10]*fUpwind[18]+0.282842712474619*fUpwind[10]*alphaSurf[18]+0.3162277660168379*alphaSurf[14]*fUpwind[17]+0.2258769757263128*alphaSurf[13]*fUpwind[17]+0.3535533905932737*alphaSurf[3]*fUpwind[17]+0.3162277660168379*fUpwind[14]*alphaSurf[17]+0.2258769757263128*fUpwind[13]*alphaSurf[17]+0.3535533905932737*fUpwind[3]*alphaSurf[17]+0.3535533905932737*alphaSurf[6]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaSurf[13]+0.2828427124746191*alphaSurf[4]*fUpwind[12]+0.2828427124746191*fUpwind[4]*alphaSurf[12]+0.3162277660168379*alphaSurf[8]*fUpwind[11]+0.2258769757263128*alphaSurf[7]*fUpwind[11]+0.3535533905932737*alphaSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaSurf[11]+0.2258769757263128*fUpwind[7]*alphaSurf[11]+0.3535533905932737*fUpwind[0]*alphaSurf[11]+0.3162277660168379*alphaSurf[5]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaSurf[10]+0.3535533905932737*alphaSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaSurf[7]+0.3162277660168379*alphaSurf[1]*fUpwind[4]+0.3162277660168379*fUpwind[1]*alphaSurf[4]; 
  Ghat[12] = 0.3162277660168379*alphaSurf[16]*fUpwind[19]+0.3162277660168379*fUpwind[16]*alphaSurf[19]+0.2258769757263128*alphaSurf[14]*fUpwind[18]+0.3162277660168379*alphaSurf[13]*fUpwind[18]+0.3535533905932737*alphaSurf[3]*fUpwind[18]+0.2258769757263128*fUpwind[14]*alphaSurf[18]+0.3162277660168379*fUpwind[13]*alphaSurf[18]+0.3535533905932737*fUpwind[3]*alphaSurf[18]+0.282842712474619*alphaSurf[10]*fUpwind[17]+0.282842712474619*fUpwind[10]*alphaSurf[17]+0.3535533905932737*alphaSurf[5]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaSurf[14]+0.2258769757263128*alphaSurf[8]*fUpwind[12]+0.3162277660168379*alphaSurf[7]*fUpwind[12]+0.3535533905932737*alphaSurf[0]*fUpwind[12]+0.2258769757263128*fUpwind[8]*alphaSurf[12]+0.3162277660168379*fUpwind[7]*alphaSurf[12]+0.3535533905932737*fUpwind[0]*alphaSurf[12]+0.2828427124746191*alphaSurf[4]*fUpwind[11]+0.2828427124746191*fUpwind[4]*alphaSurf[11]+0.3162277660168379*alphaSurf[6]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaSurf[10]+0.3535533905932737*alphaSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaSurf[8]+0.3162277660168379*alphaSurf[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alphaSurf[4]; 
  Ghat[13] = 0.282842712474619*alphaSurf[10]*fUpwind[19]+0.282842712474619*fUpwind[10]*alphaSurf[19]+0.3162277660168379*alphaSurf[12]*fUpwind[18]+0.3162277660168379*fUpwind[12]*alphaSurf[18]+0.3162277660168379*alphaSurf[16]*fUpwind[17]+0.2258769757263128*alphaSurf[11]*fUpwind[17]+0.3535533905932737*alphaSurf[2]*fUpwind[17]+0.3162277660168379*fUpwind[16]*alphaSurf[17]+0.2258769757263128*fUpwind[11]*alphaSurf[17]+0.3535533905932737*fUpwind[2]*alphaSurf[17]+0.2828427124746191*alphaSurf[5]*fUpwind[15]+0.2828427124746191*fUpwind[5]*alphaSurf[15]+0.3162277660168379*alphaSurf[9]*fUpwind[13]+0.2258769757263128*alphaSurf[7]*fUpwind[13]+0.3535533905932737*alphaSurf[0]*fUpwind[13]+0.3162277660168379*fUpwind[9]*alphaSurf[13]+0.2258769757263128*fUpwind[7]*alphaSurf[13]+0.3535533905932737*fUpwind[0]*alphaSurf[13]+0.3535533905932737*alphaSurf[6]*fUpwind[11]+0.3535533905932737*fUpwind[6]*alphaSurf[11]+0.3162277660168379*alphaSurf[4]*fUpwind[10]+0.3162277660168379*fUpwind[4]*alphaSurf[10]+0.3535533905932737*alphaSurf[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alphaSurf[7]+0.3162277660168379*alphaSurf[1]*fUpwind[5]+0.3162277660168379*fUpwind[1]*alphaSurf[5]; 
  Ghat[14] = 0.282842712474619*alphaSurf[10]*fUpwind[19]+0.282842712474619*fUpwind[10]*alphaSurf[19]+0.3162277660168379*alphaSurf[15]*fUpwind[18]+0.2258769757263128*alphaSurf[12]*fUpwind[18]+0.3535533905932737*alphaSurf[1]*fUpwind[18]+0.3162277660168379*fUpwind[15]*alphaSurf[18]+0.2258769757263128*fUpwind[12]*alphaSurf[18]+0.3535533905932737*fUpwind[1]*alphaSurf[18]+0.3162277660168379*alphaSurf[11]*fUpwind[17]+0.3162277660168379*fUpwind[11]*alphaSurf[17]+0.2828427124746191*alphaSurf[6]*fUpwind[16]+0.2828427124746191*fUpwind[6]*alphaSurf[16]+0.3162277660168379*alphaSurf[9]*fUpwind[14]+0.2258769757263128*alphaSurf[8]*fUpwind[14]+0.3535533905932737*alphaSurf[0]*fUpwind[14]+0.3162277660168379*fUpwind[9]*alphaSurf[14]+0.2258769757263128*fUpwind[8]*alphaSurf[14]+0.3535533905932737*fUpwind[0]*alphaSurf[14]+0.3535533905932737*alphaSurf[5]*fUpwind[12]+0.3535533905932737*fUpwind[5]*alphaSurf[12]+0.3162277660168379*alphaSurf[4]*fUpwind[10]+0.3162277660168379*fUpwind[4]*alphaSurf[10]+0.3535533905932737*alphaSurf[3]*fUpwind[8]+0.3535533905932737*fUpwind[3]*alphaSurf[8]+0.3162277660168379*alphaSurf[2]*fUpwind[6]+0.3162277660168379*fUpwind[2]*alphaSurf[6]; 
  Ghat[15] = 0.2258769757263128*alphaSurf[16]*fUpwind[19]+0.3162277660168379*alphaSurf[11]*fUpwind[19]+0.3535533905932737*alphaSurf[2]*fUpwind[19]+0.2258769757263128*fUpwind[16]*alphaSurf[19]+0.3162277660168379*fUpwind[11]*alphaSurf[19]+0.3535533905932737*fUpwind[2]*alphaSurf[19]+0.3162277660168379*alphaSurf[14]*fUpwind[18]+0.3162277660168379*fUpwind[14]*alphaSurf[18]+0.282842712474619*alphaSurf[10]*fUpwind[17]+0.282842712474619*fUpwind[10]*alphaSurf[17]+0.3535533905932737*alphaSurf[4]*fUpwind[16]+0.3535533905932737*fUpwind[4]*alphaSurf[16]+0.2258769757263128*alphaSurf[9]*fUpwind[15]+0.3162277660168379*alphaSurf[7]*fUpwind[15]+0.3535533905932737*alphaSurf[0]*fUpwind[15]+0.2258769757263128*fUpwind[9]*alphaSurf[15]+0.3162277660168379*fUpwind[7]*alphaSurf[15]+0.3535533905932737*fUpwind[0]*alphaSurf[15]+0.2828427124746191*alphaSurf[5]*fUpwind[13]+0.2828427124746191*fUpwind[5]*alphaSurf[13]+0.3162277660168379*alphaSurf[6]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaSurf[10]+0.3535533905932737*alphaSurf[1]*fUpwind[9]+0.3535533905932737*fUpwind[1]*alphaSurf[9]+0.3162277660168379*alphaSurf[3]*fUpwind[5]+0.3162277660168379*fUpwind[3]*alphaSurf[5]; 
  Ghat[16] = 0.2258769757263128*alphaSurf[15]*fUpwind[19]+0.3162277660168379*alphaSurf[12]*fUpwind[19]+0.3535533905932737*alphaSurf[1]*fUpwind[19]+0.2258769757263128*fUpwind[15]*alphaSurf[19]+0.3162277660168379*fUpwind[12]*alphaSurf[19]+0.3535533905932737*fUpwind[1]*alphaSurf[19]+0.282842712474619*alphaSurf[10]*fUpwind[18]+0.282842712474619*fUpwind[10]*alphaSurf[18]+0.3162277660168379*alphaSurf[13]*fUpwind[17]+0.3162277660168379*fUpwind[13]*alphaSurf[17]+0.2258769757263128*alphaSurf[9]*fUpwind[16]+0.3162277660168379*alphaSurf[8]*fUpwind[16]+0.3535533905932737*alphaSurf[0]*fUpwind[16]+0.2258769757263128*fUpwind[9]*alphaSurf[16]+0.3162277660168379*fUpwind[8]*alphaSurf[16]+0.3535533905932737*fUpwind[0]*alphaSurf[16]+0.3535533905932737*alphaSurf[4]*fUpwind[15]+0.3535533905932737*fUpwind[4]*alphaSurf[15]+0.2828427124746191*alphaSurf[6]*fUpwind[14]+0.2828427124746191*fUpwind[6]*alphaSurf[14]+0.3162277660168379*alphaSurf[5]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaSurf[10]+0.3535533905932737*alphaSurf[2]*fUpwind[9]+0.3535533905932737*fUpwind[2]*alphaSurf[9]+0.3162277660168379*alphaSurf[3]*fUpwind[6]+0.3162277660168379*fUpwind[3]*alphaSurf[6]; 
  Ghat[17] = 0.2529822128134704*alphaSurf[18]*fUpwind[19]+0.2828427124746191*alphaSurf[5]*fUpwind[19]+0.2529822128134704*fUpwind[18]*alphaSurf[19]+0.2828427124746191*fUpwind[5]*alphaSurf[19]+0.2828427124746191*alphaSurf[4]*fUpwind[18]+0.2828427124746191*fUpwind[4]*alphaSurf[18]+0.3162277660168379*alphaSurf[9]*fUpwind[17]+0.3162277660168379*alphaSurf[8]*fUpwind[17]+0.2258769757263128*alphaSurf[7]*fUpwind[17]+0.3535533905932737*alphaSurf[0]*fUpwind[17]+0.3162277660168379*fUpwind[9]*alphaSurf[17]+0.3162277660168379*fUpwind[8]*alphaSurf[17]+0.2258769757263128*fUpwind[7]*alphaSurf[17]+0.3535533905932737*fUpwind[0]*alphaSurf[17]+0.3162277660168379*alphaSurf[13]*fUpwind[16]+0.3162277660168379*fUpwind[13]*alphaSurf[16]+0.282842712474619*alphaSurf[10]*fUpwind[15]+0.282842712474619*fUpwind[10]*alphaSurf[15]+0.3162277660168379*alphaSurf[11]*fUpwind[14]+0.3162277660168379*fUpwind[11]*alphaSurf[14]+0.2258769757263128*alphaSurf[11]*fUpwind[13]+0.3535533905932737*alphaSurf[2]*fUpwind[13]+0.2258769757263128*fUpwind[11]*alphaSurf[13]+0.3535533905932737*fUpwind[2]*alphaSurf[13]+0.282842712474619*alphaSurf[10]*fUpwind[12]+0.282842712474619*fUpwind[10]*alphaSurf[12]+0.3535533905932737*alphaSurf[3]*fUpwind[11]+0.3535533905932737*fUpwind[3]*alphaSurf[11]+0.3162277660168379*alphaSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[1]*alphaSurf[10]+0.3535533905932737*alphaSurf[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alphaSurf[7]+0.3162277660168379*alphaSurf[4]*fUpwind[5]+0.3162277660168379*fUpwind[4]*alphaSurf[5]; 
  Ghat[18] = 0.2529822128134704*alphaSurf[17]*fUpwind[19]+0.2828427124746191*alphaSurf[6]*fUpwind[19]+0.2529822128134704*fUpwind[17]*alphaSurf[19]+0.2828427124746191*fUpwind[6]*alphaSurf[19]+0.3162277660168379*alphaSurf[9]*fUpwind[18]+0.2258769757263128*alphaSurf[8]*fUpwind[18]+0.3162277660168379*alphaSurf[7]*fUpwind[18]+0.3535533905932737*alphaSurf[0]*fUpwind[18]+0.3162277660168379*fUpwind[9]*alphaSurf[18]+0.2258769757263128*fUpwind[8]*alphaSurf[18]+0.3162277660168379*fUpwind[7]*alphaSurf[18]+0.3535533905932737*fUpwind[0]*alphaSurf[18]+0.2828427124746191*alphaSurf[4]*fUpwind[17]+0.2828427124746191*fUpwind[4]*alphaSurf[17]+0.282842712474619*alphaSurf[10]*fUpwind[16]+0.282842712474619*fUpwind[10]*alphaSurf[16]+0.3162277660168379*alphaSurf[14]*fUpwind[15]+0.3162277660168379*fUpwind[14]*alphaSurf[15]+0.2258769757263128*alphaSurf[12]*fUpwind[14]+0.3535533905932737*alphaSurf[1]*fUpwind[14]+0.2258769757263128*fUpwind[12]*alphaSurf[14]+0.3535533905932737*fUpwind[1]*alphaSurf[14]+0.3162277660168379*alphaSurf[12]*fUpwind[13]+0.3162277660168379*fUpwind[12]*alphaSurf[13]+0.3535533905932737*alphaSurf[3]*fUpwind[12]+0.3535533905932737*fUpwind[3]*alphaSurf[12]+0.282842712474619*alphaSurf[10]*fUpwind[11]+0.282842712474619*fUpwind[10]*alphaSurf[11]+0.3162277660168379*alphaSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[2]*alphaSurf[10]+0.3535533905932737*alphaSurf[5]*fUpwind[8]+0.3535533905932737*fUpwind[5]*alphaSurf[8]+0.3162277660168379*alphaSurf[4]*fUpwind[6]+0.3162277660168379*fUpwind[4]*alphaSurf[6]; 
  Ghat[19] = 0.2258769757263128*alphaSurf[9]*fUpwind[19]+0.3162277660168379*alphaSurf[8]*fUpwind[19]+0.3162277660168379*alphaSurf[7]*fUpwind[19]+0.3535533905932737*alphaSurf[0]*fUpwind[19]+0.2258769757263128*fUpwind[9]*alphaSurf[19]+0.3162277660168379*fUpwind[8]*alphaSurf[19]+0.3162277660168379*fUpwind[7]*alphaSurf[19]+0.3535533905932737*fUpwind[0]*alphaSurf[19]+0.2529822128134704*alphaSurf[17]*fUpwind[18]+0.2828427124746191*alphaSurf[6]*fUpwind[18]+0.2529822128134704*fUpwind[17]*alphaSurf[18]+0.2828427124746191*fUpwind[6]*alphaSurf[18]+0.2828427124746191*alphaSurf[5]*fUpwind[17]+0.2828427124746191*fUpwind[5]*alphaSurf[17]+0.2258769757263128*alphaSurf[15]*fUpwind[16]+0.3162277660168379*alphaSurf[12]*fUpwind[16]+0.3535533905932737*alphaSurf[1]*fUpwind[16]+0.2258769757263128*fUpwind[15]*alphaSurf[16]+0.3162277660168379*fUpwind[12]*alphaSurf[16]+0.3535533905932737*fUpwind[1]*alphaSurf[16]+0.3162277660168379*alphaSurf[11]*fUpwind[15]+0.3535533905932737*alphaSurf[2]*fUpwind[15]+0.3162277660168379*fUpwind[11]*alphaSurf[15]+0.3535533905932737*fUpwind[2]*alphaSurf[15]+0.282842712474619*alphaSurf[10]*fUpwind[14]+0.282842712474619*fUpwind[10]*alphaSurf[14]+0.282842712474619*alphaSurf[10]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaSurf[13]+0.3162277660168379*alphaSurf[3]*fUpwind[10]+0.3162277660168379*fUpwind[3]*alphaSurf[10]+0.3535533905932737*alphaSurf[4]*fUpwind[9]+0.3535533905932737*fUpwind[4]*alphaSurf[9]+0.3162277660168379*alphaSurf[5]*fUpwind[6]+0.3162277660168379*fUpwind[5]*alphaSurf[6]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += 0.7071067811865475*Ghat[3]*dv1par; 
  out[4] += -1.224744871391589*Ghat[0]*dv1par; 
  out[5] += 0.7071067811865475*Ghat[4]*dv1par; 
  out[6] += 0.7071067811865475*Ghat[5]*dv1par; 
  out[7] += 0.7071067811865475*Ghat[6]*dv1par; 
  out[8] += -1.224744871391589*Ghat[1]*dv1par; 
  out[9] += -1.224744871391589*Ghat[2]*dv1par; 
  out[10] += -1.224744871391589*Ghat[3]*dv1par; 
  out[11] += 0.7071067811865475*Ghat[7]*dv1par; 
  out[12] += 0.7071067811865475*Ghat[8]*dv1par; 
  out[13] += 0.7071067811865475*Ghat[9]*dv1par; 
  out[14] += 1.58113883008419*Ghat[0]*dv1par; 
  out[15] += 0.7071067811865475*Ghat[10]*dv1par; 
  out[16] += -1.224744871391589*Ghat[4]*dv1par; 
  out[17] += -1.224744871391589*Ghat[5]*dv1par; 
  out[18] += -1.224744871391589*Ghat[6]*dv1par; 
  out[19] += 0.7071067811865475*Ghat[11]*dv1par; 
  out[20] += 0.7071067811865475*Ghat[12]*dv1par; 
  out[21] += 0.7071067811865475*Ghat[13]*dv1par; 
  out[22] += 0.7071067811865475*Ghat[14]*dv1par; 
  out[23] += 0.7071067811865475*Ghat[15]*dv1par; 
  out[24] += 0.7071067811865475*Ghat[16]*dv1par; 
  out[25] += -1.224744871391589*Ghat[7]*dv1par; 
  out[26] += -1.224744871391589*Ghat[8]*dv1par; 
  out[27] += -1.224744871391589*Ghat[9]*dv1par; 
  out[28] += 1.58113883008419*Ghat[1]*dv1par; 
  out[29] += 1.58113883008419*Ghat[2]*dv1par; 
  out[30] += 1.58113883008419*Ghat[3]*dv1par; 
  out[31] += -1.224744871391589*Ghat[10]*dv1par; 
  out[32] += 0.7071067811865475*Ghat[17]*dv1par; 
  out[33] += 0.7071067811865475*Ghat[18]*dv1par; 
  out[34] += 0.7071067811865475*Ghat[19]*dv1par; 
  out[35] += -1.224744871391589*Ghat[11]*dv1par; 
  out[36] += -1.224744871391589*Ghat[12]*dv1par; 
  out[37] += -1.224744871391589*Ghat[13]*dv1par; 
  out[38] += -1.224744871391589*Ghat[14]*dv1par; 
  out[39] += -1.224744871391589*Ghat[15]*dv1par; 
  out[40] += -1.224744871391589*Ghat[16]*dv1par; 
  out[41] += 1.58113883008419*Ghat[4]*dv1par; 
  out[42] += 1.58113883008419*Ghat[5]*dv1par; 
  out[43] += 1.58113883008419*Ghat[6]*dv1par; 
  out[44] += -1.224744871391589*Ghat[17]*dv1par; 
  out[45] += -1.224744871391589*Ghat[18]*dv1par; 
  out[46] += -1.224744871391589*Ghat[19]*dv1par; 
  out[47] += 1.58113883008419*Ghat[10]*dv1par; 

  } 
} 
