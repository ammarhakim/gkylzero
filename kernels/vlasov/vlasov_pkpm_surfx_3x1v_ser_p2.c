#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_3x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // bvarl/bvarc/bvarr:  Input magnetic field unit vector in left/center/right cells.
  // u_il/u_ic/u_ir:  Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *ul = &u_il[0]; 
  const double *uc = &u_ic[0]; 
  const double *ur = &u_ir[0]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  double alpha_l[48] = {0.0}; 
  double alpha_c[48] = {0.0}; 
  double alpha_r[48] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar+1.414213562373095*ul[0]; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar+1.414213562373095*ul[1]; 
  alpha_l[2] = 1.414213562373095*bl[2]*wvpar+1.414213562373095*ul[2]; 
  alpha_l[3] = 1.414213562373095*bl[3]*wvpar+1.414213562373095*ul[3]; 
  alpha_l[4] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[5] = 1.414213562373095*bl[4]*wvpar+1.414213562373095*ul[4]; 
  alpha_l[6] = 1.414213562373095*bl[5]*wvpar+1.414213562373095*ul[5]; 
  alpha_l[7] = 1.414213562373095*bl[6]*wvpar+1.414213562373095*ul[6]; 
  alpha_l[8] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[9] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[10] = 0.408248290463863*bl[3]*dvpar; 
  alpha_l[11] = 1.414213562373095*bl[7]*wvpar+1.414213562373095*ul[7]; 
  alpha_l[12] = 1.414213562373095*bl[8]*wvpar+1.414213562373095*ul[8]; 
  alpha_l[13] = 1.414213562373095*bl[9]*wvpar+1.414213562373095*ul[9]; 
  alpha_l[15] = 1.414213562373095*bl[10]*wvpar+1.414213562373095*ul[10]; 
  alpha_l[16] = 0.408248290463863*bl[4]*dvpar; 
  alpha_l[17] = 0.408248290463863*bl[5]*dvpar; 
  alpha_l[18] = 0.408248290463863*bl[6]*dvpar; 
  alpha_l[19] = 1.414213562373095*bl[11]*wvpar+1.414213562373095*ul[11]; 
  alpha_l[20] = 1.414213562373095*bl[12]*wvpar+1.414213562373095*ul[12]; 
  alpha_l[21] = 1.414213562373095*bl[13]*wvpar+1.414213562373095*ul[13]; 
  alpha_l[22] = 1.414213562373095*bl[14]*wvpar+1.414213562373095*ul[14]; 
  alpha_l[23] = 1.414213562373095*bl[15]*wvpar+1.414213562373095*ul[15]; 
  alpha_l[24] = 1.414213562373095*bl[16]*wvpar+1.414213562373095*ul[16]; 
  alpha_l[25] = 0.408248290463863*bl[7]*dvpar; 
  alpha_l[26] = 0.408248290463863*bl[8]*dvpar; 
  alpha_l[27] = 0.408248290463863*bl[9]*dvpar; 
  alpha_l[31] = 0.408248290463863*bl[10]*dvpar; 
  alpha_l[32] = 1.414213562373095*bl[17]*wvpar+1.414213562373095*ul[17]; 
  alpha_l[33] = 1.414213562373095*bl[18]*wvpar+1.414213562373095*ul[18]; 
  alpha_l[34] = 1.414213562373095*bl[19]*wvpar+1.414213562373095*ul[19]; 
  alpha_l[35] = 0.408248290463863*bl[11]*dvpar; 
  alpha_l[36] = 0.408248290463863*bl[12]*dvpar; 
  alpha_l[37] = 0.408248290463863*bl[13]*dvpar; 
  alpha_l[38] = 0.408248290463863*bl[14]*dvpar; 
  alpha_l[39] = 0.408248290463863*bl[15]*dvpar; 
  alpha_l[40] = 0.408248290463863*bl[16]*dvpar; 
  alpha_l[44] = 0.408248290463863*bl[17]*dvpar; 
  alpha_l[45] = 0.408248290463863*bl[18]*dvpar; 
  alpha_l[46] = 0.408248290463863*bl[19]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar+1.414213562373095*uc[0]; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar+1.414213562373095*uc[1]; 
  alpha_c[2] = 1.414213562373095*bc[2]*wvpar+1.414213562373095*uc[2]; 
  alpha_c[3] = 1.414213562373095*bc[3]*wvpar+1.414213562373095*uc[3]; 
  alpha_c[4] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[5] = 1.414213562373095*bc[4]*wvpar+1.414213562373095*uc[4]; 
  alpha_c[6] = 1.414213562373095*bc[5]*wvpar+1.414213562373095*uc[5]; 
  alpha_c[7] = 1.414213562373095*bc[6]*wvpar+1.414213562373095*uc[6]; 
  alpha_c[8] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[9] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[10] = 0.408248290463863*bc[3]*dvpar; 
  alpha_c[11] = 1.414213562373095*bc[7]*wvpar+1.414213562373095*uc[7]; 
  alpha_c[12] = 1.414213562373095*bc[8]*wvpar+1.414213562373095*uc[8]; 
  alpha_c[13] = 1.414213562373095*bc[9]*wvpar+1.414213562373095*uc[9]; 
  alpha_c[15] = 1.414213562373095*bc[10]*wvpar+1.414213562373095*uc[10]; 
  alpha_c[16] = 0.408248290463863*bc[4]*dvpar; 
  alpha_c[17] = 0.408248290463863*bc[5]*dvpar; 
  alpha_c[18] = 0.408248290463863*bc[6]*dvpar; 
  alpha_c[19] = 1.414213562373095*bc[11]*wvpar+1.414213562373095*uc[11]; 
  alpha_c[20] = 1.414213562373095*bc[12]*wvpar+1.414213562373095*uc[12]; 
  alpha_c[21] = 1.414213562373095*bc[13]*wvpar+1.414213562373095*uc[13]; 
  alpha_c[22] = 1.414213562373095*bc[14]*wvpar+1.414213562373095*uc[14]; 
  alpha_c[23] = 1.414213562373095*bc[15]*wvpar+1.414213562373095*uc[15]; 
  alpha_c[24] = 1.414213562373095*bc[16]*wvpar+1.414213562373095*uc[16]; 
  alpha_c[25] = 0.408248290463863*bc[7]*dvpar; 
  alpha_c[26] = 0.408248290463863*bc[8]*dvpar; 
  alpha_c[27] = 0.408248290463863*bc[9]*dvpar; 
  alpha_c[31] = 0.408248290463863*bc[10]*dvpar; 
  alpha_c[32] = 1.414213562373095*bc[17]*wvpar+1.414213562373095*uc[17]; 
  alpha_c[33] = 1.414213562373095*bc[18]*wvpar+1.414213562373095*uc[18]; 
  alpha_c[34] = 1.414213562373095*bc[19]*wvpar+1.414213562373095*uc[19]; 
  alpha_c[35] = 0.408248290463863*bc[11]*dvpar; 
  alpha_c[36] = 0.408248290463863*bc[12]*dvpar; 
  alpha_c[37] = 0.408248290463863*bc[13]*dvpar; 
  alpha_c[38] = 0.408248290463863*bc[14]*dvpar; 
  alpha_c[39] = 0.408248290463863*bc[15]*dvpar; 
  alpha_c[40] = 0.408248290463863*bc[16]*dvpar; 
  alpha_c[44] = 0.408248290463863*bc[17]*dvpar; 
  alpha_c[45] = 0.408248290463863*bc[18]*dvpar; 
  alpha_c[46] = 0.408248290463863*bc[19]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar+1.414213562373095*ur[0]; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar+1.414213562373095*ur[1]; 
  alpha_r[2] = 1.414213562373095*br[2]*wvpar+1.414213562373095*ur[2]; 
  alpha_r[3] = 1.414213562373095*br[3]*wvpar+1.414213562373095*ur[3]; 
  alpha_r[4] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[5] = 1.414213562373095*br[4]*wvpar+1.414213562373095*ur[4]; 
  alpha_r[6] = 1.414213562373095*br[5]*wvpar+1.414213562373095*ur[5]; 
  alpha_r[7] = 1.414213562373095*br[6]*wvpar+1.414213562373095*ur[6]; 
  alpha_r[8] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[9] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[10] = 0.408248290463863*br[3]*dvpar; 
  alpha_r[11] = 1.414213562373095*br[7]*wvpar+1.414213562373095*ur[7]; 
  alpha_r[12] = 1.414213562373095*br[8]*wvpar+1.414213562373095*ur[8]; 
  alpha_r[13] = 1.414213562373095*br[9]*wvpar+1.414213562373095*ur[9]; 
  alpha_r[15] = 1.414213562373095*br[10]*wvpar+1.414213562373095*ur[10]; 
  alpha_r[16] = 0.408248290463863*br[4]*dvpar; 
  alpha_r[17] = 0.408248290463863*br[5]*dvpar; 
  alpha_r[18] = 0.408248290463863*br[6]*dvpar; 
  alpha_r[19] = 1.414213562373095*br[11]*wvpar+1.414213562373095*ur[11]; 
  alpha_r[20] = 1.414213562373095*br[12]*wvpar+1.414213562373095*ur[12]; 
  alpha_r[21] = 1.414213562373095*br[13]*wvpar+1.414213562373095*ur[13]; 
  alpha_r[22] = 1.414213562373095*br[14]*wvpar+1.414213562373095*ur[14]; 
  alpha_r[23] = 1.414213562373095*br[15]*wvpar+1.414213562373095*ur[15]; 
  alpha_r[24] = 1.414213562373095*br[16]*wvpar+1.414213562373095*ur[16]; 
  alpha_r[25] = 0.408248290463863*br[7]*dvpar; 
  alpha_r[26] = 0.408248290463863*br[8]*dvpar; 
  alpha_r[27] = 0.408248290463863*br[9]*dvpar; 
  alpha_r[31] = 0.408248290463863*br[10]*dvpar; 
  alpha_r[32] = 1.414213562373095*br[17]*wvpar+1.414213562373095*ur[17]; 
  alpha_r[33] = 1.414213562373095*br[18]*wvpar+1.414213562373095*ur[18]; 
  alpha_r[34] = 1.414213562373095*br[19]*wvpar+1.414213562373095*ur[19]; 
  alpha_r[35] = 0.408248290463863*br[11]*dvpar; 
  alpha_r[36] = 0.408248290463863*br[12]*dvpar; 
  alpha_r[37] = 0.408248290463863*br[13]*dvpar; 
  alpha_r[38] = 0.408248290463863*br[14]*dvpar; 
  alpha_r[39] = 0.408248290463863*br[15]*dvpar; 
  alpha_r[40] = 0.408248290463863*br[16]*dvpar; 
  alpha_r[44] = 0.408248290463863*br[17]*dvpar; 
  alpha_r[45] = 0.408248290463863*br[18]*dvpar; 
  alpha_r[46] = 0.408248290463863*br[19]*dvpar; 

  double alphaSurf_l[20] = {0.0}; 
  alphaSurf_l[0] = 0.3458741190809163*alpha_l[11]+0.3458741190809163*alpha_c[11]+0.4975526040028326*alpha_l[1]-0.4975526040028326*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.3458741190809163*alpha_l[19]+0.3458741190809163*alpha_c[19]+0.4975526040028326*alpha_l[5]-0.4975526040028326*alpha_c[5]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_l[2] = 0.3458741190809163*alpha_l[21]+0.3458741190809163*alpha_c[21]+0.4975526040028326*alpha_l[6]-0.4975526040028326*alpha_c[6]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.3458741190809163*alpha_l[25]+0.3458741190809163*alpha_c[25]+0.4975526040028326*alpha_l[8]-0.4975526040028326*alpha_c[8]+0.3535533905932737*alpha_l[4]+0.3535533905932737*alpha_c[4]; 
  alphaSurf_l[4] = 0.3458741190809163*alpha_l[32]+0.3458741190809163*alpha_c[32]+0.4975526040028326*alpha_l[15]-0.4975526040028326*alpha_c[15]+0.3535533905932737*alpha_l[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_l[5] = 0.3458741190809163*alpha_l[35]+0.3458741190809163*alpha_c[35]+0.4975526040028326*alpha_l[16]-0.4975526040028326*alpha_c[16]+0.3535533905932737*alpha_l[9]+0.3535533905932737*alpha_c[9]; 
  alphaSurf_l[6] = 0.3458741190809163*alpha_l[37]+0.3458741190809163*alpha_c[37]+0.4975526040028326*alpha_l[17]-0.4975526040028326*alpha_c[17]+0.3535533905932737*alpha_l[10]+0.3535533905932737*alpha_c[10]; 
  alphaSurf_l[7] = 0.4975526040028326*alpha_l[20]-0.4975526040028326*alpha_c[20]+0.3535533905932737*alpha_l[12]+0.3535533905932737*alpha_c[12]; 
  alphaSurf_l[8] = 0.4975526040028326*alpha_l[23]-0.4975526040028326*alpha_c[23]+0.3535533905932737*alpha_l[13]+0.3535533905932737*alpha_c[13]; 
  alphaSurf_l[10] = 0.3458741190809163*alpha_l[44]+0.3458741190809163*alpha_c[44]+0.4975526040028326*alpha_l[31]-0.4975526040028326*alpha_c[31]+0.3535533905932737*alpha_l[18]+0.3535533905932737*alpha_c[18]; 
  alphaSurf_l[11] = 0.4975526040028326*alpha_l[33]-0.4975526040028326*alpha_c[33]+0.3535533905932737*alpha_l[22]+0.3535533905932737*alpha_c[22]; 
  alphaSurf_l[12] = 0.4975526040028326*alpha_l[34]-0.4975526040028326*alpha_c[34]+0.3535533905932737*alpha_l[24]+0.3535533905932737*alpha_c[24]; 
  alphaSurf_l[13] = 0.4975526040028326*alpha_l[36]-0.4975526040028326*alpha_c[36]+0.3535533905932737*alpha_l[26]+0.3535533905932737*alpha_c[26]; 
  alphaSurf_l[14] = 0.4975526040028326*alpha_l[39]-0.4975526040028326*alpha_c[39]+0.3535533905932737*alpha_l[27]+0.3535533905932737*alpha_c[27]; 
  alphaSurf_l[17] = 0.4975526040028326*alpha_l[45]-0.4975526040028326*alpha_c[45]+0.3535533905932737*alpha_l[38]+0.3535533905932737*alpha_c[38]; 
  alphaSurf_l[18] = 0.4975526040028326*alpha_l[46]-0.4975526040028326*alpha_c[46]+0.3535533905932737*alpha_l[40]+0.3535533905932737*alpha_c[40]; 

  double alphaSurf_r[20] = {0.0}; 
  alphaSurf_r[0] = 0.3458741190809163*alpha_r[11]+0.3458741190809163*alpha_c[11]-0.4975526040028326*alpha_r[1]+0.4975526040028326*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = 0.3458741190809163*alpha_r[19]+0.3458741190809163*alpha_c[19]-0.4975526040028326*alpha_r[5]+0.4975526040028326*alpha_c[5]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_r[2] = 0.3458741190809163*alpha_r[21]+0.3458741190809163*alpha_c[21]-0.4975526040028326*alpha_r[6]+0.4975526040028326*alpha_c[6]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = 0.3458741190809163*alpha_r[25]+0.3458741190809163*alpha_c[25]-0.4975526040028326*alpha_r[8]+0.4975526040028326*alpha_c[8]+0.3535533905932737*alpha_r[4]+0.3535533905932737*alpha_c[4]; 
  alphaSurf_r[4] = 0.3458741190809163*alpha_r[32]+0.3458741190809163*alpha_c[32]-0.4975526040028326*alpha_r[15]+0.4975526040028326*alpha_c[15]+0.3535533905932737*alpha_r[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_r[5] = 0.3458741190809163*alpha_r[35]+0.3458741190809163*alpha_c[35]-0.4975526040028326*alpha_r[16]+0.4975526040028326*alpha_c[16]+0.3535533905932737*alpha_r[9]+0.3535533905932737*alpha_c[9]; 
  alphaSurf_r[6] = 0.3458741190809163*alpha_r[37]+0.3458741190809163*alpha_c[37]-0.4975526040028326*alpha_r[17]+0.4975526040028326*alpha_c[17]+0.3535533905932737*alpha_r[10]+0.3535533905932737*alpha_c[10]; 
  alphaSurf_r[7] = (-0.4975526040028326*alpha_r[20])+0.4975526040028326*alpha_c[20]+0.3535533905932737*alpha_r[12]+0.3535533905932737*alpha_c[12]; 
  alphaSurf_r[8] = (-0.4975526040028326*alpha_r[23])+0.4975526040028326*alpha_c[23]+0.3535533905932737*alpha_r[13]+0.3535533905932737*alpha_c[13]; 
  alphaSurf_r[10] = 0.3458741190809163*alpha_r[44]+0.3458741190809163*alpha_c[44]-0.4975526040028326*alpha_r[31]+0.4975526040028326*alpha_c[31]+0.3535533905932737*alpha_r[18]+0.3535533905932737*alpha_c[18]; 
  alphaSurf_r[11] = (-0.4975526040028326*alpha_r[33])+0.4975526040028326*alpha_c[33]+0.3535533905932737*alpha_r[22]+0.3535533905932737*alpha_c[22]; 
  alphaSurf_r[12] = (-0.4975526040028326*alpha_r[34])+0.4975526040028326*alpha_c[34]+0.3535533905932737*alpha_r[24]+0.3535533905932737*alpha_c[24]; 
  alphaSurf_r[13] = (-0.4975526040028326*alpha_r[36])+0.4975526040028326*alpha_c[36]+0.3535533905932737*alpha_r[26]+0.3535533905932737*alpha_c[26]; 
  alphaSurf_r[14] = (-0.4975526040028326*alpha_r[39])+0.4975526040028326*alpha_c[39]+0.3535533905932737*alpha_r[27]+0.3535533905932737*alpha_c[27]; 
  alphaSurf_r[17] = (-0.4975526040028326*alpha_r[45])+0.4975526040028326*alpha_c[45]+0.3535533905932737*alpha_r[38]+0.3535533905932737*alpha_c[38]; 
  alphaSurf_r[18] = (-0.4975526040028326*alpha_r[46])+0.4975526040028326*alpha_c[46]+0.3535533905932737*alpha_r[40]+0.3535533905932737*alpha_c[40]; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17])-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5]+alphaSurf_l[4])-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fc); 
  } 
  if (0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17])-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5]+alphaSurf_r[4])-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if ((-0.4242640687119281*alphaSurf_l[12])-0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fc); 
  } 
  if ((-0.4242640687119281*alphaSurf_r[12])-0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if ((-0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17]))+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5])+0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fc); 
  } 
  if ((-0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17]))+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5])+0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[18])+0.5303300858899104*alphaSurf_l[14]-0.4242640687119285*alphaSurf_l[13]+0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[5]-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[18])+0.5303300858899104*alphaSurf_r[14]-0.4242640687119285*alphaSurf_r[13]+0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[5]-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fc); 
  } 
  if (0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[18]-0.5303300858899104*alphaSurf_l[14]+0.4242640687119285*alphaSurf_l[13]+0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[5]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[18]-0.5303300858899104*alphaSurf_r[14]+0.4242640687119285*alphaSurf_r[13]+0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[5]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  if (0.5692099788303082*alphaSurf_l[18]-0.5692099788303082*alphaSurf_l[17]-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[6]+0.6363961030678926*alphaSurf_l[5]-0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*alphaSurf_l[2]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fc); 
  } 
  if (0.5692099788303082*alphaSurf_r[18]-0.5692099788303082*alphaSurf_r[17]-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[6]+0.6363961030678926*alphaSurf_r[5]-0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*alphaSurf_r[2]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  if ((-0.4242640687119281*alphaSurf_l[12])+0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[2]-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fc); 
  } 
  if ((-0.4242640687119281*alphaSurf_r[12])+0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[2]-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  if ((-0.5692099788303082*alphaSurf_l[18])+0.5692099788303082*alphaSurf_l[17]+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])-0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[6]-0.6363961030678926*(alphaSurf_l[5]+alphaSurf_l[4])+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])-0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fc); 
  } 
  if ((-0.5692099788303082*alphaSurf_r[18])+0.5692099788303082*alphaSurf_r[17]+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])-0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[6]-0.6363961030678926*(alphaSurf_r[5]+alphaSurf_r[4])+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])-0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[17])-0.4242640687119285*alphaSurf_l[14]+0.5303300858899104*(alphaSurf_l[13]+alphaSurf_l[11])+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[6]-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[17])-0.4242640687119285*alphaSurf_r[14]+0.5303300858899104*(alphaSurf_r[13]+alphaSurf_r[11])+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[6]-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fc); 
  } 
  if (0.5303300858899104*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[17]+0.4242640687119285*alphaSurf_l[14]-0.5303300858899104*alphaSurf_l[13]+0.5303300858899104*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[6]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[17]+0.4242640687119285*alphaSurf_r[14]-0.5303300858899104*alphaSurf_r[13]+0.5303300858899104*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[6]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*(alphaSurf_l[14]+alphaSurf_l[13])-0.3952847075210473*(alphaSurf_l[8]+alphaSurf_l[7])-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fc); 
  } 
  if (0.5303300858899104*(alphaSurf_r[14]+alphaSurf_r[13])-0.3952847075210473*(alphaSurf_r[8]+alphaSurf_r[7])-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alphaSurf_l[0]-0.3952847075210473*(alphaSurf_l[8]+alphaSurf_l[7]) > 0) { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fc); 
  } 
  if (0.3535533905932737*alphaSurf_r[0]-0.3952847075210473*(alphaSurf_r[8]+alphaSurf_r[7]) > 0) { 
    fUpwindQuad_r[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fr); 
  } 
  if ((-0.5303300858899104*(alphaSurf_l[14]+alphaSurf_l[13]))-0.3952847075210473*(alphaSurf_l[8]+alphaSurf_l[7])+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fc); 
  } 
  if ((-0.5303300858899104*(alphaSurf_r[14]+alphaSurf_r[13]))-0.3952847075210473*(alphaSurf_r[8]+alphaSurf_r[7])+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[17]-0.4242640687119285*alphaSurf_l[14]+0.5303300858899104*alphaSurf_l[13]-0.5303300858899104*alphaSurf_l[11]+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[6]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[17]-0.4242640687119285*alphaSurf_r[14]+0.5303300858899104*alphaSurf_r[13]-0.5303300858899104*alphaSurf_r[11]+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[6]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fr); 
  } 
  if ((-0.5303300858899104*alphaSurf_l[11])+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]+0.4743416490252568*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fc); 
  } 
  if ((-0.5303300858899104*alphaSurf_r[11])+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]+0.4743416490252568*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_r[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[17])+0.4242640687119285*alphaSurf_l[14]-0.5303300858899104*(alphaSurf_l[13]+alphaSurf_l[11])+0.3162277660168379*alphaSurf_l[8]-0.3952847075210473*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[6]+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[17])+0.4242640687119285*alphaSurf_r[14]-0.5303300858899104*(alphaSurf_r[13]+alphaSurf_r[11])+0.3162277660168379*alphaSurf_r[8]-0.3952847075210473*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[6]+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fr); 
  } 
  if ((-0.5692099788303082*alphaSurf_l[18])+0.5692099788303082*alphaSurf_l[17]-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[6]-0.6363961030678926*(alphaSurf_l[5]+alphaSurf_l[4])-0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2])+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fl); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fc); 
  } 
  if ((-0.5692099788303082*alphaSurf_r[18])+0.5692099788303082*alphaSurf_r[17]-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[6]-0.6363961030678926*(alphaSurf_r[5]+alphaSurf_r[4])-0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2])+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_r[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fr); 
  } 
  if (0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[2]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fc); 
  } 
  if (0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[2]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_r[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fr); 
  } 
  if (0.5692099788303082*alphaSurf_l[18]-0.5692099788303082*alphaSurf_l[17]+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]-0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*alphaSurf_l[6]+0.6363961030678926*alphaSurf_l[5]-0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.4743416490252568*alphaSurf_l[2]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fl); 
  } else { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fc); 
  } 
  if (0.5692099788303082*alphaSurf_r[18]-0.5692099788303082*alphaSurf_r[17]+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]-0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*alphaSurf_r[6]+0.6363961030678926*alphaSurf_r[5]-0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.4743416490252568*alphaSurf_r[2]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_r[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fr); 
  } 
  if (0.711512473537885*alphaSurf_l[18]+0.5303300858899104*alphaSurf_l[14]-0.4242640687119285*alphaSurf_l[13]-0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]-0.6363961030678926*alphaSurf_l[5]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fl); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fc); 
  } 
  if (0.711512473537885*alphaSurf_r[18]+0.5303300858899104*alphaSurf_r[14]-0.4242640687119285*alphaSurf_r[13]-0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]-0.6363961030678926*alphaSurf_r[5]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_r[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fr); 
  } 
  if ((-0.5303300858899104*alphaSurf_l[12])-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]+0.4743416490252568*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fl); 
  } else { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fc); 
  } 
  if ((-0.5303300858899104*alphaSurf_r[12])-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]+0.4743416490252568*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_r[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fr); 
  } 
  if ((-0.711512473537885*alphaSurf_l[18])-0.5303300858899104*alphaSurf_l[14]+0.4242640687119285*alphaSurf_l[13]-0.5303300858899104*alphaSurf_l[12]-0.3952847075210473*alphaSurf_l[8]+0.3162277660168379*alphaSurf_l[7]+0.6363961030678926*alphaSurf_l[5]+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fc); 
  } 
  if ((-0.711512473537885*alphaSurf_r[18])-0.5303300858899104*alphaSurf_r[14]+0.4242640687119285*alphaSurf_r[13]-0.5303300858899104*alphaSurf_r[12]-0.3952847075210473*alphaSurf_r[8]+0.3162277660168379*alphaSurf_r[7]+0.6363961030678926*alphaSurf_r[5]+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_r[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fr); 
  } 
  if ((-0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17]))-0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]-0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])-0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5])+0.6363961030678926*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fl); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fc); 
  } 
  if ((-0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17]))-0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]-0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])-0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5])+0.6363961030678926*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_r[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fr); 
  } 
  if (0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*alphaSurf_l[4]+0.4743416490252568*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fl); 
  } else { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fc); 
  } 
  if (0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*alphaSurf_r[4]+0.4743416490252568*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_r[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fr); 
  } 
  if (0.5692099788303082*(alphaSurf_l[18]+alphaSurf_l[17])+0.4242640687119285*(alphaSurf_l[14]+alphaSurf_l[13])+0.4242640687119281*alphaSurf_l[12]+0.4242640687119285*alphaSurf_l[11]+0.853814968245462*alphaSurf_l[10]+0.3162277660168379*(alphaSurf_l[8]+alphaSurf_l[7])+0.6363961030678926*(alphaSurf_l[6]+alphaSurf_l[5]+alphaSurf_l[4])+0.4743416490252568*(alphaSurf_l[3]+alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fc); 
  } 
  if (0.5692099788303082*(alphaSurf_r[18]+alphaSurf_r[17])+0.4242640687119285*(alphaSurf_r[14]+alphaSurf_r[13])+0.4242640687119281*alphaSurf_r[12]+0.4242640687119285*alphaSurf_r[11]+0.853814968245462*alphaSurf_r[10]+0.3162277660168379*(alphaSurf_r[8]+alphaSurf_r[7])+0.6363961030678926*(alphaSurf_r[6]+alphaSurf_r[5]+alphaSurf_r[4])+0.4743416490252568*(alphaSurf_r[3]+alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_r[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fr); 
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
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx1; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx1; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx1; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx1; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx1; 
  out[11] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx1; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx1; 
  out[13] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx1; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dx1; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx1; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx1; 
  out[17] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx1; 
  out[18] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dx1; 
  out[19] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx1; 
  out[20] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx1; 
  out[21] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dx1; 
  out[22] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dx1; 
  out[23] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx1; 
  out[24] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dx1; 
  out[25] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dx1; 
  out[26] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dx1; 
  out[27] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dx1; 
  out[28] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dx1; 
  out[29] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dx1; 
  out[30] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dx1; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dx1; 
  out[32] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dx1; 
  out[33] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dx1; 
  out[34] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dx1; 
  out[35] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dx1; 
  out[36] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dx1; 
  out[37] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dx1; 
  out[38] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dx1; 
  out[39] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dx1; 
  out[40] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dx1; 
  out[41] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dx1; 
  out[42] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dx1; 
  out[43] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dx1; 
  out[44] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dx1; 
  out[45] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dx1; 
  out[46] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dx1; 
  out[47] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dx1; 

} 
