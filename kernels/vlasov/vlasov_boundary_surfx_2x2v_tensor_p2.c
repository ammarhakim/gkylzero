#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_2x2v_tensor_p2(const double *w, const double *dxv, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // alpha_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell (used by general geometry version).
  // alpha_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell (used by general geometry version).
  // sgn_alpha_surf_edge: sign(alpha_surf_edge) at quadrature points (used by general geometry version).
  // sgn_alpha_surf_skin: sign(alpha_surf_skin) at quadrature points (used by general geometry version).
  // const_sgn_alpha_edge: Boolean array true if sign(alpha_surf_edge) is only one sign, either +1 or -1 (used by general geometry version).
  // const_sgn_alpha_skin: Boolean array true if sign(alpha_surf_skin) is only one sign, either +1 or -1 (used by general geometry version).
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[2], wv = w[2]; 
  double Ghat[27]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fskin[11]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*wv+(0.4564354645876384*fskin[21]+0.3535533905932737*fskin[6]+0.2041241452319315*fskin[3])*dv; 
  Ghat[1] = (1.58113883008419*fskin[19]+1.224744871391589*fskin[5]+0.7071067811865475*fskin[2])*wv+(0.4564354645876384*fskin[32]+0.3535533905932737*fskin[15]+0.2041241452319315*fskin[7])*dv; 
  Ghat[2] = (1.58113883008419*fskin[21]+1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])*wv+(0.408248290463863*fskin[45]+0.3162277660168379*fskin[23]+0.1825741858350554*fskin[13]+0.4564354645876384*fskin[11]+0.3535533905932737*fskin[1]+0.2041241452319315*fskin[0])*dv; 
  Ghat[3] = (1.58113883008419*fskin[25]+1.224744871391589*fskin[8]+0.7071067811865475*fskin[4])*wv+(0.4564354645876384*fskin[37]+0.3535533905932737*fskin[17]+0.2041241452319315*fskin[10])*dv; 
  Ghat[4] = (1.58113883008419*fskin[32]+1.224744871391589*fskin[15]+0.7071067811865475*fskin[7])*wv+(0.408248290463863*fskin[55]+0.3162277660168379*fskin[34]+0.1825741858350554*fskin[24]+0.4564354645876384*fskin[19]+0.3535533905932737*fskin[5]+0.2041241452319315*fskin[2])*dv; 
  Ghat[5] = (1.58113883008419*fskin[35]+1.224744871391589*fskin[16]+0.7071067811865475*fskin[9])*wv+(0.4564354645876384*fskin[50]+0.3535533905932737*fskin[31]+0.2041241452319315*fskin[18])*dv; 
  Ghat[6] = (1.58113883008419*fskin[37]+1.224744871391589*fskin[17]+0.7071067811865475*fskin[10])*wv+(0.408248290463863*fskin[58]+0.3162277660168379*fskin[39]+0.1825741858350554*fskin[27]+0.4564354645876384*fskin[25]+0.3535533905932737*fskin[8]+0.2041241452319315*fskin[4])*dv; 
  Ghat[7] = (1.58113883008419*fskin[44]+1.224744871391589*fskin[20]+0.7071067811865475*fskin[12])*wv+(0.4564354645876384*fskin[54]+0.3535533905932737*fskin[33]+0.2041241452319315*fskin[22])*dv; 
  Ghat[8] = (1.58113883008419*fskin[45]+1.224744871391589*fskin[23]+0.7071067811865475*fskin[13])*wv+(0.408248290463863*fskin[21]+0.3162277660168379*fskin[6]+0.1825741858350554*fskin[3])*dv; 
  Ghat[9] = (1.58113883008419*fskin[47]+1.224744871391589*fskin[28]+0.7071067811865475*fskin[14])*wv+(0.4564354645876384*fskin[62]+0.3535533905932737*fskin[42]+0.2041241452319315*fskin[30])*dv; 
  Ghat[10] = (1.58113883008419*fskin[50]+1.224744871391589*fskin[31]+0.7071067811865475*fskin[18])*wv+(0.408248290463863*fskin[67]+0.3162277660168379*fskin[52]+0.1825741858350554*fskin[40]+0.4564354645876384*fskin[35]+0.3535533905932737*fskin[16]+0.2041241452319315*fskin[9])*dv; 
  Ghat[11] = (1.58113883008419*fskin[54]+1.224744871391589*fskin[33]+0.7071067811865475*fskin[22])*wv+(0.408248290463863*fskin[72]+0.3162277660168379*fskin[56]+0.1825741858350554*fskin[46]+0.4564354645876384*fskin[44]+0.3535533905932737*fskin[20]+0.2041241452319315*fskin[12])*dv; 
  Ghat[12] = (1.58113883008419*fskin[55]+1.224744871391589*fskin[34]+0.7071067811865475*fskin[24])*wv+(0.408248290463863*fskin[32]+0.3162277660168379*fskin[15]+0.1825741858350554*fskin[7])*dv; 
  Ghat[13] = (1.58113883008419*fskin[57]+1.224744871391589*fskin[36]+0.7071067811865475*fskin[26])*wv+(0.4564354645876384*fskin[66]+0.3535533905932737*fskin[51]+0.2041241452319315*fskin[38])*dv; 
  Ghat[14] = (1.58113883008419*fskin[58]+1.224744871391589*fskin[39]+0.7071067811865475*fskin[27])*wv+(0.408248290463863*fskin[37]+0.3162277660168379*fskin[17]+0.1825741858350554*fskin[10])*dv; 
  Ghat[15] = (1.58113883008419*fskin[60]+1.224744871391589*fskin[41]+0.7071067811865475*fskin[29])*wv+(0.4564354645876384*fskin[69]+0.3535533905932737*fskin[53]+0.2041241452319315*fskin[43])*dv; 
  Ghat[16] = (1.58113883008419*fskin[62]+1.224744871391589*fskin[42]+0.7071067811865475*fskin[30])*wv+(0.408248290463863*fskin[74]+0.3162277660168379*fskin[64]+0.1825741858350554*fskin[49]+0.4564354645876384*fskin[47]+0.3535533905932737*fskin[28]+0.2041241452319315*fskin[14])*dv; 
  Ghat[17] = (1.58113883008419*fskin[66]+1.224744871391589*fskin[51]+0.7071067811865475*fskin[38])*wv+(0.408248290463863*fskin[76]+0.3162277660168379*fskin[68]+0.1825741858350554*fskin[59]+0.4564354645876384*fskin[57]+0.3535533905932737*fskin[36]+0.2041241452319315*fskin[26])*dv; 
  Ghat[18] = (1.58113883008419*fskin[67]+1.224744871391589*fskin[52]+0.7071067811865475*fskin[40])*wv+(0.408248290463863*fskin[50]+0.3162277660168379*fskin[31]+0.1825741858350554*fskin[18])*dv; 
  Ghat[19] = (1.58113883008419*fskin[69]+1.224744871391589*fskin[53]+0.7071067811865475*fskin[43])*wv+(0.408248290463863*fskin[78]+0.3162277660168379*fskin[71]+0.1825741858350554*fskin[65]+0.4564354645876384*fskin[60]+0.3535533905932737*fskin[41]+0.2041241452319315*fskin[29])*dv; 
  Ghat[20] = (1.58113883008419*fskin[72]+1.224744871391589*fskin[56]+0.7071067811865475*fskin[46])*wv+(0.408248290463863*fskin[54]+0.3162277660168379*fskin[33]+0.1825741858350554*fskin[22])*dv; 
  Ghat[21] = (1.58113883008419*fskin[73]+1.224744871391589*fskin[61]+0.7071067811865475*fskin[48])*wv+(0.4564354645876384*fskin[77]+0.3535533905932737*fskin[70]+0.2041241452319315*fskin[63])*dv; 
  Ghat[22] = (1.58113883008419*fskin[74]+1.224744871391589*fskin[64]+0.7071067811865475*fskin[49])*wv+(0.408248290463863*fskin[62]+0.3162277660168379*fskin[42]+0.1825741858350554*fskin[30])*dv; 
  Ghat[23] = (1.58113883008419*fskin[76]+1.224744871391589*fskin[68]+0.7071067811865475*fskin[59])*wv+(0.408248290463863*fskin[66]+0.3162277660168379*fskin[51]+0.1825741858350554*fskin[38])*dv; 
  Ghat[24] = (1.58113883008419*fskin[77]+1.224744871391589*fskin[70]+0.7071067811865475*fskin[63])*wv+(0.408248290463863*fskin[80]+0.3162277660168379*fskin[79]+0.1825741858350554*fskin[75]+0.4564354645876384*fskin[73]+0.3535533905932737*fskin[61]+0.2041241452319315*fskin[48])*dv; 
  Ghat[25] = (1.58113883008419*fskin[78]+1.224744871391589*fskin[71]+0.7071067811865475*fskin[65])*wv+(0.408248290463863*fskin[69]+0.3162277660168379*fskin[53]+0.1825741858350554*fskin[43])*dv; 
  Ghat[26] = (1.58113883008419*fskin[80]+1.224744871391589*fskin[79]+0.7071067811865475*fskin[75])*wv+(0.408248290463863*fskin[77]+0.3162277660168379*fskin[70]+0.1825741858350554*fskin[63])*dv; 

  } else { 

  Ghat[0] = (1.58113883008419*fedge[11]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*wv+(0.4564354645876384*fedge[21]-0.3535533905932737*fedge[6]+0.2041241452319315*fedge[3])*dv; 
  Ghat[1] = (1.58113883008419*fedge[19]-1.224744871391589*fedge[5]+0.7071067811865475*fedge[2])*wv+(0.4564354645876384*fedge[32]-0.3535533905932737*fedge[15]+0.2041241452319315*fedge[7])*dv; 
  Ghat[2] = (1.58113883008419*fedge[21]-1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])*wv+(0.408248290463863*fedge[45]-0.3162277660168379*fedge[23]+0.1825741858350554*fedge[13]+0.4564354645876384*fedge[11]-0.3535533905932737*fedge[1]+0.2041241452319315*fedge[0])*dv; 
  Ghat[3] = (1.58113883008419*fedge[25]-1.224744871391589*fedge[8]+0.7071067811865475*fedge[4])*wv+(0.4564354645876384*fedge[37]-0.3535533905932737*fedge[17]+0.2041241452319315*fedge[10])*dv; 
  Ghat[4] = (1.58113883008419*fedge[32]-1.224744871391589*fedge[15]+0.7071067811865475*fedge[7])*wv+(0.408248290463863*fedge[55]-0.3162277660168379*fedge[34]+0.1825741858350554*fedge[24]+0.4564354645876384*fedge[19]-0.3535533905932737*fedge[5]+0.2041241452319315*fedge[2])*dv; 
  Ghat[5] = (1.58113883008419*fedge[35]-1.224744871391589*fedge[16]+0.7071067811865475*fedge[9])*wv+(0.4564354645876384*fedge[50]-0.3535533905932737*fedge[31]+0.2041241452319315*fedge[18])*dv; 
  Ghat[6] = (1.58113883008419*fedge[37]-1.224744871391589*fedge[17]+0.7071067811865475*fedge[10])*wv+(0.408248290463863*fedge[58]-0.3162277660168379*fedge[39]+0.1825741858350554*fedge[27]+0.4564354645876384*fedge[25]-0.3535533905932737*fedge[8]+0.2041241452319315*fedge[4])*dv; 
  Ghat[7] = (1.58113883008419*fedge[44]-1.224744871391589*fedge[20]+0.7071067811865475*fedge[12])*wv+(0.4564354645876384*fedge[54]-0.3535533905932737*fedge[33]+0.2041241452319315*fedge[22])*dv; 
  Ghat[8] = (1.58113883008419*fedge[45]-1.224744871391589*fedge[23]+0.7071067811865475*fedge[13])*wv+(0.408248290463863*fedge[21]-0.3162277660168379*fedge[6]+0.1825741858350554*fedge[3])*dv; 
  Ghat[9] = (1.58113883008419*fedge[47]-1.224744871391589*fedge[28]+0.7071067811865475*fedge[14])*wv+(0.4564354645876384*fedge[62]-0.3535533905932737*fedge[42]+0.2041241452319315*fedge[30])*dv; 
  Ghat[10] = (1.58113883008419*fedge[50]-1.224744871391589*fedge[31]+0.7071067811865475*fedge[18])*wv+(0.408248290463863*fedge[67]-0.3162277660168379*fedge[52]+0.1825741858350554*fedge[40]+0.4564354645876384*fedge[35]-0.3535533905932737*fedge[16]+0.2041241452319315*fedge[9])*dv; 
  Ghat[11] = (1.58113883008419*fedge[54]-1.224744871391589*fedge[33]+0.7071067811865475*fedge[22])*wv+(0.408248290463863*fedge[72]-0.3162277660168379*fedge[56]+0.1825741858350554*fedge[46]+0.4564354645876384*fedge[44]-0.3535533905932737*fedge[20]+0.2041241452319315*fedge[12])*dv; 
  Ghat[12] = (1.58113883008419*fedge[55]-1.224744871391589*fedge[34]+0.7071067811865475*fedge[24])*wv+(0.408248290463863*fedge[32]-0.3162277660168379*fedge[15]+0.1825741858350554*fedge[7])*dv; 
  Ghat[13] = (1.58113883008419*fedge[57]-1.224744871391589*fedge[36]+0.7071067811865475*fedge[26])*wv+(0.4564354645876384*fedge[66]-0.3535533905932737*fedge[51]+0.2041241452319315*fedge[38])*dv; 
  Ghat[14] = (1.58113883008419*fedge[58]-1.224744871391589*fedge[39]+0.7071067811865475*fedge[27])*wv+(0.408248290463863*fedge[37]-0.3162277660168379*fedge[17]+0.1825741858350554*fedge[10])*dv; 
  Ghat[15] = (1.58113883008419*fedge[60]-1.224744871391589*fedge[41]+0.7071067811865475*fedge[29])*wv+(0.4564354645876384*fedge[69]-0.3535533905932737*fedge[53]+0.2041241452319315*fedge[43])*dv; 
  Ghat[16] = (1.58113883008419*fedge[62]-1.224744871391589*fedge[42]+0.7071067811865475*fedge[30])*wv+(0.408248290463863*fedge[74]-0.3162277660168379*fedge[64]+0.1825741858350554*fedge[49]+0.4564354645876384*fedge[47]-0.3535533905932737*fedge[28]+0.2041241452319315*fedge[14])*dv; 
  Ghat[17] = (1.58113883008419*fedge[66]-1.224744871391589*fedge[51]+0.7071067811865475*fedge[38])*wv+(0.408248290463863*fedge[76]-0.3162277660168379*fedge[68]+0.1825741858350554*fedge[59]+0.4564354645876384*fedge[57]-0.3535533905932737*fedge[36]+0.2041241452319315*fedge[26])*dv; 
  Ghat[18] = (1.58113883008419*fedge[67]-1.224744871391589*fedge[52]+0.7071067811865475*fedge[40])*wv+(0.408248290463863*fedge[50]-0.3162277660168379*fedge[31]+0.1825741858350554*fedge[18])*dv; 
  Ghat[19] = (1.58113883008419*fedge[69]-1.224744871391589*fedge[53]+0.7071067811865475*fedge[43])*wv+(0.408248290463863*fedge[78]-0.3162277660168379*fedge[71]+0.1825741858350554*fedge[65]+0.4564354645876384*fedge[60]-0.3535533905932737*fedge[41]+0.2041241452319315*fedge[29])*dv; 
  Ghat[20] = (1.58113883008419*fedge[72]-1.224744871391589*fedge[56]+0.7071067811865475*fedge[46])*wv+(0.408248290463863*fedge[54]-0.3162277660168379*fedge[33]+0.1825741858350554*fedge[22])*dv; 
  Ghat[21] = (1.58113883008419*fedge[73]-1.224744871391589*fedge[61]+0.7071067811865475*fedge[48])*wv+(0.4564354645876384*fedge[77]-0.3535533905932737*fedge[70]+0.2041241452319315*fedge[63])*dv; 
  Ghat[22] = (1.58113883008419*fedge[74]-1.224744871391589*fedge[64]+0.7071067811865475*fedge[49])*wv+(0.408248290463863*fedge[62]-0.3162277660168379*fedge[42]+0.1825741858350554*fedge[30])*dv; 
  Ghat[23] = (1.58113883008419*fedge[76]-1.224744871391589*fedge[68]+0.7071067811865475*fedge[59])*wv+(0.408248290463863*fedge[66]-0.3162277660168379*fedge[51]+0.1825741858350554*fedge[38])*dv; 
  Ghat[24] = (1.58113883008419*fedge[77]-1.224744871391589*fedge[70]+0.7071067811865475*fedge[63])*wv+(0.408248290463863*fedge[80]-0.3162277660168379*fedge[79]+0.1825741858350554*fedge[75]+0.4564354645876384*fedge[73]-0.3535533905932737*fedge[61]+0.2041241452319315*fedge[48])*dv; 
  Ghat[25] = (1.58113883008419*fedge[78]-1.224744871391589*fedge[71]+0.7071067811865475*fedge[65])*wv+(0.408248290463863*fedge[69]-0.3162277660168379*fedge[53]+0.1825741858350554*fedge[43])*dv; 
  Ghat[26] = (1.58113883008419*fedge[80]-1.224744871391589*fedge[79]+0.7071067811865475*fedge[75])*wv+(0.408248290463863*fedge[77]-0.3162277660168379*fedge[70]+0.1825741858350554*fedge[63])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -1.224744871391589*Ghat[1]*dx10; 
  out[6] += -1.224744871391589*Ghat[2]*dx10; 
  out[7] += -0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -1.224744871391589*Ghat[3]*dx10; 
  out[9] += -0.7071067811865475*Ghat[5]*dx10; 
  out[10] += -0.7071067811865475*Ghat[6]*dx10; 
  out[11] += -1.58113883008419*Ghat[0]*dx10; 
  out[12] += -0.7071067811865475*Ghat[7]*dx10; 
  out[13] += -0.7071067811865475*Ghat[8]*dx10; 
  out[14] += -0.7071067811865475*Ghat[9]*dx10; 
  out[15] += -1.224744871391589*Ghat[4]*dx10; 
  out[16] += -1.224744871391589*Ghat[5]*dx10; 
  out[17] += -1.224744871391589*Ghat[6]*dx10; 
  out[18] += -0.7071067811865475*Ghat[10]*dx10; 
  out[19] += -1.58113883008419*Ghat[1]*dx10; 
  out[20] += -1.224744871391589*Ghat[7]*dx10; 
  out[21] += -1.58113883008419*Ghat[2]*dx10; 
  out[22] += -0.7071067811865475*Ghat[11]*dx10; 
  out[23] += -1.224744871391589*Ghat[8]*dx10; 
  out[24] += -0.7071067811865475*Ghat[12]*dx10; 
  out[25] += -1.58113883008419*Ghat[3]*dx10; 
  out[26] += -0.7071067811865475*Ghat[13]*dx10; 
  out[27] += -0.7071067811865475*Ghat[14]*dx10; 
  out[28] += -1.224744871391589*Ghat[9]*dx10; 
  out[29] += -0.7071067811865475*Ghat[15]*dx10; 
  out[30] += -0.7071067811865475*Ghat[16]*dx10; 
  out[31] += -1.224744871391589*Ghat[10]*dx10; 
  out[32] += -1.58113883008419*Ghat[4]*dx10; 
  out[33] += -1.224744871391589*Ghat[11]*dx10; 
  out[34] += -1.224744871391589*Ghat[12]*dx10; 
  out[35] += -1.58113883008419*Ghat[5]*dx10; 
  out[36] += -1.224744871391589*Ghat[13]*dx10; 
  out[37] += -1.58113883008419*Ghat[6]*dx10; 
  out[38] += -0.7071067811865475*Ghat[17]*dx10; 
  out[39] += -1.224744871391589*Ghat[14]*dx10; 
  out[40] += -0.7071067811865475*Ghat[18]*dx10; 
  out[41] += -1.224744871391589*Ghat[15]*dx10; 
  out[42] += -1.224744871391589*Ghat[16]*dx10; 
  out[43] += -0.7071067811865475*Ghat[19]*dx10; 
  out[44] += -1.58113883008419*Ghat[7]*dx10; 
  out[45] += -1.58113883008419*Ghat[8]*dx10; 
  out[46] += -0.7071067811865475*Ghat[20]*dx10; 
  out[47] += -1.58113883008419*Ghat[9]*dx10; 
  out[48] += -0.7071067811865475*Ghat[21]*dx10; 
  out[49] += -0.7071067811865475*Ghat[22]*dx10; 
  out[50] += -1.58113883008419*Ghat[10]*dx10; 
  out[51] += -1.224744871391589*Ghat[17]*dx10; 
  out[52] += -1.224744871391589*Ghat[18]*dx10; 
  out[53] += -1.224744871391589*Ghat[19]*dx10; 
  out[54] += -1.58113883008419*Ghat[11]*dx10; 
  out[55] += -1.58113883008419*Ghat[12]*dx10; 
  out[56] += -1.224744871391589*Ghat[20]*dx10; 
  out[57] += -1.58113883008419*Ghat[13]*dx10; 
  out[58] += -1.58113883008419*Ghat[14]*dx10; 
  out[59] += -0.7071067811865475*Ghat[23]*dx10; 
  out[60] += -1.58113883008419*Ghat[15]*dx10; 
  out[61] += -1.224744871391589*Ghat[21]*dx10; 
  out[62] += -1.58113883008419*Ghat[16]*dx10; 
  out[63] += -0.7071067811865475*Ghat[24]*dx10; 
  out[64] += -1.224744871391589*Ghat[22]*dx10; 
  out[65] += -0.7071067811865475*Ghat[25]*dx10; 
  out[66] += -1.58113883008419*Ghat[17]*dx10; 
  out[67] += -1.58113883008419*Ghat[18]*dx10; 
  out[68] += -1.224744871391589*Ghat[23]*dx10; 
  out[69] += -1.58113883008419*Ghat[19]*dx10; 
  out[70] += -1.224744871391589*Ghat[24]*dx10; 
  out[71] += -1.224744871391589*Ghat[25]*dx10; 
  out[72] += -1.58113883008419*Ghat[20]*dx10; 
  out[73] += -1.58113883008419*Ghat[21]*dx10; 
  out[74] += -1.58113883008419*Ghat[22]*dx10; 
  out[75] += -0.7071067811865475*Ghat[26]*dx10; 
  out[76] += -1.58113883008419*Ghat[23]*dx10; 
  out[77] += -1.58113883008419*Ghat[24]*dx10; 
  out[78] += -1.58113883008419*Ghat[25]*dx10; 
  out[79] += -1.224744871391589*Ghat[26]*dx10; 
  out[80] += -1.58113883008419*Ghat[26]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fedge[11]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*wv+(0.4564354645876384*fedge[21]+0.3535533905932737*fedge[6]+0.2041241452319315*fedge[3])*dv; 
  Ghat[1] = (1.58113883008419*fedge[19]+1.224744871391589*fedge[5]+0.7071067811865475*fedge[2])*wv+(0.4564354645876384*fedge[32]+0.3535533905932737*fedge[15]+0.2041241452319315*fedge[7])*dv; 
  Ghat[2] = (1.58113883008419*fedge[21]+1.224744871391589*fedge[6]+0.7071067811865475*fedge[3])*wv+(0.408248290463863*fedge[45]+0.3162277660168379*fedge[23]+0.1825741858350554*fedge[13]+0.4564354645876384*fedge[11]+0.3535533905932737*fedge[1]+0.2041241452319315*fedge[0])*dv; 
  Ghat[3] = (1.58113883008419*fedge[25]+1.224744871391589*fedge[8]+0.7071067811865475*fedge[4])*wv+(0.4564354645876384*fedge[37]+0.3535533905932737*fedge[17]+0.2041241452319315*fedge[10])*dv; 
  Ghat[4] = (1.58113883008419*fedge[32]+1.224744871391589*fedge[15]+0.7071067811865475*fedge[7])*wv+(0.408248290463863*fedge[55]+0.3162277660168379*fedge[34]+0.1825741858350554*fedge[24]+0.4564354645876384*fedge[19]+0.3535533905932737*fedge[5]+0.2041241452319315*fedge[2])*dv; 
  Ghat[5] = (1.58113883008419*fedge[35]+1.224744871391589*fedge[16]+0.7071067811865475*fedge[9])*wv+(0.4564354645876384*fedge[50]+0.3535533905932737*fedge[31]+0.2041241452319315*fedge[18])*dv; 
  Ghat[6] = (1.58113883008419*fedge[37]+1.224744871391589*fedge[17]+0.7071067811865475*fedge[10])*wv+(0.408248290463863*fedge[58]+0.3162277660168379*fedge[39]+0.1825741858350554*fedge[27]+0.4564354645876384*fedge[25]+0.3535533905932737*fedge[8]+0.2041241452319315*fedge[4])*dv; 
  Ghat[7] = (1.58113883008419*fedge[44]+1.224744871391589*fedge[20]+0.7071067811865475*fedge[12])*wv+(0.4564354645876384*fedge[54]+0.3535533905932737*fedge[33]+0.2041241452319315*fedge[22])*dv; 
  Ghat[8] = (1.58113883008419*fedge[45]+1.224744871391589*fedge[23]+0.7071067811865475*fedge[13])*wv+(0.408248290463863*fedge[21]+0.3162277660168379*fedge[6]+0.1825741858350554*fedge[3])*dv; 
  Ghat[9] = (1.58113883008419*fedge[47]+1.224744871391589*fedge[28]+0.7071067811865475*fedge[14])*wv+(0.4564354645876384*fedge[62]+0.3535533905932737*fedge[42]+0.2041241452319315*fedge[30])*dv; 
  Ghat[10] = (1.58113883008419*fedge[50]+1.224744871391589*fedge[31]+0.7071067811865475*fedge[18])*wv+(0.408248290463863*fedge[67]+0.3162277660168379*fedge[52]+0.1825741858350554*fedge[40]+0.4564354645876384*fedge[35]+0.3535533905932737*fedge[16]+0.2041241452319315*fedge[9])*dv; 
  Ghat[11] = (1.58113883008419*fedge[54]+1.224744871391589*fedge[33]+0.7071067811865475*fedge[22])*wv+(0.408248290463863*fedge[72]+0.3162277660168379*fedge[56]+0.1825741858350554*fedge[46]+0.4564354645876384*fedge[44]+0.3535533905932737*fedge[20]+0.2041241452319315*fedge[12])*dv; 
  Ghat[12] = (1.58113883008419*fedge[55]+1.224744871391589*fedge[34]+0.7071067811865475*fedge[24])*wv+(0.408248290463863*fedge[32]+0.3162277660168379*fedge[15]+0.1825741858350554*fedge[7])*dv; 
  Ghat[13] = (1.58113883008419*fedge[57]+1.224744871391589*fedge[36]+0.7071067811865475*fedge[26])*wv+(0.4564354645876384*fedge[66]+0.3535533905932737*fedge[51]+0.2041241452319315*fedge[38])*dv; 
  Ghat[14] = (1.58113883008419*fedge[58]+1.224744871391589*fedge[39]+0.7071067811865475*fedge[27])*wv+(0.408248290463863*fedge[37]+0.3162277660168379*fedge[17]+0.1825741858350554*fedge[10])*dv; 
  Ghat[15] = (1.58113883008419*fedge[60]+1.224744871391589*fedge[41]+0.7071067811865475*fedge[29])*wv+(0.4564354645876384*fedge[69]+0.3535533905932737*fedge[53]+0.2041241452319315*fedge[43])*dv; 
  Ghat[16] = (1.58113883008419*fedge[62]+1.224744871391589*fedge[42]+0.7071067811865475*fedge[30])*wv+(0.408248290463863*fedge[74]+0.3162277660168379*fedge[64]+0.1825741858350554*fedge[49]+0.4564354645876384*fedge[47]+0.3535533905932737*fedge[28]+0.2041241452319315*fedge[14])*dv; 
  Ghat[17] = (1.58113883008419*fedge[66]+1.224744871391589*fedge[51]+0.7071067811865475*fedge[38])*wv+(0.408248290463863*fedge[76]+0.3162277660168379*fedge[68]+0.1825741858350554*fedge[59]+0.4564354645876384*fedge[57]+0.3535533905932737*fedge[36]+0.2041241452319315*fedge[26])*dv; 
  Ghat[18] = (1.58113883008419*fedge[67]+1.224744871391589*fedge[52]+0.7071067811865475*fedge[40])*wv+(0.408248290463863*fedge[50]+0.3162277660168379*fedge[31]+0.1825741858350554*fedge[18])*dv; 
  Ghat[19] = (1.58113883008419*fedge[69]+1.224744871391589*fedge[53]+0.7071067811865475*fedge[43])*wv+(0.408248290463863*fedge[78]+0.3162277660168379*fedge[71]+0.1825741858350554*fedge[65]+0.4564354645876384*fedge[60]+0.3535533905932737*fedge[41]+0.2041241452319315*fedge[29])*dv; 
  Ghat[20] = (1.58113883008419*fedge[72]+1.224744871391589*fedge[56]+0.7071067811865475*fedge[46])*wv+(0.408248290463863*fedge[54]+0.3162277660168379*fedge[33]+0.1825741858350554*fedge[22])*dv; 
  Ghat[21] = (1.58113883008419*fedge[73]+1.224744871391589*fedge[61]+0.7071067811865475*fedge[48])*wv+(0.4564354645876384*fedge[77]+0.3535533905932737*fedge[70]+0.2041241452319315*fedge[63])*dv; 
  Ghat[22] = (1.58113883008419*fedge[74]+1.224744871391589*fedge[64]+0.7071067811865475*fedge[49])*wv+(0.408248290463863*fedge[62]+0.3162277660168379*fedge[42]+0.1825741858350554*fedge[30])*dv; 
  Ghat[23] = (1.58113883008419*fedge[76]+1.224744871391589*fedge[68]+0.7071067811865475*fedge[59])*wv+(0.408248290463863*fedge[66]+0.3162277660168379*fedge[51]+0.1825741858350554*fedge[38])*dv; 
  Ghat[24] = (1.58113883008419*fedge[77]+1.224744871391589*fedge[70]+0.7071067811865475*fedge[63])*wv+(0.408248290463863*fedge[80]+0.3162277660168379*fedge[79]+0.1825741858350554*fedge[75]+0.4564354645876384*fedge[73]+0.3535533905932737*fedge[61]+0.2041241452319315*fedge[48])*dv; 
  Ghat[25] = (1.58113883008419*fedge[78]+1.224744871391589*fedge[71]+0.7071067811865475*fedge[65])*wv+(0.408248290463863*fedge[69]+0.3162277660168379*fedge[53]+0.1825741858350554*fedge[43])*dv; 
  Ghat[26] = (1.58113883008419*fedge[80]+1.224744871391589*fedge[79]+0.7071067811865475*fedge[75])*wv+(0.408248290463863*fedge[77]+0.3162277660168379*fedge[70]+0.1825741858350554*fedge[63])*dv; 

  } else { 

  Ghat[0] = (1.58113883008419*fskin[11]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*wv+(0.4564354645876384*fskin[21]-0.3535533905932737*fskin[6]+0.2041241452319315*fskin[3])*dv; 
  Ghat[1] = (1.58113883008419*fskin[19]-1.224744871391589*fskin[5]+0.7071067811865475*fskin[2])*wv+(0.4564354645876384*fskin[32]-0.3535533905932737*fskin[15]+0.2041241452319315*fskin[7])*dv; 
  Ghat[2] = (1.58113883008419*fskin[21]-1.224744871391589*fskin[6]+0.7071067811865475*fskin[3])*wv+(0.408248290463863*fskin[45]-0.3162277660168379*fskin[23]+0.1825741858350554*fskin[13]+0.4564354645876384*fskin[11]-0.3535533905932737*fskin[1]+0.2041241452319315*fskin[0])*dv; 
  Ghat[3] = (1.58113883008419*fskin[25]-1.224744871391589*fskin[8]+0.7071067811865475*fskin[4])*wv+(0.4564354645876384*fskin[37]-0.3535533905932737*fskin[17]+0.2041241452319315*fskin[10])*dv; 
  Ghat[4] = (1.58113883008419*fskin[32]-1.224744871391589*fskin[15]+0.7071067811865475*fskin[7])*wv+(0.408248290463863*fskin[55]-0.3162277660168379*fskin[34]+0.1825741858350554*fskin[24]+0.4564354645876384*fskin[19]-0.3535533905932737*fskin[5]+0.2041241452319315*fskin[2])*dv; 
  Ghat[5] = (1.58113883008419*fskin[35]-1.224744871391589*fskin[16]+0.7071067811865475*fskin[9])*wv+(0.4564354645876384*fskin[50]-0.3535533905932737*fskin[31]+0.2041241452319315*fskin[18])*dv; 
  Ghat[6] = (1.58113883008419*fskin[37]-1.224744871391589*fskin[17]+0.7071067811865475*fskin[10])*wv+(0.408248290463863*fskin[58]-0.3162277660168379*fskin[39]+0.1825741858350554*fskin[27]+0.4564354645876384*fskin[25]-0.3535533905932737*fskin[8]+0.2041241452319315*fskin[4])*dv; 
  Ghat[7] = (1.58113883008419*fskin[44]-1.224744871391589*fskin[20]+0.7071067811865475*fskin[12])*wv+(0.4564354645876384*fskin[54]-0.3535533905932737*fskin[33]+0.2041241452319315*fskin[22])*dv; 
  Ghat[8] = (1.58113883008419*fskin[45]-1.224744871391589*fskin[23]+0.7071067811865475*fskin[13])*wv+(0.408248290463863*fskin[21]-0.3162277660168379*fskin[6]+0.1825741858350554*fskin[3])*dv; 
  Ghat[9] = (1.58113883008419*fskin[47]-1.224744871391589*fskin[28]+0.7071067811865475*fskin[14])*wv+(0.4564354645876384*fskin[62]-0.3535533905932737*fskin[42]+0.2041241452319315*fskin[30])*dv; 
  Ghat[10] = (1.58113883008419*fskin[50]-1.224744871391589*fskin[31]+0.7071067811865475*fskin[18])*wv+(0.408248290463863*fskin[67]-0.3162277660168379*fskin[52]+0.1825741858350554*fskin[40]+0.4564354645876384*fskin[35]-0.3535533905932737*fskin[16]+0.2041241452319315*fskin[9])*dv; 
  Ghat[11] = (1.58113883008419*fskin[54]-1.224744871391589*fskin[33]+0.7071067811865475*fskin[22])*wv+(0.408248290463863*fskin[72]-0.3162277660168379*fskin[56]+0.1825741858350554*fskin[46]+0.4564354645876384*fskin[44]-0.3535533905932737*fskin[20]+0.2041241452319315*fskin[12])*dv; 
  Ghat[12] = (1.58113883008419*fskin[55]-1.224744871391589*fskin[34]+0.7071067811865475*fskin[24])*wv+(0.408248290463863*fskin[32]-0.3162277660168379*fskin[15]+0.1825741858350554*fskin[7])*dv; 
  Ghat[13] = (1.58113883008419*fskin[57]-1.224744871391589*fskin[36]+0.7071067811865475*fskin[26])*wv+(0.4564354645876384*fskin[66]-0.3535533905932737*fskin[51]+0.2041241452319315*fskin[38])*dv; 
  Ghat[14] = (1.58113883008419*fskin[58]-1.224744871391589*fskin[39]+0.7071067811865475*fskin[27])*wv+(0.408248290463863*fskin[37]-0.3162277660168379*fskin[17]+0.1825741858350554*fskin[10])*dv; 
  Ghat[15] = (1.58113883008419*fskin[60]-1.224744871391589*fskin[41]+0.7071067811865475*fskin[29])*wv+(0.4564354645876384*fskin[69]-0.3535533905932737*fskin[53]+0.2041241452319315*fskin[43])*dv; 
  Ghat[16] = (1.58113883008419*fskin[62]-1.224744871391589*fskin[42]+0.7071067811865475*fskin[30])*wv+(0.408248290463863*fskin[74]-0.3162277660168379*fskin[64]+0.1825741858350554*fskin[49]+0.4564354645876384*fskin[47]-0.3535533905932737*fskin[28]+0.2041241452319315*fskin[14])*dv; 
  Ghat[17] = (1.58113883008419*fskin[66]-1.224744871391589*fskin[51]+0.7071067811865475*fskin[38])*wv+(0.408248290463863*fskin[76]-0.3162277660168379*fskin[68]+0.1825741858350554*fskin[59]+0.4564354645876384*fskin[57]-0.3535533905932737*fskin[36]+0.2041241452319315*fskin[26])*dv; 
  Ghat[18] = (1.58113883008419*fskin[67]-1.224744871391589*fskin[52]+0.7071067811865475*fskin[40])*wv+(0.408248290463863*fskin[50]-0.3162277660168379*fskin[31]+0.1825741858350554*fskin[18])*dv; 
  Ghat[19] = (1.58113883008419*fskin[69]-1.224744871391589*fskin[53]+0.7071067811865475*fskin[43])*wv+(0.408248290463863*fskin[78]-0.3162277660168379*fskin[71]+0.1825741858350554*fskin[65]+0.4564354645876384*fskin[60]-0.3535533905932737*fskin[41]+0.2041241452319315*fskin[29])*dv; 
  Ghat[20] = (1.58113883008419*fskin[72]-1.224744871391589*fskin[56]+0.7071067811865475*fskin[46])*wv+(0.408248290463863*fskin[54]-0.3162277660168379*fskin[33]+0.1825741858350554*fskin[22])*dv; 
  Ghat[21] = (1.58113883008419*fskin[73]-1.224744871391589*fskin[61]+0.7071067811865475*fskin[48])*wv+(0.4564354645876384*fskin[77]-0.3535533905932737*fskin[70]+0.2041241452319315*fskin[63])*dv; 
  Ghat[22] = (1.58113883008419*fskin[74]-1.224744871391589*fskin[64]+0.7071067811865475*fskin[49])*wv+(0.408248290463863*fskin[62]-0.3162277660168379*fskin[42]+0.1825741858350554*fskin[30])*dv; 
  Ghat[23] = (1.58113883008419*fskin[76]-1.224744871391589*fskin[68]+0.7071067811865475*fskin[59])*wv+(0.408248290463863*fskin[66]-0.3162277660168379*fskin[51]+0.1825741858350554*fskin[38])*dv; 
  Ghat[24] = (1.58113883008419*fskin[77]-1.224744871391589*fskin[70]+0.7071067811865475*fskin[63])*wv+(0.408248290463863*fskin[80]-0.3162277660168379*fskin[79]+0.1825741858350554*fskin[75]+0.4564354645876384*fskin[73]-0.3535533905932737*fskin[61]+0.2041241452319315*fskin[48])*dv; 
  Ghat[25] = (1.58113883008419*fskin[78]-1.224744871391589*fskin[71]+0.7071067811865475*fskin[65])*wv+(0.408248290463863*fskin[69]-0.3162277660168379*fskin[53]+0.1825741858350554*fskin[43])*dv; 
  Ghat[26] = (1.58113883008419*fskin[80]-1.224744871391589*fskin[79]+0.7071067811865475*fskin[75])*wv+(0.408248290463863*fskin[77]-0.3162277660168379*fskin[70]+0.1825741858350554*fskin[63])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dx10; 
  out[5] += -1.224744871391589*Ghat[1]*dx10; 
  out[6] += -1.224744871391589*Ghat[2]*dx10; 
  out[7] += 0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -1.224744871391589*Ghat[3]*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dx10; 
  out[11] += 1.58113883008419*Ghat[0]*dx10; 
  out[12] += 0.7071067811865475*Ghat[7]*dx10; 
  out[13] += 0.7071067811865475*Ghat[8]*dx10; 
  out[14] += 0.7071067811865475*Ghat[9]*dx10; 
  out[15] += -1.224744871391589*Ghat[4]*dx10; 
  out[16] += -1.224744871391589*Ghat[5]*dx10; 
  out[17] += -1.224744871391589*Ghat[6]*dx10; 
  out[18] += 0.7071067811865475*Ghat[10]*dx10; 
  out[19] += 1.58113883008419*Ghat[1]*dx10; 
  out[20] += -1.224744871391589*Ghat[7]*dx10; 
  out[21] += 1.58113883008419*Ghat[2]*dx10; 
  out[22] += 0.7071067811865475*Ghat[11]*dx10; 
  out[23] += -1.224744871391589*Ghat[8]*dx10; 
  out[24] += 0.7071067811865475*Ghat[12]*dx10; 
  out[25] += 1.58113883008419*Ghat[3]*dx10; 
  out[26] += 0.7071067811865475*Ghat[13]*dx10; 
  out[27] += 0.7071067811865475*Ghat[14]*dx10; 
  out[28] += -1.224744871391589*Ghat[9]*dx10; 
  out[29] += 0.7071067811865475*Ghat[15]*dx10; 
  out[30] += 0.7071067811865475*Ghat[16]*dx10; 
  out[31] += -1.224744871391589*Ghat[10]*dx10; 
  out[32] += 1.58113883008419*Ghat[4]*dx10; 
  out[33] += -1.224744871391589*Ghat[11]*dx10; 
  out[34] += -1.224744871391589*Ghat[12]*dx10; 
  out[35] += 1.58113883008419*Ghat[5]*dx10; 
  out[36] += -1.224744871391589*Ghat[13]*dx10; 
  out[37] += 1.58113883008419*Ghat[6]*dx10; 
  out[38] += 0.7071067811865475*Ghat[17]*dx10; 
  out[39] += -1.224744871391589*Ghat[14]*dx10; 
  out[40] += 0.7071067811865475*Ghat[18]*dx10; 
  out[41] += -1.224744871391589*Ghat[15]*dx10; 
  out[42] += -1.224744871391589*Ghat[16]*dx10; 
  out[43] += 0.7071067811865475*Ghat[19]*dx10; 
  out[44] += 1.58113883008419*Ghat[7]*dx10; 
  out[45] += 1.58113883008419*Ghat[8]*dx10; 
  out[46] += 0.7071067811865475*Ghat[20]*dx10; 
  out[47] += 1.58113883008419*Ghat[9]*dx10; 
  out[48] += 0.7071067811865475*Ghat[21]*dx10; 
  out[49] += 0.7071067811865475*Ghat[22]*dx10; 
  out[50] += 1.58113883008419*Ghat[10]*dx10; 
  out[51] += -1.224744871391589*Ghat[17]*dx10; 
  out[52] += -1.224744871391589*Ghat[18]*dx10; 
  out[53] += -1.224744871391589*Ghat[19]*dx10; 
  out[54] += 1.58113883008419*Ghat[11]*dx10; 
  out[55] += 1.58113883008419*Ghat[12]*dx10; 
  out[56] += -1.224744871391589*Ghat[20]*dx10; 
  out[57] += 1.58113883008419*Ghat[13]*dx10; 
  out[58] += 1.58113883008419*Ghat[14]*dx10; 
  out[59] += 0.7071067811865475*Ghat[23]*dx10; 
  out[60] += 1.58113883008419*Ghat[15]*dx10; 
  out[61] += -1.224744871391589*Ghat[21]*dx10; 
  out[62] += 1.58113883008419*Ghat[16]*dx10; 
  out[63] += 0.7071067811865475*Ghat[24]*dx10; 
  out[64] += -1.224744871391589*Ghat[22]*dx10; 
  out[65] += 0.7071067811865475*Ghat[25]*dx10; 
  out[66] += 1.58113883008419*Ghat[17]*dx10; 
  out[67] += 1.58113883008419*Ghat[18]*dx10; 
  out[68] += -1.224744871391589*Ghat[23]*dx10; 
  out[69] += 1.58113883008419*Ghat[19]*dx10; 
  out[70] += -1.224744871391589*Ghat[24]*dx10; 
  out[71] += -1.224744871391589*Ghat[25]*dx10; 
  out[72] += 1.58113883008419*Ghat[20]*dx10; 
  out[73] += 1.58113883008419*Ghat[21]*dx10; 
  out[74] += 1.58113883008419*Ghat[22]*dx10; 
  out[75] += 0.7071067811865475*Ghat[26]*dx10; 
  out[76] += 1.58113883008419*Ghat[23]*dx10; 
  out[77] += 1.58113883008419*Ghat[24]*dx10; 
  out[78] += 1.58113883008419*Ghat[25]*dx10; 
  out[79] += -1.224744871391589*Ghat[26]*dx10; 
  out[80] += 1.58113883008419*Ghat[26]*dx10; 

  } 
  return 0.;

} 