#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_tensor_p2(const double *w, const double *dxv, 
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
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[9]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fskin[7]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*wv+(0.4564354645876384*fskin[11]+0.3535533905932737*fskin[4]+0.2041241452319315*fskin[2])*dv; 
  Ghat[1] = (1.58113883008419*fskin[11]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[2])*wv+(0.408248290463863*fskin[20]+0.3162277660168379*fskin[12]+0.1825741858350554*fskin[8]+0.4564354645876384*fskin[7]+0.3535533905932737*fskin[1]+0.2041241452319315*fskin[0])*dv; 
  Ghat[2] = (1.58113883008419*fskin[13]+1.224744871391589*fskin[5]+0.7071067811865475*fskin[3])*wv+(0.4564354645876384*fskin[17]+0.3535533905932737*fskin[10]+0.2041241452319315*fskin[6])*dv; 
  Ghat[3] = (1.58113883008419*fskin[17]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[6])*wv+(0.408248290463863*fskin[23]+0.3162277660168379*fskin[18]+0.1825741858350554*fskin[14]+0.4564354645876384*fskin[13]+0.3535533905932737*fskin[5]+0.2041241452319315*fskin[3])*dv; 
  Ghat[4] = (1.58113883008419*fskin[20]+1.224744871391589*fskin[12]+0.7071067811865475*fskin[8])*wv+(0.408248290463863*fskin[11]+0.3162277660168379*fskin[4]+0.1825741858350554*fskin[2])*dv; 
  Ghat[5] = (1.58113883008419*fskin[21]+1.224744871391589*fskin[15]+0.7071067811865475*fskin[9])*wv+(0.4564354645876384*fskin[24]+0.3535533905932737*fskin[19]+0.2041241452319315*fskin[16])*dv; 
  Ghat[6] = (1.58113883008419*fskin[23]+1.224744871391589*fskin[18]+0.7071067811865475*fskin[14])*wv+(0.408248290463863*fskin[17]+0.3162277660168379*fskin[10]+0.1825741858350554*fskin[6])*dv; 
  Ghat[7] = (1.58113883008419*fskin[24]+1.224744871391589*fskin[19]+0.7071067811865475*fskin[16])*wv+(0.408248290463863*fskin[26]+0.3162277660168379*fskin[25]+0.1825741858350554*fskin[22]+0.4564354645876384*fskin[21]+0.3535533905932737*fskin[15]+0.2041241452319315*fskin[9])*dv; 
  Ghat[8] = (1.58113883008419*fskin[26]+1.224744871391589*fskin[25]+0.7071067811865475*fskin[22])*wv+(0.408248290463863*fskin[24]+0.3162277660168379*fskin[19]+0.1825741858350554*fskin[16])*dv; 

  } else { 

  Ghat[0] = (1.58113883008419*fedge[7]-1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*wv+(0.4564354645876384*fedge[11]-0.3535533905932737*fedge[4]+0.2041241452319315*fedge[2])*dv; 
  Ghat[1] = (1.58113883008419*fedge[11]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[2])*wv+(0.408248290463863*fedge[20]-0.3162277660168379*fedge[12]+0.1825741858350554*fedge[8]+0.4564354645876384*fedge[7]-0.3535533905932737*fedge[1]+0.2041241452319315*fedge[0])*dv; 
  Ghat[2] = (1.58113883008419*fedge[13]-1.224744871391589*fedge[5]+0.7071067811865475*fedge[3])*wv+(0.4564354645876384*fedge[17]-0.3535533905932737*fedge[10]+0.2041241452319315*fedge[6])*dv; 
  Ghat[3] = (1.58113883008419*fedge[17]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[6])*wv+(0.408248290463863*fedge[23]-0.3162277660168379*fedge[18]+0.1825741858350554*fedge[14]+0.4564354645876384*fedge[13]-0.3535533905932737*fedge[5]+0.2041241452319315*fedge[3])*dv; 
  Ghat[4] = (1.58113883008419*fedge[20]-1.224744871391589*fedge[12]+0.7071067811865475*fedge[8])*wv+(0.408248290463863*fedge[11]-0.3162277660168379*fedge[4]+0.1825741858350554*fedge[2])*dv; 
  Ghat[5] = (1.58113883008419*fedge[21]-1.224744871391589*fedge[15]+0.7071067811865475*fedge[9])*wv+(0.4564354645876384*fedge[24]-0.3535533905932737*fedge[19]+0.2041241452319315*fedge[16])*dv; 
  Ghat[6] = (1.58113883008419*fedge[23]-1.224744871391589*fedge[18]+0.7071067811865475*fedge[14])*wv+(0.408248290463863*fedge[17]-0.3162277660168379*fedge[10]+0.1825741858350554*fedge[6])*dv; 
  Ghat[7] = (1.58113883008419*fedge[24]-1.224744871391589*fedge[19]+0.7071067811865475*fedge[16])*wv+(0.408248290463863*fedge[26]-0.3162277660168379*fedge[25]+0.1825741858350554*fedge[22]+0.4564354645876384*fedge[21]-0.3535533905932737*fedge[15]+0.2041241452319315*fedge[9])*dv; 
  Ghat[8] = (1.58113883008419*fedge[26]-1.224744871391589*fedge[25]+0.7071067811865475*fedge[22])*wv+(0.408248290463863*fedge[24]-0.3162277660168379*fedge[19]+0.1825741858350554*fedge[16])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += -0.7071067811865475*Ghat[3]*dx10; 
  out[7] += -1.58113883008419*Ghat[0]*dx10; 
  out[8] += -0.7071067811865475*Ghat[4]*dx10; 
  out[9] += -0.7071067811865475*Ghat[5]*dx10; 
  out[10] += -1.224744871391589*Ghat[3]*dx10; 
  out[11] += -1.58113883008419*Ghat[1]*dx10; 
  out[12] += -1.224744871391589*Ghat[4]*dx10; 
  out[13] += -1.58113883008419*Ghat[2]*dx10; 
  out[14] += -0.7071067811865475*Ghat[6]*dx10; 
  out[15] += -1.224744871391589*Ghat[5]*dx10; 
  out[16] += -0.7071067811865475*Ghat[7]*dx10; 
  out[17] += -1.58113883008419*Ghat[3]*dx10; 
  out[18] += -1.224744871391589*Ghat[6]*dx10; 
  out[19] += -1.224744871391589*Ghat[7]*dx10; 
  out[20] += -1.58113883008419*Ghat[4]*dx10; 
  out[21] += -1.58113883008419*Ghat[5]*dx10; 
  out[22] += -0.7071067811865475*Ghat[8]*dx10; 
  out[23] += -1.58113883008419*Ghat[6]*dx10; 
  out[24] += -1.58113883008419*Ghat[7]*dx10; 
  out[25] += -1.224744871391589*Ghat[8]*dx10; 
  out[26] += -1.58113883008419*Ghat[8]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fedge[7]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*wv+(0.4564354645876384*fedge[11]+0.3535533905932737*fedge[4]+0.2041241452319315*fedge[2])*dv; 
  Ghat[1] = (1.58113883008419*fedge[11]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[2])*wv+(0.408248290463863*fedge[20]+0.3162277660168379*fedge[12]+0.1825741858350554*fedge[8]+0.4564354645876384*fedge[7]+0.3535533905932737*fedge[1]+0.2041241452319315*fedge[0])*dv; 
  Ghat[2] = (1.58113883008419*fedge[13]+1.224744871391589*fedge[5]+0.7071067811865475*fedge[3])*wv+(0.4564354645876384*fedge[17]+0.3535533905932737*fedge[10]+0.2041241452319315*fedge[6])*dv; 
  Ghat[3] = (1.58113883008419*fedge[17]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[6])*wv+(0.408248290463863*fedge[23]+0.3162277660168379*fedge[18]+0.1825741858350554*fedge[14]+0.4564354645876384*fedge[13]+0.3535533905932737*fedge[5]+0.2041241452319315*fedge[3])*dv; 
  Ghat[4] = (1.58113883008419*fedge[20]+1.224744871391589*fedge[12]+0.7071067811865475*fedge[8])*wv+(0.408248290463863*fedge[11]+0.3162277660168379*fedge[4]+0.1825741858350554*fedge[2])*dv; 
  Ghat[5] = (1.58113883008419*fedge[21]+1.224744871391589*fedge[15]+0.7071067811865475*fedge[9])*wv+(0.4564354645876384*fedge[24]+0.3535533905932737*fedge[19]+0.2041241452319315*fedge[16])*dv; 
  Ghat[6] = (1.58113883008419*fedge[23]+1.224744871391589*fedge[18]+0.7071067811865475*fedge[14])*wv+(0.408248290463863*fedge[17]+0.3162277660168379*fedge[10]+0.1825741858350554*fedge[6])*dv; 
  Ghat[7] = (1.58113883008419*fedge[24]+1.224744871391589*fedge[19]+0.7071067811865475*fedge[16])*wv+(0.408248290463863*fedge[26]+0.3162277660168379*fedge[25]+0.1825741858350554*fedge[22]+0.4564354645876384*fedge[21]+0.3535533905932737*fedge[15]+0.2041241452319315*fedge[9])*dv; 
  Ghat[8] = (1.58113883008419*fedge[26]+1.224744871391589*fedge[25]+0.7071067811865475*fedge[22])*wv+(0.408248290463863*fedge[24]+0.3162277660168379*fedge[19]+0.1825741858350554*fedge[16])*dv; 

  } else { 

  Ghat[0] = (1.58113883008419*fskin[7]-1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*wv+(0.4564354645876384*fskin[11]-0.3535533905932737*fskin[4]+0.2041241452319315*fskin[2])*dv; 
  Ghat[1] = (1.58113883008419*fskin[11]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[2])*wv+(0.408248290463863*fskin[20]-0.3162277660168379*fskin[12]+0.1825741858350554*fskin[8]+0.4564354645876384*fskin[7]-0.3535533905932737*fskin[1]+0.2041241452319315*fskin[0])*dv; 
  Ghat[2] = (1.58113883008419*fskin[13]-1.224744871391589*fskin[5]+0.7071067811865475*fskin[3])*wv+(0.4564354645876384*fskin[17]-0.3535533905932737*fskin[10]+0.2041241452319315*fskin[6])*dv; 
  Ghat[3] = (1.58113883008419*fskin[17]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[6])*wv+(0.408248290463863*fskin[23]-0.3162277660168379*fskin[18]+0.1825741858350554*fskin[14]+0.4564354645876384*fskin[13]-0.3535533905932737*fskin[5]+0.2041241452319315*fskin[3])*dv; 
  Ghat[4] = (1.58113883008419*fskin[20]-1.224744871391589*fskin[12]+0.7071067811865475*fskin[8])*wv+(0.408248290463863*fskin[11]-0.3162277660168379*fskin[4]+0.1825741858350554*fskin[2])*dv; 
  Ghat[5] = (1.58113883008419*fskin[21]-1.224744871391589*fskin[15]+0.7071067811865475*fskin[9])*wv+(0.4564354645876384*fskin[24]-0.3535533905932737*fskin[19]+0.2041241452319315*fskin[16])*dv; 
  Ghat[6] = (1.58113883008419*fskin[23]-1.224744871391589*fskin[18]+0.7071067811865475*fskin[14])*wv+(0.408248290463863*fskin[17]-0.3162277660168379*fskin[10]+0.1825741858350554*fskin[6])*dv; 
  Ghat[7] = (1.58113883008419*fskin[24]-1.224744871391589*fskin[19]+0.7071067811865475*fskin[16])*wv+(0.408248290463863*fskin[26]-0.3162277660168379*fskin[25]+0.1825741858350554*fskin[22]+0.4564354645876384*fskin[21]-0.3535533905932737*fskin[15]+0.2041241452319315*fskin[9])*dv; 
  Ghat[8] = (1.58113883008419*fskin[26]-1.224744871391589*fskin[25]+0.7071067811865475*fskin[22])*wv+(0.408248290463863*fskin[24]-0.3162277660168379*fskin[19]+0.1825741858350554*fskin[16])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += 0.7071067811865475*Ghat[3]*dx10; 
  out[7] += 1.58113883008419*Ghat[0]*dx10; 
  out[8] += 0.7071067811865475*Ghat[4]*dx10; 
  out[9] += 0.7071067811865475*Ghat[5]*dx10; 
  out[10] += -1.224744871391589*Ghat[3]*dx10; 
  out[11] += 1.58113883008419*Ghat[1]*dx10; 
  out[12] += -1.224744871391589*Ghat[4]*dx10; 
  out[13] += 1.58113883008419*Ghat[2]*dx10; 
  out[14] += 0.7071067811865475*Ghat[6]*dx10; 
  out[15] += -1.224744871391589*Ghat[5]*dx10; 
  out[16] += 0.7071067811865475*Ghat[7]*dx10; 
  out[17] += 1.58113883008419*Ghat[3]*dx10; 
  out[18] += -1.224744871391589*Ghat[6]*dx10; 
  out[19] += -1.224744871391589*Ghat[7]*dx10; 
  out[20] += 1.58113883008419*Ghat[4]*dx10; 
  out[21] += 1.58113883008419*Ghat[5]*dx10; 
  out[22] += 0.7071067811865475*Ghat[8]*dx10; 
  out[23] += 1.58113883008419*Ghat[6]*dx10; 
  out[24] += 1.58113883008419*Ghat[7]*dx10; 
  out[25] += -1.224744871391589*Ghat[8]*dx10; 
  out[26] += 1.58113883008419*Ghat[8]*dx10; 

  } 
  return 0.;

} 