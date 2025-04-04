#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv,
    const double *vmap_prime_edge, const double *vmap_prime_skin,
    const double *alpha_surf_edge, const double *alpha_surf_skin, 
    const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
    const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // alpha_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // alpha_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // sgn_alpha_surf_edge: sign(alpha_surf_edge) at quadrature points.
  // sgn_alpha_surf_skin: sign(alpha_surf_skin) at quadrature points.
  // const_sgn_alpha_edge: Boolean array true if sign(alpha_surf_edge) is only one sign, either +1 or -1.
  // const_sgn_alpha_skin: Boolean array true if sign(alpha_surf_skin) is only one sign, either +1 or -1.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double rdvpar2 = 2.0/dxv[3];

  const double *alphaL = &alpha_surf_skin[72];
  const double *alphaR = &alpha_surf_edge[72];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[72];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[72];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[3];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[3];

  if (edge == -1) { 

  double fUpR[16] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = (1.58113883008419*fskin[32]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  fUpR[1] = (1.58113883008419*fskin[33]+1.224744871391589*fskin[9]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  fUpR[2] = (1.58113883008419*fskin[34]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[2])/vmap_prime_skin[0]; 
  fUpR[3] = (1.58113883008419*fskin[35]+1.224744871391589*fskin[11]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  fUpR[4] = (1.58113883008419*fskin[36]+1.224744871391589*fskin[15]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 
  fUpR[5] = (1.58113883008419*fskin[37]+1.224744871391589*fskin[17]+0.7071067811865475*fskin[6])/vmap_prime_skin[0]; 
  fUpR[6] = (1.58113883008419*fskin[38]+1.224744871391589*fskin[18]+0.7071067811865475*fskin[7])/vmap_prime_skin[0]; 
  fUpR[7] = (1.58113883008419*fskin[39]+1.224744871391589*fskin[19]+0.7071067811865475*fskin[8])/vmap_prime_skin[0]; 
  fUpR[8] = (1.58113883008419*fskin[40]+1.224744871391589*fskin[23]+0.7071067811865475*fskin[12])/vmap_prime_skin[0]; 
  fUpR[9] = (1.58113883008419*fskin[41]+1.224744871391589*fskin[24]+0.7071067811865475*fskin[13])/vmap_prime_skin[0]; 
  fUpR[10] = (1.58113883008419*fskin[42]+1.224744871391589*fskin[25]+0.7071067811865475*fskin[14])/vmap_prime_skin[0]; 
  fUpR[11] = (1.58113883008419*fskin[43]+1.224744871391589*fskin[26]+0.7071067811865475*fskin[16])/vmap_prime_skin[0]; 
  fUpR[12] = (1.58113883008419*fskin[44]+1.224744871391589*fskin[28]+0.7071067811865475*fskin[20])/vmap_prime_skin[0]; 
  fUpR[13] = (1.58113883008419*fskin[45]+1.224744871391589*fskin[29]+0.7071067811865475*fskin[21])/vmap_prime_skin[0]; 
  fUpR[14] = (1.58113883008419*fskin[46]+1.224744871391589*fskin[30]+0.7071067811865475*fskin[22])/vmap_prime_skin[0]; 
  fUpR[15] = (1.58113883008419*fskin[47]+1.224744871391589*fskin[31]+0.7071067811865475*fskin[27])/vmap_prime_skin[0]; 
    } else { 
  fUpR[0] = (1.58113883008419*fedge[32]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  fUpR[1] = (1.58113883008419*fedge[33]-1.224744871391589*fedge[9]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  fUpR[2] = (1.58113883008419*fedge[34]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[2])/vmap_prime_edge[0]; 
  fUpR[3] = (1.58113883008419*fedge[35]-1.224744871391589*fedge[11]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  fUpR[4] = (1.58113883008419*fedge[36]-1.224744871391589*fedge[15]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 
  fUpR[5] = (1.58113883008419*fedge[37]-1.224744871391589*fedge[17]+0.7071067811865475*fedge[6])/vmap_prime_edge[0]; 
  fUpR[6] = (1.58113883008419*fedge[38]-1.224744871391589*fedge[18]+0.7071067811865475*fedge[7])/vmap_prime_edge[0]; 
  fUpR[7] = (1.58113883008419*fedge[39]-1.224744871391589*fedge[19]+0.7071067811865475*fedge[8])/vmap_prime_edge[0]; 
  fUpR[8] = (1.58113883008419*fedge[40]-1.224744871391589*fedge[23]+0.7071067811865475*fedge[12])/vmap_prime_edge[0]; 
  fUpR[9] = (1.58113883008419*fedge[41]-1.224744871391589*fedge[24]+0.7071067811865475*fedge[13])/vmap_prime_edge[0]; 
  fUpR[10] = (1.58113883008419*fedge[42]-1.224744871391589*fedge[25]+0.7071067811865475*fedge[14])/vmap_prime_edge[0]; 
  fUpR[11] = (1.58113883008419*fedge[43]-1.224744871391589*fedge[26]+0.7071067811865475*fedge[16])/vmap_prime_edge[0]; 
  fUpR[12] = (1.58113883008419*fedge[44]-1.224744871391589*fedge[28]+0.7071067811865475*fedge[20])/vmap_prime_edge[0]; 
  fUpR[13] = (1.58113883008419*fedge[45]-1.224744871391589*fedge[29]+0.7071067811865475*fedge[21])/vmap_prime_edge[0]; 
  fUpR[14] = (1.58113883008419*fedge[46]-1.224744871391589*fedge[30]+0.7071067811865475*fedge[22])/vmap_prime_edge[0]; 
  fUpR[15] = (1.58113883008419*fedge[47]-1.224744871391589*fedge[31]+0.7071067811865475*fedge[27])/vmap_prime_edge[0]; 
    } 
  } else { 
  double f_cr[16] = {0.};
  double f_rl[16] = {0.};
  double sgn_alphaUpR[16] = {0.};
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = (1.58113883008419*fskin[32]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  f_cr[1] = (1.58113883008419*fskin[33]+1.224744871391589*fskin[9]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  f_cr[2] = (1.58113883008419*fskin[34]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[2])/vmap_prime_skin[0]; 
  f_cr[3] = (1.58113883008419*fskin[35]+1.224744871391589*fskin[11]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  f_cr[4] = (1.58113883008419*fskin[36]+1.224744871391589*fskin[15]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 
  f_cr[5] = (1.58113883008419*fskin[37]+1.224744871391589*fskin[17]+0.7071067811865475*fskin[6])/vmap_prime_skin[0]; 
  f_cr[6] = (1.58113883008419*fskin[38]+1.224744871391589*fskin[18]+0.7071067811865475*fskin[7])/vmap_prime_skin[0]; 
  f_cr[7] = (1.58113883008419*fskin[39]+1.224744871391589*fskin[19]+0.7071067811865475*fskin[8])/vmap_prime_skin[0]; 
  f_cr[8] = (1.58113883008419*fskin[40]+1.224744871391589*fskin[23]+0.7071067811865475*fskin[12])/vmap_prime_skin[0]; 
  f_cr[9] = (1.58113883008419*fskin[41]+1.224744871391589*fskin[24]+0.7071067811865475*fskin[13])/vmap_prime_skin[0]; 
  f_cr[10] = (1.58113883008419*fskin[42]+1.224744871391589*fskin[25]+0.7071067811865475*fskin[14])/vmap_prime_skin[0]; 
  f_cr[11] = (1.58113883008419*fskin[43]+1.224744871391589*fskin[26]+0.7071067811865475*fskin[16])/vmap_prime_skin[0]; 
  f_cr[12] = (1.58113883008419*fskin[44]+1.224744871391589*fskin[28]+0.7071067811865475*fskin[20])/vmap_prime_skin[0]; 
  f_cr[13] = (1.58113883008419*fskin[45]+1.224744871391589*fskin[29]+0.7071067811865475*fskin[21])/vmap_prime_skin[0]; 
  f_cr[14] = (1.58113883008419*fskin[46]+1.224744871391589*fskin[30]+0.7071067811865475*fskin[22])/vmap_prime_skin[0]; 
  f_cr[15] = (1.58113883008419*fskin[47]+1.224744871391589*fskin[31]+0.7071067811865475*fskin[27])/vmap_prime_skin[0]; 

  f_rl[0] = (1.58113883008419*fedge[32]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  f_rl[1] = (1.58113883008419*fedge[33]-1.224744871391589*fedge[9]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  f_rl[2] = (1.58113883008419*fedge[34]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[2])/vmap_prime_edge[0]; 
  f_rl[3] = (1.58113883008419*fedge[35]-1.224744871391589*fedge[11]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  f_rl[4] = (1.58113883008419*fedge[36]-1.224744871391589*fedge[15]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 
  f_rl[5] = (1.58113883008419*fedge[37]-1.224744871391589*fedge[17]+0.7071067811865475*fedge[6])/vmap_prime_edge[0]; 
  f_rl[6] = (1.58113883008419*fedge[38]-1.224744871391589*fedge[18]+0.7071067811865475*fedge[7])/vmap_prime_edge[0]; 
  f_rl[7] = (1.58113883008419*fedge[39]-1.224744871391589*fedge[19]+0.7071067811865475*fedge[8])/vmap_prime_edge[0]; 
  f_rl[8] = (1.58113883008419*fedge[40]-1.224744871391589*fedge[23]+0.7071067811865475*fedge[12])/vmap_prime_edge[0]; 
  f_rl[9] = (1.58113883008419*fedge[41]-1.224744871391589*fedge[24]+0.7071067811865475*fedge[13])/vmap_prime_edge[0]; 
  f_rl[10] = (1.58113883008419*fedge[42]-1.224744871391589*fedge[25]+0.7071067811865475*fedge[14])/vmap_prime_edge[0]; 
  f_rl[11] = (1.58113883008419*fedge[43]-1.224744871391589*fedge[26]+0.7071067811865475*fedge[16])/vmap_prime_edge[0]; 
  f_rl[12] = (1.58113883008419*fedge[44]-1.224744871391589*fedge[28]+0.7071067811865475*fedge[20])/vmap_prime_edge[0]; 
  f_rl[13] = (1.58113883008419*fedge[45]-1.224744871391589*fedge[29]+0.7071067811865475*fedge[21])/vmap_prime_edge[0]; 
  f_rl[14] = (1.58113883008419*fedge[46]-1.224744871391589*fedge[30]+0.7071067811865475*fedge[22])/vmap_prime_edge[0]; 
  f_rl[15] = (1.58113883008419*fedge[47]-1.224744871391589*fedge[31]+0.7071067811865475*fedge[27])/vmap_prime_edge[0]; 

  fUpR[0] = (0.125*f_cr[15]-0.125*f_rl[15])*sgn_alphaUpR[15]+(0.125*f_cr[14]-0.125*f_rl[14])*sgn_alphaUpR[14]+(0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[13]+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[12]+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[11]+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[10]+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[9]+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[8]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[7]+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[6]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[5]+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[4]+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[3]+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[2]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[1]+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.125*f_cr[14]-0.125*f_rl[14])*sgn_alphaUpR[15]+sgn_alphaUpR[14]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[13]+sgn_alphaUpR[10]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[12]+sgn_alphaUpR[9]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[11]+sgn_alphaUpR[7]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[8]+sgn_alphaUpR[4]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[6]+sgn_alphaUpR[3]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[5]+sgn_alphaUpR[2]*(0.125*f_cr[5]-0.125*f_rl[5])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[1]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[15]+sgn_alphaUpR[13]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[14]+sgn_alphaUpR[10]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[12]+sgn_alphaUpR[8]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[11]+sgn_alphaUpR[6]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[9]+sgn_alphaUpR[4]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[7]+sgn_alphaUpR[3]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[5]+sgn_alphaUpR[1]*(0.125*f_cr[5]-0.125*f_rl[5])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[2]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[15]+sgn_alphaUpR[12]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[14]+sgn_alphaUpR[9]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[13]+sgn_alphaUpR[8]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[11]+sgn_alphaUpR[5]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[10]+sgn_alphaUpR[4]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[7]+sgn_alphaUpR[2]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[6]+sgn_alphaUpR[1]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[3]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[3]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[3]; 
  fUpR[4] = (0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[15]+sgn_alphaUpR[11]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[14]+sgn_alphaUpR[7]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[13]+sgn_alphaUpR[6]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[12]+sgn_alphaUpR[5]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[10]+sgn_alphaUpR[3]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[9]+sgn_alphaUpR[2]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[8]+sgn_alphaUpR[1]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[4]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[4]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[4]; 
  fUpR[5] = (0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[15]+sgn_alphaUpR[10]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[14]+sgn_alphaUpR[13]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[12]+sgn_alphaUpR[4]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[11]+sgn_alphaUpR[3]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[9]+sgn_alphaUpR[8]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[7]+sgn_alphaUpR[6]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[5]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[5]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.125*f_cr[2]-0.125*f_rl[2]); 
  fUpR[6] = (0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[15]+sgn_alphaUpR[9]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[14]+sgn_alphaUpR[12]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[13]+sgn_alphaUpR[4]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[11]+sgn_alphaUpR[2]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[10]+sgn_alphaUpR[8]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[7]+sgn_alphaUpR[5]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[6]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[6]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[6]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.125*f_cr[3]-0.125*f_rl[3]); 
  fUpR[7] = (0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[15]+sgn_alphaUpR[8]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[14]+sgn_alphaUpR[4]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[13]+sgn_alphaUpR[12]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[11]+sgn_alphaUpR[1]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[10]+sgn_alphaUpR[9]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[7]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[7]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[7]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[6]+sgn_alphaUpR[5]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.125*f_cr[3]-0.125*f_rl[3]); 
  fUpR[8] = (0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[15]+sgn_alphaUpR[7]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[14]+sgn_alphaUpR[11]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[13]+sgn_alphaUpR[3]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[12]+sgn_alphaUpR[2]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[10]+sgn_alphaUpR[6]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[9]+sgn_alphaUpR[5]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[8]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[8]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[8]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[9] = (0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[15]+sgn_alphaUpR[6]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[14]+sgn_alphaUpR[3]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[13]+sgn_alphaUpR[11]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[12]+sgn_alphaUpR[1]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[10]+sgn_alphaUpR[7]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[9]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[9]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[9]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[8]+sgn_alphaUpR[5]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[10] = (0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[15]+sgn_alphaUpR[5]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[14]+sgn_alphaUpR[2]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[13]+sgn_alphaUpR[1]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[12]+sgn_alphaUpR[11]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[10]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[10]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[10]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[9]+sgn_alphaUpR[7]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[8]+sgn_alphaUpR[6]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[11] = (0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[15]+sgn_alphaUpR[4]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[14]+sgn_alphaUpR[8]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[13]+sgn_alphaUpR[9]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[12]+sgn_alphaUpR[10]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[11]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[11]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[11]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[7]+sgn_alphaUpR[1]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[6]+sgn_alphaUpR[2]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[3]*(0.125*f_cr[5]-0.125*f_rl[5]); 
  fUpR[12] = (0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[15]+sgn_alphaUpR[3]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[14]+sgn_alphaUpR[6]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[13]+sgn_alphaUpR[7]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[12]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[12]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[12]+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[11]+sgn_alphaUpR[10]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[9]+sgn_alphaUpR[1]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[8]+sgn_alphaUpR[2]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.125*f_cr[5]-0.125*f_rl[5]); 
  fUpR[13] = (0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[15]+sgn_alphaUpR[2]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[14]+sgn_alphaUpR[5]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[13]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[13]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[13]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[12]+sgn_alphaUpR[7]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[11]+sgn_alphaUpR[9]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[10]+sgn_alphaUpR[1]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[8]+sgn_alphaUpR[3]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[6]+sgn_alphaUpR[4]*(0.125*f_cr[6]-0.125*f_rl[6]); 
  fUpR[14] = (0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[15]+sgn_alphaUpR[1]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[14]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[14]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[14]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[13]+sgn_alphaUpR[5]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[12]+sgn_alphaUpR[6]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[11]+sgn_alphaUpR[8]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[10]+sgn_alphaUpR[2]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[9]+sgn_alphaUpR[3]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[7]+sgn_alphaUpR[4]*(0.125*f_cr[7]-0.125*f_rl[7]); 
  fUpR[15] = (0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[15]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[15]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[15]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[14]+sgn_alphaUpR[1]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[13]+sgn_alphaUpR[2]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[12]+sgn_alphaUpR[3]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[11]+sgn_alphaUpR[4]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[10]+sgn_alphaUpR[5]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[9]+sgn_alphaUpR[6]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[8]+sgn_alphaUpR[7]*(0.125*f_cr[8]-0.125*f_rl[8]); 

  } 
  double GhatR[16] = {0.};
  GhatR[0] = 0.25*(alphaR[13]*fUpR[13]+alphaR[11]*fUpR[11]+alphaR[10]*fUpR[10]+alphaR[8]*fUpR[8]+alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.25*(alphaR[10]*fUpR[13]+fUpR[10]*alphaR[13]+alphaR[7]*fUpR[11]+fUpR[7]*alphaR[11]+alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]+alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.25*(alphaR[13]*fUpR[15]+alphaR[10]*fUpR[14]+alphaR[8]*fUpR[12]+alphaR[6]*fUpR[11]+fUpR[6]*alphaR[11]+alphaR[4]*fUpR[9]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.25*(alphaR[8]*fUpR[13]+fUpR[8]*alphaR[13]+alphaR[5]*fUpR[11]+fUpR[5]*alphaR[11]+alphaR[4]*fUpR[10]+fUpR[4]*alphaR[10]+alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.25*(alphaR[11]*fUpR[15]+alphaR[7]*fUpR[14]+alphaR[6]*fUpR[13]+fUpR[6]*alphaR[13]+alphaR[5]*fUpR[12]+alphaR[3]*fUpR[10]+fUpR[3]*alphaR[10]+alphaR[2]*fUpR[9]+alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*(alphaR[10]*fUpR[15]+alphaR[13]*fUpR[14]+alphaR[4]*fUpR[12]+alphaR[3]*fUpR[11]+fUpR[3]*alphaR[11]+alphaR[8]*fUpR[9]+alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.25*(alphaR[4]*fUpR[13]+fUpR[4]*alphaR[13]+alphaR[2]*fUpR[11]+fUpR[2]*alphaR[11]+alphaR[8]*fUpR[10]+fUpR[8]*alphaR[10]+alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]+alphaR[0]*fUpR[6]+fUpR[0]*alphaR[6]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.25*(alphaR[8]*fUpR[15]+alphaR[4]*fUpR[14]+fUpR[12]*alphaR[13]+alphaR[1]*fUpR[11]+fUpR[1]*alphaR[11]+fUpR[9]*alphaR[10]+alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[5]*fUpR[6]+fUpR[5]*alphaR[6]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[7]*fUpR[15]+alphaR[11]*fUpR[14]+alphaR[3]*fUpR[13]+fUpR[3]*alphaR[13]+alphaR[2]*fUpR[12]+alphaR[6]*fUpR[10]+fUpR[6]*alphaR[10]+alphaR[5]*fUpR[9]+alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.25*(alphaR[6]*fUpR[15]+alphaR[3]*fUpR[14]+alphaR[11]*fUpR[13]+fUpR[11]*alphaR[13]+alphaR[1]*fUpR[12]+alphaR[7]*fUpR[10]+fUpR[7]*alphaR[10]+alphaR[0]*fUpR[9]+alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.25*(alphaR[5]*fUpR[15]+alphaR[2]*fUpR[14]+alphaR[1]*fUpR[13]+fUpR[1]*alphaR[13]+alphaR[11]*fUpR[12]+alphaR[0]*fUpR[10]+fUpR[0]*alphaR[10]+alphaR[7]*fUpR[9]+alphaR[6]*fUpR[8]+fUpR[6]*alphaR[8]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.25*(alphaR[4]*fUpR[15]+alphaR[8]*fUpR[14]+fUpR[9]*alphaR[13]+alphaR[10]*fUpR[12]+alphaR[0]*fUpR[11]+fUpR[0]*alphaR[11]+alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6]+fUpR[2]*alphaR[6]+alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]); 
  GhatR[12] = 0.25*(alphaR[3]*fUpR[15]+alphaR[6]*fUpR[14]+alphaR[7]*fUpR[13]+fUpR[7]*alphaR[13]+alphaR[0]*fUpR[12]+alphaR[10]*fUpR[11]+fUpR[10]*alphaR[11]+alphaR[1]*fUpR[9]+alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[13] = 0.25*(alphaR[2]*fUpR[15]+alphaR[5]*fUpR[14]+alphaR[0]*fUpR[13]+fUpR[0]*alphaR[13]+alphaR[7]*fUpR[12]+fUpR[9]*alphaR[11]+alphaR[1]*fUpR[10]+fUpR[1]*alphaR[10]+alphaR[3]*fUpR[8]+fUpR[3]*alphaR[8]+alphaR[4]*fUpR[6]+fUpR[4]*alphaR[6]); 
  GhatR[14] = 0.25*(alphaR[1]*fUpR[15]+alphaR[0]*fUpR[14]+alphaR[5]*fUpR[13]+fUpR[5]*alphaR[13]+alphaR[6]*fUpR[12]+alphaR[8]*fUpR[11]+fUpR[8]*alphaR[11]+alphaR[2]*fUpR[10]+fUpR[2]*alphaR[10]+alphaR[3]*fUpR[9]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]); 
  GhatR[15] = 0.25*(alphaR[0]*fUpR[15]+alphaR[1]*fUpR[14]+alphaR[2]*fUpR[13]+fUpR[2]*alphaR[13]+alphaR[3]*fUpR[12]+alphaR[4]*fUpR[11]+fUpR[4]*alphaR[11]+alphaR[5]*fUpR[10]+fUpR[5]*alphaR[10]+alphaR[6]*fUpR[9]+alphaR[7]*fUpR[8]+fUpR[7]*alphaR[8]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[7] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[8] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdvpar2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdvpar2; 
  out[14] += -0.7071067811865475*GhatR[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[16] += -0.7071067811865475*GhatR[11]*rdvpar2; 
  out[17] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[18] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[19] += -1.224744871391589*GhatR[7]*rdvpar2; 
  out[20] += -0.7071067811865475*GhatR[12]*rdvpar2; 
  out[21] += -0.7071067811865475*GhatR[13]*rdvpar2; 
  out[22] += -0.7071067811865475*GhatR[14]*rdvpar2; 
  out[23] += -1.224744871391589*GhatR[8]*rdvpar2; 
  out[24] += -1.224744871391589*GhatR[9]*rdvpar2; 
  out[25] += -1.224744871391589*GhatR[10]*rdvpar2; 
  out[26] += -1.224744871391589*GhatR[11]*rdvpar2; 
  out[27] += -0.7071067811865475*GhatR[15]*rdvpar2; 
  out[28] += -1.224744871391589*GhatR[12]*rdvpar2; 
  out[29] += -1.224744871391589*GhatR[13]*rdvpar2; 
  out[30] += -1.224744871391589*GhatR[14]*rdvpar2; 
  out[31] += -1.224744871391589*GhatR[15]*rdvpar2; 
  out[32] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[33] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[34] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[35] += -1.58113883008419*GhatR[3]*rdvpar2; 
  out[36] += -1.58113883008419*GhatR[4]*rdvpar2; 
  out[37] += -1.58113883008419*GhatR[5]*rdvpar2; 
  out[38] += -1.58113883008419*GhatR[6]*rdvpar2; 
  out[39] += -1.58113883008419*GhatR[7]*rdvpar2; 
  out[40] += -1.58113883008419*GhatR[8]*rdvpar2; 
  out[41] += -1.58113883008419*GhatR[9]*rdvpar2; 
  out[42] += -1.58113883008419*GhatR[10]*rdvpar2; 
  out[43] += -1.58113883008419*GhatR[11]*rdvpar2; 
  out[44] += -1.58113883008419*GhatR[12]*rdvpar2; 
  out[45] += -1.58113883008419*GhatR[13]*rdvpar2; 
  out[46] += -1.58113883008419*GhatR[14]*rdvpar2; 
  out[47] += -1.58113883008419*GhatR[15]*rdvpar2; 

  } else { 

  double fUpL[16] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = (1.58113883008419*fedge[32]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  fUpL[1] = (1.58113883008419*fedge[33]+1.224744871391589*fedge[9]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  fUpL[2] = (1.58113883008419*fedge[34]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[2])/vmap_prime_edge[0]; 
  fUpL[3] = (1.58113883008419*fedge[35]+1.224744871391589*fedge[11]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  fUpL[4] = (1.58113883008419*fedge[36]+1.224744871391589*fedge[15]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 
  fUpL[5] = (1.58113883008419*fedge[37]+1.224744871391589*fedge[17]+0.7071067811865475*fedge[6])/vmap_prime_edge[0]; 
  fUpL[6] = (1.58113883008419*fedge[38]+1.224744871391589*fedge[18]+0.7071067811865475*fedge[7])/vmap_prime_edge[0]; 
  fUpL[7] = (1.58113883008419*fedge[39]+1.224744871391589*fedge[19]+0.7071067811865475*fedge[8])/vmap_prime_edge[0]; 
  fUpL[8] = (1.58113883008419*fedge[40]+1.224744871391589*fedge[23]+0.7071067811865475*fedge[12])/vmap_prime_edge[0]; 
  fUpL[9] = (1.58113883008419*fedge[41]+1.224744871391589*fedge[24]+0.7071067811865475*fedge[13])/vmap_prime_edge[0]; 
  fUpL[10] = (1.58113883008419*fedge[42]+1.224744871391589*fedge[25]+0.7071067811865475*fedge[14])/vmap_prime_edge[0]; 
  fUpL[11] = (1.58113883008419*fedge[43]+1.224744871391589*fedge[26]+0.7071067811865475*fedge[16])/vmap_prime_edge[0]; 
  fUpL[12] = (1.58113883008419*fedge[44]+1.224744871391589*fedge[28]+0.7071067811865475*fedge[20])/vmap_prime_edge[0]; 
  fUpL[13] = (1.58113883008419*fedge[45]+1.224744871391589*fedge[29]+0.7071067811865475*fedge[21])/vmap_prime_edge[0]; 
  fUpL[14] = (1.58113883008419*fedge[46]+1.224744871391589*fedge[30]+0.7071067811865475*fedge[22])/vmap_prime_edge[0]; 
  fUpL[15] = (1.58113883008419*fedge[47]+1.224744871391589*fedge[31]+0.7071067811865475*fedge[27])/vmap_prime_edge[0]; 
    } else { 
  fUpL[0] = (1.58113883008419*fskin[32]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  fUpL[1] = (1.58113883008419*fskin[33]-1.224744871391589*fskin[9]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  fUpL[2] = (1.58113883008419*fskin[34]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[2])/vmap_prime_skin[0]; 
  fUpL[3] = (1.58113883008419*fskin[35]-1.224744871391589*fskin[11]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  fUpL[4] = (1.58113883008419*fskin[36]-1.224744871391589*fskin[15]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 
  fUpL[5] = (1.58113883008419*fskin[37]-1.224744871391589*fskin[17]+0.7071067811865475*fskin[6])/vmap_prime_skin[0]; 
  fUpL[6] = (1.58113883008419*fskin[38]-1.224744871391589*fskin[18]+0.7071067811865475*fskin[7])/vmap_prime_skin[0]; 
  fUpL[7] = (1.58113883008419*fskin[39]-1.224744871391589*fskin[19]+0.7071067811865475*fskin[8])/vmap_prime_skin[0]; 
  fUpL[8] = (1.58113883008419*fskin[40]-1.224744871391589*fskin[23]+0.7071067811865475*fskin[12])/vmap_prime_skin[0]; 
  fUpL[9] = (1.58113883008419*fskin[41]-1.224744871391589*fskin[24]+0.7071067811865475*fskin[13])/vmap_prime_skin[0]; 
  fUpL[10] = (1.58113883008419*fskin[42]-1.224744871391589*fskin[25]+0.7071067811865475*fskin[14])/vmap_prime_skin[0]; 
  fUpL[11] = (1.58113883008419*fskin[43]-1.224744871391589*fskin[26]+0.7071067811865475*fskin[16])/vmap_prime_skin[0]; 
  fUpL[12] = (1.58113883008419*fskin[44]-1.224744871391589*fskin[28]+0.7071067811865475*fskin[20])/vmap_prime_skin[0]; 
  fUpL[13] = (1.58113883008419*fskin[45]-1.224744871391589*fskin[29]+0.7071067811865475*fskin[21])/vmap_prime_skin[0]; 
  fUpL[14] = (1.58113883008419*fskin[46]-1.224744871391589*fskin[30]+0.7071067811865475*fskin[22])/vmap_prime_skin[0]; 
  fUpL[15] = (1.58113883008419*fskin[47]-1.224744871391589*fskin[31]+0.7071067811865475*fskin[27])/vmap_prime_skin[0]; 
    } 
  } else { 
  double f_lr[16] = {0.};
  double f_cl[16] = {0.};
  double sgn_alphaUpL[16] = {0.};
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = (1.58113883008419*fedge[32]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[0])/vmap_prime_edge[0]; 
  f_lr[1] = (1.58113883008419*fedge[33]+1.224744871391589*fedge[9]+0.7071067811865475*fedge[1])/vmap_prime_edge[0]; 
  f_lr[2] = (1.58113883008419*fedge[34]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[2])/vmap_prime_edge[0]; 
  f_lr[3] = (1.58113883008419*fedge[35]+1.224744871391589*fedge[11]+0.7071067811865475*fedge[3])/vmap_prime_edge[0]; 
  f_lr[4] = (1.58113883008419*fedge[36]+1.224744871391589*fedge[15]+0.7071067811865475*fedge[5])/vmap_prime_edge[0]; 
  f_lr[5] = (1.58113883008419*fedge[37]+1.224744871391589*fedge[17]+0.7071067811865475*fedge[6])/vmap_prime_edge[0]; 
  f_lr[6] = (1.58113883008419*fedge[38]+1.224744871391589*fedge[18]+0.7071067811865475*fedge[7])/vmap_prime_edge[0]; 
  f_lr[7] = (1.58113883008419*fedge[39]+1.224744871391589*fedge[19]+0.7071067811865475*fedge[8])/vmap_prime_edge[0]; 
  f_lr[8] = (1.58113883008419*fedge[40]+1.224744871391589*fedge[23]+0.7071067811865475*fedge[12])/vmap_prime_edge[0]; 
  f_lr[9] = (1.58113883008419*fedge[41]+1.224744871391589*fedge[24]+0.7071067811865475*fedge[13])/vmap_prime_edge[0]; 
  f_lr[10] = (1.58113883008419*fedge[42]+1.224744871391589*fedge[25]+0.7071067811865475*fedge[14])/vmap_prime_edge[0]; 
  f_lr[11] = (1.58113883008419*fedge[43]+1.224744871391589*fedge[26]+0.7071067811865475*fedge[16])/vmap_prime_edge[0]; 
  f_lr[12] = (1.58113883008419*fedge[44]+1.224744871391589*fedge[28]+0.7071067811865475*fedge[20])/vmap_prime_edge[0]; 
  f_lr[13] = (1.58113883008419*fedge[45]+1.224744871391589*fedge[29]+0.7071067811865475*fedge[21])/vmap_prime_edge[0]; 
  f_lr[14] = (1.58113883008419*fedge[46]+1.224744871391589*fedge[30]+0.7071067811865475*fedge[22])/vmap_prime_edge[0]; 
  f_lr[15] = (1.58113883008419*fedge[47]+1.224744871391589*fedge[31]+0.7071067811865475*fedge[27])/vmap_prime_edge[0]; 

  f_cl[0] = (1.58113883008419*fskin[32]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[0])/vmap_prime_skin[0]; 
  f_cl[1] = (1.58113883008419*fskin[33]-1.224744871391589*fskin[9]+0.7071067811865475*fskin[1])/vmap_prime_skin[0]; 
  f_cl[2] = (1.58113883008419*fskin[34]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[2])/vmap_prime_skin[0]; 
  f_cl[3] = (1.58113883008419*fskin[35]-1.224744871391589*fskin[11]+0.7071067811865475*fskin[3])/vmap_prime_skin[0]; 
  f_cl[4] = (1.58113883008419*fskin[36]-1.224744871391589*fskin[15]+0.7071067811865475*fskin[5])/vmap_prime_skin[0]; 
  f_cl[5] = (1.58113883008419*fskin[37]-1.224744871391589*fskin[17]+0.7071067811865475*fskin[6])/vmap_prime_skin[0]; 
  f_cl[6] = (1.58113883008419*fskin[38]-1.224744871391589*fskin[18]+0.7071067811865475*fskin[7])/vmap_prime_skin[0]; 
  f_cl[7] = (1.58113883008419*fskin[39]-1.224744871391589*fskin[19]+0.7071067811865475*fskin[8])/vmap_prime_skin[0]; 
  f_cl[8] = (1.58113883008419*fskin[40]-1.224744871391589*fskin[23]+0.7071067811865475*fskin[12])/vmap_prime_skin[0]; 
  f_cl[9] = (1.58113883008419*fskin[41]-1.224744871391589*fskin[24]+0.7071067811865475*fskin[13])/vmap_prime_skin[0]; 
  f_cl[10] = (1.58113883008419*fskin[42]-1.224744871391589*fskin[25]+0.7071067811865475*fskin[14])/vmap_prime_skin[0]; 
  f_cl[11] = (1.58113883008419*fskin[43]-1.224744871391589*fskin[26]+0.7071067811865475*fskin[16])/vmap_prime_skin[0]; 
  f_cl[12] = (1.58113883008419*fskin[44]-1.224744871391589*fskin[28]+0.7071067811865475*fskin[20])/vmap_prime_skin[0]; 
  f_cl[13] = (1.58113883008419*fskin[45]-1.224744871391589*fskin[29]+0.7071067811865475*fskin[21])/vmap_prime_skin[0]; 
  f_cl[14] = (1.58113883008419*fskin[46]-1.224744871391589*fskin[30]+0.7071067811865475*fskin[22])/vmap_prime_skin[0]; 
  f_cl[15] = (1.58113883008419*fskin[47]-1.224744871391589*fskin[31]+0.7071067811865475*fskin[27])/vmap_prime_skin[0]; 

  fUpL[0] = (0.125*f_lr[15]-0.125*f_cl[15])*sgn_alphaUpL[15]+(0.125*f_lr[14]-0.125*f_cl[14])*sgn_alphaUpL[14]+(0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[13]+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[12]+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[11]+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[10]+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[9]+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[8]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[7]+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[6]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[5]+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[4]+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[3]+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[2]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[1]+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.125*f_lr[14]-0.125*f_cl[14])*sgn_alphaUpL[15]+sgn_alphaUpL[14]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[13]+sgn_alphaUpL[10]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[12]+sgn_alphaUpL[9]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[11]+sgn_alphaUpL[7]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[8]+sgn_alphaUpL[4]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[6]+sgn_alphaUpL[3]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[2]*(0.125*f_lr[5]-0.125*f_cl[5])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[1]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[15]+sgn_alphaUpL[13]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[14]+sgn_alphaUpL[10]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[12]+sgn_alphaUpL[8]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[11]+sgn_alphaUpL[6]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[9]+sgn_alphaUpL[4]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[7]+sgn_alphaUpL[3]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[5]+sgn_alphaUpL[1]*(0.125*f_lr[5]-0.125*f_cl[5])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[2]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[15]+sgn_alphaUpL[12]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[14]+sgn_alphaUpL[9]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[13]+sgn_alphaUpL[8]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[11]+sgn_alphaUpL[5]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[10]+sgn_alphaUpL[4]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[7]+sgn_alphaUpL[2]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[6]+sgn_alphaUpL[1]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[3]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[3]; 
  fUpL[4] = (0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[15]+sgn_alphaUpL[11]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[14]+sgn_alphaUpL[7]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[13]+sgn_alphaUpL[6]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[12]+sgn_alphaUpL[5]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[10]+sgn_alphaUpL[3]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[9]+sgn_alphaUpL[2]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[8]+sgn_alphaUpL[1]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[4]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[4]; 
  fUpL[5] = (0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[15]+sgn_alphaUpL[10]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[14]+sgn_alphaUpL[13]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[12]+sgn_alphaUpL[4]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[11]+sgn_alphaUpL[3]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[9]+sgn_alphaUpL[8]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[7]+sgn_alphaUpL[6]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[5]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[5]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.125*f_lr[2]-0.125*f_cl[2]); 
  fUpL[6] = (0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[15]+sgn_alphaUpL[9]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[14]+sgn_alphaUpL[12]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[13]+sgn_alphaUpL[4]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[11]+sgn_alphaUpL[2]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[10]+sgn_alphaUpL[8]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[7]+sgn_alphaUpL[5]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[6]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[6]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[6]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.125*f_lr[3]-0.125*f_cl[3]); 
  fUpL[7] = (0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[15]+sgn_alphaUpL[8]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[14]+sgn_alphaUpL[4]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[13]+sgn_alphaUpL[12]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[11]+sgn_alphaUpL[1]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[10]+sgn_alphaUpL[9]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[7]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[7]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[7]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[6]+sgn_alphaUpL[5]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.125*f_lr[3]-0.125*f_cl[3]); 
  fUpL[8] = (0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[15]+sgn_alphaUpL[7]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[14]+sgn_alphaUpL[11]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[13]+sgn_alphaUpL[3]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[12]+sgn_alphaUpL[2]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[10]+sgn_alphaUpL[6]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[9]+sgn_alphaUpL[5]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[8]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[8]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[8]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[9] = (0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[15]+sgn_alphaUpL[6]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[14]+sgn_alphaUpL[3]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[13]+sgn_alphaUpL[11]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[12]+sgn_alphaUpL[1]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[10]+sgn_alphaUpL[7]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[9]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[9]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[9]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[8]+sgn_alphaUpL[5]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[10] = (0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[15]+sgn_alphaUpL[5]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[14]+sgn_alphaUpL[2]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[13]+sgn_alphaUpL[1]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[12]+sgn_alphaUpL[11]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[10]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[10]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[10]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[9]+sgn_alphaUpL[7]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[8]+sgn_alphaUpL[6]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[11] = (0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[15]+sgn_alphaUpL[4]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[14]+sgn_alphaUpL[8]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[13]+sgn_alphaUpL[9]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[12]+sgn_alphaUpL[10]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[11]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[11]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[11]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[7]+sgn_alphaUpL[1]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[6]+sgn_alphaUpL[2]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[3]*(0.125*f_lr[5]-0.125*f_cl[5]); 
  fUpL[12] = (0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[15]+sgn_alphaUpL[3]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[14]+sgn_alphaUpL[6]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[13]+sgn_alphaUpL[7]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[12]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[12]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[12]+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[11]+sgn_alphaUpL[10]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[9]+sgn_alphaUpL[1]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[8]+sgn_alphaUpL[2]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.125*f_lr[5]-0.125*f_cl[5]); 
  fUpL[13] = (0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[15]+sgn_alphaUpL[2]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[14]+sgn_alphaUpL[5]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[13]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[13]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[13]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[12]+sgn_alphaUpL[7]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[11]+sgn_alphaUpL[9]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[10]+sgn_alphaUpL[1]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[8]+sgn_alphaUpL[3]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[6]+sgn_alphaUpL[4]*(0.125*f_lr[6]-0.125*f_cl[6]); 
  fUpL[14] = (0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[15]+sgn_alphaUpL[1]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[14]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[14]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[14]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[13]+sgn_alphaUpL[5]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[12]+sgn_alphaUpL[6]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[11]+sgn_alphaUpL[8]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[10]+sgn_alphaUpL[2]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[9]+sgn_alphaUpL[3]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[7]+sgn_alphaUpL[4]*(0.125*f_lr[7]-0.125*f_cl[7]); 
  fUpL[15] = (0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[15]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[15]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[15]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[14]+sgn_alphaUpL[1]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[13]+sgn_alphaUpL[2]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[12]+sgn_alphaUpL[3]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[11]+sgn_alphaUpL[4]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[10]+sgn_alphaUpL[5]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[9]+sgn_alphaUpL[6]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[8]+sgn_alphaUpL[7]*(0.125*f_lr[8]-0.125*f_cl[8]); 

  } 
  double GhatL[16] = {0.};
  GhatL[0] = 0.25*(alphaL[13]*fUpL[13]+alphaL[11]*fUpL[11]+alphaL[10]*fUpL[10]+alphaL[8]*fUpL[8]+alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.25*(alphaL[10]*fUpL[13]+fUpL[10]*alphaL[13]+alphaL[7]*fUpL[11]+fUpL[7]*alphaL[11]+alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]+alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.25*(alphaL[13]*fUpL[15]+alphaL[10]*fUpL[14]+alphaL[8]*fUpL[12]+alphaL[6]*fUpL[11]+fUpL[6]*alphaL[11]+alphaL[4]*fUpL[9]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.25*(alphaL[8]*fUpL[13]+fUpL[8]*alphaL[13]+alphaL[5]*fUpL[11]+fUpL[5]*alphaL[11]+alphaL[4]*fUpL[10]+fUpL[4]*alphaL[10]+alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.25*(alphaL[11]*fUpL[15]+alphaL[7]*fUpL[14]+alphaL[6]*fUpL[13]+fUpL[6]*alphaL[13]+alphaL[5]*fUpL[12]+alphaL[3]*fUpL[10]+fUpL[3]*alphaL[10]+alphaL[2]*fUpL[9]+alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*(alphaL[10]*fUpL[15]+alphaL[13]*fUpL[14]+alphaL[4]*fUpL[12]+alphaL[3]*fUpL[11]+fUpL[3]*alphaL[11]+alphaL[8]*fUpL[9]+alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.25*(alphaL[4]*fUpL[13]+fUpL[4]*alphaL[13]+alphaL[2]*fUpL[11]+fUpL[2]*alphaL[11]+alphaL[8]*fUpL[10]+fUpL[8]*alphaL[10]+alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]+alphaL[0]*fUpL[6]+fUpL[0]*alphaL[6]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.25*(alphaL[8]*fUpL[15]+alphaL[4]*fUpL[14]+fUpL[12]*alphaL[13]+alphaL[1]*fUpL[11]+fUpL[1]*alphaL[11]+fUpL[9]*alphaL[10]+alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[5]*fUpL[6]+fUpL[5]*alphaL[6]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[7]*fUpL[15]+alphaL[11]*fUpL[14]+alphaL[3]*fUpL[13]+fUpL[3]*alphaL[13]+alphaL[2]*fUpL[12]+alphaL[6]*fUpL[10]+fUpL[6]*alphaL[10]+alphaL[5]*fUpL[9]+alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.25*(alphaL[6]*fUpL[15]+alphaL[3]*fUpL[14]+alphaL[11]*fUpL[13]+fUpL[11]*alphaL[13]+alphaL[1]*fUpL[12]+alphaL[7]*fUpL[10]+fUpL[7]*alphaL[10]+alphaL[0]*fUpL[9]+alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.25*(alphaL[5]*fUpL[15]+alphaL[2]*fUpL[14]+alphaL[1]*fUpL[13]+fUpL[1]*alphaL[13]+alphaL[11]*fUpL[12]+alphaL[0]*fUpL[10]+fUpL[0]*alphaL[10]+alphaL[7]*fUpL[9]+alphaL[6]*fUpL[8]+fUpL[6]*alphaL[8]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.25*(alphaL[4]*fUpL[15]+alphaL[8]*fUpL[14]+fUpL[9]*alphaL[13]+alphaL[10]*fUpL[12]+alphaL[0]*fUpL[11]+fUpL[0]*alphaL[11]+alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6]+fUpL[2]*alphaL[6]+alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]); 
  GhatL[12] = 0.25*(alphaL[3]*fUpL[15]+alphaL[6]*fUpL[14]+alphaL[7]*fUpL[13]+fUpL[7]*alphaL[13]+alphaL[0]*fUpL[12]+alphaL[10]*fUpL[11]+fUpL[10]*alphaL[11]+alphaL[1]*fUpL[9]+alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[13] = 0.25*(alphaL[2]*fUpL[15]+alphaL[5]*fUpL[14]+alphaL[0]*fUpL[13]+fUpL[0]*alphaL[13]+alphaL[7]*fUpL[12]+fUpL[9]*alphaL[11]+alphaL[1]*fUpL[10]+fUpL[1]*alphaL[10]+alphaL[3]*fUpL[8]+fUpL[3]*alphaL[8]+alphaL[4]*fUpL[6]+fUpL[4]*alphaL[6]); 
  GhatL[14] = 0.25*(alphaL[1]*fUpL[15]+alphaL[0]*fUpL[14]+alphaL[5]*fUpL[13]+fUpL[5]*alphaL[13]+alphaL[6]*fUpL[12]+alphaL[8]*fUpL[11]+fUpL[8]*alphaL[11]+alphaL[2]*fUpL[10]+fUpL[2]*alphaL[10]+alphaL[3]*fUpL[9]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]); 
  GhatL[15] = 0.25*(alphaL[0]*fUpL[15]+alphaL[1]*fUpL[14]+alphaL[2]*fUpL[13]+fUpL[2]*alphaL[13]+alphaL[3]*fUpL[12]+alphaL[4]*fUpL[11]+fUpL[4]*alphaL[11]+alphaL[5]*fUpL[10]+fUpL[5]*alphaL[10]+alphaL[6]*fUpL[9]+alphaL[7]*fUpL[8]+fUpL[7]*alphaL[8]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[7] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[8] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdvpar2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdvpar2; 
  out[14] += 0.7071067811865475*GhatL[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[16] += 0.7071067811865475*GhatL[11]*rdvpar2; 
  out[17] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[18] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[19] += -1.224744871391589*GhatL[7]*rdvpar2; 
  out[20] += 0.7071067811865475*GhatL[12]*rdvpar2; 
  out[21] += 0.7071067811865475*GhatL[13]*rdvpar2; 
  out[22] += 0.7071067811865475*GhatL[14]*rdvpar2; 
  out[23] += -1.224744871391589*GhatL[8]*rdvpar2; 
  out[24] += -1.224744871391589*GhatL[9]*rdvpar2; 
  out[25] += -1.224744871391589*GhatL[10]*rdvpar2; 
  out[26] += -1.224744871391589*GhatL[11]*rdvpar2; 
  out[27] += 0.7071067811865475*GhatL[15]*rdvpar2; 
  out[28] += -1.224744871391589*GhatL[12]*rdvpar2; 
  out[29] += -1.224744871391589*GhatL[13]*rdvpar2; 
  out[30] += -1.224744871391589*GhatL[14]*rdvpar2; 
  out[31] += -1.224744871391589*GhatL[15]*rdvpar2; 
  out[32] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[33] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[34] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[35] += 1.58113883008419*GhatL[3]*rdvpar2; 
  out[36] += 1.58113883008419*GhatL[4]*rdvpar2; 
  out[37] += 1.58113883008419*GhatL[5]*rdvpar2; 
  out[38] += 1.58113883008419*GhatL[6]*rdvpar2; 
  out[39] += 1.58113883008419*GhatL[7]*rdvpar2; 
  out[40] += 1.58113883008419*GhatL[8]*rdvpar2; 
  out[41] += 1.58113883008419*GhatL[9]*rdvpar2; 
  out[42] += 1.58113883008419*GhatL[10]*rdvpar2; 
  out[43] += 1.58113883008419*GhatL[11]*rdvpar2; 
  out[44] += 1.58113883008419*GhatL[12]*rdvpar2; 
  out[45] += 1.58113883008419*GhatL[13]*rdvpar2; 
  out[46] += 1.58113883008419*GhatL[14]*rdvpar2; 
  out[47] += 1.58113883008419*GhatL[15]*rdvpar2; 

  } 

  double vmap_prime_min = fmin(fabs(vmap_prime_edge[0]),fabs(vmap_prime_skin[0]));
  double cflFreq = fmax(fabs(alphaL[0]/vmap_prime_min), fabs(alphaR[0]/vmap_prime_min)); 
  return 0.625*rdvpar2*cflFreq; 

} 
