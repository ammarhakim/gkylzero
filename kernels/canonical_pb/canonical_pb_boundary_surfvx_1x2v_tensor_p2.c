#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x2v_tensor_p2(const double *w, const double *dxv, const double *hamil, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // alpha_surf_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // sgn_alpha_surf_edge: sign(alpha_surf_edge) at quadrature points.
  // sgn_alpha_surf_skin: sign(alpha_surf_skin) at quadrature points.
  // const_sgn_alpha_edge: Boolean array true if sign(alpha_surf_edge) is only one sign, either +1 or -1.
  // const_sgn_alpha_skin: Boolean array true if sign(alpha_surf_skin) is only one sign, either +1 or -1.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double rdvx2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_skin[9];
  const double *alphaR = &alpha_surf_edge[9];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_skin[9];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_edge[9];
  const int *const_sgn_alphaL = &const_sgn_alpha_skin[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_edge[1];

  if (edge == -1) { 

  double fUpR[9] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.58113883008419*fskin[8]+1.224744871391589*fskin[2]+0.7071067811865475*fskin[0]; 
  fUpR[1] = 1.58113883008419*fskin[12]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[1]; 
  fUpR[2] = 1.58113883008419*fskin[14]+1.224744871391589*fskin[6]+0.7071067811865475*fskin[3]; 
  fUpR[3] = 1.58113883008419*fskin[18]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[5]; 
  fUpR[4] = 1.58113883008419*fskin[20]+1.224744871391589*fskin[11]+0.7071067811865475*fskin[7]; 
  fUpR[5] = 1.58113883008419*fskin[22]+1.224744871391589*fskin[16]+0.7071067811865475*fskin[9]; 
  fUpR[6] = 1.58113883008419*fskin[23]+1.224744871391589*fskin[17]+0.7071067811865475*fskin[13]; 
  fUpR[7] = 1.58113883008419*fskin[25]+1.224744871391589*fskin[19]+0.7071067811865475*fskin[15]; 
  fUpR[8] = 1.58113883008419*fskin[26]+1.224744871391589*fskin[24]+0.7071067811865475*fskin[21]; 
    } else { 
  fUpR[0] = 1.58113883008419*fedge[8]-1.224744871391589*fedge[2]+0.7071067811865475*fedge[0]; 
  fUpR[1] = 1.58113883008419*fedge[12]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[1]; 
  fUpR[2] = 1.58113883008419*fedge[14]-1.224744871391589*fedge[6]+0.7071067811865475*fedge[3]; 
  fUpR[3] = 1.58113883008419*fedge[18]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[5]; 
  fUpR[4] = 1.58113883008419*fedge[20]-1.224744871391589*fedge[11]+0.7071067811865475*fedge[7]; 
  fUpR[5] = 1.58113883008419*fedge[22]-1.224744871391589*fedge[16]+0.7071067811865475*fedge[9]; 
  fUpR[6] = 1.58113883008419*fedge[23]-1.224744871391589*fedge[17]+0.7071067811865475*fedge[13]; 
  fUpR[7] = 1.58113883008419*fedge[25]-1.224744871391589*fedge[19]+0.7071067811865475*fedge[15]; 
  fUpR[8] = 1.58113883008419*fedge[26]-1.224744871391589*fedge[24]+0.7071067811865475*fedge[21]; 
    } 
  } else { 
  double f_cr[9] = {0.};
  double f_rl[9] = {0.};
  double sgn_alphaUpR[9] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.58113883008419*fskin[8]+1.224744871391589*fskin[2]+0.7071067811865475*fskin[0]; 
  f_cr[1] = 1.58113883008419*fskin[12]+1.224744871391589*fskin[4]+0.7071067811865475*fskin[1]; 
  f_cr[2] = 1.58113883008419*fskin[14]+1.224744871391589*fskin[6]+0.7071067811865475*fskin[3]; 
  f_cr[3] = 1.58113883008419*fskin[18]+1.224744871391589*fskin[10]+0.7071067811865475*fskin[5]; 
  f_cr[4] = 1.58113883008419*fskin[20]+1.224744871391589*fskin[11]+0.7071067811865475*fskin[7]; 
  f_cr[5] = 1.58113883008419*fskin[22]+1.224744871391589*fskin[16]+0.7071067811865475*fskin[9]; 
  f_cr[6] = 1.58113883008419*fskin[23]+1.224744871391589*fskin[17]+0.7071067811865475*fskin[13]; 
  f_cr[7] = 1.58113883008419*fskin[25]+1.224744871391589*fskin[19]+0.7071067811865475*fskin[15]; 
  f_cr[8] = 1.58113883008419*fskin[26]+1.224744871391589*fskin[24]+0.7071067811865475*fskin[21]; 

  f_rl[0] = 1.58113883008419*fedge[8]-1.224744871391589*fedge[2]+0.7071067811865475*fedge[0]; 
  f_rl[1] = 1.58113883008419*fedge[12]-1.224744871391589*fedge[4]+0.7071067811865475*fedge[1]; 
  f_rl[2] = 1.58113883008419*fedge[14]-1.224744871391589*fedge[6]+0.7071067811865475*fedge[3]; 
  f_rl[3] = 1.58113883008419*fedge[18]-1.224744871391589*fedge[10]+0.7071067811865475*fedge[5]; 
  f_rl[4] = 1.58113883008419*fedge[20]-1.224744871391589*fedge[11]+0.7071067811865475*fedge[7]; 
  f_rl[5] = 1.58113883008419*fedge[22]-1.224744871391589*fedge[16]+0.7071067811865475*fedge[9]; 
  f_rl[6] = 1.58113883008419*fedge[23]-1.224744871391589*fedge[17]+0.7071067811865475*fedge[13]; 
  f_rl[7] = 1.58113883008419*fedge[25]-1.224744871391589*fedge[19]+0.7071067811865475*fedge[15]; 
  f_rl[8] = 1.58113883008419*fedge[26]-1.224744871391589*fedge[24]+0.7071067811865475*fedge[21]; 

  fUpR[0] = (0.25*f_cr[8]-0.25*f_rl[8])*sgn_alphaUpR[8]+(0.25*f_cr[7]-0.25*f_rl[7])*sgn_alphaUpR[7]+(0.25*f_cr[6]-0.25*f_rl[6])*sgn_alphaUpR[6]+(0.25*f_cr[5]-0.25*f_rl[5])*sgn_alphaUpR[5]+(0.25*f_cr[4]-0.25*f_rl[4])*sgn_alphaUpR[4]+(0.25*f_cr[3]-0.25*f_rl[3])*sgn_alphaUpR[3]+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[2]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[1]+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])*sgn_alphaUpR[8]+sgn_alphaUpR[7]*((-0.223606797749979*f_rl[8])+0.223606797749979*f_cr[8]-0.2500000000000001*f_rl[5]+0.2500000000000001*f_cr[5])+sgn_alphaUpR[5]*(0.2500000000000001*f_cr[7]-0.2500000000000001*f_rl[7])+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[6]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[4]-0.223606797749979*f_rl[4])+(0.25*f_cr[2]-0.25*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[1]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])*sgn_alphaUpR[8]+sgn_alphaUpR[6]*(0.223606797749979*f_cr[8]-0.223606797749979*f_rl[8])+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[7]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])+(0.2500000000000001*f_cr[4]-0.2500000000000001*f_rl[4])*sgn_alphaUpR[6]+sgn_alphaUpR[4]*(0.2500000000000001*f_cr[6]-0.2500000000000001*f_rl[6])+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[5]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[5]-0.223606797749979*f_rl[5])+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.25*f_cr[3]-0.25*f_rl[3])+(0.25*f_cr[0]-0.25*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[2]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.2*f_cr[3]-0.2*f_rl[3])*sgn_alphaUpR[8]+sgn_alphaUpR[3]*(0.2*f_cr[8]-0.2*f_rl[8])+((-0.2*f_rl[6])+0.2*f_cr[6]-0.223606797749979*f_rl[2]+0.223606797749979*f_cr[2])*sgn_alphaUpR[7]+((-0.2*sgn_alphaUpR[6])-0.223606797749979*sgn_alphaUpR[2])*f_rl[7]+(0.2*sgn_alphaUpR[6]+0.223606797749979*sgn_alphaUpR[2])*f_cr[7]+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[6]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[3]*(0.223606797749979*f_cr[5]-0.223606797749979*f_rl[5])+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*((-0.223606797749979*f_rl[4])+0.223606797749979*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])+(0.5-0.25*sgn_alphaUpR[0])*f_rl[3]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[3]+(0.25*f_cr[1]-0.25*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.25*f_cr[2]-0.25*f_rl[2]); 
  fUpR[4] = ((-0.159719141249985*f_rl[8])+0.159719141249985*f_cr[8]-0.25*f_rl[5]+0.25*f_cr[5])*sgn_alphaUpR[8]+sgn_alphaUpR[5]*(0.25*f_cr[8]-0.25*f_rl[8])+(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])*sgn_alphaUpR[7]+((-0.159719141249985*f_rl[6])+0.159719141249985*f_cr[6]-0.2500000000000001*f_rl[2]+0.2500000000000001*f_cr[2])*sgn_alphaUpR[6]+sgn_alphaUpR[2]*(0.2500000000000001*f_cr[6]-0.2500000000000001*f_rl[6])+((-0.159719141249985*f_rl[4])+0.159719141249985*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[4]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[4]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[4]+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[3]+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[1]; 
  fUpR[5] = ((-0.159719141249985*f_rl[8])+0.159719141249985*f_cr[8]-0.25*f_rl[4]+0.25*f_cr[4])*sgn_alphaUpR[8]+sgn_alphaUpR[4]*(0.25*f_cr[8]-0.25*f_rl[8])+((-0.159719141249985*f_rl[7])+0.159719141249985*f_cr[7]-0.2500000000000001*f_rl[1]+0.2500000000000001*f_cr[1])*sgn_alphaUpR[7]+sgn_alphaUpR[1]*(0.2500000000000001*f_cr[7]-0.2500000000000001*f_rl[7])+(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])*sgn_alphaUpR[6]+((-0.159719141249985*f_rl[5])+0.159719141249985*f_cr[5]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[5]+(0.5-0.25*sgn_alphaUpR[0])*f_rl[5]+(0.25*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3])*sgn_alphaUpR[3]+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[2]; 
  fUpR[6] = ((-0.1428571428571428*f_rl[6])+0.1428571428571428*f_cr[6]-0.223606797749979*f_rl[2]+0.223606797749979*f_cr[2])*sgn_alphaUpR[8]+((-0.1428571428571428*sgn_alphaUpR[6])-0.223606797749979*sgn_alphaUpR[2])*f_rl[8]+(0.1428571428571428*sgn_alphaUpR[6]+0.223606797749979*sgn_alphaUpR[2])*f_cr[8]+(0.2*f_cr[3]-0.2*f_rl[3])*sgn_alphaUpR[7]+sgn_alphaUpR[3]*(0.2*f_cr[7]-0.2*f_rl[7])+((-0.223606797749979*f_rl[5])+0.223606797749979*f_cr[5]-0.159719141249985*f_rl[4]+0.159719141249985*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[6]+((-0.223606797749979*sgn_alphaUpR[5])-0.159719141249985*sgn_alphaUpR[4]-0.25*sgn_alphaUpR[0]+0.5)*f_rl[6]+(0.223606797749979*sgn_alphaUpR[5]+0.159719141249985*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[6]+(0.2500000000000001*f_cr[2]-0.2500000000000001*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.2500000000000001*f_cr[4]-0.2500000000000001*f_rl[4])+(0.223606797749979*f_cr[1]-0.223606797749979*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3]); 
  fUpR[7] = ((-0.1428571428571428*f_rl[7])+0.1428571428571428*f_cr[7]-0.223606797749979*f_rl[1]+0.223606797749979*f_cr[1])*sgn_alphaUpR[8]+((-0.1428571428571428*sgn_alphaUpR[7])-0.223606797749979*sgn_alphaUpR[1])*f_rl[8]+(0.1428571428571428*sgn_alphaUpR[7]+0.223606797749979*sgn_alphaUpR[1])*f_cr[8]+((-0.159719141249985*f_rl[5])+0.159719141249985*f_cr[5]-0.223606797749979*f_rl[4]+0.223606797749979*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[7]+((-0.159719141249985*sgn_alphaUpR[5])-0.223606797749979*sgn_alphaUpR[4]-0.25*sgn_alphaUpR[0]+0.5)*f_rl[7]+(0.159719141249985*sgn_alphaUpR[5]+0.223606797749979*sgn_alphaUpR[4]+0.25*sgn_alphaUpR[0]+0.5)*f_cr[7]+(0.2*f_cr[3]-0.2*f_rl[3])*sgn_alphaUpR[6]+sgn_alphaUpR[3]*(0.2*f_cr[6]-0.2*f_rl[6])+(0.2500000000000001*f_cr[1]-0.2500000000000001*f_rl[1])*sgn_alphaUpR[5]+sgn_alphaUpR[1]*(0.2500000000000001*f_cr[5]-0.2500000000000001*f_rl[5])+(0.223606797749979*f_cr[2]-0.223606797749979*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[3]-0.223606797749979*f_rl[3]); 
  fUpR[8] = ((-0.1020408163265306*f_rl[8])+0.1020408163265306*f_cr[8]-0.159719141249985*f_rl[5]+0.159719141249985*f_cr[5]-0.159719141249985*f_rl[4]+0.159719141249985*f_cr[4]-0.25*f_rl[0]+0.25*f_cr[0])*sgn_alphaUpR[8]+((-0.159719141249985*(sgn_alphaUpR[5]+sgn_alphaUpR[4]))-0.25*sgn_alphaUpR[0]+0.5)*f_rl[8]+(0.159719141249985*(sgn_alphaUpR[5]+sgn_alphaUpR[4])+0.25*sgn_alphaUpR[0]+0.5)*f_cr[8]+((-0.1428571428571428*f_rl[7])+0.1428571428571428*f_cr[7]-0.223606797749979*f_rl[1]+0.223606797749979*f_cr[1])*sgn_alphaUpR[7]+sgn_alphaUpR[1]*(0.223606797749979*f_cr[7]-0.223606797749979*f_rl[7])+((-0.1428571428571428*f_rl[6])+0.1428571428571428*f_cr[6]-0.223606797749979*f_rl[2]+0.223606797749979*f_cr[2])*sgn_alphaUpR[6]+sgn_alphaUpR[2]*(0.223606797749979*f_cr[6]-0.223606797749979*f_rl[6])+(0.25*f_cr[4]-0.25*f_rl[4])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.25*f_cr[5]-0.25*f_rl[5])+(0.2*f_cr[3]-0.2*f_rl[3])*sgn_alphaUpR[3]; 

  } 
  double GhatR[9] = {0.};
  GhatR[0] = 0.5*(alphaR[7]*fUpR[7]+alphaR[5]*fUpR[5]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.447213595499958*alphaR[7]*fUpR[8]+0.5000000000000001*(alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7])+0.447213595499958*alphaR[3]*fUpR[6]+0.4472135954999579*alphaR[1]*fUpR[4]+0.5*(alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.447213595499958*(alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7])+0.4472135954999579*(alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5])+0.5*(alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.4*alphaR[3]*fUpR[8]+0.447213595499958*alphaR[2]*fUpR[7]+0.4*fUpR[6]*alphaR[7]+0.447213595499958*(fUpR[2]*alphaR[7]+alphaR[1]*fUpR[6])+0.4472135954999579*(alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]+alphaR[3]*fUpR[4])+0.5*(alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[4] = 0.5*alphaR[5]*fUpR[8]+0.4472135954999579*alphaR[7]*fUpR[7]+0.5000000000000001*alphaR[2]*fUpR[6]+0.5*alphaR[0]*fUpR[4]+0.4472135954999579*(alphaR[3]*fUpR[3]+alphaR[1]*fUpR[1]); 
  GhatR[5] = 0.31943828249997*alphaR[7]*fUpR[7]+0.5000000000000001*(alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7])+0.31943828249997*alphaR[5]*fUpR[5]+0.5*(alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5])+0.4472135954999579*(alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]); 
  GhatR[6] = 0.447213595499958*alphaR[2]*fUpR[8]+0.4*(alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7])+(0.4472135954999579*alphaR[5]+0.5*alphaR[0])*fUpR[6]+0.5000000000000001*alphaR[2]*fUpR[4]+0.447213595499958*(alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = (0.2857142857142857*alphaR[7]+0.447213595499958*alphaR[1])*fUpR[8]+(0.31943828249997*alphaR[5]+0.5*alphaR[0])*fUpR[7]+(0.31943828249997*fUpR[5]+0.4472135954999579*fUpR[4]+0.5*fUpR[0])*alphaR[7]+0.4*alphaR[3]*fUpR[6]+0.5000000000000001*(alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5])+0.447213595499958*(alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = (0.31943828249997*alphaR[5]+0.5*alphaR[0])*fUpR[8]+0.2857142857142857*alphaR[7]*fUpR[7]+0.447213595499958*(alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6])+0.5*fUpR[4]*alphaR[5]+0.4*alphaR[3]*fUpR[3]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdvx2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvx2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdvx2; 
  out[4] += -1.224744871391589*GhatR[1]*rdvx2; 
  out[5] += -0.7071067811865475*GhatR[3]*rdvx2; 
  out[6] += -1.224744871391589*GhatR[2]*rdvx2; 
  out[7] += -0.7071067811865475*GhatR[4]*rdvx2; 
  out[8] += -1.58113883008419*GhatR[0]*rdvx2; 
  out[9] += -0.7071067811865475*GhatR[5]*rdvx2; 
  out[10] += -1.224744871391589*GhatR[3]*rdvx2; 
  out[11] += -1.224744871391589*GhatR[4]*rdvx2; 
  out[12] += -1.58113883008419*GhatR[1]*rdvx2; 
  out[13] += -0.7071067811865475*GhatR[6]*rdvx2; 
  out[14] += -1.58113883008419*GhatR[2]*rdvx2; 
  out[15] += -0.7071067811865475*GhatR[7]*rdvx2; 
  out[16] += -1.224744871391589*GhatR[5]*rdvx2; 
  out[17] += -1.224744871391589*GhatR[6]*rdvx2; 
  out[18] += -1.58113883008419*GhatR[3]*rdvx2; 
  out[19] += -1.224744871391589*GhatR[7]*rdvx2; 
  out[20] += -1.58113883008419*GhatR[4]*rdvx2; 
  out[21] += -0.7071067811865475*GhatR[8]*rdvx2; 
  out[22] += -1.58113883008419*GhatR[5]*rdvx2; 
  out[23] += -1.58113883008419*GhatR[6]*rdvx2; 
  out[24] += -1.224744871391589*GhatR[8]*rdvx2; 
  out[25] += -1.58113883008419*GhatR[7]*rdvx2; 
  out[26] += -1.58113883008419*GhatR[8]*rdvx2; 

  } else { 

  double fUpL[9] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.58113883008419*fedge[8]+1.224744871391589*fedge[2]+0.7071067811865475*fedge[0]; 
  fUpL[1] = 1.58113883008419*fedge[12]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[1]; 
  fUpL[2] = 1.58113883008419*fedge[14]+1.224744871391589*fedge[6]+0.7071067811865475*fedge[3]; 
  fUpL[3] = 1.58113883008419*fedge[18]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[5]; 
  fUpL[4] = 1.58113883008419*fedge[20]+1.224744871391589*fedge[11]+0.7071067811865475*fedge[7]; 
  fUpL[5] = 1.58113883008419*fedge[22]+1.224744871391589*fedge[16]+0.7071067811865475*fedge[9]; 
  fUpL[6] = 1.58113883008419*fedge[23]+1.224744871391589*fedge[17]+0.7071067811865475*fedge[13]; 
  fUpL[7] = 1.58113883008419*fedge[25]+1.224744871391589*fedge[19]+0.7071067811865475*fedge[15]; 
  fUpL[8] = 1.58113883008419*fedge[26]+1.224744871391589*fedge[24]+0.7071067811865475*fedge[21]; 
    } else { 
  fUpL[0] = 1.58113883008419*fskin[8]-1.224744871391589*fskin[2]+0.7071067811865475*fskin[0]; 
  fUpL[1] = 1.58113883008419*fskin[12]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[1]; 
  fUpL[2] = 1.58113883008419*fskin[14]-1.224744871391589*fskin[6]+0.7071067811865475*fskin[3]; 
  fUpL[3] = 1.58113883008419*fskin[18]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[5]; 
  fUpL[4] = 1.58113883008419*fskin[20]-1.224744871391589*fskin[11]+0.7071067811865475*fskin[7]; 
  fUpL[5] = 1.58113883008419*fskin[22]-1.224744871391589*fskin[16]+0.7071067811865475*fskin[9]; 
  fUpL[6] = 1.58113883008419*fskin[23]-1.224744871391589*fskin[17]+0.7071067811865475*fskin[13]; 
  fUpL[7] = 1.58113883008419*fskin[25]-1.224744871391589*fskin[19]+0.7071067811865475*fskin[15]; 
  fUpL[8] = 1.58113883008419*fskin[26]-1.224744871391589*fskin[24]+0.7071067811865475*fskin[21]; 
    } 
  } else { 
  double f_lr[9] = {0.};
  double f_cl[9] = {0.};
  double sgn_alphaUpL[9] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.58113883008419*fedge[8]+1.224744871391589*fedge[2]+0.7071067811865475*fedge[0]; 
  f_lr[1] = 1.58113883008419*fedge[12]+1.224744871391589*fedge[4]+0.7071067811865475*fedge[1]; 
  f_lr[2] = 1.58113883008419*fedge[14]+1.224744871391589*fedge[6]+0.7071067811865475*fedge[3]; 
  f_lr[3] = 1.58113883008419*fedge[18]+1.224744871391589*fedge[10]+0.7071067811865475*fedge[5]; 
  f_lr[4] = 1.58113883008419*fedge[20]+1.224744871391589*fedge[11]+0.7071067811865475*fedge[7]; 
  f_lr[5] = 1.58113883008419*fedge[22]+1.224744871391589*fedge[16]+0.7071067811865475*fedge[9]; 
  f_lr[6] = 1.58113883008419*fedge[23]+1.224744871391589*fedge[17]+0.7071067811865475*fedge[13]; 
  f_lr[7] = 1.58113883008419*fedge[25]+1.224744871391589*fedge[19]+0.7071067811865475*fedge[15]; 
  f_lr[8] = 1.58113883008419*fedge[26]+1.224744871391589*fedge[24]+0.7071067811865475*fedge[21]; 

  f_cl[0] = 1.58113883008419*fskin[8]-1.224744871391589*fskin[2]+0.7071067811865475*fskin[0]; 
  f_cl[1] = 1.58113883008419*fskin[12]-1.224744871391589*fskin[4]+0.7071067811865475*fskin[1]; 
  f_cl[2] = 1.58113883008419*fskin[14]-1.224744871391589*fskin[6]+0.7071067811865475*fskin[3]; 
  f_cl[3] = 1.58113883008419*fskin[18]-1.224744871391589*fskin[10]+0.7071067811865475*fskin[5]; 
  f_cl[4] = 1.58113883008419*fskin[20]-1.224744871391589*fskin[11]+0.7071067811865475*fskin[7]; 
  f_cl[5] = 1.58113883008419*fskin[22]-1.224744871391589*fskin[16]+0.7071067811865475*fskin[9]; 
  f_cl[6] = 1.58113883008419*fskin[23]-1.224744871391589*fskin[17]+0.7071067811865475*fskin[13]; 
  f_cl[7] = 1.58113883008419*fskin[25]-1.224744871391589*fskin[19]+0.7071067811865475*fskin[15]; 
  f_cl[8] = 1.58113883008419*fskin[26]-1.224744871391589*fskin[24]+0.7071067811865475*fskin[21]; 

  fUpL[0] = (0.25*f_lr[8]-0.25*f_cl[8])*sgn_alphaUpL[8]+(0.25*f_lr[7]-0.25*f_cl[7])*sgn_alphaUpL[7]+(0.25*f_lr[6]-0.25*f_cl[6])*sgn_alphaUpL[6]+(0.25*f_lr[5]-0.25*f_cl[5])*sgn_alphaUpL[5]+(0.25*f_lr[4]-0.25*f_cl[4])*sgn_alphaUpL[4]+(0.25*f_lr[3]-0.25*f_cl[3])*sgn_alphaUpL[3]+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[2]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[1]+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])*sgn_alphaUpL[8]+sgn_alphaUpL[7]*(0.223606797749979*f_lr[8]-0.223606797749979*f_cl[8]+0.2500000000000001*f_lr[5]-0.2500000000000001*f_cl[5])+sgn_alphaUpL[5]*(0.2500000000000001*f_lr[7]-0.2500000000000001*f_cl[7])+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[6]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4])+(0.25*f_lr[2]-0.25*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[1]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])*sgn_alphaUpL[8]+sgn_alphaUpL[6]*(0.223606797749979*f_lr[8]-0.223606797749979*f_cl[8])+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[7]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])+(0.2500000000000001*f_lr[4]-0.2500000000000001*f_cl[4])*sgn_alphaUpL[6]+sgn_alphaUpL[4]*(0.2500000000000001*f_lr[6]-0.2500000000000001*f_cl[6])+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[5]-0.223606797749979*f_cl[5])+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.25*f_lr[3]-0.25*f_cl[3])+(0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[2]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.2*f_lr[3]-0.2*f_cl[3])*sgn_alphaUpL[8]+sgn_alphaUpL[3]*(0.2*f_lr[8]-0.2*f_cl[8])+(0.2*f_lr[6]-0.2*f_cl[6]+0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[7]+(0.2*sgn_alphaUpL[6]+0.223606797749979*sgn_alphaUpL[2])*f_lr[7]+((-0.2*sgn_alphaUpL[6])-0.223606797749979*sgn_alphaUpL[2])*f_cl[7]+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[6]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[5]-0.223606797749979*f_cl[5])+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[3]+(0.25*f_lr[1]-0.25*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.25*f_lr[2]-0.25*f_cl[2]); 
  fUpL[4] = (0.159719141249985*f_lr[8]-0.159719141249985*f_cl[8]+0.25*f_lr[5]-0.25*f_cl[5])*sgn_alphaUpL[8]+sgn_alphaUpL[5]*(0.25*f_lr[8]-0.25*f_cl[8])+(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])*sgn_alphaUpL[7]+(0.159719141249985*f_lr[6]-0.159719141249985*f_cl[6]+0.2500000000000001*f_lr[2]-0.2500000000000001*f_cl[2])*sgn_alphaUpL[6]+sgn_alphaUpL[2]*(0.2500000000000001*f_lr[6]-0.2500000000000001*f_cl[6])+(0.159719141249985*f_lr[4]-0.159719141249985*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[4]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[4]+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[3]+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[1]; 
  fUpL[5] = (0.159719141249985*f_lr[8]-0.159719141249985*f_cl[8]+0.25*f_lr[4]-0.25*f_cl[4])*sgn_alphaUpL[8]+sgn_alphaUpL[4]*(0.25*f_lr[8]-0.25*f_cl[8])+(0.159719141249985*f_lr[7]-0.159719141249985*f_cl[7]+0.2500000000000001*f_lr[1]-0.2500000000000001*f_cl[1])*sgn_alphaUpL[7]+sgn_alphaUpL[1]*(0.2500000000000001*f_lr[7]-0.2500000000000001*f_cl[7])+(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])*sgn_alphaUpL[6]+(0.159719141249985*f_lr[5]-0.159719141249985*f_cl[5]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[5]+(0.25*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.25*sgn_alphaUpL[0])*f_cl[5]+(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3])*sgn_alphaUpL[3]+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[2]; 
  fUpL[6] = (0.1428571428571428*f_lr[6]-0.1428571428571428*f_cl[6]+0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[8]+(0.1428571428571428*sgn_alphaUpL[6]+0.223606797749979*sgn_alphaUpL[2])*f_lr[8]+((-0.1428571428571428*sgn_alphaUpL[6])-0.223606797749979*sgn_alphaUpL[2])*f_cl[8]+(0.2*f_lr[3]-0.2*f_cl[3])*sgn_alphaUpL[7]+sgn_alphaUpL[3]*(0.2*f_lr[7]-0.2*f_cl[7])+(0.223606797749979*f_lr[5]-0.223606797749979*f_cl[5]+0.159719141249985*f_lr[4]-0.159719141249985*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[6]+(0.223606797749979*sgn_alphaUpL[5]+0.159719141249985*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[6]+((-0.223606797749979*sgn_alphaUpL[5])-0.159719141249985*sgn_alphaUpL[4]-0.25*sgn_alphaUpL[0]+0.5)*f_cl[6]+(0.2500000000000001*f_lr[2]-0.2500000000000001*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.2500000000000001*f_lr[4]-0.2500000000000001*f_cl[4])+(0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3]); 
  fUpL[7] = (0.1428571428571428*f_lr[7]-0.1428571428571428*f_cl[7]+0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[8]+(0.1428571428571428*sgn_alphaUpL[7]+0.223606797749979*sgn_alphaUpL[1])*f_lr[8]+((-0.1428571428571428*sgn_alphaUpL[7])-0.223606797749979*sgn_alphaUpL[1])*f_cl[8]+(0.159719141249985*f_lr[5]-0.159719141249985*f_cl[5]+0.223606797749979*f_lr[4]-0.223606797749979*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[7]+(0.159719141249985*sgn_alphaUpL[5]+0.223606797749979*sgn_alphaUpL[4]+0.25*sgn_alphaUpL[0]+0.5)*f_lr[7]+((-0.159719141249985*sgn_alphaUpL[5])-0.223606797749979*sgn_alphaUpL[4]-0.25*sgn_alphaUpL[0]+0.5)*f_cl[7]+(0.2*f_lr[3]-0.2*f_cl[3])*sgn_alphaUpL[6]+sgn_alphaUpL[3]*(0.2*f_lr[6]-0.2*f_cl[6])+(0.2500000000000001*f_lr[1]-0.2500000000000001*f_cl[1])*sgn_alphaUpL[5]+sgn_alphaUpL[1]*(0.2500000000000001*f_lr[5]-0.2500000000000001*f_cl[5])+(0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[3]-0.223606797749979*f_cl[3]); 
  fUpL[8] = (0.1020408163265306*f_lr[8]-0.1020408163265306*f_cl[8]+0.159719141249985*f_lr[5]-0.159719141249985*f_cl[5]+0.159719141249985*f_lr[4]-0.159719141249985*f_cl[4]+0.25*f_lr[0]-0.25*f_cl[0])*sgn_alphaUpL[8]+(0.159719141249985*(sgn_alphaUpL[5]+sgn_alphaUpL[4])+0.25*sgn_alphaUpL[0]+0.5)*f_lr[8]+((-0.159719141249985*(sgn_alphaUpL[5]+sgn_alphaUpL[4]))-0.25*sgn_alphaUpL[0]+0.5)*f_cl[8]+(0.1428571428571428*f_lr[7]-0.1428571428571428*f_cl[7]+0.223606797749979*f_lr[1]-0.223606797749979*f_cl[1])*sgn_alphaUpL[7]+sgn_alphaUpL[1]*(0.223606797749979*f_lr[7]-0.223606797749979*f_cl[7])+(0.1428571428571428*f_lr[6]-0.1428571428571428*f_cl[6]+0.223606797749979*f_lr[2]-0.223606797749979*f_cl[2])*sgn_alphaUpL[6]+sgn_alphaUpL[2]*(0.223606797749979*f_lr[6]-0.223606797749979*f_cl[6])+(0.25*f_lr[4]-0.25*f_cl[4])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.25*f_lr[5]-0.25*f_cl[5])+(0.2*f_lr[3]-0.2*f_cl[3])*sgn_alphaUpL[3]; 

  } 
  double GhatL[9] = {0.};
  GhatL[0] = 0.5*(alphaL[7]*fUpL[7]+alphaL[5]*fUpL[5]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.447213595499958*alphaL[7]*fUpL[8]+0.5000000000000001*(alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7])+0.447213595499958*alphaL[3]*fUpL[6]+0.4472135954999579*alphaL[1]*fUpL[4]+0.5*(alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.447213595499958*(alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7])+0.4472135954999579*(alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5])+0.5*(alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.4*alphaL[3]*fUpL[8]+0.447213595499958*alphaL[2]*fUpL[7]+0.4*fUpL[6]*alphaL[7]+0.447213595499958*(fUpL[2]*alphaL[7]+alphaL[1]*fUpL[6])+0.4472135954999579*(alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]+alphaL[3]*fUpL[4])+0.5*(alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[4] = 0.5*alphaL[5]*fUpL[8]+0.4472135954999579*alphaL[7]*fUpL[7]+0.5000000000000001*alphaL[2]*fUpL[6]+0.5*alphaL[0]*fUpL[4]+0.4472135954999579*(alphaL[3]*fUpL[3]+alphaL[1]*fUpL[1]); 
  GhatL[5] = 0.31943828249997*alphaL[7]*fUpL[7]+0.5000000000000001*(alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7])+0.31943828249997*alphaL[5]*fUpL[5]+0.5*(alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5])+0.4472135954999579*(alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]); 
  GhatL[6] = 0.447213595499958*alphaL[2]*fUpL[8]+0.4*(alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7])+(0.4472135954999579*alphaL[5]+0.5*alphaL[0])*fUpL[6]+0.5000000000000001*alphaL[2]*fUpL[4]+0.447213595499958*(alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = (0.2857142857142857*alphaL[7]+0.447213595499958*alphaL[1])*fUpL[8]+(0.31943828249997*alphaL[5]+0.5*alphaL[0])*fUpL[7]+(0.31943828249997*fUpL[5]+0.4472135954999579*fUpL[4]+0.5*fUpL[0])*alphaL[7]+0.4*alphaL[3]*fUpL[6]+0.5000000000000001*(alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5])+0.447213595499958*(alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = (0.31943828249997*alphaL[5]+0.5*alphaL[0])*fUpL[8]+0.2857142857142857*alphaL[7]*fUpL[7]+0.447213595499958*(alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6])+0.5*fUpL[4]*alphaL[5]+0.4*alphaL[3]*fUpL[3]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdvx2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvx2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdvx2; 
  out[4] += -1.224744871391589*GhatL[1]*rdvx2; 
  out[5] += 0.7071067811865475*GhatL[3]*rdvx2; 
  out[6] += -1.224744871391589*GhatL[2]*rdvx2; 
  out[7] += 0.7071067811865475*GhatL[4]*rdvx2; 
  out[8] += 1.58113883008419*GhatL[0]*rdvx2; 
  out[9] += 0.7071067811865475*GhatL[5]*rdvx2; 
  out[10] += -1.224744871391589*GhatL[3]*rdvx2; 
  out[11] += -1.224744871391589*GhatL[4]*rdvx2; 
  out[12] += 1.58113883008419*GhatL[1]*rdvx2; 
  out[13] += 0.7071067811865475*GhatL[6]*rdvx2; 
  out[14] += 1.58113883008419*GhatL[2]*rdvx2; 
  out[15] += 0.7071067811865475*GhatL[7]*rdvx2; 
  out[16] += -1.224744871391589*GhatL[5]*rdvx2; 
  out[17] += -1.224744871391589*GhatL[6]*rdvx2; 
  out[18] += 1.58113883008419*GhatL[3]*rdvx2; 
  out[19] += -1.224744871391589*GhatL[7]*rdvx2; 
  out[20] += 1.58113883008419*GhatL[4]*rdvx2; 
  out[21] += 0.7071067811865475*GhatL[8]*rdvx2; 
  out[22] += 1.58113883008419*GhatL[5]*rdvx2; 
  out[23] += 1.58113883008419*GhatL[6]*rdvx2; 
  out[24] += -1.224744871391589*GhatL[8]*rdvx2; 
  out[25] += 1.58113883008419*GhatL[7]*rdvx2; 
  out[26] += 1.58113883008419*GhatL[8]*rdvx2; 

  } 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 1.25*rdvx2*cflFreq; 

} 
