#include <gkyl_rad_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:     cell-center coordinates. 
  // dxv[4]:   cell spacing. 
  // nI:        ion density. 
  // nuField:       2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  double rdv2 = 2.0/dxv[2]; 
  double alphaDrSurf[8] = {0.0}; 
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.1178511301977579*(6.708203932499369*nI[3]*nuField[20]+6.708203932499369*(nI[2]*nuField[18]+nI[1]*nuField[17])+6.708203932499369*nI[0]*nuField[16]+5.196152422706631*(nI[3]*nuField[11]+nI[2]*nuField[7]+nI[1]*nuField[6])+3.0*nI[3]*nuField[5]+5.196152422706631*nI[0]*nuField[3]+3.0*(nI[2]*nuField[2]+nI[1]*nuField[1]+nI[0]*nuField[0])); 
  alphaDrSurf[1] = 0.1178511301977579*(6.708203932499369*nI[2]*nuField[20]+6.708203932499369*(nI[3]*nuField[18]+nI[0]*nuField[17])+6.708203932499369*nI[1]*nuField[16]+5.196152422706631*(nI[2]*nuField[11]+nI[3]*nuField[7]+nI[0]*nuField[6])+3.0*nI[2]*nuField[5]+5.196152422706631*nI[1]*nuField[3]+3.0*(nuField[2]*nI[3]+nI[0]*nuField[1]+nuField[0]*nI[1])); 
  alphaDrSurf[2] = 0.1178511301977579*(6.708203932499369*nI[1]*nuField[20]+6.708203932499369*(nI[0]*nuField[18]+nI[3]*nuField[17])+6.708203932499369*nI[2]*nuField[16]+5.196152422706631*(nI[1]*nuField[11]+nI[0]*nuField[7]+nI[3]*nuField[6])+3.0*nI[1]*nuField[5]+5.196152422706631*nI[2]*nuField[3]+3.0*(nuField[1]*nI[3]+nI[0]*nuField[2]+nuField[0]*nI[2])); 
  alphaDrSurf[3] = 0.1178511301977579*(6.708203932499369*nI[3]*nuField[23]+6.708203932499369*(nI[2]*nuField[22]+nI[1]*nuField[21])+6.708203932499369*nI[0]*nuField[19]+5.196152422706631*(nI[3]*nuField[15]+nI[2]*nuField[14]+nI[1]*nuField[13])+3.0*nI[3]*nuField[12]+5.196152422706631*nI[0]*nuField[10]+3.0*(nI[2]*nuField[9]+nI[1]*nuField[8]+nI[0]*nuField[4])); 
  alphaDrSurf[4] = 0.1178511301977579*(6.708203932499369*nI[0]*nuField[20]+6.708203932499369*(nI[1]*nuField[18]+nI[2]*nuField[17])+6.708203932499369*nI[3]*nuField[16]+5.196152422706631*(nI[0]*nuField[11]+nI[1]*nuField[7]+nI[2]*nuField[6])+3.0*nI[0]*nuField[5]+5.196152422706631*nI[3]*nuField[3]+3.0*(nuField[0]*nI[3]+nI[1]*nuField[2]+nuField[1]*nI[2])); 
  alphaDrSurf[5] = 0.1178511301977579*(6.708203932499369*nI[2]*nuField[23]+6.708203932499369*(nI[3]*nuField[22]+nI[0]*nuField[21])+6.708203932499369*nI[1]*nuField[19]+5.196152422706631*(nI[2]*nuField[15]+nI[3]*nuField[14]+nI[0]*nuField[13])+3.0*nI[2]*nuField[12]+5.196152422706631*nI[1]*nuField[10]+3.0*(nI[3]*nuField[9]+nI[0]*nuField[8]+nI[1]*nuField[4])); 
  alphaDrSurf[6] = 0.1178511301977579*(6.708203932499369*nI[1]*nuField[23]+6.708203932499369*(nI[0]*nuField[22]+nI[3]*nuField[21])+6.708203932499369*nI[2]*nuField[19]+5.196152422706631*(nI[1]*nuField[15]+nI[0]*nuField[14]+nI[3]*nuField[13])+3.0*nI[1]*nuField[12]+5.196152422706631*nI[2]*nuField[10]+3.0*(nI[0]*nuField[9]+nI[3]*nuField[8]+nI[2]*nuField[4])); 
  alphaDrSurf[7] = 0.1178511301977579*(6.708203932499369*nI[0]*nuField[23]+6.708203932499369*(nI[1]*nuField[22]+nI[2]*nuField[21])+6.708203932499369*nI[3]*nuField[19]+5.196152422706631*(nI[0]*nuField[15]+nI[1]*nuField[14]+nI[2]*nuField[13])+3.0*nI[0]*nuField[12]+5.196152422706631*nI[3]*nuField[10]+3.0*(nI[1]*nuField[9]+nI[2]*nuField[8]+nI[3]*nuField[4])); 

  Ghat[0] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[7]*fEdge[23]+6.708203932499369*(alphaDrSurf[6]*fEdge[22]+alphaDrSurf[5]*fEdge[21]+alphaDrSurf[4]*fEdge[20])+6.708203932499369*(alphaDrSurf[3]*fEdge[19]+alphaDrSurf[2]*fEdge[18]+alphaDrSurf[1]*fEdge[17])+6.708203932499369*alphaDrSurf[0]*fEdge[16]-5.196152422706631*(alphaDrSurf[7]*fEdge[15]+alphaDrSurf[6]*fEdge[14]+alphaDrSurf[5]*fEdge[13])+3.0*alphaDrSurf[7]*fEdge[12]-5.196152422706631*(alphaDrSurf[4]*fEdge[11]+alphaDrSurf[3]*fEdge[10])+3.0*(alphaDrSurf[6]*fEdge[9]+alphaDrSurf[5]*fEdge[8])-5.196152422706631*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[1]*fEdge[6])+3.0*(alphaDrSurf[4]*fEdge[5]+alphaDrSurf[3]*fEdge[4])-5.196152422706631*alphaDrSurf[0]*fEdge[3]+3.0*(alphaDrSurf[2]*fEdge[2]+alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[6]*fEdge[23]+6.708203932499369*(alphaDrSurf[7]*fEdge[22]+alphaDrSurf[3]*fEdge[21]+alphaDrSurf[2]*fEdge[20])+6.708203932499369*(alphaDrSurf[5]*fEdge[19]+alphaDrSurf[4]*fEdge[18]+alphaDrSurf[0]*fEdge[17])+6.708203932499369*alphaDrSurf[1]*fEdge[16]-5.196152422706631*(alphaDrSurf[6]*fEdge[15]+alphaDrSurf[7]*fEdge[14]+alphaDrSurf[3]*fEdge[13])+3.0*alphaDrSurf[6]*fEdge[12]-5.196152422706631*(alphaDrSurf[2]*fEdge[11]+alphaDrSurf[5]*fEdge[10])+3.0*(alphaDrSurf[7]*fEdge[9]+alphaDrSurf[3]*fEdge[8])-5.196152422706631*(alphaDrSurf[4]*fEdge[7]+alphaDrSurf[0]*fEdge[6])+3.0*(alphaDrSurf[2]*fEdge[5]+fEdge[4]*alphaDrSurf[5]+fEdge[2]*alphaDrSurf[4])-5.196152422706631*alphaDrSurf[1]*fEdge[3]+3.0*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 
  Ghat[2] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[5]*fEdge[23]+6.708203932499369*(alphaDrSurf[3]*fEdge[22]+alphaDrSurf[7]*fEdge[21]+alphaDrSurf[1]*fEdge[20])+6.708203932499369*(alphaDrSurf[6]*fEdge[19]+alphaDrSurf[0]*fEdge[18]+alphaDrSurf[4]*fEdge[17])+6.708203932499369*alphaDrSurf[2]*fEdge[16]-5.196152422706631*(alphaDrSurf[5]*fEdge[15]+alphaDrSurf[3]*fEdge[14]+alphaDrSurf[7]*fEdge[13])+3.0*alphaDrSurf[5]*fEdge[12]-5.196152422706631*(alphaDrSurf[1]*fEdge[11]+alphaDrSurf[6]*fEdge[10])+3.0*(alphaDrSurf[3]*fEdge[9]+alphaDrSurf[7]*fEdge[8])-5.196152422706631*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[4]*fEdge[6])+3.0*(fEdge[4]*alphaDrSurf[6]+alphaDrSurf[1]*fEdge[5]+fEdge[1]*alphaDrSurf[4])-5.196152422706631*alphaDrSurf[2]*fEdge[3]+3.0*(alphaDrSurf[0]*fEdge[2]+fEdge[0]*alphaDrSurf[2])); 
  Ghat[3] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[4]*fEdge[23]+6.708203932499369*(alphaDrSurf[2]*fEdge[22]+alphaDrSurf[1]*fEdge[21]+alphaDrSurf[7]*fEdge[20])+6.708203932499369*(alphaDrSurf[0]*fEdge[19]+alphaDrSurf[6]*fEdge[18]+alphaDrSurf[5]*fEdge[17])+6.708203932499369*alphaDrSurf[3]*fEdge[16]-5.196152422706631*(alphaDrSurf[4]*fEdge[15]+alphaDrSurf[2]*fEdge[14]+alphaDrSurf[1]*fEdge[13])+3.0*alphaDrSurf[4]*fEdge[12]-5.196152422706631*(alphaDrSurf[7]*fEdge[11]+alphaDrSurf[0]*fEdge[10])+3.0*(alphaDrSurf[2]*fEdge[9]+alphaDrSurf[1]*fEdge[8])-5.196152422706631*alphaDrSurf[6]*fEdge[7]+3.0*fEdge[5]*alphaDrSurf[7]-5.196152422706631*alphaDrSurf[5]*fEdge[6]+3.0*(fEdge[2]*alphaDrSurf[6]+fEdge[1]*alphaDrSurf[5]+alphaDrSurf[0]*fEdge[4])+alphaDrSurf[3]*(3.0*fEdge[0]-5.196152422706631*fEdge[3])); 
  Ghat[4] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[3]*fEdge[23]+6.708203932499369*(alphaDrSurf[5]*fEdge[22]+alphaDrSurf[6]*fEdge[21]+alphaDrSurf[0]*fEdge[20])+6.708203932499369*(alphaDrSurf[7]*fEdge[19]+alphaDrSurf[1]*fEdge[18]+alphaDrSurf[2]*fEdge[17])+6.708203932499369*alphaDrSurf[4]*fEdge[16]-5.196152422706631*(alphaDrSurf[3]*fEdge[15]+alphaDrSurf[5]*fEdge[14]+alphaDrSurf[6]*fEdge[13])+3.0*alphaDrSurf[3]*fEdge[12]-5.196152422706631*(alphaDrSurf[0]*fEdge[11]+alphaDrSurf[7]*fEdge[10])+3.0*(alphaDrSurf[5]*fEdge[9]+alphaDrSurf[6]*fEdge[8])-5.196152422706631*alphaDrSurf[1]*fEdge[7]+3.0*fEdge[4]*alphaDrSurf[7]-5.196152422706631*alphaDrSurf[2]*fEdge[6]+3.0*alphaDrSurf[0]*fEdge[5]+(3.0*fEdge[0]-5.196152422706631*fEdge[3])*alphaDrSurf[4]+3.0*(alphaDrSurf[1]*fEdge[2]+fEdge[1]*alphaDrSurf[2])); 
  Ghat[5] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[2]*fEdge[23]+6.708203932499369*(alphaDrSurf[4]*fEdge[22]+alphaDrSurf[0]*fEdge[21]+alphaDrSurf[6]*fEdge[20])+6.708203932499369*(alphaDrSurf[1]*fEdge[19]+alphaDrSurf[7]*fEdge[18]+alphaDrSurf[3]*fEdge[17])+6.708203932499369*alphaDrSurf[5]*fEdge[16]-5.196152422706631*(alphaDrSurf[2]*fEdge[15]+alphaDrSurf[4]*fEdge[14]+alphaDrSurf[0]*fEdge[13])+3.0*alphaDrSurf[2]*fEdge[12]-5.196152422706631*(alphaDrSurf[6]*fEdge[11]+alphaDrSurf[1]*fEdge[10])+3.0*(alphaDrSurf[4]*fEdge[9]+alphaDrSurf[0]*fEdge[8])+alphaDrSurf[7]*(3.0*fEdge[2]-5.196152422706631*fEdge[7])-5.196152422706631*alphaDrSurf[3]*fEdge[6]+3.0*fEdge[5]*alphaDrSurf[6]+(3.0*fEdge[0]-5.196152422706631*fEdge[3])*alphaDrSurf[5]+3.0*(alphaDrSurf[1]*fEdge[4]+fEdge[1]*alphaDrSurf[3])); 
  Ghat[6] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[1]*fEdge[23]+6.708203932499369*(alphaDrSurf[0]*fEdge[22]+alphaDrSurf[4]*fEdge[21]+alphaDrSurf[5]*fEdge[20])+6.708203932499369*(alphaDrSurf[2]*fEdge[19]+alphaDrSurf[3]*fEdge[18]+alphaDrSurf[7]*fEdge[17])+6.708203932499369*alphaDrSurf[6]*fEdge[16]-5.196152422706631*(alphaDrSurf[1]*fEdge[15]+alphaDrSurf[0]*fEdge[14]+alphaDrSurf[4]*fEdge[13])+3.0*alphaDrSurf[1]*fEdge[12]-5.196152422706631*(alphaDrSurf[5]*fEdge[11]+alphaDrSurf[2]*fEdge[10])+3.0*(alphaDrSurf[0]*fEdge[9]+alphaDrSurf[4]*fEdge[8])-5.196152422706631*alphaDrSurf[3]*fEdge[7]+(3.0*fEdge[1]-5.196152422706631*fEdge[6])*alphaDrSurf[7]+(3.0*fEdge[0]-5.196152422706631*fEdge[3])*alphaDrSurf[6]+3.0*(alphaDrSurf[5]*fEdge[5]+alphaDrSurf[2]*fEdge[4]+fEdge[2]*alphaDrSurf[3])); 
  Ghat[7] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[0]*fEdge[23]+6.708203932499369*(alphaDrSurf[1]*fEdge[22]+alphaDrSurf[2]*fEdge[21]+alphaDrSurf[3]*fEdge[20])+6.708203932499369*(alphaDrSurf[4]*fEdge[19]+alphaDrSurf[5]*fEdge[18]+alphaDrSurf[6]*fEdge[17])+6.708203932499369*alphaDrSurf[7]*fEdge[16]-5.196152422706631*(alphaDrSurf[0]*fEdge[15]+alphaDrSurf[1]*fEdge[14]+alphaDrSurf[2]*fEdge[13])+3.0*alphaDrSurf[0]*fEdge[12]-5.196152422706631*(alphaDrSurf[3]*fEdge[11]+alphaDrSurf[4]*fEdge[10])+3.0*(alphaDrSurf[1]*fEdge[9]+alphaDrSurf[2]*fEdge[8])-5.196152422706631*alphaDrSurf[5]*fEdge[7]+(3.0*fEdge[0]-5.196152422706631*fEdge[3])*alphaDrSurf[7]-5.196152422706631*alphaDrSurf[6]*fEdge[6]+3.0*(fEdge[1]*alphaDrSurf[6]+alphaDrSurf[3]*fEdge[5]+fEdge[2]*alphaDrSurf[5]+alphaDrSurf[4]*fEdge[4])); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 1.224744871391589*Ghat[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 1.224744871391589*Ghat[4]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 
  out[16] += 1.58113883008419*Ghat[0]*rdv2; 
  out[17] += 1.58113883008419*Ghat[1]*rdv2; 
  out[18] += 1.58113883008419*Ghat[2]*rdv2; 
  out[19] += 1.58113883008419*Ghat[3]*rdv2; 
  out[20] += 1.58113883008419*Ghat[4]*rdv2; 
  out[21] += 1.58113883008419*Ghat[5]*rdv2; 
  out[22] += 1.58113883008419*Ghat[6]*rdv2; 
  out[23] += 1.58113883008419*Ghat[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.1178511301977579*(6.708203932499369*nI[3]*nuField[20]+6.708203932499369*(nI[2]*nuField[18]+nI[1]*nuField[17])+6.708203932499369*nI[0]*nuField[16]-5.196152422706631*(nI[3]*nuField[11]+nI[2]*nuField[7]+nI[1]*nuField[6])+3.0*nI[3]*nuField[5]-5.196152422706631*nI[0]*nuField[3]+3.0*(nI[2]*nuField[2]+nI[1]*nuField[1]+nI[0]*nuField[0])); 
  alphaDrSurf[1] = 0.1178511301977579*(6.708203932499369*nI[2]*nuField[20]+6.708203932499369*(nI[3]*nuField[18]+nI[0]*nuField[17])+6.708203932499369*nI[1]*nuField[16]-5.196152422706631*(nI[2]*nuField[11]+nI[3]*nuField[7]+nI[0]*nuField[6])+3.0*nI[2]*nuField[5]-5.196152422706631*nI[1]*nuField[3]+3.0*(nuField[2]*nI[3]+nI[0]*nuField[1]+nuField[0]*nI[1])); 
  alphaDrSurf[2] = 0.1178511301977579*(6.708203932499369*nI[1]*nuField[20]+6.708203932499369*(nI[0]*nuField[18]+nI[3]*nuField[17])+6.708203932499369*nI[2]*nuField[16]-5.196152422706631*(nI[1]*nuField[11]+nI[0]*nuField[7]+nI[3]*nuField[6])+3.0*nI[1]*nuField[5]-5.196152422706631*nI[2]*nuField[3]+3.0*(nuField[1]*nI[3]+nI[0]*nuField[2]+nuField[0]*nI[2])); 
  alphaDrSurf[3] = 0.1178511301977579*(6.708203932499369*nI[3]*nuField[23]+6.708203932499369*(nI[2]*nuField[22]+nI[1]*nuField[21])+6.708203932499369*nI[0]*nuField[19]-5.196152422706631*(nI[3]*nuField[15]+nI[2]*nuField[14]+nI[1]*nuField[13])+3.0*nI[3]*nuField[12]-5.196152422706631*nI[0]*nuField[10]+3.0*(nI[2]*nuField[9]+nI[1]*nuField[8]+nI[0]*nuField[4])); 
  alphaDrSurf[4] = 0.1178511301977579*(6.708203932499369*nI[0]*nuField[20]+6.708203932499369*(nI[1]*nuField[18]+nI[2]*nuField[17])+6.708203932499369*nI[3]*nuField[16]-5.196152422706631*(nI[0]*nuField[11]+nI[1]*nuField[7]+nI[2]*nuField[6])+3.0*nI[0]*nuField[5]-5.196152422706631*nI[3]*nuField[3]+3.0*(nuField[0]*nI[3]+nI[1]*nuField[2]+nuField[1]*nI[2])); 
  alphaDrSurf[5] = 0.1178511301977579*(6.708203932499369*nI[2]*nuField[23]+6.708203932499369*(nI[3]*nuField[22]+nI[0]*nuField[21])+6.708203932499369*nI[1]*nuField[19]-5.196152422706631*(nI[2]*nuField[15]+nI[3]*nuField[14]+nI[0]*nuField[13])+3.0*nI[2]*nuField[12]-5.196152422706631*nI[1]*nuField[10]+3.0*(nI[3]*nuField[9]+nI[0]*nuField[8]+nI[1]*nuField[4])); 
  alphaDrSurf[6] = 0.1178511301977579*(6.708203932499369*nI[1]*nuField[23]+6.708203932499369*(nI[0]*nuField[22]+nI[3]*nuField[21])+6.708203932499369*nI[2]*nuField[19]-5.196152422706631*(nI[1]*nuField[15]+nI[0]*nuField[14]+nI[3]*nuField[13])+3.0*nI[1]*nuField[12]-5.196152422706631*nI[2]*nuField[10]+3.0*(nI[0]*nuField[9]+nI[3]*nuField[8]+nI[2]*nuField[4])); 
  alphaDrSurf[7] = 0.1178511301977579*(6.708203932499369*nI[0]*nuField[23]+6.708203932499369*(nI[1]*nuField[22]+nI[2]*nuField[21])+6.708203932499369*nI[3]*nuField[19]-5.196152422706631*(nI[0]*nuField[15]+nI[1]*nuField[14]+nI[2]*nuField[13])+3.0*nI[0]*nuField[12]-5.196152422706631*nI[3]*nuField[10]+3.0*(nI[1]*nuField[9]+nI[2]*nuField[8]+nI[3]*nuField[4])); 

  Ghat[0] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[7]*fSkin[23]+6.708203932499369*(alphaDrSurf[6]*fSkin[22]+alphaDrSurf[5]*fSkin[21]+alphaDrSurf[4]*fSkin[20])+6.708203932499369*(alphaDrSurf[3]*fSkin[19]+alphaDrSurf[2]*fSkin[18]+alphaDrSurf[1]*fSkin[17])+6.708203932499369*alphaDrSurf[0]*fSkin[16]-5.196152422706631*(alphaDrSurf[7]*fSkin[15]+alphaDrSurf[6]*fSkin[14]+alphaDrSurf[5]*fSkin[13])+3.0*alphaDrSurf[7]*fSkin[12]-5.196152422706631*(alphaDrSurf[4]*fSkin[11]+alphaDrSurf[3]*fSkin[10])+3.0*(alphaDrSurf[6]*fSkin[9]+alphaDrSurf[5]*fSkin[8])-5.196152422706631*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[1]*fSkin[6])+3.0*(alphaDrSurf[4]*fSkin[5]+alphaDrSurf[3]*fSkin[4])-5.196152422706631*alphaDrSurf[0]*fSkin[3]+3.0*(alphaDrSurf[2]*fSkin[2]+alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[6]*fSkin[23]+6.708203932499369*(alphaDrSurf[7]*fSkin[22]+alphaDrSurf[3]*fSkin[21]+alphaDrSurf[2]*fSkin[20])+6.708203932499369*(alphaDrSurf[5]*fSkin[19]+alphaDrSurf[4]*fSkin[18]+alphaDrSurf[0]*fSkin[17])+6.708203932499369*alphaDrSurf[1]*fSkin[16]-5.196152422706631*(alphaDrSurf[6]*fSkin[15]+alphaDrSurf[7]*fSkin[14]+alphaDrSurf[3]*fSkin[13])+3.0*alphaDrSurf[6]*fSkin[12]-5.196152422706631*(alphaDrSurf[2]*fSkin[11]+alphaDrSurf[5]*fSkin[10])+3.0*(alphaDrSurf[7]*fSkin[9]+alphaDrSurf[3]*fSkin[8])-5.196152422706631*(alphaDrSurf[4]*fSkin[7]+alphaDrSurf[0]*fSkin[6])+3.0*(alphaDrSurf[2]*fSkin[5]+fSkin[4]*alphaDrSurf[5]+fSkin[2]*alphaDrSurf[4])-5.196152422706631*alphaDrSurf[1]*fSkin[3]+3.0*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 
  Ghat[2] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[5]*fSkin[23]+6.708203932499369*(alphaDrSurf[3]*fSkin[22]+alphaDrSurf[7]*fSkin[21]+alphaDrSurf[1]*fSkin[20])+6.708203932499369*(alphaDrSurf[6]*fSkin[19]+alphaDrSurf[0]*fSkin[18]+alphaDrSurf[4]*fSkin[17])+6.708203932499369*alphaDrSurf[2]*fSkin[16]-5.196152422706631*(alphaDrSurf[5]*fSkin[15]+alphaDrSurf[3]*fSkin[14]+alphaDrSurf[7]*fSkin[13])+3.0*alphaDrSurf[5]*fSkin[12]-5.196152422706631*(alphaDrSurf[1]*fSkin[11]+alphaDrSurf[6]*fSkin[10])+3.0*(alphaDrSurf[3]*fSkin[9]+alphaDrSurf[7]*fSkin[8])-5.196152422706631*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[4]*fSkin[6])+3.0*(fSkin[4]*alphaDrSurf[6]+alphaDrSurf[1]*fSkin[5]+fSkin[1]*alphaDrSurf[4])-5.196152422706631*alphaDrSurf[2]*fSkin[3]+3.0*(alphaDrSurf[0]*fSkin[2]+fSkin[0]*alphaDrSurf[2])); 
  Ghat[3] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[4]*fSkin[23]+6.708203932499369*(alphaDrSurf[2]*fSkin[22]+alphaDrSurf[1]*fSkin[21]+alphaDrSurf[7]*fSkin[20])+6.708203932499369*(alphaDrSurf[0]*fSkin[19]+alphaDrSurf[6]*fSkin[18]+alphaDrSurf[5]*fSkin[17])+6.708203932499369*alphaDrSurf[3]*fSkin[16]-5.196152422706631*(alphaDrSurf[4]*fSkin[15]+alphaDrSurf[2]*fSkin[14]+alphaDrSurf[1]*fSkin[13])+3.0*alphaDrSurf[4]*fSkin[12]-5.196152422706631*(alphaDrSurf[7]*fSkin[11]+alphaDrSurf[0]*fSkin[10])+3.0*(alphaDrSurf[2]*fSkin[9]+alphaDrSurf[1]*fSkin[8])-5.196152422706631*alphaDrSurf[6]*fSkin[7]+3.0*fSkin[5]*alphaDrSurf[7]-5.196152422706631*alphaDrSurf[5]*fSkin[6]+3.0*(fSkin[2]*alphaDrSurf[6]+fSkin[1]*alphaDrSurf[5]+alphaDrSurf[0]*fSkin[4])+alphaDrSurf[3]*(3.0*fSkin[0]-5.196152422706631*fSkin[3])); 
  Ghat[4] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[3]*fSkin[23]+6.708203932499369*(alphaDrSurf[5]*fSkin[22]+alphaDrSurf[6]*fSkin[21]+alphaDrSurf[0]*fSkin[20])+6.708203932499369*(alphaDrSurf[7]*fSkin[19]+alphaDrSurf[1]*fSkin[18]+alphaDrSurf[2]*fSkin[17])+6.708203932499369*alphaDrSurf[4]*fSkin[16]-5.196152422706631*(alphaDrSurf[3]*fSkin[15]+alphaDrSurf[5]*fSkin[14]+alphaDrSurf[6]*fSkin[13])+3.0*alphaDrSurf[3]*fSkin[12]-5.196152422706631*(alphaDrSurf[0]*fSkin[11]+alphaDrSurf[7]*fSkin[10])+3.0*(alphaDrSurf[5]*fSkin[9]+alphaDrSurf[6]*fSkin[8])-5.196152422706631*alphaDrSurf[1]*fSkin[7]+3.0*fSkin[4]*alphaDrSurf[7]-5.196152422706631*alphaDrSurf[2]*fSkin[6]+3.0*alphaDrSurf[0]*fSkin[5]+(3.0*fSkin[0]-5.196152422706631*fSkin[3])*alphaDrSurf[4]+3.0*(alphaDrSurf[1]*fSkin[2]+fSkin[1]*alphaDrSurf[2])); 
  Ghat[5] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[2]*fSkin[23]+6.708203932499369*(alphaDrSurf[4]*fSkin[22]+alphaDrSurf[0]*fSkin[21]+alphaDrSurf[6]*fSkin[20])+6.708203932499369*(alphaDrSurf[1]*fSkin[19]+alphaDrSurf[7]*fSkin[18]+alphaDrSurf[3]*fSkin[17])+6.708203932499369*alphaDrSurf[5]*fSkin[16]-5.196152422706631*(alphaDrSurf[2]*fSkin[15]+alphaDrSurf[4]*fSkin[14]+alphaDrSurf[0]*fSkin[13])+3.0*alphaDrSurf[2]*fSkin[12]-5.196152422706631*(alphaDrSurf[6]*fSkin[11]+alphaDrSurf[1]*fSkin[10])+3.0*(alphaDrSurf[4]*fSkin[9]+alphaDrSurf[0]*fSkin[8])+alphaDrSurf[7]*(3.0*fSkin[2]-5.196152422706631*fSkin[7])-5.196152422706631*alphaDrSurf[3]*fSkin[6]+3.0*fSkin[5]*alphaDrSurf[6]+(3.0*fSkin[0]-5.196152422706631*fSkin[3])*alphaDrSurf[5]+3.0*(alphaDrSurf[1]*fSkin[4]+fSkin[1]*alphaDrSurf[3])); 
  Ghat[6] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[1]*fSkin[23]+6.708203932499369*(alphaDrSurf[0]*fSkin[22]+alphaDrSurf[4]*fSkin[21]+alphaDrSurf[5]*fSkin[20])+6.708203932499369*(alphaDrSurf[2]*fSkin[19]+alphaDrSurf[3]*fSkin[18]+alphaDrSurf[7]*fSkin[17])+6.708203932499369*alphaDrSurf[6]*fSkin[16]-5.196152422706631*(alphaDrSurf[1]*fSkin[15]+alphaDrSurf[0]*fSkin[14]+alphaDrSurf[4]*fSkin[13])+3.0*alphaDrSurf[1]*fSkin[12]-5.196152422706631*(alphaDrSurf[5]*fSkin[11]+alphaDrSurf[2]*fSkin[10])+3.0*(alphaDrSurf[0]*fSkin[9]+alphaDrSurf[4]*fSkin[8])-5.196152422706631*alphaDrSurf[3]*fSkin[7]+(3.0*fSkin[1]-5.196152422706631*fSkin[6])*alphaDrSurf[7]+(3.0*fSkin[0]-5.196152422706631*fSkin[3])*alphaDrSurf[6]+3.0*(alphaDrSurf[5]*fSkin[5]+alphaDrSurf[2]*fSkin[4]+fSkin[2]*alphaDrSurf[3])); 
  Ghat[7] = 0.08333333333333333*(6.708203932499369*alphaDrSurf[0]*fSkin[23]+6.708203932499369*(alphaDrSurf[1]*fSkin[22]+alphaDrSurf[2]*fSkin[21]+alphaDrSurf[3]*fSkin[20])+6.708203932499369*(alphaDrSurf[4]*fSkin[19]+alphaDrSurf[5]*fSkin[18]+alphaDrSurf[6]*fSkin[17])+6.708203932499369*alphaDrSurf[7]*fSkin[16]-5.196152422706631*(alphaDrSurf[0]*fSkin[15]+alphaDrSurf[1]*fSkin[14]+alphaDrSurf[2]*fSkin[13])+3.0*alphaDrSurf[0]*fSkin[12]-5.196152422706631*(alphaDrSurf[3]*fSkin[11]+alphaDrSurf[4]*fSkin[10])+3.0*(alphaDrSurf[1]*fSkin[9]+alphaDrSurf[2]*fSkin[8])-5.196152422706631*alphaDrSurf[5]*fSkin[7]+(3.0*fSkin[0]-5.196152422706631*fSkin[3])*alphaDrSurf[7]-5.196152422706631*alphaDrSurf[6]*fSkin[6]+3.0*(fSkin[1]*alphaDrSurf[6]+alphaDrSurf[3]*fSkin[5]+fSkin[2]*alphaDrSurf[5]+alphaDrSurf[4]*fSkin[4])); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 1.224744871391589*Ghat[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat[2]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 1.224744871391589*Ghat[4]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 
  out[16] += -1.58113883008419*Ghat[0]*rdv2; 
  out[17] += -1.58113883008419*Ghat[1]*rdv2; 
  out[18] += -1.58113883008419*Ghat[2]*rdv2; 
  out[19] += -1.58113883008419*Ghat[3]*rdv2; 
  out[20] += -1.58113883008419*Ghat[4]*rdv2; 
  out[21] += -1.58113883008419*Ghat[5]*rdv2; 
  out[22] += -1.58113883008419*Ghat[6]*rdv2; 
  out[23] += -1.58113883008419*Ghat[7]*rdv2; 

  } 

  return 0.;

} 
