#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, const double *nu, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nu:         collisionality. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 

  double alphaDrSurf[8] = {0.0}; 
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 
  alphaDrSurf[3] = nu[3]*wvpar+0.5*nu[3]*dvpar; 
  alphaDrSurf[4] = nu[4]*wvpar+0.5*nu[4]*dvpar; 
  alphaDrSurf[5] = nu[5]*wvpar+0.5*nu[5]*dvpar; 
  alphaDrSurf[6] = nu[6]*wvpar+0.5*nu[6]*dvpar; 
  alphaDrSurf[7] = nu[7]*wvpar+0.5*nu[7]*dvpar; 

  Ghat[0] = 0.5590169943749476*alphaDrSurf[7]*fSkin[23]+0.5590169943749475*alphaDrSurf[6]*fSkin[22]+0.5590169943749475*alphaDrSurf[5]*fSkin[21]+0.5590169943749475*alphaDrSurf[4]*fSkin[20]+0.5590169943749476*alphaDrSurf[3]*fSkin[19]+0.5590169943749476*alphaDrSurf[2]*fSkin[18]+0.5590169943749476*alphaDrSurf[1]*fSkin[17]+0.5590169943749475*alphaDrSurf[0]*fSkin[16]+0.4330127018922193*alphaDrSurf[7]*fSkin[15]+0.4330127018922193*alphaDrSurf[6]*fSkin[14]+0.4330127018922193*alphaDrSurf[5]*fSkin[13]+0.4330127018922193*alphaDrSurf[4]*fSkin[12]+0.25*alphaDrSurf[7]*fSkin[11]+0.4330127018922193*alphaDrSurf[3]*fSkin[10]+0.4330127018922193*alphaDrSurf[2]*fSkin[9]+0.4330127018922193*alphaDrSurf[1]*fSkin[8]+0.25*alphaDrSurf[6]*fSkin[7]+0.25*alphaDrSurf[5]*fSkin[6]+0.25*alphaDrSurf[4]*fSkin[5]+0.4330127018922193*alphaDrSurf[0]*fSkin[4]+0.25*alphaDrSurf[3]*fSkin[3]+0.25*alphaDrSurf[2]*fSkin[2]+0.25*alphaDrSurf[1]*fSkin[1]+0.25*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 0.5590169943749476*alphaDrSurf[6]*fSkin[23]+0.5590169943749475*alphaDrSurf[7]*fSkin[22]+0.5590169943749475*alphaDrSurf[3]*fSkin[21]+0.5590169943749475*alphaDrSurf[2]*fSkin[20]+0.5590169943749476*alphaDrSurf[5]*fSkin[19]+0.5590169943749476*alphaDrSurf[4]*fSkin[18]+0.5590169943749476*alphaDrSurf[0]*fSkin[17]+0.5590169943749475*alphaDrSurf[1]*fSkin[16]+0.4330127018922193*alphaDrSurf[6]*fSkin[15]+0.4330127018922193*alphaDrSurf[7]*fSkin[14]+0.4330127018922193*alphaDrSurf[3]*fSkin[13]+0.4330127018922193*alphaDrSurf[2]*fSkin[12]+0.25*alphaDrSurf[6]*fSkin[11]+0.4330127018922193*alphaDrSurf[5]*fSkin[10]+0.4330127018922193*alphaDrSurf[4]*fSkin[9]+0.4330127018922193*alphaDrSurf[0]*fSkin[8]+0.25*alphaDrSurf[7]*fSkin[7]+0.25*alphaDrSurf[3]*fSkin[6]+0.25*alphaDrSurf[2]*fSkin[5]+0.25*fSkin[3]*alphaDrSurf[5]+0.4330127018922193*alphaDrSurf[1]*fSkin[4]+0.25*fSkin[2]*alphaDrSurf[4]+0.25*alphaDrSurf[0]*fSkin[1]+0.25*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.5590169943749476*alphaDrSurf[5]*fSkin[23]+0.5590169943749475*alphaDrSurf[3]*fSkin[22]+0.5590169943749475*alphaDrSurf[7]*fSkin[21]+0.5590169943749475*alphaDrSurf[1]*fSkin[20]+0.5590169943749476*alphaDrSurf[6]*fSkin[19]+0.5590169943749476*alphaDrSurf[0]*fSkin[18]+0.5590169943749476*alphaDrSurf[4]*fSkin[17]+0.5590169943749475*alphaDrSurf[2]*fSkin[16]+0.4330127018922193*alphaDrSurf[5]*fSkin[15]+0.4330127018922193*alphaDrSurf[3]*fSkin[14]+0.4330127018922193*alphaDrSurf[7]*fSkin[13]+0.4330127018922193*alphaDrSurf[1]*fSkin[12]+0.25*alphaDrSurf[5]*fSkin[11]+0.4330127018922193*alphaDrSurf[6]*fSkin[10]+0.4330127018922193*alphaDrSurf[0]*fSkin[9]+0.4330127018922193*alphaDrSurf[4]*fSkin[8]+0.25*alphaDrSurf[3]*fSkin[7]+0.25*fSkin[6]*alphaDrSurf[7]+0.25*fSkin[3]*alphaDrSurf[6]+0.25*alphaDrSurf[1]*fSkin[5]+0.4330127018922193*alphaDrSurf[2]*fSkin[4]+0.25*fSkin[1]*alphaDrSurf[4]+0.25*alphaDrSurf[0]*fSkin[2]+0.25*fSkin[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.5590169943749476*alphaDrSurf[4]*fSkin[23]+0.5590169943749475*alphaDrSurf[2]*fSkin[22]+0.5590169943749475*alphaDrSurf[1]*fSkin[21]+0.5590169943749475*alphaDrSurf[7]*fSkin[20]+0.5590169943749476*alphaDrSurf[0]*fSkin[19]+0.5590169943749476*alphaDrSurf[6]*fSkin[18]+0.5590169943749476*alphaDrSurf[5]*fSkin[17]+0.5590169943749475*alphaDrSurf[3]*fSkin[16]+0.4330127018922193*alphaDrSurf[4]*fSkin[15]+0.4330127018922193*alphaDrSurf[2]*fSkin[14]+0.4330127018922193*alphaDrSurf[1]*fSkin[13]+0.4330127018922193*alphaDrSurf[7]*fSkin[12]+0.25*alphaDrSurf[4]*fSkin[11]+0.4330127018922193*alphaDrSurf[0]*fSkin[10]+0.4330127018922193*alphaDrSurf[6]*fSkin[9]+0.4330127018922193*alphaDrSurf[5]*fSkin[8]+0.25*alphaDrSurf[2]*fSkin[7]+0.25*fSkin[5]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fSkin[6]+0.25*fSkin[2]*alphaDrSurf[6]+0.25*fSkin[1]*alphaDrSurf[5]+0.4330127018922193*alphaDrSurf[3]*fSkin[4]+0.25*alphaDrSurf[0]*fSkin[3]+0.25*fSkin[0]*alphaDrSurf[3]; 
  Ghat[4] = 0.5590169943749476*alphaDrSurf[3]*fSkin[23]+0.5590169943749475*alphaDrSurf[5]*fSkin[22]+0.5590169943749475*alphaDrSurf[6]*fSkin[21]+0.5590169943749475*alphaDrSurf[0]*fSkin[20]+0.5590169943749476*alphaDrSurf[7]*fSkin[19]+0.5590169943749476*alphaDrSurf[1]*fSkin[18]+0.5590169943749476*alphaDrSurf[2]*fSkin[17]+0.5590169943749475*alphaDrSurf[4]*fSkin[16]+0.4330127018922193*alphaDrSurf[3]*fSkin[15]+0.4330127018922193*alphaDrSurf[5]*fSkin[14]+0.4330127018922193*alphaDrSurf[6]*fSkin[13]+0.4330127018922193*alphaDrSurf[0]*fSkin[12]+0.25*alphaDrSurf[3]*fSkin[11]+0.4330127018922193*alphaDrSurf[7]*fSkin[10]+0.4330127018922193*alphaDrSurf[1]*fSkin[9]+0.4330127018922193*alphaDrSurf[2]*fSkin[8]+0.25*alphaDrSurf[5]*fSkin[7]+0.25*fSkin[3]*alphaDrSurf[7]+0.25*alphaDrSurf[6]*fSkin[6]+0.25*alphaDrSurf[0]*fSkin[5]+0.4330127018922193*alphaDrSurf[4]*fSkin[4]+0.25*fSkin[0]*alphaDrSurf[4]+0.25*alphaDrSurf[1]*fSkin[2]+0.25*fSkin[1]*alphaDrSurf[2]; 
  Ghat[5] = 0.5590169943749476*alphaDrSurf[2]*fSkin[23]+0.5590169943749475*alphaDrSurf[4]*fSkin[22]+0.5590169943749475*alphaDrSurf[0]*fSkin[21]+0.5590169943749475*alphaDrSurf[6]*fSkin[20]+0.5590169943749476*alphaDrSurf[1]*fSkin[19]+0.5590169943749476*alphaDrSurf[7]*fSkin[18]+0.5590169943749476*alphaDrSurf[3]*fSkin[17]+0.5590169943749475*alphaDrSurf[5]*fSkin[16]+0.4330127018922193*alphaDrSurf[2]*fSkin[15]+0.4330127018922193*alphaDrSurf[4]*fSkin[14]+0.4330127018922193*alphaDrSurf[0]*fSkin[13]+0.4330127018922193*alphaDrSurf[6]*fSkin[12]+0.25*alphaDrSurf[2]*fSkin[11]+0.4330127018922193*alphaDrSurf[1]*fSkin[10]+0.4330127018922193*alphaDrSurf[7]*fSkin[9]+0.4330127018922193*alphaDrSurf[3]*fSkin[8]+0.25*alphaDrSurf[4]*fSkin[7]+0.25*fSkin[2]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fSkin[6]+0.25*fSkin[5]*alphaDrSurf[6]+0.4330127018922193*fSkin[4]*alphaDrSurf[5]+0.25*fSkin[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fSkin[3]+0.25*fSkin[1]*alphaDrSurf[3]; 
  Ghat[6] = 0.5590169943749476*alphaDrSurf[1]*fSkin[23]+0.5590169943749475*alphaDrSurf[0]*fSkin[22]+0.5590169943749475*alphaDrSurf[4]*fSkin[21]+0.5590169943749475*alphaDrSurf[5]*fSkin[20]+0.5590169943749476*alphaDrSurf[2]*fSkin[19]+0.5590169943749476*alphaDrSurf[3]*fSkin[18]+0.5590169943749476*alphaDrSurf[7]*fSkin[17]+0.5590169943749475*alphaDrSurf[6]*fSkin[16]+0.4330127018922193*alphaDrSurf[1]*fSkin[15]+0.4330127018922193*alphaDrSurf[0]*fSkin[14]+0.4330127018922193*alphaDrSurf[4]*fSkin[13]+0.4330127018922193*alphaDrSurf[5]*fSkin[12]+0.25*alphaDrSurf[1]*fSkin[11]+0.4330127018922193*alphaDrSurf[2]*fSkin[10]+0.4330127018922193*alphaDrSurf[3]*fSkin[9]+0.4330127018922193*alphaDrSurf[7]*fSkin[8]+0.25*alphaDrSurf[0]*fSkin[7]+0.25*fSkin[1]*alphaDrSurf[7]+0.25*alphaDrSurf[4]*fSkin[6]+0.4330127018922193*fSkin[4]*alphaDrSurf[6]+0.25*fSkin[0]*alphaDrSurf[6]+0.25*alphaDrSurf[5]*fSkin[5]+0.25*alphaDrSurf[2]*fSkin[3]+0.25*fSkin[2]*alphaDrSurf[3]; 
  Ghat[7] = 0.5590169943749476*alphaDrSurf[0]*fSkin[23]+0.5590169943749475*alphaDrSurf[1]*fSkin[22]+0.5590169943749475*alphaDrSurf[2]*fSkin[21]+0.5590169943749475*alphaDrSurf[3]*fSkin[20]+0.5590169943749476*alphaDrSurf[4]*fSkin[19]+0.5590169943749476*alphaDrSurf[5]*fSkin[18]+0.5590169943749476*alphaDrSurf[6]*fSkin[17]+0.5590169943749475*alphaDrSurf[7]*fSkin[16]+0.4330127018922193*alphaDrSurf[0]*fSkin[15]+0.4330127018922193*alphaDrSurf[1]*fSkin[14]+0.4330127018922193*alphaDrSurf[2]*fSkin[13]+0.4330127018922193*alphaDrSurf[3]*fSkin[12]+0.25*alphaDrSurf[0]*fSkin[11]+0.4330127018922193*alphaDrSurf[4]*fSkin[10]+0.4330127018922193*alphaDrSurf[5]*fSkin[9]+0.4330127018922193*alphaDrSurf[6]*fSkin[8]+0.25*alphaDrSurf[1]*fSkin[7]+0.4330127018922193*fSkin[4]*alphaDrSurf[7]+0.25*fSkin[0]*alphaDrSurf[7]+0.25*alphaDrSurf[2]*fSkin[6]+0.25*fSkin[1]*alphaDrSurf[6]+0.25*alphaDrSurf[3]*fSkin[5]+0.25*fSkin[2]*alphaDrSurf[5]+0.25*fSkin[3]*alphaDrSurf[4]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += 0.7071067811865475*Ghat[3]*dv1par; 
  out[4] += 1.224744871391589*Ghat[0]*dv1par; 
  out[5] += 0.7071067811865475*Ghat[4]*dv1par; 
  out[6] += 0.7071067811865475*Ghat[5]*dv1par; 
  out[7] += 0.7071067811865475*Ghat[6]*dv1par; 
  out[8] += 1.224744871391589*Ghat[1]*dv1par; 
  out[9] += 1.224744871391589*Ghat[2]*dv1par; 
  out[10] += 1.224744871391589*Ghat[3]*dv1par; 
  out[11] += 0.7071067811865475*Ghat[7]*dv1par; 
  out[12] += 1.224744871391589*Ghat[4]*dv1par; 
  out[13] += 1.224744871391589*Ghat[5]*dv1par; 
  out[14] += 1.224744871391589*Ghat[6]*dv1par; 
  out[15] += 1.224744871391589*Ghat[7]*dv1par; 
  out[16] += 1.58113883008419*Ghat[0]*dv1par; 
  out[17] += 1.58113883008419*Ghat[1]*dv1par; 
  out[18] += 1.58113883008419*Ghat[2]*dv1par; 
  out[19] += 1.58113883008419*Ghat[3]*dv1par; 
  out[20] += 1.58113883008419*Ghat[4]*dv1par; 
  out[21] += 1.58113883008419*Ghat[5]*dv1par; 
  out[22] += 1.58113883008419*Ghat[6]*dv1par; 
  out[23] += 1.58113883008419*Ghat[7]*dv1par; 

  } else { 

  alphaDrSurf[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 
  alphaDrSurf[3] = nu[3]*wvpar-0.5*nu[3]*dvpar; 
  alphaDrSurf[4] = nu[4]*wvpar-0.5*nu[4]*dvpar; 
  alphaDrSurf[5] = nu[5]*wvpar-0.5*nu[5]*dvpar; 
  alphaDrSurf[6] = nu[6]*wvpar-0.5*nu[6]*dvpar; 
  alphaDrSurf[7] = nu[7]*wvpar-0.5*nu[7]*dvpar; 

  Ghat[0] = 0.5590169943749476*alphaDrSurf[7]*fSkin[23]+0.5590169943749475*alphaDrSurf[6]*fSkin[22]+0.5590169943749475*alphaDrSurf[5]*fSkin[21]+0.5590169943749475*alphaDrSurf[4]*fSkin[20]+0.5590169943749476*alphaDrSurf[3]*fSkin[19]+0.5590169943749476*alphaDrSurf[2]*fSkin[18]+0.5590169943749476*alphaDrSurf[1]*fSkin[17]+0.5590169943749475*alphaDrSurf[0]*fSkin[16]-0.4330127018922193*alphaDrSurf[7]*fSkin[15]-0.4330127018922193*alphaDrSurf[6]*fSkin[14]-0.4330127018922193*alphaDrSurf[5]*fSkin[13]-0.4330127018922193*alphaDrSurf[4]*fSkin[12]+0.25*alphaDrSurf[7]*fSkin[11]-0.4330127018922193*alphaDrSurf[3]*fSkin[10]-0.4330127018922193*alphaDrSurf[2]*fSkin[9]-0.4330127018922193*alphaDrSurf[1]*fSkin[8]+0.25*alphaDrSurf[6]*fSkin[7]+0.25*alphaDrSurf[5]*fSkin[6]+0.25*alphaDrSurf[4]*fSkin[5]-0.4330127018922193*alphaDrSurf[0]*fSkin[4]+0.25*alphaDrSurf[3]*fSkin[3]+0.25*alphaDrSurf[2]*fSkin[2]+0.25*alphaDrSurf[1]*fSkin[1]+0.25*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 0.5590169943749476*alphaDrSurf[6]*fSkin[23]+0.5590169943749475*alphaDrSurf[7]*fSkin[22]+0.5590169943749475*alphaDrSurf[3]*fSkin[21]+0.5590169943749475*alphaDrSurf[2]*fSkin[20]+0.5590169943749476*alphaDrSurf[5]*fSkin[19]+0.5590169943749476*alphaDrSurf[4]*fSkin[18]+0.5590169943749476*alphaDrSurf[0]*fSkin[17]+0.5590169943749475*alphaDrSurf[1]*fSkin[16]-0.4330127018922193*alphaDrSurf[6]*fSkin[15]-0.4330127018922193*alphaDrSurf[7]*fSkin[14]-0.4330127018922193*alphaDrSurf[3]*fSkin[13]-0.4330127018922193*alphaDrSurf[2]*fSkin[12]+0.25*alphaDrSurf[6]*fSkin[11]-0.4330127018922193*alphaDrSurf[5]*fSkin[10]-0.4330127018922193*alphaDrSurf[4]*fSkin[9]-0.4330127018922193*alphaDrSurf[0]*fSkin[8]+0.25*alphaDrSurf[7]*fSkin[7]+0.25*alphaDrSurf[3]*fSkin[6]+0.25*alphaDrSurf[2]*fSkin[5]+0.25*fSkin[3]*alphaDrSurf[5]-0.4330127018922193*alphaDrSurf[1]*fSkin[4]+0.25*fSkin[2]*alphaDrSurf[4]+0.25*alphaDrSurf[0]*fSkin[1]+0.25*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.5590169943749476*alphaDrSurf[5]*fSkin[23]+0.5590169943749475*alphaDrSurf[3]*fSkin[22]+0.5590169943749475*alphaDrSurf[7]*fSkin[21]+0.5590169943749475*alphaDrSurf[1]*fSkin[20]+0.5590169943749476*alphaDrSurf[6]*fSkin[19]+0.5590169943749476*alphaDrSurf[0]*fSkin[18]+0.5590169943749476*alphaDrSurf[4]*fSkin[17]+0.5590169943749475*alphaDrSurf[2]*fSkin[16]-0.4330127018922193*alphaDrSurf[5]*fSkin[15]-0.4330127018922193*alphaDrSurf[3]*fSkin[14]-0.4330127018922193*alphaDrSurf[7]*fSkin[13]-0.4330127018922193*alphaDrSurf[1]*fSkin[12]+0.25*alphaDrSurf[5]*fSkin[11]-0.4330127018922193*alphaDrSurf[6]*fSkin[10]-0.4330127018922193*alphaDrSurf[0]*fSkin[9]-0.4330127018922193*alphaDrSurf[4]*fSkin[8]+0.25*alphaDrSurf[3]*fSkin[7]+0.25*fSkin[6]*alphaDrSurf[7]+0.25*fSkin[3]*alphaDrSurf[6]+0.25*alphaDrSurf[1]*fSkin[5]-0.4330127018922193*alphaDrSurf[2]*fSkin[4]+0.25*fSkin[1]*alphaDrSurf[4]+0.25*alphaDrSurf[0]*fSkin[2]+0.25*fSkin[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.5590169943749476*alphaDrSurf[4]*fSkin[23]+0.5590169943749475*alphaDrSurf[2]*fSkin[22]+0.5590169943749475*alphaDrSurf[1]*fSkin[21]+0.5590169943749475*alphaDrSurf[7]*fSkin[20]+0.5590169943749476*alphaDrSurf[0]*fSkin[19]+0.5590169943749476*alphaDrSurf[6]*fSkin[18]+0.5590169943749476*alphaDrSurf[5]*fSkin[17]+0.5590169943749475*alphaDrSurf[3]*fSkin[16]-0.4330127018922193*alphaDrSurf[4]*fSkin[15]-0.4330127018922193*alphaDrSurf[2]*fSkin[14]-0.4330127018922193*alphaDrSurf[1]*fSkin[13]-0.4330127018922193*alphaDrSurf[7]*fSkin[12]+0.25*alphaDrSurf[4]*fSkin[11]-0.4330127018922193*alphaDrSurf[0]*fSkin[10]-0.4330127018922193*alphaDrSurf[6]*fSkin[9]-0.4330127018922193*alphaDrSurf[5]*fSkin[8]+0.25*alphaDrSurf[2]*fSkin[7]+0.25*fSkin[5]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fSkin[6]+0.25*fSkin[2]*alphaDrSurf[6]+0.25*fSkin[1]*alphaDrSurf[5]-0.4330127018922193*alphaDrSurf[3]*fSkin[4]+0.25*alphaDrSurf[0]*fSkin[3]+0.25*fSkin[0]*alphaDrSurf[3]; 
  Ghat[4] = 0.5590169943749476*alphaDrSurf[3]*fSkin[23]+0.5590169943749475*alphaDrSurf[5]*fSkin[22]+0.5590169943749475*alphaDrSurf[6]*fSkin[21]+0.5590169943749475*alphaDrSurf[0]*fSkin[20]+0.5590169943749476*alphaDrSurf[7]*fSkin[19]+0.5590169943749476*alphaDrSurf[1]*fSkin[18]+0.5590169943749476*alphaDrSurf[2]*fSkin[17]+0.5590169943749475*alphaDrSurf[4]*fSkin[16]-0.4330127018922193*alphaDrSurf[3]*fSkin[15]-0.4330127018922193*alphaDrSurf[5]*fSkin[14]-0.4330127018922193*alphaDrSurf[6]*fSkin[13]-0.4330127018922193*alphaDrSurf[0]*fSkin[12]+0.25*alphaDrSurf[3]*fSkin[11]-0.4330127018922193*alphaDrSurf[7]*fSkin[10]-0.4330127018922193*alphaDrSurf[1]*fSkin[9]-0.4330127018922193*alphaDrSurf[2]*fSkin[8]+0.25*alphaDrSurf[5]*fSkin[7]+0.25*fSkin[3]*alphaDrSurf[7]+0.25*alphaDrSurf[6]*fSkin[6]+0.25*alphaDrSurf[0]*fSkin[5]-0.4330127018922193*alphaDrSurf[4]*fSkin[4]+0.25*fSkin[0]*alphaDrSurf[4]+0.25*alphaDrSurf[1]*fSkin[2]+0.25*fSkin[1]*alphaDrSurf[2]; 
  Ghat[5] = 0.5590169943749476*alphaDrSurf[2]*fSkin[23]+0.5590169943749475*alphaDrSurf[4]*fSkin[22]+0.5590169943749475*alphaDrSurf[0]*fSkin[21]+0.5590169943749475*alphaDrSurf[6]*fSkin[20]+0.5590169943749476*alphaDrSurf[1]*fSkin[19]+0.5590169943749476*alphaDrSurf[7]*fSkin[18]+0.5590169943749476*alphaDrSurf[3]*fSkin[17]+0.5590169943749475*alphaDrSurf[5]*fSkin[16]-0.4330127018922193*alphaDrSurf[2]*fSkin[15]-0.4330127018922193*alphaDrSurf[4]*fSkin[14]-0.4330127018922193*alphaDrSurf[0]*fSkin[13]-0.4330127018922193*alphaDrSurf[6]*fSkin[12]+0.25*alphaDrSurf[2]*fSkin[11]-0.4330127018922193*alphaDrSurf[1]*fSkin[10]-0.4330127018922193*alphaDrSurf[7]*fSkin[9]-0.4330127018922193*alphaDrSurf[3]*fSkin[8]+0.25*alphaDrSurf[4]*fSkin[7]+0.25*fSkin[2]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fSkin[6]+0.25*fSkin[5]*alphaDrSurf[6]-0.4330127018922193*fSkin[4]*alphaDrSurf[5]+0.25*fSkin[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fSkin[3]+0.25*fSkin[1]*alphaDrSurf[3]; 
  Ghat[6] = 0.5590169943749476*alphaDrSurf[1]*fSkin[23]+0.5590169943749475*alphaDrSurf[0]*fSkin[22]+0.5590169943749475*alphaDrSurf[4]*fSkin[21]+0.5590169943749475*alphaDrSurf[5]*fSkin[20]+0.5590169943749476*alphaDrSurf[2]*fSkin[19]+0.5590169943749476*alphaDrSurf[3]*fSkin[18]+0.5590169943749476*alphaDrSurf[7]*fSkin[17]+0.5590169943749475*alphaDrSurf[6]*fSkin[16]-0.4330127018922193*alphaDrSurf[1]*fSkin[15]-0.4330127018922193*alphaDrSurf[0]*fSkin[14]-0.4330127018922193*alphaDrSurf[4]*fSkin[13]-0.4330127018922193*alphaDrSurf[5]*fSkin[12]+0.25*alphaDrSurf[1]*fSkin[11]-0.4330127018922193*alphaDrSurf[2]*fSkin[10]-0.4330127018922193*alphaDrSurf[3]*fSkin[9]-0.4330127018922193*alphaDrSurf[7]*fSkin[8]+0.25*alphaDrSurf[0]*fSkin[7]+0.25*fSkin[1]*alphaDrSurf[7]+0.25*alphaDrSurf[4]*fSkin[6]-0.4330127018922193*fSkin[4]*alphaDrSurf[6]+0.25*fSkin[0]*alphaDrSurf[6]+0.25*alphaDrSurf[5]*fSkin[5]+0.25*alphaDrSurf[2]*fSkin[3]+0.25*fSkin[2]*alphaDrSurf[3]; 
  Ghat[7] = 0.5590169943749476*alphaDrSurf[0]*fSkin[23]+0.5590169943749475*alphaDrSurf[1]*fSkin[22]+0.5590169943749475*alphaDrSurf[2]*fSkin[21]+0.5590169943749475*alphaDrSurf[3]*fSkin[20]+0.5590169943749476*alphaDrSurf[4]*fSkin[19]+0.5590169943749476*alphaDrSurf[5]*fSkin[18]+0.5590169943749476*alphaDrSurf[6]*fSkin[17]+0.5590169943749475*alphaDrSurf[7]*fSkin[16]-0.4330127018922193*alphaDrSurf[0]*fSkin[15]-0.4330127018922193*alphaDrSurf[1]*fSkin[14]-0.4330127018922193*alphaDrSurf[2]*fSkin[13]-0.4330127018922193*alphaDrSurf[3]*fSkin[12]+0.25*alphaDrSurf[0]*fSkin[11]-0.4330127018922193*alphaDrSurf[4]*fSkin[10]-0.4330127018922193*alphaDrSurf[5]*fSkin[9]-0.4330127018922193*alphaDrSurf[6]*fSkin[8]+0.25*alphaDrSurf[1]*fSkin[7]-0.4330127018922193*fSkin[4]*alphaDrSurf[7]+0.25*fSkin[0]*alphaDrSurf[7]+0.25*alphaDrSurf[2]*fSkin[6]+0.25*fSkin[1]*alphaDrSurf[6]+0.25*alphaDrSurf[3]*fSkin[5]+0.25*fSkin[2]*alphaDrSurf[5]+0.25*fSkin[3]*alphaDrSurf[4]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += -0.7071067811865475*Ghat[3]*dv1par; 
  out[4] += 1.224744871391589*Ghat[0]*dv1par; 
  out[5] += -0.7071067811865475*Ghat[4]*dv1par; 
  out[6] += -0.7071067811865475*Ghat[5]*dv1par; 
  out[7] += -0.7071067811865475*Ghat[6]*dv1par; 
  out[8] += 1.224744871391589*Ghat[1]*dv1par; 
  out[9] += 1.224744871391589*Ghat[2]*dv1par; 
  out[10] += 1.224744871391589*Ghat[3]*dv1par; 
  out[11] += -0.7071067811865475*Ghat[7]*dv1par; 
  out[12] += 1.224744871391589*Ghat[4]*dv1par; 
  out[13] += 1.224744871391589*Ghat[5]*dv1par; 
  out[14] += 1.224744871391589*Ghat[6]*dv1par; 
  out[15] += 1.224744871391589*Ghat[7]*dv1par; 
  out[16] += -1.58113883008419*Ghat[0]*dv1par; 
  out[17] += -1.58113883008419*Ghat[1]*dv1par; 
  out[18] += -1.58113883008419*Ghat[2]*dv1par; 
  out[19] += -1.58113883008419*Ghat[3]*dv1par; 
  out[20] += -1.58113883008419*Ghat[4]*dv1par; 
  out[21] += -1.58113883008419*Ghat[5]*dv1par; 
  out[22] += -1.58113883008419*Ghat[6]*dv1par; 
  out[23] += -1.58113883008419*Ghat[7]*dv1par; 

  } 
} 
