#include <gkyl_rad_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:     cell-center coordinates. 
  // dxv[2]:   cell spacing. 
  // nI:        ion density. 
  // nuField:       2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 
  double alphaDrSurf[2] = {0.0}; 
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  if (w[1]>0) { 

  Ghat[0] = 1.118033988749895*(alphaDrSurf[1]*fEdge[5]+alphaDrSurf[0]*fEdge[4])-0.8660254037844386*(alphaDrSurf[1]*fEdge[3]+alphaDrSurf[0]*fEdge[2])+0.5*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0]); 
  Ghat[1] = 1.118033988749895*(alphaDrSurf[0]*fEdge[5]+alphaDrSurf[1]*fEdge[4])-0.8660254037844386*(alphaDrSurf[0]*fEdge[3]+alphaDrSurf[1]*fEdge[2])+0.5*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1]); 

  } else { 

  Ghat[0] = 0.1666666666666667*(6.708203932499369*alphaDrSurf[1]*fSkin[5]+6.708203932499369*alphaDrSurf[0]*fSkin[4]+5.196152422706631*(alphaDrSurf[1]*fSkin[3]+alphaDrSurf[0]*fSkin[2])+3.0*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = 0.1666666666666667*(6.708203932499369*alphaDrSurf[0]*fSkin[5]+6.708203932499369*alphaDrSurf[1]*fSkin[4]+5.196152422706631*(alphaDrSurf[0]*fSkin[3]+alphaDrSurf[1]*fSkin[2])+3.0*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -1.58113883008419*Ghat[0]*rdv2; 
  out[5] += -1.58113883008419*Ghat[1]*rdv2; 

  } else { 

  if (w[1]>0) { 

  Ghat[0] = 1.118033988749895*(alphaDrSurf[1]*fSkin[5]+alphaDrSurf[0]*fSkin[4])-0.8660254037844386*(alphaDrSurf[1]*fSkin[3]+alphaDrSurf[0]*fSkin[2])+0.5*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0]); 
  Ghat[1] = 1.118033988749895*(alphaDrSurf[0]*fSkin[5]+alphaDrSurf[1]*fSkin[4])-0.8660254037844386*(alphaDrSurf[0]*fSkin[3]+alphaDrSurf[1]*fSkin[2])+0.5*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1]); 

  } else { 

  Ghat[0] = 0.1666666666666667*(6.708203932499369*alphaDrSurf[1]*fEdge[5]+6.708203932499369*alphaDrSurf[0]*fEdge[4]+5.196152422706631*(alphaDrSurf[1]*fEdge[3]+alphaDrSurf[0]*fEdge[2])+3.0*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = 0.1666666666666667*(6.708203932499369*alphaDrSurf[0]*fEdge[5]+6.708203932499369*alphaDrSurf[1]*fEdge[4]+5.196152422706631*(alphaDrSurf[0]*fEdge[3]+alphaDrSurf[1]*fEdge[2])+3.0*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 1.58113883008419*Ghat[0]*rdv2; 
  out[5] += 1.58113883008419*Ghat[1]*rdv2; 

  }
  //  printf("w[1]=%e, rdv= %f, out: 0=%e, 1=%e, 2=%e, 3=%e, 4=%e, 5=%e\n",w[1], rdv2,out[0],out[1],out[2],out[3],out[4],out[5]);
  return 0.;

} 
