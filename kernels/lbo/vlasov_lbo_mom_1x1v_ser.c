#include <gkyl_vlasov_lbo_mom_kernels.h> 
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_f_p1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 1.0; 
 
  if (atLower) {
 
    out[0] += (1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_vf_p1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 1.0; 
 
  if (atLower) {
 
    out[0] += (1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary+(0.6123724356957944*dxv[1]*fIn[2]-0.3535533905932737*fIn[0]*dxv[1])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary+(0.6123724356957944*dxv[1]*fIn[3]-0.3535533905932737*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary+((-0.6123724356957944*dxv[1]*fIn[2])-0.3535533905932737*fIn[0]*dxv[1])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary+((-0.6123724356957944*dxv[1]*fIn[3])-0.3535533905932737*dxv[1]*fIn[1])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_f_p2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 1.0; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x1v_ser_vf_p2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 1.0; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS*vBoundary; 
 
  } 
 
}
 
