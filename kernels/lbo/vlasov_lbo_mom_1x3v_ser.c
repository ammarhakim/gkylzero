#include <gkyl_vlasov_lbo_mom_kernels.h> 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_f_p1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[2] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[3] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
 
  } else {
 
    out[2] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[3] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
 
  } else {
 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_vf_p1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[1]*fIn[2]-0.7071067811865475*fIn[0]*dxv[1])*dS; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[1]*fIn[5]-0.7071067811865475*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[1]*fIn[2])-0.7071067811865475*fIn[0]*dxv[1])*dS; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[1]*fIn[5])-0.7071067811865475*dxv[1]*fIn[1])*dS; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[6]-0.7071067811865475*fIn[1]*dxv[2])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[6])-0.7071067811865475*fIn[1]*dxv[2])*dS; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[8]-0.7071067811865475*fIn[1]*dxv[3])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[8])-0.7071067811865475*fIn[1]*dxv[3])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_f_p2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[3] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[4] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
 
  } else {
 
    out[3] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[4] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[6] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[7] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
 
  } else {
 
    out[6] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[7] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x3v_ser_vf_p2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } 
 
  dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } 
 
}
 
