#include <gkyl_vlasov_lbo_mom_kernels.h> 
GKYL_CU_DH void vlasov_boundary_integral_2x3v_ser_f_p1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]-1.414213562373095*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
  dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS; 
    out[7] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[6])*dS; 
 
  } else {
 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS; 
    out[7] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
  dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[8] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[20]-1.414213562373095*fIn[6])*dS; 
 
  } else {
 
    out[8] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[20]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_2x3v_ser_vf_p1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS; 
 
  dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[7]-0.7071067811865475*fIn[1]*dxv[2])*dS; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[8]-0.7071067811865475*dxv[2]*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]-1.414213562373095*fIn[6])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[16]-0.7071067811865475*dxv[2]*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[7])-0.7071067811865475*fIn[1]*dxv[2])*dS; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[8])-0.7071067811865475*dxv[2]*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]+1.414213562373095*fIn[6])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[16])-0.7071067811865475*dxv[2]*fIn[6])*dS; 
 
  } 
 
  dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[9]-0.7071067811865475*fIn[1]*dxv[3])*dS; 
    out[2] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[10]-0.7071067811865475*fIn[2]*dxv[3])*dS; 
    out[3] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[6])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[17]-0.7071067811865475*dxv[3]*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[9])-0.7071067811865475*fIn[1]*dxv[3])*dS; 
    out[2] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[10])-0.7071067811865475*fIn[2]*dxv[3])*dS; 
    out[3] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[6])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[17])-0.7071067811865475*dxv[3]*fIn[6])*dS; 
 
  } 
 
  dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[5]-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[12]-0.7071067811865475*fIn[1]*dxv[4])*dS; 
    out[2] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[13]-0.7071067811865475*fIn[2]*dxv[4])*dS; 
    out[3] += (2.449489742783178*fIn[20]-1.414213562373095*fIn[6])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[20]-0.7071067811865475*dxv[4]*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[5])-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[12])-0.7071067811865475*fIn[1]*dxv[4])*dS; 
    out[2] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[13])-0.7071067811865475*fIn[2]*dxv[4])*dS; 
    out[3] += (2.449489742783178*fIn[20]+1.414213562373095*fIn[6])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[20])-0.7071067811865475*dxv[4]*fIn[6])*dS; 
 
  } 
 
}
 
