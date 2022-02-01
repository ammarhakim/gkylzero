#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x3v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[3]*dxv[4]; 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]-1.414213562373095*fIn[6])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[3]*dxv[4]; 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[4]; 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS; 
    out[7] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[6])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[4]; 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS; 
    out[7] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[8] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[20]-1.414213562373095*fIn[6])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
    out[8] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[20]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_2x3v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.25*dxv[3]*dxv[4]; 
      out[0] += (1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS+vBoundary[0]*(2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
      out[1] += (1.224744871391589*dxv[2]*fIn[7]-0.7071067811865475*fIn[1]*dxv[2])*dS+vBoundary[0]*(2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS; 
      out[2] += (1.224744871391589*dxv[2]*fIn[8]-0.7071067811865475*dxv[2]*fIn[2])*dS+vBoundary[0]*(2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS; 
      out[3] += (1.224744871391589*dxv[2]*fIn[16]-0.7071067811865475*dxv[2]*fIn[6])*dS+vBoundary[0]*(2.449489742783178*fIn[16]-1.414213562373095*fIn[6])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.25*dxv[3]*dxv[4]; 
      out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*vBoundary[3]*dS+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
      out[1] += ((-1.224744871391589*dxv[2]*fIn[7])-0.7071067811865475*fIn[1]*dxv[2])*dS+vBoundary[3]*(2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS; 
      out[2] += ((-1.224744871391589*dxv[2]*fIn[8])-0.7071067811865475*dxv[2]*fIn[2])*dS+vBoundary[3]*(2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS; 
      out[3] += ((-1.224744871391589*dxv[2]*fIn[16])-0.7071067811865475*dxv[2]*fIn[6])*dS+vBoundary[3]*(2.449489742783178*fIn[16]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[4]; 
      out[0] += (1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS+vBoundary[1]*(2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
      out[1] += (1.224744871391589*dxv[3]*fIn[9]-0.7071067811865475*fIn[1]*dxv[3])*dS+vBoundary[1]*(2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS; 
      out[2] += (1.224744871391589*dxv[3]*fIn[10]-0.7071067811865475*fIn[2]*dxv[3])*dS+vBoundary[1]*(2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS; 
      out[3] += (1.224744871391589*dxv[3]*fIn[17]-0.7071067811865475*dxv[3]*fIn[6])*dS+vBoundary[1]*(2.449489742783178*fIn[17]-1.414213562373095*fIn[6])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[4]; 
      out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*vBoundary[4]*dS+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
      out[1] += ((-1.224744871391589*dxv[3]*fIn[9])-0.7071067811865475*fIn[1]*dxv[3])*dS+vBoundary[4]*(2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS; 
      out[2] += ((-1.224744871391589*dxv[3]*fIn[10])-0.7071067811865475*fIn[2]*dxv[3])*dS+vBoundary[4]*(2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS; 
      out[3] += ((-1.224744871391589*dxv[3]*fIn[17])-0.7071067811865475*dxv[3]*fIn[6])*dS+vBoundary[4]*(2.449489742783178*fIn[17]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
  if (edge == GKYL_VZ_LOWER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
      out[0] += (1.224744871391589*dxv[4]*fIn[5]-0.7071067811865475*fIn[0]*dxv[4])*dS+vBoundary[2]*(2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
      out[1] += (1.224744871391589*dxv[4]*fIn[12]-0.7071067811865475*fIn[1]*dxv[4])*dS+vBoundary[2]*(2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS; 
      out[2] += (1.224744871391589*dxv[4]*fIn[13]-0.7071067811865475*fIn[2]*dxv[4])*dS+vBoundary[2]*(2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS; 
      out[3] += (1.224744871391589*dxv[4]*fIn[20]-0.7071067811865475*dxv[4]*fIn[6])*dS+vBoundary[2]*(2.449489742783178*fIn[20]-1.414213562373095*fIn[6])*dS; 
 
  } else if (edge == GKYL_VZ_UPPER) {
 
    dS = 0.25*dxv[2]*dxv[3]; 
      out[0] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*vBoundary[5]*dS+((-1.224744871391589*dxv[4]*fIn[5])-0.7071067811865475*fIn[0]*dxv[4])*dS; 
      out[1] += ((-1.224744871391589*dxv[4]*fIn[12])-0.7071067811865475*fIn[1]*dxv[4])*dS+vBoundary[5]*(2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS; 
      out[2] += ((-1.224744871391589*dxv[4]*fIn[13])-0.7071067811865475*fIn[2]*dxv[4])*dS+vBoundary[5]*(2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS; 
      out[3] += ((-1.224744871391589*dxv[4]*fIn[20])-0.7071067811865475*dxv[4]*fIn[6])*dS+vBoundary[5]*(2.449489742783178*fIn[20]+1.414213562373095*fIn[6])*dS; 
 
  } 
 
}
 
