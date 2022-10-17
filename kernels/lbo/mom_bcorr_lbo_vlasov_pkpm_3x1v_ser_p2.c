#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_pkpm_3x1v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direction. 
  // mass:      species mass. 
  // fIn[48]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += vBoundary[0]*((-1.58113883008419*fIn[14])+1.224744871391589*fIn[4]-0.7071067811865475*fIn[0])*mass; 
    out[1] += vBoundary[0]*((-1.58113883008419*fIn[28])+1.224744871391589*fIn[8]-0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[0]*((-1.58113883008419*fIn[29])+1.224744871391589*fIn[9]-0.7071067811865475*fIn[2])*mass; 
    out[3] += vBoundary[0]*((-1.58113883008419*fIn[30])+1.224744871391589*fIn[10]-0.7071067811865475*fIn[3])*mass; 
    out[4] += vBoundary[0]*((-1.58113883008419*fIn[41])+1.224744871391589*fIn[16]-0.7071067811865475*fIn[5])*mass; 
    out[5] += vBoundary[0]*((-1.58113883008419*fIn[42])+1.224744871391589*fIn[17]-0.7071067811865475*fIn[6])*mass; 
    out[6] += vBoundary[0]*((-1.58113883008419*fIn[43])+1.224744871391589*fIn[18]-0.7071067811865475*fIn[7])*mass; 
    out[7] += vBoundary[0]*(1.224744871391589*fIn[25]-0.7071067811865475*fIn[11])*mass; 
    out[8] += vBoundary[0]*(1.224744871391589*fIn[26]-0.7071067811865475*fIn[12])*mass; 
    out[9] += vBoundary[0]*(1.224744871391589*fIn[27]-0.7071067811865475*fIn[13])*mass; 
    out[10] += vBoundary[0]*((-1.58113883008419*fIn[47])+1.224744871391589*fIn[31]-0.7071067811865475*fIn[15])*mass; 
    out[11] += vBoundary[0]*(1.224744871391589*fIn[35]-0.7071067811865475*fIn[19])*mass; 
    out[12] += vBoundary[0]*(1.224744871391589*fIn[36]-0.7071067811865475*fIn[20])*mass; 
    out[13] += vBoundary[0]*(1.224744871391589*fIn[37]-0.7071067811865475*fIn[21])*mass; 
    out[14] += vBoundary[0]*(1.224744871391589*fIn[38]-0.7071067811865475*fIn[22])*mass; 
    out[15] += vBoundary[0]*(1.224744871391589*fIn[39]-0.7071067811865475*fIn[23])*mass; 
    out[16] += vBoundary[0]*(1.224744871391589*fIn[40]-0.7071067811865475*fIn[24])*mass; 
    out[17] += vBoundary[0]*(1.224744871391589*fIn[44]-0.7071067811865475*fIn[32])*mass; 
    out[18] += vBoundary[0]*(1.224744871391589*fIn[45]-0.7071067811865475*fIn[33])*mass; 
    out[19] += vBoundary[0]*(1.224744871391589*fIn[46]-0.7071067811865475*fIn[34])*mass; 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += vBoundary[1]*(1.58113883008419*fIn[14]+1.224744871391589*fIn[4]+0.7071067811865475*fIn[0])*mass; 
    out[1] += vBoundary[1]*(1.58113883008419*fIn[28]+1.224744871391589*fIn[8]+0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[1]*(1.58113883008419*fIn[29]+1.224744871391589*fIn[9]+0.7071067811865475*fIn[2])*mass; 
    out[3] += vBoundary[1]*(1.58113883008419*fIn[30]+1.224744871391589*fIn[10]+0.7071067811865475*fIn[3])*mass; 
    out[4] += vBoundary[1]*(1.58113883008419*fIn[41]+1.224744871391589*fIn[16]+0.7071067811865475*fIn[5])*mass; 
    out[5] += vBoundary[1]*(1.58113883008419*fIn[42]+1.224744871391589*fIn[17]+0.7071067811865475*fIn[6])*mass; 
    out[6] += vBoundary[1]*(1.58113883008419*fIn[43]+1.224744871391589*fIn[18]+0.7071067811865475*fIn[7])*mass; 
    out[7] += vBoundary[1]*(1.224744871391589*fIn[25]+0.7071067811865475*fIn[11])*mass; 
    out[8] += vBoundary[1]*(1.224744871391589*fIn[26]+0.7071067811865475*fIn[12])*mass; 
    out[9] += vBoundary[1]*(1.224744871391589*fIn[27]+0.7071067811865475*fIn[13])*mass; 
    out[10] += vBoundary[1]*(1.58113883008419*fIn[47]+1.224744871391589*fIn[31]+0.7071067811865475*fIn[15])*mass; 
    out[11] += vBoundary[1]*(1.224744871391589*fIn[35]+0.7071067811865475*fIn[19])*mass; 
    out[12] += vBoundary[1]*(1.224744871391589*fIn[36]+0.7071067811865475*fIn[20])*mass; 
    out[13] += vBoundary[1]*(1.224744871391589*fIn[37]+0.7071067811865475*fIn[21])*mass; 
    out[14] += vBoundary[1]*(1.224744871391589*fIn[38]+0.7071067811865475*fIn[22])*mass; 
    out[15] += vBoundary[1]*(1.224744871391589*fIn[39]+0.7071067811865475*fIn[23])*mass; 
    out[16] += vBoundary[1]*(1.224744871391589*fIn[40]+0.7071067811865475*fIn[24])*mass; 
    out[17] += vBoundary[1]*(1.224744871391589*fIn[44]+0.7071067811865475*fIn[32])*mass; 
    out[18] += vBoundary[1]*(1.224744871391589*fIn[45]+0.7071067811865475*fIn[33])*mass; 
    out[19] += vBoundary[1]*(1.224744871391589*fIn[46]+0.7071067811865475*fIn[34])*mass; 
  }
 
} 
 
