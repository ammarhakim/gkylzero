#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_pkpm_2x1v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direction. 
  // fIn[20]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += vBoundary[0]*((-1.58113883008419*fIn[9])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[0]); 
    out[1] += vBoundary[0]*((-1.58113883008419*fIn[15])+1.224744871391589*fIn[5]-0.7071067811865475*fIn[1]); 
    out[2] += vBoundary[0]*((-1.58113883008419*fIn[16])+1.224744871391589*fIn[6]-0.7071067811865475*fIn[2]); 
    out[3] += vBoundary[0]*((-1.58113883008419*fIn[19])+1.224744871391589*fIn[10]-0.7071067811865475*fIn[4]); 
    out[4] += vBoundary[0]*(1.224744871391589*fIn[13]-0.7071067811865475*fIn[7]); 
    out[5] += vBoundary[0]*(1.224744871391589*fIn[14]-0.7071067811865475*fIn[8]); 
    out[6] += vBoundary[0]*(1.224744871391589*fIn[17]-0.7071067811865475*fIn[11]); 
    out[7] += vBoundary[0]*(1.224744871391589*fIn[18]-0.7071067811865475*fIn[12]); 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += vBoundary[1]*(1.58113883008419*fIn[9]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[0]); 
    out[1] += vBoundary[1]*(1.58113883008419*fIn[15]+1.224744871391589*fIn[5]+0.7071067811865475*fIn[1]); 
    out[2] += vBoundary[1]*(1.58113883008419*fIn[16]+1.224744871391589*fIn[6]+0.7071067811865475*fIn[2]); 
    out[3] += vBoundary[1]*(1.58113883008419*fIn[19]+1.224744871391589*fIn[10]+0.7071067811865475*fIn[4]); 
    out[4] += vBoundary[1]*(1.224744871391589*fIn[13]+0.7071067811865475*fIn[7]); 
    out[5] += vBoundary[1]*(1.224744871391589*fIn[14]+0.7071067811865475*fIn[8]); 
    out[6] += vBoundary[1]*(1.224744871391589*fIn[17]+0.7071067811865475*fIn[11]); 
    out[7] += vBoundary[1]*(1.224744871391589*fIn[18]+0.7071067811865475*fIn[12]); 
  }
 
} 
 
