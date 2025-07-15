#include <gkyl_mom_bcorr_lbo_pkpm_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_pkpm_2x1v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direction. 
  // mass:      species mass. 
  // fIn[27]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f^(vmax)_(vmin) and vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += ((-1.58113883008419*fIn[9])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[0])*mass; 
    out[1] += ((-1.58113883008419*fIn[15])+1.224744871391589*fIn[5]-0.7071067811865475*fIn[1])*mass; 
    out[2] += ((-1.58113883008419*fIn[16])+1.224744871391589*fIn[6]-0.7071067811865475*fIn[2])*mass; 
    out[3] += ((-1.58113883008419*fIn[19])+1.224744871391589*fIn[10]-0.7071067811865475*fIn[4])*mass; 
    out[4] += ((-1.58113883008419*fIn[21])+1.224744871391589*fIn[13]-0.7071067811865475*fIn[7])*mass; 
    out[5] += ((-1.58113883008419*fIn[22])+1.224744871391589*fIn[14]-0.7071067811865475*fIn[8])*mass; 
    out[6] += ((-1.58113883008419*fIn[24])+1.224744871391589*fIn[17]-0.7071067811865475*fIn[11])*mass; 
    out[7] += ((-1.58113883008419*fIn[25])+1.224744871391589*fIn[18]-0.7071067811865475*fIn[12])*mass; 
    out[8] += ((-1.58113883008419*fIn[26])+1.224744871391589*fIn[23]-0.7071067811865475*fIn[20])*mass; 
    out[9] += out[0]*vBoundary[0]; 
    out[10] += vBoundary[0]*out[1]; 
    out[11] += vBoundary[0]*out[2]; 
    out[12] += vBoundary[0]*out[3]; 
    out[13] += vBoundary[0]*out[4]; 
    out[14] += vBoundary[0]*out[5]; 
    out[15] += vBoundary[0]*out[6]; 
    out[16] += vBoundary[0]*out[7]; 
    out[17] += vBoundary[0]*out[8]; 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += (1.58113883008419*fIn[9]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[0])*mass; 
    out[1] += (1.58113883008419*fIn[15]+1.224744871391589*fIn[5]+0.7071067811865475*fIn[1])*mass; 
    out[2] += (1.58113883008419*fIn[16]+1.224744871391589*fIn[6]+0.7071067811865475*fIn[2])*mass; 
    out[3] += (1.58113883008419*fIn[19]+1.224744871391589*fIn[10]+0.7071067811865475*fIn[4])*mass; 
    out[4] += (1.58113883008419*fIn[21]+1.224744871391589*fIn[13]+0.7071067811865475*fIn[7])*mass; 
    out[5] += (1.58113883008419*fIn[22]+1.224744871391589*fIn[14]+0.7071067811865475*fIn[8])*mass; 
    out[6] += (1.58113883008419*fIn[24]+1.224744871391589*fIn[17]+0.7071067811865475*fIn[11])*mass; 
    out[7] += (1.58113883008419*fIn[25]+1.224744871391589*fIn[18]+0.7071067811865475*fIn[12])*mass; 
    out[8] += (1.58113883008419*fIn[26]+1.224744871391589*fIn[23]+0.7071067811865475*fIn[20])*mass; 
    out[9] += out[0]*vBoundary[1]; 
    out[10] += out[1]*vBoundary[1]; 
    out[11] += vBoundary[1]*out[2]; 
    out[12] += vBoundary[1]*out[3]; 
    out[13] += vBoundary[1]*out[4]; 
    out[14] += vBoundary[1]*out[5]; 
    out[15] += vBoundary[1]*out[6]; 
    out[16] += vBoundary[1]*out[7]; 
    out[17] += vBoundary[1]*out[8]; 
  }
 
} 
 
