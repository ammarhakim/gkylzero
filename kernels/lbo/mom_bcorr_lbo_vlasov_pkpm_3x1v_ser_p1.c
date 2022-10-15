#include <gkyl_mom_bcorr_lbo_vlasov_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_vlasov_pkpm_3x1v_ser_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double mass, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direction. 
  // mass:      species mass. 
  // fIn[24]:    distribution function at lower/upper velocity boundaries. 
  // out:       int dS of vf^(vmax)_(vmin). 
 
  if (edge == GKYL_VX_LOWER) {
    out[0] += vBoundary[0]*((-1.58113883008419*fIn[16])+1.224744871391589*fIn[4]-0.7071067811865475*fIn[0])*mass; 
    out[1] += vBoundary[0]*((-1.58113883008419*fIn[17])+1.224744871391589*fIn[8]-0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[0]*((-1.58113883008419*fIn[18])+1.224744871391589*fIn[9]-0.7071067811865475*fIn[2])*mass; 
    out[3] += vBoundary[0]*((-1.58113883008419*fIn[19])+1.224744871391589*fIn[10]-0.7071067811865475*fIn[3])*mass; 
    out[4] += vBoundary[0]*((-1.58113883008419*fIn[20])+1.224744871391589*fIn[12]-0.7071067811865475*fIn[5])*mass; 
    out[5] += vBoundary[0]*((-1.58113883008419*fIn[21])+1.224744871391589*fIn[13]-0.7071067811865475*fIn[6])*mass; 
    out[6] += vBoundary[0]*((-1.58113883008419*fIn[22])+1.224744871391589*fIn[14]-0.7071067811865475*fIn[7])*mass; 
    out[7] += vBoundary[0]*((-1.58113883008419*fIn[23])+1.224744871391589*fIn[15]-0.7071067811865475*fIn[11])*mass; 
  } 
  else if (edge == GKYL_VX_UPPER) {
    out[0] += vBoundary[1]*(1.58113883008419*fIn[16]+1.224744871391589*fIn[4]+0.7071067811865475*fIn[0])*mass; 
    out[1] += vBoundary[1]*(1.58113883008419*fIn[17]+1.224744871391589*fIn[8]+0.7071067811865475*fIn[1])*mass; 
    out[2] += vBoundary[1]*(1.58113883008419*fIn[18]+1.224744871391589*fIn[9]+0.7071067811865475*fIn[2])*mass; 
    out[3] += vBoundary[1]*(1.58113883008419*fIn[19]+1.224744871391589*fIn[10]+0.7071067811865475*fIn[3])*mass; 
    out[4] += vBoundary[1]*(1.58113883008419*fIn[20]+1.224744871391589*fIn[12]+0.7071067811865475*fIn[5])*mass; 
    out[5] += vBoundary[1]*(1.58113883008419*fIn[21]+1.224744871391589*fIn[13]+0.7071067811865475*fIn[6])*mass; 
    out[6] += vBoundary[1]*(1.58113883008419*fIn[22]+1.224744871391589*fIn[14]+0.7071067811865475*fIn[7])*mass; 
    out[7] += vBoundary[1]*(1.58113883008419*fIn[23]+1.224744871391589*fIn[15]+0.7071067811865475*fIn[11])*mass; 
  }
 
} 
 
