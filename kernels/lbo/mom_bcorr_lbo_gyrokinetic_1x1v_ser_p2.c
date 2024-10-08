#include <gkyl_mom_bcorr_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_1x1v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double _m, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direction. 
  // _m:        species mass. 
  // fIn[8]:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  dS = 1.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
    out[3] += vBoundary[0]*((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[4] += vBoundary[0]*((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[5] += vBoundary[0]*(1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
    out[3] += vBoundary[1]*(1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[4] += vBoundary[1]*(1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[5] += vBoundary[1]*(1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
  }
 
} 
 
