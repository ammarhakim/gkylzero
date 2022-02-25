#include <gkyl_mom_bcorr_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_2x2v_tensor_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double _m, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direction. 
  // _m:        species mass. 
  // fIn[81]:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 1.0*dxv[3]*M_PI/_m; 
    out[0] += ((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += ((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += ((-2.23606797749979*fIn[34])+1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
    out[4] += ((-2.23606797749979*fIn[45])+1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
    out[5] += ((-2.23606797749979*fIn[46])+1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
    out[6] += ((-2.23606797749979*fIn[55])+1.732050807568877*fIn[32]-1.0*fIn[19])*dS; 
    out[7] += ((-2.23606797749979*fIn[56])+1.732050807568877*fIn[33]-1.0*fIn[20])*dS; 
    out[8] += ((-2.23606797749979*fIn[72])+1.732050807568877*fIn[54]-1.0*fIn[44])*dS; 
    out[9] += vBoundary[0]*((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[10] += vBoundary[0]*((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[11] += vBoundary[0]*((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[12] += vBoundary[0]*((-2.23606797749979*fIn[34])+1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
    out[13] += vBoundary[0]*((-2.23606797749979*fIn[45])+1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
    out[14] += vBoundary[0]*((-2.23606797749979*fIn[46])+1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
    out[15] += vBoundary[0]*((-2.23606797749979*fIn[55])+1.732050807568877*fIn[32]-1.0*fIn[19])*dS; 
    out[16] += vBoundary[0]*((-2.23606797749979*fIn[56])+1.732050807568877*fIn[33]-1.0*fIn[20])*dS; 
    out[17] += vBoundary[0]*((-2.23606797749979*fIn[72])+1.732050807568877*fIn[54]-1.0*fIn[44])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 1.0*dxv[3]*M_PI/_m; 
    out[0] += (2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (2.23606797749979*fIn[34]+1.732050807568877*fIn[15]+fIn[5])*dS; 
    out[4] += (2.23606797749979*fIn[45]+1.732050807568877*fIn[21]+fIn[11])*dS; 
    out[5] += (2.23606797749979*fIn[46]+1.732050807568877*fIn[22]+fIn[12])*dS; 
    out[6] += (2.23606797749979*fIn[55]+1.732050807568877*fIn[32]+fIn[19])*dS; 
    out[7] += (2.23606797749979*fIn[56]+1.732050807568877*fIn[33]+fIn[20])*dS; 
    out[8] += (2.23606797749979*fIn[72]+1.732050807568877*fIn[54]+fIn[44])*dS; 
    out[9] += vBoundary[2]*(2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[10] += vBoundary[2]*(2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[11] += vBoundary[2]*(2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[12] += vBoundary[2]*(2.23606797749979*fIn[34]+1.732050807568877*fIn[15]+fIn[5])*dS; 
    out[13] += vBoundary[2]*(2.23606797749979*fIn[45]+1.732050807568877*fIn[21]+fIn[11])*dS; 
    out[14] += vBoundary[2]*(2.23606797749979*fIn[46]+1.732050807568877*fIn[22]+fIn[12])*dS; 
    out[15] += vBoundary[2]*(2.23606797749979*fIn[55]+1.732050807568877*fIn[32]+fIn[19])*dS; 
    out[16] += vBoundary[2]*(2.23606797749979*fIn[56]+1.732050807568877*fIn[33]+fIn[20])*dS; 
    out[17] += vBoundary[2]*(2.23606797749979*fIn[72]+1.732050807568877*fIn[54]+fIn[44])*dS; 
  }
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 2.0*dxv[2]*M_PI/_m; 
    out[9] += vBoundary[1]*((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[10] += vBoundary[1]*((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[11] += vBoundary[1]*((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[12] += vBoundary[1]*((-2.23606797749979*fIn[41])+1.732050807568877*fIn[16]-1.0*fIn[5])*dS; 
    out[13] += vBoundary[1]*((-2.23606797749979*fIn[47])+1.732050807568877*fIn[25]-1.0*fIn[11])*dS; 
    out[14] += vBoundary[1]*((-2.23606797749979*fIn[48])+1.732050807568877*fIn[26]-1.0*fIn[12])*dS; 
    out[15] += vBoundary[1]*((-2.23606797749979*fIn[60])+1.732050807568877*fIn[35]-1.0*fIn[19])*dS; 
    out[16] += vBoundary[1]*((-2.23606797749979*fIn[61])+1.732050807568877*fIn[36]-1.0*fIn[20])*dS; 
    out[17] += vBoundary[1]*((-2.23606797749979*fIn[73])+1.732050807568877*fIn[57]-1.0*fIn[44])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 2.0*dxv[2]*M_PI/_m; 
    out[9] += vBoundary[3]*(2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[10] += vBoundary[3]*(2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[11] += vBoundary[3]*(2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[12] += vBoundary[3]*(2.23606797749979*fIn[41]+1.732050807568877*fIn[16]+fIn[5])*dS; 
    out[13] += vBoundary[3]*(2.23606797749979*fIn[47]+1.732050807568877*fIn[25]+fIn[11])*dS; 
    out[14] += vBoundary[3]*(2.23606797749979*fIn[48]+1.732050807568877*fIn[26]+fIn[12])*dS; 
    out[15] += vBoundary[3]*(2.23606797749979*fIn[60]+1.732050807568877*fIn[35]+fIn[19])*dS; 
    out[16] += vBoundary[3]*(2.23606797749979*fIn[61]+1.732050807568877*fIn[36]+fIn[20])*dS; 
    out[17] += vBoundary[3]*(2.23606797749979*fIn[73]+1.732050807568877*fIn[57]+fIn[44])*dS; 
  }
 
} 
 
