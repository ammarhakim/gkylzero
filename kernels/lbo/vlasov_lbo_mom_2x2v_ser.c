#include <gkyl_vlasov_lbo_mom_kernels.h> 
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_f_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.5*dxv[3]; 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[5])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.5*dxv[3]; 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[5])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.5*dxv[2]; 
    out[4] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[5] += (1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[6] += (1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[7] += (1.732050807568877*fIn[12]-1.0*fIn[5])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.5*dxv[2]; 
    out[4] += (1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[5] += (1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[6] += (1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[7] += (1.732050807568877*fIn[12]+fIn[5])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_vf_p1(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.5*dxv[3]; 
      out[0] += (0.8660254037844386*dxv[2]*fIn[3]-0.5*fIn[0]*dxv[2])*dS+vBoundary[0]*(1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
      out[1] += (0.8660254037844386*dxv[2]*fIn[6]-0.5*fIn[1]*dxv[2])*dS+vBoundary[0]*(1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
      out[2] += (0.8660254037844386*dxv[2]*fIn[7]-0.5*dxv[2]*fIn[2])*dS+vBoundary[0]*(1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
      out[3] += (0.8660254037844386*dxv[2]*fIn[11]-0.5*dxv[2]*fIn[5])*dS+vBoundary[0]*(1.732050807568877*fIn[11]-1.0*fIn[5])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.5*dxv[3]; 
      out[0] += ((-0.8660254037844386*dxv[2]*fIn[3])-0.5*fIn[0]*dxv[2])*dS+vBoundary[2]*(1.732050807568877*fIn[3]+fIn[0])*dS; 
      out[1] += ((-0.8660254037844386*dxv[2]*fIn[6])-0.5*fIn[1]*dxv[2])*dS+vBoundary[2]*(1.732050807568877*fIn[6]+fIn[1])*dS; 
      out[2] += ((-0.8660254037844386*dxv[2]*fIn[7])-0.5*dxv[2]*fIn[2])*dS+vBoundary[2]*(1.732050807568877*fIn[7]+fIn[2])*dS; 
      out[3] += ((-0.8660254037844386*dxv[2]*fIn[11])-0.5*dxv[2]*fIn[5])*dS+vBoundary[2]*(1.732050807568877*fIn[11]+fIn[5])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.5*dxv[2]; 
      out[0] += (0.8660254037844386*dxv[3]*fIn[4]-0.5*fIn[0]*dxv[3])*dS+vBoundary[1]*(1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
      out[1] += (0.8660254037844386*dxv[3]*fIn[8]-0.5*fIn[1]*dxv[3])*dS+vBoundary[1]*(1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
      out[2] += (0.8660254037844386*dxv[3]*fIn[9]-0.5*fIn[2]*dxv[3])*dS+vBoundary[1]*(1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
      out[3] += (0.8660254037844386*dxv[3]*fIn[12]-0.5*dxv[3]*fIn[5])*dS+vBoundary[1]*(1.732050807568877*fIn[12]-1.0*fIn[5])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.5*dxv[2]; 
      out[0] += ((-0.8660254037844386*dxv[3]*fIn[4])-0.5*fIn[0]*dxv[3])*dS+vBoundary[3]*(1.732050807568877*fIn[4]+fIn[0])*dS; 
      out[1] += ((-0.8660254037844386*dxv[3]*fIn[8])-0.5*fIn[1]*dxv[3])*dS+vBoundary[3]*(1.732050807568877*fIn[8]+fIn[1])*dS; 
      out[2] += ((-0.8660254037844386*dxv[3]*fIn[9])-0.5*fIn[2]*dxv[3])*dS+vBoundary[3]*(1.732050807568877*fIn[9]+fIn[2])*dS; 
      out[3] += ((-0.8660254037844386*dxv[3]*fIn[12])-0.5*dxv[3]*fIn[5])*dS+vBoundary[3]*(1.732050807568877*fIn[12]+fIn[5])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_f_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.5*dxv[3]; 
    out[0] += ((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += ((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += ((-2.23606797749979*fIn[34])+1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
    out[6] += (1.732050807568877*fIn[32]-1.0*fIn[19])*dS; 
    out[7] += (1.732050807568877*fIn[33]-1.0*fIn[20])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.5*dxv[3]; 
    out[0] += (2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (2.23606797749979*fIn[34]+1.732050807568877*fIn[15]+fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS; 
    out[6] += (1.732050807568877*fIn[32]+fIn[19])*dS; 
    out[7] += (1.732050807568877*fIn[33]+fIn[20])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.5*dxv[2]; 
    out[8] += ((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[9] += ((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[10] += ((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[11] += ((-2.23606797749979*fIn[41])+1.732050807568877*fIn[16]-1.0*fIn[5])*dS; 
    out[12] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS; 
    out[13] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS; 
    out[14] += (1.732050807568877*fIn[35]-1.0*fIn[19])*dS; 
    out[15] += (1.732050807568877*fIn[36]-1.0*fIn[20])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.5*dxv[2]; 
    out[8] += (2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[9] += (2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[10] += (2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[11] += (2.23606797749979*fIn[41]+1.732050807568877*fIn[16]+fIn[5])*dS; 
    out[12] += (1.732050807568877*fIn[25]+fIn[11])*dS; 
    out[13] += (1.732050807568877*fIn[26]+fIn[12])*dS; 
    out[14] += (1.732050807568877*fIn[35]+fIn[19])*dS; 
    out[15] += (1.732050807568877*fIn[36]+fIn[20])*dS; 
 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_2x2v_ser_vf_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:   indicator of which velocity grid edge is being considered. 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
 
  if (edge == GKYL_VX_LOWER) {
 
    dS = 0.5*dxv[3]; 
      out[0] += vBoundary[0]*((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
      out[1] += vBoundary[0]*((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
      out[2] += vBoundary[0]*((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
      out[3] += vBoundary[0]*((-2.23606797749979*fIn[34])+1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
      out[4] += vBoundary[0]*(1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
      out[5] += vBoundary[0]*(1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
      out[6] += vBoundary[0]*(1.732050807568877*fIn[32]-1.0*fIn[19])*dS; 
      out[7] += vBoundary[0]*(1.732050807568877*fIn[33]-1.0*fIn[20])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    dS = 0.5*dxv[3]; 
      out[0] += vBoundary[2]*(2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
      out[1] += vBoundary[2]*(2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
      out[2] += vBoundary[2]*(2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
      out[3] += vBoundary[2]*(2.23606797749979*fIn[34]+1.732050807568877*fIn[15]+fIn[5])*dS; 
      out[4] += vBoundary[2]*(1.732050807568877*fIn[21]+fIn[11])*dS; 
      out[5] += vBoundary[2]*(1.732050807568877*fIn[22]+fIn[12])*dS; 
      out[6] += vBoundary[2]*(1.732050807568877*fIn[32]+fIn[19])*dS; 
      out[7] += vBoundary[2]*(1.732050807568877*fIn[33]+fIn[20])*dS; 
 
  } 
 
  if (edge == GKYL_VY_LOWER) {
 
    dS = 0.5*dxv[2]; 
      out[0] += vBoundary[1]*((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
      out[1] += vBoundary[1]*((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
      out[2] += vBoundary[1]*((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
      out[3] += vBoundary[1]*((-2.23606797749979*fIn[41])+1.732050807568877*fIn[16]-1.0*fIn[5])*dS; 
      out[4] += vBoundary[1]*(1.732050807568877*fIn[25]-1.0*fIn[11])*dS; 
      out[5] += vBoundary[1]*(1.732050807568877*fIn[26]-1.0*fIn[12])*dS; 
      out[6] += vBoundary[1]*(1.732050807568877*fIn[35]-1.0*fIn[19])*dS; 
      out[7] += vBoundary[1]*(1.732050807568877*fIn[36]-1.0*fIn[20])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    dS = 0.5*dxv[2]; 
      out[0] += vBoundary[3]*(2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
      out[1] += vBoundary[3]*(2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS; 
      out[2] += vBoundary[3]*(2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS; 
      out[3] += vBoundary[3]*(2.23606797749979*fIn[41]+1.732050807568877*fIn[16]+fIn[5])*dS; 
      out[4] += vBoundary[3]*(1.732050807568877*fIn[25]+fIn[11])*dS; 
      out[5] += vBoundary[3]*(1.732050807568877*fIn[26]+fIn[12])*dS; 
      out[6] += vBoundary[3]*(1.732050807568877*fIn[35]+fIn[19])*dS; 
      out[7] += vBoundary[3]*(1.732050807568877*fIn[36]+fIn[20])*dS; 
 
  } 
 
}
 
