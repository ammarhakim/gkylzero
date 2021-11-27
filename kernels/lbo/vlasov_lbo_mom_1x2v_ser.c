#include <gkyl_vlasov_lbo_mom_kernels.h> 
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_f_p1(const int *idx, const int *atLower, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
  int eval = 0; 
 
  if (idx[3] == 0) { 
    dS = 0.5*dxv[2]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[0]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[0] += (1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
      out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
 
    } else if (eval == 1) {
 
      out[0] += (1.732050807568877*fIn[2]+fIn[0])*dS; 
      out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS; 
 
    } 
  } 
 
  if (idx[3] == 1) { 
    dS = 0.5*dxv[1]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[1]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[2] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
      out[3] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
 
    } else if (eval == 1) {
 
      out[2] += (1.732050807568877*fIn[3]+fIn[0])*dS; 
      out[3] += (1.732050807568877*fIn[5]+fIn[1])*dS; 
 
    } 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_vf_p1(const int *idx, const int *atLower, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
  int eval = 0; 
 
  if (idx[3] == 0) { 
    dS = 0.5*dxv[2]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[0]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[0] += (0.8660254037844386*dxv[1]*fIn[2]-0.5*fIn[0]*dxv[1])*dS+vBoundary[0]*(1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
      out[1] += (0.8660254037844386*dxv[1]*fIn[4]-0.5*dxv[1]*fIn[1])*dS+vBoundary[0]*(1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
 
    } else if (eval == 1) {
 
      out[0] += (1.732050807568877*fIn[2]+fIn[0])*vBoundary[2]*dS+((-0.8660254037844386*dxv[1]*fIn[2])-0.5*fIn[0]*dxv[1])*dS; 
      out[1] += ((-0.8660254037844386*dxv[1]*fIn[4])-0.5*dxv[1]*fIn[1])*dS+vBoundary[2]*(1.732050807568877*fIn[4]+fIn[1])*dS; 
 
    } 
  } 
 
  if (idx[3] == 1) { 
    dS = 0.5*dxv[1]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[1]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[0] += (0.8660254037844386*dxv[2]*fIn[3]-0.5*fIn[0]*dxv[2])*dS+vBoundary[1]*(1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
      out[1] += (0.8660254037844386*dxv[2]*fIn[5]-0.5*fIn[1]*dxv[2])*dS+vBoundary[1]*(1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
 
    } else if (eval == 1) {
 
      out[0] += (1.732050807568877*fIn[3]+fIn[0])*vBoundary[3]*dS+((-0.8660254037844386*dxv[2]*fIn[3])-0.5*fIn[0]*dxv[2])*dS; 
      out[1] += ((-0.8660254037844386*dxv[2]*fIn[5])-0.5*fIn[1]*dxv[2])*dS+vBoundary[3]*(1.732050807568877*fIn[5]+fIn[1])*dS; 
 
    } 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_f_p2(const int *idx, const int *atLower, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
  int eval = 0; 
 
  if (idx[3] == 0) { 
    dS = 0.5*dxv[2]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[0]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
      out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
      out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
 
    } else if (eval == 1) {
 
      out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
      out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
      out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS; 
 
    } 
  } 
 
  if (idx[3] == 1) { 
    dS = 0.5*dxv[1]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[1]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[3] += ((-2.23606797749979*fIn[9])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
      out[4] += ((-2.23606797749979*fIn[15])+1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
      out[5] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS; 
 
    } else if (eval == 1) {
 
      out[3] += (2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS; 
      out[4] += (2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS; 
      out[5] += (1.732050807568877*fIn[13]+fIn[7])*dS; 
 
    } 
  } 
 
}
 
GKYL_CU_DH void vlasov_boundary_integral_1x2v_ser_vf_p2(const int *idx, const int *atLower, const double *vBoundary, const double *dxv, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:   cell length in each direction. 
  // fIn[NP]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  double dS = 0.0; 
  int eval = 0; 
 
  if (idx[3] == 0) { 
    dS = 0.5*dxv[2]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[0]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[0] += vBoundary[0]*((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
      out[1] += vBoundary[0]*((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
      out[2] += vBoundary[0]*(1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
 
    } else if (eval == 1) {
 
      out[0] += vBoundary[2]*(2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
      out[1] += vBoundary[2]*(2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
      out[2] += vBoundary[2]*(1.732050807568877*fIn[11]+fIn[7])*dS; 
 
    } 
  } 
 
  if (idx[3] == 1) { 
    dS = 0.5*dxv[1]; 
    eval = ((idx[2] == 0) ? -1 : ((idx[2] == atLower[1]) ? 1 : 0)); 
 
    if (eval == -1) {
 
      out[0] += vBoundary[1]*((-2.23606797749979*fIn[9])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
      out[1] += vBoundary[1]*((-2.23606797749979*fIn[15])+1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
      out[2] += vBoundary[1]*(1.732050807568877*fIn[13]-1.0*fIn[7])*dS; 
 
    } else if (eval == 1) {
 
      out[0] += vBoundary[3]*(2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS; 
      out[1] += vBoundary[3]*(2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS; 
      out[2] += vBoundary[3]*(1.732050807568877*fIn[13]+fIn[7])*dS; 
 
    } 
  } 
 
}
 
