#include <gkyl_mom_bcorr_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void mom_bcorr_lbo_gyrokinetic_3x2v_ser_p2(const int *idx, enum gkyl_vel_edge edge, const double *vBoundary, const double *dxv, double _m, const double *fIn, double* GKYL_RESTRICT out) 
{ 
  // edge:      indicator of which velocity grid edge is being considered (VX/VPAR, VY/MU, VZ). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direction. 
  // _m:        species mass. 
  // fIn[112]:   distribution function at lower/upper velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) [first vdim components] and vf^(vmax)_(vmin) [last component]. 
 
  double dS = 0.0; 
 
  dS = 3.141592653589793*dxv[4]/_m; 
 
  if (edge == GKYL_VX_LOWER) {
 
    out[0] += ((-2.23606797749979*fIn[19])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[40])+1.732050807568877*fIn[9]-1.0*fIn[1])*dS; 
    out[2] += ((-2.23606797749979*fIn[41])+1.732050807568877*fIn[10]-1.0*fIn[2])*dS; 
    out[3] += ((-2.23606797749979*fIn[42])+1.732050807568877*fIn[11]-1.0*fIn[3])*dS; 
    out[4] += ((-2.23606797749979*fIn[65])+1.732050807568877*fIn[22]-1.0*fIn[6])*dS; 
    out[5] += ((-2.23606797749979*fIn[66])+1.732050807568877*fIn[23]-1.0*fIn[7])*dS; 
    out[6] += ((-2.23606797749979*fIn[67])+1.732050807568877*fIn[24]-1.0*fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[37]-1.0*fIn[16])*dS; 
    out[8] += (1.732050807568877*fIn[38]-1.0*fIn[17])*dS; 
    out[9] += (1.732050807568877*fIn[39]-1.0*fIn[18])*dS; 
    out[10] += ((-2.23606797749979*fIn[90])+1.732050807568877*fIn[51]-1.0*fIn[21])*dS; 
    out[11] += (1.732050807568877*fIn[59]-1.0*fIn[31])*dS; 
    out[12] += (1.732050807568877*fIn[60]-1.0*fIn[32])*dS; 
    out[13] += (1.732050807568877*fIn[61]-1.0*fIn[33])*dS; 
    out[14] += (1.732050807568877*fIn[62]-1.0*fIn[34])*dS; 
    out[15] += (1.732050807568877*fIn[63]-1.0*fIn[35])*dS; 
    out[16] += (1.732050807568877*fIn[64]-1.0*fIn[36])*dS; 
    out[17] += (1.732050807568877*fIn[87]-1.0*fIn[56])*dS; 
    out[18] += (1.732050807568877*fIn[88]-1.0*fIn[57])*dS; 
    out[19] += (1.732050807568877*fIn[89]-1.0*fIn[58])*dS; 
    out[20] += vBoundary[0]*((-2.23606797749979*fIn[19])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[21] += vBoundary[0]*((-2.23606797749979*fIn[40])+1.732050807568877*fIn[9]-1.0*fIn[1])*dS; 
    out[22] += vBoundary[0]*((-2.23606797749979*fIn[41])+1.732050807568877*fIn[10]-1.0*fIn[2])*dS; 
    out[23] += vBoundary[0]*((-2.23606797749979*fIn[42])+1.732050807568877*fIn[11]-1.0*fIn[3])*dS; 
    out[24] += vBoundary[0]*((-2.23606797749979*fIn[65])+1.732050807568877*fIn[22]-1.0*fIn[6])*dS; 
    out[25] += vBoundary[0]*((-2.23606797749979*fIn[66])+1.732050807568877*fIn[23]-1.0*fIn[7])*dS; 
    out[26] += vBoundary[0]*((-2.23606797749979*fIn[67])+1.732050807568877*fIn[24]-1.0*fIn[8])*dS; 
    out[27] += vBoundary[0]*(1.732050807568877*fIn[37]-1.0*fIn[16])*dS; 
    out[28] += vBoundary[0]*(1.732050807568877*fIn[38]-1.0*fIn[17])*dS; 
    out[29] += vBoundary[0]*(1.732050807568877*fIn[39]-1.0*fIn[18])*dS; 
    out[30] += vBoundary[0]*((-2.23606797749979*fIn[90])+1.732050807568877*fIn[51]-1.0*fIn[21])*dS; 
    out[31] += vBoundary[0]*(1.732050807568877*fIn[59]-1.0*fIn[31])*dS; 
    out[32] += vBoundary[0]*(1.732050807568877*fIn[60]-1.0*fIn[32])*dS; 
    out[33] += vBoundary[0]*(1.732050807568877*fIn[61]-1.0*fIn[33])*dS; 
    out[34] += vBoundary[0]*(1.732050807568877*fIn[62]-1.0*fIn[34])*dS; 
    out[35] += vBoundary[0]*(1.732050807568877*fIn[63]-1.0*fIn[35])*dS; 
    out[36] += vBoundary[0]*(1.732050807568877*fIn[64]-1.0*fIn[36])*dS; 
    out[37] += vBoundary[0]*(1.732050807568877*fIn[87]-1.0*fIn[56])*dS; 
    out[38] += vBoundary[0]*(1.732050807568877*fIn[88]-1.0*fIn[57])*dS; 
    out[39] += vBoundary[0]*(1.732050807568877*fIn[89]-1.0*fIn[58])*dS; 
 
  } else if (edge == GKYL_VX_UPPER) {
 
    out[0] += (2.23606797749979*fIn[19]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[40]+1.732050807568877*fIn[9]+fIn[1])*dS; 
    out[2] += (2.23606797749979*fIn[41]+1.732050807568877*fIn[10]+fIn[2])*dS; 
    out[3] += (2.23606797749979*fIn[42]+1.732050807568877*fIn[11]+fIn[3])*dS; 
    out[4] += (2.23606797749979*fIn[65]+1.732050807568877*fIn[22]+fIn[6])*dS; 
    out[5] += (2.23606797749979*fIn[66]+1.732050807568877*fIn[23]+fIn[7])*dS; 
    out[6] += (2.23606797749979*fIn[67]+1.732050807568877*fIn[24]+fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[37]+fIn[16])*dS; 
    out[8] += (1.732050807568877*fIn[38]+fIn[17])*dS; 
    out[9] += (1.732050807568877*fIn[39]+fIn[18])*dS; 
    out[10] += (2.23606797749979*fIn[90]+1.732050807568877*fIn[51]+fIn[21])*dS; 
    out[11] += (1.732050807568877*fIn[59]+fIn[31])*dS; 
    out[12] += (1.732050807568877*fIn[60]+fIn[32])*dS; 
    out[13] += (1.732050807568877*fIn[61]+fIn[33])*dS; 
    out[14] += (1.732050807568877*fIn[62]+fIn[34])*dS; 
    out[15] += (1.732050807568877*fIn[63]+fIn[35])*dS; 
    out[16] += (1.732050807568877*fIn[64]+fIn[36])*dS; 
    out[17] += (1.732050807568877*fIn[87]+fIn[56])*dS; 
    out[18] += (1.732050807568877*fIn[88]+fIn[57])*dS; 
    out[19] += (1.732050807568877*fIn[89]+fIn[58])*dS; 
    out[20] += vBoundary[2]*(2.23606797749979*fIn[19]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[21] += vBoundary[2]*(2.23606797749979*fIn[40]+1.732050807568877*fIn[9]+fIn[1])*dS; 
    out[22] += vBoundary[2]*(2.23606797749979*fIn[41]+1.732050807568877*fIn[10]+fIn[2])*dS; 
    out[23] += vBoundary[2]*(2.23606797749979*fIn[42]+1.732050807568877*fIn[11]+fIn[3])*dS; 
    out[24] += vBoundary[2]*(2.23606797749979*fIn[65]+1.732050807568877*fIn[22]+fIn[6])*dS; 
    out[25] += vBoundary[2]*(2.23606797749979*fIn[66]+1.732050807568877*fIn[23]+fIn[7])*dS; 
    out[26] += vBoundary[2]*(2.23606797749979*fIn[67]+1.732050807568877*fIn[24]+fIn[8])*dS; 
    out[27] += vBoundary[2]*(1.732050807568877*fIn[37]+fIn[16])*dS; 
    out[28] += vBoundary[2]*(1.732050807568877*fIn[38]+fIn[17])*dS; 
    out[29] += vBoundary[2]*(1.732050807568877*fIn[39]+fIn[18])*dS; 
    out[30] += vBoundary[2]*(2.23606797749979*fIn[90]+1.732050807568877*fIn[51]+fIn[21])*dS; 
    out[31] += vBoundary[2]*(1.732050807568877*fIn[59]+fIn[31])*dS; 
    out[32] += vBoundary[2]*(1.732050807568877*fIn[60]+fIn[32])*dS; 
    out[33] += vBoundary[2]*(1.732050807568877*fIn[61]+fIn[33])*dS; 
    out[34] += vBoundary[2]*(1.732050807568877*fIn[62]+fIn[34])*dS; 
    out[35] += vBoundary[2]*(1.732050807568877*fIn[63]+fIn[35])*dS; 
    out[36] += vBoundary[2]*(1.732050807568877*fIn[64]+fIn[36])*dS; 
    out[37] += vBoundary[2]*(1.732050807568877*fIn[87]+fIn[56])*dS; 
    out[38] += vBoundary[2]*(1.732050807568877*fIn[88]+fIn[57])*dS; 
    out[39] += vBoundary[2]*(1.732050807568877*fIn[89]+fIn[58])*dS; 
  }
 
  dS = 6.283185307179586*dxv[3]/_m; 
 
  if (edge == GKYL_VY_LOWER) {
 
    out[20] += vBoundary[1]*((-2.23606797749979*fIn[20])+1.732050807568877*fIn[5]-1.0*fIn[0])*dS; 
    out[21] += vBoundary[1]*((-2.23606797749979*fIn[47])+1.732050807568877*fIn[12]-1.0*fIn[1])*dS; 
    out[22] += vBoundary[1]*((-2.23606797749979*fIn[48])+1.732050807568877*fIn[13]-1.0*fIn[2])*dS; 
    out[23] += vBoundary[1]*((-2.23606797749979*fIn[49])+1.732050807568877*fIn[14]-1.0*fIn[3])*dS; 
    out[24] += vBoundary[1]*((-2.23606797749979*fIn[80])+1.732050807568877*fIn[25]-1.0*fIn[6])*dS; 
    out[25] += vBoundary[1]*((-2.23606797749979*fIn[81])+1.732050807568877*fIn[26]-1.0*fIn[7])*dS; 
    out[26] += vBoundary[1]*((-2.23606797749979*fIn[82])+1.732050807568877*fIn[27]-1.0*fIn[8])*dS; 
    out[27] += vBoundary[1]*(1.732050807568877*fIn[43]-1.0*fIn[16])*dS; 
    out[28] += vBoundary[1]*(1.732050807568877*fIn[44]-1.0*fIn[17])*dS; 
    out[29] += vBoundary[1]*(1.732050807568877*fIn[45]-1.0*fIn[18])*dS; 
    out[30] += vBoundary[1]*((-2.23606797749979*fIn[103])+1.732050807568877*fIn[52]-1.0*fIn[21])*dS; 
    out[31] += vBoundary[1]*(1.732050807568877*fIn[68]-1.0*fIn[31])*dS; 
    out[32] += vBoundary[1]*(1.732050807568877*fIn[69]-1.0*fIn[32])*dS; 
    out[33] += vBoundary[1]*(1.732050807568877*fIn[70]-1.0*fIn[33])*dS; 
    out[34] += vBoundary[1]*(1.732050807568877*fIn[71]-1.0*fIn[34])*dS; 
    out[35] += vBoundary[1]*(1.732050807568877*fIn[72]-1.0*fIn[35])*dS; 
    out[36] += vBoundary[1]*(1.732050807568877*fIn[73]-1.0*fIn[36])*dS; 
    out[37] += vBoundary[1]*(1.732050807568877*fIn[91]-1.0*fIn[56])*dS; 
    out[38] += vBoundary[1]*(1.732050807568877*fIn[92]-1.0*fIn[57])*dS; 
    out[39] += vBoundary[1]*(1.732050807568877*fIn[93]-1.0*fIn[58])*dS; 
 
  } else if (edge == GKYL_VY_UPPER) {
 
    out[20] += vBoundary[3]*(2.23606797749979*fIn[20]+1.732050807568877*fIn[5]+fIn[0])*dS; 
    out[21] += vBoundary[3]*(2.23606797749979*fIn[47]+1.732050807568877*fIn[12]+fIn[1])*dS; 
    out[22] += vBoundary[3]*(2.23606797749979*fIn[48]+1.732050807568877*fIn[13]+fIn[2])*dS; 
    out[23] += vBoundary[3]*(2.23606797749979*fIn[49]+1.732050807568877*fIn[14]+fIn[3])*dS; 
    out[24] += vBoundary[3]*(2.23606797749979*fIn[80]+1.732050807568877*fIn[25]+fIn[6])*dS; 
    out[25] += vBoundary[3]*(2.23606797749979*fIn[81]+1.732050807568877*fIn[26]+fIn[7])*dS; 
    out[26] += vBoundary[3]*(2.23606797749979*fIn[82]+1.732050807568877*fIn[27]+fIn[8])*dS; 
    out[27] += vBoundary[3]*(1.732050807568877*fIn[43]+fIn[16])*dS; 
    out[28] += vBoundary[3]*(1.732050807568877*fIn[44]+fIn[17])*dS; 
    out[29] += vBoundary[3]*(1.732050807568877*fIn[45]+fIn[18])*dS; 
    out[30] += vBoundary[3]*(2.23606797749979*fIn[103]+1.732050807568877*fIn[52]+fIn[21])*dS; 
    out[31] += vBoundary[3]*(1.732050807568877*fIn[68]+fIn[31])*dS; 
    out[32] += vBoundary[3]*(1.732050807568877*fIn[69]+fIn[32])*dS; 
    out[33] += vBoundary[3]*(1.732050807568877*fIn[70]+fIn[33])*dS; 
    out[34] += vBoundary[3]*(1.732050807568877*fIn[71]+fIn[34])*dS; 
    out[35] += vBoundary[3]*(1.732050807568877*fIn[72]+fIn[35])*dS; 
    out[36] += vBoundary[3]*(1.732050807568877*fIn[73]+fIn[36])*dS; 
    out[37] += vBoundary[3]*(1.732050807568877*fIn[91]+fIn[56])*dS; 
    out[38] += vBoundary[3]*(1.732050807568877*fIn[92]+fIn[57])*dS; 
    out[39] += vBoundary[3]*(1.732050807568877*fIn[93]+fIn[58])*dS; 
  }
 
} 
 
