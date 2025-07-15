#include <gkyl_array_integrate_kernels.h>
GKYL_CU_DH void gkyl_array_integrate_op_none_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += 1.4142135623730951*fIn[2*c]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_none_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += 1.4142135623730951*fIn[3*c]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_none_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += 2.0*fIn[4*c]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_none_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += 2.0*fIn[8*c]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_none_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += 2.8284271247461907*fIn[8*c]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_none_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += 2.8284271247461907*fIn[20*c]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(1.4142135623730951*fIn[2*c])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(1.4142135623730951*fIn[3*c])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(2.0*fIn[4*c])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(2.0*fIn[8*c])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(2.8284271247461907*fIn[8*c])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(2.8284271247461907*fIn[20*c])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c) {
    out[c] += ((fIn[2*c+1]*fIn[2*c+1])+(fIn[2*c]*fIn[2*c]))*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c) {
    out[c] += ((fIn[3*c+2]*fIn[3*c+2])+(fIn[3*c+1]*fIn[3*c+1])+(fIn[3*c]*fIn[3*c]))*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c) {
    out[c] += ((fIn[4*c+3]*fIn[4*c+3])+(fIn[4*c+2]*fIn[4*c+2])+(fIn[4*c+1]*fIn[4*c+1])+(fIn[4*c]*fIn[4*c]))*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c) {
    out[c] += ((fIn[8*c+7]*fIn[8*c+7])+(fIn[8*c+6]*fIn[8*c+6])+(fIn[8*c+5]*fIn[8*c+5])+(fIn[8*c+4]*fIn[8*c+4])+(fIn[8*c+3]*fIn[8*c+3])+(fIn[8*c+2]*fIn[8*c+2])+(fIn[8*c+1]*fIn[8*c+1])+(fIn[8*c]*fIn[8*c]))*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c) {
    out[c] += ((fIn[8*c+7]*fIn[8*c+7])+(fIn[8*c+6]*fIn[8*c+6])+(fIn[8*c+5]*fIn[8*c+5])+(fIn[8*c+4]*fIn[8*c+4])+(fIn[8*c+3]*fIn[8*c+3])+(fIn[8*c+2]*fIn[8*c+2])+(fIn[8*c+1]*fIn[8*c+1])+(fIn[8*c]*fIn[8*c]))*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 
  for (unsigned c=0; c<num_comp; ++c) {
    out[c] += ((fIn[20*c+19]*fIn[20*c+19])+(fIn[20*c+18]*fIn[20*c+18])+(fIn[20*c+17]*fIn[20*c+17])+(fIn[20*c+16]*fIn[20*c+16])+(fIn[20*c+15]*fIn[20*c+15])+(fIn[20*c+14]*fIn[20*c+14])+(fIn[20*c+13]*fIn[20*c+13])+(fIn[20*c+12]*fIn[20*c+12])+(fIn[20*c+11]*fIn[20*c+11])+(fIn[20*c+10]*fIn[20*c+10])+(fIn[20*c+9]*fIn[20*c+9])+(fIn[20*c+8]*fIn[20*c+8])+(fIn[20*c+7]*fIn[20*c+7])+(fIn[20*c+6]*fIn[20*c+6])+(fIn[20*c+5]*fIn[20*c+5])+(fIn[20*c+4]*fIn[20*c+4])+(fIn[20*c+3]*fIn[20*c+3])+(fIn[20*c+2]*fIn[20*c+2])+(fIn[20*c+1]*fIn[20*c+1])+(fIn[20*c]*fIn[20*c]))*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 0.7071067811865475*vol;

  out[0] += (2.0*fIn[0]*fIn[1]*weight[1]+weight[0]*(fIn[1]*fIn[1])+(fIn[0]*fIn[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 0.020203050891044214*vol;

  out[0] += ((22.3606797749979*(fIn[2]*fIn[2])+70.0*fIn[0]*fIn[2]+31.304951684997057*(fIn[1]*fIn[1]))*weight[2]+35.0*weight[0]*(fIn[2]*fIn[2])+62.609903369994115*fIn[1]*weight[1]*fIn[2]+70.0*fIn[0]*fIn[1]*weight[1]+35.0*weight[0]*(fIn[1]*fIn[1])+35.0*(fIn[0]*fIn[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 0.5*vol;

  out[0] += ((2.0*fIn[0]*fIn[3]+2.0*fIn[1]*fIn[2])*weight[3]+weight[0]*(fIn[3]*fIn[3])+(2.0*fIn[1]*weight[2]+2.0*weight[1]*fIn[2])*fIn[3]+2.0*fIn[0]*fIn[2]*weight[2]+weight[0]*(fIn[2]*fIn[2])+2.0*fIn[0]*fIn[1]*weight[1]+weight[0]*(fIn[1]*fIn[1])+(fIn[0]*fIn[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 0.004761904761904762*vol;

  out[0] += (((134.1640786499874*fIn[5]+187.82971010998233*fIn[4]+210.0*fIn[0])*fIn[7]+168.0*fIn[3]*fIn[6]+210.00000000000003*fIn[1]*fIn[5]+187.82971010998233*fIn[2]*fIn[3])*weight[7]+(67.0820393249937*weight[5]+93.91485505499116*weight[4]+105.0*weight[0])*(fIn[7]*fIn[7])+(168.0*fIn[3]*weight[6]+168.0*weight[3]*fIn[6]+210.00000000000003*fIn[1]*weight[5]+210.00000000000003*weight[1]*fIn[5]+187.82971010998233*fIn[2]*weight[3]+187.82971010998233*weight[2]*fIn[3])*fIn[7]+((187.82971010998233*fIn[5]+134.1640786499874*fIn[4]+210.0*fIn[0])*fIn[6]+210.00000000000003*fIn[2]*fIn[4]+187.82971010998233*fIn[1]*fIn[3])*weight[6]+(93.91485505499116*weight[5]+67.0820393249937*weight[4]+105.0*weight[0])*(fIn[6]*fIn[6])+(210.00000000000003*fIn[2]*weight[4]+210.00000000000003*weight[2]*fIn[4]+187.82971010998233*fIn[1]*weight[3]+187.82971010998233*weight[1]*fIn[3])*fIn[6]+(67.0820393249937*(fIn[5]*fIn[5])+210.0*fIn[0]*fIn[5]+93.91485505499116*(fIn[3]*fIn[3])+93.91485505499116*(fIn[2]*fIn[2]))*weight[5]+105.0*weight[0]*(fIn[5]*fIn[5])+(187.82971010998233*fIn[3]*weight[3]+187.82971010998233*fIn[2]*weight[2])*fIn[5]+(67.0820393249937*(fIn[4]*fIn[4])+210.0*fIn[0]*fIn[4]+93.91485505499116*(fIn[3]*fIn[3])+93.91485505499116*(fIn[1]*fIn[1]))*weight[4]+105.0*weight[0]*(fIn[4]*fIn[4])+(187.82971010998233*fIn[3]*weight[3]+187.82971010998233*fIn[1]*weight[1])*fIn[4]+(210.0*fIn[0]*fIn[3]+210.0*fIn[1]*fIn[2])*weight[3]+105.0*weight[0]*(fIn[3]*fIn[3])+(210.0*fIn[1]*weight[2]+210.0*weight[1]*fIn[2])*fIn[3]+210.0*fIn[0]*fIn[2]*weight[2]+105.0*weight[0]*(fIn[2]*fIn[2])+210.0*fIn[0]*fIn[1]*weight[1]+105.0*weight[0]*(fIn[1]*fIn[1])+105.0*(fIn[0]*fIn[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 0.3535533905932737*vol;

  out[0] += ((2.0*fIn[0]*fIn[7]+2.0*fIn[1]*fIn[6]+2.0*fIn[2]*fIn[5]+2.0*fIn[3]*fIn[4])*weight[7]+weight[0]*(fIn[7]*fIn[7])+(2.0*fIn[1]*weight[6]+2.0*weight[1]*fIn[6]+2.0*fIn[2]*weight[5]+2.0*weight[2]*fIn[5]+2.0*fIn[3]*weight[4]+2.0*weight[3]*fIn[4])*fIn[7]+(2.0*fIn[0]*fIn[6]+2.0*fIn[4]*fIn[5]+2.0*fIn[2]*fIn[3])*weight[6]+weight[0]*(fIn[6]*fIn[6])+(2.0*fIn[4]*weight[5]+2.0*weight[4]*fIn[5]+2.0*fIn[2]*weight[3]+2.0*weight[2]*fIn[3])*fIn[6]+(2.0*fIn[0]*fIn[5]+2.0*fIn[1]*fIn[3])*weight[5]+weight[0]*(fIn[5]*fIn[5])+(2.0*fIn[1]*weight[3]+2.0*weight[1]*fIn[3])*fIn[5]+(2.0*fIn[0]*fIn[4]+2.0*fIn[1]*fIn[2])*weight[4]+weight[0]*(fIn[4]*fIn[4])+(2.0*fIn[1]*weight[2]+2.0*weight[1]*fIn[2])*fIn[4]+2.0*fIn[0]*fIn[3]*weight[3]+weight[0]*(fIn[3]*fIn[3])+2.0*fIn[0]*fIn[2]*weight[2]+weight[0]*(fIn[2]*fIn[2])+2.0*fIn[0]*fIn[1]*weight[1]+weight[0]*(fIn[1]*fIn[1])+(fIn[0]*fIn[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 6.734350297014737e-4*vol;

  out[0] += (((670.8203932499371*fIn[9]+939.1485505499119*fIn[8]+939.1485505499119*fIn[7]+1050.0*fIn[0])*fIn[19]+(751.3188404399293*fIn[17]+840.0*fIn[6])*fIn[18]+840.0*fIn[5]*fIn[17]+(670.8203932499371*fIn[15]+939.1485505499119*fIn[12]+1050.0000000000002*fIn[1])*fIn[16]+(939.1485505499119*fIn[11]+1050.0000000000002*fIn[2])*fIn[15]+840.0000000000001*fIn[10]*fIn[14]+840.0000000000001*fIn[10]*fIn[13]+939.1485505499119*fIn[3]*fIn[10]+1050.0*fIn[4]*fIn[9]+939.1485505499119*fIn[5]*fIn[6])*weight[19]+(335.41019662496853*weight[9]+469.57427527495594*weight[8]+469.57427527495594*weight[7]+525.0*weight[0])*(fIn[19]*fIn[19])+((751.3188404399293*fIn[17]+840.0*fIn[6])*weight[18]+(751.3188404399293*weight[17]+840.0*weight[6])*fIn[18]+840.0*fIn[5]*weight[17]+840.0*weight[5]*fIn[17]+(670.8203932499371*fIn[15]+939.1485505499119*fIn[12]+1050.0000000000002*fIn[1])*weight[16]+(670.8203932499371*weight[15]+939.1485505499119*weight[12]+1050.0000000000002*weight[1])*fIn[16]+(939.1485505499119*fIn[11]+1050.0000000000002*fIn[2])*weight[15]+(939.1485505499119*weight[11]+1050.0000000000002*weight[2])*fIn[15]+840.0000000000001*fIn[10]*weight[14]+840.0000000000001*weight[10]*fIn[14]+840.0000000000001*fIn[10]*weight[13]+840.0000000000001*weight[10]*fIn[13]+939.1485505499119*fIn[3]*weight[10]+939.1485505499119*weight[3]*fIn[10]+1050.0*fIn[4]*weight[9]+1050.0*weight[4]*fIn[9]+939.1485505499119*fIn[5]*weight[6]+939.1485505499119*weight[5]*fIn[6])*fIn[19]+((939.1485505499119*fIn[9]+670.8203932499371*fIn[8]+939.1485505499119*fIn[7]+1050.0*fIn[0])*fIn[18]+840.0*fIn[4]*fIn[17]+840.0000000000001*fIn[10]*fIn[16]+939.1485505499119*fIn[14]*fIn[15]+(670.8203932499371*fIn[12]+1050.0000000000002*fIn[1])*fIn[14]+939.1485505499119*fIn[12]*fIn[13]+1050.0000000000002*fIn[3]*fIn[12]+840.0000000000001*fIn[10]*fIn[11]+939.1485505499119*fIn[2]*fIn[10]+1050.0*fIn[5]*fIn[8]+939.1485505499119*fIn[4]*fIn[6])*weight[18]+(469.57427527495594*weight[9]+335.41019662496853*weight[8]+469.57427527495594*weight[7]+525.0*weight[0])*(fIn[18]*fIn[18])+(840.0*fIn[4]*weight[17]+840.0*weight[4]*fIn[17]+840.0000000000001*fIn[10]*weight[16]+840.0000000000001*weight[10]*fIn[16]+939.1485505499119*fIn[14]*weight[15]+939.1485505499119*weight[14]*fIn[15]+(670.8203932499371*fIn[12]+1050.0000000000002*fIn[1])*weight[14]+(670.8203932499371*weight[12]+1050.0000000000002*weight[1])*fIn[14]+939.1485505499119*fIn[12]*weight[13]+939.1485505499119*weight[12]*fIn[13]+1050.0000000000002*fIn[3]*weight[12]+1050.0000000000002*weight[3]*fIn[12]+840.0000000000001*fIn[10]*weight[11]+840.0000000000001*weight[10]*fIn[11]+939.1485505499119*fIn[2]*weight[10]+939.1485505499119*weight[2]*fIn[10]+1050.0*fIn[5]*weight[8]+1050.0*weight[5]*fIn[8]+939.1485505499119*fIn[4]*weight[6]+939.1485505499119*weight[4]*fIn[6])*fIn[18]+((939.1485505499119*fIn[9]+939.1485505499119*fIn[8]+670.8203932499371*fIn[7]+1050.0*fIn[0])*fIn[17]+939.1485505499119*fIn[13]*fIn[16]+840.0000000000001*fIn[10]*fIn[15]+939.1485505499119*fIn[11]*fIn[14]+(670.8203932499371*fIn[11]+1050.0000000000002*fIn[2])*fIn[13]+840.0000000000001*fIn[10]*fIn[12]+1050.0000000000002*fIn[3]*fIn[11]+939.1485505499119*fIn[1]*fIn[10]+1050.0*fIn[6]*fIn[7]+939.1485505499119*fIn[4]*fIn[5])*weight[17]+(469.57427527495594*weight[9]+469.57427527495594*weight[8]+335.41019662496853*weight[7]+525.0*weight[0])*(fIn[17]*fIn[17])+(939.1485505499119*fIn[13]*weight[16]+939.1485505499119*weight[13]*fIn[16]+840.0000000000001*fIn[10]*weight[15]+840.0000000000001*weight[10]*fIn[15]+939.1485505499119*fIn[11]*weight[14]+939.1485505499119*weight[11]*fIn[14]+(670.8203932499371*fIn[11]+1050.0000000000002*fIn[2])*weight[13]+(670.8203932499371*weight[11]+1050.0000000000002*weight[2])*fIn[13]+840.0000000000001*fIn[10]*weight[12]+840.0000000000001*weight[10]*fIn[12]+1050.0000000000002*fIn[3]*weight[11]+1050.0000000000002*weight[3]*fIn[11]+939.1485505499119*fIn[1]*weight[10]+939.1485505499119*weight[1]*fIn[10]+1050.0*fIn[6]*weight[7]+1050.0*weight[6]*fIn[7]+939.1485505499119*fIn[4]*weight[5]+939.1485505499119*weight[4]*fIn[5])*fIn[17]+((670.8203932499371*fIn[9]+939.1485505499119*fIn[8]+1050.0*fIn[0])*fIn[16]+1050.0*fIn[4]*fIn[15]+840.0*fIn[6]*fIn[14]+939.1485505499116*fIn[5]*fIn[10]+1050.0000000000002*fIn[2]*fIn[9]+939.1485505499116*fIn[3]*fIn[6])*weight[16]+(335.41019662496853*weight[9]+469.57427527495594*weight[8]+525.0*weight[0])*(fIn[16]*fIn[16])+(1050.0*fIn[4]*weight[15]+1050.0*weight[4]*fIn[15]+840.0*fIn[6]*weight[14]+840.0*weight[6]*fIn[14]+939.1485505499116*fIn[5]*weight[10]+939.1485505499116*weight[5]*fIn[10]+1050.0000000000002*fIn[2]*weight[9]+1050.0000000000002*weight[2]*fIn[9]+939.1485505499116*fIn[3]*weight[6]+939.1485505499116*weight[3]*fIn[6])*fIn[16]+((670.8203932499371*fIn[9]+939.1485505499119*fIn[7]+1050.0*fIn[0])*fIn[15]+840.0*fIn[5]*fIn[13]+939.1485505499116*fIn[6]*fIn[10]+1050.0000000000002*fIn[1]*fIn[9]+939.1485505499116*fIn[3]*fIn[5])*weight[15]+(335.41019662496853*weight[9]+469.57427527495594*weight[7]+525.0*weight[0])*(fIn[15]*fIn[15])+(840.0*fIn[5]*weight[13]+840.0*weight[5]*fIn[13]+939.1485505499116*fIn[6]*weight[10]+939.1485505499116*weight[6]*fIn[10]+1050.0000000000002*fIn[1]*weight[9]+1050.0000000000002*weight[1]*fIn[9]+939.1485505499116*fIn[3]*weight[5]+939.1485505499116*weight[3]*fIn[5])*fIn[15]+((939.1485505499119*fIn[9]+670.8203932499371*fIn[8]+1050.0*fIn[0])*fIn[14]+1050.0*fIn[5]*fIn[12]+939.1485505499116*fIn[4]*fIn[10]+1050.0000000000002*fIn[3]*fIn[8]+939.1485505499116*fIn[2]*fIn[6])*weight[14]+(469.57427527495594*weight[9]+335.41019662496853*weight[8]+525.0*weight[0])*(fIn[14]*fIn[14])+(1050.0*fIn[5]*weight[12]+1050.0*weight[5]*fIn[12]+939.1485505499116*fIn[4]*weight[10]+939.1485505499116*weight[4]*fIn[10]+1050.0000000000002*fIn[3]*weight[8]+1050.0000000000002*weight[3]*fIn[8]+939.1485505499116*fIn[2]*weight[6]+939.1485505499116*weight[2]*fIn[6])*fIn[14]+((939.1485505499119*fIn[9]+670.8203932499371*fIn[7]+1050.0*fIn[0])*fIn[13]+1050.0*fIn[6]*fIn[11]+939.1485505499116*fIn[4]*fIn[10]+1050.0000000000002*fIn[3]*fIn[7]+939.1485505499116*fIn[1]*fIn[5])*weight[13]+(469.57427527495594*weight[9]+335.41019662496853*weight[7]+525.0*weight[0])*(fIn[13]*fIn[13])+(1050.0*fIn[6]*weight[11]+1050.0*weight[6]*fIn[11]+939.1485505499116*fIn[4]*weight[10]+939.1485505499116*weight[4]*fIn[10]+1050.0000000000002*fIn[3]*weight[7]+1050.0000000000002*weight[3]*fIn[7]+939.1485505499116*fIn[1]*weight[5]+939.1485505499116*weight[1]*fIn[5])*fIn[13]+((670.8203932499371*fIn[8]+939.1485505499119*fIn[7]+1050.0*fIn[0])*fIn[12]+840.0*fIn[4]*fIn[11]+939.1485505499116*fIn[6]*fIn[10]+1050.0000000000002*fIn[1]*fIn[8]+939.1485505499116*fIn[2]*fIn[4])*weight[12]+(335.41019662496853*weight[8]+469.57427527495594*weight[7]+525.0*weight[0])*(fIn[12]*fIn[12])+(840.0*fIn[4]*weight[11]+840.0*weight[4]*fIn[11]+939.1485505499116*fIn[6]*weight[10]+939.1485505499116*weight[6]*fIn[10]+1050.0000000000002*fIn[1]*weight[8]+1050.0000000000002*weight[1]*fIn[8]+939.1485505499116*fIn[2]*weight[4]+939.1485505499116*weight[2]*fIn[4])*fIn[12]+((939.1485505499119*fIn[8]+670.8203932499371*fIn[7]+1050.0*fIn[0])*fIn[11]+939.1485505499116*fIn[5]*fIn[10]+1050.0000000000002*fIn[2]*fIn[7]+939.1485505499116*fIn[1]*fIn[4])*weight[11]+(469.57427527495594*weight[8]+335.41019662496853*weight[7]+525.0*weight[0])*(fIn[11]*fIn[11])+(939.1485505499116*fIn[5]*weight[10]+939.1485505499116*weight[5]*fIn[10]+1050.0000000000002*fIn[2]*weight[7]+1050.0000000000002*weight[2]*fIn[7]+939.1485505499116*fIn[1]*weight[4]+939.1485505499116*weight[1]*fIn[4])*fIn[11]+((939.1485505499119*fIn[9]+939.1485505499119*fIn[8]+939.1485505499119*fIn[7]+1050.0*fIn[0])*fIn[10]+1050.0*fIn[1]*fIn[6]+1050.0*fIn[2]*fIn[5]+1050.0*fIn[3]*fIn[4])*weight[10]+(469.57427527495594*weight[9]+469.57427527495594*weight[8]+469.57427527495594*weight[7]+525.0*weight[0])*(fIn[10]*fIn[10])+(1050.0*fIn[1]*weight[6]+1050.0*weight[1]*fIn[6]+1050.0*fIn[2]*weight[5]+1050.0*weight[2]*fIn[5]+1050.0*fIn[3]*weight[4]+1050.0*weight[3]*fIn[4])*fIn[10]+(335.41019662496853*(fIn[9]*fIn[9])+1050.0*fIn[0]*fIn[9]+469.57427527495594*(fIn[6]*fIn[6])+469.57427527495594*(fIn[5]*fIn[5])+469.57427527495594*(fIn[3]*fIn[3]))*weight[9]+525.0*weight[0]*(fIn[9]*fIn[9])+(939.1485505499119*fIn[6]*weight[6]+939.1485505499119*fIn[5]*weight[5]+939.1485505499119*fIn[3]*weight[3])*fIn[9]+(335.41019662496853*(fIn[8]*fIn[8])+1050.0*fIn[0]*fIn[8]+469.57427527495594*(fIn[6]*fIn[6])+469.57427527495594*(fIn[4]*fIn[4])+469.57427527495594*(fIn[2]*fIn[2]))*weight[8]+525.0*weight[0]*(fIn[8]*fIn[8])+(939.1485505499119*fIn[6]*weight[6]+939.1485505499119*fIn[4]*weight[4]+939.1485505499119*fIn[2]*weight[2])*fIn[8]+(335.41019662496853*(fIn[7]*fIn[7])+1050.0*fIn[0]*fIn[7]+469.57427527495594*(fIn[5]*fIn[5])+469.57427527495594*(fIn[4]*fIn[4])+469.57427527495594*(fIn[1]*fIn[1]))*weight[7]+525.0*weight[0]*(fIn[7]*fIn[7])+(939.1485505499119*fIn[5]*weight[5]+939.1485505499119*fIn[4]*weight[4]+939.1485505499119*fIn[1]*weight[1])*fIn[7]+(1050.0*fIn[0]*fIn[6]+1050.0*fIn[4]*fIn[5]+1050.0*fIn[2]*fIn[3])*weight[6]+525.0*weight[0]*(fIn[6]*fIn[6])+(1050.0*fIn[4]*weight[5]+1050.0*weight[4]*fIn[5]+1050.0*fIn[2]*weight[3]+1050.0*weight[2]*fIn[3])*fIn[6]+(1050.0*fIn[0]*fIn[5]+1050.0*fIn[1]*fIn[3])*weight[5]+525.0*weight[0]*(fIn[5]*fIn[5])+(1050.0*fIn[1]*weight[3]+1050.0*weight[1]*fIn[3])*fIn[5]+(1050.0*fIn[0]*fIn[4]+1050.0*fIn[1]*fIn[2])*weight[4]+525.0*weight[0]*(fIn[4]*fIn[4])+(1050.0*fIn[1]*weight[2]+1050.0*weight[1]*fIn[2])*fIn[4]+1050.0*fIn[0]*fIn[3]*weight[3]+525.0*weight[0]*(fIn[3]*fIn[3])+1050.0*fIn[0]*fIn[2]*weight[2]+525.0*weight[0]*(fIn[2]*fIn[2])+1050.0*fIn[0]*fIn[1]*weight[1]+525.0*weight[0]*(fIn[1]*fIn[1])+525.0*(fIn[0]*fIn[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_1x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/dxSq[0];

  out[0] += (fIn[1]*fIn[1])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_1x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/dxSq[0];

  out[0] += (5.0*(fIn[2]*fIn[2])+(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]);

  out[0] += ((dxSq[1]+dxSq[0])*(fIn[3]*fIn[3])+dxSq[0]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]);

  out[0] += ((dxSq[1]+5.0*dxSq[0])*(fIn[7]*fIn[7])+(5.0*dxSq[1]+dxSq[0])*(fIn[6]*fIn[6])+5.0*dxSq[0]*(fIn[5]*fIn[5])+5.0*dxSq[1]*(fIn[4]*fIn[4])+(dxSq[1]+dxSq[0])*(fIn[3]*fIn[3])+dxSq[0]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]*dxSq[2]);

  out[0] += (((dxSq[1]+dxSq[0])*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[7]*fIn[7])+(dxSq[0]*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[6]*fIn[6])+(dxSq[1]*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[5]*fIn[5])+(dxSq[1]+dxSq[0])*dxSq[2]*(fIn[4]*fIn[4])+dxSq[0]*dxSq[1]*(fIn[3]*fIn[3])+dxSq[0]*dxSq[2]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1])*dxSq[2])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]*dxSq[2]);

  out[0] += (((dxSq[1]+dxSq[0])*dxSq[2]+5.0*dxSq[0]*dxSq[1])*(fIn[19]*fIn[19])+((dxSq[1]+5.0*dxSq[0])*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[18]*fIn[18])+((5.0*dxSq[1]+dxSq[0])*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[17]*fIn[17])+(dxSq[0]*dxSq[2]+5.0*dxSq[0]*dxSq[1])*(fIn[16]*fIn[16])+(dxSq[1]*dxSq[2]+5.0*dxSq[0]*dxSq[1])*(fIn[15]*fIn[15])+(5.0*dxSq[0]*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[14]*fIn[14])+(5.0*dxSq[1]*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[13]*fIn[13])+(dxSq[1]+5.0*dxSq[0])*dxSq[2]*(fIn[12]*fIn[12])+(5.0*dxSq[1]+dxSq[0])*dxSq[2]*(fIn[11]*fIn[11])+((dxSq[1]+dxSq[0])*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[10]*fIn[10])+5.0*dxSq[0]*dxSq[1]*(fIn[9]*fIn[9])+5.0*dxSq[0]*dxSq[2]*(fIn[8]*fIn[8])+5.0*dxSq[1]*dxSq[2]*(fIn[7]*fIn[7])+(dxSq[0]*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[6]*fIn[6])+(dxSq[1]*dxSq[2]+dxSq[0]*dxSq[1])*(fIn[5]*fIn[5])+(dxSq[1]+dxSq[0])*dxSq[2]*(fIn[4]*fIn[4])+dxSq[0]*dxSq[1]*(fIn[3]*fIn[3])+dxSq[0]*dxSq[2]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1])*dxSq[2])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]);

  out[0] += ((dxSq[1]+dxSq[0])*(fIn[3]*fIn[3])+dxSq[0]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]);

  out[0] += ((dxSq[1]+5.0*dxSq[0])*(fIn[7]*fIn[7])+(5.0*dxSq[1]+dxSq[0])*(fIn[6]*fIn[6])+5.0*dxSq[0]*(fIn[5]*fIn[5])+5.0*dxSq[1]*(fIn[4]*fIn[4])+(dxSq[1]+dxSq[0])*(fIn[3]*fIn[3])+dxSq[0]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]);

  out[0] += ((dxSq[1]+dxSq[0])*(fIn[7]*fIn[7])+dxSq[0]*(fIn[6]*fIn[6])+dxSq[1]*(fIn[5]*fIn[5])+(dxSq[1]+dxSq[0])*(fIn[4]*fIn[4])+dxSq[0]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = (12.0*vol)/(dxSq[0]*dxSq[1]);

  out[0] += ((dxSq[1]+dxSq[0])*(fIn[19]*fIn[19])+(dxSq[1]+5.0*dxSq[0])*(fIn[18]*fIn[18])+(5.0*dxSq[1]+dxSq[0])*(fIn[17]*fIn[17])+dxSq[0]*(fIn[16]*fIn[16])+dxSq[1]*(fIn[15]*fIn[15])+5.0*dxSq[0]*(fIn[14]*fIn[14])+5.0*dxSq[1]*(fIn[13]*fIn[13])+(dxSq[1]+5.0*dxSq[0])*(fIn[12]*fIn[12])+(5.0*dxSq[1]+dxSq[0])*(fIn[11]*fIn[11])+(dxSq[1]+dxSq[0])*(fIn[10]*fIn[10])+5.0*dxSq[0]*(fIn[8]*fIn[8])+5.0*dxSq[1]*(fIn[7]*fIn[7])+dxSq[0]*(fIn[6]*fIn[6])+dxSq[1]*(fIn[5]*fIn[5])+(dxSq[1]+dxSq[0])*(fIn[4]*fIn[4])+dxSq[0]*(fIn[2]*fIn[2])+dxSq[1]*(fIn[1]*fIn[1]))*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double *epsxx = &weight[0];

  double rdx00 = 4.0/dxSq[0];

  out[0] += 0.5*(3.0*epsxx[0]*(fIn[3]*fIn[3])+6.0*fIn[1]*epsxx[2]*fIn[3]+3.0*epsxx[0]*(fIn[1]*fIn[1]))*rdx00*vol; 
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double *epsxx = &weight[0];

  double rdx00 = 4.0/dxSq[0];

  out[0] += 0.014285714285714285*((67.0820393249937*epsxx[5]+105.0*epsxx[0])*(fIn[7]*fIn[7])+(469.57427527495594*fIn[4]*epsxx[7]+420.0*epsxx[3]*fIn[6]+210.00000000000003*fIn[1]*epsxx[5]+187.82971010998233*epsxx[2]*fIn[3])*fIn[7]+420.0*fIn[3]*fIn[6]*epsxx[7]+(469.57427527495594*epsxx[5]+469.57427527495594*epsxx[4]+525.0*epsxx[0])*(fIn[6]*fIn[6])+(939.1485505499119*fIn[4]*epsxx[6]+1050.0000000000002*epsxx[2]*fIn[4]+469.5742752749558*epsxx[1]*fIn[3]+469.5742752749558*fIn[1]*epsxx[3])*fIn[6]+93.91485505499116*(fIn[3]*fIn[3])*epsxx[5]+(469.57427527495594*epsxx[4]+525.0*epsxx[0])*(fIn[4]*fIn[4])+(469.57427527495594*epsxx[3]*fIn[3]+469.57427527495594*epsxx[1]*fIn[1])*fIn[4]+105.0*epsxx[0]*(fIn[3]*fIn[3])+210.0*fIn[1]*epsxx[2]*fIn[3]+105.0*epsxx[0]*(fIn[1]*fIn[1]))*rdx00*vol; 
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double *epsxx = &weight[0];
  const double *epsxy = &weight[8];
  const double *epsyy = &weight[16];

  double rdx00 = 4.0/dxSq[0];
  double rdx01 = 4.0/sqrt(dxSq[0]*dxSq[1]);
  double rdx11 = 4.0/dxSq[1];

  out[0] += 0.25*((4.242640687119286*epsyy[0]*(fIn[7]*fIn[7])+(8.485281374238571*epsyy[1]*fIn[6]+8.485281374238571*fIn[2]*epsyy[5]+8.485281374238571*epsyy[3]*fIn[4])*fIn[7]+4.242640687119286*epsyy[0]*(fIn[6]*fIn[6])+(8.485281374238571*fIn[4]*epsyy[5]+8.485281374238571*fIn[2]*epsyy[3])*fIn[6]+4.242640687119286*epsyy[0]*(fIn[4]*fIn[4])+8.485281374238571*epsyy[1]*fIn[2]*fIn[4]+4.242640687119286*epsyy[0]*(fIn[2]*fIn[2]))*rdx11+(8.485281374238571*epsxy[4]*(fIn[7]*fIn[7])+(16.970562748477146*fIn[4]*epsxy[7]+8.485281374238571*epsxy[2]*fIn[6]+8.485281374238571*fIn[2]*epsxy[6]+8.485281374238571*epsxy[1]*fIn[5]+8.485281374238571*fIn[1]*epsxy[5])*fIn[7]+(8.485281374238571*fIn[4]*epsxy[6]+8.485281374238571*epsxy[0]*fIn[5]+8.485281374238571*fIn[1]*epsxy[3])*fIn[6]+(8.485281374238571*fIn[4]*epsxy[5]+8.485281374238571*fIn[2]*epsxy[3])*fIn[5]+8.485281374238571*epsxy[4]*(fIn[4]*fIn[4])+(8.485281374238571*epsxy[2]*fIn[2]+8.485281374238571*epsxy[1]*fIn[1])*fIn[4]+8.485281374238571*epsxy[0]*fIn[1]*fIn[2])*rdx01+(4.242640687119286*epsxx[0]*(fIn[7]*fIn[7])+(8.485281374238571*fIn[1]*epsxx[6]+8.485281374238571*epsxx[2]*fIn[5]+8.485281374238571*epsxx[3]*fIn[4])*fIn[7]+8.485281374238571*fIn[4]*fIn[5]*epsxx[6]+4.242640687119286*epsxx[0]*(fIn[5]*fIn[5])+8.485281374238571*fIn[1]*epsxx[3]*fIn[5]+4.242640687119286*epsxx[0]*(fIn[4]*fIn[4])+8.485281374238571*fIn[1]*epsxx[2]*fIn[4]+4.242640687119286*epsxx[0]*(fIn[1]*fIn[1]))*rdx00)*vol; 
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double *epsxx = &weight[0];
  const double *epsxy = &weight[20];
  const double *epsyy = &weight[40];

  double rdx00 = 4.0/dxSq[0];
  double rdx01 = 4.0/sqrt(dxSq[0]*dxSq[1]);
  double rdx11 = 4.0/dxSq[1];

  out[0] += 0.0014285714285714286*(((474.3416490252571*epsyy[9]+664.0783086353599*epsyy[7]+742.462120245875*epsyy[0])*(fIn[19]*fIn[19])+(3320.3915431767996*fIn[8]*epsyy[19]+(2656.3132345414406*epsyy[17]+2969.848480983501*epsyy[6])*fIn[18]+1187.9393923934006*epsyy[5]*fIn[17]+(948.6832980505144*epsyy[15]+1484.9242404917502*epsyy[1])*fIn[16]+3320.3915431767996*fIn[12]*epsyy[16]+(1328.15661727072*fIn[11]+1484.9242404917502*fIn[2])*epsyy[15]+2969.8484809835013*epsyy[10]*fIn[14]+1187.9393923934006*fIn[10]*epsyy[13]+1328.15661727072*epsyy[3]*fIn[10]+1484.9242404917502*fIn[4]*epsyy[9]+1328.15661727072*epsyy[5]*fIn[6])*fIn[19]+((2656.3132345414406*fIn[17]+2969.848480983501*fIn[6])*fIn[18]+3320.3915431767996*fIn[12]*fIn[16]+2969.8484809835013*fIn[10]*fIn[14])*epsyy[19]+(3320.3915431767996*epsyy[9]+3320.3915431767996*epsyy[8]+3320.3915431767996*epsyy[7]+3712.3106012293747*epsyy[0])*(fIn[18]*fIn[18])+(6640.783086353601*fIn[8]*epsyy[18]+2969.848480983501*epsyy[4]*fIn[17]+2969.848480983501*fIn[4]*epsyy[17]+2969.8484809835013*epsyy[10]*fIn[16]+2969.8484809835013*fIn[10]*epsyy[16]+6640.783086353601*fIn[14]*epsyy[15]+(6640.783086353601*epsyy[12]+7424.621202458751*epsyy[1])*fIn[14]+6640.783086353601*fIn[12]*epsyy[14]+6640.783086353601*fIn[12]*epsyy[13]+7424.621202458751*epsyy[3]*fIn[12]+2969.8484809835013*epsyy[10]*fIn[11]+2969.8484809835013*fIn[10]*epsyy[11]+3320.3915431767996*epsyy[2]*fIn[10]+3320.3915431767996*fIn[2]*epsyy[10]+7424.62120245875*epsyy[5]*fIn[8]+3320.3915431767996*epsyy[4]*fIn[6]+3320.3915431767996*fIn[4]*epsyy[6])*fIn[18]+6640.783086353601*fIn[12]*fIn[14]*epsyy[18]+(664.0783086353599*epsyy[9]+474.3416490252571*epsyy[7]+742.462120245875*epsyy[0])*(fIn[17]*fIn[17])+(3320.3915431767996*fIn[8]*epsyy[17]+1328.15661727072*epsyy[13]*fIn[16]+1187.9393923934006*fIn[10]*epsyy[15]+3320.3915431767996*epsyy[11]*fIn[14]+(948.6832980505144*fIn[11]+1484.9242404917502*fIn[2])*epsyy[13]+2969.8484809835013*epsyy[10]*fIn[12]+1484.9242404917502*epsyy[3]*fIn[11]+1328.15661727072*epsyy[1]*fIn[10]+1484.9242404917502*fIn[6]*epsyy[7]+1328.15661727072*fIn[4]*epsyy[5])*fIn[17]+(3320.3915431767996*fIn[11]*fIn[14]+2969.8484809835013*fIn[10]*fIn[12])*epsyy[17]+(474.3416490252571*epsyy[9]+742.462120245875*epsyy[0])*(fIn[16]*fIn[16])+(3320.3915431767996*fIn[8]*epsyy[16]+1484.9242404917502*fIn[4]*epsyy[15]+2969.848480983501*epsyy[6]*fIn[14]+1328.1566172707196*epsyy[5]*fIn[10]+1484.9242404917502*fIn[2]*epsyy[9]+1328.1566172707196*epsyy[3]*fIn[6])*fIn[16]+2969.848480983501*fIn[6]*fIn[14]*epsyy[16]+1328.1566172707196*fIn[6]*fIn[10]*epsyy[15]+(3320.3915431767996*epsyy[9]+3320.3915431767996*epsyy[8]+3712.3106012293747*epsyy[0])*(fIn[14]*fIn[14])+(6640.783086353601*fIn[8]*epsyy[14]+7424.62120245875*epsyy[5]*fIn[12]+3320.391543176799*epsyy[4]*fIn[10]+3320.391543176799*fIn[4]*epsyy[10]+7424.621202458751*epsyy[3]*fIn[8]+3320.391543176799*epsyy[2]*fIn[6]+3320.391543176799*fIn[2]*epsyy[6])*fIn[14]+(1484.9242404917502*fIn[6]*fIn[11]+1328.1566172707196*fIn[4]*fIn[10])*epsyy[13]+(3320.3915431767996*epsyy[8]+3320.3915431767996*epsyy[7]+3712.3106012293747*epsyy[0])*(fIn[12]*fIn[12])+(6640.783086353601*fIn[8]*epsyy[12]+2969.848480983501*epsyy[4]*fIn[11]+2969.848480983501*fIn[4]*epsyy[11]+3320.391543176799*epsyy[6]*fIn[10]+3320.391543176799*fIn[6]*epsyy[10]+7424.621202458751*epsyy[1]*fIn[8]+3320.391543176799*epsyy[2]*fIn[4]+3320.391543176799*fIn[2]*epsyy[4])*fIn[12]+(474.3416490252571*epsyy[7]+742.462120245875*epsyy[0])*(fIn[11]*fIn[11])+(3320.3915431767996*fIn[8]*epsyy[11]+1328.1566172707196*epsyy[5]*fIn[10]+1484.9242404917502*fIn[2]*epsyy[7]+1328.1566172707196*epsyy[1]*fIn[4])*fIn[11]+(664.0783086353599*epsyy[9]+664.0783086353599*epsyy[7]+742.462120245875*epsyy[0])*(fIn[10]*fIn[10])+(3320.3915431767996*fIn[8]*epsyy[10]+1484.9242404917502*epsyy[1]*fIn[6]+1484.9242404917502*fIn[2]*epsyy[5]+1484.9242404917502*epsyy[3]*fIn[4])*fIn[10]+664.0783086353599*(fIn[6]*fIn[6])*epsyy[9]+(3320.3915431767996*epsyy[8]+3712.3106012293747*epsyy[0])*(fIn[8]*fIn[8])+(3320.3915431767996*epsyy[6]*fIn[6]+3320.3915431767996*epsyy[4]*fIn[4]+3320.3915431767996*epsyy[2]*fIn[2])*fIn[8]+664.0783086353599*(fIn[4]*fIn[4])*epsyy[7]+742.462120245875*epsyy[0]*(fIn[6]*fIn[6])+(1484.9242404917502*fIn[4]*epsyy[5]+1484.9242404917502*fIn[2]*epsyy[3])*fIn[6]+742.462120245875*epsyy[0]*(fIn[4]*fIn[4])+1484.9242404917502*epsyy[1]*fIn[2]*fIn[4]+742.462120245875*epsyy[0]*(fIn[2]*fIn[2]))*rdx11+((948.6832980505144*epsxy[19]+1484.9242404917502*epsxy[4])*(fIn[19]*fIn[19])+(2969.848480983501*fIn[4]*epsxy[19]+(3984.4698518121604*epsxy[18]+2969.848480983501*epsxy[5])*fIn[18]+(3984.4698518121604*epsxy[17]+2969.848480983501*epsxy[6])*fIn[17]+(948.6832980505144*epsxy[16]+1484.9242404917502*epsxy[2])*fIn[16]+(3320.3915431767996*fIn[11]+1484.9242404917502*fIn[2])*epsxy[16]+(948.6832980505144*epsxy[15]+1484.9242404917502*epsxy[1])*fIn[15]+(3320.3915431767996*fIn[12]+1484.9242404917502*fIn[1])*epsxy[15]+(2656.3132345414406*epsxy[14]+2969.8484809835013*epsxy[3])*fIn[14]+(2656.3132345414406*epsxy[13]+2969.8484809835013*epsxy[3])*fIn[13]+2656.3132345414406*epsxy[10]*fIn[10]+(3320.3915431767996*fIn[8]+3320.3915431767996*fIn[7])*epsxy[9]+1328.15661727072*epsxy[6]*fIn[6]+1328.15661727072*epsxy[5]*fIn[5])*fIn[19]+(2656.3132345414406*(fIn[18]*fIn[18])+2969.848480983501*fIn[5]*fIn[18]+2656.3132345414406*(fIn[17]*fIn[17])+2969.848480983501*fIn[6]*fIn[17]+3320.3915431767996*fIn[11]*fIn[16]+3320.3915431767996*fIn[12]*fIn[15]+6640.783086353601*fIn[13]*fIn[14]+1328.15661727072*(fIn[10]*fIn[10]))*epsxy[19]+2969.848480983501*epsxy[4]*(fIn[18]*fIn[18])+(4454.77272147525*fIn[4]*epsxy[18]+(6640.783086353601*epsxy[9]+6640.783086353601*epsxy[8]+6640.783086353601*epsxy[7]+7424.62120245875*epsxy[0])*fIn[17]+6640.783086353601*fIn[7]*epsxy[17]+1328.15661727072*epsxy[14]*fIn[16]+(2656.3132345414406*fIn[14]+6640.783086353601*fIn[13])*epsxy[16]+2969.8484809835013*epsxy[10]*fIn[15]+2969.8484809835013*fIn[10]*epsxy[15]+2969.8484809835013*epsxy[2]*fIn[14]+(6640.783086353601*fIn[11]+1484.9242404917502*fIn[2])*epsxy[14]+(6640.783086353601*epsxy[11]+7424.621202458751*epsxy[2])*fIn[13]+6640.783086353601*fIn[11]*epsxy[13]+5939.6969619670035*epsxy[10]*fIn[12]+4454.77272147525*fIn[10]*epsxy[12]+7424.621202458751*epsxy[3]*fIn[11]+3320.3915431767996*epsxy[1]*fIn[10]+3320.3915431767996*fIn[1]*epsxy[10]+2969.848480983501*epsxy[6]*fIn[8]+1484.9242404917502*fIn[6]*epsxy[8]+7424.62120245875*epsxy[6]*fIn[7]+3320.3915431767996*epsxy[4]*fIn[5]+3320.3915431767996*fIn[4]*epsxy[5])*fIn[18]+(6640.783086353601*fIn[8]*fIn[17]+6640.783086353601*fIn[11]*fIn[14]+4454.77272147525*fIn[10]*fIn[12])*epsxy[18]+2969.848480983501*epsxy[4]*(fIn[17]*fIn[17])+(4454.77272147525*fIn[4]*epsxy[17]+2969.8484809835013*epsxy[10]*fIn[16]+2969.8484809835013*fIn[10]*epsxy[16]+1328.15661727072*epsxy[13]*fIn[15]+(6640.783086353601*fIn[14]+2656.3132345414406*fIn[13])*epsxy[15]+(6640.783086353601*epsxy[12]+7424.621202458751*epsxy[1])*fIn[14]+6640.783086353601*fIn[12]*epsxy[14]+2969.8484809835013*epsxy[1]*fIn[13]+(6640.783086353601*fIn[12]+1484.9242404917502*fIn[1])*epsxy[13]+7424.621202458751*epsxy[3]*fIn[12]+5939.6969619670035*epsxy[10]*fIn[11]+4454.77272147525*fIn[10]*epsxy[11]+3320.3915431767996*epsxy[2]*fIn[10]+3320.3915431767996*fIn[2]*epsxy[10]+7424.62120245875*epsxy[5]*fIn[8]+2969.848480983501*epsxy[5]*fIn[7]+1484.9242404917502*fIn[5]*epsxy[7]+3320.3915431767996*epsxy[4]*fIn[6]+3320.3915431767996*fIn[4]*epsxy[6])*fIn[17]+(6640.783086353601*fIn[12]*fIn[13]+4454.77272147525*fIn[10]*fIn[11])*epsxy[17]+(1484.9242404917502*fIn[4]*epsxy[16]+(948.6832980505144*epsxy[9]+1484.9242404917502*epsxy[0])*fIn[15]+3320.3915431767996*fIn[7]*epsxy[15]+2969.848480983501*epsxy[5]*fIn[13]+1328.1566172707196*epsxy[6]*fIn[10]+1484.9242404917502*fIn[1]*epsxy[9]+1328.1566172707196*epsxy[3]*fIn[5])*fIn[16]+(3320.3915431767996*fIn[8]*fIn[15]+2969.848480983501*fIn[5]*fIn[14]+1328.1566172707196*fIn[6]*fIn[10])*epsxy[16]+(1484.9242404917502*fIn[4]*epsxy[15]+2969.848480983501*epsxy[6]*fIn[14]+1328.1566172707196*epsxy[5]*fIn[10]+1484.9242404917502*fIn[2]*epsxy[9]+1328.1566172707196*epsxy[3]*fIn[6])*fIn[15]+(2969.848480983501*fIn[6]*fIn[13]+1328.1566172707196*fIn[5]*fIn[10])*epsxy[15]+(2969.848480983501*fIn[4]*epsxy[14]+7424.62120245875*epsxy[4]*fIn[13]+2969.848480983501*epsxy[6]*fIn[12]+7424.62120245875*epsxy[5]*fIn[11]+(2969.8484809835013*epsxy[9]+2969.8484809835013*epsxy[8]+3320.391543176799*epsxy[0])*fIn[10]+7424.621202458751*fIn[7]*epsxy[10]+3320.391543176799*fIn[1]*epsxy[6]+3320.391543176799*epsxy[2]*fIn[5]+3320.391543176799*epsxy[3]*fIn[4])*fIn[14]+(1484.9242404917502*fIn[6]*fIn[12]+2969.8484809835013*fIn[8]*fIn[10])*epsxy[14]+(2969.848480983501*fIn[4]*epsxy[13]+7424.62120245875*epsxy[6]*fIn[12]+2969.848480983501*epsxy[5]*fIn[11]+(2969.8484809835013*epsxy[9]+2969.8484809835013*epsxy[7]+3320.391543176799*epsxy[0])*fIn[10]+7424.621202458751*fIn[8]*epsxy[10]+3320.391543176799*epsxy[1]*fIn[6]+3320.391543176799*fIn[2]*epsxy[5]+3320.391543176799*epsxy[3]*fIn[4])*fIn[13]+(1484.9242404917502*fIn[5]*fIn[11]+2969.8484809835013*fIn[7]*fIn[10])*epsxy[13]+2969.848480983501*epsxy[4]*(fIn[12]*fIn[12])+(4454.77272147525*fIn[4]*epsxy[12]+(6640.783086353601*epsxy[8]+6640.783086353601*epsxy[7]+7424.62120245875*epsxy[0])*fIn[11]+6640.783086353601*fIn[7]*epsxy[11]+3320.391543176799*epsxy[5]*fIn[10]+3320.391543176799*fIn[5]*epsxy[10]+2969.8484809835013*epsxy[2]*fIn[8]+1484.9242404917502*fIn[2]*epsxy[8]+7424.621202458751*epsxy[2]*fIn[7]+3320.391543176799*epsxy[1]*fIn[4]+3320.391543176799*fIn[1]*epsxy[4])*fIn[12]+6640.783086353601*fIn[8]*fIn[11]*epsxy[12]+2969.848480983501*epsxy[4]*(fIn[11]*fIn[11])+(4454.77272147525*fIn[4]*epsxy[11]+3320.391543176799*epsxy[6]*fIn[10]+3320.391543176799*fIn[6]*epsxy[10]+7424.621202458751*epsxy[1]*fIn[8]+2969.8484809835013*epsxy[1]*fIn[7]+1484.9242404917502*fIn[1]*epsxy[7]+3320.391543176799*epsxy[2]*fIn[4]+3320.391543176799*fIn[2]*epsxy[4])*fIn[11]+1484.9242404917502*epsxy[4]*(fIn[10]*fIn[10])+(2969.848480983501*fIn[4]*epsxy[10]+3320.3915431767996*epsxy[3]*fIn[8]+3320.3915431767996*epsxy[3]*fIn[7]+1484.9242404917502*epsxy[2]*fIn[6]+1484.9242404917502*fIn[2]*epsxy[6]+1484.9242404917502*epsxy[1]*fIn[5]+1484.9242404917502*fIn[1]*epsxy[5])*fIn[10]+1328.15661727072*fIn[5]*fIn[6]*epsxy[9]+(2969.848480983501*fIn[4]*epsxy[8]+7424.62120245875*epsxy[4]*fIn[7]+3320.3915431767996*fIn[5]*epsxy[6]+3320.3915431767996*epsxy[0]*fIn[4]+3320.3915431767996*fIn[1]*epsxy[2])*fIn[8]+(2969.848480983501*fIn[4]*epsxy[7]+3320.3915431767996*epsxy[5]*fIn[6]+3320.3915431767996*epsxy[0]*fIn[4]+3320.3915431767996*epsxy[1]*fIn[2])*fIn[7]+(1484.9242404917502*fIn[4]*epsxy[6]+1484.9242404917502*epsxy[0]*fIn[5]+1484.9242404917502*fIn[1]*epsxy[3])*fIn[6]+(1484.9242404917502*fIn[4]*epsxy[5]+1484.9242404917502*fIn[2]*epsxy[3])*fIn[5]+1484.9242404917502*epsxy[4]*(fIn[4]*fIn[4])+(1484.9242404917502*epsxy[2]*fIn[2]+1484.9242404917502*epsxy[1]*fIn[1])*fIn[4]+1484.9242404917502*epsxy[0]*fIn[1]*fIn[2])*rdx01+((474.3416490252571*epsxx[9]+664.0783086353599*epsxx[8]+742.462120245875*epsxx[0])*(fIn[19]*fIn[19])+(3320.3915431767996*fIn[7]*epsxx[19]+1187.9393923934006*epsxx[6]*fIn[18]+2656.3132345414406*fIn[17]*epsxx[18]+2969.848480983501*epsxx[5]*fIn[17]+(948.6832980505144*fIn[15]+1328.15661727072*fIn[12]+1484.9242404917502*fIn[1])*epsxx[16]+1484.9242404917502*epsxx[2]*fIn[15]+3320.3915431767996*fIn[11]*epsxx[15]+1187.9393923934006*fIn[10]*epsxx[14]+2969.8484809835013*epsxx[10]*fIn[13]+1328.15661727072*epsxx[3]*fIn[10]+1484.9242404917502*fIn[4]*epsxx[9]+1328.15661727072*fIn[5]*epsxx[6])*fIn[19]+(2656.3132345414406*fIn[17]*fIn[18]+2969.848480983501*fIn[5]*fIn[17]+3320.3915431767996*fIn[11]*fIn[15]+2969.8484809835013*fIn[10]*fIn[13])*epsxx[19]+(664.0783086353599*epsxx[9]+474.3416490252571*epsxx[8]+742.462120245875*epsxx[0])*(fIn[18]*fIn[18])+(3320.3915431767996*fIn[7]*epsxx[18]+2969.848480983501*epsxx[4]*fIn[17]+1187.9393923934006*fIn[10]*epsxx[16]+1328.15661727072*epsxx[14]*fIn[15]+(948.6832980505144*fIn[12]+1484.9242404917502*fIn[1])*epsxx[14]+3320.3915431767996*epsxx[12]*fIn[13]+1484.9242404917502*epsxx[3]*fIn[12]+2969.8484809835013*epsxx[10]*fIn[11]+1328.15661727072*epsxx[2]*fIn[10]+1484.9242404917502*fIn[5]*epsxx[8]+1328.15661727072*fIn[4]*epsxx[6])*fIn[18]+(2969.848480983501*fIn[4]*fIn[17]+3320.3915431767996*fIn[12]*fIn[13]+2969.8484809835013*fIn[10]*fIn[11])*epsxx[18]+(3320.3915431767996*epsxx[9]+3320.3915431767996*epsxx[8]+3320.3915431767996*epsxx[7]+3712.3106012293747*epsxx[0])*(fIn[17]*fIn[17])+(6640.783086353601*fIn[7]*epsxx[17]+6640.783086353601*fIn[13]*epsxx[16]+2969.8484809835013*epsxx[10]*fIn[15]+2969.8484809835013*fIn[10]*epsxx[15]+6640.783086353601*fIn[11]*epsxx[14]+(6640.783086353601*epsxx[11]+7424.621202458751*epsxx[2])*fIn[13]+6640.783086353601*fIn[11]*epsxx[13]+2969.8484809835013*epsxx[10]*fIn[12]+2969.8484809835013*fIn[10]*epsxx[12]+7424.621202458751*epsxx[3]*fIn[11]+3320.3915431767996*epsxx[1]*fIn[10]+3320.3915431767996*fIn[1]*epsxx[10]+7424.62120245875*epsxx[6]*fIn[7]+3320.3915431767996*epsxx[4]*fIn[5]+3320.3915431767996*fIn[4]*epsxx[5])*fIn[17]+6640.783086353601*fIn[11]*fIn[13]*epsxx[17]+(1484.9242404917502*fIn[4]*fIn[15]+1328.1566172707196*fIn[5]*fIn[10])*epsxx[16]+(474.3416490252571*epsxx[9]+742.462120245875*epsxx[0])*(fIn[15]*fIn[15])+(3320.3915431767996*fIn[7]*epsxx[15]+2969.848480983501*epsxx[5]*fIn[13]+1328.1566172707196*epsxx[6]*fIn[10]+1484.9242404917502*fIn[1]*epsxx[9]+1328.1566172707196*epsxx[3]*fIn[5])*fIn[15]+2969.848480983501*fIn[5]*fIn[13]*epsxx[15]+(1484.9242404917502*fIn[5]*fIn[12]+1328.1566172707196*fIn[4]*fIn[10])*epsxx[14]+(3320.3915431767996*epsxx[9]+3320.3915431767996*epsxx[7]+3712.3106012293747*epsxx[0])*(fIn[13]*fIn[13])+(6640.783086353601*fIn[7]*epsxx[13]+7424.62120245875*epsxx[6]*fIn[11]+3320.391543176799*epsxx[4]*fIn[10]+3320.391543176799*fIn[4]*epsxx[10]+7424.621202458751*epsxx[3]*fIn[7]+3320.391543176799*epsxx[1]*fIn[5]+3320.391543176799*fIn[1]*epsxx[5])*fIn[13]+(474.3416490252571*epsxx[8]+742.462120245875*epsxx[0])*(fIn[12]*fIn[12])+(3320.3915431767996*fIn[7]*epsxx[12]+2969.848480983501*epsxx[4]*fIn[11]+1328.1566172707196*epsxx[6]*fIn[10]+1484.9242404917502*fIn[1]*epsxx[8]+1328.1566172707196*epsxx[2]*fIn[4])*fIn[12]+2969.848480983501*fIn[4]*fIn[11]*epsxx[12]+(3320.3915431767996*epsxx[8]+3320.3915431767996*epsxx[7]+3712.3106012293747*epsxx[0])*(fIn[11]*fIn[11])+(6640.783086353601*fIn[7]*epsxx[11]+3320.391543176799*epsxx[5]*fIn[10]+3320.391543176799*fIn[5]*epsxx[10]+7424.621202458751*epsxx[2]*fIn[7]+3320.391543176799*epsxx[1]*fIn[4]+3320.391543176799*fIn[1]*epsxx[4])*fIn[11]+(664.0783086353599*epsxx[9]+664.0783086353599*epsxx[8]+742.462120245875*epsxx[0])*(fIn[10]*fIn[10])+(3320.3915431767996*fIn[7]*epsxx[10]+1484.9242404917502*fIn[1]*epsxx[6]+1484.9242404917502*epsxx[2]*fIn[5]+1484.9242404917502*epsxx[3]*fIn[4])*fIn[10]+664.0783086353599*(fIn[5]*fIn[5])*epsxx[9]+664.0783086353599*(fIn[4]*fIn[4])*epsxx[8]+(3320.3915431767996*epsxx[7]+3712.3106012293747*epsxx[0])*(fIn[7]*fIn[7])+(3320.3915431767996*epsxx[5]*fIn[5]+3320.3915431767996*epsxx[4]*fIn[4]+3320.3915431767996*epsxx[1]*fIn[1])*fIn[7]+1484.9242404917502*fIn[4]*fIn[5]*epsxx[6]+742.462120245875*epsxx[0]*(fIn[5]*fIn[5])+1484.9242404917502*fIn[1]*epsxx[3]*fIn[5]+742.462120245875*epsxx[0]*(fIn[4]*fIn[4])+1484.9242404917502*fIn[1]*epsxx[2]*fIn[4]+742.462120245875*epsxx[0]*(fIn[1]*fIn[1]))*rdx00)*vol; 
}

