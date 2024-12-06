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

  double dfdx0Sq[4] = {0.};
  dfdx0Sq[0] = (6.0*(fIn[3]*fIn[3])+6.0*(fIn[1]*fIn[1]))/dxSq[0]; 
  dfdx0Sq[2] = (12.0*fIn[1]*fIn[3])/dxSq[0]; 

  double dfdx1Sq[4] = {0.};
  dfdx1Sq[0] = (6.0*(fIn[3]*fIn[3])+6.0*(fIn[2]*fIn[2]))/dxSq[1]; 
  dfdx1Sq[1] = (12.0*fIn[2]*fIn[3])/dxSq[1]; 

  const double volFac = vol;

  out[0] += (dfdx0Sq[2]*weight[2]+dfdx1Sq[1]*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  double dfdx0Sq[8] = {0.};
  dfdx0Sq[0] = (6.0*(fIn[7]*fIn[7])+30.0*(fIn[6]*fIn[6])+30.0*(fIn[4]*fIn[4])+6.0*(fIn[3]*fIn[3])+6.0*(fIn[1]*fIn[1]))/dxSq[0]; 
  dfdx0Sq[1] = (26.832815729997474*fIn[3]*fIn[6]+26.832815729997478*fIn[1]*fIn[4])/dxSq[0]; 
  dfdx0Sq[2] = (0.2*(53.66563145999495*fIn[3]*fIn[7]+300.00000000000006*fIn[4]*fIn[6]+60.0*fIn[1]*fIn[3]))/dxSq[0]; 
  dfdx0Sq[3] = (24.0*fIn[6]*fIn[7]+26.832815729997474*fIn[1]*fIn[6]+26.832815729997478*fIn[3]*fIn[4])/dxSq[0]; 
  dfdx0Sq[4] = (26.832815729997478*(fIn[6]*fIn[6])+26.832815729997478*(fIn[4]*fIn[4]))/dxSq[0]; 
  dfdx0Sq[5] = (0.02857142857142857*(134.1640786499874*(fIn[7]*fIn[7])+420.00000000000006*fIn[1]*fIn[7]+939.1485505499119*(fIn[6]*fIn[6])+187.82971010998233*(fIn[3]*fIn[3])))/dxSq[0]; 
  dfdx0Sq[6] = (53.665631459994955*fIn[4]*fIn[6])/dxSq[0]; 
  dfdx0Sq[7] = (26.832815729997478*fIn[4]*fIn[7]+24.0*fIn[3]*fIn[6])/dxSq[0]; 

  double dfdx1Sq[8] = {0.};
  dfdx1Sq[0] = (30.0*(fIn[7]*fIn[7])+6.0*(fIn[6]*fIn[6])+30.0*(fIn[5]*fIn[5])+6.0*(fIn[3]*fIn[3])+6.0*(fIn[2]*fIn[2]))/dxSq[1]; 
  dfdx1Sq[1] = (0.2*(300.00000000000006*fIn[5]*fIn[7]+53.66563145999495*fIn[3]*fIn[6]+60.0*fIn[2]*fIn[3]))/dxSq[1]; 
  dfdx1Sq[2] = (26.832815729997474*fIn[3]*fIn[7]+26.832815729997478*fIn[2]*fIn[5])/dxSq[1]; 
  dfdx1Sq[3] = ((24.0*fIn[6]+26.832815729997474*fIn[2])*fIn[7]+26.832815729997478*fIn[3]*fIn[5])/dxSq[1]; 
  dfdx1Sq[4] = (0.02857142857142857*(939.1485505499119*(fIn[7]*fIn[7])+134.1640786499874*(fIn[6]*fIn[6])+420.00000000000006*fIn[2]*fIn[6]+187.82971010998233*(fIn[3]*fIn[3])))/dxSq[1]; 
  dfdx1Sq[5] = (26.832815729997478*(fIn[7]*fIn[7])+26.832815729997478*(fIn[5]*fIn[5]))/dxSq[1]; 
  dfdx1Sq[6] = (24.0*fIn[3]*fIn[7]+26.832815729997478*fIn[5]*fIn[6])/dxSq[1]; 
  dfdx1Sq[7] = (53.665631459994955*fIn[5]*fIn[7])/dxSq[1]; 

  const double volFac = vol;

  out[0] += ((dfdx1Sq[7]+dfdx0Sq[7])*weight[7]+(dfdx1Sq[6]+dfdx0Sq[6])*weight[6]+(dfdx1Sq[5]+dfdx0Sq[5])*weight[5]+(dfdx1Sq[4]+dfdx0Sq[4])*weight[4]+(dfdx1Sq[3]+dfdx0Sq[3])*weight[3]+(dfdx1Sq[2]+dfdx0Sq[2])*weight[2]+(dfdx1Sq[1]+dfdx0Sq[1])*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  double dfdx0Sq[8] = {0.};
  dfdx0Sq[0] = (4.242640687119286*(fIn[7]*fIn[7])+4.242640687119286*(fIn[5]*fIn[5])+4.242640687119286*(fIn[4]*fIn[4])+4.242640687119286*(fIn[1]*fIn[1]))/dxSq[0]; 
  dfdx0Sq[2] = (8.485281374238571*fIn[5]*fIn[7]+8.485281374238571*fIn[1]*fIn[4])/dxSq[0]; 
  dfdx0Sq[3] = (8.485281374238571*fIn[4]*fIn[7]+8.485281374238571*fIn[1]*fIn[5])/dxSq[0]; 
  dfdx0Sq[6] = (8.485281374238571*fIn[1]*fIn[7]+8.485281374238571*fIn[4]*fIn[5])/dxSq[0]; 

  double dfdx1Sq[8] = {0.};
  dfdx1Sq[0] = (4.242640687119286*(fIn[7]*fIn[7])+4.242640687119286*(fIn[6]*fIn[6])+4.242640687119286*(fIn[4]*fIn[4])+4.242640687119286*(fIn[2]*fIn[2]))/dxSq[1]; 
  dfdx1Sq[1] = (8.485281374238571*fIn[6]*fIn[7]+8.485281374238571*fIn[2]*fIn[4])/dxSq[1]; 
  dfdx1Sq[3] = (8.485281374238571*fIn[4]*fIn[7]+8.485281374238571*fIn[2]*fIn[6])/dxSq[1]; 
  dfdx1Sq[5] = (8.485281374238571*fIn[2]*fIn[7]+8.485281374238571*fIn[4]*fIn[6])/dxSq[1]; 

  const double volFac = vol;

  out[0] += (dfdx0Sq[6]*weight[6]+dfdx1Sq[5]*weight[5]+(dfdx1Sq[3]+dfdx0Sq[3])*weight[3]+dfdx0Sq[2]*weight[2]+dfdx1Sq[1]*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*volFac;
}

GKYL_CU_DH void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p2(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  double dfdx0Sq[20] = {0.};
  dfdx0Sq[0] = (4.242640687119286*(fIn[19]*fIn[19])+4.242640687119286*(fIn[18]*fIn[18])+21.213203435596427*(fIn[17]*fIn[17])+4.242640687119286*(fIn[15]*fIn[15])+21.213203435596427*(fIn[13]*fIn[13])+4.242640687119286*(fIn[12]*fIn[12])+21.213203435596427*(fIn[11]*fIn[11])+4.242640687119286*(fIn[10]*fIn[10])+21.213203435596427*(fIn[7]*fIn[7])+4.242640687119286*(fIn[5]*fIn[5])+4.242640687119286*(fIn[4]*fIn[4])+4.242640687119286*(fIn[1]*fIn[1]))/dxSq[0]; 
  dfdx0Sq[1] = (18.97366596101028*fIn[10]*fIn[17]+18.97366596101028*fIn[5]*fIn[13]+18.97366596101028*fIn[4]*fIn[11]+18.97366596101028*fIn[1]*fIn[7])/dxSq[0]; 
  dfdx0Sq[2] = (0.2*(42.42640687119286*fIn[15]*fIn[19]+37.94733192202057*fIn[10]*fIn[18]+212.13203435596432*fIn[13]*fIn[17]+37.94733192202057*fIn[4]*fIn[12]+212.13203435596432*fIn[7]*fIn[11]+42.42640687119286*fIn[5]*fIn[10]+42.42640687119286*fIn[1]*fIn[4]))/dxSq[0]; 
  dfdx0Sq[3] = (0.2*(37.94733192202057*fIn[10]*fIn[19]+42.42640687119286*fIn[12]*fIn[18]+212.13203435596432*fIn[11]*fIn[17]+37.94733192202057*fIn[5]*fIn[15]+212.13203435596432*fIn[7]*fIn[13]+42.42640687119286*fIn[4]*fIn[10]+42.42640687119286*fIn[1]*fIn[5]))/dxSq[0]; 
  dfdx0Sq[4] = (16.970562748477146*fIn[17]*fIn[18]+18.97366596101028*fIn[5]*fIn[17]+18.97366596101028*fIn[10]*fIn[13]+16.970562748477146*fIn[11]*fIn[12]+18.97366596101028*fIn[1]*fIn[11]+18.97366596101028*fIn[4]*fIn[7])/dxSq[0]; 
  dfdx0Sq[5] = (16.970562748477146*fIn[17]*fIn[19]+18.97366596101028*fIn[4]*fIn[17]+16.970562748477146*fIn[13]*fIn[15]+18.97366596101028*fIn[1]*fIn[13]+18.97366596101028*fIn[10]*fIn[11]+18.97366596101028*fIn[5]*fIn[7])/dxSq[0]; 
  dfdx0Sq[6] = (0.2*((33.9411254969543*fIn[18]+37.94733192202057*fIn[5])*fIn[19]+37.94733192202057*fIn[4]*fIn[18]+212.1320343559643*fIn[7]*fIn[17]+37.94733192202057*fIn[10]*fIn[15]+212.1320343559643*fIn[11]*fIn[13]+37.94733192202057*fIn[10]*fIn[12]+42.42640687119286*fIn[1]*fIn[10]+42.42640687119286*fIn[4]*fIn[5]))/dxSq[0]; 
  dfdx0Sq[7] = (0.7071067811865475*(26.832815729997478*(fIn[17]*fIn[17])+26.832815729997478*(fIn[13]*fIn[13])+26.832815729997478*(fIn[11]*fIn[11])+26.832815729997478*(fIn[7]*fIn[7])))/dxSq[0]; 
  dfdx0Sq[8] = (0.020203050891044214*(187.82971010998233*(fIn[19]*fIn[19])+134.1640786499874*(fIn[18]*fIn[18])+420.0*fIn[5]*fIn[18]+939.1485505499119*(fIn[17]*fIn[17])+134.1640786499874*(fIn[12]*fIn[12])+420.00000000000006*fIn[1]*fIn[12]+939.1485505499119*(fIn[11]*fIn[11])+187.82971010998233*(fIn[10]*fIn[10])+187.82971010998233*(fIn[4]*fIn[4])))/dxSq[0]; 
  dfdx0Sq[9] = (0.020203050891044214*(134.1640786499874*(fIn[19]*fIn[19])+420.0*fIn[4]*fIn[19]+187.82971010998233*(fIn[18]*fIn[18])+939.1485505499119*(fIn[17]*fIn[17])+134.1640786499874*(fIn[15]*fIn[15])+420.00000000000006*fIn[1]*fIn[15]+939.1485505499119*(fIn[13]*fIn[13])+187.82971010998233*(fIn[10]*fIn[10])+187.82971010998233*(fIn[5]*fIn[5])))/dxSq[0]; 
  dfdx0Sq[10] = (0.2*(84.85281374238573*fIn[13]*fIn[19]+84.85281374238573*fIn[11]*fIn[18]+(84.85281374238573*fIn[15]+84.85281374238573*fIn[12]+94.86832980505142*fIn[1])*fIn[17]+94.8683298050514*fIn[4]*fIn[13]+94.8683298050514*fIn[5]*fIn[11]+94.86832980505142*fIn[7]*fIn[10]))/dxSq[0]; 
  dfdx0Sq[11] = (0.7071067811865475*(53.665631459994955*fIn[13]*fIn[17]+53.665631459994955*fIn[7]*fIn[11]))/dxSq[0]; 
  dfdx0Sq[12] = (0.1414213562373095*(134.1640786499874*fIn[13]*fIn[18]+120.00000000000001*fIn[10]*fIn[17]+134.1640786499874*fIn[7]*fIn[12]+120.0*fIn[4]*fIn[11]))/dxSq[0]; 
  dfdx0Sq[13] = (0.7071067811865475*(53.665631459994955*fIn[11]*fIn[17]+53.665631459994955*fIn[7]*fIn[13]))/dxSq[0]; 
  dfdx0Sq[14] = (0.004040610178208843*(1680.0000000000002*fIn[10]*fIn[19]+(1878.2971010998237*fIn[15]+1341.6407864998741*fIn[12]+2100.0000000000005*fIn[1])*fIn[18]+9391.485505499119*fIn[11]*fIn[17]+2100.0*fIn[5]*fIn[12]+1878.2971010998233*fIn[4]*fIn[10]))/dxSq[0]; 
  dfdx0Sq[15] = (0.1414213562373095*(134.1640786499874*fIn[11]*fIn[19]+120.00000000000001*fIn[10]*fIn[17]+134.1640786499874*fIn[7]*fIn[15]+120.0*fIn[5]*fIn[13]))/dxSq[0]; 
  dfdx0Sq[16] = (0.004040610178208843*((1341.6407864998741*fIn[15]+1878.2971010998237*fIn[12]+2100.0000000000005*fIn[1])*fIn[19]+1680.0000000000002*fIn[10]*fIn[18]+9391.485505499119*fIn[13]*fIn[17]+2100.0*fIn[4]*fIn[15]+1878.2971010998233*fIn[5]*fIn[10]))/dxSq[0]; 
  dfdx0Sq[17] = (0.7071067811865475*(53.665631459994955*fIn[7]*fIn[17]+53.665631459994955*fIn[11]*fIn[13]))/dxSq[0]; 
  dfdx0Sq[18] = (0.1414213562373095*(107.33126291998991*fIn[17]*fIn[19]+134.1640786499874*fIn[7]*fIn[18]+120.0*fIn[4]*fIn[17]+134.1640786499874*fIn[12]*fIn[13]+120.00000000000001*fIn[10]*fIn[11]))/dxSq[0]; 
  dfdx0Sq[19] = (0.1414213562373095*(134.1640786499874*fIn[7]*fIn[19]+107.33126291998991*fIn[17]*fIn[18]+120.0*fIn[5]*fIn[17]+134.1640786499874*fIn[11]*fIn[15]+120.00000000000001*fIn[10]*fIn[13]))/dxSq[0]; 

  double dfdx1Sq[20] = {0.};
  dfdx1Sq[0] = (4.242640687119286*(fIn[19]*fIn[19])+21.213203435596427*(fIn[18]*fIn[18])+4.242640687119286*(fIn[17]*fIn[17])+4.242640687119286*(fIn[16]*fIn[16])+21.213203435596427*(fIn[14]*fIn[14])+21.213203435596427*(fIn[12]*fIn[12])+4.242640687119286*(fIn[11]*fIn[11])+4.242640687119286*(fIn[10]*fIn[10])+21.213203435596427*(fIn[8]*fIn[8])+4.242640687119286*(fIn[6]*fIn[6])+4.242640687119286*(fIn[4]*fIn[4])+4.242640687119286*(fIn[2]*fIn[2]))/dxSq[1]; 
  dfdx1Sq[1] = (0.2*(42.42640687119286*fIn[16]*fIn[19]+212.13203435596432*fIn[14]*fIn[18]+37.94733192202057*fIn[10]*fIn[17]+212.13203435596432*fIn[8]*fIn[12]+37.94733192202057*fIn[4]*fIn[11]+42.42640687119286*fIn[6]*fIn[10]+42.42640687119286*fIn[2]*fIn[4]))/dxSq[1]; 
  dfdx1Sq[2] = (18.97366596101028*fIn[10]*fIn[18]+18.97366596101028*fIn[6]*fIn[14]+18.97366596101028*fIn[4]*fIn[12]+18.97366596101028*fIn[2]*fIn[8])/dxSq[1]; 
  dfdx1Sq[3] = (0.2*(37.94733192202057*fIn[10]*fIn[19]+212.13203435596432*fIn[12]*fIn[18]+42.42640687119286*fIn[11]*fIn[17]+37.94733192202057*fIn[6]*fIn[16]+212.13203435596432*fIn[8]*fIn[14]+42.42640687119286*fIn[4]*fIn[10]+42.42640687119286*fIn[2]*fIn[6]))/dxSq[1]; 
  dfdx1Sq[4] = ((16.970562748477146*fIn[17]+18.97366596101028*fIn[6])*fIn[18]+18.97366596101028*fIn[10]*fIn[14]+(16.970562748477146*fIn[11]+18.97366596101028*fIn[2])*fIn[12]+18.97366596101028*fIn[4]*fIn[8])/dxSq[1]; 
  dfdx1Sq[5] = (0.2*((33.9411254969543*fIn[17]+37.94733192202057*fIn[6])*fIn[19]+212.1320343559643*fIn[8]*fIn[18]+37.94733192202057*fIn[4]*fIn[17]+37.94733192202057*fIn[10]*fIn[16]+212.1320343559643*fIn[12]*fIn[14]+37.94733192202057*fIn[10]*fIn[11]+42.42640687119286*fIn[2]*fIn[10]+42.42640687119286*fIn[4]*fIn[6]))/dxSq[1]; 
  dfdx1Sq[6] = (16.970562748477146*fIn[18]*fIn[19]+18.97366596101028*fIn[4]*fIn[18]+16.970562748477146*fIn[14]*fIn[16]+18.97366596101028*fIn[2]*fIn[14]+18.97366596101028*fIn[10]*fIn[12]+18.97366596101028*fIn[6]*fIn[8])/dxSq[1]; 
  dfdx1Sq[7] = (0.020203050891044214*(187.82971010998233*(fIn[19]*fIn[19])+939.1485505499119*(fIn[18]*fIn[18])+134.1640786499874*(fIn[17]*fIn[17])+420.0*fIn[6]*fIn[17]+939.1485505499119*(fIn[12]*fIn[12])+134.1640786499874*(fIn[11]*fIn[11])+420.00000000000006*fIn[2]*fIn[11]+187.82971010998233*(fIn[10]*fIn[10])+187.82971010998233*(fIn[4]*fIn[4])))/dxSq[1]; 
  dfdx1Sq[8] = (0.7071067811865475*(26.832815729997478*(fIn[18]*fIn[18])+26.832815729997478*(fIn[14]*fIn[14])+26.832815729997478*(fIn[12]*fIn[12])+26.832815729997478*(fIn[8]*fIn[8])))/dxSq[1]; 
  dfdx1Sq[9] = (0.020203050891044214*(134.1640786499874*(fIn[19]*fIn[19])+420.0*fIn[4]*fIn[19]+939.1485505499119*(fIn[18]*fIn[18])+187.82971010998233*(fIn[17]*fIn[17])+134.1640786499874*(fIn[16]*fIn[16])+420.00000000000006*fIn[2]*fIn[16]+939.1485505499119*(fIn[14]*fIn[14])+187.82971010998233*(fIn[10]*fIn[10])+187.82971010998233*(fIn[6]*fIn[6])))/dxSq[1]; 
  dfdx1Sq[10] = (0.2*(84.85281374238573*fIn[14]*fIn[19]+(84.85281374238573*fIn[16]+84.85281374238573*fIn[11]+94.86832980505142*fIn[2])*fIn[18]+84.85281374238573*fIn[12]*fIn[17]+94.8683298050514*fIn[4]*fIn[14]+94.8683298050514*fIn[6]*fIn[12]+94.86832980505142*fIn[8]*fIn[10]))/dxSq[1]; 
  dfdx1Sq[11] = (0.1414213562373095*(120.00000000000001*fIn[10]*fIn[18]+134.1640786499874*fIn[14]*fIn[17]+120.0*fIn[4]*fIn[12]+134.1640786499874*fIn[8]*fIn[11]))/dxSq[1]; 
  dfdx1Sq[12] = (0.7071067811865475*(53.665631459994955*fIn[14]*fIn[18]+53.665631459994955*fIn[8]*fIn[12]))/dxSq[1]; 
  dfdx1Sq[13] = (0.004040610178208843*(1680.0000000000002*fIn[10]*fIn[19]+9391.485505499119*fIn[12]*fIn[18]+(1878.2971010998237*fIn[16]+1341.6407864998741*fIn[11]+2100.0000000000005*fIn[2])*fIn[17]+2100.0*fIn[6]*fIn[11]+1878.2971010998233*fIn[4]*fIn[10]))/dxSq[1]; 
  dfdx1Sq[14] = (0.7071067811865475*(53.665631459994955*fIn[12]*fIn[18]+53.665631459994955*fIn[8]*fIn[14]))/dxSq[1]; 
  dfdx1Sq[15] = (0.004040610178208843*((1341.6407864998741*fIn[16]+1878.2971010998237*fIn[11]+2100.0000000000005*fIn[2])*fIn[19]+9391.485505499119*fIn[14]*fIn[18]+1680.0000000000002*fIn[10]*fIn[17]+2100.0*fIn[4]*fIn[16]+1878.2971010998233*fIn[6]*fIn[10]))/dxSq[1]; 
  dfdx1Sq[16] = (0.1414213562373095*(134.1640786499874*fIn[12]*fIn[19]+120.00000000000001*fIn[10]*fIn[18]+134.1640786499874*fIn[8]*fIn[16]+120.0*fIn[6]*fIn[14]))/dxSq[1]; 
  dfdx1Sq[17] = (0.1414213562373095*(107.33126291998991*fIn[18]*fIn[19]+120.0*fIn[4]*fIn[18]+134.1640786499874*fIn[8]*fIn[17]+134.1640786499874*fIn[11]*fIn[14]+120.00000000000001*fIn[10]*fIn[12]))/dxSq[1]; 
  dfdx1Sq[18] = (0.7071067811865475*(53.665631459994955*fIn[8]*fIn[18]+53.665631459994955*fIn[12]*fIn[14]))/dxSq[1]; 
  dfdx1Sq[19] = (0.1414213562373095*(134.1640786499874*fIn[8]*fIn[19]+(107.33126291998991*fIn[17]+120.0*fIn[6])*fIn[18]+134.1640786499874*fIn[12]*fIn[16]+120.00000000000001*fIn[10]*fIn[14]))/dxSq[1]; 

  const double volFac = vol;

  out[0] += ((dfdx1Sq[19]+dfdx0Sq[19])*weight[19]+(dfdx1Sq[18]+dfdx0Sq[18])*weight[18]+(dfdx1Sq[17]+dfdx0Sq[17])*weight[17]+(dfdx1Sq[16]+dfdx0Sq[16])*weight[16]+(dfdx1Sq[15]+dfdx0Sq[15])*weight[15]+(dfdx1Sq[14]+dfdx0Sq[14])*weight[14]+(dfdx1Sq[13]+dfdx0Sq[13])*weight[13]+(dfdx1Sq[12]+dfdx0Sq[12])*weight[12]+(dfdx1Sq[11]+dfdx0Sq[11])*weight[11]+(dfdx1Sq[10]+dfdx0Sq[10])*weight[10]+(dfdx1Sq[9]+dfdx0Sq[9])*weight[9]+(dfdx1Sq[8]+dfdx0Sq[8])*weight[8]+(dfdx1Sq[7]+dfdx0Sq[7])*weight[7]+(dfdx1Sq[6]+dfdx0Sq[6])*weight[6]+(dfdx1Sq[5]+dfdx0Sq[5])*weight[5]+(dfdx1Sq[4]+dfdx0Sq[4])*weight[4]+(dfdx1Sq[3]+dfdx0Sq[3])*weight[3]+(dfdx1Sq[2]+dfdx0Sq[2])*weight[2]+(dfdx1Sq[1]+dfdx0Sq[1])*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*volFac;
}

