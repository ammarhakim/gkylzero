#include <gkyl_array_integrate_kernels.h>

GKYL_CU_DH void gkyl_array_integrate_op_none_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fIn[c*num_basis]*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_abs_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  for (unsigned c=0; c<num_comp; ++c)
    out[c] += fabs(fIn[c*num_basis])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_sq_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  for (unsigned c=0; c<num_comp; ++c) {
    for (unsigned b=0; b<num_basis; ++b)
      out[c] += pow(fIn[c*num_basis+b],2)*vol;
  }
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_1x_ser_p1_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += (fIn[1]*fIn[1])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_1x_ser_p2_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += (5.0*fIn[2]*fIn[2]+fIn[1]*fIn[1])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_2x_ser_p1_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += ((dxSq[1]+dxSq[0])*fIn[3]*fIn[3]+dxSq[0]*fIn[2]*fIn[2]+dxSq[1]*fIn[1]*fIn[1])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_2x_ser_p2_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += ((dxSq[1]+5*dxSq[0])*fIn[7]*fIn[7]+(5*dxSq[1]+dxSq[0])*fIn[6]*fIn[6]+5*dxSq[0]*fIn[5]*fIn[5]+5*dxSq[1]*fIn[4]*fIn[4]+(dxSq[1]+dxSq[0])*fIn[3]*fIn[3]+dxSq[0]*fIn[2]*fIn[2]+dxSq[1]*fIn[1]*fIn[1])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_grad_sq_3x_ser_p1_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += (((dxSq[1]+dxSq[0])*dxSq[2]+dxSq[0]*dxSq[1])*fIn[7]*fIn[7]+(dxSq[0]*dxSq[2]+dxSq[0]*dxSq[1])*fIn[6]*fIn[6]+(dxSq[1]*dxSq[2]+dxSq[0]*dxSq[1])*fIn[5]*fIn[5]+(dxSq[1]+dxSq[0])*dxSq[2]*fIn[4]*fIn[4]+dxSq[0]*dxSq[1]*fIn[3]*fIn[3]+dxSq[0]*dxSq[2]*fIn[2]*fIn[2]+dxSq[1]*fIn[1]*fIn[1]*dxSq[2])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_3x_ser_p1_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += ((dxSq[1]+dxSq[0])*fIn[7]*fIn[7]+dxSq[0]*fIn[6]*fIn[6]+dxSq[1]*fIn[5]*fIn[5]+(dxSq[1]+dxSq[0])*fIn[4]*fIn[4]+dxSq[0]*fIn[2]*fIn[2]+dxSq[1]*fIn[1]*fIn[1])*vol;
}

GKYL_CU_DH void gkyl_array_integrate_op_gradperp_sq_3x_ser_p2_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  out[0] += ((dxSq[1]+dxSq[0])*fIn[19]*fIn[19]+(dxSq[1]+5*dxSq[0])*fIn[18]*fIn[18]+(5*dxSq[1]+dxSq[0])*fIn[17]*fIn[17]+dxSq[0]*fIn[16]*fIn[16]+dxSq[1]*fIn[15]*fIn[15]+5*dxSq[0]*fIn[14]*fIn[14]+5*dxSq[1]*fIn[13]*fIn[13]+(dxSq[1]+5*dxSq[0])*fIn[12]*fIn[12]+(5*dxSq[1]+dxSq[0])*fIn[11]*fIn[11]+(dxSq[1]+dxSq[0])*fIn[10]*fIn[10]+5*dxSq[0]*fIn[8]*fIn[8]+5*dxSq[1]*fIn[7]*fIn[7]+dxSq[0]*fIn[6]*fIn[6]+dxSq[1]*fIn[5]*fIn[5]+(dxSq[1]+dxSq[0])*fIn[4]*fIn[4]+dxSq[0]*fIn[2]*fIn[2]+dxSq[1]*fIn[1]*fIn[1])*vol;
}

void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p1_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  double dfdx0Sq[4] = {0.};
  dfdx0Sq[0] = (6.0*(fIn[3]*fIn[3])+6.0*(fIn[1]*fIn[1]))/dxSq[0];
  dfdx0Sq[2] = (12.0*fIn[1]*fIn[3])/dxSq[0];

  double dfdx1Sq[4] = {0.};
  dfdx1Sq[0] = (6.0*(fIn[3]*fIn[3])+6.0*(fIn[2]*fIn[2]))/dxSq[1];
  dfdx1Sq[1] = (12.0*fIn[2]*fIn[3])/dxSq[1];

  out[0] += (dfdx0Sq[2]*weight[2]+dfdx1Sq[1]*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*vol;
}

void gkyl_array_integrate_op_eps_gradperp_sq_2x_ser_p2_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  double dfdx0Sq[8] = {0.};
  dfdx0Sq[0] = (6.0*(fIn[7]*fIn[7])+30.0*(fIn[6]*fIn[6])+30.0*(fIn[4]*fIn[4])+6.0*(fIn[3]*fIn[3])+6.0*(fIn[1]*fIn[1]))/dxSq[0];
  dfdx0Sq[1] = (26.83281572999747*fIn[3]*fIn[6]+26.83281572999748*fIn[1]*fIn[4])/dxSq[0];
  dfdx0Sq[2] = (0.2*(53.66563145999495*fIn[3]*fIn[7]+300.0000000000001*fIn[4]*fIn[6]+60.0*fIn[1]*fIn[3]))/dxSq[0];
  dfdx0Sq[3] = (24.0*fIn[6]*fIn[7]+26.83281572999747*fIn[1]*fIn[6]+26.83281572999748*fIn[3]*fIn[4])/dxSq[0];
  dfdx0Sq[4] = (26.83281572999748*(fIn[6]*fIn[6])+26.83281572999748*(fIn[4]*fIn[4]))/dxSq[0];
  dfdx0Sq[5] = (0.02857142857142857*(134.1640786499874*(fIn[7]*fIn[7])+420.0000000000001*fIn[1]*fIn[7]+939.1485505499119*(fIn[6]*fIn[6])+187.8297101099823*(fIn[3]*fIn[3])))/dxSq[0];
  dfdx0Sq[6] = (53.66563145999496*fIn[4]*fIn[6])/dxSq[0];
  dfdx0Sq[7] = (26.83281572999748*fIn[4]*fIn[7]+24.0*fIn[3]*fIn[6])/dxSq[0];

  double dfdx1Sq[8] = {0.};
  dfdx1Sq[0] = (30.0*(fIn[7]*fIn[7])+6.0*(fIn[6]*fIn[6])+30.0*(fIn[5]*fIn[5])+6.0*(fIn[3]*fIn[3])+6.0*(fIn[2]*fIn[2]))/dxSq[1];
  dfdx1Sq[1] = (0.2*(300.0000000000001*fIn[5]*fIn[7]+53.66563145999495*fIn[3]*fIn[6]+60.0*fIn[2]*fIn[3]))/dxSq[1];
  dfdx1Sq[2] = (26.83281572999747*fIn[3]*fIn[7]+26.83281572999748*fIn[2]*fIn[5])/dxSq[1];
  dfdx1Sq[3] = ((24.0*fIn[6]+26.83281572999747*fIn[2])*fIn[7]+26.83281572999748*fIn[3]*fIn[5])/dxSq[1];
  dfdx1Sq[4] = (0.02857142857142857*(939.1485505499119*(fIn[7]*fIn[7])+134.1640786499874*(fIn[6]*fIn[6])+420.0000000000001*fIn[2]*fIn[6]+187.8297101099823*(fIn[3]*fIn[3])))/dxSq[1];
  dfdx1Sq[5] = (26.83281572999748*(fIn[7]*fIn[7])+26.83281572999748*(fIn[5]*fIn[5]))/dxSq[1];
  dfdx1Sq[6] = (24.0*fIn[3]*fIn[7]+26.83281572999748*fIn[5]*fIn[6])/dxSq[1];
  dfdx1Sq[7] = (53.66563145999496*fIn[5]*fIn[7])/dxSq[1];

  out[0] += ((dfdx1Sq[7]+dfdx0Sq[7])*weight[7]+(dfdx1Sq[6]+dfdx0Sq[6])*weight[6]+(dfdx1Sq[5]+dfdx0Sq[5])*weight[5]+(dfdx1Sq[4]+dfdx0Sq[4])*weight[4]+(dfdx1Sq[3]+dfdx0Sq[3])*weight[3]+(dfdx1Sq[2]+dfdx0Sq[2])*weight[2]+(dfdx1Sq[1]+dfdx0Sq[1])*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*vol;
}

void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p1_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
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

  out[0] += (dfdx0Sq[6]*weight[6]+dfdx1Sq[5]*weight[5]+(dfdx1Sq[3]+dfdx0Sq[3])*weight[3]+dfdx0Sq[2]*weight[2]+dfdx1Sq[1]*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*vol;
}

void gkyl_array_integrate_op_eps_gradperp_sq_3x_ser_p2_ker(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out)
{
  double dfdx0Sq[20] = {0.};
  dfdx0Sq[0] = (4.242640687119286*(fIn[19]*fIn[19])+4.242640687119286*(fIn[18]*fIn[18])+21.21320343559643*(fIn[17]*fIn[17])+4.242640687119286*(fIn[15]*fIn[15])+21.21320343559643*(fIn[13]*fIn[13])+4.242640687119286*(fIn[12]*fIn[12])+21.21320343559643*(fIn[11]*fIn[11])+4.242640687119286*(fIn[10]*fIn[10])+21.21320343559643*(fIn[7]*fIn[7])+4.242640687119286*(fIn[5]*fIn[5])+4.242640687119286*(fIn[4]*fIn[4])+4.242640687119286*(fIn[1]*fIn[1]))/dxSq[0];
  dfdx0Sq[1] = (18.97366596101028*fIn[10]*fIn[17]+18.97366596101028*fIn[5]*fIn[13]+18.97366596101028*fIn[4]*fIn[11]+18.97366596101028*fIn[1]*fIn[7])/dxSq[0];
  dfdx0Sq[2] = (0.2*(42.42640687119286*fIn[15]*fIn[19]+37.94733192202057*fIn[10]*fIn[18]+212.1320343559643*fIn[13]*fIn[17]+37.94733192202057*fIn[4]*fIn[12]+212.1320343559643*fIn[7]*fIn[11]+42.42640687119286*fIn[5]*fIn[10]+42.42640687119286*fIn[1]*fIn[4]))/dxSq[0];
  dfdx0Sq[3] = (0.2*(37.94733192202057*fIn[10]*fIn[19]+42.42640687119286*fIn[12]*fIn[18]+212.1320343559643*fIn[11]*fIn[17]+37.94733192202057*fIn[5]*fIn[15]+212.1320343559643*fIn[7]*fIn[13]+42.42640687119286*fIn[4]*fIn[10]+42.42640687119286*fIn[1]*fIn[5]))/dxSq[0];
  dfdx0Sq[4] = (16.97056274847715*fIn[17]*fIn[18]+18.97366596101028*fIn[5]*fIn[17]+18.97366596101028*fIn[10]*fIn[13]+16.97056274847715*fIn[11]*fIn[12]+18.97366596101028*fIn[1]*fIn[11]+18.97366596101028*fIn[4]*fIn[7])/dxSq[0];
  dfdx0Sq[5] = (16.97056274847715*fIn[17]*fIn[19]+18.97366596101028*fIn[4]*fIn[17]+16.97056274847715*fIn[13]*fIn[15]+18.97366596101028*fIn[1]*fIn[13]+18.97366596101028*fIn[10]*fIn[11]+18.97366596101028*fIn[5]*fIn[7])/dxSq[0];
  dfdx0Sq[6] = (0.2*((33.9411254969543*fIn[18]+37.94733192202057*fIn[5])*fIn[19]+37.94733192202057*fIn[4]*fIn[18]+212.1320343559643*fIn[7]*fIn[17]+37.94733192202057*fIn[10]*fIn[15]+212.1320343559643*fIn[11]*fIn[13]+37.94733192202057*fIn[10]*fIn[12]+42.42640687119286*fIn[1]*fIn[10]+42.42640687119286*fIn[4]*fIn[5]))/dxSq[0];
  dfdx0Sq[7] = (0.7071067811865475*(26.83281572999748*(fIn[17]*fIn[17])+26.83281572999748*(fIn[13]*fIn[13])+26.83281572999748*(fIn[11]*fIn[11])+26.83281572999748*(fIn[7]*fIn[7])))/dxSq[0];
  dfdx0Sq[8] = (0.02020305089104421*(187.8297101099823*(fIn[19]*fIn[19])+134.1640786499874*(fIn[18]*fIn[18])+420.0*fIn[5]*fIn[18]+939.1485505499119*(fIn[17]*fIn[17])+134.1640786499874*(fIn[12]*fIn[12])+420.0000000000001*fIn[1]*fIn[12]+939.1485505499119*(fIn[11]*fIn[11])+187.8297101099823*(fIn[10]*fIn[10])+187.8297101099823*(fIn[4]*fIn[4])))/dxSq[0];
  dfdx0Sq[9] = (0.02020305089104421*(134.1640786499874*(fIn[19]*fIn[19])+420.0*fIn[4]*fIn[19]+187.8297101099823*(fIn[18]*fIn[18])+939.1485505499119*(fIn[17]*fIn[17])+134.1640786499874*(fIn[15]*fIn[15])+420.0000000000001*fIn[1]*fIn[15]+939.1485505499119*(fIn[13]*fIn[13])+187.8297101099823*(fIn[10]*fIn[10])+187.8297101099823*(fIn[5]*fIn[5])))/dxSq[0];
  dfdx0Sq[10] = (0.2*(84.85281374238573*fIn[13]*fIn[19]+84.85281374238573*fIn[11]*fIn[18]+(84.85281374238573*fIn[15]+84.85281374238573*fIn[12]+94.86832980505142*fIn[1])*fIn[17]+94.8683298050514*fIn[4]*fIn[13]+94.8683298050514*fIn[5]*fIn[11]+94.86832980505142*fIn[7]*fIn[10]))/dxSq[0];
  dfdx0Sq[11] = (0.7071067811865475*(53.66563145999496*fIn[13]*fIn[17]+53.66563145999496*fIn[7]*fIn[11]))/dxSq[0];
  dfdx0Sq[12] = (0.1414213562373095*(134.1640786499874*fIn[13]*fIn[18]+120.0*fIn[10]*fIn[17]+134.1640786499874*fIn[7]*fIn[12]+120.0*fIn[4]*fIn[11]))/dxSq[0];
  dfdx0Sq[13] = (0.7071067811865475*(53.66563145999496*fIn[11]*fIn[17]+53.66563145999496*fIn[7]*fIn[13]))/dxSq[0];
  dfdx0Sq[14] = (0.004040610178208843*(1680.0*fIn[10]*fIn[19]+(1878.297101099824*fIn[15]+1341.640786499874*fIn[12]+2100.0*fIn[1])*fIn[18]+9391.485505499119*fIn[11]*fIn[17]+2100.0*fIn[5]*fIn[12]+1878.297101099823*fIn[4]*fIn[10]))/dxSq[0];
  dfdx0Sq[15] = (0.1414213562373095*(134.1640786499874*fIn[11]*fIn[19]+120.0*fIn[10]*fIn[17]+134.1640786499874*fIn[7]*fIn[15]+120.0*fIn[5]*fIn[13]))/dxSq[0];
  dfdx0Sq[16] = (0.004040610178208843*((1341.640786499874*fIn[15]+1878.297101099824*fIn[12]+2100.0*fIn[1])*fIn[19]+1680.0*fIn[10]*fIn[18]+9391.485505499119*fIn[13]*fIn[17]+2100.0*fIn[4]*fIn[15]+1878.297101099823*fIn[5]*fIn[10]))/dxSq[0];
  dfdx0Sq[17] = (0.7071067811865475*(53.66563145999496*fIn[7]*fIn[17]+53.66563145999496*fIn[11]*fIn[13]))/dxSq[0];
  dfdx0Sq[18] = (0.1414213562373095*(107.3312629199899*fIn[17]*fIn[19]+134.1640786499874*fIn[7]*fIn[18]+120.0*fIn[4]*fIn[17]+134.1640786499874*fIn[12]*fIn[13]+120.0*fIn[10]*fIn[11]))/dxSq[0];
  dfdx0Sq[19] = (0.1414213562373095*(134.1640786499874*fIn[7]*fIn[19]+107.3312629199899*fIn[17]*fIn[18]+120.0*fIn[5]*fIn[17]+134.1640786499874*fIn[11]*fIn[15]+120.0*fIn[10]*fIn[13]))/dxSq[0];

  double dfdx1Sq[20] = {0.};
  dfdx1Sq[0] = (4.242640687119286*(fIn[19]*fIn[19])+21.21320343559643*(fIn[18]*fIn[18])+4.242640687119286*(fIn[17]*fIn[17])+4.242640687119286*(fIn[16]*fIn[16])+21.21320343559643*(fIn[14]*fIn[14])+21.21320343559643*(fIn[12]*fIn[12])+4.242640687119286*(fIn[11]*fIn[11])+4.242640687119286*(fIn[10]*fIn[10])+21.21320343559643*(fIn[8]*fIn[8])+4.242640687119286*(fIn[6]*fIn[6])+4.242640687119286*(fIn[4]*fIn[4])+4.242640687119286*(fIn[2]*fIn[2]))/dxSq[1];
  dfdx1Sq[1] = (0.2*(42.42640687119286*fIn[16]*fIn[19]+212.1320343559643*fIn[14]*fIn[18]+37.94733192202057*fIn[10]*fIn[17]+212.1320343559643*fIn[8]*fIn[12]+37.94733192202057*fIn[4]*fIn[11]+42.42640687119286*fIn[6]*fIn[10]+42.42640687119286*fIn[2]*fIn[4]))/dxSq[1];
  dfdx1Sq[2] = (18.97366596101028*fIn[10]*fIn[18]+18.97366596101028*fIn[6]*fIn[14]+18.97366596101028*fIn[4]*fIn[12]+18.97366596101028*fIn[2]*fIn[8])/dxSq[1];
  dfdx1Sq[3] = (0.2*(37.94733192202057*fIn[10]*fIn[19]+212.1320343559643*fIn[12]*fIn[18]+42.42640687119286*fIn[11]*fIn[17]+37.94733192202057*fIn[6]*fIn[16]+212.1320343559643*fIn[8]*fIn[14]+42.42640687119286*fIn[4]*fIn[10]+42.42640687119286*fIn[2]*fIn[6]))/dxSq[1];
  dfdx1Sq[4] = ((16.97056274847715*fIn[17]+18.97366596101028*fIn[6])*fIn[18]+18.97366596101028*fIn[10]*fIn[14]+(16.97056274847715*fIn[11]+18.97366596101028*fIn[2])*fIn[12]+18.97366596101028*fIn[4]*fIn[8])/dxSq[1];
  dfdx1Sq[5] = (0.2*((33.9411254969543*fIn[17]+37.94733192202057*fIn[6])*fIn[19]+212.1320343559643*fIn[8]*fIn[18]+37.94733192202057*fIn[4]*fIn[17]+37.94733192202057*fIn[10]*fIn[16]+212.1320343559643*fIn[12]*fIn[14]+37.94733192202057*fIn[10]*fIn[11]+42.42640687119286*fIn[2]*fIn[10]+42.42640687119286*fIn[4]*fIn[6]))/dxSq[1];
  dfdx1Sq[6] = (16.97056274847715*fIn[18]*fIn[19]+18.97366596101028*fIn[4]*fIn[18]+16.97056274847715*fIn[14]*fIn[16]+18.97366596101028*fIn[2]*fIn[14]+18.97366596101028*fIn[10]*fIn[12]+18.97366596101028*fIn[6]*fIn[8])/dxSq[1];
  dfdx1Sq[7] = (0.02020305089104421*(187.8297101099823*(fIn[19]*fIn[19])+939.1485505499119*(fIn[18]*fIn[18])+134.1640786499874*(fIn[17]*fIn[17])+420.0*fIn[6]*fIn[17]+939.1485505499119*(fIn[12]*fIn[12])+134.1640786499874*(fIn[11]*fIn[11])+420.0000000000001*fIn[2]*fIn[11]+187.8297101099823*(fIn[10]*fIn[10])+187.8297101099823*(fIn[4]*fIn[4])))/dxSq[1];
  dfdx1Sq[8] = (0.7071067811865475*(26.83281572999748*(fIn[18]*fIn[18])+26.83281572999748*(fIn[14]*fIn[14])+26.83281572999748*(fIn[12]*fIn[12])+26.83281572999748*(fIn[8]*fIn[8])))/dxSq[1];
  dfdx1Sq[9] = (0.02020305089104421*(134.1640786499874*(fIn[19]*fIn[19])+420.0*fIn[4]*fIn[19]+939.1485505499119*(fIn[18]*fIn[18])+187.8297101099823*(fIn[17]*fIn[17])+134.1640786499874*(fIn[16]*fIn[16])+420.0000000000001*fIn[2]*fIn[16]+939.1485505499119*(fIn[14]*fIn[14])+187.8297101099823*(fIn[10]*fIn[10])+187.8297101099823*(fIn[6]*fIn[6])))/dxSq[1];
  dfdx1Sq[10] = (0.2*(84.85281374238573*fIn[14]*fIn[19]+(84.85281374238573*fIn[16]+84.85281374238573*fIn[11]+94.86832980505142*fIn[2])*fIn[18]+84.85281374238573*fIn[12]*fIn[17]+94.8683298050514*fIn[4]*fIn[14]+94.8683298050514*fIn[6]*fIn[12]+94.86832980505142*fIn[8]*fIn[10]))/dxSq[1];
  dfdx1Sq[11] = (0.1414213562373095*(120.0*fIn[10]*fIn[18]+134.1640786499874*fIn[14]*fIn[17]+120.0*fIn[4]*fIn[12]+134.1640786499874*fIn[8]*fIn[11]))/dxSq[1];
  dfdx1Sq[12] = (0.7071067811865475*(53.66563145999496*fIn[14]*fIn[18]+53.66563145999496*fIn[8]*fIn[12]))/dxSq[1];
  dfdx1Sq[13] = (0.004040610178208843*(1680.0*fIn[10]*fIn[19]+9391.485505499119*fIn[12]*fIn[18]+(1878.297101099824*fIn[16]+1341.640786499874*fIn[11]+2100.0*fIn[2])*fIn[17]+2100.0*fIn[6]*fIn[11]+1878.297101099823*fIn[4]*fIn[10]))/dxSq[1];
  dfdx1Sq[14] = (0.7071067811865475*(53.66563145999496*fIn[12]*fIn[18]+53.66563145999496*fIn[8]*fIn[14]))/dxSq[1];
  dfdx1Sq[15] = (0.004040610178208843*((1341.640786499874*fIn[16]+1878.297101099824*fIn[11]+2100.0*fIn[2])*fIn[19]+9391.485505499119*fIn[14]*fIn[18]+1680.0*fIn[10]*fIn[17]+2100.0*fIn[4]*fIn[16]+1878.297101099823*fIn[6]*fIn[10]))/dxSq[1];
  dfdx1Sq[16] = (0.1414213562373095*(134.1640786499874*fIn[12]*fIn[19]+120.0*fIn[10]*fIn[18]+134.1640786499874*fIn[8]*fIn[16]+120.0*fIn[6]*fIn[14]))/dxSq[1];
  dfdx1Sq[17] = (0.1414213562373095*(107.3312629199899*fIn[18]*fIn[19]+120.0*fIn[4]*fIn[18]+134.1640786499874*fIn[8]*fIn[17]+134.1640786499874*fIn[11]*fIn[14]+120.0*fIn[10]*fIn[12]))/dxSq[1];
  dfdx1Sq[18] = (0.7071067811865475*(53.66563145999496*fIn[8]*fIn[18]+53.66563145999496*fIn[12]*fIn[14]))/dxSq[1];
  dfdx1Sq[19] = (0.1414213562373095*(134.1640786499874*fIn[8]*fIn[19]+(107.3312629199899*fIn[17]+120.0*fIn[6])*fIn[18]+134.1640786499874*fIn[12]*fIn[16]+120.0*fIn[10]*fIn[14]))/dxSq[1];

  out[0] += ((dfdx1Sq[19]+dfdx0Sq[19])*weight[19]+(dfdx1Sq[18]+dfdx0Sq[18])*weight[18]+(dfdx1Sq[17]+dfdx0Sq[17])*weight[17]+(dfdx1Sq[16]+dfdx0Sq[16])*weight[16]+(dfdx1Sq[15]+dfdx0Sq[15])*weight[15]+(dfdx1Sq[14]+dfdx0Sq[14])*weight[14]+(dfdx1Sq[13]+dfdx0Sq[13])*weight[13]+(dfdx1Sq[12]+dfdx0Sq[12])*weight[12]+(dfdx1Sq[11]+dfdx0Sq[11])*weight[11]+(dfdx1Sq[10]+dfdx0Sq[10])*weight[10]+(dfdx1Sq[9]+dfdx0Sq[9])*weight[9]+(dfdx1Sq[8]+dfdx0Sq[8])*weight[8]+(dfdx1Sq[7]+dfdx0Sq[7])*weight[7]+(dfdx1Sq[6]+dfdx0Sq[6])*weight[6]+(dfdx1Sq[5]+dfdx0Sq[5])*weight[5]+(dfdx1Sq[4]+dfdx0Sq[4])*weight[4]+(dfdx1Sq[3]+dfdx0Sq[3])*weight[3]+(dfdx1Sq[2]+dfdx0Sq[2])*weight[2]+(dfdx1Sq[1]+dfdx0Sq[1])*weight[1]+(dfdx1Sq[0]+dfdx0Sq[0])*weight[0])*vol;
}
