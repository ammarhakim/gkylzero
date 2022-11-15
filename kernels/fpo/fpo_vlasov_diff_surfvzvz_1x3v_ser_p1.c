#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* w, const double* dx,
  const double* D[], const double* f[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: 
  // f: 
  // out: Incremented output

  const double Jvzvz = 4/dx[3]/dx[3];

  const double* Dlll = D[0] + 80;
  const double* flll = f[0];
  const double* Dllc = D[1] + 80;
  const double* fllc = f[1];
  const double* Dllu = D[2] + 80;
  const double* fllu = f[2];
  const double* Dlcl = D[3] + 80;
  const double* flcl = f[3];
  const double* Dlcc = D[4] + 80;
  const double* flcc = f[4];
  const double* Dlcu = D[5] + 80;
  const double* flcu = f[5];
  const double* Dlul = D[6] + 80;
  const double* flul = f[6];
  const double* Dluc = D[7] + 80;
  const double* fluc = f[7];
  const double* Dluu = D[8] + 80;
  const double* fluu = f[8];
  const double* Dcll = D[9] + 80;
  const double* fcll = f[9];
  const double* Dclc = D[10] + 80;
  const double* fclc = f[10];
  const double* Dclu = D[11] + 80;
  const double* fclu = f[11];
  const double* Dccl = D[12] + 80;
  const double* fccl = f[12];
  const double* Dccc = D[13] + 80;
  const double* fccc = f[13];
  const double* Dccu = D[14] + 80;
  const double* fccu = f[14];
  const double* Dcul = D[15] + 80;
  const double* fcul = f[15];
  const double* Dcuc = D[16] + 80;
  const double* fcuc = f[16];
  const double* Dcuu = D[17] + 80;
  const double* fcuu = f[17];
  const double* Dull = D[18] + 80;
  const double* full = f[18];
  const double* Dulc = D[19] + 80;
  const double* fulc = f[19];
  const double* Dulu = D[20] + 80;
  const double* fulu = f[20];
  const double* Ducl = D[21] + 80;
  const double* fucl = f[21];
  const double* Ducc = D[22] + 80;
  const double* fucc = f[22];
  const double* Ducu = D[23] + 80;
  const double* fucu = f[23];
  const double* Duul = D[24] + 80;
  const double* fuul = f[24];
  const double* Duuc = D[25] + 80;
  const double* fuuc = f[25];
  const double* Duuu = D[26] + 80;
  const double* fuuu = f[26];

  double D_proj1_l[8];
  D_proj1_l[0] = 0.408248290463863*Dccl[4]-0.408248290463863*Dccc[4]+0.3535533905932737*Dccl[0]+0.3535533905932737*Dccc[0];
  D_proj1_l[1] = 0.408248290463863*Dccl[8]-0.408248290463863*Dccc[8]+0.3535533905932737*Dccl[1]+0.3535533905932737*Dccc[1];
  D_proj1_l[2] = 0.408248290463863*Dccl[9]-0.408248290463863*Dccc[9]+0.3535533905932737*Dccl[2]+0.3535533905932737*Dccc[2];
  D_proj1_l[3] = 0.408248290463863*Dccl[10]-0.408248290463863*Dccc[10]+0.3535533905932737*Dccl[3]+0.3535533905932737*Dccc[3];
  D_proj1_l[4] = 0.408248290463863*Dccl[12]-0.408248290463863*Dccc[12]+0.3535533905932737*Dccl[5]+0.3535533905932737*Dccc[5];
  D_proj1_l[5] = 0.408248290463863*Dccl[13]-0.408248290463863*Dccc[13]+0.3535533905932737*Dccl[6]+0.3535533905932737*Dccc[6];
  D_proj1_l[6] = 0.408248290463863*Dccl[14]-0.408248290463863*Dccc[14]+0.3535533905932737*Dccl[7]+0.3535533905932737*Dccc[7];
  D_proj1_l[7] = 0.408248290463863*Dccl[15]-0.408248290463863*Dccc[15]+0.3535533905932737*Dccl[11]+0.3535533905932737*Dccc[11];

  double D_proj1_u[8];
  D_proj1_u[0] = (-0.408248290463863*Dccu[4])+0.408248290463863*Dccc[4]+0.3535533905932737*Dccu[0]+0.3535533905932737*Dccc[0];
  D_proj1_u[1] = (-0.408248290463863*Dccu[8])+0.408248290463863*Dccc[8]+0.3535533905932737*Dccu[1]+0.3535533905932737*Dccc[1];
  D_proj1_u[2] = (-0.408248290463863*Dccu[9])+0.408248290463863*Dccc[9]+0.3535533905932737*Dccu[2]+0.3535533905932737*Dccc[2];
  D_proj1_u[3] = (-0.408248290463863*Dccu[10])+0.408248290463863*Dccc[10]+0.3535533905932737*Dccu[3]+0.3535533905932737*Dccc[3];
  D_proj1_u[4] = (-0.408248290463863*Dccu[12])+0.408248290463863*Dccc[12]+0.3535533905932737*Dccu[5]+0.3535533905932737*Dccc[5];
  D_proj1_u[5] = (-0.408248290463863*Dccu[13])+0.408248290463863*Dccc[13]+0.3535533905932737*Dccu[6]+0.3535533905932737*Dccc[6];
  D_proj1_u[6] = (-0.408248290463863*Dccu[14])+0.408248290463863*Dccc[14]+0.3535533905932737*Dccu[7]+0.3535533905932737*Dccc[7];
  D_proj1_u[7] = (-0.408248290463863*Dccu[15])+0.408248290463863*Dccc[15]+0.3535533905932737*Dccu[11]+0.3535533905932737*Dccc[11];

  double df_proj1_l[8];
  df_proj1_l[0] = (-0.7654655446197428*fccl[4])-0.7654655446197428*fccc[4]-0.7954951288348656*fccl[0]+0.7954951288348656*fccc[0];
  df_proj1_l[1] = (-0.7654655446197428*fccl[8])-0.7654655446197428*fccc[8]-0.7954951288348656*fccl[1]+0.7954951288348656*fccc[1];
  df_proj1_l[2] = (-0.7654655446197428*fccl[9])-0.7654655446197428*fccc[9]-0.7954951288348656*fccl[2]+0.7954951288348656*fccc[2];
  df_proj1_l[3] = (-0.7654655446197428*fccl[10])-0.7654655446197428*fccc[10]-0.7954951288348656*fccl[3]+0.7954951288348656*fccc[3];
  df_proj1_l[4] = (-0.7654655446197428*fccl[12])-0.7654655446197428*fccc[12]-0.7954951288348656*fccl[5]+0.7954951288348656*fccc[5];
  df_proj1_l[5] = (-0.7654655446197428*fccl[13])-0.7654655446197428*fccc[13]-0.7954951288348656*fccl[6]+0.7954951288348656*fccc[6];
  df_proj1_l[6] = (-0.7654655446197428*fccl[14])-0.7654655446197428*fccc[14]-0.7954951288348656*fccl[7]+0.7954951288348656*fccc[7];
  df_proj1_l[7] = (-0.7654655446197428*fccl[15])-0.7654655446197428*fccc[15]-0.7954951288348656*fccl[11]+0.7954951288348656*fccc[11];

  double df_proj1_u[8];
  df_proj1_u[0] = (-0.7654655446197428*fccu[4])-0.7654655446197428*fccc[4]+0.7954951288348656*fccu[0]-0.7954951288348656*fccc[0];
  df_proj1_u[1] = (-0.7654655446197428*fccu[8])-0.7654655446197428*fccc[8]+0.7954951288348656*fccu[1]-0.7954951288348656*fccc[1];
  df_proj1_u[2] = (-0.7654655446197428*fccu[9])-0.7654655446197428*fccc[9]+0.7954951288348656*fccu[2]-0.7954951288348656*fccc[2];
  df_proj1_u[3] = (-0.7654655446197428*fccu[10])-0.7654655446197428*fccc[10]+0.7954951288348656*fccu[3]-0.7954951288348656*fccc[3];
  df_proj1_u[4] = (-0.7654655446197428*fccu[12])-0.7654655446197428*fccc[12]+0.7954951288348656*fccu[5]-0.7954951288348656*fccc[5];
  df_proj1_u[5] = (-0.7654655446197428*fccu[13])-0.7654655446197428*fccc[13]+0.7954951288348656*fccu[6]-0.7954951288348656*fccc[6];
  df_proj1_u[6] = (-0.7654655446197428*fccu[14])-0.7654655446197428*fccc[14]+0.7954951288348656*fccu[7]-0.7954951288348656*fccc[7];
  df_proj1_u[7] = (-0.7654655446197428*fccu[15])-0.7654655446197428*fccc[15]+0.7954951288348656*fccu[11]-0.7954951288348656*fccc[11];

  double f_proj2_l[8];
  f_proj2_l[0] = 0.408248290463863*fccl[4]-0.408248290463863*fccc[4]+0.3535533905932737*fccl[0]+0.3535533905932737*fccc[0];
  f_proj2_l[1] = 0.408248290463863*fccl[8]-0.408248290463863*fccc[8]+0.3535533905932737*fccl[1]+0.3535533905932737*fccc[1];
  f_proj2_l[2] = 0.408248290463863*fccl[9]-0.408248290463863*fccc[9]+0.3535533905932737*fccl[2]+0.3535533905932737*fccc[2];
  f_proj2_l[3] = 0.408248290463863*fccl[10]-0.408248290463863*fccc[10]+0.3535533905932737*fccl[3]+0.3535533905932737*fccc[3];
  f_proj2_l[4] = 0.408248290463863*fccl[12]-0.408248290463863*fccc[12]+0.3535533905932737*fccl[5]+0.3535533905932737*fccc[5];
  f_proj2_l[5] = 0.408248290463863*fccl[13]-0.408248290463863*fccc[13]+0.3535533905932737*fccl[6]+0.3535533905932737*fccc[6];
  f_proj2_l[6] = 0.408248290463863*fccl[14]-0.408248290463863*fccc[14]+0.3535533905932737*fccl[7]+0.3535533905932737*fccc[7];
  f_proj2_l[7] = 0.408248290463863*fccl[15]-0.408248290463863*fccc[15]+0.3535533905932737*fccl[11]+0.3535533905932737*fccc[11];

  double f_proj2_u[8];
  f_proj2_u[0] = (-0.408248290463863*fccu[4])+0.408248290463863*fccc[4]+0.3535533905932737*fccu[0]+0.3535533905932737*fccc[0];
  f_proj2_u[1] = (-0.408248290463863*fccu[8])+0.408248290463863*fccc[8]+0.3535533905932737*fccu[1]+0.3535533905932737*fccc[1];
  f_proj2_u[2] = (-0.408248290463863*fccu[9])+0.408248290463863*fccc[9]+0.3535533905932737*fccu[2]+0.3535533905932737*fccc[2];
  f_proj2_u[3] = (-0.408248290463863*fccu[10])+0.408248290463863*fccc[10]+0.3535533905932737*fccu[3]+0.3535533905932737*fccc[3];
  f_proj2_u[4] = (-0.408248290463863*fccu[12])+0.408248290463863*fccc[12]+0.3535533905932737*fccu[5]+0.3535533905932737*fccc[5];
  f_proj2_u[5] = (-0.408248290463863*fccu[13])+0.408248290463863*fccc[13]+0.3535533905932737*fccu[6]+0.3535533905932737*fccc[6];
  f_proj2_u[6] = (-0.408248290463863*fccu[14])+0.408248290463863*fccc[14]+0.3535533905932737*fccu[7]+0.3535533905932737*fccc[7];
  f_proj2_u[7] = (-0.408248290463863*fccu[15])+0.408248290463863*fccc[15]+0.3535533905932737*fccu[11]+0.3535533905932737*fccc[11];

  out[0] += Jvzvz*(0.125*D_proj1_u[7]*df_proj1_u[7]-0.125*D_proj1_l[7]*df_proj1_l[7]+0.125*D_proj1_u[6]*df_proj1_u[6]-0.125*D_proj1_l[6]*df_proj1_l[6]+0.125*D_proj1_u[5]*df_proj1_u[5]-0.125*D_proj1_l[5]*df_proj1_l[5]+0.125*D_proj1_u[4]*df_proj1_u[4]-0.125*D_proj1_l[4]*df_proj1_l[4]+0.125*D_proj1_u[3]*df_proj1_u[3]-0.125*D_proj1_l[3]*df_proj1_l[3]+0.125*D_proj1_u[2]*df_proj1_u[2]-0.125*D_proj1_l[2]*df_proj1_l[2]+0.125*D_proj1_u[1]*df_proj1_u[1]-0.125*D_proj1_l[1]*df_proj1_l[1]+0.125*D_proj1_u[0]*df_proj1_u[0]-0.125*D_proj1_l[0]*df_proj1_l[0]);
  out[1] += Jvzvz*(0.125*D_proj1_u[6]*df_proj1_u[7]-0.125*D_proj1_l[6]*df_proj1_l[7]+0.125*df_proj1_u[6]*D_proj1_u[7]-0.125*df_proj1_l[6]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[5]-0.125*D_proj1_l[3]*df_proj1_l[5]+0.125*df_proj1_u[3]*D_proj1_u[5]-0.125*df_proj1_l[3]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[4]-0.125*D_proj1_l[2]*df_proj1_l[4]+0.125*df_proj1_u[2]*D_proj1_u[4]-0.125*df_proj1_l[2]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[1]-0.125*D_proj1_l[0]*df_proj1_l[1]+0.125*df_proj1_u[0]*D_proj1_u[1]-0.125*df_proj1_l[0]*D_proj1_l[1]);
  out[2] += Jvzvz*(0.125*D_proj1_u[5]*df_proj1_u[7]-0.125*D_proj1_l[5]*df_proj1_l[7]+0.125*df_proj1_u[5]*D_proj1_u[7]-0.125*df_proj1_l[5]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[6]-0.125*D_proj1_l[3]*df_proj1_l[6]+0.125*df_proj1_u[3]*D_proj1_u[6]-0.125*df_proj1_l[3]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[4]-0.125*D_proj1_l[1]*df_proj1_l[4]+0.125*df_proj1_u[1]*D_proj1_u[4]-0.125*df_proj1_l[1]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[2]-0.125*D_proj1_l[0]*df_proj1_l[2]+0.125*df_proj1_u[0]*D_proj1_u[2]-0.125*df_proj1_l[0]*D_proj1_l[2]);
  out[3] += Jvzvz*(0.125*D_proj1_u[4]*df_proj1_u[7]-0.125*D_proj1_l[4]*df_proj1_l[7]+0.125*df_proj1_u[4]*D_proj1_u[7]-0.125*df_proj1_l[4]*D_proj1_l[7]+0.125*D_proj1_u[2]*df_proj1_u[6]-0.125*D_proj1_l[2]*df_proj1_l[6]+0.125*df_proj1_u[2]*D_proj1_u[6]-0.125*df_proj1_l[2]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[5]-0.125*D_proj1_l[1]*df_proj1_l[5]+0.125*df_proj1_u[1]*D_proj1_u[5]-0.125*df_proj1_l[1]*D_proj1_l[5]+0.125*D_proj1_u[0]*df_proj1_u[3]-0.125*D_proj1_l[0]*df_proj1_l[3]+0.125*df_proj1_u[0]*D_proj1_u[3]-0.125*df_proj1_l[0]*D_proj1_l[3]);
  out[4] += Jvzvz*((-0.2165063509461096*D_proj1_u[7]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[7]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[7]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[7]*df_proj1_l[7]-0.2165063509461096*D_proj1_u[6]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[6]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[6]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[5]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[4]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[3]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[2]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[1]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[0]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[0]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[0]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[0]);
  out[5] += Jvzvz*(0.125*D_proj1_u[3]*df_proj1_u[7]-0.125*D_proj1_l[3]*df_proj1_l[7]+0.125*df_proj1_u[3]*D_proj1_u[7]-0.125*df_proj1_l[3]*D_proj1_l[7]+0.125*D_proj1_u[5]*df_proj1_u[6]-0.125*D_proj1_l[5]*df_proj1_l[6]+0.125*df_proj1_u[5]*D_proj1_u[6]-0.125*df_proj1_l[5]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[4]-0.125*D_proj1_l[0]*df_proj1_l[4]+0.125*df_proj1_u[0]*D_proj1_u[4]-0.125*df_proj1_l[0]*D_proj1_l[4]+0.125*D_proj1_u[1]*df_proj1_u[2]-0.125*D_proj1_l[1]*df_proj1_l[2]+0.125*df_proj1_u[1]*D_proj1_u[2]-0.125*df_proj1_l[1]*D_proj1_l[2]);
  out[6] += Jvzvz*(0.125*D_proj1_u[2]*df_proj1_u[7]-0.125*D_proj1_l[2]*df_proj1_l[7]+0.125*df_proj1_u[2]*D_proj1_u[7]-0.125*df_proj1_l[2]*D_proj1_l[7]+0.125*D_proj1_u[4]*df_proj1_u[6]-0.125*D_proj1_l[4]*df_proj1_l[6]+0.125*df_proj1_u[4]*D_proj1_u[6]-0.125*df_proj1_l[4]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[5]-0.125*D_proj1_l[0]*df_proj1_l[5]+0.125*df_proj1_u[0]*D_proj1_u[5]-0.125*df_proj1_l[0]*D_proj1_l[5]+0.125*D_proj1_u[1]*df_proj1_u[3]-0.125*D_proj1_l[1]*df_proj1_l[3]+0.125*df_proj1_u[1]*D_proj1_u[3]-0.125*df_proj1_l[1]*D_proj1_l[3]);
  out[7] += Jvzvz*(0.125*D_proj1_u[1]*df_proj1_u[7]-0.125*D_proj1_l[1]*df_proj1_l[7]+0.125*df_proj1_u[1]*D_proj1_u[7]-0.125*df_proj1_l[1]*D_proj1_l[7]+0.125*D_proj1_u[0]*df_proj1_u[6]-0.125*D_proj1_l[0]*df_proj1_l[6]+0.125*df_proj1_u[0]*D_proj1_u[6]-0.125*df_proj1_l[0]*D_proj1_l[6]+0.125*D_proj1_u[4]*df_proj1_u[5]-0.125*D_proj1_l[4]*df_proj1_l[5]+0.125*df_proj1_u[4]*D_proj1_u[5]-0.125*df_proj1_l[4]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[3]-0.125*D_proj1_l[2]*df_proj1_l[3]+0.125*df_proj1_u[2]*D_proj1_u[3]-0.125*df_proj1_l[2]*D_proj1_l[3]);
  out[8] += Jvzvz*((-0.2165063509461096*D_proj1_u[6]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[6]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[6]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[6]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[6]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[6]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[1]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[1]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[1]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[1]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[1]);
  out[9] += Jvzvz*((-0.2165063509461096*D_proj1_u[5]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[5]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[2]);
  out[10] += Jvzvz*((-0.2165063509461096*D_proj1_u[4]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[4]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[3]);
  out[11] += Jvzvz*(0.125*D_proj1_u[0]*df_proj1_u[7]-0.125*D_proj1_l[0]*df_proj1_l[7]+0.125*df_proj1_u[0]*D_proj1_u[7]-0.125*df_proj1_l[0]*D_proj1_l[7]+0.125*D_proj1_u[1]*df_proj1_u[6]-0.125*D_proj1_l[1]*df_proj1_l[6]+0.125*df_proj1_u[1]*D_proj1_u[6]-0.125*df_proj1_l[1]*D_proj1_l[6]+0.125*D_proj1_u[2]*df_proj1_u[5]-0.125*D_proj1_l[2]*df_proj1_l[5]+0.125*df_proj1_u[2]*D_proj1_u[5]-0.125*df_proj1_l[2]*D_proj1_l[5]+0.125*D_proj1_u[3]*df_proj1_u[4]-0.125*D_proj1_l[3]*df_proj1_l[4]+0.125*df_proj1_u[3]*D_proj1_u[4]-0.125*df_proj1_l[3]*D_proj1_l[4]);
  out[12] += Jvzvz*((-0.2165063509461096*D_proj1_u[3]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[3]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[2]);
  out[13] += Jvzvz*((-0.2165063509461096*D_proj1_u[2]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[2]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[3]);
  out[14] += Jvzvz*((-0.2165063509461096*D_proj1_u[1]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[1]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[3]);
  out[15] += Jvzvz*((-0.2165063509461096*D_proj1_u[0]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[0]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[4]);
}

