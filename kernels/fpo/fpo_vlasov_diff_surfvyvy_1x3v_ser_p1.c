#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_diff_surfvyvy_1x3v_ser_p1(const double* w, const double* dx,
  const double* D[], const double* f[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: 
  // f: 
  // out: Incremented output

  const double Jvyvy = 4/dx[2]/dx[2];

  const double* Dlll = D[0] + 48;
  const double* flll = f[0];
  const double* Dllc = D[1] + 48;
  const double* fllc = f[1];
  const double* Dllu = D[2] + 48;
  const double* fllu = f[2];
  const double* Dlcl = D[3] + 48;
  const double* flcl = f[3];
  const double* Dlcc = D[4] + 48;
  const double* flcc = f[4];
  const double* Dlcu = D[5] + 48;
  const double* flcu = f[5];
  const double* Dlul = D[6] + 48;
  const double* flul = f[6];
  const double* Dluc = D[7] + 48;
  const double* fluc = f[7];
  const double* Dluu = D[8] + 48;
  const double* fluu = f[8];
  const double* Dcll = D[9] + 48;
  const double* fcll = f[9];
  const double* Dclc = D[10] + 48;
  const double* fclc = f[10];
  const double* Dclu = D[11] + 48;
  const double* fclu = f[11];
  const double* Dccl = D[12] + 48;
  const double* fccl = f[12];
  const double* Dccc = D[13] + 48;
  const double* fccc = f[13];
  const double* Dccu = D[14] + 48;
  const double* fccu = f[14];
  const double* Dcul = D[15] + 48;
  const double* fcul = f[15];
  const double* Dcuc = D[16] + 48;
  const double* fcuc = f[16];
  const double* Dcuu = D[17] + 48;
  const double* fcuu = f[17];
  const double* Dull = D[18] + 48;
  const double* full = f[18];
  const double* Dulc = D[19] + 48;
  const double* fulc = f[19];
  const double* Dulu = D[20] + 48;
  const double* fulu = f[20];
  const double* Ducl = D[21] + 48;
  const double* fucl = f[21];
  const double* Ducc = D[22] + 48;
  const double* fucc = f[22];
  const double* Ducu = D[23] + 48;
  const double* fucu = f[23];
  const double* Duul = D[24] + 48;
  const double* fuul = f[24];
  const double* Duuc = D[25] + 48;
  const double* fuuc = f[25];
  const double* Duuu = D[26] + 48;
  const double* fuuu = f[26];

  double D_proj1_l[8];
  D_proj1_l[0] = 0.408248290463863*Dclc[3]-0.408248290463863*Dccc[3]+0.3535533905932737*Dclc[0]+0.3535533905932737*Dccc[0];
  D_proj1_l[1] = 0.408248290463863*Dclc[6]-0.408248290463863*Dccc[6]+0.3535533905932737*Dclc[1]+0.3535533905932737*Dccc[1];
  D_proj1_l[2] = 0.408248290463863*Dclc[7]-0.408248290463863*Dccc[7]+0.3535533905932737*Dclc[2]+0.3535533905932737*Dccc[2];
  D_proj1_l[3] = 0.408248290463863*Dclc[10]-0.408248290463863*Dccc[10]+0.3535533905932737*Dclc[4]+0.3535533905932737*Dccc[4];
  D_proj1_l[4] = 0.408248290463863*Dclc[11]-0.408248290463863*Dccc[11]+0.3535533905932737*Dclc[5]+0.3535533905932737*Dccc[5];
  D_proj1_l[5] = 0.408248290463863*Dclc[13]-0.408248290463863*Dccc[13]+0.3535533905932737*Dclc[8]+0.3535533905932737*Dccc[8];
  D_proj1_l[6] = 0.408248290463863*Dclc[14]-0.408248290463863*Dccc[14]+0.3535533905932737*Dclc[9]+0.3535533905932737*Dccc[9];
  D_proj1_l[7] = 0.408248290463863*Dclc[15]-0.408248290463863*Dccc[15]+0.3535533905932737*Dclc[12]+0.3535533905932737*Dccc[12];

  double D_proj1_u[8];
  D_proj1_u[0] = (-0.408248290463863*Dcuc[3])+0.408248290463863*Dccc[3]+0.3535533905932737*Dcuc[0]+0.3535533905932737*Dccc[0];
  D_proj1_u[1] = (-0.408248290463863*Dcuc[6])+0.408248290463863*Dccc[6]+0.3535533905932737*Dcuc[1]+0.3535533905932737*Dccc[1];
  D_proj1_u[2] = (-0.408248290463863*Dcuc[7])+0.408248290463863*Dccc[7]+0.3535533905932737*Dcuc[2]+0.3535533905932737*Dccc[2];
  D_proj1_u[3] = (-0.408248290463863*Dcuc[10])+0.408248290463863*Dccc[10]+0.3535533905932737*Dcuc[4]+0.3535533905932737*Dccc[4];
  D_proj1_u[4] = (-0.408248290463863*Dcuc[11])+0.408248290463863*Dccc[11]+0.3535533905932737*Dcuc[5]+0.3535533905932737*Dccc[5];
  D_proj1_u[5] = (-0.408248290463863*Dcuc[13])+0.408248290463863*Dccc[13]+0.3535533905932737*Dcuc[8]+0.3535533905932737*Dccc[8];
  D_proj1_u[6] = (-0.408248290463863*Dcuc[14])+0.408248290463863*Dccc[14]+0.3535533905932737*Dcuc[9]+0.3535533905932737*Dccc[9];
  D_proj1_u[7] = (-0.408248290463863*Dcuc[15])+0.408248290463863*Dccc[15]+0.3535533905932737*Dcuc[12]+0.3535533905932737*Dccc[12];

  double df_proj1_l[8];
  df_proj1_l[0] = (-0.7654655446197428*fclc[3])-0.7654655446197428*fccc[3]-0.7954951288348656*fclc[0]+0.7954951288348656*fccc[0];
  df_proj1_l[1] = (-0.7654655446197428*fclc[6])-0.7654655446197428*fccc[6]-0.7954951288348656*fclc[1]+0.7954951288348656*fccc[1];
  df_proj1_l[2] = (-0.7654655446197428*fclc[7])-0.7654655446197428*fccc[7]-0.7954951288348656*fclc[2]+0.7954951288348656*fccc[2];
  df_proj1_l[3] = (-0.7654655446197428*fclc[10])-0.7654655446197428*fccc[10]-0.7954951288348656*fclc[4]+0.7954951288348656*fccc[4];
  df_proj1_l[4] = (-0.7654655446197428*fclc[11])-0.7654655446197428*fccc[11]-0.7954951288348656*fclc[5]+0.7954951288348656*fccc[5];
  df_proj1_l[5] = (-0.7654655446197428*fclc[13])-0.7654655446197428*fccc[13]-0.7954951288348656*fclc[8]+0.7954951288348656*fccc[8];
  df_proj1_l[6] = (-0.7654655446197428*fclc[14])-0.7654655446197428*fccc[14]-0.7954951288348656*fclc[9]+0.7954951288348656*fccc[9];
  df_proj1_l[7] = (-0.7654655446197428*fclc[15])-0.7654655446197428*fccc[15]-0.7954951288348656*fclc[12]+0.7954951288348656*fccc[12];

  double df_proj1_u[8];
  df_proj1_u[0] = (-0.7654655446197428*fcuc[3])-0.7654655446197428*fccc[3]+0.7954951288348656*fcuc[0]-0.7954951288348656*fccc[0];
  df_proj1_u[1] = (-0.7654655446197428*fcuc[6])-0.7654655446197428*fccc[6]+0.7954951288348656*fcuc[1]-0.7954951288348656*fccc[1];
  df_proj1_u[2] = (-0.7654655446197428*fcuc[7])-0.7654655446197428*fccc[7]+0.7954951288348656*fcuc[2]-0.7954951288348656*fccc[2];
  df_proj1_u[3] = (-0.7654655446197428*fcuc[10])-0.7654655446197428*fccc[10]+0.7954951288348656*fcuc[4]-0.7954951288348656*fccc[4];
  df_proj1_u[4] = (-0.7654655446197428*fcuc[11])-0.7654655446197428*fccc[11]+0.7954951288348656*fcuc[5]-0.7954951288348656*fccc[5];
  df_proj1_u[5] = (-0.7654655446197428*fcuc[13])-0.7654655446197428*fccc[13]+0.7954951288348656*fcuc[8]-0.7954951288348656*fccc[8];
  df_proj1_u[6] = (-0.7654655446197428*fcuc[14])-0.7654655446197428*fccc[14]+0.7954951288348656*fcuc[9]-0.7954951288348656*fccc[9];
  df_proj1_u[7] = (-0.7654655446197428*fcuc[15])-0.7654655446197428*fccc[15]+0.7954951288348656*fcuc[12]-0.7954951288348656*fccc[12];

  double f_proj2_l[8];
  f_proj2_l[0] = 0.408248290463863*fclc[3]-0.408248290463863*fccc[3]+0.3535533905932737*fclc[0]+0.3535533905932737*fccc[0];
  f_proj2_l[1] = 0.408248290463863*fclc[6]-0.408248290463863*fccc[6]+0.3535533905932737*fclc[1]+0.3535533905932737*fccc[1];
  f_proj2_l[2] = 0.408248290463863*fclc[7]-0.408248290463863*fccc[7]+0.3535533905932737*fclc[2]+0.3535533905932737*fccc[2];
  f_proj2_l[3] = 0.408248290463863*fclc[10]-0.408248290463863*fccc[10]+0.3535533905932737*fclc[4]+0.3535533905932737*fccc[4];
  f_proj2_l[4] = 0.408248290463863*fclc[11]-0.408248290463863*fccc[11]+0.3535533905932737*fclc[5]+0.3535533905932737*fccc[5];
  f_proj2_l[5] = 0.408248290463863*fclc[13]-0.408248290463863*fccc[13]+0.3535533905932737*fclc[8]+0.3535533905932737*fccc[8];
  f_proj2_l[6] = 0.408248290463863*fclc[14]-0.408248290463863*fccc[14]+0.3535533905932737*fclc[9]+0.3535533905932737*fccc[9];
  f_proj2_l[7] = 0.408248290463863*fclc[15]-0.408248290463863*fccc[15]+0.3535533905932737*fclc[12]+0.3535533905932737*fccc[12];

  double f_proj2_u[8];
  f_proj2_u[0] = (-0.408248290463863*fcuc[3])+0.408248290463863*fccc[3]+0.3535533905932737*fcuc[0]+0.3535533905932737*fccc[0];
  f_proj2_u[1] = (-0.408248290463863*fcuc[6])+0.408248290463863*fccc[6]+0.3535533905932737*fcuc[1]+0.3535533905932737*fccc[1];
  f_proj2_u[2] = (-0.408248290463863*fcuc[7])+0.408248290463863*fccc[7]+0.3535533905932737*fcuc[2]+0.3535533905932737*fccc[2];
  f_proj2_u[3] = (-0.408248290463863*fcuc[10])+0.408248290463863*fccc[10]+0.3535533905932737*fcuc[4]+0.3535533905932737*fccc[4];
  f_proj2_u[4] = (-0.408248290463863*fcuc[11])+0.408248290463863*fccc[11]+0.3535533905932737*fcuc[5]+0.3535533905932737*fccc[5];
  f_proj2_u[5] = (-0.408248290463863*fcuc[13])+0.408248290463863*fccc[13]+0.3535533905932737*fcuc[8]+0.3535533905932737*fccc[8];
  f_proj2_u[6] = (-0.408248290463863*fcuc[14])+0.408248290463863*fccc[14]+0.3535533905932737*fcuc[9]+0.3535533905932737*fccc[9];
  f_proj2_u[7] = (-0.408248290463863*fcuc[15])+0.408248290463863*fccc[15]+0.3535533905932737*fcuc[12]+0.3535533905932737*fccc[12];

  out[0] += Jvyvy*(0.125*D_proj1_u[7]*df_proj1_u[7]-0.125*D_proj1_l[7]*df_proj1_l[7]+0.125*D_proj1_u[6]*df_proj1_u[6]-0.125*D_proj1_l[6]*df_proj1_l[6]+0.125*D_proj1_u[5]*df_proj1_u[5]-0.125*D_proj1_l[5]*df_proj1_l[5]+0.125*D_proj1_u[4]*df_proj1_u[4]-0.125*D_proj1_l[4]*df_proj1_l[4]+0.125*D_proj1_u[3]*df_proj1_u[3]-0.125*D_proj1_l[3]*df_proj1_l[3]+0.125*D_proj1_u[2]*df_proj1_u[2]-0.125*D_proj1_l[2]*df_proj1_l[2]+0.125*D_proj1_u[1]*df_proj1_u[1]-0.125*D_proj1_l[1]*df_proj1_l[1]+0.125*D_proj1_u[0]*df_proj1_u[0]-0.125*D_proj1_l[0]*df_proj1_l[0]);
  out[1] += Jvyvy*(0.125*D_proj1_u[6]*df_proj1_u[7]-0.125*D_proj1_l[6]*df_proj1_l[7]+0.125*df_proj1_u[6]*D_proj1_u[7]-0.125*df_proj1_l[6]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[5]-0.125*D_proj1_l[3]*df_proj1_l[5]+0.125*df_proj1_u[3]*D_proj1_u[5]-0.125*df_proj1_l[3]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[4]-0.125*D_proj1_l[2]*df_proj1_l[4]+0.125*df_proj1_u[2]*D_proj1_u[4]-0.125*df_proj1_l[2]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[1]-0.125*D_proj1_l[0]*df_proj1_l[1]+0.125*df_proj1_u[0]*D_proj1_u[1]-0.125*df_proj1_l[0]*D_proj1_l[1]);
  out[2] += Jvyvy*(0.125*D_proj1_u[5]*df_proj1_u[7]-0.125*D_proj1_l[5]*df_proj1_l[7]+0.125*df_proj1_u[5]*D_proj1_u[7]-0.125*df_proj1_l[5]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[6]-0.125*D_proj1_l[3]*df_proj1_l[6]+0.125*df_proj1_u[3]*D_proj1_u[6]-0.125*df_proj1_l[3]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[4]-0.125*D_proj1_l[1]*df_proj1_l[4]+0.125*df_proj1_u[1]*D_proj1_u[4]-0.125*df_proj1_l[1]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[2]-0.125*D_proj1_l[0]*df_proj1_l[2]+0.125*df_proj1_u[0]*D_proj1_u[2]-0.125*df_proj1_l[0]*D_proj1_l[2]);
  out[3] += Jvyvy*((-0.2165063509461096*D_proj1_u[7]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[7]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[7]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[7]*df_proj1_l[7]-0.2165063509461096*D_proj1_u[6]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[6]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[6]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[5]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[4]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[3]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[2]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[1]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[0]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[0]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[0]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[0]);
  out[4] += Jvyvy*(0.125*D_proj1_u[4]*df_proj1_u[7]-0.125*D_proj1_l[4]*df_proj1_l[7]+0.125*df_proj1_u[4]*D_proj1_u[7]-0.125*df_proj1_l[4]*D_proj1_l[7]+0.125*D_proj1_u[2]*df_proj1_u[6]-0.125*D_proj1_l[2]*df_proj1_l[6]+0.125*df_proj1_u[2]*D_proj1_u[6]-0.125*df_proj1_l[2]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[5]-0.125*D_proj1_l[1]*df_proj1_l[5]+0.125*df_proj1_u[1]*D_proj1_u[5]-0.125*df_proj1_l[1]*D_proj1_l[5]+0.125*D_proj1_u[0]*df_proj1_u[3]-0.125*D_proj1_l[0]*df_proj1_l[3]+0.125*df_proj1_u[0]*D_proj1_u[3]-0.125*df_proj1_l[0]*D_proj1_l[3]);
  out[5] += Jvyvy*(0.125*D_proj1_u[3]*df_proj1_u[7]-0.125*D_proj1_l[3]*df_proj1_l[7]+0.125*df_proj1_u[3]*D_proj1_u[7]-0.125*df_proj1_l[3]*D_proj1_l[7]+0.125*D_proj1_u[5]*df_proj1_u[6]-0.125*D_proj1_l[5]*df_proj1_l[6]+0.125*df_proj1_u[5]*D_proj1_u[6]-0.125*df_proj1_l[5]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[4]-0.125*D_proj1_l[0]*df_proj1_l[4]+0.125*df_proj1_u[0]*D_proj1_u[4]-0.125*df_proj1_l[0]*D_proj1_l[4]+0.125*D_proj1_u[1]*df_proj1_u[2]-0.125*D_proj1_l[1]*df_proj1_l[2]+0.125*df_proj1_u[1]*D_proj1_u[2]-0.125*df_proj1_l[1]*D_proj1_l[2]);
  out[6] += Jvyvy*((-0.2165063509461096*D_proj1_u[6]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[6]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[6]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[6]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[6]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[6]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[1]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[1]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[1]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[1]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[1]);
  out[7] += Jvyvy*((-0.2165063509461096*D_proj1_u[5]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[5]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[2]);
  out[8] += Jvyvy*(0.125*D_proj1_u[2]*df_proj1_u[7]-0.125*D_proj1_l[2]*df_proj1_l[7]+0.125*df_proj1_u[2]*D_proj1_u[7]-0.125*df_proj1_l[2]*D_proj1_l[7]+0.125*D_proj1_u[4]*df_proj1_u[6]-0.125*D_proj1_l[4]*df_proj1_l[6]+0.125*df_proj1_u[4]*D_proj1_u[6]-0.125*df_proj1_l[4]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[5]-0.125*D_proj1_l[0]*df_proj1_l[5]+0.125*df_proj1_u[0]*D_proj1_u[5]-0.125*df_proj1_l[0]*D_proj1_l[5]+0.125*D_proj1_u[1]*df_proj1_u[3]-0.125*D_proj1_l[1]*df_proj1_l[3]+0.125*df_proj1_u[1]*D_proj1_u[3]-0.125*df_proj1_l[1]*D_proj1_l[3]);
  out[9] += Jvyvy*(0.125*D_proj1_u[1]*df_proj1_u[7]-0.125*D_proj1_l[1]*df_proj1_l[7]+0.125*df_proj1_u[1]*D_proj1_u[7]-0.125*df_proj1_l[1]*D_proj1_l[7]+0.125*D_proj1_u[0]*df_proj1_u[6]-0.125*D_proj1_l[0]*df_proj1_l[6]+0.125*df_proj1_u[0]*D_proj1_u[6]-0.125*df_proj1_l[0]*D_proj1_l[6]+0.125*D_proj1_u[4]*df_proj1_u[5]-0.125*D_proj1_l[4]*df_proj1_l[5]+0.125*df_proj1_u[4]*D_proj1_u[5]-0.125*df_proj1_l[4]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[3]-0.125*D_proj1_l[2]*df_proj1_l[3]+0.125*df_proj1_u[2]*D_proj1_u[3]-0.125*df_proj1_l[2]*D_proj1_l[3]);
  out[10] += Jvyvy*((-0.2165063509461096*D_proj1_u[4]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[4]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[3]);
  out[11] += Jvyvy*((-0.2165063509461096*D_proj1_u[3]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[3]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[2]);
  out[12] += Jvyvy*(0.125*D_proj1_u[0]*df_proj1_u[7]-0.125*D_proj1_l[0]*df_proj1_l[7]+0.125*df_proj1_u[0]*D_proj1_u[7]-0.125*df_proj1_l[0]*D_proj1_l[7]+0.125*D_proj1_u[1]*df_proj1_u[6]-0.125*D_proj1_l[1]*df_proj1_l[6]+0.125*df_proj1_u[1]*D_proj1_u[6]-0.125*df_proj1_l[1]*D_proj1_l[6]+0.125*D_proj1_u[2]*df_proj1_u[5]-0.125*D_proj1_l[2]*df_proj1_l[5]+0.125*df_proj1_u[2]*D_proj1_u[5]-0.125*df_proj1_l[2]*D_proj1_l[5]+0.125*D_proj1_u[3]*df_proj1_u[4]-0.125*D_proj1_l[3]*df_proj1_l[4]+0.125*df_proj1_u[3]*D_proj1_u[4]-0.125*df_proj1_l[3]*D_proj1_l[4]);
  out[13] += Jvyvy*((-0.2165063509461096*D_proj1_u[2]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[2]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[3]);
  out[14] += Jvyvy*((-0.2165063509461096*D_proj1_u[1]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[1]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[3]);
  out[15] += Jvyvy*((-0.2165063509461096*D_proj1_u[0]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[0]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[4]);
}

