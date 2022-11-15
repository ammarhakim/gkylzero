#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_Dvxvx_surf_1x3v_ser_p1(const double* w, const double* dx,
  const double* g[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // g: 
  // out: Incremented output

  const double Jvxvx = 4/dx[1]/dx[1];

  const double* glll = g[0];
  const double* gllc = g[1];
  const double* gllu = g[2];
  const double* glcl = g[3];
  const double* glcc = g[4];
  const double* glcu = g[5];
  const double* glul = g[6];
  const double* gluc = g[7];
  const double* gluu = g[8];
  const double* gcll = g[9];
  const double* gclc = g[10];
  const double* gclu = g[11];
  const double* gccl = g[12];
  const double* gccc = g[13];
  const double* gccu = g[14];
  const double* gcul = g[15];
  const double* gcuc = g[16];
  const double* gcuu = g[17];
  const double* gull = g[18];
  const double* gulc = g[19];
  const double* gulu = g[20];
  const double* gucl = g[21];
  const double* gucc = g[22];
  const double* gucu = g[23];
  const double* guul = g[24];
  const double* guuc = g[25];
  const double* guuu = g[26];

  out[0] += Jvxvx*((-0.5412658773652741*gucc[2])+0.5412658773652741*glcc[2]+0.5625*gucc[0]+0.5625*glcc[0]-1.125*gccc[0]);
  out[1] += Jvxvx*((-0.5412658773652741*gucc[5])+0.5412658773652741*glcc[5]+0.5625*gucc[1]+0.5625*glcc[1]-1.125*gccc[1]);
  out[2] += Jvxvx*((-0.4375*gucc[2])-0.4375*glcc[2]-2.875*gccc[2]+0.5412658773652739*gucc[0]-0.5412658773652739*glcc[0]);
  out[3] += Jvxvx*((-0.5412658773652741*gucc[7])+0.5412658773652741*glcc[7]+0.5625*gucc[3]+0.5625*glcc[3]-1.125*gccc[3]);
  out[4] += Jvxvx*((-0.5412658773652741*gucc[9])+0.5412658773652741*glcc[9]+0.5625*gucc[4]+0.5625*glcc[4]-1.125*gccc[4]);
  out[5] += Jvxvx*((-0.4375*gucc[5])-0.4375*glcc[5]-2.875*gccc[5]+0.5412658773652739*gucc[1]-0.5412658773652739*glcc[1]);
  out[6] += Jvxvx*((-0.5412658773652741*gucc[11])+0.5412658773652741*glcc[11]+0.5625*gucc[6]+0.5625*glcc[6]-1.125*gccc[6]);
  out[7] += Jvxvx*((-0.4375*gucc[7])-0.4375*glcc[7]-2.875*gccc[7]+0.5412658773652739*gucc[3]-0.5412658773652739*glcc[3]);
  out[8] += Jvxvx*((-0.5412658773652741*gucc[12])+0.5412658773652741*glcc[12]+0.5625*gucc[8]+0.5625*glcc[8]-1.125*gccc[8]);
  out[9] += Jvxvx*((-0.4375*gucc[9])-0.4375*glcc[9]-2.875*gccc[9]+0.5412658773652739*gucc[4]-0.5412658773652739*glcc[4]);
  out[10] += Jvxvx*((-0.5412658773652741*gucc[14])+0.5412658773652741*glcc[14]+0.5625*gucc[10]+0.5625*glcc[10]-1.125*gccc[10]);
  out[11] += Jvxvx*((-0.4375*gucc[11])-0.4375*glcc[11]-2.875*gccc[11]+0.5412658773652739*gucc[6]-0.5412658773652739*glcc[6]);
  out[12] += Jvxvx*((-0.4375*gucc[12])-0.4375*glcc[12]-2.875*gccc[12]+0.5412658773652739*gucc[8]-0.5412658773652739*glcc[8]);
  out[13] += Jvxvx*((-0.5412658773652741*gucc[15])+0.5412658773652741*glcc[15]+0.5625*gucc[13]+0.5625*glcc[13]-1.125*gccc[13]);
  out[14] += Jvxvx*((-0.4375*gucc[14])-0.4375*glcc[14]-2.875*gccc[14]+0.5412658773652739*gucc[10]-0.5412658773652739*glcc[10]);
  out[15] += Jvxvx*((-0.4375*gucc[15])-0.4375*glcc[15]-2.875*gccc[15]+0.5412658773652739*gucc[13]-0.5412658773652739*glcc[13]);
}

