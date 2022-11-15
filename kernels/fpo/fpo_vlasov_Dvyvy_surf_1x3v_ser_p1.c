#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_Dvyvy_surf_1x3v_ser_p1(const double* w, const double* dx,
  const double* g[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // g: 
  // out: Incremented output

  const double Jvyvy = 4/dx[2]/dx[2];

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

  out[48] += Jvyvy*((-0.5412658773652741*gcuc[3])+0.5412658773652741*gclc[3]+0.5625*gcuc[0]+0.5625*gclc[0]-1.125*gccc[0]);
  out[49] += Jvyvy*((-0.5412658773652741*gcuc[6])+0.5412658773652741*gclc[6]+0.5625*gcuc[1]+0.5625*gclc[1]-1.125*gccc[1]);
  out[50] += Jvyvy*((-0.5412658773652741*gcuc[7])+0.5412658773652741*gclc[7]+0.5625*gcuc[2]+0.5625*gclc[2]-1.125*gccc[2]);
  out[51] += Jvyvy*((-0.4375*gcuc[3])-0.4375*gclc[3]-2.875*gccc[3]+0.5412658773652739*gcuc[0]-0.5412658773652739*gclc[0]);
  out[52] += Jvyvy*((-0.5412658773652741*gcuc[10])+0.5412658773652741*gclc[10]+0.5625*gcuc[4]+0.5625*gclc[4]-1.125*gccc[4]);
  out[53] += Jvyvy*((-0.5412658773652741*gcuc[11])+0.5412658773652741*gclc[11]+0.5625*gcuc[5]+0.5625*gclc[5]-1.125*gccc[5]);
  out[54] += Jvyvy*((-0.4375*gcuc[6])-0.4375*gclc[6]-2.875*gccc[6]+0.5412658773652739*gcuc[1]-0.5412658773652739*gclc[1]);
  out[55] += Jvyvy*((-0.4375*gcuc[7])-0.4375*gclc[7]-2.875*gccc[7]+0.5412658773652739*gcuc[2]-0.5412658773652739*gclc[2]);
  out[56] += Jvyvy*((-0.5412658773652741*gcuc[13])+0.5412658773652741*gclc[13]+0.5625*gcuc[8]+0.5625*gclc[8]-1.125*gccc[8]);
  out[57] += Jvyvy*((-0.5412658773652741*gcuc[14])+0.5412658773652741*gclc[14]+0.5625*gcuc[9]+0.5625*gclc[9]-1.125*gccc[9]);
  out[58] += Jvyvy*((-0.4375*gcuc[10])-0.4375*gclc[10]-2.875*gccc[10]+0.5412658773652739*gcuc[4]-0.5412658773652739*gclc[4]);
  out[59] += Jvyvy*((-0.4375*gcuc[11])-0.4375*gclc[11]-2.875*gccc[11]+0.5412658773652739*gcuc[5]-0.5412658773652739*gclc[5]);
  out[60] += Jvyvy*((-0.5412658773652741*gcuc[15])+0.5412658773652741*gclc[15]+0.5625*gcuc[12]+0.5625*gclc[12]-1.125*gccc[12]);
  out[61] += Jvyvy*((-0.4375*gcuc[13])-0.4375*gclc[13]-2.875*gccc[13]+0.5412658773652739*gcuc[8]-0.5412658773652739*gclc[8]);
  out[62] += Jvyvy*((-0.4375*gcuc[14])-0.4375*gclc[14]-2.875*gccc[14]+0.5412658773652739*gcuc[9]-0.5412658773652739*gclc[9]);
  out[63] += Jvyvy*((-0.4375*gcuc[15])-0.4375*gclc[15]-2.875*gccc[15]+0.5412658773652739*gcuc[12]-0.5412658773652739*gclc[12]);
}

