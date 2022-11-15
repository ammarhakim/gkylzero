#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_Dvzvz_surf_1x3v_ser_p1(const double* w, const double* dx,
  const double* g[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dx[NDIM]: Cell spacing
  // g: 
  // out: Incremented output

  const double Jvzvz = 4/dx[3]/dx[3];

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

  out[80] += Jvzvz*((-0.5412658773652741*gccu[4])+0.5412658773652741*gccl[4]+0.5625*gccu[0]+0.5625*gccl[0]-1.125*gccc[0]);
  out[81] += Jvzvz*((-0.5412658773652741*gccu[8])+0.5412658773652741*gccl[8]+0.5625*gccu[1]+0.5625*gccl[1]-1.125*gccc[1]);
  out[82] += Jvzvz*((-0.5412658773652741*gccu[9])+0.5412658773652741*gccl[9]+0.5625*gccu[2]+0.5625*gccl[2]-1.125*gccc[2]);
  out[83] += Jvzvz*((-0.5412658773652741*gccu[10])+0.5412658773652741*gccl[10]+0.5625*gccu[3]+0.5625*gccl[3]-1.125*gccc[3]);
  out[84] += Jvzvz*((-0.4375*gccu[4])-0.4375*gccl[4]-2.875*gccc[4]+0.5412658773652739*gccu[0]-0.5412658773652739*gccl[0]);
  out[85] += Jvzvz*((-0.5412658773652741*gccu[12])+0.5412658773652741*gccl[12]+0.5625*gccu[5]+0.5625*gccl[5]-1.125*gccc[5]);
  out[86] += Jvzvz*((-0.5412658773652741*gccu[13])+0.5412658773652741*gccl[13]+0.5625*gccu[6]+0.5625*gccl[6]-1.125*gccc[6]);
  out[87] += Jvzvz*((-0.5412658773652741*gccu[14])+0.5412658773652741*gccl[14]+0.5625*gccu[7]+0.5625*gccl[7]-1.125*gccc[7]);
  out[88] += Jvzvz*((-0.4375*gccu[8])-0.4375*gccl[8]-2.875*gccc[8]+0.5412658773652739*gccu[1]-0.5412658773652739*gccl[1]);
  out[89] += Jvzvz*((-0.4375*gccu[9])-0.4375*gccl[9]-2.875*gccc[9]+0.5412658773652739*gccu[2]-0.5412658773652739*gccl[2]);
  out[90] += Jvzvz*((-0.4375*gccu[10])-0.4375*gccl[10]-2.875*gccc[10]+0.5412658773652739*gccu[3]-0.5412658773652739*gccl[3]);
  out[91] += Jvzvz*((-0.5412658773652741*gccu[15])+0.5412658773652741*gccl[15]+0.5625*gccu[11]+0.5625*gccl[11]-1.125*gccc[11]);
  out[92] += Jvzvz*((-0.4375*gccu[12])-0.4375*gccl[12]-2.875*gccc[12]+0.5412658773652739*gccu[5]-0.5412658773652739*gccl[5]);
  out[93] += Jvzvz*((-0.4375*gccu[13])-0.4375*gccl[13]-2.875*gccc[13]+0.5412658773652739*gccu[6]-0.5412658773652739*gccl[6]);
  out[94] += Jvzvz*((-0.4375*gccu[14])-0.4375*gccl[14]-2.875*gccc[14]+0.5412658773652739*gccu[7]-0.5412658773652739*gccl[7]);
  out[95] += Jvzvz*((-0.4375*gccu[15])-0.4375*gccl[15]-2.875*gccc[15]+0.5412658773652739*gccu[11]-0.5412658773652739*gccl[11]);
}

