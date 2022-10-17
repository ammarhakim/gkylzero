#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_2x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2.0; 
  const double wvpar = w[2], dvpar = dxv[2]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  out[0] += 1.414213562373095*f[0]*mass*volFact; 
  out[1] += 1.414213562373095*f[1]*mass*volFact; 
  out[2] += 1.414213562373095*f[2]*mass*volFact; 
  out[3] += 1.414213562373095*f[4]*mass*volFact; 
  out[4] += 1.414213562373095*f[7]*mass*volFact; 
  out[5] += 1.414213562373095*f[8]*mass*volFact; 
  out[6] += 1.414213562373095*f[11]*mass*volFact; 
  out[7] += 1.414213562373095*f[12]*mass*volFact; 
  out[8] += mass*volFact*(1.414213562373095*f[0]*wvpar_sq+0.8164965809277261*f[3]*dvpar*wvpar+0.105409255338946*f[9]*dvpar_sq+0.1178511301977579*f[0]*dvpar_sq); 
  out[9] += mass*volFact*(1.414213562373095*f[1]*wvpar_sq+0.8164965809277261*f[5]*dvpar*wvpar+0.105409255338946*f[15]*dvpar_sq+0.1178511301977579*f[1]*dvpar_sq); 
  out[10] += mass*volFact*(1.414213562373095*f[2]*wvpar_sq+0.8164965809277261*f[6]*dvpar*wvpar+0.105409255338946*f[16]*dvpar_sq+0.1178511301977579*f[2]*dvpar_sq); 
  out[11] += mass*volFact*(1.414213562373095*f[4]*wvpar_sq+0.8164965809277261*f[10]*dvpar*wvpar+0.105409255338946*f[19]*dvpar_sq+0.1178511301977579*f[4]*dvpar_sq); 
  out[12] += mass*volFact*(1.414213562373095*f[7]*wvpar_sq+0.816496580927726*f[13]*dvpar*wvpar+0.1178511301977579*f[7]*dvpar_sq); 
  out[13] += mass*volFact*(1.414213562373095*f[8]*wvpar_sq+0.816496580927726*f[14]*dvpar*wvpar+0.1178511301977579*f[8]*dvpar_sq); 
  out[14] += mass*volFact*(1.414213562373095*f[11]*wvpar_sq+0.816496580927726*f[17]*dvpar*wvpar+0.1178511301977579*f[11]*dvpar_sq); 
  out[15] += mass*volFact*(1.414213562373095*f[12]*wvpar_sq+0.816496580927726*f[18]*dvpar*wvpar+0.1178511301977579*f[12]*dvpar_sq); 
  out[16] += mass*volFact*(0.6123724356957944*f[3]*dvpar*wvpar_sq+0.7071067811865475*f[0]*wvpar_cu+0.1581138830084189*f[9]*dvpar_sq*wvpar+0.1767766952966368*f[0]*dvpar_sq*wvpar+0.03061862178478971*f[3]*dvpar_cu); 
  out[17] += mass*volFact*(0.6123724356957944*f[5]*dvpar*wvpar_sq+0.7071067811865475*f[1]*wvpar_cu+0.1581138830084189*f[15]*dvpar_sq*wvpar+0.1767766952966368*f[1]*dvpar_sq*wvpar+0.03061862178478971*f[5]*dvpar_cu); 
  out[18] += mass*volFact*(0.6123724356957944*f[6]*dvpar*wvpar_sq+0.7071067811865475*f[2]*wvpar_cu+0.1581138830084189*f[16]*dvpar_sq*wvpar+0.1767766952966368*f[2]*dvpar_sq*wvpar+0.03061862178478971*f[6]*dvpar_cu); 
  out[19] += mass*volFact*(0.6123724356957944*f[10]*dvpar*wvpar_sq+0.7071067811865475*f[4]*wvpar_cu+0.1581138830084189*f[19]*dvpar_sq*wvpar+0.1767766952966368*f[4]*dvpar_sq*wvpar+0.03061862178478971*f[10]*dvpar_cu); 
  out[20] += mass*volFact*(0.6123724356957944*f[13]*dvpar*wvpar_sq+0.7071067811865475*f[7]*wvpar_cu+0.1767766952966368*f[7]*dvpar_sq*wvpar+0.03061862178478971*f[13]*dvpar_cu); 
  out[21] += mass*volFact*(0.6123724356957944*f[14]*dvpar*wvpar_sq+0.7071067811865475*f[8]*wvpar_cu+0.1767766952966368*f[8]*dvpar_sq*wvpar+0.03061862178478971*f[14]*dvpar_cu); 
  out[22] += mass*volFact*(0.6123724356957944*f[17]*dvpar*wvpar_sq+0.7071067811865475*f[11]*wvpar_cu+0.1767766952966368*f[11]*dvpar_sq*wvpar+0.03061862178478971*f[17]*dvpar_cu); 
  out[23] += mass*volFact*(0.6123724356957944*f[18]*dvpar*wvpar_sq+0.7071067811865475*f[12]*wvpar_cu+0.1767766952966368*f[12]*dvpar_sq*wvpar+0.03061862178478971*f[18]*dvpar_cu); 
} 
