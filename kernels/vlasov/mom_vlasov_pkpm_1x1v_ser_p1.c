#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2.0; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  out[0] += 1.414213562373095*f[0]*mass*volFact; 
  out[1] += 1.414213562373095*f[1]*mass*volFact; 
  out[2] += mass*volFact*(1.414213562373095*f[0]*wvpar_sq+0.8164965809277261*f[2]*dvpar*wvpar+0.105409255338946*f[4]*dvpar_sq+0.1178511301977579*f[0]*dvpar_sq); 
  out[3] += mass*volFact*(1.414213562373095*f[1]*wvpar_sq+0.8164965809277261*f[3]*dvpar*wvpar+0.105409255338946*f[5]*dvpar_sq+0.1178511301977579*f[1]*dvpar_sq); 
  out[4] += mass*volFact*(0.6123724356957944*f[2]*dvpar*wvpar_sq+0.7071067811865475*f[0]*wvpar_cu+0.1581138830084189*f[4]*dvpar_sq*wvpar+0.1767766952966368*f[0]*dvpar_sq*wvpar+0.03061862178478971*f[2]*dvpar_cu); 
  out[5] += mass*volFact*(0.6123724356957944*f[3]*dvpar*wvpar_sq+0.7071067811865475*f[1]*wvpar_cu+0.1581138830084189*f[5]*dvpar_sq*wvpar+0.1767766952966368*f[1]*dvpar_sq*wvpar+0.03061862178478971*f[3]*dvpar_cu); 
} 
