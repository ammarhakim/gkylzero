#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2.0; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[8]; 
  out[0] += 1.414213562373095*F_0[0]*mass*volFact; 
  out[1] += 1.414213562373095*F_0[1]*mass*volFact; 
  out[2] += 1.414213562373095*F_0[4]*mass*volFact; 
  out[3] += mass*volFact*(1.414213562373095*F_0[0]*wvpar_sq+0.8164965809277261*F_0[2]*dvpar*wvpar+0.105409255338946*F_0[5]*dvpar_sq+0.1178511301977579*F_0[0]*dvpar_sq); 
  out[4] += mass*volFact*(1.414213562373095*F_0[1]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[7]*dvpar_sq+0.1178511301977579*F_0[1]*dvpar_sq); 
  out[5] += mass*volFact*(1.414213562373095*F_0[4]*wvpar_sq+0.816496580927726*F_0[6]*dvpar*wvpar+0.1178511301977579*F_0[4]*dvpar_sq); 
  out[6] += 1.414213562373095*G_1[0]*mass*volFact; 
  out[7] += 1.414213562373095*G_1[1]*mass*volFact; 
  out[8] += 1.414213562373095*G_1[4]*mass*volFact; 
} 
