#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_2x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2.0; 
  const double wvpar = w[2], dvpar = dxv[2]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[12]; 
  out[0] += 1.414213562373095*F_0[0]*mass*volFact; 
  out[1] += 1.414213562373095*F_0[1]*mass*volFact; 
  out[2] += 1.414213562373095*F_0[2]*mass*volFact; 
  out[3] += 1.414213562373095*F_0[4]*mass*volFact; 
  out[4] += mass*volFact*(1.414213562373095*F_0[0]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[8]*dvpar_sq+0.1178511301977579*F_0[0]*dvpar_sq); 
  out[5] += mass*volFact*(1.414213562373095*F_0[1]*wvpar_sq+0.8164965809277261*F_0[5]*dvpar*wvpar+0.105409255338946*F_0[9]*dvpar_sq+0.1178511301977579*F_0[1]*dvpar_sq); 
  out[6] += mass*volFact*(1.414213562373095*F_0[2]*wvpar_sq+0.8164965809277261*F_0[6]*dvpar*wvpar+0.105409255338946*F_0[10]*dvpar_sq+0.1178511301977579*F_0[2]*dvpar_sq); 
  out[7] += mass*volFact*(1.414213562373095*F_0[4]*wvpar_sq+0.8164965809277261*F_0[7]*dvpar*wvpar+0.105409255338946*F_0[11]*dvpar_sq+0.1178511301977579*F_0[4]*dvpar_sq); 
  out[8] += 1.414213562373095*G_1[0]*mass*volFact; 
  out[9] += 1.414213562373095*G_1[1]*mass*volFact; 
  out[10] += 1.414213562373095*G_1[2]*mass*volFact; 
  out[11] += 1.414213562373095*G_1[4]*mass*volFact; 
} 
