#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_3x1v_ser_p1(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[3]/2.0; 
  const double wvpar = w[3], dvpar = dxv[3]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[24]; 
  out[0] += 1.414213562373095*F_0[0]*mass*volFact; 
  out[1] += 1.414213562373095*F_0[1]*mass*volFact; 
  out[2] += 1.414213562373095*F_0[2]*mass*volFact; 
  out[3] += 1.414213562373095*F_0[3]*mass*volFact; 
  out[4] += 1.414213562373095*F_0[5]*mass*volFact; 
  out[5] += 1.414213562373095*F_0[6]*mass*volFact; 
  out[6] += 1.414213562373095*F_0[7]*mass*volFact; 
  out[7] += 1.414213562373095*F_0[11]*mass*volFact; 
  out[8] += mass*volFact*(1.414213562373095*F_0[0]*wvpar_sq+0.8164965809277261*F_0[4]*dvpar*wvpar+0.105409255338946*F_0[16]*dvpar_sq+0.1178511301977579*F_0[0]*dvpar_sq); 
  out[9] += mass*volFact*(1.414213562373095*F_0[1]*wvpar_sq+0.8164965809277261*F_0[8]*dvpar*wvpar+0.105409255338946*F_0[17]*dvpar_sq+0.1178511301977579*F_0[1]*dvpar_sq); 
  out[10] += mass*volFact*(1.414213562373095*F_0[2]*wvpar_sq+0.8164965809277261*F_0[9]*dvpar*wvpar+0.105409255338946*F_0[18]*dvpar_sq+0.1178511301977579*F_0[2]*dvpar_sq); 
  out[11] += mass*volFact*(1.414213562373095*F_0[3]*wvpar_sq+0.8164965809277261*F_0[10]*dvpar*wvpar+0.105409255338946*F_0[19]*dvpar_sq+0.1178511301977579*F_0[3]*dvpar_sq); 
  out[12] += mass*volFact*(1.414213562373095*F_0[5]*wvpar_sq+0.8164965809277261*F_0[12]*dvpar*wvpar+0.105409255338946*F_0[20]*dvpar_sq+0.1178511301977579*F_0[5]*dvpar_sq); 
  out[13] += mass*volFact*(1.414213562373095*F_0[6]*wvpar_sq+0.8164965809277261*F_0[13]*dvpar*wvpar+0.105409255338946*F_0[21]*dvpar_sq+0.1178511301977579*F_0[6]*dvpar_sq); 
  out[14] += mass*volFact*(1.414213562373095*F_0[7]*wvpar_sq+0.8164965809277261*F_0[14]*dvpar*wvpar+0.105409255338946*F_0[22]*dvpar_sq+0.1178511301977579*F_0[7]*dvpar_sq); 
  out[15] += mass*volFact*(1.414213562373095*F_0[11]*wvpar_sq+0.8164965809277261*F_0[15]*dvpar*wvpar+0.105409255338946*F_0[23]*dvpar_sq+0.1178511301977579*F_0[11]*dvpar_sq); 
  out[16] += 1.414213562373095*G_1[0]*mass*volFact; 
  out[17] += 1.414213562373095*G_1[1]*mass*volFact; 
  out[18] += 1.414213562373095*G_1[2]*mass*volFact; 
  out[19] += 1.414213562373095*G_1[3]*mass*volFact; 
  out[20] += 1.414213562373095*G_1[5]*mass*volFact; 
  out[21] += 1.414213562373095*G_1[6]*mass*volFact; 
  out[22] += 1.414213562373095*G_1[7]*mass*volFact; 
  out[23] += 1.414213562373095*G_1[11]*mass*volFact; 
} 
