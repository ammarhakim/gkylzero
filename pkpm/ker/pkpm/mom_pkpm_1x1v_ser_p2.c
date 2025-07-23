#include <gkyl_mom_pkpm_kernels.h> 
GKYL_CU_DH void mom_pkpm_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
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
  out[9] += mass*volFact*(1.414213562373095*F_0[0]*wvpar+0.408248290463863*F_0[2]*dvpar); 
  out[10] += mass*volFact*(1.414213562373095*F_0[1]*wvpar+0.408248290463863*F_0[3]*dvpar); 
  out[11] += mass*volFact*(1.414213562373095*F_0[4]*wvpar+0.408248290463863*F_0[6]*dvpar); 
} 
GKYL_CU_DH void mom_pkpm_diag_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2.0; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 
  const double wvpar_qu = wvpar*wvpar*wvpar*wvpar, dvpar_qu = dvpar*dvpar*dvpar*dvpar; 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[8]; 
  out[0] += 1.414213562373095*F_0[0]*mass*volFact; 
  out[1] += 1.414213562373095*F_0[1]*mass*volFact; 
  out[2] += 1.414213562373095*F_0[4]*mass*volFact; 
  out[3] += mass*volFact*(1.414213562373095*F_0[0]*wvpar+0.408248290463863*F_0[2]*dvpar); 
  out[4] += mass*volFact*(1.414213562373095*F_0[1]*wvpar+0.408248290463863*F_0[3]*dvpar); 
  out[5] += mass*volFact*(1.414213562373095*F_0[4]*wvpar+0.408248290463863*F_0[6]*dvpar); 
  out[6] += mass*volFact*(1.414213562373095*F_0[0]*wvpar_sq+0.8164965809277261*F_0[2]*dvpar*wvpar+0.105409255338946*F_0[5]*dvpar_sq+0.1178511301977579*F_0[0]*dvpar_sq); 
  out[7] += mass*volFact*(1.414213562373095*F_0[1]*wvpar_sq+0.8164965809277261*F_0[3]*dvpar*wvpar+0.105409255338946*F_0[7]*dvpar_sq+0.1178511301977579*F_0[1]*dvpar_sq); 
  out[8] += mass*volFact*(1.414213562373095*F_0[4]*wvpar_sq+0.816496580927726*F_0[6]*dvpar*wvpar+0.1178511301977579*F_0[4]*dvpar_sq); 
  out[9] += 1.414213562373095*G_1[0]*mass*volFact; 
  out[10] += 1.414213562373095*G_1[1]*mass*volFact; 
  out[11] += 1.414213562373095*G_1[4]*mass*volFact; 
  out[12] += mass*volFact*(1.224744871391589*F_0[2]*dvpar*wvpar_sq+1.414213562373095*F_0[0]*wvpar_cu+0.3162277660168379*F_0[5]*dvpar_sq*wvpar+0.3535533905932737*F_0[0]*dvpar_sq*wvpar+0.06123724356957942*F_0[2]*dvpar_cu); 
  out[13] += mass*volFact*(1.224744871391589*F_0[3]*dvpar*wvpar_sq+1.414213562373095*F_0[1]*wvpar_cu+0.3162277660168379*F_0[7]*dvpar_sq*wvpar+0.3535533905932737*F_0[1]*dvpar_sq*wvpar+0.06123724356957942*F_0[3]*dvpar_cu); 
  out[14] += mass*volFact*(1.224744871391589*F_0[6]*dvpar*wvpar_sq+1.414213562373095*F_0[4]*wvpar_cu+0.3535533905932737*F_0[4]*dvpar_sq*wvpar+0.06123724356957942*F_0[6]*dvpar_cu); 
  out[15] += mass*volFact*(1.414213562373095*G_1[0]*wvpar+0.408248290463863*G_1[2]*dvpar); 
  out[16] += mass*volFact*(1.414213562373095*G_1[1]*wvpar+0.408248290463863*G_1[3]*dvpar); 
  out[17] += mass*volFact*(1.414213562373095*G_1[4]*wvpar+0.408248290463863*G_1[6]*dvpar); 
  out[18] += mass*volFact*(0.6324555320336759*F_0[5]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[0]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[0]*wvpar_qu+1.632993161855453*F_0[2]*dvpar*wvpar_cu+0.2449489742783178*F_0[2]*dvpar_cu*wvpar+0.02258769757263127*F_0[5]*dvpar_qu+0.01767766952966368*F_0[0]*dvpar_qu); 
  out[19] += mass*volFact*(0.632455532033676*F_0[7]*dvpar_sq*wvpar_sq+0.7071067811865475*F_0[1]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[1]*wvpar_qu+1.632993161855453*F_0[3]*dvpar*wvpar_cu+0.2449489742783178*F_0[3]*dvpar_cu*wvpar+0.02258769757263128*F_0[7]*dvpar_qu+0.01767766952966368*F_0[1]*dvpar_qu); 
  out[20] += mass*volFact*(0.7071067811865475*F_0[4]*dvpar_sq*wvpar_sq+1.414213562373095*F_0[4]*wvpar_qu+1.632993161855453*F_0[6]*dvpar*wvpar_cu+0.2449489742783177*F_0[6]*dvpar_cu*wvpar+0.01767766952966368*F_0[4]*dvpar_qu); 
  out[21] += mass*volFact*(1.414213562373095*G_1[0]*wvpar_sq+0.8164965809277261*G_1[2]*dvpar*wvpar+0.105409255338946*G_1[5]*dvpar_sq+0.1178511301977579*G_1[0]*dvpar_sq); 
  out[22] += mass*volFact*(1.414213562373095*G_1[1]*wvpar_sq+0.8164965809277261*G_1[3]*dvpar*wvpar+0.105409255338946*G_1[7]*dvpar_sq+0.1178511301977579*G_1[1]*dvpar_sq); 
  out[23] += mass*volFact*(1.414213562373095*G_1[4]*wvpar_sq+0.816496580927726*G_1[6]*dvpar*wvpar+0.1178511301977579*G_1[4]*dvpar_sq); 
} 
