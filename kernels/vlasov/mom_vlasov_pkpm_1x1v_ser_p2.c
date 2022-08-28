#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void mom_vlasov_pkpm_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double mass, const double *bvar, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wvpar = w[1], dvpar = dxv[1]; 
  const double wvpar_sq = wvpar*wvpar, dvpar_sq = dvpar*dvpar; 
  const double wvpar_cu = wvpar*wvpar*wvpar, dvpar_cu = dvpar*dvpar*dvpar; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[3]; 
  const double *bz = &bvar[6]; 

  out[0] += 1.414213562373095*f[0]*mass*volFact; 
  out[1] += 1.414213562373095*f[1]*mass*volFact; 
  out[2] += 1.414213562373095*f[4]*mass*volFact; 
  out[3] += mass*volFact*(1.414213562373095*f[0]*wvpar_sq+0.8164965809277261*f[2]*dvpar*wvpar+0.105409255338946*f[5]*dvpar_sq+0.1178511301977579*f[0]*dvpar_sq); 
  out[4] += mass*volFact*(1.414213562373095*f[1]*wvpar_sq+0.8164965809277261*f[3]*dvpar*wvpar+0.105409255338946*f[7]*dvpar_sq+0.1178511301977579*f[1]*dvpar_sq); 
  out[5] += mass*volFact*(1.414213562373095*f[4]*wvpar_sq+0.816496580927726*f[6]*dvpar*wvpar+0.1178511301977579*f[4]*dvpar_sq); 
  out[6] += mass*volFact*(0.4330127018922194*bx[2]*f[6]*dvpar*wvpar_sq+0.4330127018922193*bx[1]*f[3]*dvpar*wvpar_sq+0.4330127018922193*bx[0]*f[2]*dvpar*wvpar_sq+0.5*bx[2]*f[4]*wvpar_cu+0.5*bx[1]*f[1]*wvpar_cu+0.5*bx[0]*f[0]*wvpar_cu+0.1118033988749895*bx[1]*f[7]*dvpar_sq*wvpar+0.1118033988749895*bx[0]*f[5]*dvpar_sq*wvpar+0.125*bx[2]*f[4]*dvpar_sq*wvpar+0.125*bx[1]*f[1]*dvpar_sq*wvpar+0.125*bx[0]*f[0]*dvpar_sq*wvpar+0.02165063509461096*bx[2]*f[6]*dvpar_cu+0.02165063509461097*bx[1]*f[3]*dvpar_cu+0.02165063509461097*bx[0]*f[2]*dvpar_cu); 
  out[7] += mass*volFact*(0.3872983346207417*bx[1]*f[6]*dvpar*wvpar_sq+0.3872983346207416*bx[2]*f[3]*dvpar*wvpar_sq+0.4330127018922193*bx[0]*f[3]*dvpar*wvpar_sq+0.4330127018922193*bx[1]*f[2]*dvpar*wvpar_sq+0.4472135954999579*bx[1]*f[4]*wvpar_cu+0.4472135954999579*f[1]*bx[2]*wvpar_cu+0.5*bx[0]*f[1]*wvpar_cu+0.5*f[0]*bx[1]*wvpar_cu+0.1*bx[2]*f[7]*dvpar_sq*wvpar+0.1118033988749895*bx[0]*f[7]*dvpar_sq*wvpar+0.1118033988749895*bx[1]*f[5]*dvpar_sq*wvpar+0.1118033988749895*bx[1]*f[4]*dvpar_sq*wvpar+0.1118033988749895*f[1]*bx[2]*dvpar_sq*wvpar+0.125*bx[0]*f[1]*dvpar_sq*wvpar+0.125*f[0]*bx[1]*dvpar_sq*wvpar+0.01936491673103708*bx[1]*f[6]*dvpar_cu+0.01936491673103708*bx[2]*f[3]*dvpar_cu+0.02165063509461097*bx[0]*f[3]*dvpar_cu+0.02165063509461097*bx[1]*f[2]*dvpar_cu); 
  out[8] += mass*volFact*(0.276641667586244*bx[2]*f[6]*dvpar*wvpar_sq+0.4330127018922194*bx[0]*f[6]*dvpar*wvpar_sq+0.3872983346207416*bx[1]*f[3]*dvpar*wvpar_sq+0.4330127018922193*bx[2]*f[2]*dvpar*wvpar_sq+0.31943828249997*bx[2]*f[4]*wvpar_cu+0.5*bx[0]*f[4]*wvpar_cu+0.5*f[0]*bx[2]*wvpar_cu+0.4472135954999579*bx[1]*f[1]*wvpar_cu+0.1*bx[1]*f[7]*dvpar_sq*wvpar+0.1118033988749895*bx[2]*f[5]*dvpar_sq*wvpar+0.07985957062499249*bx[2]*f[4]*dvpar_sq*wvpar+0.125*bx[0]*f[4]*dvpar_sq*wvpar+0.125*f[0]*bx[2]*dvpar_sq*wvpar+0.1118033988749895*bx[1]*f[1]*dvpar_sq*wvpar+0.0138320833793122*bx[2]*f[6]*dvpar_cu+0.02165063509461096*bx[0]*f[6]*dvpar_cu+0.01936491673103708*bx[1]*f[3]*dvpar_cu+0.02165063509461097*bx[2]*f[2]*dvpar_cu); 
} 
