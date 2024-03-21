#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_1x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: cell spacing. 
  // vmapl,vmapc,vmapr: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double fl_over_jacv[12], fc_over_jacv[12], fr_over_jacv[12];
  fl_over_jacv[0] = fl[0]/jacobvel[0]; 
  fl_over_jacv[1] = fl[1]/jacobvel[0]; 
  fl_over_jacv[2] = fl[2]/jacobvel[0]; 
  fl_over_jacv[3] = fl[3]/jacobvel[0]; 
  fl_over_jacv[4] = fl[4]/jacobvel[0]; 
  fl_over_jacv[5] = fl[5]/jacobvel[0]; 
  fl_over_jacv[6] = fl[6]/jacobvel[0]; 
  fl_over_jacv[7] = fl[7]/jacobvel[0]; 
  fl_over_jacv[8] = fl[8]/jacobvel[0]; 
  fl_over_jacv[9] = fl[9]/jacobvel[0]; 
  fl_over_jacv[10] = fl[10]/jacobvel[0]; 
  fl_over_jacv[11] = fl[11]/jacobvel[0]; 

  fc_over_jacv[0] = fc[0]/jacobvel[0]; 
  fc_over_jacv[1] = fc[1]/jacobvel[0]; 
  fc_over_jacv[2] = fc[2]/jacobvel[0]; 
  fc_over_jacv[3] = fc[3]/jacobvel[0]; 
  fc_over_jacv[4] = fc[4]/jacobvel[0]; 
  fc_over_jacv[5] = fc[5]/jacobvel[0]; 
  fc_over_jacv[6] = fc[6]/jacobvel[0]; 
  fc_over_jacv[7] = fc[7]/jacobvel[0]; 
  fc_over_jacv[8] = fc[8]/jacobvel[0]; 
  fc_over_jacv[9] = fc[9]/jacobvel[0]; 
  fc_over_jacv[10] = fc[10]/jacobvel[0]; 
  fc_over_jacv[11] = fc[11]/jacobvel[0]; 

  fr_over_jacv[0] = fr[0]/jacobvel[0]; 
  fr_over_jacv[1] = fr[1]/jacobvel[0]; 
  fr_over_jacv[2] = fr[2]/jacobvel[0]; 
  fr_over_jacv[3] = fr[3]/jacobvel[0]; 
  fr_over_jacv[4] = fr[4]/jacobvel[0]; 
  fr_over_jacv[5] = fr[5]/jacobvel[0]; 
  fr_over_jacv[6] = fr[6]/jacobvel[0]; 
  fr_over_jacv[7] = fr[7]/jacobvel[0]; 
  fr_over_jacv[8] = fr[8]/jacobvel[0]; 
  fr_over_jacv[9] = fr[9]/jacobvel[0]; 
  fr_over_jacv[10] = fr[10]/jacobvel[0]; 
  fr_over_jacv[11] = fr[11]/jacobvel[0]; 

  double dvl = 2.449489742783178*vmapl[3];
  double dvc = 2.449489742783178*vmapc[3];
  double dvr = 2.449489742783178*vmapr[3];

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double diffFac[2] = {0.}; 
  diffFac[0] = 1.414213562373095*bmag_inv[1]*nuVtSqSum[1]*vmap_prime[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*vmap_prime[1]*m_; 
  diffFac[1] = 1.414213562373095*bmag_inv[0]*nuVtSqSum[1]*vmap_prime[1]*m_+1.414213562373095*nuVtSqSum[0]*bmag_inv[1]*vmap_prime[1]*m_; 

  double fprimel[12] = {0.0}; 
  const double dvlR2 = pow(dvl,2);
  const double dvcR2 = pow(dvc,2);

  fprimel[0] = (5.929270612815709*vmapc[3]*fl_over_jacv[10]-3.423265984407287*vmapc[2]*fl_over_jacv[10]+6.161878771933115*vmapc[3]*fl_over_jacv[8]-3.557562367689425*vmapc[2]*fl_over_jacv[8])/dvlR2+(5.929270612815709*vmapc[3]*fc_over_jacv[10]-3.423265984407287*vmapc[2]*fc_over_jacv[10]-6.161878771933115*vmapc[3]*fc_over_jacv[8]+3.557562367689425*vmapc[2]*fc_over_jacv[8])/dvcR2-1.482317653203927*vmapc[3]*fl_over_jacv[10]+0.8558164961018215*vmapc[2]*fl_over_jacv[10]-1.482317653203927*vmapc[3]*fc_over_jacv[10]+0.8558164961018215*vmapc[2]*fc_over_jacv[10]-1.540469692983278*vmapc[3]*fl_over_jacv[8]+0.8893905919223563*vmapc[2]*fl_over_jacv[8]+1.540469692983278*vmapc[3]*fc_over_jacv[8]-0.8893905919223563*vmapc[2]*fc_over_jacv[8]+1.325825214724776*fl_over_jacv[3]*vmapc[3]+1.325825214724776*fc_over_jacv[3]*vmapc[3]+1.377837980315537*fl_over_jacv[0]*vmapc[3]-1.377837980315537*fc_over_jacv[0]*vmapc[3]-0.7654655446197428*vmapc[2]*fl_over_jacv[3]-0.7654655446197428*vmapc[2]*fc_over_jacv[3]-0.7954951288348656*fl_over_jacv[0]*vmapc[2]+0.7954951288348656*fc_over_jacv[0]*vmapc[2]; 
  fprimel[1] = (5.929270612815711*vmapc[3]*fl_over_jacv[11]-3.423265984407287*vmapc[2]*fl_over_jacv[11]+6.161878771933116*vmapc[3]*fl_over_jacv[9]-3.557562367689425*vmapc[2]*fl_over_jacv[9])/dvlR2+(5.929270612815711*vmapc[3]*fc_over_jacv[11]-3.423265984407287*vmapc[2]*fc_over_jacv[11]-6.161878771933116*vmapc[3]*fc_over_jacv[9]+3.557562367689425*vmapc[2]*fc_over_jacv[9])/dvcR2-1.482317653203927*vmapc[3]*fl_over_jacv[11]+0.8558164961018216*vmapc[2]*fl_over_jacv[11]-1.482317653203927*vmapc[3]*fc_over_jacv[11]+0.8558164961018216*vmapc[2]*fc_over_jacv[11]-1.540469692983279*vmapc[3]*fl_over_jacv[9]+0.8893905919223561*vmapc[2]*fl_over_jacv[9]+1.540469692983279*vmapc[3]*fc_over_jacv[9]-0.8893905919223561*vmapc[2]*fc_over_jacv[9]+1.325825214724776*vmapc[3]*fl_over_jacv[5]-0.7654655446197428*vmapc[2]*fl_over_jacv[5]+1.325825214724776*vmapc[3]*fc_over_jacv[5]-0.7654655446197428*vmapc[2]*fc_over_jacv[5]+1.377837980315537*fl_over_jacv[1]*vmapc[3]-1.377837980315537*fc_over_jacv[1]*vmapc[3]-0.7954951288348656*fl_over_jacv[1]*vmapc[2]+0.7954951288348656*fc_over_jacv[1]*vmapc[2]; 
  fprimel[2] = (2.651650429449552*vmapc[3]*fl_over_jacv[6]-1.530931089239486*vmapc[2]*fl_over_jacv[6]+2.755675960631073*fl_over_jacv[2]*vmapc[3]-1.590990257669731*fl_over_jacv[2]*vmapc[2])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[6]-1.530931089239486*vmapc[2]*fc_over_jacv[6]-2.755675960631073*fc_over_jacv[2]*vmapc[3]+1.590990257669731*fc_over_jacv[2]*vmapc[2])/dvc; 
  fprimel[4] = (2.651650429449552*vmapc[3]*fl_over_jacv[7]-1.530931089239486*vmapc[2]*fl_over_jacv[7]+2.755675960631073*vmapc[3]*fl_over_jacv[4]-1.590990257669731*vmapc[2]*fl_over_jacv[4])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[7]-1.530931089239486*vmapc[2]*fc_over_jacv[7]-2.755675960631073*vmapc[3]*fc_over_jacv[4]+1.590990257669731*vmapc[2]*fc_over_jacv[4])/dvc; 
  fprimel[8] = (5.303300858899106*vmapc[3]*fl_over_jacv[10]-3.061862178478972*vmapc[2]*fl_over_jacv[10]+5.511351921262148*vmapc[3]*fl_over_jacv[8]-3.181980515339463*vmapc[2]*fl_over_jacv[8])/dvlR2+(5.303300858899106*vmapc[3]*fc_over_jacv[10]-3.061862178478972*vmapc[2]*fc_over_jacv[10]-5.511351921262148*vmapc[3]*fc_over_jacv[8]+3.181980515339463*vmapc[2]*fc_over_jacv[8])/dvcR2; 
  fprimel[9] = (5.303300858899106*vmapc[3]*fl_over_jacv[11]-3.061862178478972*vmapc[2]*fl_over_jacv[11]+5.511351921262148*vmapc[3]*fl_over_jacv[9]-3.181980515339463*vmapc[2]*fl_over_jacv[9])/dvlR2+(5.303300858899106*vmapc[3]*fc_over_jacv[11]-3.061862178478972*vmapc[2]*fc_over_jacv[11]-5.511351921262148*vmapc[3]*fc_over_jacv[9]+3.181980515339463*vmapc[2]*fc_over_jacv[9])/dvcR2; 

  double fprimer[12] = {0.0}; 
  const double dvrR2 = pow(dvr,2);

  fprimer[0] = ((-5.929270612815709*vmapc[3]*fr_over_jacv[10])-3.423265984407287*vmapc[2]*fr_over_jacv[10]+6.161878771933115*vmapc[3]*fr_over_jacv[8]+3.557562367689425*vmapc[2]*fr_over_jacv[8])/dvrR2+((-5.929270612815709*vmapc[3]*fc_over_jacv[10])-3.423265984407287*vmapc[2]*fc_over_jacv[10]-6.161878771933115*vmapc[3]*fc_over_jacv[8]-3.557562367689425*vmapc[2]*fc_over_jacv[8])/dvcR2+1.482317653203927*vmapc[3]*fr_over_jacv[10]+0.8558164961018215*vmapc[2]*fr_over_jacv[10]+1.482317653203927*vmapc[3]*fc_over_jacv[10]+0.8558164961018215*vmapc[2]*fc_over_jacv[10]-1.540469692983278*vmapc[3]*fr_over_jacv[8]-0.8893905919223563*vmapc[2]*fr_over_jacv[8]+1.540469692983278*vmapc[3]*fc_over_jacv[8]+0.8893905919223563*vmapc[2]*fc_over_jacv[8]-1.325825214724776*fr_over_jacv[3]*vmapc[3]-1.325825214724776*fc_over_jacv[3]*vmapc[3]+1.377837980315537*fr_over_jacv[0]*vmapc[3]-1.377837980315537*fc_over_jacv[0]*vmapc[3]-0.7654655446197428*vmapc[2]*fr_over_jacv[3]-0.7654655446197428*vmapc[2]*fc_over_jacv[3]+0.7954951288348656*fr_over_jacv[0]*vmapc[2]-0.7954951288348656*fc_over_jacv[0]*vmapc[2]; 
  fprimer[1] = ((-5.929270612815711*vmapc[3]*fr_over_jacv[11])-3.423265984407287*vmapc[2]*fr_over_jacv[11]+6.161878771933116*vmapc[3]*fr_over_jacv[9]+3.557562367689425*vmapc[2]*fr_over_jacv[9])/dvrR2+((-5.929270612815711*vmapc[3]*fc_over_jacv[11])-3.423265984407287*vmapc[2]*fc_over_jacv[11]-6.161878771933116*vmapc[3]*fc_over_jacv[9]-3.557562367689425*vmapc[2]*fc_over_jacv[9])/dvcR2+1.482317653203927*vmapc[3]*fr_over_jacv[11]+0.8558164961018216*vmapc[2]*fr_over_jacv[11]+1.482317653203927*vmapc[3]*fc_over_jacv[11]+0.8558164961018216*vmapc[2]*fc_over_jacv[11]-1.540469692983279*vmapc[3]*fr_over_jacv[9]-0.8893905919223561*vmapc[2]*fr_over_jacv[9]+1.540469692983279*vmapc[3]*fc_over_jacv[9]+0.8893905919223561*vmapc[2]*fc_over_jacv[9]-1.325825214724776*vmapc[3]*fr_over_jacv[5]-0.7654655446197428*vmapc[2]*fr_over_jacv[5]-1.325825214724776*vmapc[3]*fc_over_jacv[5]-0.7654655446197428*vmapc[2]*fc_over_jacv[5]+1.377837980315537*fr_over_jacv[1]*vmapc[3]-1.377837980315537*fc_over_jacv[1]*vmapc[3]+0.7954951288348656*fr_over_jacv[1]*vmapc[2]-0.7954951288348656*fc_over_jacv[1]*vmapc[2]; 
  fprimer[2] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[6])-1.530931089239486*vmapc[2]*fr_over_jacv[6]+2.755675960631073*fr_over_jacv[2]*vmapc[3]+1.590990257669731*fr_over_jacv[2]*vmapc[2])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[6])-1.530931089239486*vmapc[2]*fc_over_jacv[6]-2.755675960631073*fc_over_jacv[2]*vmapc[3]-1.590990257669731*fc_over_jacv[2]*vmapc[2])/dvc; 
  fprimer[4] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[7])-1.530931089239486*vmapc[2]*fr_over_jacv[7]+2.755675960631073*vmapc[3]*fr_over_jacv[4]+1.590990257669731*vmapc[2]*fr_over_jacv[4])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[7])-1.530931089239486*vmapc[2]*fc_over_jacv[7]-2.755675960631073*vmapc[3]*fc_over_jacv[4]-1.590990257669731*vmapc[2]*fc_over_jacv[4])/dvc; 
  fprimer[8] = ((-5.303300858899106*vmapc[3]*fr_over_jacv[10])-3.061862178478972*vmapc[2]*fr_over_jacv[10]+5.511351921262148*vmapc[3]*fr_over_jacv[8]+3.181980515339463*vmapc[2]*fr_over_jacv[8])/dvrR2+((-5.303300858899106*vmapc[3]*fc_over_jacv[10])-3.061862178478972*vmapc[2]*fc_over_jacv[10]-5.511351921262148*vmapc[3]*fc_over_jacv[8]-3.181980515339463*vmapc[2]*fc_over_jacv[8])/dvcR2; 
  fprimer[9] = ((-5.303300858899106*vmapc[3]*fr_over_jacv[11])-3.061862178478972*vmapc[2]*fr_over_jacv[11]+5.511351921262148*vmapc[3]*fr_over_jacv[9]+3.181980515339463*vmapc[2]*fr_over_jacv[9])/dvrR2+((-5.303300858899106*vmapc[3]*fc_over_jacv[11])-3.061862178478972*vmapc[2]*fc_over_jacv[11]-5.511351921262148*vmapc[3]*fc_over_jacv[9]-3.181980515339463*vmapc[2]*fc_over_jacv[9])/dvcR2; 

  double incrl[12] = {0.0}; 
  incrl[0] = 0.7071067811865475*diffFac[1]*fprimel[1]+0.7071067811865475*diffFac[0]*fprimel[0]; 
  incrl[1] = 0.7071067811865475*diffFac[0]*fprimel[1]+0.7071067811865475*fprimel[0]*diffFac[1]; 
  incrl[2] = 0.7071067811865475*diffFac[1]*fprimel[4]+0.7071067811865475*diffFac[0]*fprimel[2]; 
  incrl[4] = 0.7071067811865475*diffFac[0]*fprimel[4]+0.7071067811865475*diffFac[1]*fprimel[2]; 
  incrl[8] = 0.7071067811865475*diffFac[1]*fprimel[9]+0.7071067811865475*diffFac[0]*fprimel[8]; 
  incrl[9] = 0.7071067811865475*diffFac[0]*fprimel[9]+0.7071067811865475*diffFac[1]*fprimel[8]; 

  double incrr[12] = {0.0}; 
  incrr[0] = 0.7071067811865475*diffFac[1]*fprimer[1]+0.7071067811865475*diffFac[0]*fprimer[0]; 
  incrr[1] = 0.7071067811865475*diffFac[0]*fprimer[1]+0.7071067811865475*fprimer[0]*diffFac[1]; 
  incrr[2] = 0.7071067811865475*diffFac[1]*fprimer[4]+0.7071067811865475*diffFac[0]*fprimer[2]; 
  incrr[4] = 0.7071067811865475*diffFac[0]*fprimer[4]+0.7071067811865475*diffFac[1]*fprimer[2]; 
  incrr[8] = 0.7071067811865475*diffFac[1]*fprimer[9]+0.7071067811865475*diffFac[0]*fprimer[8]; 
  incrr[9] = 0.7071067811865475*diffFac[0]*fprimer[9]+0.7071067811865475*diffFac[1]*fprimer[8]; 

  out[0] += (incrr[0]-1.0*incrl[0])*rdvSq4; 
  out[1] += (incrr[1]-1.0*incrl[1])*rdvSq4; 
  out[2] += (incrr[2]-1.0*incrl[2])*rdvSq4; 
  out[4] += (incrr[4]-1.0*incrl[4])*rdvSq4; 
  out[8] += (incrr[8]-1.0*incrl[8])*rdvSq4; 
  out[9] += (incrr[9]-1.0*incrl[9])*rdvSq4; 

  return 0.;

} 
