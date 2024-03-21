#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_2x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: cell spacing. 
  // vmapl,vmapc,vmapr: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double fl_over_jacv[24], fc_over_jacv[24], fr_over_jacv[24];
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
  fl_over_jacv[12] = fl[12]/jacobvel[0]; 
  fl_over_jacv[13] = fl[13]/jacobvel[0]; 
  fl_over_jacv[14] = fl[14]/jacobvel[0]; 
  fl_over_jacv[15] = fl[15]/jacobvel[0]; 
  fl_over_jacv[16] = fl[16]/jacobvel[0]; 
  fl_over_jacv[17] = fl[17]/jacobvel[0]; 
  fl_over_jacv[18] = fl[18]/jacobvel[0]; 
  fl_over_jacv[19] = fl[19]/jacobvel[0]; 
  fl_over_jacv[20] = fl[20]/jacobvel[0]; 
  fl_over_jacv[21] = fl[21]/jacobvel[0]; 
  fl_over_jacv[22] = fl[22]/jacobvel[0]; 
  fl_over_jacv[23] = fl[23]/jacobvel[0]; 

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
  fc_over_jacv[12] = fc[12]/jacobvel[0]; 
  fc_over_jacv[13] = fc[13]/jacobvel[0]; 
  fc_over_jacv[14] = fc[14]/jacobvel[0]; 
  fc_over_jacv[15] = fc[15]/jacobvel[0]; 
  fc_over_jacv[16] = fc[16]/jacobvel[0]; 
  fc_over_jacv[17] = fc[17]/jacobvel[0]; 
  fc_over_jacv[18] = fc[18]/jacobvel[0]; 
  fc_over_jacv[19] = fc[19]/jacobvel[0]; 
  fc_over_jacv[20] = fc[20]/jacobvel[0]; 
  fc_over_jacv[21] = fc[21]/jacobvel[0]; 
  fc_over_jacv[22] = fc[22]/jacobvel[0]; 
  fc_over_jacv[23] = fc[23]/jacobvel[0]; 

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
  fr_over_jacv[12] = fr[12]/jacobvel[0]; 
  fr_over_jacv[13] = fr[13]/jacobvel[0]; 
  fr_over_jacv[14] = fr[14]/jacobvel[0]; 
  fr_over_jacv[15] = fr[15]/jacobvel[0]; 
  fr_over_jacv[16] = fr[16]/jacobvel[0]; 
  fr_over_jacv[17] = fr[17]/jacobvel[0]; 
  fr_over_jacv[18] = fr[18]/jacobvel[0]; 
  fr_over_jacv[19] = fr[19]/jacobvel[0]; 
  fr_over_jacv[20] = fr[20]/jacobvel[0]; 
  fr_over_jacv[21] = fr[21]/jacobvel[0]; 
  fr_over_jacv[22] = fr[22]/jacobvel[0]; 
  fr_over_jacv[23] = fr[23]/jacobvel[0]; 

  double dvl = 2.449489742783178*vmapl[3];
  double dvc = 2.449489742783178*vmapc[3];
  double dvr = 2.449489742783178*vmapr[3];

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double diffFac[4] = {0.}; 
  diffFac[0] = bmag_inv[1]*nuVtSqSum[1]*vmap_prime[1]*m_+bmag_inv[0]*nuVtSqSum[0]*vmap_prime[1]*m_; 
  diffFac[1] = bmag_inv[0]*nuVtSqSum[1]*vmap_prime[1]*m_+nuVtSqSum[0]*bmag_inv[1]*vmap_prime[1]*m_; 
  diffFac[2] = bmag_inv[1]*vmap_prime[1]*nuVtSqSum[3]*m_+bmag_inv[0]*vmap_prime[1]*nuVtSqSum[2]*m_; 
  diffFac[3] = bmag_inv[0]*vmap_prime[1]*nuVtSqSum[3]*m_+bmag_inv[1]*vmap_prime[1]*nuVtSqSum[2]*m_; 

  double fprimel[24] = {0.0}; 

  fprimel[0] = 1.325825214724776*vmapc[3]*fl_over_jacv[4]-0.7654655446197428*vmapc[2]*fl_over_jacv[4]+1.325825214724776*vmapc[3]*fc_over_jacv[4]-0.7654655446197428*vmapc[2]*fc_over_jacv[4]+1.377837980315537*fl_over_jacv[0]*vmapc[3]-1.377837980315537*fc_over_jacv[0]*vmapc[3]-0.7954951288348656*fl_over_jacv[0]*vmapc[2]+0.7954951288348656*fc_over_jacv[0]*vmapc[2]; 
  fprimel[1] = 1.325825214724776*vmapc[3]*fl_over_jacv[8]-0.7654655446197428*vmapc[2]*fl_over_jacv[8]+1.325825214724776*vmapc[3]*fc_over_jacv[8]-0.7654655446197428*vmapc[2]*fc_over_jacv[8]+1.377837980315537*fl_over_jacv[1]*vmapc[3]-1.377837980315537*fc_over_jacv[1]*vmapc[3]-0.7954951288348656*fl_over_jacv[1]*vmapc[2]+0.7954951288348656*fc_over_jacv[1]*vmapc[2]; 
  fprimel[2] = (2.651650429449552*vmapc[3]*fl_over_jacv[9]-1.530931089239486*vmapc[2]*fl_over_jacv[9]+2.755675960631073*fl_over_jacv[2]*vmapc[3]-1.590990257669731*fl_over_jacv[2]*vmapc[2])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[9]-1.530931089239486*vmapc[2]*fc_over_jacv[9]-2.755675960631073*fc_over_jacv[2]*vmapc[3]+1.590990257669731*fc_over_jacv[2]*vmapc[2])/dvc; 
  fprimel[3] = 1.325825214724776*vmapc[3]*fl_over_jacv[10]-0.7654655446197428*vmapc[2]*fl_over_jacv[10]+1.325825214724776*vmapc[3]*fc_over_jacv[10]-0.7654655446197428*vmapc[2]*fc_over_jacv[10]+1.377837980315537*fl_over_jacv[3]*vmapc[3]-1.377837980315537*fc_over_jacv[3]*vmapc[3]-0.7954951288348656*vmapc[2]*fl_over_jacv[3]+0.7954951288348656*vmapc[2]*fc_over_jacv[3]; 
  fprimel[5] = (2.651650429449552*vmapc[3]*fl_over_jacv[12]-1.530931089239486*vmapc[2]*fl_over_jacv[12]+2.755675960631073*vmapc[3]*fl_over_jacv[5]-1.590990257669731*vmapc[2]*fl_over_jacv[5])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[12]-1.530931089239486*vmapc[2]*fc_over_jacv[12]-2.755675960631073*vmapc[3]*fc_over_jacv[5]+1.590990257669731*vmapc[2]*fc_over_jacv[5])/dvc; 
  fprimel[6] = 1.325825214724776*vmapc[3]*fl_over_jacv[13]-0.7654655446197428*vmapc[2]*fl_over_jacv[13]+1.325825214724776*vmapc[3]*fc_over_jacv[13]-0.7654655446197428*vmapc[2]*fc_over_jacv[13]+1.377837980315537*vmapc[3]*fl_over_jacv[6]-0.7954951288348656*vmapc[2]*fl_over_jacv[6]-1.377837980315537*vmapc[3]*fc_over_jacv[6]+0.7954951288348656*vmapc[2]*fc_over_jacv[6]; 
  fprimel[7] = (2.651650429449552*vmapc[3]*fl_over_jacv[14]-1.530931089239486*vmapc[2]*fl_over_jacv[14]+2.755675960631073*vmapc[3]*fl_over_jacv[7]-1.590990257669731*vmapc[2]*fl_over_jacv[7])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[14]-1.530931089239486*vmapc[2]*fc_over_jacv[14]-2.755675960631073*vmapc[3]*fc_over_jacv[7]+1.590990257669731*vmapc[2]*fc_over_jacv[7])/dvc; 
  fprimel[11] = (2.651650429449552*vmapc[3]*fl_over_jacv[15]-1.530931089239486*vmapc[2]*fl_over_jacv[15]+2.755675960631073*vmapc[3]*fl_over_jacv[11]-1.590990257669731*vmapc[2]*fl_over_jacv[11])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[15]-1.530931089239486*vmapc[2]*fc_over_jacv[15]-2.755675960631073*vmapc[3]*fc_over_jacv[11]+1.590990257669731*vmapc[2]*fc_over_jacv[11])/dvc; 
  fprimel[16] = 1.325825214724776*vmapc[3]*fl_over_jacv[19]-0.7654655446197429*vmapc[2]*fl_over_jacv[19]+1.325825214724776*vmapc[3]*fc_over_jacv[19]-0.7654655446197429*vmapc[2]*fc_over_jacv[19]+1.377837980315537*vmapc[3]*fl_over_jacv[16]-0.7954951288348656*vmapc[2]*fl_over_jacv[16]-1.377837980315537*vmapc[3]*fc_over_jacv[16]+0.7954951288348656*vmapc[2]*fc_over_jacv[16]; 
  fprimel[17] = 1.325825214724776*vmapc[3]*fl_over_jacv[21]-0.7654655446197429*vmapc[2]*fl_over_jacv[21]+1.325825214724776*vmapc[3]*fc_over_jacv[21]-0.7654655446197429*vmapc[2]*fc_over_jacv[21]+1.377837980315537*vmapc[3]*fl_over_jacv[17]-0.7954951288348656*vmapc[2]*fl_over_jacv[17]-1.377837980315537*vmapc[3]*fc_over_jacv[17]+0.7954951288348656*vmapc[2]*fc_over_jacv[17]; 
  fprimel[18] = (2.651650429449552*vmapc[3]*fl_over_jacv[22]-1.530931089239486*vmapc[2]*fl_over_jacv[22]+2.755675960631073*vmapc[3]*fl_over_jacv[18]-1.590990257669731*vmapc[2]*fl_over_jacv[18])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[22]-1.530931089239486*vmapc[2]*fc_over_jacv[22]-2.755675960631073*vmapc[3]*fc_over_jacv[18]+1.590990257669731*vmapc[2]*fc_over_jacv[18])/dvc; 
  fprimel[20] = (2.651650429449552*vmapc[3]*fl_over_jacv[23]-1.530931089239486*vmapc[2]*fl_over_jacv[23]+2.755675960631073*vmapc[3]*fl_over_jacv[20]-1.590990257669731*vmapc[2]*fl_over_jacv[20])/dvl+(2.651650429449552*vmapc[3]*fc_over_jacv[23]-1.530931089239486*vmapc[2]*fc_over_jacv[23]-2.755675960631073*vmapc[3]*fc_over_jacv[20]+1.590990257669731*vmapc[2]*fc_over_jacv[20])/dvc; 

  double fprimer[24] = {0.0}; 

  fprimer[0] = (-1.325825214724776*vmapc[3]*fr_over_jacv[4])-0.7654655446197428*vmapc[2]*fr_over_jacv[4]-1.325825214724776*vmapc[3]*fc_over_jacv[4]-0.7654655446197428*vmapc[2]*fc_over_jacv[4]+1.377837980315537*fr_over_jacv[0]*vmapc[3]-1.377837980315537*fc_over_jacv[0]*vmapc[3]+0.7954951288348656*fr_over_jacv[0]*vmapc[2]-0.7954951288348656*fc_over_jacv[0]*vmapc[2]; 
  fprimer[1] = (-1.325825214724776*vmapc[3]*fr_over_jacv[8])-0.7654655446197428*vmapc[2]*fr_over_jacv[8]-1.325825214724776*vmapc[3]*fc_over_jacv[8]-0.7654655446197428*vmapc[2]*fc_over_jacv[8]+1.377837980315537*fr_over_jacv[1]*vmapc[3]-1.377837980315537*fc_over_jacv[1]*vmapc[3]+0.7954951288348656*fr_over_jacv[1]*vmapc[2]-0.7954951288348656*fc_over_jacv[1]*vmapc[2]; 
  fprimer[2] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[9])-1.530931089239486*vmapc[2]*fr_over_jacv[9]+2.755675960631073*fr_over_jacv[2]*vmapc[3]+1.590990257669731*fr_over_jacv[2]*vmapc[2])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[9])-1.530931089239486*vmapc[2]*fc_over_jacv[9]-2.755675960631073*fc_over_jacv[2]*vmapc[3]-1.590990257669731*fc_over_jacv[2]*vmapc[2])/dvc; 
  fprimer[3] = (-1.325825214724776*vmapc[3]*fr_over_jacv[10])-0.7654655446197428*vmapc[2]*fr_over_jacv[10]-1.325825214724776*vmapc[3]*fc_over_jacv[10]-0.7654655446197428*vmapc[2]*fc_over_jacv[10]+1.377837980315537*fr_over_jacv[3]*vmapc[3]-1.377837980315537*fc_over_jacv[3]*vmapc[3]+0.7954951288348656*vmapc[2]*fr_over_jacv[3]-0.7954951288348656*vmapc[2]*fc_over_jacv[3]; 
  fprimer[5] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[12])-1.530931089239486*vmapc[2]*fr_over_jacv[12]+2.755675960631073*vmapc[3]*fr_over_jacv[5]+1.590990257669731*vmapc[2]*fr_over_jacv[5])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[12])-1.530931089239486*vmapc[2]*fc_over_jacv[12]-2.755675960631073*vmapc[3]*fc_over_jacv[5]-1.590990257669731*vmapc[2]*fc_over_jacv[5])/dvc; 
  fprimer[6] = (-1.325825214724776*vmapc[3]*fr_over_jacv[13])-0.7654655446197428*vmapc[2]*fr_over_jacv[13]-1.325825214724776*vmapc[3]*fc_over_jacv[13]-0.7654655446197428*vmapc[2]*fc_over_jacv[13]+1.377837980315537*vmapc[3]*fr_over_jacv[6]+0.7954951288348656*vmapc[2]*fr_over_jacv[6]-1.377837980315537*vmapc[3]*fc_over_jacv[6]-0.7954951288348656*vmapc[2]*fc_over_jacv[6]; 
  fprimer[7] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[14])-1.530931089239486*vmapc[2]*fr_over_jacv[14]+2.755675960631073*vmapc[3]*fr_over_jacv[7]+1.590990257669731*vmapc[2]*fr_over_jacv[7])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[14])-1.530931089239486*vmapc[2]*fc_over_jacv[14]-2.755675960631073*vmapc[3]*fc_over_jacv[7]-1.590990257669731*vmapc[2]*fc_over_jacv[7])/dvc; 
  fprimer[11] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[15])-1.530931089239486*vmapc[2]*fr_over_jacv[15]+2.755675960631073*vmapc[3]*fr_over_jacv[11]+1.590990257669731*vmapc[2]*fr_over_jacv[11])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[15])-1.530931089239486*vmapc[2]*fc_over_jacv[15]-2.755675960631073*vmapc[3]*fc_over_jacv[11]-1.590990257669731*vmapc[2]*fc_over_jacv[11])/dvc; 
  fprimer[16] = (-1.325825214724776*vmapc[3]*fr_over_jacv[19])-0.7654655446197429*vmapc[2]*fr_over_jacv[19]-1.325825214724776*vmapc[3]*fc_over_jacv[19]-0.7654655446197429*vmapc[2]*fc_over_jacv[19]+1.377837980315537*vmapc[3]*fr_over_jacv[16]+0.7954951288348656*vmapc[2]*fr_over_jacv[16]-1.377837980315537*vmapc[3]*fc_over_jacv[16]-0.7954951288348656*vmapc[2]*fc_over_jacv[16]; 
  fprimer[17] = (-1.325825214724776*vmapc[3]*fr_over_jacv[21])-0.7654655446197429*vmapc[2]*fr_over_jacv[21]-1.325825214724776*vmapc[3]*fc_over_jacv[21]-0.7654655446197429*vmapc[2]*fc_over_jacv[21]+1.377837980315537*vmapc[3]*fr_over_jacv[17]+0.7954951288348656*vmapc[2]*fr_over_jacv[17]-1.377837980315537*vmapc[3]*fc_over_jacv[17]-0.7954951288348656*vmapc[2]*fc_over_jacv[17]; 
  fprimer[18] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[22])-1.530931089239486*vmapc[2]*fr_over_jacv[22]+2.755675960631073*vmapc[3]*fr_over_jacv[18]+1.590990257669731*vmapc[2]*fr_over_jacv[18])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[22])-1.530931089239486*vmapc[2]*fc_over_jacv[22]-2.755675960631073*vmapc[3]*fc_over_jacv[18]-1.590990257669731*vmapc[2]*fc_over_jacv[18])/dvc; 
  fprimer[20] = ((-2.651650429449552*vmapc[3]*fr_over_jacv[23])-1.530931089239486*vmapc[2]*fr_over_jacv[23]+2.755675960631073*vmapc[3]*fr_over_jacv[20]+1.590990257669731*vmapc[2]*fr_over_jacv[20])/dvr+((-2.651650429449552*vmapc[3]*fc_over_jacv[23])-1.530931089239486*vmapc[2]*fc_over_jacv[23]-2.755675960631073*vmapc[3]*fc_over_jacv[20]-1.590990257669731*vmapc[2]*fc_over_jacv[20])/dvc; 

  double incrl[24] = {0.0}; 
  incrl[0] = 0.5*diffFac[3]*fprimel[5]+0.5*diffFac[2]*fprimel[2]+0.5*diffFac[1]*fprimel[1]+0.5*diffFac[0]*fprimel[0]; 
  incrl[1] = 0.5*diffFac[2]*fprimel[5]+0.5*fprimel[2]*diffFac[3]+0.5*diffFac[0]*fprimel[1]+0.5*fprimel[0]*diffFac[1]; 
  incrl[2] = 0.5*diffFac[1]*fprimel[5]+0.5*fprimel[1]*diffFac[3]+0.5*diffFac[0]*fprimel[2]+0.5*fprimel[0]*diffFac[2]; 
  incrl[3] = 0.5*diffFac[3]*fprimel[11]+0.5*diffFac[2]*fprimel[7]+0.5*diffFac[1]*fprimel[6]+0.5*diffFac[0]*fprimel[3]; 
  incrl[5] = 0.5*diffFac[0]*fprimel[5]+0.5*fprimel[0]*diffFac[3]+0.5*diffFac[1]*fprimel[2]+0.5*fprimel[1]*diffFac[2]; 
  incrl[6] = 0.5*diffFac[2]*fprimel[11]+0.5*diffFac[3]*fprimel[7]+0.5*diffFac[0]*fprimel[6]+0.5*diffFac[1]*fprimel[3]; 
  incrl[7] = 0.5*diffFac[1]*fprimel[11]+0.5*diffFac[0]*fprimel[7]+0.5*diffFac[3]*fprimel[6]+0.5*diffFac[2]*fprimel[3]; 
  incrl[11] = 0.5*diffFac[0]*fprimel[11]+0.5*diffFac[1]*fprimel[7]+0.5*diffFac[2]*fprimel[6]+0.5*diffFac[3]*fprimel[3]; 
  incrl[16] = 0.5*diffFac[3]*fprimel[20]+0.5000000000000001*diffFac[2]*fprimel[18]+0.5000000000000001*diffFac[1]*fprimel[17]+0.5*diffFac[0]*fprimel[16]; 
  incrl[17] = 0.5000000000000001*diffFac[2]*fprimel[20]+0.5*diffFac[3]*fprimel[18]+0.5*diffFac[0]*fprimel[17]+0.5000000000000001*diffFac[1]*fprimel[16]; 
  incrl[18] = 0.5000000000000001*diffFac[1]*fprimel[20]+0.5*diffFac[0]*fprimel[18]+0.5*diffFac[3]*fprimel[17]+0.5000000000000001*diffFac[2]*fprimel[16]; 
  incrl[20] = 0.5*diffFac[0]*fprimel[20]+0.5000000000000001*diffFac[1]*fprimel[18]+0.5000000000000001*diffFac[2]*fprimel[17]+0.5*diffFac[3]*fprimel[16]; 

  double incrr[24] = {0.0}; 
  incrr[0] = 0.5*diffFac[3]*fprimer[5]+0.5*diffFac[2]*fprimer[2]+0.5*diffFac[1]*fprimer[1]+0.5*diffFac[0]*fprimer[0]; 
  incrr[1] = 0.5*diffFac[2]*fprimer[5]+0.5*fprimer[2]*diffFac[3]+0.5*diffFac[0]*fprimer[1]+0.5*fprimer[0]*diffFac[1]; 
  incrr[2] = 0.5*diffFac[1]*fprimer[5]+0.5*fprimer[1]*diffFac[3]+0.5*diffFac[0]*fprimer[2]+0.5*fprimer[0]*diffFac[2]; 
  incrr[3] = 0.5*diffFac[3]*fprimer[11]+0.5*diffFac[2]*fprimer[7]+0.5*diffFac[1]*fprimer[6]+0.5*diffFac[0]*fprimer[3]; 
  incrr[5] = 0.5*diffFac[0]*fprimer[5]+0.5*fprimer[0]*diffFac[3]+0.5*diffFac[1]*fprimer[2]+0.5*fprimer[1]*diffFac[2]; 
  incrr[6] = 0.5*diffFac[2]*fprimer[11]+0.5*diffFac[3]*fprimer[7]+0.5*diffFac[0]*fprimer[6]+0.5*diffFac[1]*fprimer[3]; 
  incrr[7] = 0.5*diffFac[1]*fprimer[11]+0.5*diffFac[0]*fprimer[7]+0.5*diffFac[3]*fprimer[6]+0.5*diffFac[2]*fprimer[3]; 
  incrr[11] = 0.5*diffFac[0]*fprimer[11]+0.5*diffFac[1]*fprimer[7]+0.5*diffFac[2]*fprimer[6]+0.5*diffFac[3]*fprimer[3]; 
  incrr[16] = 0.5*diffFac[3]*fprimer[20]+0.5000000000000001*diffFac[2]*fprimer[18]+0.5000000000000001*diffFac[1]*fprimer[17]+0.5*diffFac[0]*fprimer[16]; 
  incrr[17] = 0.5000000000000001*diffFac[2]*fprimer[20]+0.5*diffFac[3]*fprimer[18]+0.5*diffFac[0]*fprimer[17]+0.5000000000000001*diffFac[1]*fprimer[16]; 
  incrr[18] = 0.5000000000000001*diffFac[1]*fprimer[20]+0.5*diffFac[0]*fprimer[18]+0.5*diffFac[3]*fprimer[17]+0.5000000000000001*diffFac[2]*fprimer[16]; 
  incrr[20] = 0.5*diffFac[0]*fprimer[20]+0.5000000000000001*diffFac[1]*fprimer[18]+0.5000000000000001*diffFac[2]*fprimer[17]+0.5*diffFac[3]*fprimer[16]; 

  out[0] += (incrr[0]-1.0*incrl[0])*rdvSq4; 
  out[1] += (incrr[1]-1.0*incrl[1])*rdvSq4; 
  out[2] += (incrr[2]-1.0*incrl[2])*rdvSq4; 
  out[3] += (incrr[3]-1.0*incrl[3])*rdvSq4; 
  out[5] += (incrr[5]-1.0*incrl[5])*rdvSq4; 
  out[6] += (incrr[6]-1.0*incrl[6])*rdvSq4; 
  out[7] += (incrr[7]-1.0*incrl[7])*rdvSq4; 
  out[11] += (incrr[11]-1.0*incrl[11])*rdvSq4; 
  out[16] += (incrr[16]-1.0*incrl[16])*rdvSq4; 
  out[17] += (incrr[17]-1.0*incrl[17])*rdvSq4; 
  out[18] += (incrr[18]-1.0*incrl[18])*rdvSq4; 
  out[20] += (incrr[20]-1.0*incrl[20])*rdvSq4; 

  return 0.;

} 
