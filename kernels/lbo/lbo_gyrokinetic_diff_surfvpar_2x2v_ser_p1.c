#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_2x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  double dvl = 2.449489742783178*vmapl[1];
  double dvc = 2.449489742783178*vmapc[1];
  double dvr = 2.449489742783178*vmapr[1];

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double fprimel[24] = {0.0}; 

  fprimel[0] = (-1.341640786499874*fl_over_jacv[16])+1.341640786499874*fc_over_jacv[16]-2.381569860407206*fl_over_jacv[3]-2.381569860407206*fc_over_jacv[3]-1.875*fl_over_jacv[0]+1.875*fc_over_jacv[0]; 
  fprimel[1] = ((-2.683281572999748*fl_over_jacv[17])-4.763139720814412*fl_over_jacv[6]-3.75*fl_over_jacv[1])/dvl+(2.683281572999748*fc_over_jacv[17]-4.763139720814412*fc_over_jacv[6]+3.75*fc_over_jacv[1])/dvc; 
  fprimel[2] = (-1.341640786499874*fl_over_jacv[18])+1.341640786499874*fc_over_jacv[18]-2.381569860407206*fl_over_jacv[7]-2.381569860407206*fc_over_jacv[7]-1.875*fl_over_jacv[2]+1.875*fc_over_jacv[2]; 
  fprimel[4] = (-1.341640786499874*fl_over_jacv[19])+1.341640786499874*fc_over_jacv[19]-2.381569860407206*fl_over_jacv[10]-2.381569860407206*fc_over_jacv[10]-1.875*fl_over_jacv[4]+1.875*fc_over_jacv[4]; 
  fprimel[5] = ((-2.683281572999748*fl_over_jacv[20])-4.763139720814412*fl_over_jacv[11]-3.75*fl_over_jacv[5])/dvl+(2.683281572999748*fc_over_jacv[20]-4.763139720814412*fc_over_jacv[11]+3.75*fc_over_jacv[5])/dvc; 
  fprimel[8] = ((-2.683281572999748*fl_over_jacv[21])-4.763139720814412*fl_over_jacv[13]-3.75*fl_over_jacv[8])/dvl+(2.683281572999748*fc_over_jacv[21]-4.763139720814412*fc_over_jacv[13]+3.75*fc_over_jacv[8])/dvc; 
  fprimel[9] = (-1.341640786499874*fl_over_jacv[22])+1.341640786499874*fc_over_jacv[22]-2.381569860407206*fl_over_jacv[14]-2.381569860407206*fc_over_jacv[14]-1.875*fl_over_jacv[9]+1.875*fc_over_jacv[9]; 
  fprimel[12] = ((-2.683281572999748*fl_over_jacv[23])-4.763139720814412*fl_over_jacv[15]-3.75*fl_over_jacv[12])/dvl+(2.683281572999748*fc_over_jacv[23]-4.763139720814412*fc_over_jacv[15]+3.75*fc_over_jacv[12])/dvc; 

  double fprimer[24] = {0.0}; 

  fprimer[0] = 1.341640786499874*fr_over_jacv[16]-1.341640786499874*fc_over_jacv[16]-2.381569860407206*fr_over_jacv[3]-2.381569860407206*fc_over_jacv[3]+1.875*fr_over_jacv[0]-1.875*fc_over_jacv[0]; 
  fprimer[1] = (2.683281572999748*fr_over_jacv[17]-4.763139720814412*fr_over_jacv[6]+3.75*fr_over_jacv[1])/dvr+((-2.683281572999748*fc_over_jacv[17])-4.763139720814412*fc_over_jacv[6]-3.75*fc_over_jacv[1])/dvc; 
  fprimer[2] = 1.341640786499874*fr_over_jacv[18]-1.341640786499874*fc_over_jacv[18]-2.381569860407206*fr_over_jacv[7]-2.381569860407206*fc_over_jacv[7]+1.875*fr_over_jacv[2]-1.875*fc_over_jacv[2]; 
  fprimer[4] = 1.341640786499874*fr_over_jacv[19]-1.341640786499874*fc_over_jacv[19]-2.381569860407206*fr_over_jacv[10]-2.381569860407206*fc_over_jacv[10]+1.875*fr_over_jacv[4]-1.875*fc_over_jacv[4]; 
  fprimer[5] = (2.683281572999748*fr_over_jacv[20]-4.763139720814412*fr_over_jacv[11]+3.75*fr_over_jacv[5])/dvr+((-2.683281572999748*fc_over_jacv[20])-4.763139720814412*fc_over_jacv[11]-3.75*fc_over_jacv[5])/dvc; 
  fprimer[8] = (2.683281572999748*fr_over_jacv[21]-4.763139720814412*fr_over_jacv[13]+3.75*fr_over_jacv[8])/dvr+((-2.683281572999748*fc_over_jacv[21])-4.763139720814412*fc_over_jacv[13]-3.75*fc_over_jacv[8])/dvc; 
  fprimer[9] = 1.341640786499874*fr_over_jacv[22]-1.341640786499874*fc_over_jacv[22]-2.381569860407206*fr_over_jacv[14]-2.381569860407206*fc_over_jacv[14]+1.875*fr_over_jacv[9]-1.875*fc_over_jacv[9]; 
  fprimer[12] = (2.683281572999748*fr_over_jacv[23]-4.763139720814412*fr_over_jacv[15]+3.75*fr_over_jacv[12])/dvr+((-2.683281572999748*fc_over_jacv[23])-4.763139720814412*fc_over_jacv[15]-3.75*fc_over_jacv[12])/dvc; 

  double incrl[24] = {0.0}; 
  incrl[0] = 0.5*vmap_prime[1]*nuVtSqSum[3]*fprimel[5]+0.5*vmap_prime[1]*fprimel[2]*nuVtSqSum[2]+0.5*fprimel[1]*nuVtSqSum[1]*vmap_prime[1]+0.5*fprimel[0]*nuVtSqSum[0]*vmap_prime[1]; 
  incrl[1] = 0.5*vmap_prime[1]*nuVtSqSum[2]*fprimel[5]+0.5*vmap_prime[1]*fprimel[2]*nuVtSqSum[3]+0.5*fprimel[0]*nuVtSqSum[1]*vmap_prime[1]+0.5*nuVtSqSum[0]*fprimel[1]*vmap_prime[1]; 
  incrl[2] = 0.5*nuVtSqSum[1]*vmap_prime[1]*fprimel[5]+0.5*fprimel[1]*vmap_prime[1]*nuVtSqSum[3]+0.5*fprimel[0]*vmap_prime[1]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimel[2]; 
  incrl[4] = 0.5*vmap_prime[1]*nuVtSqSum[3]*fprimel[12]+0.5*vmap_prime[1]*nuVtSqSum[2]*fprimel[9]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimel[8]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimel[4]; 
  incrl[5] = 0.5*nuVtSqSum[0]*vmap_prime[1]*fprimel[5]+0.5*fprimel[0]*vmap_prime[1]*nuVtSqSum[3]+0.5*fprimel[1]*vmap_prime[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimel[2]; 
  incrl[8] = 0.5*vmap_prime[1]*nuVtSqSum[2]*fprimel[12]+0.5*vmap_prime[1]*nuVtSqSum[3]*fprimel[9]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimel[8]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimel[4]; 
  incrl[9] = 0.5*nuVtSqSum[1]*vmap_prime[1]*fprimel[12]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimel[9]+0.5*vmap_prime[1]*nuVtSqSum[3]*fprimel[8]+0.5*vmap_prime[1]*nuVtSqSum[2]*fprimel[4]; 
  incrl[12] = 0.5*nuVtSqSum[0]*vmap_prime[1]*fprimel[12]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimel[9]+0.5*vmap_prime[1]*nuVtSqSum[2]*fprimel[8]+0.5*vmap_prime[1]*nuVtSqSum[3]*fprimel[4]; 

  double incrr[24] = {0.0}; 
  incrr[0] = 0.5*vmap_prime[1]*nuVtSqSum[3]*fprimer[5]+0.5*vmap_prime[1]*fprimer[2]*nuVtSqSum[2]+0.5*fprimer[1]*nuVtSqSum[1]*vmap_prime[1]+0.5*fprimer[0]*nuVtSqSum[0]*vmap_prime[1]; 
  incrr[1] = 0.5*vmap_prime[1]*nuVtSqSum[2]*fprimer[5]+0.5*vmap_prime[1]*fprimer[2]*nuVtSqSum[3]+0.5*fprimer[0]*nuVtSqSum[1]*vmap_prime[1]+0.5*nuVtSqSum[0]*fprimer[1]*vmap_prime[1]; 
  incrr[2] = 0.5*nuVtSqSum[1]*vmap_prime[1]*fprimer[5]+0.5*fprimer[1]*vmap_prime[1]*nuVtSqSum[3]+0.5*fprimer[0]*vmap_prime[1]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimer[2]; 
  incrr[4] = 0.5*vmap_prime[1]*nuVtSqSum[3]*fprimer[12]+0.5*vmap_prime[1]*nuVtSqSum[2]*fprimer[9]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimer[8]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimer[4]; 
  incrr[5] = 0.5*nuVtSqSum[0]*vmap_prime[1]*fprimer[5]+0.5*fprimer[0]*vmap_prime[1]*nuVtSqSum[3]+0.5*fprimer[1]*vmap_prime[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimer[2]; 
  incrr[8] = 0.5*vmap_prime[1]*nuVtSqSum[2]*fprimer[12]+0.5*vmap_prime[1]*nuVtSqSum[3]*fprimer[9]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimer[8]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimer[4]; 
  incrr[9] = 0.5*nuVtSqSum[1]*vmap_prime[1]*fprimer[12]+0.5*nuVtSqSum[0]*vmap_prime[1]*fprimer[9]+0.5*vmap_prime[1]*nuVtSqSum[3]*fprimer[8]+0.5*vmap_prime[1]*nuVtSqSum[2]*fprimer[4]; 
  incrr[12] = 0.5*nuVtSqSum[0]*vmap_prime[1]*fprimer[12]+0.5*nuVtSqSum[1]*vmap_prime[1]*fprimer[9]+0.5*vmap_prime[1]*nuVtSqSum[2]*fprimer[8]+0.5*vmap_prime[1]*nuVtSqSum[3]*fprimer[4]; 

  out[0] += (incrr[0]-1.0*incrl[0])*rdvSq4; 
  out[1] += (incrr[1]-1.0*incrl[1])*rdvSq4; 
  out[2] += (incrr[2]-1.0*incrl[2])*rdvSq4; 
  out[4] += (incrr[4]-1.0*incrl[4])*rdvSq4; 
  out[5] += (incrr[5]-1.0*incrl[5])*rdvSq4; 
  out[8] += (incrr[8]-1.0*incrl[8])*rdvSq4; 
  out[9] += (incrr[9]-1.0*incrl[9])*rdvSq4; 
  out[12] += (incrr[12]-1.0*incrl[12])*rdvSq4; 

  return 0.;

} 
