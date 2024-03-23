#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_3x2v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[5]: cell spacing. 
  // vmapl,vmapc,vmapr: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double fl_over_jacv[48], fc_over_jacv[48], fr_over_jacv[48];
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
  fl_over_jacv[24] = fl[24]/jacobvel[0]; 
  fl_over_jacv[25] = fl[25]/jacobvel[0]; 
  fl_over_jacv[26] = fl[26]/jacobvel[0]; 
  fl_over_jacv[27] = fl[27]/jacobvel[0]; 
  fl_over_jacv[28] = fl[28]/jacobvel[0]; 
  fl_over_jacv[29] = fl[29]/jacobvel[0]; 
  fl_over_jacv[30] = fl[30]/jacobvel[0]; 
  fl_over_jacv[31] = fl[31]/jacobvel[0]; 
  fl_over_jacv[32] = fl[32]/jacobvel[0]; 
  fl_over_jacv[33] = fl[33]/jacobvel[0]; 
  fl_over_jacv[34] = fl[34]/jacobvel[0]; 
  fl_over_jacv[35] = fl[35]/jacobvel[0]; 
  fl_over_jacv[36] = fl[36]/jacobvel[0]; 
  fl_over_jacv[37] = fl[37]/jacobvel[0]; 
  fl_over_jacv[38] = fl[38]/jacobvel[0]; 
  fl_over_jacv[39] = fl[39]/jacobvel[0]; 
  fl_over_jacv[40] = fl[40]/jacobvel[0]; 
  fl_over_jacv[41] = fl[41]/jacobvel[0]; 
  fl_over_jacv[42] = fl[42]/jacobvel[0]; 
  fl_over_jacv[43] = fl[43]/jacobvel[0]; 
  fl_over_jacv[44] = fl[44]/jacobvel[0]; 
  fl_over_jacv[45] = fl[45]/jacobvel[0]; 
  fl_over_jacv[46] = fl[46]/jacobvel[0]; 
  fl_over_jacv[47] = fl[47]/jacobvel[0]; 

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
  fc_over_jacv[24] = fc[24]/jacobvel[0]; 
  fc_over_jacv[25] = fc[25]/jacobvel[0]; 
  fc_over_jacv[26] = fc[26]/jacobvel[0]; 
  fc_over_jacv[27] = fc[27]/jacobvel[0]; 
  fc_over_jacv[28] = fc[28]/jacobvel[0]; 
  fc_over_jacv[29] = fc[29]/jacobvel[0]; 
  fc_over_jacv[30] = fc[30]/jacobvel[0]; 
  fc_over_jacv[31] = fc[31]/jacobvel[0]; 
  fc_over_jacv[32] = fc[32]/jacobvel[0]; 
  fc_over_jacv[33] = fc[33]/jacobvel[0]; 
  fc_over_jacv[34] = fc[34]/jacobvel[0]; 
  fc_over_jacv[35] = fc[35]/jacobvel[0]; 
  fc_over_jacv[36] = fc[36]/jacobvel[0]; 
  fc_over_jacv[37] = fc[37]/jacobvel[0]; 
  fc_over_jacv[38] = fc[38]/jacobvel[0]; 
  fc_over_jacv[39] = fc[39]/jacobvel[0]; 
  fc_over_jacv[40] = fc[40]/jacobvel[0]; 
  fc_over_jacv[41] = fc[41]/jacobvel[0]; 
  fc_over_jacv[42] = fc[42]/jacobvel[0]; 
  fc_over_jacv[43] = fc[43]/jacobvel[0]; 
  fc_over_jacv[44] = fc[44]/jacobvel[0]; 
  fc_over_jacv[45] = fc[45]/jacobvel[0]; 
  fc_over_jacv[46] = fc[46]/jacobvel[0]; 
  fc_over_jacv[47] = fc[47]/jacobvel[0]; 

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
  fr_over_jacv[24] = fr[24]/jacobvel[0]; 
  fr_over_jacv[25] = fr[25]/jacobvel[0]; 
  fr_over_jacv[26] = fr[26]/jacobvel[0]; 
  fr_over_jacv[27] = fr[27]/jacobvel[0]; 
  fr_over_jacv[28] = fr[28]/jacobvel[0]; 
  fr_over_jacv[29] = fr[29]/jacobvel[0]; 
  fr_over_jacv[30] = fr[30]/jacobvel[0]; 
  fr_over_jacv[31] = fr[31]/jacobvel[0]; 
  fr_over_jacv[32] = fr[32]/jacobvel[0]; 
  fr_over_jacv[33] = fr[33]/jacobvel[0]; 
  fr_over_jacv[34] = fr[34]/jacobvel[0]; 
  fr_over_jacv[35] = fr[35]/jacobvel[0]; 
  fr_over_jacv[36] = fr[36]/jacobvel[0]; 
  fr_over_jacv[37] = fr[37]/jacobvel[0]; 
  fr_over_jacv[38] = fr[38]/jacobvel[0]; 
  fr_over_jacv[39] = fr[39]/jacobvel[0]; 
  fr_over_jacv[40] = fr[40]/jacobvel[0]; 
  fr_over_jacv[41] = fr[41]/jacobvel[0]; 
  fr_over_jacv[42] = fr[42]/jacobvel[0]; 
  fr_over_jacv[43] = fr[43]/jacobvel[0]; 
  fr_over_jacv[44] = fr[44]/jacobvel[0]; 
  fr_over_jacv[45] = fr[45]/jacobvel[0]; 
  fr_over_jacv[46] = fr[46]/jacobvel[0]; 
  fr_over_jacv[47] = fr[47]/jacobvel[0]; 

  double dvl = 2.449489742783178*vmapl[1];
  double dvc = 2.449489742783178*vmapc[1];
  double dvr = 2.449489742783178*vmapr[1];

  const double *nuVtSqSum = &nuPrimMomsSum[8];

  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double fprimel[48] = {0.0}; 

  fprimel[0] = (-1.341640786499874*fl_over_jacv[32])+1.341640786499874*fc_over_jacv[32]-2.381569860407206*fl_over_jacv[4]-2.381569860407206*fc_over_jacv[4]-1.875*fl_over_jacv[0]+1.875*fc_over_jacv[0]; 
  fprimel[1] = ((-2.683281572999748*fl_over_jacv[33])-4.763139720814412*fl_over_jacv[9]-3.75*fl_over_jacv[1])/dvl+(2.683281572999748*fc_over_jacv[33]-4.763139720814412*fc_over_jacv[9]+3.75*fc_over_jacv[1])/dvc; 
  fprimel[2] = (-1.341640786499874*fl_over_jacv[34])+1.341640786499874*fc_over_jacv[34]-2.381569860407206*fl_over_jacv[10]-2.381569860407206*fc_over_jacv[10]-1.875*fl_over_jacv[2]+1.875*fc_over_jacv[2]; 
  fprimel[3] = (-1.341640786499874*fl_over_jacv[35])+1.341640786499874*fc_over_jacv[35]-2.381569860407206*fl_over_jacv[11]-2.381569860407206*fc_over_jacv[11]-1.875*fl_over_jacv[3]+1.875*fc_over_jacv[3]; 
  fprimel[5] = (-1.341640786499874*fl_over_jacv[36])+1.341640786499874*fc_over_jacv[36]-2.381569860407206*fl_over_jacv[15]-2.381569860407206*fc_over_jacv[15]-1.875*fl_over_jacv[5]+1.875*fc_over_jacv[5]; 
  fprimel[6] = ((-2.683281572999748*fl_over_jacv[37])-4.763139720814412*fl_over_jacv[17]-3.75*fl_over_jacv[6])/dvl+(2.683281572999748*fc_over_jacv[37]-4.763139720814412*fc_over_jacv[17]+3.75*fc_over_jacv[6])/dvc; 
  fprimel[7] = ((-2.683281572999748*fl_over_jacv[38])-4.763139720814412*fl_over_jacv[18]-3.75*fl_over_jacv[7])/dvl+(2.683281572999748*fc_over_jacv[38]-4.763139720814412*fc_over_jacv[18]+3.75*fc_over_jacv[7])/dvc; 
  fprimel[8] = (-1.341640786499874*fl_over_jacv[39])+1.341640786499874*fc_over_jacv[39]-2.381569860407206*fl_over_jacv[19]-2.381569860407206*fc_over_jacv[19]-1.875*fl_over_jacv[8]+1.875*fc_over_jacv[8]; 
  fprimel[12] = ((-2.683281572999748*fl_over_jacv[40])-4.763139720814412*fl_over_jacv[23]-3.75*fl_over_jacv[12])/dvl+(2.683281572999748*fc_over_jacv[40]-4.763139720814412*fc_over_jacv[23]+3.75*fc_over_jacv[12])/dvc; 
  fprimel[13] = (-1.341640786499874*fl_over_jacv[41])+1.341640786499874*fc_over_jacv[41]-2.381569860407206*fl_over_jacv[24]-2.381569860407206*fc_over_jacv[24]-1.875*fl_over_jacv[13]+1.875*fc_over_jacv[13]; 
  fprimel[14] = (-1.341640786499874*fl_over_jacv[42])+1.341640786499874*fc_over_jacv[42]-2.381569860407206*fl_over_jacv[25]-2.381569860407206*fc_over_jacv[25]-1.875*fl_over_jacv[14]+1.875*fc_over_jacv[14]; 
  fprimel[16] = ((-2.683281572999748*fl_over_jacv[43])-4.763139720814412*fl_over_jacv[26]-3.75*fl_over_jacv[16])/dvl+(2.683281572999748*fc_over_jacv[43]-4.763139720814412*fc_over_jacv[26]+3.75*fc_over_jacv[16])/dvc; 
  fprimel[20] = ((-2.683281572999748*fl_over_jacv[44])-4.763139720814412*fl_over_jacv[28]-3.75*fl_over_jacv[20])/dvl+(2.683281572999748*fc_over_jacv[44]-4.763139720814412*fc_over_jacv[28]+3.75*fc_over_jacv[20])/dvc; 
  fprimel[21] = ((-2.683281572999748*fl_over_jacv[45])-4.763139720814412*fl_over_jacv[29]-3.75*fl_over_jacv[21])/dvl+(2.683281572999748*fc_over_jacv[45]-4.763139720814412*fc_over_jacv[29]+3.75*fc_over_jacv[21])/dvc; 
  fprimel[22] = (-1.341640786499874*fl_over_jacv[46])+1.341640786499874*fc_over_jacv[46]-2.381569860407206*fl_over_jacv[30]-2.381569860407206*fc_over_jacv[30]-1.875*fl_over_jacv[22]+1.875*fc_over_jacv[22]; 
  fprimel[27] = ((-2.683281572999748*fl_over_jacv[47])-4.763139720814412*fl_over_jacv[31]-3.75*fl_over_jacv[27])/dvl+(2.683281572999748*fc_over_jacv[47]-4.763139720814412*fc_over_jacv[31]+3.75*fc_over_jacv[27])/dvc; 

  double fprimer[48] = {0.0}; 

  fprimer[0] = 1.341640786499874*fr_over_jacv[32]-1.341640786499874*fc_over_jacv[32]-2.381569860407206*fr_over_jacv[4]-2.381569860407206*fc_over_jacv[4]+1.875*fr_over_jacv[0]-1.875*fc_over_jacv[0]; 
  fprimer[1] = (2.683281572999748*fr_over_jacv[33]-4.763139720814412*fr_over_jacv[9]+3.75*fr_over_jacv[1])/dvr+((-2.683281572999748*fc_over_jacv[33])-4.763139720814412*fc_over_jacv[9]-3.75*fc_over_jacv[1])/dvc; 
  fprimer[2] = 1.341640786499874*fr_over_jacv[34]-1.341640786499874*fc_over_jacv[34]-2.381569860407206*fr_over_jacv[10]-2.381569860407206*fc_over_jacv[10]+1.875*fr_over_jacv[2]-1.875*fc_over_jacv[2]; 
  fprimer[3] = 1.341640786499874*fr_over_jacv[35]-1.341640786499874*fc_over_jacv[35]-2.381569860407206*fr_over_jacv[11]-2.381569860407206*fc_over_jacv[11]+1.875*fr_over_jacv[3]-1.875*fc_over_jacv[3]; 
  fprimer[5] = 1.341640786499874*fr_over_jacv[36]-1.341640786499874*fc_over_jacv[36]-2.381569860407206*fr_over_jacv[15]-2.381569860407206*fc_over_jacv[15]+1.875*fr_over_jacv[5]-1.875*fc_over_jacv[5]; 
  fprimer[6] = (2.683281572999748*fr_over_jacv[37]-4.763139720814412*fr_over_jacv[17]+3.75*fr_over_jacv[6])/dvr+((-2.683281572999748*fc_over_jacv[37])-4.763139720814412*fc_over_jacv[17]-3.75*fc_over_jacv[6])/dvc; 
  fprimer[7] = (2.683281572999748*fr_over_jacv[38]-4.763139720814412*fr_over_jacv[18]+3.75*fr_over_jacv[7])/dvr+((-2.683281572999748*fc_over_jacv[38])-4.763139720814412*fc_over_jacv[18]-3.75*fc_over_jacv[7])/dvc; 
  fprimer[8] = 1.341640786499874*fr_over_jacv[39]-1.341640786499874*fc_over_jacv[39]-2.381569860407206*fr_over_jacv[19]-2.381569860407206*fc_over_jacv[19]+1.875*fr_over_jacv[8]-1.875*fc_over_jacv[8]; 
  fprimer[12] = (2.683281572999748*fr_over_jacv[40]-4.763139720814412*fr_over_jacv[23]+3.75*fr_over_jacv[12])/dvr+((-2.683281572999748*fc_over_jacv[40])-4.763139720814412*fc_over_jacv[23]-3.75*fc_over_jacv[12])/dvc; 
  fprimer[13] = 1.341640786499874*fr_over_jacv[41]-1.341640786499874*fc_over_jacv[41]-2.381569860407206*fr_over_jacv[24]-2.381569860407206*fc_over_jacv[24]+1.875*fr_over_jacv[13]-1.875*fc_over_jacv[13]; 
  fprimer[14] = 1.341640786499874*fr_over_jacv[42]-1.341640786499874*fc_over_jacv[42]-2.381569860407206*fr_over_jacv[25]-2.381569860407206*fc_over_jacv[25]+1.875*fr_over_jacv[14]-1.875*fc_over_jacv[14]; 
  fprimer[16] = (2.683281572999748*fr_over_jacv[43]-4.763139720814412*fr_over_jacv[26]+3.75*fr_over_jacv[16])/dvr+((-2.683281572999748*fc_over_jacv[43])-4.763139720814412*fc_over_jacv[26]-3.75*fc_over_jacv[16])/dvc; 
  fprimer[20] = (2.683281572999748*fr_over_jacv[44]-4.763139720814412*fr_over_jacv[28]+3.75*fr_over_jacv[20])/dvr+((-2.683281572999748*fc_over_jacv[44])-4.763139720814412*fc_over_jacv[28]-3.75*fc_over_jacv[20])/dvc; 
  fprimer[21] = (2.683281572999748*fr_over_jacv[45]-4.763139720814412*fr_over_jacv[29]+3.75*fr_over_jacv[21])/dvr+((-2.683281572999748*fc_over_jacv[45])-4.763139720814412*fc_over_jacv[29]-3.75*fc_over_jacv[21])/dvc; 
  fprimer[22] = 1.341640786499874*fr_over_jacv[46]-1.341640786499874*fc_over_jacv[46]-2.381569860407206*fr_over_jacv[30]-2.381569860407206*fc_over_jacv[30]+1.875*fr_over_jacv[22]-1.875*fc_over_jacv[22]; 
  fprimer[27] = (2.683281572999748*fr_over_jacv[47]-4.763139720814412*fr_over_jacv[31]+3.75*fr_over_jacv[27])/dvr+((-2.683281572999748*fc_over_jacv[47])-4.763139720814412*fc_over_jacv[31]-3.75*fc_over_jacv[27])/dvc; 

  double incrl[48] = {0.0}; 
  incrl[0] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[8]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[6]+0.3535533905932737*vmap_prime[1]*fprimel[3]*nuVtSqSum[3]+0.3535533905932737*vmap_prime[1]*fprimel[2]*nuVtSqSum[2]+0.3535533905932737*fprimel[1]*nuVtSqSum[1]*vmap_prime[1]+0.3535533905932737*fprimel[0]*nuVtSqSum[0]*vmap_prime[1]; 
  incrl[1] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[8]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[6]+0.3535533905932737*vmap_prime[1]*fprimel[3]*nuVtSqSum[5]+0.3535533905932737*vmap_prime[1]*fprimel[2]*nuVtSqSum[4]+0.3535533905932737*fprimel[0]*nuVtSqSum[1]*vmap_prime[1]+0.3535533905932737*nuVtSqSum[0]*fprimel[1]*vmap_prime[1]; 
  incrl[2] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[8]+0.3535533905932737*vmap_prime[1]*fprimel[7]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*fprimel[3]*nuVtSqSum[6]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[6]+0.3535533905932737*fprimel[1]*vmap_prime[1]*nuVtSqSum[4]+0.3535533905932737*fprimel[0]*vmap_prime[1]*nuVtSqSum[2]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[2]; 
  incrl[3] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[8]+0.3535533905932737*vmap_prime[1]*fprimel[6]*nuVtSqSum[7]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[7]+0.3535533905932737*vmap_prime[1]*fprimel[2]*nuVtSqSum[6]+0.3535533905932737*fprimel[1]*vmap_prime[1]*nuVtSqSum[5]+0.3535533905932737*fprimel[0]*vmap_prime[1]*nuVtSqSum[3]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[3]; 
  incrl[5] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[13]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[12]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[5]; 
  incrl[6] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[8]+0.3535533905932737*vmap_prime[1]*fprimel[3]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[7]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[6]+0.3535533905932737*fprimel[0]*vmap_prime[1]*nuVtSqSum[4]+0.3535533905932737*fprimel[1]*vmap_prime[1]*nuVtSqSum[2]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[2]; 
  incrl[7] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[8]+0.3535533905932737*vmap_prime[1]*fprimel[2]*nuVtSqSum[7]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[7]+0.3535533905932737*vmap_prime[1]*fprimel[6]*nuVtSqSum[6]+0.3535533905932737*fprimel[0]*vmap_prime[1]*nuVtSqSum[5]+0.3535533905932737*fprimel[1]*vmap_prime[1]*nuVtSqSum[3]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[3]; 
  incrl[8] = 0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[16]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[8]+0.3535533905932737*fprimel[1]*vmap_prime[1]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[7]+0.3535533905932737*fprimel[0]*vmap_prime[1]*nuVtSqSum[6]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[6]+0.3535533905932737*vmap_prime[1]*fprimel[2]*nuVtSqSum[3]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[3]; 
  incrl[12] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[13]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[12]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[5]; 
  incrl[13] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[21]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[14]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[12]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[5]; 
  incrl[14] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[22]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[20]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[12]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[5]; 
  incrl[16] = 0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[16]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[8]+0.3535533905932737*fprimel[0]*vmap_prime[1]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[7]+0.3535533905932737*fprimel[1]*vmap_prime[1]*nuVtSqSum[6]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[6]+0.3535533905932737*vmap_prime[1]*fprimel[2]*nuVtSqSum[5]+0.3535533905932737*vmap_prime[1]*fprimel[3]*nuVtSqSum[4]; 
  incrl[20] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[21]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[14]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[12]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[5]; 
  incrl[21] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[22]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[20]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[12]+0.3535533905932737*vmap_prime[1]*fprimel[5]*nuVtSqSum[5]; 
  incrl[22] = 0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[27]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimel[12]+0.3535533905932737*vmap_prime[1]*fprimel[5]*nuVtSqSum[6]; 
  incrl[27] = 0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimel[27]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimel[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimel[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimel[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimel[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimel[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimel[12]+0.3535533905932737*vmap_prime[1]*fprimel[5]*nuVtSqSum[7]; 

  double incrr[48] = {0.0}; 
  incrr[0] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[8]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[6]+0.3535533905932737*vmap_prime[1]*fprimer[3]*nuVtSqSum[3]+0.3535533905932737*vmap_prime[1]*fprimer[2]*nuVtSqSum[2]+0.3535533905932737*fprimer[1]*nuVtSqSum[1]*vmap_prime[1]+0.3535533905932737*fprimer[0]*nuVtSqSum[0]*vmap_prime[1]; 
  incrr[1] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[8]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[6]+0.3535533905932737*vmap_prime[1]*fprimer[3]*nuVtSqSum[5]+0.3535533905932737*vmap_prime[1]*fprimer[2]*nuVtSqSum[4]+0.3535533905932737*fprimer[0]*nuVtSqSum[1]*vmap_prime[1]+0.3535533905932737*nuVtSqSum[0]*fprimer[1]*vmap_prime[1]; 
  incrr[2] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[8]+0.3535533905932737*vmap_prime[1]*fprimer[7]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*fprimer[3]*nuVtSqSum[6]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[6]+0.3535533905932737*fprimer[1]*vmap_prime[1]*nuVtSqSum[4]+0.3535533905932737*fprimer[0]*vmap_prime[1]*nuVtSqSum[2]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[2]; 
  incrr[3] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[8]+0.3535533905932737*vmap_prime[1]*fprimer[6]*nuVtSqSum[7]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[7]+0.3535533905932737*vmap_prime[1]*fprimer[2]*nuVtSqSum[6]+0.3535533905932737*fprimer[1]*vmap_prime[1]*nuVtSqSum[5]+0.3535533905932737*fprimer[0]*vmap_prime[1]*nuVtSqSum[3]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[3]; 
  incrr[5] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[13]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[12]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[5]; 
  incrr[6] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[8]+0.3535533905932737*vmap_prime[1]*fprimer[3]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[7]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[6]+0.3535533905932737*fprimer[0]*vmap_prime[1]*nuVtSqSum[4]+0.3535533905932737*fprimer[1]*vmap_prime[1]*nuVtSqSum[2]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[2]; 
  incrr[7] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[16]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[8]+0.3535533905932737*vmap_prime[1]*fprimer[2]*nuVtSqSum[7]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[7]+0.3535533905932737*vmap_prime[1]*fprimer[6]*nuVtSqSum[6]+0.3535533905932737*fprimer[0]*vmap_prime[1]*nuVtSqSum[5]+0.3535533905932737*fprimer[1]*vmap_prime[1]*nuVtSqSum[3]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[3]; 
  incrr[8] = 0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[16]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[8]+0.3535533905932737*fprimer[1]*vmap_prime[1]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[7]+0.3535533905932737*fprimer[0]*vmap_prime[1]*nuVtSqSum[6]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[6]+0.3535533905932737*vmap_prime[1]*fprimer[2]*nuVtSqSum[3]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[3]; 
  incrr[12] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[13]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[12]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[5]; 
  incrr[13] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[21]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[14]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[12]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[5]; 
  incrr[14] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[22]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[20]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[12]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[5]; 
  incrr[16] = 0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[16]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[8]+0.3535533905932737*fprimer[0]*vmap_prime[1]*nuVtSqSum[7]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[7]+0.3535533905932737*fprimer[1]*vmap_prime[1]*nuVtSqSum[6]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[6]+0.3535533905932737*vmap_prime[1]*fprimer[2]*nuVtSqSum[5]+0.3535533905932737*vmap_prime[1]*fprimer[3]*nuVtSqSum[4]; 
  incrr[20] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[21]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[14]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[12]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[5]; 
  incrr[21] = 0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[27]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[22]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[20]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[12]+0.3535533905932737*vmap_prime[1]*fprimer[5]*nuVtSqSum[5]; 
  incrr[22] = 0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[27]+0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[7]*fprimer[12]+0.3535533905932737*vmap_prime[1]*fprimer[5]*nuVtSqSum[6]; 
  incrr[27] = 0.3535533905932737*nuVtSqSum[0]*vmap_prime[1]*fprimer[27]+0.3535533905932737*nuVtSqSum[1]*vmap_prime[1]*fprimer[22]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[2]*fprimer[21]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[3]*fprimer[20]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[4]*fprimer[14]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[5]*fprimer[13]+0.3535533905932737*vmap_prime[1]*nuVtSqSum[6]*fprimer[12]+0.3535533905932737*vmap_prime[1]*fprimer[5]*nuVtSqSum[7]; 

  out[0] += (incrr[0]-1.0*incrl[0])*rdvSq4; 
  out[1] += (incrr[1]-1.0*incrl[1])*rdvSq4; 
  out[2] += (incrr[2]-1.0*incrl[2])*rdvSq4; 
  out[3] += (incrr[3]-1.0*incrl[3])*rdvSq4; 
  out[5] += (incrr[5]-1.0*incrl[5])*rdvSq4; 
  out[6] += (incrr[6]-1.0*incrl[6])*rdvSq4; 
  out[7] += (incrr[7]-1.0*incrl[7])*rdvSq4; 
  out[8] += (incrr[8]-1.0*incrl[8])*rdvSq4; 
  out[12] += (incrr[12]-1.0*incrl[12])*rdvSq4; 
  out[13] += (incrr[13]-1.0*incrl[13])*rdvSq4; 
  out[14] += (incrr[14]-1.0*incrl[14])*rdvSq4; 
  out[16] += (incrr[16]-1.0*incrl[16])*rdvSq4; 
  out[20] += (incrr[20]-1.0*incrl[20])*rdvSq4; 
  out[21] += (incrr[21]-1.0*incrl[21])*rdvSq4; 
  out[22] += (incrr[22]-1.0*incrl[22])*rdvSq4; 
  out[27] += (incrr[27]-1.0*incrl[27])*rdvSq4; 

  return 0.;

} 
