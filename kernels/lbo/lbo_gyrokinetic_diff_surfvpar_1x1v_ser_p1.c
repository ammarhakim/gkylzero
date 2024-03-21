#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_1x1v_ser_p1(const double *dxv, const double *vmapl, const double *vmapc, const double *vmapr, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[2]: cell spacing. 
  // vmapl,vmapc,vmapr: velocity space mapping.
  // vmap_prime: velocity space mapping derivative.
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double fl_over_jacv[6], fc_over_jacv[6], fr_over_jacv[6];
  fl_over_jacv[0] = fl[0]/jacobvel[0]; 
  fl_over_jacv[1] = fl[1]/jacobvel[0]; 
  fl_over_jacv[2] = fl[2]/jacobvel[0]; 
  fl_over_jacv[3] = fl[3]/jacobvel[0]; 
  fl_over_jacv[4] = fl[4]/jacobvel[0]; 
  fl_over_jacv[5] = fl[5]/jacobvel[0]; 

  fc_over_jacv[0] = fc[0]/jacobvel[0]; 
  fc_over_jacv[1] = fc[1]/jacobvel[0]; 
  fc_over_jacv[2] = fc[2]/jacobvel[0]; 
  fc_over_jacv[3] = fc[3]/jacobvel[0]; 
  fc_over_jacv[4] = fc[4]/jacobvel[0]; 
  fc_over_jacv[5] = fc[5]/jacobvel[0]; 

  fr_over_jacv[0] = fr[0]/jacobvel[0]; 
  fr_over_jacv[1] = fr[1]/jacobvel[0]; 
  fr_over_jacv[2] = fr[2]/jacobvel[0]; 
  fr_over_jacv[3] = fr[3]/jacobvel[0]; 
  fr_over_jacv[4] = fr[4]/jacobvel[0]; 
  fr_over_jacv[5] = fr[5]/jacobvel[0]; 

  double dvl = 2.449489742783178*vmapl[1];
  double dvc = 2.449489742783178*vmapc[1];
  double dvr = 2.449489742783178*vmapr[1];

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double fprimel[6] = {0.0}; 

  fprimel[0] = (-1.341640786499874*fl_over_jacv[4])+1.341640786499874*fc_over_jacv[4]-2.381569860407206*fl_over_jacv[2]-2.381569860407206*fc_over_jacv[2]-1.875*fl_over_jacv[0]+1.875*fc_over_jacv[0]; 
  fprimel[1] = ((-2.683281572999748*fl_over_jacv[5])-4.763139720814412*fl_over_jacv[3]-3.75*fl_over_jacv[1])/dvl+(2.683281572999748*fc_over_jacv[5]-4.763139720814412*fc_over_jacv[3]+3.75*fc_over_jacv[1])/dvc; 

  double fprimer[6] = {0.0}; 

  fprimer[0] = 1.341640786499874*fr_over_jacv[4]-1.341640786499874*fc_over_jacv[4]-2.381569860407206*fr_over_jacv[2]-2.381569860407206*fc_over_jacv[2]+1.875*fr_over_jacv[0]-1.875*fc_over_jacv[0]; 
  fprimer[1] = (2.683281572999748*fr_over_jacv[5]-4.763139720814412*fr_over_jacv[3]+3.75*fr_over_jacv[1])/dvr+((-2.683281572999748*fc_over_jacv[5])-4.763139720814412*fc_over_jacv[3]-3.75*fc_over_jacv[1])/dvc; 

  double incrl[6] = {0.0}; 
  incrl[0] = 0.7071067811865475*fprimel[1]*nuVtSqSum[1]*vmap_prime[2]+0.7071067811865475*fprimel[0]*nuVtSqSum[0]*vmap_prime[2]; 
  incrl[1] = 0.7071067811865475*fprimel[0]*nuVtSqSum[1]*vmap_prime[2]+0.7071067811865475*nuVtSqSum[0]*fprimel[1]*vmap_prime[2]; 

  double incrr[6] = {0.0}; 
  incrr[0] = 0.7071067811865475*fprimer[1]*nuVtSqSum[1]*vmap_prime[2]+0.7071067811865475*fprimer[0]*nuVtSqSum[0]*vmap_prime[2]; 
  incrr[1] = 0.7071067811865475*fprimer[0]*nuVtSqSum[1]*vmap_prime[2]+0.7071067811865475*nuVtSqSum[0]*fprimer[1]*vmap_prime[2]; 

  out[0] += (incrr[0]-1.0*incrl[0])*rdvSq4; 
  out[1] += (incrr[1]-1.0*incrl[1])*rdvSq4; 

  return 0.;

} 
