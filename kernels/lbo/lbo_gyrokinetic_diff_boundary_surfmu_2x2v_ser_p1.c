#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_2x2v_ser_p1(const double *dxv, const double *vmap_edge, const double *vmap_skin, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: Cell spacing. 
  // vmap_edge,vmap_skin: velocity space mapping.
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  double fedge_over_jacv[24], fskin_over_jacv[24];
  fedge_over_jacv[0] = fedge[0]/jacobvel[0]; 
  fedge_over_jacv[1] = fedge[1]/jacobvel[0]; 
  fedge_over_jacv[2] = fedge[2]/jacobvel[0]; 
  fedge_over_jacv[3] = fedge[3]/jacobvel[0]; 
  fedge_over_jacv[4] = fedge[4]/jacobvel[0]; 
  fedge_over_jacv[5] = fedge[5]/jacobvel[0]; 
  fedge_over_jacv[6] = fedge[6]/jacobvel[0]; 
  fedge_over_jacv[7] = fedge[7]/jacobvel[0]; 
  fedge_over_jacv[8] = fedge[8]/jacobvel[0]; 
  fedge_over_jacv[9] = fedge[9]/jacobvel[0]; 
  fedge_over_jacv[10] = fedge[10]/jacobvel[0]; 
  fedge_over_jacv[11] = fedge[11]/jacobvel[0]; 
  fedge_over_jacv[12] = fedge[12]/jacobvel[0]; 
  fedge_over_jacv[13] = fedge[13]/jacobvel[0]; 
  fedge_over_jacv[14] = fedge[14]/jacobvel[0]; 
  fedge_over_jacv[15] = fedge[15]/jacobvel[0]; 
  fedge_over_jacv[16] = fedge[16]/jacobvel[0]; 
  fedge_over_jacv[17] = fedge[17]/jacobvel[0]; 
  fedge_over_jacv[18] = fedge[18]/jacobvel[0]; 
  fedge_over_jacv[19] = fedge[19]/jacobvel[0]; 
  fedge_over_jacv[20] = fedge[20]/jacobvel[0]; 
  fedge_over_jacv[21] = fedge[21]/jacobvel[0]; 
  fedge_over_jacv[22] = fedge[22]/jacobvel[0]; 
  fedge_over_jacv[23] = fedge[23]/jacobvel[0]; 

  fskin_over_jacv[0] = fskin[0]/jacobvel[0]; 
  fskin_over_jacv[1] = fskin[1]/jacobvel[0]; 
  fskin_over_jacv[2] = fskin[2]/jacobvel[0]; 
  fskin_over_jacv[3] = fskin[3]/jacobvel[0]; 
  fskin_over_jacv[4] = fskin[4]/jacobvel[0]; 
  fskin_over_jacv[5] = fskin[5]/jacobvel[0]; 
  fskin_over_jacv[6] = fskin[6]/jacobvel[0]; 
  fskin_over_jacv[7] = fskin[7]/jacobvel[0]; 
  fskin_over_jacv[8] = fskin[8]/jacobvel[0]; 
  fskin_over_jacv[9] = fskin[9]/jacobvel[0]; 
  fskin_over_jacv[10] = fskin[10]/jacobvel[0]; 
  fskin_over_jacv[11] = fskin[11]/jacobvel[0]; 
  fskin_over_jacv[12] = fskin[12]/jacobvel[0]; 
  fskin_over_jacv[13] = fskin[13]/jacobvel[0]; 
  fskin_over_jacv[14] = fskin[14]/jacobvel[0]; 
  fskin_over_jacv[15] = fskin[15]/jacobvel[0]; 
  fskin_over_jacv[16] = fskin[16]/jacobvel[0]; 
  fskin_over_jacv[17] = fskin[17]/jacobvel[0]; 
  fskin_over_jacv[18] = fskin[18]/jacobvel[0]; 
  fskin_over_jacv[19] = fskin[19]/jacobvel[0]; 
  fskin_over_jacv[20] = fskin[20]/jacobvel[0]; 
  fskin_over_jacv[21] = fskin[21]/jacobvel[0]; 
  fskin_over_jacv[22] = fskin[22]/jacobvel[0]; 
  fskin_over_jacv[23] = fskin[23]/jacobvel[0]; 

  double dv_edge = 2.449489742783178*vmap_edge[3];
  double dv_skin = 2.449489742783178*vmap_skin[3];

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double facDiff[4]; 
  facDiff[0] = vmap_prime[0]*bmag_inv[1]*nuVtSqSum[1]*m_+bmag_inv[0]*nuVtSqSum[0]*vmap_prime[0]*m_; 
  facDiff[1] = bmag_inv[0]*vmap_prime[0]*nuVtSqSum[1]*m_+nuVtSqSum[0]*vmap_prime[0]*bmag_inv[1]*m_; 
  facDiff[2] = vmap_prime[0]*bmag_inv[1]*nuVtSqSum[3]*m_+bmag_inv[0]*vmap_prime[0]*nuVtSqSum[2]*m_; 
  facDiff[3] = bmag_inv[0]*vmap_prime[0]*nuVtSqSum[3]*m_+vmap_prime[0]*bmag_inv[1]*nuVtSqSum[2]*m_; 

  double surfVar_l = 0.7071067811865475*vmap_skin[2]-1.224744871391589*vmap_skin[3];
  double surfVar_r = 1.224744871391589*vmap_skin[3]+0.7071067811865475*vmap_skin[2];

  double edgeSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[24] = {0.0}; 

  edgeSurf[0] = ((-0.5412658773652741*fskin_over_jacv[4])-0.5412658773652741*fedge_over_jacv[4]-0.5625*fskin_over_jacv[0]+0.5625*fedge_over_jacv[0])*surfVar_r; 
  edgeSurf[1] = ((-0.5412658773652741*fskin_over_jacv[8])-0.5412658773652741*fedge_over_jacv[8]-0.5625*fskin_over_jacv[1]+0.5625*fedge_over_jacv[1])*surfVar_r; 
  edgeSurf[2] = (((-1.082531754730548*fskin_over_jacv[9])-1.125*fskin_over_jacv[2])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[2]-1.082531754730548*fedge_over_jacv[9])*surfVar_r)/dv_edge; 
  edgeSurf[3] = ((-0.5412658773652741*fskin_over_jacv[10])-0.5412658773652741*fedge_over_jacv[10]-0.5625*fskin_over_jacv[3]+0.5625*fedge_over_jacv[3])*surfVar_r; 
  edgeSurf[4] = ((-0.9375*fskin_over_jacv[4])-0.9375*fedge_over_jacv[4]-0.9742785792574932*fskin_over_jacv[0]+0.9742785792574932*fedge_over_jacv[0])*surfVar_r; 
  edgeSurf[5] = (((-1.082531754730548*fskin_over_jacv[12])-1.125*fskin_over_jacv[5])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[5]-1.082531754730548*fedge_over_jacv[12])*surfVar_r)/dv_edge; 
  edgeSurf[6] = ((-0.5412658773652741*fskin_over_jacv[13])-0.5412658773652741*fedge_over_jacv[13]-0.5625*fskin_over_jacv[6]+0.5625*fedge_over_jacv[6])*surfVar_r; 
  edgeSurf[7] = (((-1.082531754730548*fskin_over_jacv[14])-1.125*fskin_over_jacv[7])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[7]-1.082531754730548*fedge_over_jacv[14])*surfVar_r)/dv_edge; 
  edgeSurf[8] = ((-0.9375*fskin_over_jacv[8])-0.9375*fedge_over_jacv[8]-0.9742785792574932*fskin_over_jacv[1]+0.9742785792574932*fedge_over_jacv[1])*surfVar_r; 
  edgeSurf[9] = (((-1.875*fskin_over_jacv[9])-1.948557158514986*fskin_over_jacv[2])*surfVar_r)/dv_skin+((1.948557158514986*fedge_over_jacv[2]-1.875*fedge_over_jacv[9])*surfVar_r)/dv_edge; 
  edgeSurf[10] = ((-0.9375*fskin_over_jacv[10])-0.9375*fedge_over_jacv[10]-0.9742785792574932*fskin_over_jacv[3]+0.9742785792574932*fedge_over_jacv[3])*surfVar_r; 
  edgeSurf[11] = (((-1.082531754730548*fskin_over_jacv[15])-1.125*fskin_over_jacv[11])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[11]-1.082531754730548*fedge_over_jacv[15])*surfVar_r)/dv_edge; 
  edgeSurf[12] = (((-1.875*fskin_over_jacv[12])-1.948557158514986*fskin_over_jacv[5])*surfVar_r)/dv_skin+((1.948557158514986*fedge_over_jacv[5]-1.875*fedge_over_jacv[12])*surfVar_r)/dv_edge; 
  edgeSurf[13] = ((-0.9375*fskin_over_jacv[13])-0.9375*fedge_over_jacv[13]-0.9742785792574932*fskin_over_jacv[6]+0.9742785792574932*fedge_over_jacv[6])*surfVar_r; 
  edgeSurf[14] = (((-1.875*fskin_over_jacv[14])-1.948557158514986*fskin_over_jacv[7])*surfVar_r)/dv_skin+((1.948557158514986*fedge_over_jacv[7]-1.875*fedge_over_jacv[14])*surfVar_r)/dv_edge; 
  edgeSurf[15] = (((-1.875*fskin_over_jacv[15])-1.948557158514986*fskin_over_jacv[11])*surfVar_r)/dv_skin+((1.948557158514986*fedge_over_jacv[11]-1.875*fedge_over_jacv[15])*surfVar_r)/dv_edge; 
  edgeSurf[16] = ((-0.5412658773652742*fskin_over_jacv[19])-0.5412658773652742*fedge_over_jacv[19]-0.5625*fskin_over_jacv[16]+0.5625*fedge_over_jacv[16])*surfVar_r; 
  edgeSurf[17] = ((-0.5412658773652742*fskin_over_jacv[21])-0.5412658773652742*fedge_over_jacv[21]-0.5625*fskin_over_jacv[17]+0.5625*fedge_over_jacv[17])*surfVar_r; 
  edgeSurf[18] = (((-1.082531754730548*fskin_over_jacv[22])-1.125*fskin_over_jacv[18])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[18]-1.082531754730548*fedge_over_jacv[22])*surfVar_r)/dv_edge; 
  edgeSurf[19] = ((-0.9375*fskin_over_jacv[19])-0.9375*fedge_over_jacv[19]-0.9742785792574935*fskin_over_jacv[16]+0.9742785792574935*fedge_over_jacv[16])*surfVar_r; 
  edgeSurf[20] = (((-1.082531754730548*fskin_over_jacv[23])-1.125*fskin_over_jacv[20])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[20]-1.082531754730548*fedge_over_jacv[23])*surfVar_r)/dv_edge; 
  edgeSurf[21] = ((-0.9375*fskin_over_jacv[21])-0.9375*fedge_over_jacv[21]-0.9742785792574935*fskin_over_jacv[17]+0.9742785792574935*fedge_over_jacv[17])*surfVar_r; 
  edgeSurf[22] = (((-1.875*fskin_over_jacv[22])-1.948557158514987*fskin_over_jacv[18])*surfVar_r)/dv_skin+((1.948557158514987*fedge_over_jacv[18]-1.875*fedge_over_jacv[22])*surfVar_r)/dv_edge; 
  edgeSurf[23] = (((-1.875*fskin_over_jacv[23])-1.948557158514987*fskin_over_jacv[20])*surfVar_r)/dv_skin+((1.948557158514987*fedge_over_jacv[20]-1.875*fedge_over_jacv[23])*surfVar_r)/dv_edge; 

  edgeSurf_incr[0] = 0.5*(facDiff[3]*edgeSurf[5]+edgeSurf[2]*facDiff[2]+edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.5*(facDiff[2]*edgeSurf[5]+edgeSurf[2]*facDiff[3]+edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.5*(facDiff[1]*edgeSurf[5]+edgeSurf[1]*facDiff[3]+edgeSurf[0]*facDiff[2]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.5*(facDiff[3]*edgeSurf[11]+facDiff[2]*edgeSurf[7]+facDiff[1]*edgeSurf[6]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.5*(facDiff[3]*edgeSurf[12]+facDiff[2]*edgeSurf[9]+facDiff[1]*edgeSurf[8]+facDiff[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.5*(facDiff[0]*edgeSurf[5]+edgeSurf[0]*facDiff[3]+edgeSurf[1]*facDiff[2]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[6] = 0.5*(facDiff[2]*edgeSurf[11]+facDiff[3]*edgeSurf[7]+facDiff[0]*edgeSurf[6]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[7] = 0.5*(facDiff[1]*edgeSurf[11]+facDiff[0]*edgeSurf[7]+facDiff[3]*edgeSurf[6]+facDiff[2]*edgeSurf[3]); 
  edgeSurf_incr[8] = 0.5*(facDiff[2]*edgeSurf[12]+facDiff[3]*edgeSurf[9]+facDiff[0]*edgeSurf[8]+facDiff[1]*edgeSurf[4]); 
  edgeSurf_incr[9] = 0.5*(facDiff[1]*edgeSurf[12]+facDiff[0]*edgeSurf[9]+facDiff[3]*edgeSurf[8]+facDiff[2]*edgeSurf[4]); 
  edgeSurf_incr[10] = 0.5*(facDiff[3]*edgeSurf[15]+facDiff[2]*edgeSurf[14]+facDiff[1]*edgeSurf[13]+facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.5*(facDiff[0]*edgeSurf[11]+facDiff[1]*edgeSurf[7]+facDiff[2]*edgeSurf[6]+edgeSurf[3]*facDiff[3]); 
  edgeSurf_incr[12] = 0.5*(facDiff[0]*edgeSurf[12]+facDiff[1]*edgeSurf[9]+facDiff[2]*edgeSurf[8]+facDiff[3]*edgeSurf[4]); 
  edgeSurf_incr[13] = 0.5*(facDiff[2]*edgeSurf[15]+facDiff[3]*edgeSurf[14]+facDiff[0]*edgeSurf[13]+facDiff[1]*edgeSurf[10]); 
  edgeSurf_incr[14] = 0.5*(facDiff[1]*edgeSurf[15]+facDiff[0]*edgeSurf[14]+facDiff[3]*edgeSurf[13]+facDiff[2]*edgeSurf[10]); 
  edgeSurf_incr[15] = 0.5*(facDiff[0]*edgeSurf[15]+facDiff[1]*edgeSurf[14]+facDiff[2]*edgeSurf[13]+facDiff[3]*edgeSurf[10]); 
  edgeSurf_incr[16] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[20]+15.0*(facDiff[2]*edgeSurf[18]+facDiff[1]*edgeSurf[17])+15.0*facDiff[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[20]+15.0*(facDiff[3]*edgeSurf[18]+facDiff[0]*edgeSurf[17])+15.0*facDiff[1]*edgeSurf[16]); 
  edgeSurf_incr[18] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[20]+15.0*(facDiff[0]*edgeSurf[18]+facDiff[3]*edgeSurf[17])+15.0*facDiff[2]*edgeSurf[16]); 
  edgeSurf_incr[19] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[23]+15.0*(facDiff[2]*edgeSurf[22]+facDiff[1]*edgeSurf[21])+15.0*facDiff[0]*edgeSurf[19]); 
  edgeSurf_incr[20] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[20]+15.0*(facDiff[1]*edgeSurf[18]+facDiff[2]*edgeSurf[17])+15.0*facDiff[3]*edgeSurf[16]); 
  edgeSurf_incr[21] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[23]+15.0*(facDiff[3]*edgeSurf[22]+facDiff[0]*edgeSurf[21])+15.0*facDiff[1]*edgeSurf[19]); 
  edgeSurf_incr[22] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[23]+15.0*(facDiff[0]*edgeSurf[22]+facDiff[3]*edgeSurf[21])+15.0*facDiff[2]*edgeSurf[19]); 
  edgeSurf_incr[23] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[23]+15.0*(facDiff[1]*edgeSurf[22]+facDiff[2]*edgeSurf[21])+15.0*facDiff[3]*edgeSurf[19]); 


  } else { 

  double edgeSurf[24] = {0.0}; 

  edgeSurf[0] = (0.5412658773652741*fskin_over_jacv[4]+0.5412658773652741*fedge_over_jacv[4]-0.5625*fskin_over_jacv[0]+0.5625*fedge_over_jacv[0])*surfVar_l; 
  edgeSurf[1] = (0.5412658773652741*fskin_over_jacv[8]+0.5412658773652741*fedge_over_jacv[8]-0.5625*fskin_over_jacv[1]+0.5625*fedge_over_jacv[1])*surfVar_l; 
  edgeSurf[2] = ((1.082531754730548*fskin_over_jacv[9]-1.125*fskin_over_jacv[2])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[9]+1.125*fedge_over_jacv[2])*surfVar_l)/dv_edge; 
  edgeSurf[3] = (0.5412658773652741*fskin_over_jacv[10]+0.5412658773652741*fedge_over_jacv[10]-0.5625*fskin_over_jacv[3]+0.5625*fedge_over_jacv[3])*surfVar_l; 
  edgeSurf[4] = ((-0.9375*fskin_over_jacv[4])-0.9375*fedge_over_jacv[4]+0.9742785792574932*fskin_over_jacv[0]-0.9742785792574932*fedge_over_jacv[0])*surfVar_l; 
  edgeSurf[5] = ((1.082531754730548*fskin_over_jacv[12]-1.125*fskin_over_jacv[5])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[12]+1.125*fedge_over_jacv[5])*surfVar_l)/dv_edge; 
  edgeSurf[6] = (0.5412658773652741*fskin_over_jacv[13]+0.5412658773652741*fedge_over_jacv[13]-0.5625*fskin_over_jacv[6]+0.5625*fedge_over_jacv[6])*surfVar_l; 
  edgeSurf[7] = ((1.082531754730548*fskin_over_jacv[14]-1.125*fskin_over_jacv[7])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[14]+1.125*fedge_over_jacv[7])*surfVar_l)/dv_edge; 
  edgeSurf[8] = ((-0.9375*fskin_over_jacv[8])-0.9375*fedge_over_jacv[8]+0.9742785792574932*fskin_over_jacv[1]-0.9742785792574932*fedge_over_jacv[1])*surfVar_l; 
  edgeSurf[9] = ((1.948557158514986*fskin_over_jacv[2]-1.875*fskin_over_jacv[9])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[9])-1.948557158514986*fedge_over_jacv[2])*surfVar_l)/dv_edge; 
  edgeSurf[10] = ((-0.9375*fskin_over_jacv[10])-0.9375*fedge_over_jacv[10]+0.9742785792574932*fskin_over_jacv[3]-0.9742785792574932*fedge_over_jacv[3])*surfVar_l; 
  edgeSurf[11] = ((1.082531754730548*fskin_over_jacv[15]-1.125*fskin_over_jacv[11])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[15]+1.125*fedge_over_jacv[11])*surfVar_l)/dv_edge; 
  edgeSurf[12] = ((1.948557158514986*fskin_over_jacv[5]-1.875*fskin_over_jacv[12])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[12])-1.948557158514986*fedge_over_jacv[5])*surfVar_l)/dv_edge; 
  edgeSurf[13] = ((-0.9375*fskin_over_jacv[13])-0.9375*fedge_over_jacv[13]+0.9742785792574932*fskin_over_jacv[6]-0.9742785792574932*fedge_over_jacv[6])*surfVar_l; 
  edgeSurf[14] = ((1.948557158514986*fskin_over_jacv[7]-1.875*fskin_over_jacv[14])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[14])-1.948557158514986*fedge_over_jacv[7])*surfVar_l)/dv_edge; 
  edgeSurf[15] = ((1.948557158514986*fskin_over_jacv[11]-1.875*fskin_over_jacv[15])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[15])-1.948557158514986*fedge_over_jacv[11])*surfVar_l)/dv_edge; 
  edgeSurf[16] = (0.5412658773652742*fskin_over_jacv[19]+0.5412658773652742*fedge_over_jacv[19]-0.5625*fskin_over_jacv[16]+0.5625*fedge_over_jacv[16])*surfVar_l; 
  edgeSurf[17] = (0.5412658773652742*fskin_over_jacv[21]+0.5412658773652742*fedge_over_jacv[21]-0.5625*fskin_over_jacv[17]+0.5625*fedge_over_jacv[17])*surfVar_l; 
  edgeSurf[18] = ((1.082531754730548*fskin_over_jacv[22]-1.125*fskin_over_jacv[18])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[22]+1.125*fedge_over_jacv[18])*surfVar_l)/dv_edge; 
  edgeSurf[19] = ((-0.9375*fskin_over_jacv[19])-0.9375*fedge_over_jacv[19]+0.9742785792574935*fskin_over_jacv[16]-0.9742785792574935*fedge_over_jacv[16])*surfVar_l; 
  edgeSurf[20] = ((1.082531754730548*fskin_over_jacv[23]-1.125*fskin_over_jacv[20])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[23]+1.125*fedge_over_jacv[20])*surfVar_l)/dv_edge; 
  edgeSurf[21] = ((-0.9375*fskin_over_jacv[21])-0.9375*fedge_over_jacv[21]+0.9742785792574935*fskin_over_jacv[17]-0.9742785792574935*fedge_over_jacv[17])*surfVar_l; 
  edgeSurf[22] = ((1.948557158514987*fskin_over_jacv[18]-1.875*fskin_over_jacv[22])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[22])-1.948557158514987*fedge_over_jacv[18])*surfVar_l)/dv_edge; 
  edgeSurf[23] = ((1.948557158514987*fskin_over_jacv[20]-1.875*fskin_over_jacv[23])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[23])-1.948557158514987*fedge_over_jacv[20])*surfVar_l)/dv_edge; 

  edgeSurf_incr[0] = 0.5*(facDiff[3]*edgeSurf[5]+edgeSurf[2]*facDiff[2]+edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.5*(facDiff[2]*edgeSurf[5]+edgeSurf[2]*facDiff[3]+edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.5*(facDiff[1]*edgeSurf[5]+edgeSurf[1]*facDiff[3]+edgeSurf[0]*facDiff[2]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.5*(facDiff[3]*edgeSurf[11]+facDiff[2]*edgeSurf[7]+facDiff[1]*edgeSurf[6]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.5*(facDiff[3]*edgeSurf[12]+facDiff[2]*edgeSurf[9]+facDiff[1]*edgeSurf[8]+facDiff[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.5*(facDiff[0]*edgeSurf[5]+edgeSurf[0]*facDiff[3]+edgeSurf[1]*facDiff[2]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[6] = 0.5*(facDiff[2]*edgeSurf[11]+facDiff[3]*edgeSurf[7]+facDiff[0]*edgeSurf[6]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[7] = 0.5*(facDiff[1]*edgeSurf[11]+facDiff[0]*edgeSurf[7]+facDiff[3]*edgeSurf[6]+facDiff[2]*edgeSurf[3]); 
  edgeSurf_incr[8] = 0.5*(facDiff[2]*edgeSurf[12]+facDiff[3]*edgeSurf[9]+facDiff[0]*edgeSurf[8]+facDiff[1]*edgeSurf[4]); 
  edgeSurf_incr[9] = 0.5*(facDiff[1]*edgeSurf[12]+facDiff[0]*edgeSurf[9]+facDiff[3]*edgeSurf[8]+facDiff[2]*edgeSurf[4]); 
  edgeSurf_incr[10] = 0.5*(facDiff[3]*edgeSurf[15]+facDiff[2]*edgeSurf[14]+facDiff[1]*edgeSurf[13]+facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.5*(facDiff[0]*edgeSurf[11]+facDiff[1]*edgeSurf[7]+facDiff[2]*edgeSurf[6]+edgeSurf[3]*facDiff[3]); 
  edgeSurf_incr[12] = 0.5*(facDiff[0]*edgeSurf[12]+facDiff[1]*edgeSurf[9]+facDiff[2]*edgeSurf[8]+facDiff[3]*edgeSurf[4]); 
  edgeSurf_incr[13] = 0.5*(facDiff[2]*edgeSurf[15]+facDiff[3]*edgeSurf[14]+facDiff[0]*edgeSurf[13]+facDiff[1]*edgeSurf[10]); 
  edgeSurf_incr[14] = 0.5*(facDiff[1]*edgeSurf[15]+facDiff[0]*edgeSurf[14]+facDiff[3]*edgeSurf[13]+facDiff[2]*edgeSurf[10]); 
  edgeSurf_incr[15] = 0.5*(facDiff[0]*edgeSurf[15]+facDiff[1]*edgeSurf[14]+facDiff[2]*edgeSurf[13]+facDiff[3]*edgeSurf[10]); 
  edgeSurf_incr[16] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[20]+15.0*(facDiff[2]*edgeSurf[18]+facDiff[1]*edgeSurf[17])+15.0*facDiff[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[20]+15.0*(facDiff[3]*edgeSurf[18]+facDiff[0]*edgeSurf[17])+15.0*facDiff[1]*edgeSurf[16]); 
  edgeSurf_incr[18] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[20]+15.0*(facDiff[0]*edgeSurf[18]+facDiff[3]*edgeSurf[17])+15.0*facDiff[2]*edgeSurf[16]); 
  edgeSurf_incr[19] = 0.03333333333333333*(15.0*facDiff[3]*edgeSurf[23]+15.0*(facDiff[2]*edgeSurf[22]+facDiff[1]*edgeSurf[21])+15.0*facDiff[0]*edgeSurf[19]); 
  edgeSurf_incr[20] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[20]+15.0*(facDiff[1]*edgeSurf[18]+facDiff[2]*edgeSurf[17])+15.0*facDiff[3]*edgeSurf[16]); 
  edgeSurf_incr[21] = 0.03333333333333333*(15.0*facDiff[2]*edgeSurf[23]+15.0*(facDiff[3]*edgeSurf[22]+facDiff[0]*edgeSurf[21])+15.0*facDiff[1]*edgeSurf[19]); 
  edgeSurf_incr[22] = 0.03333333333333333*(15.0*facDiff[1]*edgeSurf[23]+15.0*(facDiff[0]*edgeSurf[22]+facDiff[3]*edgeSurf[21])+15.0*facDiff[2]*edgeSurf[19]); 
  edgeSurf_incr[23] = 0.03333333333333333*(15.0*facDiff[0]*edgeSurf[23]+15.0*(facDiff[1]*edgeSurf[22]+facDiff[2]*edgeSurf[21])+15.0*facDiff[3]*edgeSurf[19]); 

  } 

  out[0] += edgeSurf_incr[0]*rdvSq4; 
  out[1] += edgeSurf_incr[1]*rdvSq4; 
  out[2] += edgeSurf_incr[2]*rdvSq4; 
  out[3] += edgeSurf_incr[3]*rdvSq4; 
  out[4] += edgeSurf_incr[4]*rdvSq4; 
  out[5] += edgeSurf_incr[5]*rdvSq4; 
  out[6] += edgeSurf_incr[6]*rdvSq4; 
  out[7] += edgeSurf_incr[7]*rdvSq4; 
  out[8] += edgeSurf_incr[8]*rdvSq4; 
  out[9] += edgeSurf_incr[9]*rdvSq4; 
  out[10] += edgeSurf_incr[10]*rdvSq4; 
  out[11] += edgeSurf_incr[11]*rdvSq4; 
  out[12] += edgeSurf_incr[12]*rdvSq4; 
  out[13] += edgeSurf_incr[13]*rdvSq4; 
  out[14] += edgeSurf_incr[14]*rdvSq4; 
  out[15] += edgeSurf_incr[15]*rdvSq4; 
  out[16] += edgeSurf_incr[16]*rdvSq4; 
  out[17] += edgeSurf_incr[17]*rdvSq4; 
  out[18] += edgeSurf_incr[18]*rdvSq4; 
  out[19] += edgeSurf_incr[19]*rdvSq4; 
  out[20] += edgeSurf_incr[20]*rdvSq4; 
  out[21] += edgeSurf_incr[21]*rdvSq4; 
  out[22] += edgeSurf_incr[22]*rdvSq4; 
  out[23] += edgeSurf_incr[23]*rdvSq4; 

  return 0.;

} 
