#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_1x2v_ser_p1(const double *dxv, const double *vmap_edge, const double *vmap_skin, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: Cell spacing. 
  // vmap_edge,vmap_skin: velocity space mapping.
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  double fedge_over_jacv[12], fskin_over_jacv[12];
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

  double dv_edge = 2.449489742783178*vmap_edge[3];
  double dv_skin = 2.449489742783178*vmap_skin[3];

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double facDiff[2]; 
  facDiff[0] = 1.414213562373095*vmap_prime[0]*bmag_inv[1]*nuVtSqSum[1]*m_+1.414213562373095*bmag_inv[0]*nuVtSqSum[0]*vmap_prime[0]*m_; 
  facDiff[1] = 1.414213562373095*bmag_inv[0]*vmap_prime[0]*nuVtSqSum[1]*m_+1.414213562373095*nuVtSqSum[0]*vmap_prime[0]*bmag_inv[1]*m_; 

  double surfVar_l = 0.7071067811865475*vmap_skin[2]-1.224744871391589*vmap_skin[3];
  double surfVar_r = 1.224744871391589*vmap_skin[3]+0.7071067811865475*vmap_skin[2];

  double edgeSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[12] = {0.0}; 
  const double dv_edgeR2 = pow(dv_edge,2);
  const double dv_skinR2 = pow(dv_skin,2);

  edgeSurf[0] = (((-2.420614591379636*fskin_over_jacv[10])-2.515576474687264*fskin_over_jacv[8])*surfVar_r)/dv_skinR2+((2.515576474687264*fedge_over_jacv[8]-2.420614591379636*fedge_over_jacv[10])*surfVar_r)/dv_edgeR2+(0.6051536478449089*fskin_over_jacv[10]+0.6051536478449089*fedge_over_jacv[10]+0.6288941186718159*fskin_over_jacv[8]-0.6288941186718159*fedge_over_jacv[8]-0.5412658773652741*fskin_over_jacv[3]-0.5412658773652741*fedge_over_jacv[3]-0.5625*fskin_over_jacv[0]+0.5625*fedge_over_jacv[0])*surfVar_r; 
  edgeSurf[1] = (((-2.420614591379636*fskin_over_jacv[11])-2.515576474687263*fskin_over_jacv[9])*surfVar_r)/dv_skinR2+((2.515576474687263*fedge_over_jacv[9]-2.420614591379636*fedge_over_jacv[11])*surfVar_r)/dv_edgeR2+(0.605153647844909*fskin_over_jacv[11]+0.605153647844909*fedge_over_jacv[11]+0.6288941186718158*fskin_over_jacv[9]-0.6288941186718158*fedge_over_jacv[9]-0.5412658773652741*fskin_over_jacv[5]-0.5412658773652741*fedge_over_jacv[5]-0.5625*fskin_over_jacv[1]+0.5625*fedge_over_jacv[1])*surfVar_r; 
  edgeSurf[2] = (((-1.082531754730548*fskin_over_jacv[6])-1.125*fskin_over_jacv[2])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[2]-1.082531754730548*fedge_over_jacv[6])*surfVar_r)/dv_edge; 
  edgeSurf[3] = (((-4.192627457812105*fskin_over_jacv[10])-4.357106264483343*fskin_over_jacv[8])*surfVar_r)/dv_skinR2+((4.357106264483343*fedge_over_jacv[8]-4.192627457812105*fedge_over_jacv[10])*surfVar_r)/dv_edgeR2+(1.048156864453026*fskin_over_jacv[10]+1.048156864453026*fedge_over_jacv[10]+1.089276566120836*fskin_over_jacv[8]-1.089276566120836*fedge_over_jacv[8]-0.9375*fskin_over_jacv[3]-0.9375*fedge_over_jacv[3]-0.9742785792574932*fskin_over_jacv[0]+0.9742785792574932*fedge_over_jacv[0])*surfVar_r; 
  edgeSurf[4] = (((-1.082531754730548*fskin_over_jacv[7])-1.125*fskin_over_jacv[4])*surfVar_r)/dv_skin+((1.125*fedge_over_jacv[4]-1.082531754730548*fedge_over_jacv[7])*surfVar_r)/dv_edge; 
  edgeSurf[5] = (((-4.192627457812106*fskin_over_jacv[11])-4.357106264483344*fskin_over_jacv[9])*surfVar_r)/dv_skinR2+((4.357106264483344*fedge_over_jacv[9]-4.192627457812106*fedge_over_jacv[11])*surfVar_r)/dv_edgeR2+(1.048156864453027*fskin_over_jacv[11]+1.048156864453027*fedge_over_jacv[11]+1.089276566120836*fskin_over_jacv[9]-1.089276566120836*fedge_over_jacv[9]-0.9375*fskin_over_jacv[5]-0.9375*fedge_over_jacv[5]-0.9742785792574932*fskin_over_jacv[1]+0.9742785792574932*fedge_over_jacv[1])*surfVar_r; 
  edgeSurf[6] = (((-1.875*fskin_over_jacv[6])-1.948557158514986*fskin_over_jacv[2])*surfVar_r)/dv_skin+((1.948557158514986*fedge_over_jacv[2]-1.875*fedge_over_jacv[6])*surfVar_r)/dv_edge; 
  edgeSurf[7] = (((-1.875*fskin_over_jacv[7])-1.948557158514986*fskin_over_jacv[4])*surfVar_r)/dv_skin+((1.948557158514986*fedge_over_jacv[4]-1.875*fedge_over_jacv[7])*surfVar_r)/dv_edge; 
  edgeSurf[8] = (((-2.165063509461097*fskin_over_jacv[10])-2.25*fskin_over_jacv[8])*surfVar_r)/dv_skinR2+((2.25*fedge_over_jacv[8]-2.165063509461097*fedge_over_jacv[10])*surfVar_r)/dv_edgeR2; 
  edgeSurf[9] = (((-2.165063509461097*fskin_over_jacv[11])-2.25*fskin_over_jacv[9])*surfVar_r)/dv_skinR2+((2.25*fedge_over_jacv[9]-2.165063509461097*fedge_over_jacv[11])*surfVar_r)/dv_edgeR2; 
  edgeSurf[10] = (((-3.75*fskin_over_jacv[10])-3.897114317029974*fskin_over_jacv[8])*surfVar_r)/dv_skinR2+((3.897114317029974*fedge_over_jacv[8]-3.75*fedge_over_jacv[10])*surfVar_r)/dv_edgeR2; 
  edgeSurf[11] = (((-3.75*fskin_over_jacv[11])-3.897114317029974*fskin_over_jacv[9])*surfVar_r)/dv_skinR2+((3.897114317029974*fedge_over_jacv[9]-3.75*fedge_over_jacv[11])*surfVar_r)/dv_edgeR2; 

  edgeSurf_incr[0] = 0.7071067811865475*(edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.7071067811865475*(edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865475*(facDiff[1]*edgeSurf[4]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865475*(facDiff[1]*edgeSurf[5]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.7071067811865475*(facDiff[0]*edgeSurf[4]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.7071067811865475*(facDiff[0]*edgeSurf[5]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(facDiff[1]*edgeSurf[7]+facDiff[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.7071067811865475*(facDiff[0]*edgeSurf[7]+facDiff[1]*edgeSurf[6]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[9]+15.0*facDiff[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[9]+15.0*facDiff[1]*edgeSurf[8]); 
  edgeSurf_incr[10] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[11]+15.0*facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[11]+15.0*facDiff[1]*edgeSurf[10]); 


  } else { 

  double edgeSurf[12] = {0.0}; 
  const double dv_edgeR2 = pow(dv_edge,2);
  const double dv_skinR2 = pow(dv_skin,2);

  edgeSurf[0] = ((2.420614591379636*fskin_over_jacv[10]-2.515576474687264*fskin_over_jacv[8])*surfVar_l)/dv_skinR2+((2.420614591379636*fedge_over_jacv[10]+2.515576474687264*fedge_over_jacv[8])*surfVar_l)/dv_edgeR2+((-0.6051536478449089*fskin_over_jacv[10])-0.6051536478449089*fedge_over_jacv[10]+0.6288941186718159*fskin_over_jacv[8]-0.6288941186718159*fedge_over_jacv[8]+0.5412658773652741*fskin_over_jacv[3]+0.5412658773652741*fedge_over_jacv[3]-0.5625*fskin_over_jacv[0]+0.5625*fedge_over_jacv[0])*surfVar_l; 
  edgeSurf[1] = ((2.420614591379636*fskin_over_jacv[11]-2.515576474687263*fskin_over_jacv[9])*surfVar_l)/dv_skinR2+((2.420614591379636*fedge_over_jacv[11]+2.515576474687263*fedge_over_jacv[9])*surfVar_l)/dv_edgeR2+((-0.605153647844909*fskin_over_jacv[11])-0.605153647844909*fedge_over_jacv[11]+0.6288941186718158*fskin_over_jacv[9]-0.6288941186718158*fedge_over_jacv[9]+0.5412658773652741*fskin_over_jacv[5]+0.5412658773652741*fedge_over_jacv[5]-0.5625*fskin_over_jacv[1]+0.5625*fedge_over_jacv[1])*surfVar_l; 
  edgeSurf[2] = ((1.082531754730548*fskin_over_jacv[6]-1.125*fskin_over_jacv[2])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[6]+1.125*fedge_over_jacv[2])*surfVar_l)/dv_edge; 
  edgeSurf[3] = ((4.357106264483343*fskin_over_jacv[8]-4.192627457812105*fskin_over_jacv[10])*surfVar_l)/dv_skinR2+(((-4.192627457812105*fedge_over_jacv[10])-4.357106264483343*fedge_over_jacv[8])*surfVar_l)/dv_edgeR2+(1.048156864453026*fskin_over_jacv[10]+1.048156864453026*fedge_over_jacv[10]-1.089276566120836*fskin_over_jacv[8]+1.089276566120836*fedge_over_jacv[8]-0.9375*fskin_over_jacv[3]-0.9375*fedge_over_jacv[3]+0.9742785792574932*fskin_over_jacv[0]-0.9742785792574932*fedge_over_jacv[0])*surfVar_l; 
  edgeSurf[4] = ((1.082531754730548*fskin_over_jacv[7]-1.125*fskin_over_jacv[4])*surfVar_l)/dv_skin+((1.082531754730548*fedge_over_jacv[7]+1.125*fedge_over_jacv[4])*surfVar_l)/dv_edge; 
  edgeSurf[5] = ((4.357106264483344*fskin_over_jacv[9]-4.192627457812106*fskin_over_jacv[11])*surfVar_l)/dv_skinR2+(((-4.192627457812106*fedge_over_jacv[11])-4.357106264483344*fedge_over_jacv[9])*surfVar_l)/dv_edgeR2+(1.048156864453027*fskin_over_jacv[11]+1.048156864453027*fedge_over_jacv[11]-1.089276566120836*fskin_over_jacv[9]+1.089276566120836*fedge_over_jacv[9]-0.9375*fskin_over_jacv[5]-0.9375*fedge_over_jacv[5]+0.9742785792574932*fskin_over_jacv[1]-0.9742785792574932*fedge_over_jacv[1])*surfVar_l; 
  edgeSurf[6] = ((1.948557158514986*fskin_over_jacv[2]-1.875*fskin_over_jacv[6])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[6])-1.948557158514986*fedge_over_jacv[2])*surfVar_l)/dv_edge; 
  edgeSurf[7] = ((1.948557158514986*fskin_over_jacv[4]-1.875*fskin_over_jacv[7])*surfVar_l)/dv_skin+(((-1.875*fedge_over_jacv[7])-1.948557158514986*fedge_over_jacv[4])*surfVar_l)/dv_edge; 
  edgeSurf[8] = ((2.165063509461097*fskin_over_jacv[10]-2.25*fskin_over_jacv[8])*surfVar_l)/dv_skinR2+((2.165063509461097*fedge_over_jacv[10]+2.25*fedge_over_jacv[8])*surfVar_l)/dv_edgeR2; 
  edgeSurf[9] = ((2.165063509461097*fskin_over_jacv[11]-2.25*fskin_over_jacv[9])*surfVar_l)/dv_skinR2+((2.165063509461097*fedge_over_jacv[11]+2.25*fedge_over_jacv[9])*surfVar_l)/dv_edgeR2; 
  edgeSurf[10] = ((3.897114317029974*fskin_over_jacv[8]-3.75*fskin_over_jacv[10])*surfVar_l)/dv_skinR2+(((-3.75*fedge_over_jacv[10])-3.897114317029974*fedge_over_jacv[8])*surfVar_l)/dv_edgeR2; 
  edgeSurf[11] = ((3.897114317029974*fskin_over_jacv[9]-3.75*fskin_over_jacv[11])*surfVar_l)/dv_skinR2+(((-3.75*fedge_over_jacv[11])-3.897114317029974*fedge_over_jacv[9])*surfVar_l)/dv_edgeR2; 

  edgeSurf_incr[0] = 0.7071067811865475*(edgeSurf[1]*facDiff[1]+edgeSurf[0]*facDiff[0]); 
  edgeSurf_incr[1] = 0.7071067811865475*(edgeSurf[0]*facDiff[1]+facDiff[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865475*(facDiff[1]*edgeSurf[4]+facDiff[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865475*(facDiff[1]*edgeSurf[5]+facDiff[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.7071067811865475*(facDiff[0]*edgeSurf[4]+facDiff[1]*edgeSurf[2]); 
  edgeSurf_incr[5] = 0.7071067811865475*(facDiff[0]*edgeSurf[5]+facDiff[1]*edgeSurf[3]); 
  edgeSurf_incr[6] = 0.7071067811865475*(facDiff[1]*edgeSurf[7]+facDiff[0]*edgeSurf[6]); 
  edgeSurf_incr[7] = 0.7071067811865475*(facDiff[0]*edgeSurf[7]+facDiff[1]*edgeSurf[6]); 
  edgeSurf_incr[8] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[9]+15.0*facDiff[0]*edgeSurf[8]); 
  edgeSurf_incr[9] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[9]+15.0*facDiff[1]*edgeSurf[8]); 
  edgeSurf_incr[10] = 0.04714045207910316*(15.0*facDiff[1]*edgeSurf[11]+15.0*facDiff[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.04714045207910316*(15.0*facDiff[0]*edgeSurf[11]+15.0*facDiff[1]*edgeSurf[10]); 

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

  return 0.;

} 
