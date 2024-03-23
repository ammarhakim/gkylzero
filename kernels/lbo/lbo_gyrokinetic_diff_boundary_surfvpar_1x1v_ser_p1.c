#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_1x1v_ser_p1(const double *dxv, const double *vmap_edge, const double *vmap_skin, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // dxv[2]: Cell spacing. 
  // vmap_edge,vmap_skin: velocity space mapping.
  // vmap_prime: velocity space mapping derivative (in the skin cell).
  // jacobvel: velocity space jacobian.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[2*NC]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fskin/edge: Distribution function in cells 
  // out: Incremented distribution function in cell 

  double fedge_over_jacv[6], fskin_over_jacv[6];
  fedge_over_jacv[0] = fedge[0]/jacobvel[0]; 
  fedge_over_jacv[1] = fedge[1]/jacobvel[0]; 
  fedge_over_jacv[2] = fedge[2]/jacobvel[0]; 
  fedge_over_jacv[3] = fedge[3]/jacobvel[0]; 
  fedge_over_jacv[4] = fedge[4]/jacobvel[0]; 
  fedge_over_jacv[5] = fedge[5]/jacobvel[0]; 

  fskin_over_jacv[0] = fskin[0]/jacobvel[0]; 
  fskin_over_jacv[1] = fskin[1]/jacobvel[0]; 
  fskin_over_jacv[2] = fskin[2]/jacobvel[0]; 
  fskin_over_jacv[3] = fskin[3]/jacobvel[0]; 
  fskin_over_jacv[4] = fskin[4]/jacobvel[0]; 
  fskin_over_jacv[5] = fskin[5]/jacobvel[0]; 

  double dv_edge = 2.449489742783178*vmap_edge[1];
  double dv_skin = 2.449489742783178*vmap_skin[1];

  const double *nuVtSqSum = &nuPrimMomsSum[2];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double edgeSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[6] = {0.0}; 

  edgeSurf[0] = (-0.6708203932499369*fskin_over_jacv[4])+0.6708203932499369*fedge_over_jacv[4]-1.190784930203603*fskin_over_jacv[2]-1.190784930203603*fedge_over_jacv[2]-0.9375*fskin_over_jacv[0]+0.9375*fedge_over_jacv[0]; 
  edgeSurf[1] = ((-1.341640786499874*fskin_over_jacv[5])-2.381569860407206*fskin_over_jacv[3])/dv_skin-(1.875*fskin_over_jacv[1])/dv_skin+(1.341640786499874*fedge_over_jacv[5]-2.381569860407206*fedge_over_jacv[3])/dv_edge+(1.875*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[2] = (-1.161895003862225*fskin_over_jacv[4])+1.161895003862225*fedge_over_jacv[4]-2.0625*fskin_over_jacv[2]-2.0625*fedge_over_jacv[2]-1.623797632095822*fskin_over_jacv[0]+1.623797632095822*fedge_over_jacv[0]; 
  edgeSurf[3] = ((-2.32379000772445*fskin_over_jacv[5])-4.125*fskin_over_jacv[3])/dv_skin-(3.247595264191645*fskin_over_jacv[1])/dv_skin+(2.32379000772445*fedge_over_jacv[5]-4.125*fedge_over_jacv[3])/dv_edge+(3.247595264191645*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[4] = (-1.5*fskin_over_jacv[4])+1.5*fedge_over_jacv[4]-2.662676050517599*fskin_over_jacv[2]-2.662676050517599*fedge_over_jacv[2]-2.096313728906053*fskin_over_jacv[0]+2.096313728906053*fedge_over_jacv[0]; 
  edgeSurf[5] = ((-3.0*fskin_over_jacv[5])-5.325352101035199*fskin_over_jacv[3])/dv_skin-(4.192627457812105*fskin_over_jacv[1])/dv_skin+(3.0*fedge_over_jacv[5]-5.325352101035199*fedge_over_jacv[3])/dv_edge+(4.192627457812105*fedge_over_jacv[1])/dv_edge; 

  edgeSurf_incr[0] = 0.7071067811865476*(edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.7071067811865476*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865476*(nuVtSqSum[1]*edgeSurf[3]+nuVtSqSum[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865476*(nuVtSqSum[0]*edgeSurf[3]+nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[4] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[1]*edgeSurf[5]+21.21320343559643*nuVtSqSum[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[0]*edgeSurf[5]+21.21320343559643*nuVtSqSum[1]*edgeSurf[4]); 


  } else { 

  double edgeSurf[6] = {0.0}; 

  edgeSurf[0] = (-0.6708203932499369*fskin_over_jacv[4])+0.6708203932499369*fedge_over_jacv[4]+1.190784930203603*fskin_over_jacv[2]+1.190784930203603*fedge_over_jacv[2]-0.9375*fskin_over_jacv[0]+0.9375*fedge_over_jacv[0]; 
  edgeSurf[1] = (2.381569860407206*fskin_over_jacv[3]-1.341640786499874*fskin_over_jacv[5])/dv_skin-(1.875*fskin_over_jacv[1])/dv_skin+(1.341640786499874*fedge_over_jacv[5]+2.381569860407206*fedge_over_jacv[3])/dv_edge+(1.875*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[2] = 1.161895003862225*fskin_over_jacv[4]-1.161895003862225*fedge_over_jacv[4]-2.0625*fskin_over_jacv[2]-2.0625*fedge_over_jacv[2]+1.623797632095822*fskin_over_jacv[0]-1.623797632095822*fedge_over_jacv[0]; 
  edgeSurf[3] = (2.32379000772445*fskin_over_jacv[5]-4.125*fskin_over_jacv[3])/dv_skin+(3.247595264191645*fskin_over_jacv[1])/dv_skin+((-2.32379000772445*fedge_over_jacv[5])-4.125*fedge_over_jacv[3])/dv_edge-(3.247595264191645*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[4] = (-1.5*fskin_over_jacv[4])+1.5*fedge_over_jacv[4]+2.662676050517599*fskin_over_jacv[2]+2.662676050517599*fedge_over_jacv[2]-2.096313728906053*fskin_over_jacv[0]+2.096313728906053*fedge_over_jacv[0]; 
  edgeSurf[5] = (5.325352101035199*fskin_over_jacv[3]-3.0*fskin_over_jacv[5])/dv_skin-(4.192627457812105*fskin_over_jacv[1])/dv_skin+(3.0*fedge_over_jacv[5]+5.325352101035199*fedge_over_jacv[3])/dv_edge+(4.192627457812105*fedge_over_jacv[1])/dv_edge; 

  edgeSurf_incr[0] = 0.7071067811865476*(edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.7071067811865476*(edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.7071067811865476*(nuVtSqSum[1]*edgeSurf[3]+nuVtSqSum[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.7071067811865476*(nuVtSqSum[0]*edgeSurf[3]+nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[4] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[1]*edgeSurf[5]+21.21320343559643*nuVtSqSum[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.03333333333333333*(21.21320343559643*nuVtSqSum[0]*edgeSurf[5]+21.21320343559643*nuVtSqSum[1]*edgeSurf[4]); 

  } 

  out[0] += edgeSurf_incr[0]*rdvSq4; 
  out[1] += edgeSurf_incr[1]*rdvSq4; 
  out[2] += edgeSurf_incr[2]*rdvSq4; 
  out[3] += edgeSurf_incr[3]*rdvSq4; 
  out[4] += edgeSurf_incr[4]*rdvSq4; 
  out[5] += edgeSurf_incr[5]*rdvSq4; 

  return 0.;

} 
