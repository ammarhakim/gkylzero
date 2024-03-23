#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_2x2v_ser_p1(const double *dxv, const double *vmap_edge, const double *vmap_skin, const double *vmap_prime, const double *jacobvel, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  double dv_edge = 2.449489742783178*vmap_edge[1];
  double dv_skin = 2.449489742783178*vmap_skin[1];

  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  double edgeSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[24] = {0.0}; 

  edgeSurf[0] = (-0.6708203932499369*fskin_over_jacv[16])+0.6708203932499369*fedge_over_jacv[16]-1.190784930203603*fskin_over_jacv[3]-1.190784930203603*fedge_over_jacv[3]-0.9375*fskin_over_jacv[0]+0.9375*fedge_over_jacv[0]; 
  edgeSurf[1] = ((-1.341640786499874*fskin_over_jacv[17])-2.381569860407206*fskin_over_jacv[6])/dv_skin-(1.875*fskin_over_jacv[1])/dv_skin+(1.341640786499874*fedge_over_jacv[17]-2.381569860407206*fedge_over_jacv[6])/dv_edge+(1.875*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[2] = (-0.6708203932499369*fskin_over_jacv[18])+0.6708203932499369*fedge_over_jacv[18]-1.190784930203603*fskin_over_jacv[7]-1.190784930203603*fedge_over_jacv[7]-0.9375*fskin_over_jacv[2]+0.9375*fedge_over_jacv[2]; 
  edgeSurf[3] = (-1.161895003862225*fskin_over_jacv[16])+1.161895003862225*fedge_over_jacv[16]-2.0625*fskin_over_jacv[3]-2.0625*fedge_over_jacv[3]-1.623797632095822*fskin_over_jacv[0]+1.623797632095822*fedge_over_jacv[0]; 
  edgeSurf[4] = (-0.6708203932499369*fskin_over_jacv[19])+0.6708203932499369*fedge_over_jacv[19]-1.190784930203603*fskin_over_jacv[10]-1.190784930203603*fedge_over_jacv[10]-0.9375*fskin_over_jacv[4]+0.9375*fedge_over_jacv[4]; 
  edgeSurf[5] = ((-1.341640786499874*fskin_over_jacv[20])-2.381569860407206*fskin_over_jacv[11]-1.875*fskin_over_jacv[5])/dv_skin+(1.341640786499874*fedge_over_jacv[20]-2.381569860407206*fedge_over_jacv[11]+1.875*fedge_over_jacv[5])/dv_edge; 
  edgeSurf[6] = ((-2.32379000772445*fskin_over_jacv[17])-4.125*fskin_over_jacv[6])/dv_skin-(3.247595264191645*fskin_over_jacv[1])/dv_skin+(2.32379000772445*fedge_over_jacv[17]-4.125*fedge_over_jacv[6])/dv_edge+(3.247595264191645*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[7] = (-1.161895003862225*fskin_over_jacv[18])+1.161895003862225*fedge_over_jacv[18]-2.0625*fskin_over_jacv[7]-2.0625*fedge_over_jacv[7]-1.623797632095822*fskin_over_jacv[2]+1.623797632095822*fedge_over_jacv[2]; 
  edgeSurf[8] = ((-1.341640786499874*fskin_over_jacv[21])-2.381569860407206*fskin_over_jacv[13]-1.875*fskin_over_jacv[8])/dv_skin+(1.341640786499874*fedge_over_jacv[21]-2.381569860407206*fedge_over_jacv[13]+1.875*fedge_over_jacv[8])/dv_edge; 
  edgeSurf[9] = (-0.6708203932499369*fskin_over_jacv[22])+0.6708203932499369*fedge_over_jacv[22]-1.190784930203603*fskin_over_jacv[14]-1.190784930203603*fedge_over_jacv[14]-0.9375*fskin_over_jacv[9]+0.9375*fedge_over_jacv[9]; 
  edgeSurf[10] = (-1.161895003862225*fskin_over_jacv[19])+1.161895003862225*fedge_over_jacv[19]-2.0625*fskin_over_jacv[10]-2.0625*fedge_over_jacv[10]-1.623797632095822*fskin_over_jacv[4]+1.623797632095822*fedge_over_jacv[4]; 
  edgeSurf[11] = ((-2.32379000772445*fskin_over_jacv[20])-4.125*fskin_over_jacv[11]-3.247595264191645*fskin_over_jacv[5])/dv_skin+(2.32379000772445*fedge_over_jacv[20]-4.125*fedge_over_jacv[11]+3.247595264191645*fedge_over_jacv[5])/dv_edge; 
  edgeSurf[12] = ((-1.341640786499874*fskin_over_jacv[23])-2.381569860407206*fskin_over_jacv[15]-1.875*fskin_over_jacv[12])/dv_skin+(1.341640786499874*fedge_over_jacv[23]-2.381569860407206*fedge_over_jacv[15]+1.875*fedge_over_jacv[12])/dv_edge; 
  edgeSurf[13] = ((-2.32379000772445*fskin_over_jacv[21])-4.125*fskin_over_jacv[13]-3.247595264191645*fskin_over_jacv[8])/dv_skin+(2.32379000772445*fedge_over_jacv[21]-4.125*fedge_over_jacv[13]+3.247595264191645*fedge_over_jacv[8])/dv_edge; 
  edgeSurf[14] = (-1.161895003862225*fskin_over_jacv[22])+1.161895003862225*fedge_over_jacv[22]-2.0625*fskin_over_jacv[14]-2.0625*fedge_over_jacv[14]-1.623797632095822*fskin_over_jacv[9]+1.623797632095822*fedge_over_jacv[9]; 
  edgeSurf[15] = ((-2.32379000772445*fskin_over_jacv[23])-4.125*fskin_over_jacv[15]-3.247595264191645*fskin_over_jacv[12])/dv_skin+(2.32379000772445*fedge_over_jacv[23]-4.125*fedge_over_jacv[15]+3.247595264191645*fedge_over_jacv[12])/dv_edge; 
  edgeSurf[16] = (-1.5*fskin_over_jacv[16])+1.5*fedge_over_jacv[16]-2.662676050517599*fskin_over_jacv[3]-2.662676050517599*fedge_over_jacv[3]-2.096313728906053*fskin_over_jacv[0]+2.096313728906053*fedge_over_jacv[0]; 
  edgeSurf[17] = ((-3.0*fskin_over_jacv[17])-5.325352101035199*fskin_over_jacv[6])/dv_skin-(4.192627457812105*fskin_over_jacv[1])/dv_skin+(3.0*fedge_over_jacv[17]-5.325352101035199*fedge_over_jacv[6])/dv_edge+(4.192627457812105*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[18] = (-1.5*fskin_over_jacv[18])+1.5*fedge_over_jacv[18]-2.662676050517599*fskin_over_jacv[7]-2.662676050517599*fedge_over_jacv[7]-2.096313728906053*fskin_over_jacv[2]+2.096313728906053*fedge_over_jacv[2]; 
  edgeSurf[19] = (-1.5*fskin_over_jacv[19])+1.5*fedge_over_jacv[19]-2.662676050517599*fskin_over_jacv[10]-2.662676050517599*fedge_over_jacv[10]-2.096313728906053*fskin_over_jacv[4]+2.096313728906053*fedge_over_jacv[4]; 
  edgeSurf[20] = ((-3.0*fskin_over_jacv[20])-5.325352101035199*fskin_over_jacv[11]-4.192627457812106*fskin_over_jacv[5])/dv_skin+(3.0*fedge_over_jacv[20]-5.325352101035199*fedge_over_jacv[11]+4.192627457812106*fedge_over_jacv[5])/dv_edge; 
  edgeSurf[21] = ((-3.0*fskin_over_jacv[21])-5.325352101035199*fskin_over_jacv[13]-4.192627457812106*fskin_over_jacv[8])/dv_skin+(3.0*fedge_over_jacv[21]-5.325352101035199*fedge_over_jacv[13]+4.192627457812106*fedge_over_jacv[8])/dv_edge; 
  edgeSurf[22] = (-1.5*fskin_over_jacv[22])+1.5*fedge_over_jacv[22]-2.662676050517599*fskin_over_jacv[14]-2.662676050517599*fedge_over_jacv[14]-2.096313728906053*fskin_over_jacv[9]+2.096313728906053*fedge_over_jacv[9]; 
  edgeSurf[23] = ((-3.0*fskin_over_jacv[23])-5.325352101035199*fskin_over_jacv[15]-4.192627457812105*fskin_over_jacv[12])/dv_skin+(3.0*fedge_over_jacv[23]-5.325352101035199*fedge_over_jacv[15]+4.192627457812105*fedge_over_jacv[12])/dv_edge; 

  edgeSurf_incr[0] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[5]+edgeSurf[2]*nuVtSqSum[2]+edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[5]+edgeSurf[2]*nuVtSqSum[3]+edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[5]+edgeSurf[1]*nuVtSqSum[3]+edgeSurf[0]*nuVtSqSum[2]+nuVtSqSum[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[11]+nuVtSqSum[2]*edgeSurf[7]+nuVtSqSum[1]*edgeSurf[6]+nuVtSqSum[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[12]+nuVtSqSum[2]*edgeSurf[9]+nuVtSqSum[1]*edgeSurf[8]+nuVtSqSum[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[5]+edgeSurf[0]*nuVtSqSum[3]+edgeSurf[1]*nuVtSqSum[2]+nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[6] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[11]+nuVtSqSum[3]*edgeSurf[7]+nuVtSqSum[0]*edgeSurf[6]+nuVtSqSum[1]*edgeSurf[3]); 
  edgeSurf_incr[7] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[11]+nuVtSqSum[0]*edgeSurf[7]+nuVtSqSum[3]*edgeSurf[6]+nuVtSqSum[2]*edgeSurf[3]); 
  edgeSurf_incr[8] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[12]+nuVtSqSum[3]*edgeSurf[9]+nuVtSqSum[0]*edgeSurf[8]+nuVtSqSum[1]*edgeSurf[4]); 
  edgeSurf_incr[9] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[12]+nuVtSqSum[0]*edgeSurf[9]+nuVtSqSum[3]*edgeSurf[8]+nuVtSqSum[2]*edgeSurf[4]); 
  edgeSurf_incr[10] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[15]+nuVtSqSum[2]*edgeSurf[14]+nuVtSqSum[1]*edgeSurf[13]+nuVtSqSum[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[11]+nuVtSqSum[1]*edgeSurf[7]+nuVtSqSum[2]*edgeSurf[6]+edgeSurf[3]*nuVtSqSum[3]); 
  edgeSurf_incr[12] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[12]+nuVtSqSum[1]*edgeSurf[9]+nuVtSqSum[2]*edgeSurf[8]+nuVtSqSum[3]*edgeSurf[4]); 
  edgeSurf_incr[13] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[15]+nuVtSqSum[3]*edgeSurf[14]+nuVtSqSum[0]*edgeSurf[13]+nuVtSqSum[1]*edgeSurf[10]); 
  edgeSurf_incr[14] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[15]+nuVtSqSum[0]*edgeSurf[14]+nuVtSqSum[3]*edgeSurf[13]+nuVtSqSum[2]*edgeSurf[10]); 
  edgeSurf_incr[15] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[15]+nuVtSqSum[1]*edgeSurf[14]+nuVtSqSum[2]*edgeSurf[13]+nuVtSqSum[3]*edgeSurf[10]); 
  edgeSurf_incr[16] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[3]*edgeSurf[20]+15.0*(nuVtSqSum[2]*edgeSurf[18]+nuVtSqSum[1]*edgeSurf[17])+15.0*nuVtSqSum[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[2]*edgeSurf[20]+15.0*(nuVtSqSum[3]*edgeSurf[18]+nuVtSqSum[0]*edgeSurf[17])+15.0*nuVtSqSum[1]*edgeSurf[16]); 
  edgeSurf_incr[18] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[1]*edgeSurf[20]+15.0*(nuVtSqSum[0]*edgeSurf[18]+nuVtSqSum[3]*edgeSurf[17])+15.0*nuVtSqSum[2]*edgeSurf[16]); 
  edgeSurf_incr[19] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[3]*edgeSurf[23]+15.0*(nuVtSqSum[2]*edgeSurf[22]+nuVtSqSum[1]*edgeSurf[21])+15.0*nuVtSqSum[0]*edgeSurf[19]); 
  edgeSurf_incr[20] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[0]*edgeSurf[20]+15.0*(nuVtSqSum[1]*edgeSurf[18]+nuVtSqSum[2]*edgeSurf[17])+15.0*nuVtSqSum[3]*edgeSurf[16]); 
  edgeSurf_incr[21] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[2]*edgeSurf[23]+15.0*(nuVtSqSum[3]*edgeSurf[22]+nuVtSqSum[0]*edgeSurf[21])+15.0*nuVtSqSum[1]*edgeSurf[19]); 
  edgeSurf_incr[22] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[1]*edgeSurf[23]+15.0*(nuVtSqSum[0]*edgeSurf[22]+nuVtSqSum[3]*edgeSurf[21])+15.0*nuVtSqSum[2]*edgeSurf[19]); 
  edgeSurf_incr[23] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[0]*edgeSurf[23]+15.0*(nuVtSqSum[1]*edgeSurf[22]+nuVtSqSum[2]*edgeSurf[21])+15.0*nuVtSqSum[3]*edgeSurf[19]); 


  } else { 

  double edgeSurf[24] = {0.0}; 

  edgeSurf[0] = (-0.6708203932499369*fskin_over_jacv[16])+0.6708203932499369*fedge_over_jacv[16]+1.190784930203603*fskin_over_jacv[3]+1.190784930203603*fedge_over_jacv[3]-0.9375*fskin_over_jacv[0]+0.9375*fedge_over_jacv[0]; 
  edgeSurf[1] = (2.381569860407206*fskin_over_jacv[6]-1.341640786499874*fskin_over_jacv[17])/dv_skin-(1.875*fskin_over_jacv[1])/dv_skin+(1.341640786499874*fedge_over_jacv[17]+2.381569860407206*fedge_over_jacv[6])/dv_edge+(1.875*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[2] = (-0.6708203932499369*fskin_over_jacv[18])+0.6708203932499369*fedge_over_jacv[18]+1.190784930203603*fskin_over_jacv[7]+1.190784930203603*fedge_over_jacv[7]-0.9375*fskin_over_jacv[2]+0.9375*fedge_over_jacv[2]; 
  edgeSurf[3] = 1.161895003862225*fskin_over_jacv[16]-1.161895003862225*fedge_over_jacv[16]-2.0625*fskin_over_jacv[3]-2.0625*fedge_over_jacv[3]+1.623797632095822*fskin_over_jacv[0]-1.623797632095822*fedge_over_jacv[0]; 
  edgeSurf[4] = (-0.6708203932499369*fskin_over_jacv[19])+0.6708203932499369*fedge_over_jacv[19]+1.190784930203603*fskin_over_jacv[10]+1.190784930203603*fedge_over_jacv[10]-0.9375*fskin_over_jacv[4]+0.9375*fedge_over_jacv[4]; 
  edgeSurf[5] = ((-1.341640786499874*fskin_over_jacv[20])+2.381569860407206*fskin_over_jacv[11]-1.875*fskin_over_jacv[5])/dv_skin+(1.341640786499874*fedge_over_jacv[20]+2.381569860407206*fedge_over_jacv[11]+1.875*fedge_over_jacv[5])/dv_edge; 
  edgeSurf[6] = (2.32379000772445*fskin_over_jacv[17]-4.125*fskin_over_jacv[6])/dv_skin+(3.247595264191645*fskin_over_jacv[1])/dv_skin+((-2.32379000772445*fedge_over_jacv[17])-4.125*fedge_over_jacv[6])/dv_edge-(3.247595264191645*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[7] = 1.161895003862225*fskin_over_jacv[18]-1.161895003862225*fedge_over_jacv[18]-2.0625*fskin_over_jacv[7]-2.0625*fedge_over_jacv[7]+1.623797632095822*fskin_over_jacv[2]-1.623797632095822*fedge_over_jacv[2]; 
  edgeSurf[8] = ((-1.341640786499874*fskin_over_jacv[21])+2.381569860407206*fskin_over_jacv[13]-1.875*fskin_over_jacv[8])/dv_skin+(1.341640786499874*fedge_over_jacv[21]+2.381569860407206*fedge_over_jacv[13]+1.875*fedge_over_jacv[8])/dv_edge; 
  edgeSurf[9] = (-0.6708203932499369*fskin_over_jacv[22])+0.6708203932499369*fedge_over_jacv[22]+1.190784930203603*fskin_over_jacv[14]+1.190784930203603*fedge_over_jacv[14]-0.9375*fskin_over_jacv[9]+0.9375*fedge_over_jacv[9]; 
  edgeSurf[10] = 1.161895003862225*fskin_over_jacv[19]-1.161895003862225*fedge_over_jacv[19]-2.0625*fskin_over_jacv[10]-2.0625*fedge_over_jacv[10]+1.623797632095822*fskin_over_jacv[4]-1.623797632095822*fedge_over_jacv[4]; 
  edgeSurf[11] = (2.32379000772445*fskin_over_jacv[20]-4.125*fskin_over_jacv[11]+3.247595264191645*fskin_over_jacv[5])/dv_skin+((-2.32379000772445*fedge_over_jacv[20])-4.125*fedge_over_jacv[11]-3.247595264191645*fedge_over_jacv[5])/dv_edge; 
  edgeSurf[12] = ((-1.341640786499874*fskin_over_jacv[23])+2.381569860407206*fskin_over_jacv[15]-1.875*fskin_over_jacv[12])/dv_skin+(1.341640786499874*fedge_over_jacv[23]+2.381569860407206*fedge_over_jacv[15]+1.875*fedge_over_jacv[12])/dv_edge; 
  edgeSurf[13] = (2.32379000772445*fskin_over_jacv[21]-4.125*fskin_over_jacv[13]+3.247595264191645*fskin_over_jacv[8])/dv_skin+((-2.32379000772445*fedge_over_jacv[21])-4.125*fedge_over_jacv[13]-3.247595264191645*fedge_over_jacv[8])/dv_edge; 
  edgeSurf[14] = 1.161895003862225*fskin_over_jacv[22]-1.161895003862225*fedge_over_jacv[22]-2.0625*fskin_over_jacv[14]-2.0625*fedge_over_jacv[14]+1.623797632095822*fskin_over_jacv[9]-1.623797632095822*fedge_over_jacv[9]; 
  edgeSurf[15] = (2.32379000772445*fskin_over_jacv[23]-4.125*fskin_over_jacv[15]+3.247595264191645*fskin_over_jacv[12])/dv_skin+((-2.32379000772445*fedge_over_jacv[23])-4.125*fedge_over_jacv[15]-3.247595264191645*fedge_over_jacv[12])/dv_edge; 
  edgeSurf[16] = (-1.5*fskin_over_jacv[16])+1.5*fedge_over_jacv[16]+2.662676050517599*fskin_over_jacv[3]+2.662676050517599*fedge_over_jacv[3]-2.096313728906053*fskin_over_jacv[0]+2.096313728906053*fedge_over_jacv[0]; 
  edgeSurf[17] = (5.325352101035199*fskin_over_jacv[6]-3.0*fskin_over_jacv[17])/dv_skin-(4.192627457812105*fskin_over_jacv[1])/dv_skin+(3.0*fedge_over_jacv[17]+5.325352101035199*fedge_over_jacv[6])/dv_edge+(4.192627457812105*fedge_over_jacv[1])/dv_edge; 
  edgeSurf[18] = (-1.5*fskin_over_jacv[18])+1.5*fedge_over_jacv[18]+2.662676050517599*fskin_over_jacv[7]+2.662676050517599*fedge_over_jacv[7]-2.096313728906053*fskin_over_jacv[2]+2.096313728906053*fedge_over_jacv[2]; 
  edgeSurf[19] = (-1.5*fskin_over_jacv[19])+1.5*fedge_over_jacv[19]+2.662676050517599*fskin_over_jacv[10]+2.662676050517599*fedge_over_jacv[10]-2.096313728906053*fskin_over_jacv[4]+2.096313728906053*fedge_over_jacv[4]; 
  edgeSurf[20] = ((-3.0*fskin_over_jacv[20])+5.325352101035199*fskin_over_jacv[11]-4.192627457812106*fskin_over_jacv[5])/dv_skin+(3.0*fedge_over_jacv[20]+5.325352101035199*fedge_over_jacv[11]+4.192627457812106*fedge_over_jacv[5])/dv_edge; 
  edgeSurf[21] = ((-3.0*fskin_over_jacv[21])+5.325352101035199*fskin_over_jacv[13]-4.192627457812106*fskin_over_jacv[8])/dv_skin+(3.0*fedge_over_jacv[21]+5.325352101035199*fedge_over_jacv[13]+4.192627457812106*fedge_over_jacv[8])/dv_edge; 
  edgeSurf[22] = (-1.5*fskin_over_jacv[22])+1.5*fedge_over_jacv[22]+2.662676050517599*fskin_over_jacv[14]+2.662676050517599*fedge_over_jacv[14]-2.096313728906053*fskin_over_jacv[9]+2.096313728906053*fedge_over_jacv[9]; 
  edgeSurf[23] = ((-3.0*fskin_over_jacv[23])+5.325352101035199*fskin_over_jacv[15]-4.192627457812105*fskin_over_jacv[12])/dv_skin+(3.0*fedge_over_jacv[23]+5.325352101035199*fedge_over_jacv[15]+4.192627457812105*fedge_over_jacv[12])/dv_edge; 

  edgeSurf_incr[0] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[5]+edgeSurf[2]*nuVtSqSum[2]+edgeSurf[1]*nuVtSqSum[1]+edgeSurf[0]*nuVtSqSum[0]); 
  edgeSurf_incr[1] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[5]+edgeSurf[2]*nuVtSqSum[3]+edgeSurf[0]*nuVtSqSum[1]+nuVtSqSum[0]*edgeSurf[1]); 
  edgeSurf_incr[2] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[5]+edgeSurf[1]*nuVtSqSum[3]+edgeSurf[0]*nuVtSqSum[2]+nuVtSqSum[0]*edgeSurf[2]); 
  edgeSurf_incr[3] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[11]+nuVtSqSum[2]*edgeSurf[7]+nuVtSqSum[1]*edgeSurf[6]+nuVtSqSum[0]*edgeSurf[3]); 
  edgeSurf_incr[4] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[12]+nuVtSqSum[2]*edgeSurf[9]+nuVtSqSum[1]*edgeSurf[8]+nuVtSqSum[0]*edgeSurf[4]); 
  edgeSurf_incr[5] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[5]+edgeSurf[0]*nuVtSqSum[3]+edgeSurf[1]*nuVtSqSum[2]+nuVtSqSum[1]*edgeSurf[2]); 
  edgeSurf_incr[6] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[11]+nuVtSqSum[3]*edgeSurf[7]+nuVtSqSum[0]*edgeSurf[6]+nuVtSqSum[1]*edgeSurf[3]); 
  edgeSurf_incr[7] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[11]+nuVtSqSum[0]*edgeSurf[7]+nuVtSqSum[3]*edgeSurf[6]+nuVtSqSum[2]*edgeSurf[3]); 
  edgeSurf_incr[8] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[12]+nuVtSqSum[3]*edgeSurf[9]+nuVtSqSum[0]*edgeSurf[8]+nuVtSqSum[1]*edgeSurf[4]); 
  edgeSurf_incr[9] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[12]+nuVtSqSum[0]*edgeSurf[9]+nuVtSqSum[3]*edgeSurf[8]+nuVtSqSum[2]*edgeSurf[4]); 
  edgeSurf_incr[10] = 0.5*vmap_prime[1]*(nuVtSqSum[3]*edgeSurf[15]+nuVtSqSum[2]*edgeSurf[14]+nuVtSqSum[1]*edgeSurf[13]+nuVtSqSum[0]*edgeSurf[10]); 
  edgeSurf_incr[11] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[11]+nuVtSqSum[1]*edgeSurf[7]+nuVtSqSum[2]*edgeSurf[6]+edgeSurf[3]*nuVtSqSum[3]); 
  edgeSurf_incr[12] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[12]+nuVtSqSum[1]*edgeSurf[9]+nuVtSqSum[2]*edgeSurf[8]+nuVtSqSum[3]*edgeSurf[4]); 
  edgeSurf_incr[13] = 0.5*vmap_prime[1]*(nuVtSqSum[2]*edgeSurf[15]+nuVtSqSum[3]*edgeSurf[14]+nuVtSqSum[0]*edgeSurf[13]+nuVtSqSum[1]*edgeSurf[10]); 
  edgeSurf_incr[14] = 0.5*vmap_prime[1]*(nuVtSqSum[1]*edgeSurf[15]+nuVtSqSum[0]*edgeSurf[14]+nuVtSqSum[3]*edgeSurf[13]+nuVtSqSum[2]*edgeSurf[10]); 
  edgeSurf_incr[15] = 0.5*vmap_prime[1]*(nuVtSqSum[0]*edgeSurf[15]+nuVtSqSum[1]*edgeSurf[14]+nuVtSqSum[2]*edgeSurf[13]+nuVtSqSum[3]*edgeSurf[10]); 
  edgeSurf_incr[16] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[3]*edgeSurf[20]+15.0*(nuVtSqSum[2]*edgeSurf[18]+nuVtSqSum[1]*edgeSurf[17])+15.0*nuVtSqSum[0]*edgeSurf[16]); 
  edgeSurf_incr[17] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[2]*edgeSurf[20]+15.0*(nuVtSqSum[3]*edgeSurf[18]+nuVtSqSum[0]*edgeSurf[17])+15.0*nuVtSqSum[1]*edgeSurf[16]); 
  edgeSurf_incr[18] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[1]*edgeSurf[20]+15.0*(nuVtSqSum[0]*edgeSurf[18]+nuVtSqSum[3]*edgeSurf[17])+15.0*nuVtSqSum[2]*edgeSurf[16]); 
  edgeSurf_incr[19] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[3]*edgeSurf[23]+15.0*(nuVtSqSum[2]*edgeSurf[22]+nuVtSqSum[1]*edgeSurf[21])+15.0*nuVtSqSum[0]*edgeSurf[19]); 
  edgeSurf_incr[20] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[0]*edgeSurf[20]+15.0*(nuVtSqSum[1]*edgeSurf[18]+nuVtSqSum[2]*edgeSurf[17])+15.0*nuVtSqSum[3]*edgeSurf[16]); 
  edgeSurf_incr[21] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[2]*edgeSurf[23]+15.0*(nuVtSqSum[3]*edgeSurf[22]+nuVtSqSum[0]*edgeSurf[21])+15.0*nuVtSqSum[1]*edgeSurf[19]); 
  edgeSurf_incr[22] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[1]*edgeSurf[23]+15.0*(nuVtSqSum[0]*edgeSurf[22]+nuVtSqSum[3]*edgeSurf[21])+15.0*nuVtSqSum[2]*edgeSurf[19]); 
  edgeSurf_incr[23] = 0.03333333333333333*vmap_prime[1]*(15.0*nuVtSqSum[0]*edgeSurf[23]+15.0*(nuVtSqSum[1]*edgeSurf[22]+nuVtSqSum[2]*edgeSurf[21])+15.0*nuVtSqSum[3]*edgeSurf[19]); 

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
