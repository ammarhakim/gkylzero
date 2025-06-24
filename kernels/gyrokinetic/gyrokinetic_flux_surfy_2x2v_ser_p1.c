#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_flux_surfy_2x2v_ser_p1(
    const struct gkyl_basis *basis, const double *w, const double *dxv,
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
    const double *bmag, const double *phi, const double *JfL, const double *JfR, 
    double* GKYL_RESTRICT flux_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential.
  // bmag_surf: bmag represented on the surface.
  // jacobtot_inv_surf: jacobtot_inv represented on the surface.
  // cmag_surf: cmag represented on the surface.
  // b_i_surf: b_i represented on the surface.
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[24] = {0.}; 
  hamil[0] = (1.414213562373095*phi[0]-2.449489742783179*phi[2])*q_+vmapSq[0]*m_+1.414213562373095*bmag[0]*vmap[2]; 
  hamil[1] = 1.414213562373095*(phi[1]*q_+bmag[1]*vmap[2])-2.449489742783179*phi[3]*q_; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

double flux_surf_nodal[12]= {0.0}; 
double cfl = 0.0; 
double bmag_quad = 0.0; 
double Jc_quad = 0.0; 
double B3_quad = 0.0; 
double normcurlbhat_quad = 0.0; 
double area_elem_quad = 0.0; 
double bhat_quad[3] = {0.0}; 
double alpha_quad = 0.0; 
double JfL_quad = 0.0; 
double JfR_quad = 0.0; 
double Jfavg_quad = 0.0; 
double Jfjump_quad = 0.0; 



bmag_quad = gkdgs[0].bmag; 
Jc_quad = gkdgs[0].Jc; 
B3_quad = gkdgs[0].B3; 
normcurlbhat_quad = gkdgs[0].normcurlbhat; 
bhat_quad[0] = gkdgs[0].bhat.x[0];  
bhat_quad[1] = gkdgs[0].bhat.x[1];  
bhat_quad[2] = gkdgs[0].bhat.x[2];  
area_elem_quad = dgs[0].area_elem; 

alpha_quad = ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21])-0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])-0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])-0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])+0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])+0.25*(1.732050807568877*JfL[12]+JfL[8])+0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])-0.25*(1.732050807568877*JfL[9]+JfL[4])-0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])-0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  (-0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21]))+0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])+0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])+0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])-0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])-0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])-0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])+0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])+0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])+0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  (-0.05590169943749474*(8.660254037844388*JfL[23]+5.0*JfL[21]))+0.05590169943749476*(8.660254037844388*JfL[22]+5.0*JfL[19])+0.05590169943749476*(8.660254037844388*JfL[20]+5.0*JfL[17])-0.05590169943749474*(8.660254037844388*JfL[18]+5.0*JfL[16])+0.25*(1.732050807568877*JfL[12]+JfL[8])-0.25*(1.732050807568877*JfL[9]+JfL[4])-0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  0.05590169943749474*(8.660254037844388*JfR[23]-5.0*JfR[21])-0.05590169943749476*(8.660254037844388*JfR[22]-5.0*JfR[19])-0.05590169943749476*(8.660254037844388*JfR[20]-5.0*JfR[17])+0.05590169943749474*(8.660254037844388*JfR[18]-5.0*JfR[16])-0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])+0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])+0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * (1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21])-0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])-0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])+0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])-0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])+0.25*(1.732050807568877*JfL[12]+JfL[8])-0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])-0.25*(1.732050807568877*JfL[9]+JfL[4])+0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])-0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  (-0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21]))+0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])+0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])-0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])+0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])-0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])+0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])+0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])-0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])+0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  (-0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21]))+0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])-0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])+0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])-0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])-0.25*(1.732050807568877*JfL[12]+JfL[8])+0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])+0.25*(1.732050807568877*JfL[9]+JfL[4])-0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])-0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21])-0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])+0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])-0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])+0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])+0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])-0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])-0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])+0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])+0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.05590169943749474*(8.660254037844388*JfL[23]+5.0*JfL[21])-0.05590169943749476*(8.660254037844388*JfL[22]+5.0*JfL[19])+0.05590169943749476*(8.660254037844388*JfL[20]+5.0*JfL[17])-0.05590169943749474*(8.660254037844388*JfL[18]+5.0*JfL[16])-0.25*(1.732050807568877*JfL[12]+JfL[8])+0.25*(1.732050807568877*JfL[9]+JfL[4])-0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  (-0.05590169943749474*(8.660254037844388*JfR[23]-5.0*JfR[21]))+0.05590169943749476*(8.660254037844388*JfR[22]-5.0*JfR[19])-0.05590169943749476*(8.660254037844388*JfR[20]-5.0*JfR[17])+0.05590169943749474*(8.660254037844388*JfR[18]-5.0*JfR[16])+0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])-0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])+0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * (1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  (-0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21]))+0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])-0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])-0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])+0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])-0.25*(1.732050807568877*JfL[12]+JfL[8])-0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])+0.25*(1.732050807568877*JfL[9]+JfL[4])+0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])-0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21])-0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])+0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])+0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])-0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])+0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])+0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])-0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])-0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])+0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;


bmag_quad = gkdgs[1].bmag; 
Jc_quad = gkdgs[1].Jc; 
B3_quad = gkdgs[1].B3; 
normcurlbhat_quad = gkdgs[1].normcurlbhat; 
bhat_quad[0] = gkdgs[1].bhat.x[0];  
bhat_quad[1] = gkdgs[1].bhat.x[1];  
bhat_quad[2] = gkdgs[1].bhat.x[2];  
area_elem_quad = dgs[1].area_elem; 

alpha_quad = ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  (-0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21]))-0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])+0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])+0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])+0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])-0.25*(1.732050807568877*JfL[12]+JfL[8])-0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])-0.25*(1.732050807568877*JfL[9]+JfL[4])-0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])+0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21])+0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])-0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])-0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])-0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])+0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])+0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])+0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])+0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])-0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[6] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.05590169943749474*(8.660254037844388*JfL[23]+5.0*JfL[21])+0.05590169943749476*(8.660254037844388*JfL[22]+5.0*JfL[19])-0.05590169943749476*(8.660254037844388*JfL[20]+5.0*JfL[17])-0.05590169943749474*(8.660254037844388*JfL[18]+5.0*JfL[16])-0.25*(1.732050807568877*JfL[12]+JfL[8])-0.25*(1.732050807568877*JfL[9]+JfL[4])+0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  (-0.05590169943749474*(8.660254037844388*JfR[23]-5.0*JfR[21]))-0.05590169943749476*(8.660254037844388*JfR[22]-5.0*JfR[19])+0.05590169943749476*(8.660254037844388*JfR[20]-5.0*JfR[17])+0.05590169943749474*(8.660254037844388*JfR[18]-5.0*JfR[16])+0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])+0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])-0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[7] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * (1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  (-0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21]))-0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])+0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])-0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])-0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])-0.25*(1.732050807568877*JfL[12]+JfL[8])+0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])-0.25*(1.732050807568877*JfL[9]+JfL[4])+0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])+0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21])+0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])-0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])+0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])+0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])+0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])-0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])+0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])-0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])-0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[8] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2]-1.837117307087383*hamil[8])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21])+0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])+0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])-0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])-0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])+0.25*(1.732050807568877*JfL[12]+JfL[8])-0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])+0.25*(1.732050807568877*JfL[9]+JfL[4])-0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])+0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  (-0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21]))-0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])-0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])+0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])+0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])-0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])+0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])-0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])+0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])-0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[9] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.5*hamil[2])/vmap[1])/m_/bmag_quad * (0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  (-0.05590169943749474*(8.660254037844388*JfL[23]+5.0*JfL[21]))-0.05590169943749476*(8.660254037844388*JfL[22]+5.0*JfL[19])-0.05590169943749476*(8.660254037844388*JfL[20]+5.0*JfL[17])-0.05590169943749474*(8.660254037844388*JfL[18]+5.0*JfL[16])+0.25*(1.732050807568877*JfL[12]+JfL[8])+0.25*(1.732050807568877*JfL[9]+JfL[4])+0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  0.05590169943749474*(8.660254037844388*JfR[23]-5.0*JfR[21])+0.05590169943749476*(8.660254037844388*JfR[22]-5.0*JfR[19])+0.05590169943749476*(8.660254037844388*JfR[20]-5.0*JfR[17])+0.05590169943749474*(8.660254037844388*JfR[18]-5.0*JfR[16])-0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])-0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])-0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[10] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1])/m_/bmag_quad * (1.837117307087383*hamil[8]+0.6123724356957944*hamil[2])/q_ * normcurlbhat_quad ;
alpha_quad += 1/q_/bmag_quad/area_elem_quad * -bhat_quad[0]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdx2) ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.04472135954999579*(8.660254037844388*JfL[23]+5.0*JfL[21])+0.04472135954999579*(8.660254037844388*JfL[22]+5.0*JfL[19])+0.04472135954999581*(8.660254037844388*JfL[20]+5.0*JfL[17])+0.04472135954999579*(8.660254037844388*JfL[18]+5.0*JfL[16])+0.3354101966249685*(1.732050807568877*JfL[15]+JfL[13])+0.3354101966249685*(1.732050807568877*JfL[14]+JfL[10])+0.25*(1.732050807568877*JfL[12]+JfL[8])+0.3354101966249685*(1.732050807568877*JfL[11]+JfL[6])+0.25*(1.732050807568877*JfL[9]+JfL[4])+0.3354101966249685*(1.732050807568877*JfL[7]+JfL[3])+0.25*(1.732050807568877*JfL[5]+JfL[1])+0.25*(1.732050807568877*JfL[2]+JfL[0]);
JfR_quad =  (-0.04472135954999579*(8.660254037844388*JfR[23]-5.0*JfR[21]))-0.04472135954999579*(8.660254037844388*JfR[22]-5.0*JfR[19])-0.04472135954999581*(8.660254037844388*JfR[20]-5.0*JfR[17])-0.04472135954999579*(8.660254037844388*JfR[18]-5.0*JfR[16])-0.3354101966249685*(1.732050807568877*JfR[15]-1.0*JfR[13])-0.3354101966249685*(1.732050807568877*JfR[14]-1.0*JfR[10])-0.25*(1.732050807568877*JfR[12]-1.0*JfR[8])-0.3354101966249685*(1.732050807568877*JfR[11]-1.0*JfR[6])-0.25*(1.732050807568877*JfR[9]-1.0*JfR[4])-0.3354101966249685*(1.732050807568877*JfR[7]-1.0*JfR[3])-0.25*(1.732050807568877*JfR[5]-1.0*JfR[1])-0.25*(1.732050807568877*JfR[2]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[11] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 0) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 1) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 2) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 3) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 4) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 5) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 6) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 7) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 8) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 9) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 10) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 11) ;


  return cfl*1.5*rdz2; 

} 
