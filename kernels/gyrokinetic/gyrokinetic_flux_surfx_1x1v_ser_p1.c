#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_flux_surfx_1x1v_ser_p1(
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
  double rdvpar2 = 2.0/dxv[1];

  double hamil[6] = {0.}; 
  hamil[0] = (phi[0]-1.732050807568877*phi[1])*q_+0.5*vmapSq[0]*m_; 
  hamil[1] = 0.5*vmapSq[1]*m_; 
  hamil[2] = 0.5*vmapSq[2]*m_; 

double flux_surf_nodal[3]= {0.0}; 
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

alpha_quad = ((0.8164965809277262*(1.224744871391589*hamil[1]-3.674234614174767*hamil[2]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(1.224744871391589*hamil[1]-3.674234614174767*hamil[2]))/vmap[1])/m_/bmag_quad * (1.224744871391589*hamil[1]-3.674234614174767*hamil[2])/q_ * normcurlbhat_quad ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.08944271909999158*(8.660254037844388*JfL[5]+5.0*JfL[4])-0.6708203932499369*(1.732050807568877*JfL[3]+JfL[2])+0.5*(1.732050807568877*JfL[1]+JfL[0]);
JfR_quad =  (-0.08944271909999158*(8.660254037844388*JfR[5]-5.0*JfR[4]))+0.6708203932499369*(1.732050807568877*JfR[3]-1.0*JfR[2])-0.5*(1.732050807568877*JfR[1]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((1.0*hamil[1])/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((1.0*hamil[1])/vmap[1])/m_/bmag_quad * (1.224744871391589*hamil[1])/q_ * normcurlbhat_quad ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.5*(1.732050807568877*JfL[1]+JfL[0])-0.1118033988749895*(8.660254037844388*JfL[5]+5.0*JfL[4]);
JfR_quad =  0.1118033988749895*(8.660254037844388*JfR[5]-5.0*JfR[4])-0.5*(1.732050807568877*JfR[1]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;
alpha_quad = ((0.8164965809277262*(3.674234614174767*hamil[2]+1.224744871391589*hamil[1]))/vmap[1])/m_/bmag_quad * B3_quad ;
alpha_quad += ((0.8164965809277262*(3.674234614174767*hamil[2]+1.224744871391589*hamil[1]))/vmap[1])/m_/bmag_quad * (3.674234614174767*hamil[2]+1.224744871391589*hamil[1])/q_ * normcurlbhat_quad ;

alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
JfL_quad =  0.08944271909999158*(8.660254037844388*JfL[5]+5.0*JfL[4])+0.6708203932499369*(1.732050807568877*JfL[3]+JfL[2])+0.5*(1.732050807568877*JfL[1]+JfL[0]);
JfR_quad =  (-0.08944271909999158*(8.660254037844388*JfR[5]-5.0*JfR[4]))-0.6708203932499369*(1.732050807568877*JfR[3]-1.0*JfR[2])-0.5*(1.732050807568877*JfR[1]-1.0*JfR[0]);
Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[0], 0) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[0], 1) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[0], 2) ;


  return cfl*1.5*rdx2; 

} 
