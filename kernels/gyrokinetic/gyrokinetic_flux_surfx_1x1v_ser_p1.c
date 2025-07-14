#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_flux_surfx_1x1v_ser_p1(
    const double *w, const double *dxv,
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

  double hamil[3] = {0.}; 
  hamil[0] = (phi[0]-1.732050807568877*phi[1])*q_+0.5*vmapSq[0]*m_; 
  hamil[1] = 0.5*vmapSq[1]*m_; 
  hamil[2] = 0.5*vmapSq[2]*m_; 

  double *flux_surf_nodal = &flux_surf[0]; 
  double cfl = 0.0; 
  double bmag_quad = 0.0; 
  double Jc_quad = 0.0; 
  double B3_quad = 0.0; 
  double normcurlbhat_quad = 0.0; 
  double area_elem_quad = 0.0; 
  double bhat_quad[3] = {0.0}; 
  double alpha_quad = 0.0; 
  double Jf_quad = 0.0; 
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277261*(1.224744871391589*hamil[1]-3.674234614174766*hamil[2]))/vmap[1]; 
  mvpar_quad[1] = (1.0*hamil[1])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(3.674234614174766*hamil[2]+1.224744871391589*hamil[1]))/vmap[1]; 
  double mvparsq_quad[3] = {0.0}; 
  mvparsq_quad[0] = mvpar_quad[0]*mvpar_quad[0]/m_; 
  mvparsq_quad[1] = mvpar_quad[1]*mvpar_quad[1]/m_; 
  mvparsq_quad[2] = mvpar_quad[2]*mvpar_quad[2]/m_; 



  bmag_quad = gkdgs[0].bmag; 
  Jc_quad = gkdgs[0].Jc; 
  B3_quad = gkdgs[0].B3; 
  normcurlbhat_quad = gkdgs[0].normcurlbhat; 
  bhat_quad[0] = gkdgs[0].bhat.x[0];  
  bhat_quad[1] = gkdgs[0].bhat.x[1];  
  bhat_quad[2] = gkdgs[0].bhat.x[2];  
  area_elem_quad = dgs[0].area_elem; 


  alpha_quad = (mvpar_quad[0]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]-1.161895003862225*JfL[3]-0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0];
    flux_surf_nodal[0] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.7745966692414834*JfR[5])+0.4472135954999579*JfR[4]+1.161895003862225*JfR[3]-0.6708203932499369*JfR[2]-0.8660254037844386*JfR[1]+0.5*JfR[0];
    flux_surf_nodal[0] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvpar_quad[1]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.9682458365518543*JfL[5])-0.5590169943749475*JfL[4]+0.8660254037844386*JfL[1]+0.5*JfL[0];
    flux_surf_nodal[1] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.9682458365518543*JfR[5]-0.5590169943749475*JfR[4]-0.8660254037844386*JfR[1]+0.5*JfR[0];
    flux_surf_nodal[1] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvpar_quad[2]*B3_quad/(m_*bmag_quad))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.7745966692414834*JfL[5]+0.4472135954999579*JfL[4]+1.161895003862225*JfL[3]+0.6708203932499369*JfL[2]+0.8660254037844386*JfL[1]+0.5*JfL[0];
    flux_surf_nodal[2] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.7745966692414834*JfR[5])+0.4472135954999579*JfR[4]-1.161895003862225*JfR[3]+0.6708203932499369*JfR[2]-0.8660254037844386*JfR[1]+0.5*JfR[0];
    flux_surf_nodal[2] = alpha_quad*Jf_quad;
  }

  return cfl*1.5*rdx2; 

} 
