#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_flux_surfvpar_1x1v_ser_p1(const struct gkyl_basis *basis, const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_,
    const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, const double *bmag, 
    const double *phi, const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf) 
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
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double hamil[6] = {0.}; 
  hamil[0] = 0.7071067811865475*vmapSq[0]*mass+1.414213562373095*phi[0]*charge; 
  hamil[1] = 1.414213562373095*phi[1]*charge; 
  hamil[2] = 0.7071067811865475*vmapSq[1]*mass; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*mass; 

double flux_surf_nodal[2]= {0.0}; 
double bmag_quad = 0.0; 
double B3_quad = 0.0; 
double Jc_quad = 0.0; 
double dualcurlbhat_quad[3] = {0.0}; 


bmag_quad = gkdgv[0].bmag; 
B3_quad = gkdgv[0].B3; 
Jc_quad = dgv[0].area_elem; 
dualcurlbhat_quad[0] = gkdgv[0].dualcurlbhat.x[0]; 
dualcurlbhat_quad[1] = gkdgv[0].dualcurlbhat.x[1]; 
dualcurlbhat_quad[2] = gkdgv[0].dualcurlbhat.x[2]; 

flux_surf_nodal[0] = (0.8660254037844386*hamil[1])/m_/bmag_quad * B3_quad ;
flux_surf_nodal[0] += (0.8660254037844386*hamil[1])/m_/bmag_quad * dualcurlbhat_quad[k];
flux_surf_nodal[0] = flux_surf_nodal[0]/2.0 * (1.677050983124842*JfR[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfR[5]*(0.1924500897298753-0.5773502691896258*vpar^2)-0.8660254037844386*JfR[3]*vpar+0.8660254037844386*JfR[2]*vpar-0.5*JfR[1]+0.5*JfR[0] + 1.677050983124842*JfL[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfL[5]*(0.1924500897298753-0.5773502691896258*vpar^2)-0.8660254037844386*JfL[3]*vpar+0.8660254037844386*JfL[2]*vpar-0.5*JfL[1]+0.5*JfL[0]) + fabs(flux_surf_nodal[0])/2.0 * (1.677050983124842*JfR[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfR[5]*(0.1924500897298753-0.5773502691896258*vpar^2)-0.8660254037844386*JfR[3]*vpar+0.8660254037844386*JfR[2]*vpar-0.5*JfR[1]+0.5*JfR[0] - (1.677050983124842*JfL[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfL[5]*(0.1924500897298753-0.5773502691896258*vpar^2)-0.8660254037844386*JfL[3]*vpar+0.8660254037844386*JfL[2]*vpar-0.5*JfL[1]+0.5*JfL[0])) ;

bmag_quad = gkdgv[1].bmag; 
B3_quad = gkdgv[1].B3; 
Jc_quad = dgv[1].area_elem; 
dualcurlbhat_quad[0] = gkdgv[1].dualcurlbhat.x[0]; 
dualcurlbhat_quad[1] = gkdgv[1].dualcurlbhat.x[1]; 
dualcurlbhat_quad[2] = gkdgv[1].dualcurlbhat.x[2]; 

flux_surf_nodal[1] = (0.8660254037844386*hamil[1])/m_/bmag_quad * B3_quad ;
flux_surf_nodal[1] += (0.8660254037844386*hamil[1])/m_/bmag_quad * dualcurlbhat_quad[k];
flux_surf_nodal[1] = flux_surf_nodal[1]/2.0 * (1.677050983124842*JfR[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfR[5]*(0.5773502691896258*vpar^2-0.1924500897298753)+0.8660254037844386*JfR[3]*vpar+0.8660254037844386*JfR[2]*vpar+0.5*JfR[1]+0.5*JfR[0] + 1.677050983124842*JfL[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfL[5]*(0.5773502691896258*vpar^2-0.1924500897298753)+0.8660254037844386*JfL[3]*vpar+0.8660254037844386*JfL[2]*vpar+0.5*JfL[1]+0.5*JfL[0]) + fabs(flux_surf_nodal[1])/2.0 * (1.677050983124842*JfR[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfR[5]*(0.5773502691896258*vpar^2-0.1924500897298753)+0.8660254037844386*JfR[3]*vpar+0.8660254037844386*JfR[2]*vpar+0.5*JfR[1]+0.5*JfR[0] - (1.677050983124842*JfL[4]*(vpar^2-0.3333333333333333)+2.904737509655563*JfL[5]*(0.5773502691896258*vpar^2-0.1924500897298753)+0.8660254037844386*JfL[3]*vpar+0.8660254037844386*JfL[2]*vpar+0.5*JfL[1]+0.5*JfL[0])) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[6], 1) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[6], 2) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[6], 3) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[6], 4) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[6], 5) ;
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[6], 6) ;

} 
