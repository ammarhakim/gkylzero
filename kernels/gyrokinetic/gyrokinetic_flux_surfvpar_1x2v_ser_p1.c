#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_flux_surfvpar_1x2v_ser_p1(
    const struct gkyl_basis *basis, const double *w, const double *dxv, 
    const double *vmap, const double *vmapSq, const double q_, const double m_,
    const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv,
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
  // flux_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.

  double rdx2 = 2.0/dxv[0];
  double rdvpar2 = 2.0/dxv[1];

  double hamil[12] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+vmapSq[0]*m_+1.414213562373095*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

double flux_surf_nodal[4]= {0.0}; 
double bmag_quad = 0.0; 
double B3_quad = 0.0; 
double Jc_quad = 0.0; 
double dualcurlbhat_quad[3] = {0.0}; 


bmag_quad = gkdgv[0].bmag; 
B3_quad = gkdgv[0].B3; 
Jc_quad = dgv[0].Jc; 
dualcurlbhat_quad[0] = gkdgv[0].dualcurlbhat.x[0]; 
dualcurlbhat_quad[1] = gkdgv[0].dualcurlbhat.x[1]; 
dualcurlbhat_quad[2] = gkdgv[0].dualcurlbhat.x[2]; 

flux_surf_nodal[0] = (0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])/m_/bmag_quad * B3_quad ;
flux_surf_nodal[0] += (0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])/m_/bmag_quad * dualcurlbhat_quad[0];
flux_surf_nodal[0] = flux_surf_nodal[0]/2.0 * (0.7905694150420947*JfR[11]-0.7905694150420948*JfR[10]-0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]-0.6123724356957944*JfR[7]+0.6123724356957944*JfR[6]+0.3535533905932737*JfR[5]+0.6123724356957944*JfR[4]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]-0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] + 0.7905694150420947*JfL[11]-0.7905694150420948*JfL[10]-0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]-0.6123724356957944*JfL[7]+0.6123724356957944*JfL[6]+0.3535533905932737*JfL[5]+0.6123724356957944*JfL[4]-0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]-0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0]) + fabs(flux_surf_nodal[0])/2.0 * (0.7905694150420947*JfR[11]-0.7905694150420948*JfR[10]-0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]-0.6123724356957944*JfR[7]+0.6123724356957944*JfR[6]+0.3535533905932737*JfR[5]+0.6123724356957944*JfR[4]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]-0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] - (0.7905694150420947*JfL[11]-0.7905694150420948*JfL[10]-0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]-0.6123724356957944*JfL[7]+0.6123724356957944*JfL[6]+0.3535533905932737*JfL[5]+0.6123724356957944*JfL[4]-0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]-0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])) ;
flux_surf_nodal[1] = (0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])/m_/bmag_quad * B3_quad ;
flux_surf_nodal[1] += (0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])/m_/bmag_quad * dualcurlbhat_quad[0];
flux_surf_nodal[1] = flux_surf_nodal[1]/2.0 * ((-0.7905694150420947*JfR[11])+0.7905694150420948*JfR[10]-0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]+0.6123724356957944*JfR[7]-0.6123724356957944*JfR[6]-0.3535533905932737*JfR[5]+0.6123724356957944*JfR[4]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]-0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] + (-0.7905694150420947*JfL[11])+0.7905694150420948*JfL[10]-0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]+0.6123724356957944*JfL[7]-0.6123724356957944*JfL[6]-0.3535533905932737*JfL[5]+0.6123724356957944*JfL[4]+0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]-0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0]) + fabs(flux_surf_nodal[1])/2.0 * ((-0.7905694150420947*JfR[11])+0.7905694150420948*JfR[10]-0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]+0.6123724356957944*JfR[7]-0.6123724356957944*JfR[6]-0.3535533905932737*JfR[5]+0.6123724356957944*JfR[4]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]-0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] - ((-0.7905694150420947*JfL[11])+0.7905694150420948*JfL[10]-0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]+0.6123724356957944*JfL[7]-0.6123724356957944*JfL[6]-0.3535533905932737*JfL[5]+0.6123724356957944*JfL[4]+0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]-0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])) ;

bmag_quad = gkdgv[1].bmag; 
B3_quad = gkdgv[1].B3; 
Jc_quad = dgv[1].Jc; 
dualcurlbhat_quad[0] = gkdgv[1].dualcurlbhat.x[0]; 
dualcurlbhat_quad[1] = gkdgv[1].dualcurlbhat.x[1]; 
dualcurlbhat_quad[2] = gkdgv[1].dualcurlbhat.x[2]; 

flux_surf_nodal[2] = (0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])/m_/bmag_quad * B3_quad ;
flux_surf_nodal[2] += (0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])/m_/bmag_quad * dualcurlbhat_quad[0];
flux_surf_nodal[2] = flux_surf_nodal[2]/2.0 * ((-0.7905694150420947*JfR[11])-0.7905694150420948*JfR[10]+0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]+0.6123724356957944*JfR[7]+0.6123724356957944*JfR[6]-0.3535533905932737*JfR[5]-0.6123724356957944*JfR[4]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]+0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] + (-0.7905694150420947*JfL[11])-0.7905694150420948*JfL[10]+0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]+0.6123724356957944*JfL[7]+0.6123724356957944*JfL[6]-0.3535533905932737*JfL[5]-0.6123724356957944*JfL[4]-0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]+0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0]) + fabs(flux_surf_nodal[2])/2.0 * ((-0.7905694150420947*JfR[11])-0.7905694150420948*JfR[10]+0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]+0.6123724356957944*JfR[7]+0.6123724356957944*JfR[6]-0.3535533905932737*JfR[5]-0.6123724356957944*JfR[4]-0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]+0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] - ((-0.7905694150420947*JfL[11])-0.7905694150420948*JfL[10]+0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]+0.6123724356957944*JfL[7]+0.6123724356957944*JfL[6]-0.3535533905932737*JfL[5]-0.6123724356957944*JfL[4]-0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]+0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])) ;
flux_surf_nodal[3] = (0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])/m_/bmag_quad * B3_quad ;
flux_surf_nodal[3] += (0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])/m_/bmag_quad * dualcurlbhat_quad[0];
flux_surf_nodal[3] = flux_surf_nodal[3]/2.0 * (0.7905694150420947*JfR[11]+0.7905694150420948*JfR[10]+0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]-0.6123724356957944*JfR[7]-0.6123724356957944*JfR[6]+0.3535533905932737*JfR[5]-0.6123724356957944*JfR[4]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]+0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] + 0.7905694150420947*JfL[11]+0.7905694150420948*JfL[10]+0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]-0.6123724356957944*JfL[7]-0.6123724356957944*JfL[6]+0.3535533905932737*JfL[5]-0.6123724356957944*JfL[4]+0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]+0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0]) + fabs(flux_surf_nodal[3])/2.0 * (0.7905694150420947*JfR[11]+0.7905694150420948*JfR[10]+0.7905694150420948*JfR[9]+0.7905694150420947*JfR[8]-0.6123724356957944*JfR[7]-0.6123724356957944*JfR[6]+0.3535533905932737*JfR[5]-0.6123724356957944*JfR[4]+0.3535533905932737*JfR[3]-0.6123724356957944*JfR[2]+0.3535533905932737*JfR[1]+0.3535533905932737*JfR[0] - (0.7905694150420947*JfL[11]+0.7905694150420948*JfL[10]+0.7905694150420948*JfL[9]+0.7905694150420947*JfL[8]-0.6123724356957944*JfL[7]-0.6123724356957944*JfL[6]+0.3535533905932737*JfL[5]-0.6123724356957944*JfL[4]+0.3535533905932737*JfL[3]-0.6123724356957944*JfL[2]+0.3535533905932737*JfL[1]+0.3535533905932737*JfL[0])) ;
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
basis->quad_nodal_to_modal(flux_surf_nodal, &flux_surf[12], 12) ;

} 
