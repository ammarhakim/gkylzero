#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_1x2v_ser_p1(
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

  double hamil[12] = {0.}; 
  hamil[0] = (2.449489742783179*phi[1]+1.414213562373095*phi[0])*q_+0.7071067811865475*vmapSq[0]*m_+(1.732050807568877*bmag[1]+bmag[0])*vmap[2]; 
  hamil[1] = 0.7071067811865475*vmapSq[1]*m_; 
  hamil[2] = (1.732050807568877*bmag[1]+bmag[0])*vmap[3]; 
  hamil[4] = 0.7071067811865475*vmapSq[2]*m_; 

  double flux_surf_nodal[6]= {0.0}; 
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
  double mvpar_quad[3] = {0.0}; 
  mvpar_quad[0] = (0.8164965809277262*(0.8660254037844386*hamil[1]-2.598076211353316*hamil[4]))/vmap[1]; 
  mvpar_quad[1] = (0.7071067811865476*hamil[1])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277262*(2.598076211353316*hamil[4]+0.8660254037844386*hamil[1]))/vmap[1]; 
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


  alpha_quad = mvpar_quad[0]*B3_quad/(m_*bmag_quad)  ;
  alpha_quad += mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) ;

  alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-0.04472135954999579*(12.24744871391589*JfL[11]+7.071067811865476*JfL[10]))+0.04472135954999579*(12.24744871391589*JfL[9]+7.071067811865476*JfL[8])+0.3354101966249685*(2.449489742783179*JfL[7]+1.414213562373095*JfL[6])-0.25*(2.449489742783179*JfL[5]+1.414213562373095*JfL[3])-0.3354101966249685*(2.449489742783179*JfL[4]+1.414213562373095*JfL[2])+0.25*(2.449489742783179*JfL[1]+1.414213562373095*JfL[0]);
  JfR_quad =  0.04472135954999579*(12.24744871391589*JfR[11]-7.071067811865476*JfR[10])-0.04472135954999579*(12.24744871391589*JfR[9]-7.071067811865476*JfR[8])-0.3354101966249685*(2.449489742783179*JfR[7]-1.414213562373095*JfR[6])+0.25*(2.449489742783179*JfR[5]-1.414213562373095*JfR[3])+0.3354101966249685*(2.449489742783179*JfR[4]-1.414213562373095*JfR[2])-0.25*(2.449489742783179*JfR[1]-1.414213562373095*JfR[0]);
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[0] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = mvpar_quad[1]*B3_quad/(m_*bmag_quad)  ;
  alpha_quad += mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) ;

  alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  0.05590169943749476*(12.24744871391589*JfL[11]+7.071067811865476*JfL[10])-0.05590169943749474*(12.24744871391589*JfL[9]+7.071067811865476*JfL[8])-0.25*(2.449489742783179*JfL[5]+1.414213562373095*JfL[3])+0.25*(2.449489742783179*JfL[1]+1.414213562373095*JfL[0]);
  JfR_quad =  (-0.05590169943749476*(12.24744871391589*JfR[11]-7.071067811865476*JfR[10]))+0.05590169943749474*(12.24744871391589*JfR[9]-7.071067811865476*JfR[8])+0.25*(2.449489742783179*JfR[5]-1.414213562373095*JfR[3])-0.25*(2.449489742783179*JfR[1]-1.414213562373095*JfR[0]);
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[1] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = mvpar_quad[2]*B3_quad/(m_*bmag_quad)  ;
  alpha_quad += mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) ;

  alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-0.04472135954999579*(12.24744871391589*JfL[11]+7.071067811865476*JfL[10]))+0.04472135954999579*(12.24744871391589*JfL[9]+7.071067811865476*JfL[8])-0.3354101966249685*(2.449489742783179*JfL[7]+1.414213562373095*JfL[6])-0.25*(2.449489742783179*JfL[5]+1.414213562373095*JfL[3])+0.3354101966249685*(2.449489742783179*JfL[4]+1.414213562373095*JfL[2])+0.25*(2.449489742783179*JfL[1]+1.414213562373095*JfL[0]);
  JfR_quad =  0.04472135954999579*(12.24744871391589*JfR[11]-7.071067811865476*JfR[10])-0.04472135954999579*(12.24744871391589*JfR[9]-7.071067811865476*JfR[8])+0.3354101966249685*(2.449489742783179*JfR[7]-1.414213562373095*JfR[6])+0.25*(2.449489742783179*JfR[5]-1.414213562373095*JfR[3])-0.3354101966249685*(2.449489742783179*JfR[4]-1.414213562373095*JfR[2])-0.25*(2.449489742783179*JfR[1]-1.414213562373095*JfR[0]);
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[2] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = mvpar_quad[0]*B3_quad/(m_*bmag_quad)  ;
  alpha_quad += mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) ;

  alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  0.04472135954999579*(12.24744871391589*JfL[11]+7.071067811865476*JfL[10])+0.04472135954999579*(12.24744871391589*JfL[9]+7.071067811865476*JfL[8])-0.3354101966249685*(2.449489742783179*JfL[7]+1.414213562373095*JfL[6])+0.25*(2.449489742783179*JfL[5]+1.414213562373095*JfL[3])-0.3354101966249685*(2.449489742783179*JfL[4]+1.414213562373095*JfL[2])+0.25*(2.449489742783179*JfL[1]+1.414213562373095*JfL[0]);
  JfR_quad =  (-0.04472135954999579*(12.24744871391589*JfR[11]-7.071067811865476*JfR[10]))-0.04472135954999579*(12.24744871391589*JfR[9]-7.071067811865476*JfR[8])+0.3354101966249685*(2.449489742783179*JfR[7]-1.414213562373095*JfR[6])-0.25*(2.449489742783179*JfR[5]-1.414213562373095*JfR[3])+0.3354101966249685*(2.449489742783179*JfR[4]-1.414213562373095*JfR[2])-0.25*(2.449489742783179*JfR[1]-1.414213562373095*JfR[0]);
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[3] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = mvpar_quad[1]*B3_quad/(m_*bmag_quad)  ;
  alpha_quad += mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) ;

  alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  (-0.05590169943749476*(12.24744871391589*JfL[11]+7.071067811865476*JfL[10]))-0.05590169943749474*(12.24744871391589*JfL[9]+7.071067811865476*JfL[8])+0.25*(2.449489742783179*JfL[5]+1.414213562373095*JfL[3])+0.25*(2.449489742783179*JfL[1]+1.414213562373095*JfL[0]);
  JfR_quad =  0.05590169943749476*(12.24744871391589*JfR[11]-7.071067811865476*JfR[10])+0.05590169943749474*(12.24744871391589*JfR[9]-7.071067811865476*JfR[8])-0.25*(2.449489742783179*JfR[5]-1.414213562373095*JfR[3])-0.25*(2.449489742783179*JfR[1]-1.414213562373095*JfR[0]);
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[4] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  alpha_quad = mvpar_quad[2]*B3_quad/(m_*bmag_quad)  ;
  alpha_quad += mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) ;

  alpha_quad = alpha_quad*area_elem_quad/Jc_quad ;
  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  JfL_quad =  0.04472135954999579*(12.24744871391589*JfL[11]+7.071067811865476*JfL[10])+0.04472135954999579*(12.24744871391589*JfL[9]+7.071067811865476*JfL[8])+0.3354101966249685*(2.449489742783179*JfL[7]+1.414213562373095*JfL[6])+0.25*(2.449489742783179*JfL[5]+1.414213562373095*JfL[3])+0.3354101966249685*(2.449489742783179*JfL[4]+1.414213562373095*JfL[2])+0.25*(2.449489742783179*JfL[1]+1.414213562373095*JfL[0]);
  JfR_quad =  (-0.04472135954999579*(12.24744871391589*JfR[11]-7.071067811865476*JfR[10]))-0.04472135954999579*(12.24744871391589*JfR[9]-7.071067811865476*JfR[8])-0.3354101966249685*(2.449489742783179*JfR[7]-1.414213562373095*JfR[6])-0.25*(2.449489742783179*JfR[5]-1.414213562373095*JfR[3])-0.3354101966249685*(2.449489742783179*JfR[4]-1.414213562373095*JfR[2])-0.25*(2.449489742783179*JfR[1]-1.414213562373095*JfR[0]);
  Jfavg_quad = (JfL_quad + JfR_quad)/2.0 ;
  Jfjump_quad = (JfR_quad - JfL_quad)/2.0 ;
  flux_surf_nodal[5] = alpha_quad*Jfavg_quad - fabs(alpha_quad)*Jfjump_quad ;

  double *fmodal = &flux_surf[0]; 
  fmodal[0] = 0.2777777777777778*flux_surf_nodal[5]+0.4444444444444444*flux_surf_nodal[4]+0.2777777777777778*flux_surf_nodal[3]+0.2777777777777778*flux_surf_nodal[2]+0.4444444444444444*flux_surf_nodal[1]+0.2777777777777778*flux_surf_nodal[0]; 
  fmodal[1] = 0.372677996249965*flux_surf_nodal[5]-0.372677996249965*flux_surf_nodal[3]+0.372677996249965*flux_surf_nodal[2]-0.372677996249965*flux_surf_nodal[0]; 
  fmodal[2] = 0.2777777777777778*flux_surf_nodal[5]+0.4444444444444444*flux_surf_nodal[4]+0.2777777777777778*flux_surf_nodal[3]-0.2777777777777778*flux_surf_nodal[2]-0.4444444444444444*flux_surf_nodal[1]-0.2777777777777778*flux_surf_nodal[0]; 
  fmodal[3] = 0.372677996249965*flux_surf_nodal[5]-0.372677996249965*flux_surf_nodal[3]-0.372677996249965*flux_surf_nodal[2]+0.372677996249965*flux_surf_nodal[0]; 
  fmodal[4] = 0.2484519974999766*flux_surf_nodal[5]-0.4969039949999533*flux_surf_nodal[4]+0.2484519974999766*flux_surf_nodal[3]+0.2484519974999766*flux_surf_nodal[2]-0.4969039949999533*flux_surf_nodal[1]+0.2484519974999766*flux_surf_nodal[0]; 
  fmodal[5] = 0.2484519974999767*flux_surf_nodal[5]-0.4969039949999534*flux_surf_nodal[4]+0.2484519974999767*flux_surf_nodal[3]-0.2484519974999767*flux_surf_nodal[2]+0.4969039949999534*flux_surf_nodal[1]-0.2484519974999767*flux_surf_nodal[0]; 


  return cfl*1.5*rdx2; 

} 
