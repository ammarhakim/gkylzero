#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_flux_edge_surfx_2x2v_ser_p1(
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
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  double hamil[12] = {0.}; 
  hamil[0] = (2.449489742783178*phi[1]+1.414213562373095*phi[0])*q_+vmapSq[0]*m_+(1.732050807568877*bmag[1]+bmag[0])*vmap[2]; 
  hamil[1] = (2.449489742783178*phi[3]+1.414213562373095*phi[2])*q_+vmap[2]*(1.732050807568877*bmag[3]+bmag[2]); 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = (1.732050807568877*bmag[1]+bmag[0])*vmap[3]; 
  hamil[5] = (1.732050807568877*bmag[3]+bmag[2])*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

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
  mvpar_quad[0] = (0.8164965809277261*(0.6123724356957944*hamil[2]-1.837117307087383*hamil[8]))/vmap[1]; 
  mvpar_quad[1] = (0.5*hamil[2])/vmap[1]; 
  mvpar_quad[2] = (0.8164965809277261*(1.837117307087383*hamil[8]+0.6123724356957944*hamil[2]))/vmap[1]; 
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


  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.3872983346207417*JfL[23]+0.223606797749979*JfL[22]-0.3872983346207416*(JfL[21]+JfL[20])-0.2236067977499789*JfL[19]-0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]-0.5809475019311124*JfL[15]-0.3354101966249685*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.3354101966249685*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]+0.3354101966249685*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]-0.3354101966249685*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[0] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.3872983346207417*JfR[23])+0.223606797749979*JfR[22]+0.3872983346207416*(JfR[21]+JfR[20])-0.2236067977499789*JfR[19]-0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]+0.5809475019311124*JfR[15]-0.3354101966249685*JfR[14]-0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]+0.3354101966249685*JfR[10]+0.25*JfR[9]+0.4330127018922193*JfR[8]+0.3354101966249685*JfR[7]+0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]-0.25*JfR[4]-0.3354101966249685*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[0] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]+0.4841229182759271*(JfL[21]+JfL[20])+0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]-0.4330127018922193*(JfL[8]+JfL[5])-0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[1] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.4841229182759271*JfR[23]-0.2795084971874737*JfR[22]-0.4841229182759271*(JfR[21]+JfR[20])+0.2795084971874738*(JfR[19]+JfR[18])+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]-0.4330127018922193*JfR[12]+0.25*JfR[9]+0.4330127018922193*(JfR[8]+JfR[5])-0.25*(JfR[4]+JfR[2])-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[1] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.3872983346207417*JfL[23]+0.223606797749979*JfL[22]-0.3872983346207416*(JfL[21]+JfL[20])-0.2236067977499789*JfL[19]-0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]+0.5809475019311124*JfL[15]+0.3354101966249685*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.3354101966249685*JfL[10]+0.25*JfL[9]-0.4330127018922193*JfL[8]-0.3354101966249685*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]-0.25*JfL[4]+0.3354101966249685*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[2] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.3872983346207417*JfR[23])+0.223606797749979*JfR[22]+0.3872983346207416*(JfR[21]+JfR[20])-0.2236067977499789*JfR[19]-0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]-0.5809475019311124*JfR[15]+0.3354101966249685*JfR[14]+0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]-0.3354101966249685*JfR[10]+0.25*JfR[9]+0.4330127018922193*JfR[8]-0.3354101966249685*JfR[7]-0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]-0.25*JfR[4]+0.3354101966249685*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[2] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.3872983346207417*JfL[23])-0.223606797749979*JfL[22]+0.3872983346207416*JfL[21]-0.3872983346207416*JfL[20]+0.2236067977499789*JfL[19]-0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]+0.5809475019311124*JfL[15]+0.3354101966249685*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.3354101966249685*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]+0.3354101966249685*JfL[7]-0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.3354101966249685*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[3] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.3872983346207417*JfR[23]-0.223606797749979*JfR[22]-0.3872983346207416*JfR[21]+0.3872983346207416*JfR[20]+0.2236067977499789*JfR[19]-0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]-0.5809475019311124*JfR[15]+0.3354101966249685*JfR[14]+0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]-0.3354101966249685*JfR[10]-0.25*JfR[9]-0.4330127018922193*JfR[8]+0.3354101966249685*JfR[7]+0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]+0.25*JfR[4]-0.3354101966249685*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[3] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]-0.4841229182759271*JfL[21]+0.4841229182759271*JfL[20]-0.2795084971874738*JfL[19]+0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.4330127018922193*JfL[5]+0.25*JfL[4]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[4] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.4841229182759271*JfR[23])+0.2795084971874737*JfR[22]+0.4841229182759271*JfR[21]-0.4841229182759271*JfR[20]-0.2795084971874738*JfR[19]+0.2795084971874738*JfR[18]+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]+0.4330127018922193*JfR[12]-0.25*JfR[9]-0.4330127018922193*JfR[8]+0.4330127018922193*JfR[5]+0.25*JfR[4]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[4] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.3872983346207417*JfL[23])-0.223606797749979*JfL[22]+0.3872983346207416*JfL[21]-0.3872983346207416*JfL[20]+0.2236067977499789*JfL[19]-0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]-0.5809475019311124*JfL[15]-0.3354101966249685*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.3354101966249685*JfL[10]-0.25*JfL[9]+0.4330127018922193*JfL[8]-0.3354101966249685*JfL[7]+0.5809475019311124*JfL[6]-0.4330127018922193*JfL[5]+0.25*JfL[4]+0.3354101966249685*JfL[3]-0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[5] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.3872983346207417*JfR[23]-0.223606797749979*JfR[22]-0.3872983346207416*JfR[21]+0.3872983346207416*JfR[20]+0.2236067977499789*JfR[19]-0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]+0.5809475019311124*JfR[15]-0.3354101966249685*JfR[14]-0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]+0.3354101966249685*JfR[10]-0.25*JfR[9]-0.4330127018922193*JfR[8]-0.3354101966249685*JfR[7]-0.5809475019311124*JfR[6]+0.4330127018922193*JfR[5]+0.25*JfR[4]+0.3354101966249685*JfR[3]-0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[5] = alpha_quad*Jf_quad;
  }


  bmag_quad = gkdgs[1].bmag; 
  Jc_quad = gkdgs[1].Jc; 
  B3_quad = gkdgs[1].B3; 
  normcurlbhat_quad = gkdgs[1].normcurlbhat; 
  bhat_quad[0] = gkdgs[1].bhat.x[0];  
  bhat_quad[1] = gkdgs[1].bhat.x[1];  
  bhat_quad[2] = gkdgs[1].bhat.x[2];  
  area_elem_quad = dgs[1].area_elem; 


  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.3872983346207417*JfL[23])-0.223606797749979*JfL[22]-0.3872983346207416*JfL[21]+0.3872983346207416*JfL[20]-0.2236067977499789*JfL[19]+0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]+0.5809475019311124*JfL[15]+0.3354101966249685*JfL[14]+0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]+0.3354101966249685*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]-0.3354101966249685*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]-0.3354101966249685*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[6] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.3872983346207417*JfR[23]-0.223606797749979*JfR[22]+0.3872983346207416*JfR[21]-0.3872983346207416*JfR[20]-0.2236067977499789*JfR[19]+0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]-0.5809475019311124*JfR[15]+0.3354101966249685*JfR[14]-0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]+0.3354101966249685*JfR[10]-0.25*JfR[9]+0.4330127018922193*JfR[8]-0.3354101966249685*JfR[7]+0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]-0.25*JfR[4]-0.3354101966249685*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[6] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.4841229182759271*JfL[23]+0.2795084971874737*JfL[22]+0.4841229182759271*JfL[21]-0.4841229182759271*JfL[20]+0.2795084971874738*JfL[19]-0.2795084971874738*JfL[18]-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]-0.4330127018922193*JfL[12]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[7] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.4841229182759271*JfR[23])+0.2795084971874737*JfR[22]-0.4841229182759271*JfR[21]+0.4841229182759271*JfR[20]+0.2795084971874738*JfR[19]-0.2795084971874738*JfR[18]+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]+0.4330127018922193*JfR[12]-0.25*JfR[9]+0.4330127018922193*JfR[8]-0.4330127018922193*JfR[5]-0.25*JfR[4]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[7] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[1]-0.6123724356957944*hamil[5])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.3872983346207417*JfL[23])-0.223606797749979*JfL[22]-0.3872983346207416*JfL[21]+0.3872983346207416*JfL[20]-0.2236067977499789*JfL[19]+0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]-0.5809475019311124*JfL[15]-0.3354101966249685*JfL[14]-0.5809475019311124*JfL[13]-0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]-0.3354101966249685*JfL[10]-0.25*JfL[9]-0.4330127018922193*JfL[8]+0.3354101966249685*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]-0.25*JfL[4]+0.3354101966249685*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[8] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.3872983346207417*JfR[23]-0.223606797749979*JfR[22]+0.3872983346207416*JfR[21]-0.3872983346207416*JfR[20]-0.2236067977499789*JfR[19]+0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]+0.5809475019311124*JfR[15]-0.3354101966249685*JfR[14]+0.5809475019311124*JfR[13]+0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]-0.3354101966249685*JfR[10]-0.25*JfR[9]+0.4330127018922193*JfR[8]+0.3354101966249685*JfR[7]-0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]-0.25*JfR[4]+0.3354101966249685*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[8] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[0]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.3872983346207417*JfL[23]+0.223606797749979*JfL[22]+0.3872983346207416*(JfL[21]+JfL[20])+0.2236067977499789*JfL[19]+0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]-0.5809475019311124*JfL[15]-0.3354101966249685*JfL[14]-0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]-0.5809475019311124*JfL[11]-0.3354101966249685*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]-0.3354101966249685*JfL[7]-0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]-0.3354101966249685*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[9] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.3872983346207417*JfR[23])+0.223606797749979*JfR[22]-0.3872983346207416*(JfR[21]+JfR[20])+0.2236067977499789*JfR[19]+0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]+0.5809475019311124*JfR[15]-0.3354101966249685*JfR[14]+0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]+0.5809475019311124*JfR[11]-0.3354101966249685*JfR[10]+0.25*JfR[9]-0.4330127018922193*JfR[8]-0.3354101966249685*JfR[7]+0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]+0.25*JfR[4]-0.3354101966249685*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[9] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[1]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  (-0.4841229182759271*JfL[23])-0.2795084971874737*JfL[22]-0.4841229182759271*(JfL[21]+JfL[20])-0.2795084971874738*(JfL[19]+JfL[18])-0.4841229182759271*JfL[17]-0.2795084971874737*JfL[16]+0.4330127018922193*JfL[12]+0.25*JfL[9]+0.4330127018922193*(JfL[8]+JfL[5])+0.25*(JfL[4]+JfL[2])+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[10] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  0.4841229182759271*JfR[23]-0.2795084971874737*JfR[22]+0.4841229182759271*(JfR[21]+JfR[20])-0.2795084971874738*(JfR[19]+JfR[18])+0.4841229182759271*JfR[17]-0.2795084971874737*JfR[16]-0.4330127018922193*JfR[12]+0.25*JfR[9]-0.4330127018922193*(JfR[8]+JfR[5])+0.25*(JfR[4]+JfR[2])-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[10] = alpha_quad*Jf_quad;
  }

  alpha_quad = (mvparsq_quad[2]*normcurlbhat_quad/(bmag_quad*q_) + 1/(q_*bmag_quad*area_elem_quad) * bhat_quad[1]*((0.6123724356957944*hamil[5]+0.6123724356957944*hamil[1])*rdz2))*area_elem_quad/Jc_quad  ;

  cfl = fmax(fabs(alpha_quad), fabs(cfl)) ;
  if (alpha_quad > 0) {
    Jf_quad =  0.3872983346207417*JfL[23]+0.223606797749979*JfL[22]+0.3872983346207416*(JfL[21]+JfL[20])+0.2236067977499789*JfL[19]+0.223606797749979*JfL[18]+0.3872983346207417*JfL[17]+0.223606797749979*JfL[16]+0.5809475019311124*JfL[15]+0.3354101966249685*JfL[14]+0.5809475019311124*JfL[13]+0.4330127018922193*JfL[12]+0.5809475019311124*JfL[11]+0.3354101966249685*JfL[10]+0.25*JfL[9]+0.4330127018922193*JfL[8]+0.3354101966249685*JfL[7]+0.5809475019311124*JfL[6]+0.4330127018922193*JfL[5]+0.25*JfL[4]+0.3354101966249685*JfL[3]+0.25*JfL[2]+0.4330127018922193*JfL[1]+0.25*JfL[0];
    flux_surf_nodal[11] = alpha_quad*Jf_quad;
  }
  else {
    Jf_quad =  (-0.3872983346207417*JfR[23])+0.223606797749979*JfR[22]-0.3872983346207416*(JfR[21]+JfR[20])+0.2236067977499789*JfR[19]+0.223606797749979*JfR[18]-0.3872983346207417*JfR[17]+0.223606797749979*JfR[16]-0.5809475019311124*JfR[15]+0.3354101966249685*JfR[14]-0.5809475019311124*JfR[13]-0.4330127018922193*JfR[12]-0.5809475019311124*JfR[11]+0.3354101966249685*JfR[10]+0.25*JfR[9]-0.4330127018922193*JfR[8]+0.3354101966249685*JfR[7]-0.5809475019311124*JfR[6]-0.4330127018922193*JfR[5]+0.25*JfR[4]+0.3354101966249685*JfR[3]+0.25*JfR[2]-0.4330127018922193*JfR[1]+0.25*JfR[0];
    flux_surf_nodal[11] = alpha_quad*Jf_quad;
  }

  return cfl*1.5*rdx2; 

} 
