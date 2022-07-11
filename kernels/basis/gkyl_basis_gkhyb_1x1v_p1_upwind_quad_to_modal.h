GKYL_CU_DH static inline void 
gkhyb_1x1v_p1_xdir_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.392837100659193*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.392837100659193*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 
  fUpwind[2] = 0.3513641844631533*fUpwindQuad[2]-0.7027283689263066*fUpwindQuad[1]+0.3513641844631533*fUpwindQuad[0]; 

} 

GKYL_CU_DH static inline void 
gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.5*fUpwindQuad[1]+0.5*fUpwindQuad[0]; 
  fUpwind[1] = 0.5*fUpwindQuad[1]-0.5*fUpwindQuad[0]; 
  fUpwind[2] = 0.8660254037844386*fUpwindQuad[1]*mu+0.8660254037844386*fUpwindQuad[0]*mu; 
  fUpwind[3] = 0.8660254037844386*fUpwindQuad[1]*mu-0.8660254037844386*fUpwindQuad[0]*mu; 

} 

