#include <gkyl_dg_vlasov_gen_geo_alpha_kernels.h> 
GKYL_CU_DH int vlasov_gen_geo_alpha_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing.
  // tvComp[72]: Components for tangent basis vectors.
  // gij[48]:    Contravariant components of metric tensor.
  // alpha_surf: output surface phase space flux, v^i, in each direction (cdim components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  const double *gix = &gij[0]; 
  const double *giy = &gij[8]; 
  const double *giz = &gij[16]; 

  const double *e_1x = &tvComp[0]; 
  const double *e_1y = &tvComp[8]; 
  const double *e_1z = &tvComp[16]; 
  const double *e_2x = &tvComp[24]; 
  const double *e_2y = &tvComp[32]; 
  const double *e_2z = &tvComp[40]; 
  const double *e_3x = &tvComp[48]; 
  const double *e_3y = &tvComp[56]; 
  const double *e_3z = &tvComp[64]; 

  double *alpha_surfL = &alpha_surf[0]; 
  double *sgn_alpha_surfL = &sgn_alpha_surf[0];

  double ei_x_surf[4] = {0.0};
  double ei_y_surf[4] = {0.0};
  double ei_z_surf[4] = {0.0};

  ei_x_surf[0] = 0.75*e_3x[7]*giz[7]-0.4330127018922193*e_3x[6]*giz[7]+0.75*e_2x[7]*giy[7]-0.4330127018922193*e_2x[6]*giy[7]+0.75*e_1x[7]*gix[7]-0.4330127018922193*e_1x[6]*gix[7]-0.4330127018922193*giz[6]*e_3x[7]-0.4330127018922193*giy[6]*e_2x[7]-0.4330127018922193*gix[6]*e_1x[7]+0.25*e_3x[6]*giz[6]+0.25*e_2x[6]*giy[6]+0.25*e_1x[6]*gix[6]+0.75*e_3x[5]*giz[5]-0.4330127018922193*e_3x[3]*giz[5]+0.75*e_2x[5]*giy[5]-0.4330127018922193*e_2x[3]*giy[5]+0.75*e_1x[5]*gix[5]-0.4330127018922193*e_1x[3]*gix[5]-0.4330127018922193*giz[3]*e_3x[5]-0.4330127018922193*giy[3]*e_2x[5]-0.4330127018922193*gix[3]*e_1x[5]+0.75*e_3x[4]*giz[4]-0.4330127018922193*e_3x[2]*giz[4]+0.75*e_2x[4]*giy[4]-0.4330127018922193*e_2x[2]*giy[4]+0.75*e_1x[4]*gix[4]-0.4330127018922193*e_1x[2]*gix[4]-0.4330127018922193*giz[2]*e_3x[4]-0.4330127018922193*giy[2]*e_2x[4]-0.4330127018922193*gix[2]*e_1x[4]+0.25*e_3x[3]*giz[3]+0.25*e_2x[3]*giy[3]+0.25*e_1x[3]*gix[3]+0.25*e_3x[2]*giz[2]+0.25*e_2x[2]*giy[2]+0.25*e_1x[2]*gix[2]+0.75*e_3x[1]*giz[1]-0.4330127018922193*e_3x[0]*giz[1]+0.75*e_2x[1]*giy[1]-0.4330127018922193*e_2x[0]*giy[1]+0.75*e_1x[1]*gix[1]-0.4330127018922193*e_1x[0]*gix[1]-0.4330127018922193*giz[0]*e_3x[1]-0.4330127018922193*giy[0]*e_2x[1]-0.4330127018922193*gix[0]*e_1x[1]+0.25*e_3x[0]*giz[0]+0.25*e_2x[0]*giy[0]+0.25*e_1x[0]*gix[0]; 
  ei_x_surf[1] = 0.75*e_3x[5]*giz[7]-0.4330127018922193*e_3x[3]*giz[7]+0.75*e_2x[5]*giy[7]-0.4330127018922193*e_2x[3]*giy[7]+0.75*e_1x[5]*gix[7]-0.4330127018922193*e_1x[3]*gix[7]+0.75*giz[5]*e_3x[7]-0.4330127018922193*giz[3]*e_3x[7]+0.75*giy[5]*e_2x[7]-0.4330127018922193*giy[3]*e_2x[7]+0.75*gix[5]*e_1x[7]-0.4330127018922193*gix[3]*e_1x[7]-0.4330127018922193*e_3x[5]*giz[6]+0.25*e_3x[3]*giz[6]-0.4330127018922193*e_2x[5]*giy[6]+0.25*e_2x[3]*giy[6]-0.4330127018922193*e_1x[5]*gix[6]+0.25*e_1x[3]*gix[6]-0.4330127018922193*giz[5]*e_3x[6]+0.25*giz[3]*e_3x[6]-0.4330127018922193*giy[5]*e_2x[6]+0.25*giy[3]*e_2x[6]-0.4330127018922193*gix[5]*e_1x[6]+0.25*gix[3]*e_1x[6]+0.75*e_3x[1]*giz[4]-0.4330127018922193*e_3x[0]*giz[4]+0.75*e_2x[1]*giy[4]-0.4330127018922193*e_2x[0]*giy[4]+0.75*e_1x[1]*gix[4]-0.4330127018922193*e_1x[0]*gix[4]+0.75*giz[1]*e_3x[4]-0.4330127018922193*giz[0]*e_3x[4]+0.75*giy[1]*e_2x[4]-0.4330127018922193*giy[0]*e_2x[4]+0.75*gix[1]*e_1x[4]-0.4330127018922193*gix[0]*e_1x[4]-0.4330127018922193*e_3x[1]*giz[2]+0.25*e_3x[0]*giz[2]-0.4330127018922193*e_2x[1]*giy[2]+0.25*e_2x[0]*giy[2]-0.4330127018922193*e_1x[1]*gix[2]+0.25*e_1x[0]*gix[2]-0.4330127018922193*giz[1]*e_3x[2]+0.25*giz[0]*e_3x[2]-0.4330127018922193*giy[1]*e_2x[2]+0.25*giy[0]*e_2x[2]-0.4330127018922193*gix[1]*e_1x[2]+0.25*gix[0]*e_1x[2]; 
  ei_x_surf[2] = 0.75*e_3x[4]*giz[7]-0.4330127018922193*e_3x[2]*giz[7]+0.75*e_2x[4]*giy[7]-0.4330127018922193*e_2x[2]*giy[7]+0.75*e_1x[4]*gix[7]-0.4330127018922193*e_1x[2]*gix[7]+0.75*giz[4]*e_3x[7]-0.4330127018922193*giz[2]*e_3x[7]+0.75*giy[4]*e_2x[7]-0.4330127018922193*giy[2]*e_2x[7]+0.75*gix[4]*e_1x[7]-0.4330127018922193*gix[2]*e_1x[7]-0.4330127018922193*e_3x[4]*giz[6]+0.25*e_3x[2]*giz[6]-0.4330127018922193*e_2x[4]*giy[6]+0.25*e_2x[2]*giy[6]-0.4330127018922193*e_1x[4]*gix[6]+0.25*e_1x[2]*gix[6]-0.4330127018922193*giz[4]*e_3x[6]+0.25*giz[2]*e_3x[6]-0.4330127018922193*giy[4]*e_2x[6]+0.25*giy[2]*e_2x[6]-0.4330127018922193*gix[4]*e_1x[6]+0.25*gix[2]*e_1x[6]+0.75*e_3x[1]*giz[5]-0.4330127018922193*e_3x[0]*giz[5]+0.75*e_2x[1]*giy[5]-0.4330127018922193*e_2x[0]*giy[5]+0.75*e_1x[1]*gix[5]-0.4330127018922193*e_1x[0]*gix[5]+0.75*giz[1]*e_3x[5]-0.4330127018922193*giz[0]*e_3x[5]+0.75*giy[1]*e_2x[5]-0.4330127018922193*giy[0]*e_2x[5]+0.75*gix[1]*e_1x[5]-0.4330127018922193*gix[0]*e_1x[5]-0.4330127018922193*e_3x[1]*giz[3]+0.25*e_3x[0]*giz[3]-0.4330127018922193*e_2x[1]*giy[3]+0.25*e_2x[0]*giy[3]-0.4330127018922193*e_1x[1]*gix[3]+0.25*e_1x[0]*gix[3]-0.4330127018922193*giz[1]*e_3x[3]+0.25*giz[0]*e_3x[3]-0.4330127018922193*giy[1]*e_2x[3]+0.25*giy[0]*e_2x[3]-0.4330127018922193*gix[1]*e_1x[3]+0.25*gix[0]*e_1x[3]; 
  ei_x_surf[3] = 0.75*e_3x[1]*giz[7]-0.4330127018922193*e_3x[0]*giz[7]+0.75*e_2x[1]*giy[7]-0.4330127018922193*e_2x[0]*giy[7]+0.75*e_1x[1]*gix[7]-0.4330127018922193*e_1x[0]*gix[7]+0.75*giz[1]*e_3x[7]-0.4330127018922193*giz[0]*e_3x[7]+0.75*giy[1]*e_2x[7]-0.4330127018922193*giy[0]*e_2x[7]+0.75*gix[1]*e_1x[7]-0.4330127018922193*gix[0]*e_1x[7]-0.4330127018922193*e_3x[1]*giz[6]+0.25*e_3x[0]*giz[6]-0.4330127018922193*e_2x[1]*giy[6]+0.25*e_2x[0]*giy[6]-0.4330127018922193*e_1x[1]*gix[6]+0.25*e_1x[0]*gix[6]-0.4330127018922193*giz[1]*e_3x[6]+0.25*giz[0]*e_3x[6]-0.4330127018922193*giy[1]*e_2x[6]+0.25*giy[0]*e_2x[6]-0.4330127018922193*gix[1]*e_1x[6]+0.25*gix[0]*e_1x[6]+0.75*e_3x[4]*giz[5]-0.4330127018922193*e_3x[2]*giz[5]+0.75*e_2x[4]*giy[5]-0.4330127018922193*e_2x[2]*giy[5]+0.75*e_1x[4]*gix[5]-0.4330127018922193*e_1x[2]*gix[5]+0.75*giz[4]*e_3x[5]-0.4330127018922193*giz[2]*e_3x[5]+0.75*giy[4]*e_2x[5]-0.4330127018922193*giy[2]*e_2x[5]+0.75*gix[4]*e_1x[5]-0.4330127018922193*gix[2]*e_1x[5]-0.4330127018922193*e_3x[3]*giz[4]-0.4330127018922193*e_2x[3]*giy[4]-0.4330127018922193*e_1x[3]*gix[4]-0.4330127018922193*giz[3]*e_3x[4]-0.4330127018922193*giy[3]*e_2x[4]-0.4330127018922193*gix[3]*e_1x[4]+0.25*e_3x[2]*giz[3]+0.25*e_2x[2]*giy[3]+0.25*e_1x[2]*gix[3]+0.25*giz[2]*e_3x[3]+0.25*giy[2]*e_2x[3]+0.25*gix[2]*e_1x[3]; 

  ei_y_surf[0] = 0.75*e_3y[7]*giz[7]-0.4330127018922193*e_3y[6]*giz[7]+0.75*e_2y[7]*giy[7]-0.4330127018922193*e_2y[6]*giy[7]+0.75*e_1y[7]*gix[7]-0.4330127018922193*e_1y[6]*gix[7]-0.4330127018922193*giz[6]*e_3y[7]-0.4330127018922193*giy[6]*e_2y[7]-0.4330127018922193*gix[6]*e_1y[7]+0.25*e_3y[6]*giz[6]+0.25*e_2y[6]*giy[6]+0.25*e_1y[6]*gix[6]+0.75*e_3y[5]*giz[5]-0.4330127018922193*e_3y[3]*giz[5]+0.75*e_2y[5]*giy[5]-0.4330127018922193*e_2y[3]*giy[5]+0.75*e_1y[5]*gix[5]-0.4330127018922193*e_1y[3]*gix[5]-0.4330127018922193*giz[3]*e_3y[5]-0.4330127018922193*giy[3]*e_2y[5]-0.4330127018922193*gix[3]*e_1y[5]+0.75*e_3y[4]*giz[4]-0.4330127018922193*e_3y[2]*giz[4]+0.75*e_2y[4]*giy[4]-0.4330127018922193*e_2y[2]*giy[4]+0.75*e_1y[4]*gix[4]-0.4330127018922193*e_1y[2]*gix[4]-0.4330127018922193*giz[2]*e_3y[4]-0.4330127018922193*giy[2]*e_2y[4]-0.4330127018922193*gix[2]*e_1y[4]+0.25*e_3y[3]*giz[3]+0.25*e_2y[3]*giy[3]+0.25*e_1y[3]*gix[3]+0.25*e_3y[2]*giz[2]+0.25*e_2y[2]*giy[2]+0.25*e_1y[2]*gix[2]+0.75*e_3y[1]*giz[1]-0.4330127018922193*e_3y[0]*giz[1]+0.75*e_2y[1]*giy[1]-0.4330127018922193*e_2y[0]*giy[1]+0.75*e_1y[1]*gix[1]-0.4330127018922193*e_1y[0]*gix[1]-0.4330127018922193*giz[0]*e_3y[1]-0.4330127018922193*giy[0]*e_2y[1]-0.4330127018922193*gix[0]*e_1y[1]+0.25*e_3y[0]*giz[0]+0.25*e_2y[0]*giy[0]+0.25*e_1y[0]*gix[0]; 
  ei_y_surf[1] = 0.75*e_3y[5]*giz[7]-0.4330127018922193*e_3y[3]*giz[7]+0.75*e_2y[5]*giy[7]-0.4330127018922193*e_2y[3]*giy[7]+0.75*e_1y[5]*gix[7]-0.4330127018922193*e_1y[3]*gix[7]+0.75*giz[5]*e_3y[7]-0.4330127018922193*giz[3]*e_3y[7]+0.75*giy[5]*e_2y[7]-0.4330127018922193*giy[3]*e_2y[7]+0.75*gix[5]*e_1y[7]-0.4330127018922193*gix[3]*e_1y[7]-0.4330127018922193*e_3y[5]*giz[6]+0.25*e_3y[3]*giz[6]-0.4330127018922193*e_2y[5]*giy[6]+0.25*e_2y[3]*giy[6]-0.4330127018922193*e_1y[5]*gix[6]+0.25*e_1y[3]*gix[6]-0.4330127018922193*giz[5]*e_3y[6]+0.25*giz[3]*e_3y[6]-0.4330127018922193*giy[5]*e_2y[6]+0.25*giy[3]*e_2y[6]-0.4330127018922193*gix[5]*e_1y[6]+0.25*gix[3]*e_1y[6]+0.75*e_3y[1]*giz[4]-0.4330127018922193*e_3y[0]*giz[4]+0.75*e_2y[1]*giy[4]-0.4330127018922193*e_2y[0]*giy[4]+0.75*e_1y[1]*gix[4]-0.4330127018922193*e_1y[0]*gix[4]+0.75*giz[1]*e_3y[4]-0.4330127018922193*giz[0]*e_3y[4]+0.75*giy[1]*e_2y[4]-0.4330127018922193*giy[0]*e_2y[4]+0.75*gix[1]*e_1y[4]-0.4330127018922193*gix[0]*e_1y[4]-0.4330127018922193*e_3y[1]*giz[2]+0.25*e_3y[0]*giz[2]-0.4330127018922193*e_2y[1]*giy[2]+0.25*e_2y[0]*giy[2]-0.4330127018922193*e_1y[1]*gix[2]+0.25*e_1y[0]*gix[2]-0.4330127018922193*giz[1]*e_3y[2]+0.25*giz[0]*e_3y[2]-0.4330127018922193*giy[1]*e_2y[2]+0.25*giy[0]*e_2y[2]-0.4330127018922193*gix[1]*e_1y[2]+0.25*gix[0]*e_1y[2]; 
  ei_y_surf[2] = 0.75*e_3y[4]*giz[7]-0.4330127018922193*e_3y[2]*giz[7]+0.75*e_2y[4]*giy[7]-0.4330127018922193*e_2y[2]*giy[7]+0.75*e_1y[4]*gix[7]-0.4330127018922193*e_1y[2]*gix[7]+0.75*giz[4]*e_3y[7]-0.4330127018922193*giz[2]*e_3y[7]+0.75*giy[4]*e_2y[7]-0.4330127018922193*giy[2]*e_2y[7]+0.75*gix[4]*e_1y[7]-0.4330127018922193*gix[2]*e_1y[7]-0.4330127018922193*e_3y[4]*giz[6]+0.25*e_3y[2]*giz[6]-0.4330127018922193*e_2y[4]*giy[6]+0.25*e_2y[2]*giy[6]-0.4330127018922193*e_1y[4]*gix[6]+0.25*e_1y[2]*gix[6]-0.4330127018922193*giz[4]*e_3y[6]+0.25*giz[2]*e_3y[6]-0.4330127018922193*giy[4]*e_2y[6]+0.25*giy[2]*e_2y[6]-0.4330127018922193*gix[4]*e_1y[6]+0.25*gix[2]*e_1y[6]+0.75*e_3y[1]*giz[5]-0.4330127018922193*e_3y[0]*giz[5]+0.75*e_2y[1]*giy[5]-0.4330127018922193*e_2y[0]*giy[5]+0.75*e_1y[1]*gix[5]-0.4330127018922193*e_1y[0]*gix[5]+0.75*giz[1]*e_3y[5]-0.4330127018922193*giz[0]*e_3y[5]+0.75*giy[1]*e_2y[5]-0.4330127018922193*giy[0]*e_2y[5]+0.75*gix[1]*e_1y[5]-0.4330127018922193*gix[0]*e_1y[5]-0.4330127018922193*e_3y[1]*giz[3]+0.25*e_3y[0]*giz[3]-0.4330127018922193*e_2y[1]*giy[3]+0.25*e_2y[0]*giy[3]-0.4330127018922193*e_1y[1]*gix[3]+0.25*e_1y[0]*gix[3]-0.4330127018922193*giz[1]*e_3y[3]+0.25*giz[0]*e_3y[3]-0.4330127018922193*giy[1]*e_2y[3]+0.25*giy[0]*e_2y[3]-0.4330127018922193*gix[1]*e_1y[3]+0.25*gix[0]*e_1y[3]; 
  ei_y_surf[3] = 0.75*e_3y[1]*giz[7]-0.4330127018922193*e_3y[0]*giz[7]+0.75*e_2y[1]*giy[7]-0.4330127018922193*e_2y[0]*giy[7]+0.75*e_1y[1]*gix[7]-0.4330127018922193*e_1y[0]*gix[7]+0.75*giz[1]*e_3y[7]-0.4330127018922193*giz[0]*e_3y[7]+0.75*giy[1]*e_2y[7]-0.4330127018922193*giy[0]*e_2y[7]+0.75*gix[1]*e_1y[7]-0.4330127018922193*gix[0]*e_1y[7]-0.4330127018922193*e_3y[1]*giz[6]+0.25*e_3y[0]*giz[6]-0.4330127018922193*e_2y[1]*giy[6]+0.25*e_2y[0]*giy[6]-0.4330127018922193*e_1y[1]*gix[6]+0.25*e_1y[0]*gix[6]-0.4330127018922193*giz[1]*e_3y[6]+0.25*giz[0]*e_3y[6]-0.4330127018922193*giy[1]*e_2y[6]+0.25*giy[0]*e_2y[6]-0.4330127018922193*gix[1]*e_1y[6]+0.25*gix[0]*e_1y[6]+0.75*e_3y[4]*giz[5]-0.4330127018922193*e_3y[2]*giz[5]+0.75*e_2y[4]*giy[5]-0.4330127018922193*e_2y[2]*giy[5]+0.75*e_1y[4]*gix[5]-0.4330127018922193*e_1y[2]*gix[5]+0.75*giz[4]*e_3y[5]-0.4330127018922193*giz[2]*e_3y[5]+0.75*giy[4]*e_2y[5]-0.4330127018922193*giy[2]*e_2y[5]+0.75*gix[4]*e_1y[5]-0.4330127018922193*gix[2]*e_1y[5]-0.4330127018922193*e_3y[3]*giz[4]-0.4330127018922193*e_2y[3]*giy[4]-0.4330127018922193*e_1y[3]*gix[4]-0.4330127018922193*giz[3]*e_3y[4]-0.4330127018922193*giy[3]*e_2y[4]-0.4330127018922193*gix[3]*e_1y[4]+0.25*e_3y[2]*giz[3]+0.25*e_2y[2]*giy[3]+0.25*e_1y[2]*gix[3]+0.25*giz[2]*e_3y[3]+0.25*giy[2]*e_2y[3]+0.25*gix[2]*e_1y[3]; 

  ei_z_surf[0] = 0.75*e_3z[7]*giz[7]-0.4330127018922193*e_3z[6]*giz[7]+0.75*e_2z[7]*giy[7]-0.4330127018922193*e_2z[6]*giy[7]+0.75*e_1z[7]*gix[7]-0.4330127018922193*e_1z[6]*gix[7]-0.4330127018922193*giz[6]*e_3z[7]-0.4330127018922193*giy[6]*e_2z[7]-0.4330127018922193*gix[6]*e_1z[7]+0.25*e_3z[6]*giz[6]+0.25*e_2z[6]*giy[6]+0.25*e_1z[6]*gix[6]+0.75*e_3z[5]*giz[5]-0.4330127018922193*e_3z[3]*giz[5]+0.75*e_2z[5]*giy[5]-0.4330127018922193*e_2z[3]*giy[5]+0.75*e_1z[5]*gix[5]-0.4330127018922193*e_1z[3]*gix[5]-0.4330127018922193*giz[3]*e_3z[5]-0.4330127018922193*giy[3]*e_2z[5]-0.4330127018922193*gix[3]*e_1z[5]+0.75*e_3z[4]*giz[4]-0.4330127018922193*e_3z[2]*giz[4]+0.75*e_2z[4]*giy[4]-0.4330127018922193*e_2z[2]*giy[4]+0.75*e_1z[4]*gix[4]-0.4330127018922193*e_1z[2]*gix[4]-0.4330127018922193*giz[2]*e_3z[4]-0.4330127018922193*giy[2]*e_2z[4]-0.4330127018922193*gix[2]*e_1z[4]+0.25*e_3z[3]*giz[3]+0.25*e_2z[3]*giy[3]+0.25*e_1z[3]*gix[3]+0.25*e_3z[2]*giz[2]+0.25*e_2z[2]*giy[2]+0.25*e_1z[2]*gix[2]+0.75*e_3z[1]*giz[1]-0.4330127018922193*e_3z[0]*giz[1]+0.75*e_2z[1]*giy[1]-0.4330127018922193*e_2z[0]*giy[1]+0.75*e_1z[1]*gix[1]-0.4330127018922193*e_1z[0]*gix[1]-0.4330127018922193*giz[0]*e_3z[1]-0.4330127018922193*giy[0]*e_2z[1]-0.4330127018922193*gix[0]*e_1z[1]+0.25*e_3z[0]*giz[0]+0.25*e_2z[0]*giy[0]+0.25*e_1z[0]*gix[0]; 
  ei_z_surf[1] = 0.75*e_3z[5]*giz[7]-0.4330127018922193*e_3z[3]*giz[7]+0.75*e_2z[5]*giy[7]-0.4330127018922193*e_2z[3]*giy[7]+0.75*e_1z[5]*gix[7]-0.4330127018922193*e_1z[3]*gix[7]+0.75*giz[5]*e_3z[7]-0.4330127018922193*giz[3]*e_3z[7]+0.75*giy[5]*e_2z[7]-0.4330127018922193*giy[3]*e_2z[7]+0.75*gix[5]*e_1z[7]-0.4330127018922193*gix[3]*e_1z[7]-0.4330127018922193*e_3z[5]*giz[6]+0.25*e_3z[3]*giz[6]-0.4330127018922193*e_2z[5]*giy[6]+0.25*e_2z[3]*giy[6]-0.4330127018922193*e_1z[5]*gix[6]+0.25*e_1z[3]*gix[6]-0.4330127018922193*giz[5]*e_3z[6]+0.25*giz[3]*e_3z[6]-0.4330127018922193*giy[5]*e_2z[6]+0.25*giy[3]*e_2z[6]-0.4330127018922193*gix[5]*e_1z[6]+0.25*gix[3]*e_1z[6]+0.75*e_3z[1]*giz[4]-0.4330127018922193*e_3z[0]*giz[4]+0.75*e_2z[1]*giy[4]-0.4330127018922193*e_2z[0]*giy[4]+0.75*e_1z[1]*gix[4]-0.4330127018922193*e_1z[0]*gix[4]+0.75*giz[1]*e_3z[4]-0.4330127018922193*giz[0]*e_3z[4]+0.75*giy[1]*e_2z[4]-0.4330127018922193*giy[0]*e_2z[4]+0.75*gix[1]*e_1z[4]-0.4330127018922193*gix[0]*e_1z[4]-0.4330127018922193*e_3z[1]*giz[2]+0.25*e_3z[0]*giz[2]-0.4330127018922193*e_2z[1]*giy[2]+0.25*e_2z[0]*giy[2]-0.4330127018922193*e_1z[1]*gix[2]+0.25*e_1z[0]*gix[2]-0.4330127018922193*giz[1]*e_3z[2]+0.25*giz[0]*e_3z[2]-0.4330127018922193*giy[1]*e_2z[2]+0.25*giy[0]*e_2z[2]-0.4330127018922193*gix[1]*e_1z[2]+0.25*gix[0]*e_1z[2]; 
  ei_z_surf[2] = 0.75*e_3z[4]*giz[7]-0.4330127018922193*e_3z[2]*giz[7]+0.75*e_2z[4]*giy[7]-0.4330127018922193*e_2z[2]*giy[7]+0.75*e_1z[4]*gix[7]-0.4330127018922193*e_1z[2]*gix[7]+0.75*giz[4]*e_3z[7]-0.4330127018922193*giz[2]*e_3z[7]+0.75*giy[4]*e_2z[7]-0.4330127018922193*giy[2]*e_2z[7]+0.75*gix[4]*e_1z[7]-0.4330127018922193*gix[2]*e_1z[7]-0.4330127018922193*e_3z[4]*giz[6]+0.25*e_3z[2]*giz[6]-0.4330127018922193*e_2z[4]*giy[6]+0.25*e_2z[2]*giy[6]-0.4330127018922193*e_1z[4]*gix[6]+0.25*e_1z[2]*gix[6]-0.4330127018922193*giz[4]*e_3z[6]+0.25*giz[2]*e_3z[6]-0.4330127018922193*giy[4]*e_2z[6]+0.25*giy[2]*e_2z[6]-0.4330127018922193*gix[4]*e_1z[6]+0.25*gix[2]*e_1z[6]+0.75*e_3z[1]*giz[5]-0.4330127018922193*e_3z[0]*giz[5]+0.75*e_2z[1]*giy[5]-0.4330127018922193*e_2z[0]*giy[5]+0.75*e_1z[1]*gix[5]-0.4330127018922193*e_1z[0]*gix[5]+0.75*giz[1]*e_3z[5]-0.4330127018922193*giz[0]*e_3z[5]+0.75*giy[1]*e_2z[5]-0.4330127018922193*giy[0]*e_2z[5]+0.75*gix[1]*e_1z[5]-0.4330127018922193*gix[0]*e_1z[5]-0.4330127018922193*e_3z[1]*giz[3]+0.25*e_3z[0]*giz[3]-0.4330127018922193*e_2z[1]*giy[3]+0.25*e_2z[0]*giy[3]-0.4330127018922193*e_1z[1]*gix[3]+0.25*e_1z[0]*gix[3]-0.4330127018922193*giz[1]*e_3z[3]+0.25*giz[0]*e_3z[3]-0.4330127018922193*giy[1]*e_2z[3]+0.25*giy[0]*e_2z[3]-0.4330127018922193*gix[1]*e_1z[3]+0.25*gix[0]*e_1z[3]; 
  ei_z_surf[3] = 0.75*e_3z[1]*giz[7]-0.4330127018922193*e_3z[0]*giz[7]+0.75*e_2z[1]*giy[7]-0.4330127018922193*e_2z[0]*giy[7]+0.75*e_1z[1]*gix[7]-0.4330127018922193*e_1z[0]*gix[7]+0.75*giz[1]*e_3z[7]-0.4330127018922193*giz[0]*e_3z[7]+0.75*giy[1]*e_2z[7]-0.4330127018922193*giy[0]*e_2z[7]+0.75*gix[1]*e_1z[7]-0.4330127018922193*gix[0]*e_1z[7]-0.4330127018922193*e_3z[1]*giz[6]+0.25*e_3z[0]*giz[6]-0.4330127018922193*e_2z[1]*giy[6]+0.25*e_2z[0]*giy[6]-0.4330127018922193*e_1z[1]*gix[6]+0.25*e_1z[0]*gix[6]-0.4330127018922193*giz[1]*e_3z[6]+0.25*giz[0]*e_3z[6]-0.4330127018922193*giy[1]*e_2z[6]+0.25*giy[0]*e_2z[6]-0.4330127018922193*gix[1]*e_1z[6]+0.25*gix[0]*e_1z[6]+0.75*e_3z[4]*giz[5]-0.4330127018922193*e_3z[2]*giz[5]+0.75*e_2z[4]*giy[5]-0.4330127018922193*e_2z[2]*giy[5]+0.75*e_1z[4]*gix[5]-0.4330127018922193*e_1z[2]*gix[5]+0.75*giz[4]*e_3z[5]-0.4330127018922193*giz[2]*e_3z[5]+0.75*giy[4]*e_2z[5]-0.4330127018922193*giy[2]*e_2z[5]+0.75*gix[4]*e_1z[5]-0.4330127018922193*gix[2]*e_1z[5]-0.4330127018922193*e_3z[3]*giz[4]-0.4330127018922193*e_2z[3]*giy[4]-0.4330127018922193*e_1z[3]*gix[4]-0.4330127018922193*giz[3]*e_3z[4]-0.4330127018922193*giy[3]*e_2z[4]-0.4330127018922193*gix[3]*e_1z[4]+0.25*e_3z[2]*giz[3]+0.25*e_2z[2]*giy[3]+0.25*e_1z[2]*gix[3]+0.25*giz[2]*e_3z[3]+0.25*giy[2]*e_2z[3]+0.25*gix[2]*e_1z[3]; 

  alpha_surfL[0] = 2.828427124746191*ei_z_surf[0]*w[5]+2.828427124746191*ei_y_surf[0]*w[4]+2.828427124746191*ei_x_surf[0]*w[3]; 
  alpha_surfL[1] = 2.828427124746191*ei_z_surf[1]*w[5]+2.828427124746191*ei_y_surf[1]*w[4]+2.828427124746191*ei_x_surf[1]*w[3]; 
  alpha_surfL[2] = 2.828427124746191*ei_z_surf[2]*w[5]+2.828427124746191*ei_y_surf[2]*w[4]+2.828427124746191*ei_x_surf[2]*w[3]; 
  alpha_surfL[3] = 0.8164965809277261*ei_x_surf[0]*dxv[3]; 
  alpha_surfL[4] = 0.8164965809277261*ei_y_surf[0]*dxv[4]; 
  alpha_surfL[5] = 0.8164965809277261*ei_z_surf[0]*dxv[5]; 
  alpha_surfL[6] = 2.828427124746191*ei_z_surf[3]*w[5]+2.828427124746191*ei_y_surf[3]*w[4]+2.828427124746191*ei_x_surf[3]*w[3]; 
  alpha_surfL[7] = 0.8164965809277261*ei_x_surf[1]*dxv[3]; 
  alpha_surfL[8] = 0.8164965809277261*ei_x_surf[2]*dxv[3]; 
  alpha_surfL[9] = 0.8164965809277261*ei_y_surf[1]*dxv[4]; 
  alpha_surfL[10] = 0.8164965809277261*ei_y_surf[2]*dxv[4]; 
  alpha_surfL[12] = 0.8164965809277261*ei_z_surf[1]*dxv[5]; 
  alpha_surfL[13] = 0.8164965809277261*ei_z_surf[2]*dxv[5]; 
  alpha_surfL[16] = 0.8164965809277261*dxv[3]*ei_x_surf[3]; 
  alpha_surfL[17] = 0.8164965809277261*ei_y_surf[3]*dxv[4]; 
  alpha_surfL[20] = 0.8164965809277261*ei_z_surf[3]*dxv[5]; 
  int const_sgn_alpha_surf = 1;  
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16]))+0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6])-0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12])+0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5])-0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*alpha_surfL[17]-0.1767766952966367*alpha_surfL[16]+0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12])-0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9])+0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6])-0.1767766952966367*alpha_surfL[5]+0.1767766952966367*alpha_surfL[4]-0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17])-0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9])+0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4])-0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]))+0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9])-0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*alpha_surfL[6]-0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4])+0.1767766952966367*alpha_surfL[3]-0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*alpha_surfL[17]+0.1767766952966367*alpha_surfL[16]-0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12])+0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9])-0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5])-0.1767766952966367*alpha_surfL[4]+0.1767766952966367*alpha_surfL[3]-0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12])-0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*alpha_surfL[6]-0.1767766952966367*alpha_surfL[5]+0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3])-0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16])-0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3])-0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1])+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16])-0.1767766952966367*alpha_surfL[13]+0.1767766952966367*alpha_surfL[12]-0.1767766952966367*alpha_surfL[10]+0.1767766952966367*alpha_surfL[9]-0.1767766952966367*alpha_surfL[8]+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3])+0.1767766952966367*alpha_surfL[2]-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13])-0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])+0.1767766952966367*alpha_surfL[9]-0.1767766952966367*alpha_surfL[8]+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*alpha_surfL[6]+0.1767766952966367*alpha_surfL[5]-0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3])+0.1767766952966367*alpha_surfL[2]-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*alpha_surfL[17]+0.1767766952966367*alpha_surfL[16]-0.1767766952966367*alpha_surfL[13]+0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])-0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5])+0.1767766952966367*alpha_surfL[4]-0.1767766952966367*alpha_surfL[3]+0.1767766952966367*alpha_surfL[2]-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]))+0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13])-0.1767766952966367*alpha_surfL[12]+0.1767766952966367*alpha_surfL[10]-0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*alpha_surfL[6]+0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4])-0.1767766952966367*alpha_surfL[3]+0.1767766952966367*alpha_surfL[2]-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17])-0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13])+0.1767766952966367*alpha_surfL[12]-0.1767766952966367*alpha_surfL[10]+0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4])+0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2])-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*alpha_surfL[17]-0.1767766952966367*alpha_surfL[16]+0.1767766952966367*alpha_surfL[13]-0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])+0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6])+0.1767766952966367*alpha_surfL[5]-0.1767766952966367*alpha_surfL[4]+0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2])-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13])+0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])-0.1767766952966367*alpha_surfL[9]+0.1767766952966367*alpha_surfL[8]-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5])+0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2])-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16]))+0.1767766952966367*alpha_surfL[13]-0.1767766952966367*alpha_surfL[12]+0.1767766952966367*alpha_surfL[10]-0.1767766952966367*alpha_surfL[9]+0.1767766952966367*alpha_surfL[8]-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6])+0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2])-0.1767766952966367*alpha_surfL[1]+0.1767766952966367*alpha_surfL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13])-0.1767766952966367*alpha_surfL[12]+0.1767766952966367*alpha_surfL[10]-0.1767766952966367*alpha_surfL[9]+0.1767766952966367*alpha_surfL[8]-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2])+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16])-0.1767766952966367*alpha_surfL[13]+0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])-0.1767766952966367*alpha_surfL[9]+0.1767766952966367*alpha_surfL[8]-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6])+0.1767766952966367*alpha_surfL[5]-0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2])+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*alpha_surfL[17]+0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13])-0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])+0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5])+0.1767766952966367*alpha_surfL[4]-0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2])+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]))+0.1767766952966367*alpha_surfL[16]-0.1767766952966367*alpha_surfL[13]+0.1767766952966367*alpha_surfL[12]-0.1767766952966367*alpha_surfL[10]+0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])-0.1767766952966367*(alpha_surfL[7]+alpha_surfL[6])+0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4])-0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2])+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17])-0.1767766952966367*alpha_surfL[16]+0.1767766952966367*alpha_surfL[13]-0.1767766952966367*alpha_surfL[12]+0.1767766952966367*alpha_surfL[10]-0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4])+0.1767766952966367*alpha_surfL[3]-0.1767766952966367*alpha_surfL[2]+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*alpha_surfL[17]-0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13])+0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])-0.1767766952966367*(alpha_surfL[9]+alpha_surfL[8])+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*alpha_surfL[6]+0.1767766952966367*alpha_surfL[5]-0.1767766952966367*alpha_surfL[4]+0.1767766952966367*alpha_surfL[3]-0.1767766952966367*alpha_surfL[2]+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16])+0.1767766952966367*alpha_surfL[13]-0.1767766952966367*(alpha_surfL[12]+alpha_surfL[10])+0.1767766952966367*alpha_surfL[9]-0.1767766952966367*alpha_surfL[8]+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5])+0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3])-0.1767766952966367*alpha_surfL[2]+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13]))+0.1767766952966367*alpha_surfL[12]-0.1767766952966367*alpha_surfL[10]+0.1767766952966367*alpha_surfL[9]-0.1767766952966367*alpha_surfL[8]+0.1767766952966367*alpha_surfL[7]-0.1767766952966367*alpha_surfL[6]+0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3])-0.1767766952966367*alpha_surfL[2]+0.1767766952966367*(alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7]))+0.1767766952966367*alpha_surfL[6]-0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3])+0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16])+0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12])-0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5])-0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3])+0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*alpha_surfL[17]-0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12])+0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9])-0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*alpha_surfL[6]-0.1767766952966367*alpha_surfL[5]+0.1767766952966367*alpha_surfL[4]-0.1767766952966367*alpha_surfL[3]+0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17])-0.1767766952966367*alpha_surfL[16]+0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9])-0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7])+0.1767766952966367*(alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4])-0.1767766952966367*alpha_surfL[3]+0.1767766952966367*(alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[27] = 1.0; 
  else  
    sgn_alpha_surfL[27] = -1.0; 
  
  if (sgn_alpha_surfL[27] == sgn_alpha_surfL[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]))+0.1767766952966367*alpha_surfL[16]-0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9])+0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6])-0.1767766952966367*(alpha_surfL[5]+alpha_surfL[4])+0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[28] = 1.0; 
  else  
    sgn_alpha_surfL[28] = -1.0; 
  
  if (sgn_alpha_surfL[28] == sgn_alpha_surfL[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*alpha_surfL[20]-0.1767766952966367*alpha_surfL[17]+0.1767766952966367*(alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12])-0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9])+0.1767766952966367*(alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5])-0.1767766952966367*alpha_surfL[4]+0.1767766952966367*(alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[29] = 1.0; 
  else  
    sgn_alpha_surfL[29] = -1.0; 
  
  if (sgn_alpha_surfL[29] == sgn_alpha_surfL[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.1767766952966367*alpha_surfL[20])+0.1767766952966367*(alpha_surfL[17]+alpha_surfL[16])-0.1767766952966367*(alpha_surfL[13]+alpha_surfL[12])+0.1767766952966367*(alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6])-0.1767766952966367*alpha_surfL[5]+0.1767766952966367*(alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[30] = 1.0; 
  else  
    sgn_alpha_surfL[30] = -1.0; 
  
  if (sgn_alpha_surfL[30] == sgn_alpha_surfL[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.1767766952966367*(alpha_surfL[20]+alpha_surfL[17]+alpha_surfL[16]+alpha_surfL[13]+alpha_surfL[12]+alpha_surfL[10]+alpha_surfL[9]+alpha_surfL[8]+alpha_surfL[7]+alpha_surfL[6]+alpha_surfL[5]+alpha_surfL[4]+alpha_surfL[3]+alpha_surfL[2]+alpha_surfL[1]+alpha_surfL[0]) > 0.) 
    sgn_alpha_surfL[31] = 1.0; 
  else  
    sgn_alpha_surfL[31] = -1.0; 
  
  if (sgn_alpha_surfL[31] == sgn_alpha_surfL[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
