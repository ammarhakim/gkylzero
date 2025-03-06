#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double vol_incr[6] = {0.0}; 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[0]*qSkin[1])-0.5412658773652741*coeff[0]*qEdge[1]-0.5625*coeff[0]*qSkin[0]+0.5625*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*qSkin[1])-0.4375*coeff[0]*qEdge[1]-1.4072912811497125*coeff[0]*qSkin[0]+0.5412658773652739*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = -(0.5412658773652741*coeff[0]*qSkin[3])-0.5412658773652741*coeff[0]*qEdge[3]-0.5625*coeff[0]*qSkin[2]+0.5625*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(1.4375*coeff[0]*qSkin[3])-0.4375*coeff[0]*qEdge[3]-1.4072912811497125*coeff[0]*qSkin[2]+0.5412658773652739*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = -(0.5412658773652742*coeff[0]*qSkin[5])-0.5412658773652742*coeff[0]*qEdge[5]-0.5625*coeff[0]*qSkin[4]+0.5625*coeff[0]*qEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*qSkin[5])-0.4375*coeff[0]*qEdge[5]-1.4072912811497127*coeff[0]*qSkin[4]+0.5412658773652742*coeff[0]*qEdge[4]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*qSkin[0]-1.0*coeff[0]*qSkin[1]; 
  boundSurf_incr[3] = 0.8660254037844386*coeff[0]*qSkin[2]-1.0*coeff[0]*qSkin[3]; 
  boundSurf_incr[5] = 0.8660254037844387*coeff[0]*qSkin[4]-1.0*coeff[0]*qSkin[5]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*qSkin[1]+0.5412658773652741*coeff[0]*qEdge[1]-0.5625*coeff[0]*qSkin[0]+0.5625*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*qSkin[1])-0.4375*coeff[0]*qEdge[1]+1.4072912811497125*coeff[0]*qSkin[0]-0.5412658773652739*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*qSkin[3]+0.5412658773652741*coeff[0]*qEdge[3]-0.5625*coeff[0]*qSkin[2]+0.5625*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(1.4375*coeff[0]*qSkin[3])-0.4375*coeff[0]*qEdge[3]+1.4072912811497125*coeff[0]*qSkin[2]-0.5412658773652739*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = 0.5412658773652742*coeff[0]*qSkin[5]+0.5412658773652742*coeff[0]*qEdge[5]-0.5625*coeff[0]*qSkin[4]+0.5625*coeff[0]*qEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*qSkin[5])-0.4375*coeff[0]*qEdge[5]+1.4072912811497127*coeff[0]*qSkin[4]-0.5412658773652742*coeff[0]*qEdge[4]; 

  boundSurf_incr[1] = -(1.0*coeff[0]*qSkin[1])-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = -(1.0*coeff[0]*qSkin[3])-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[5] = -(1.0*coeff[0]*qSkin[5])-0.8660254037844387*coeff[0]*qSkin[4]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fSkin[6];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.7071067811865476*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1]); 
  fSkin[2] = 0.7071067811865476*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2]); 
  fSkin[3] = 0.7071067811865476*(jacobgeo_inv[0]*qSkin[3]+jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[5]+21.213203435596427*jacobgeo_inv[0]*qSkin[4]); 
  fSkin[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_inv[0]*qSkin[5]+21.21320343559643*jacobgeo_inv[1]*qSkin[4]); 

  double fEdge[6];
  fEdge[0] = 0.7071067811865476*(jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.7071067811865476*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1]); 
  fEdge[2] = 0.7071067811865476*(jacobgeo_inv[1]*qEdge[3]+jacobgeo_inv[0]*qEdge[2]); 
  fEdge[3] = 0.7071067811865476*(jacobgeo_inv[0]*qEdge[3]+jacobgeo_inv[1]*qEdge[2]); 
  fEdge[4] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qEdge[5]+21.213203435596427*jacobgeo_inv[0]*qEdge[4]); 
  fEdge[5] = 0.03333333333333333*(21.213203435596427*jacobgeo_inv[0]*qEdge[5]+21.21320343559643*jacobgeo_inv[1]*qEdge[4]); 

  double vol_incr[6] = {0.0}; 
  vol_incr[1] = 1.5*fSkin[4]*coeff[5]+1.5*fSkin[2]*coeff[3]+1.5*fSkin[0]*coeff[1]; 
  vol_incr[3] = 1.3416407864998738*fSkin[2]*coeff[5]+1.3416407864998738*coeff[3]*fSkin[4]+1.5*fSkin[0]*coeff[3]+1.5*coeff[1]*fSkin[2]; 
  vol_incr[5] = 0.9583148474999099*fSkin[4]*coeff[5]+1.5*fSkin[0]*coeff[5]+1.5*coeff[1]*fSkin[4]+1.3416407864998738*fSkin[2]*coeff[3]; 

  double edgeSurf_incr[6] = {0.0}; 
  double boundSurf_incr[6] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.2706329386826371*coeff[4]*fSkin[5])-0.2706329386826371*coeff[4]*fEdge[5]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]-0.27063293868263705*coeff[2]*fSkin[3]-0.27063293868263705*coeff[2]*fEdge[3]-0.28125*coeff[2]*fSkin[2]+0.28125*coeff[2]*fEdge[2]-0.27063293868263705*coeff[0]*fSkin[1]-0.27063293868263705*coeff[0]*fEdge[1]-0.28125*coeff[0]*fSkin[0]+0.28125*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(0.4330127018922193*coeff[5]*fSkin[5])-0.7187500000000001*coeff[4]*fSkin[5]+0.4330127018922193*coeff[5]*fEdge[5]-0.21875*coeff[4]*fEdge[5]-0.375*fSkin[4]*coeff[5]-0.375*fEdge[4]*coeff[5]-0.7036456405748562*coeff[4]*fSkin[4]+0.27063293868263694*coeff[4]*fEdge[4]-0.4330127018922193*coeff[3]*fSkin[3]-0.71875*coeff[2]*fSkin[3]+0.4330127018922193*coeff[3]*fEdge[3]-0.21875*coeff[2]*fEdge[3]-0.375*fSkin[2]*coeff[3]-0.375*fEdge[2]*coeff[3]-0.7036456405748562*coeff[2]*fSkin[2]+0.27063293868263694*coeff[2]*fEdge[2]-0.4330127018922193*coeff[1]*fSkin[1]-0.71875*coeff[0]*fSkin[1]+0.4330127018922193*coeff[1]*fEdge[1]-0.21875*coeff[0]*fEdge[1]-0.375*fSkin[0]*coeff[1]-0.375*fEdge[0]*coeff[1]-0.7036456405748562*coeff[0]*fSkin[0]+0.27063293868263694*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.24206145913796356*coeff[2]*fSkin[5])-0.24206145913796356*coeff[2]*fEdge[5]-0.25155764746872633*coeff[2]*fSkin[4]+0.25155764746872633*coeff[2]*fEdge[4]-0.24206145913796356*fSkin[3]*coeff[4]-0.24206145913796356*fEdge[3]*coeff[4]-0.25155764746872633*fSkin[2]*coeff[4]+0.25155764746872633*fEdge[2]*coeff[4]-0.27063293868263705*coeff[0]*fSkin[3]-0.27063293868263705*coeff[0]*fEdge[3]-0.28125*coeff[0]*fSkin[2]+0.28125*coeff[0]*fEdge[2]-0.27063293868263705*fSkin[1]*coeff[2]-0.27063293868263705*fEdge[1]*coeff[2]-0.28125*fSkin[0]*coeff[2]+0.28125*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = -(0.3872983346207417*coeff[3]*fSkin[5])-0.6428695435311895*coeff[2]*fSkin[5]+0.3872983346207417*coeff[3]*fEdge[5]-0.1956559480312315*coeff[2]*fEdge[5]-0.3872983346207417*fSkin[3]*coeff[5]+0.3872983346207417*fEdge[3]*coeff[5]-0.33541019662496846*fSkin[2]*coeff[5]-0.33541019662496846*fEdge[2]*coeff[5]-0.33541019662496846*coeff[3]*fSkin[4]-0.6293597937587051*coeff[2]*fSkin[4]-0.33541019662496846*coeff[3]*fEdge[4]+0.24206145913796343*coeff[2]*fEdge[4]-0.6428695435311895*fSkin[3]*coeff[4]-0.19565594803123162*fEdge[3]*coeff[4]-0.6293597937587051*fSkin[2]*coeff[4]+0.24206145913796343*fEdge[2]*coeff[4]-0.4330127018922193*coeff[1]*fSkin[3]-0.71875*coeff[0]*fSkin[3]+0.4330127018922193*coeff[1]*fEdge[3]-0.21875*coeff[0]*fEdge[3]-0.4330127018922193*fSkin[1]*coeff[3]+0.4330127018922193*fEdge[1]*coeff[3]-0.375*fSkin[0]*coeff[3]-0.375*fEdge[0]*coeff[3]-0.375*coeff[1]*fSkin[2]-0.7036456405748562*coeff[0]*fSkin[2]-0.375*coeff[1]*fEdge[2]+0.27063293868263694*coeff[0]*fEdge[2]-0.71875*fSkin[1]*coeff[2]-0.21875*fEdge[1]*coeff[2]-0.7036456405748562*fSkin[0]*coeff[2]+0.27063293868263694*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = -(0.17290104224140254*coeff[4]*fSkin[5])-0.2706329386826371*coeff[0]*fSkin[5]-0.17290104224140254*coeff[4]*fEdge[5]-0.2706329386826371*coeff[0]*fEdge[5]-0.1796840339062331*coeff[4]*fSkin[4]-0.28125*coeff[0]*fSkin[4]+0.1796840339062331*coeff[4]*fEdge[4]+0.28125*coeff[0]*fEdge[4]-0.27063293868263705*fSkin[1]*coeff[4]-0.27063293868263705*fEdge[1]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]-0.24206145913796356*coeff[2]*fSkin[3]-0.24206145913796356*coeff[2]*fEdge[3]-0.25155764746872633*coeff[2]*fSkin[2]+0.25155764746872633*coeff[2]*fEdge[2]; 
  edgeSurf_incr[5] = -(0.27664166758624403*coeff[5]*fSkin[5])-0.4591925310937069*coeff[4]*fSkin[5]-0.4330127018922193*coeff[1]*fSkin[5]-0.71875*coeff[0]*fSkin[5]+0.27664166758624403*coeff[5]*fEdge[5]-0.13975424859373692*coeff[4]*fEdge[5]+0.4330127018922193*coeff[1]*fEdge[5]-0.21875*coeff[0]*fEdge[5]-0.23957871187497748*fSkin[4]*coeff[5]-0.23957871187497748*fEdge[4]*coeff[5]-0.4330127018922193*fSkin[1]*coeff[5]+0.4330127018922193*fEdge[1]*coeff[5]-0.375*fSkin[0]*coeff[5]-0.375*fEdge[0]*coeff[5]-0.44954270982764666*coeff[4]*fSkin[4]-0.375*coeff[1]*fSkin[4]-0.7036456405748563*coeff[0]*fSkin[4]+0.1729010422414026*coeff[4]*fEdge[4]-0.375*coeff[1]*fEdge[4]+0.2706329386826371*coeff[0]*fEdge[4]-0.7187500000000001*fSkin[1]*coeff[4]-0.21875*fEdge[1]*coeff[4]-0.7036456405748563*fSkin[0]*coeff[4]+0.2706329386826371*fEdge[0]*coeff[4]-0.3872983346207417*coeff[3]*fSkin[3]-0.6428695435311895*coeff[2]*fSkin[3]+0.3872983346207417*coeff[3]*fEdge[3]-0.1956559480312315*coeff[2]*fEdge[3]-0.33541019662496846*fSkin[2]*coeff[3]-0.33541019662496846*fEdge[2]*coeff[3]-0.6293597937587053*coeff[2]*fSkin[2]+0.24206145913796356*coeff[2]*fEdge[2]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[5]*fSkin[5]-0.5000000000000001*coeff[4]*fSkin[5]-0.75*fSkin[4]*coeff[5]+0.4330127018922193*coeff[4]*fSkin[4]+0.8660254037844386*coeff[3]*fSkin[3]-0.5*coeff[2]*fSkin[3]-0.75*fSkin[2]*coeff[3]+0.4330127018922193*coeff[2]*fSkin[2]+0.8660254037844386*coeff[1]*fSkin[1]-0.5*coeff[0]*fSkin[1]-0.75*fSkin[0]*coeff[1]+0.4330127018922193*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = 0.7745966692414834*coeff[3]*fSkin[5]-0.44721359549995804*coeff[2]*fSkin[5]+0.7745966692414834*fSkin[3]*coeff[5]-0.6708203932499369*fSkin[2]*coeff[5]-0.6708203932499369*coeff[3]*fSkin[4]+0.38729833462074165*coeff[2]*fSkin[4]-0.4472135954999579*fSkin[3]*coeff[4]+0.38729833462074165*fSkin[2]*coeff[4]+0.8660254037844386*coeff[1]*fSkin[3]-0.5*coeff[0]*fSkin[3]+0.8660254037844386*fSkin[1]*coeff[3]-0.75*fSkin[0]*coeff[3]-0.75*coeff[1]*fSkin[2]+0.4330127018922193*coeff[0]*fSkin[2]-0.5*fSkin[1]*coeff[2]+0.4330127018922193*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = 0.5532833351724881*coeff[5]*fSkin[5]-0.31943828249996997*coeff[4]*fSkin[5]+0.8660254037844386*coeff[1]*fSkin[5]-0.5*coeff[0]*fSkin[5]-0.47915742374995496*fSkin[4]*coeff[5]+0.8660254037844386*fSkin[1]*coeff[5]-0.75*fSkin[0]*coeff[5]+0.27664166758624403*coeff[4]*fSkin[4]-0.75*coeff[1]*fSkin[4]+0.43301270189221935*coeff[0]*fSkin[4]-0.5000000000000001*fSkin[1]*coeff[4]+0.43301270189221935*fSkin[0]*coeff[4]+0.7745966692414834*coeff[3]*fSkin[3]-0.44721359549995804*coeff[2]*fSkin[3]-0.6708203932499369*fSkin[2]*coeff[3]+0.3872983346207417*coeff[2]*fSkin[2]; 

  } else { 

  edgeSurf_incr[0] = 0.2706329386826371*coeff[4]*fSkin[5]+0.2706329386826371*coeff[4]*fEdge[5]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]+0.27063293868263705*coeff[2]*fSkin[3]+0.27063293868263705*coeff[2]*fEdge[3]-0.28125*coeff[2]*fSkin[2]+0.28125*coeff[2]*fEdge[2]+0.27063293868263705*coeff[0]*fSkin[1]+0.27063293868263705*coeff[0]*fEdge[1]-0.28125*coeff[0]*fSkin[0]+0.28125*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.4330127018922193*coeff[5]*fSkin[5]-0.7187500000000001*coeff[4]*fSkin[5]-0.4330127018922193*coeff[5]*fEdge[5]-0.21875*coeff[4]*fEdge[5]-0.375*fSkin[4]*coeff[5]-0.375*fEdge[4]*coeff[5]+0.7036456405748562*coeff[4]*fSkin[4]-0.27063293868263694*coeff[4]*fEdge[4]+0.4330127018922193*coeff[3]*fSkin[3]-0.71875*coeff[2]*fSkin[3]-0.4330127018922193*coeff[3]*fEdge[3]-0.21875*coeff[2]*fEdge[3]-0.375*fSkin[2]*coeff[3]-0.375*fEdge[2]*coeff[3]+0.7036456405748562*coeff[2]*fSkin[2]-0.27063293868263694*coeff[2]*fEdge[2]+0.4330127018922193*coeff[1]*fSkin[1]-0.71875*coeff[0]*fSkin[1]-0.4330127018922193*coeff[1]*fEdge[1]-0.21875*coeff[0]*fEdge[1]-0.375*fSkin[0]*coeff[1]-0.375*fEdge[0]*coeff[1]+0.7036456405748562*coeff[0]*fSkin[0]-0.27063293868263694*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.24206145913796356*coeff[2]*fSkin[5]+0.24206145913796356*coeff[2]*fEdge[5]-0.25155764746872633*coeff[2]*fSkin[4]+0.25155764746872633*coeff[2]*fEdge[4]+0.24206145913796356*fSkin[3]*coeff[4]+0.24206145913796356*fEdge[3]*coeff[4]-0.25155764746872633*fSkin[2]*coeff[4]+0.25155764746872633*fEdge[2]*coeff[4]+0.27063293868263705*coeff[0]*fSkin[3]+0.27063293868263705*coeff[0]*fEdge[3]-0.28125*coeff[0]*fSkin[2]+0.28125*coeff[0]*fEdge[2]+0.27063293868263705*fSkin[1]*coeff[2]+0.27063293868263705*fEdge[1]*coeff[2]-0.28125*fSkin[0]*coeff[2]+0.28125*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = 0.3872983346207417*coeff[3]*fSkin[5]-0.6428695435311895*coeff[2]*fSkin[5]-0.3872983346207417*coeff[3]*fEdge[5]-0.1956559480312315*coeff[2]*fEdge[5]+0.3872983346207417*fSkin[3]*coeff[5]-0.3872983346207417*fEdge[3]*coeff[5]-0.33541019662496846*fSkin[2]*coeff[5]-0.33541019662496846*fEdge[2]*coeff[5]-0.33541019662496846*coeff[3]*fSkin[4]+0.6293597937587051*coeff[2]*fSkin[4]-0.33541019662496846*coeff[3]*fEdge[4]-0.24206145913796343*coeff[2]*fEdge[4]-0.6428695435311895*fSkin[3]*coeff[4]-0.19565594803123162*fEdge[3]*coeff[4]+0.6293597937587051*fSkin[2]*coeff[4]-0.24206145913796343*fEdge[2]*coeff[4]+0.4330127018922193*coeff[1]*fSkin[3]-0.71875*coeff[0]*fSkin[3]-0.4330127018922193*coeff[1]*fEdge[3]-0.21875*coeff[0]*fEdge[3]+0.4330127018922193*fSkin[1]*coeff[3]-0.4330127018922193*fEdge[1]*coeff[3]-0.375*fSkin[0]*coeff[3]-0.375*fEdge[0]*coeff[3]-0.375*coeff[1]*fSkin[2]+0.7036456405748562*coeff[0]*fSkin[2]-0.375*coeff[1]*fEdge[2]-0.27063293868263694*coeff[0]*fEdge[2]-0.71875*fSkin[1]*coeff[2]-0.21875*fEdge[1]*coeff[2]+0.7036456405748562*fSkin[0]*coeff[2]-0.27063293868263694*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = 0.17290104224140254*coeff[4]*fSkin[5]+0.2706329386826371*coeff[0]*fSkin[5]+0.17290104224140254*coeff[4]*fEdge[5]+0.2706329386826371*coeff[0]*fEdge[5]-0.1796840339062331*coeff[4]*fSkin[4]-0.28125*coeff[0]*fSkin[4]+0.1796840339062331*coeff[4]*fEdge[4]+0.28125*coeff[0]*fEdge[4]+0.27063293868263705*fSkin[1]*coeff[4]+0.27063293868263705*fEdge[1]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]+0.24206145913796356*coeff[2]*fSkin[3]+0.24206145913796356*coeff[2]*fEdge[3]-0.25155764746872633*coeff[2]*fSkin[2]+0.25155764746872633*coeff[2]*fEdge[2]; 
  edgeSurf_incr[5] = 0.27664166758624403*coeff[5]*fSkin[5]-0.4591925310937069*coeff[4]*fSkin[5]+0.4330127018922193*coeff[1]*fSkin[5]-0.71875*coeff[0]*fSkin[5]-0.27664166758624403*coeff[5]*fEdge[5]-0.13975424859373692*coeff[4]*fEdge[5]-0.4330127018922193*coeff[1]*fEdge[5]-0.21875*coeff[0]*fEdge[5]-0.23957871187497748*fSkin[4]*coeff[5]-0.23957871187497748*fEdge[4]*coeff[5]+0.4330127018922193*fSkin[1]*coeff[5]-0.4330127018922193*fEdge[1]*coeff[5]-0.375*fSkin[0]*coeff[5]-0.375*fEdge[0]*coeff[5]+0.44954270982764666*coeff[4]*fSkin[4]-0.375*coeff[1]*fSkin[4]+0.7036456405748563*coeff[0]*fSkin[4]-0.1729010422414026*coeff[4]*fEdge[4]-0.375*coeff[1]*fEdge[4]-0.2706329386826371*coeff[0]*fEdge[4]-0.7187500000000001*fSkin[1]*coeff[4]-0.21875*fEdge[1]*coeff[4]+0.7036456405748563*fSkin[0]*coeff[4]-0.2706329386826371*fEdge[0]*coeff[4]+0.3872983346207417*coeff[3]*fSkin[3]-0.6428695435311895*coeff[2]*fSkin[3]-0.3872983346207417*coeff[3]*fEdge[3]-0.1956559480312315*coeff[2]*fEdge[3]-0.33541019662496846*fSkin[2]*coeff[3]-0.33541019662496846*fEdge[2]*coeff[3]+0.6293597937587053*coeff[2]*fSkin[2]-0.24206145913796356*coeff[2]*fEdge[2]; 

  boundSurf_incr[1] = -(0.8660254037844386*coeff[5]*fSkin[5])-0.5000000000000001*coeff[4]*fSkin[5]-0.75*fSkin[4]*coeff[5]-0.4330127018922193*coeff[4]*fSkin[4]-0.8660254037844386*coeff[3]*fSkin[3]-0.5*coeff[2]*fSkin[3]-0.75*fSkin[2]*coeff[3]-0.4330127018922193*coeff[2]*fSkin[2]-0.8660254037844386*coeff[1]*fSkin[1]-0.5*coeff[0]*fSkin[1]-0.75*fSkin[0]*coeff[1]-0.4330127018922193*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(0.7745966692414834*coeff[3]*fSkin[5])-0.44721359549995804*coeff[2]*fSkin[5]-0.7745966692414834*fSkin[3]*coeff[5]-0.6708203932499369*fSkin[2]*coeff[5]-0.6708203932499369*coeff[3]*fSkin[4]-0.38729833462074165*coeff[2]*fSkin[4]-0.4472135954999579*fSkin[3]*coeff[4]-0.38729833462074165*fSkin[2]*coeff[4]-0.8660254037844386*coeff[1]*fSkin[3]-0.5*coeff[0]*fSkin[3]-0.8660254037844386*fSkin[1]*coeff[3]-0.75*fSkin[0]*coeff[3]-0.75*coeff[1]*fSkin[2]-0.4330127018922193*coeff[0]*fSkin[2]-0.5*fSkin[1]*coeff[2]-0.4330127018922193*fSkin[0]*coeff[2]; 
  boundSurf_incr[5] = -(0.5532833351724881*coeff[5]*fSkin[5])-0.31943828249996997*coeff[4]*fSkin[5]-0.8660254037844386*coeff[1]*fSkin[5]-0.5*coeff[0]*fSkin[5]-0.47915742374995496*fSkin[4]*coeff[5]-0.8660254037844386*fSkin[1]*coeff[5]-0.75*fSkin[0]*coeff[5]-0.27664166758624403*coeff[4]*fSkin[4]-0.75*coeff[1]*fSkin[4]-0.43301270189221935*coeff[0]*fSkin[4]-0.5000000000000001*fSkin[1]*coeff[4]-0.43301270189221935*fSkin[0]*coeff[4]-0.7745966692414834*coeff[3]*fSkin[3]-0.44721359549995804*coeff[2]*fSkin[3]-0.6708203932499369*fSkin[2]*coeff[3]-0.3872983346207417*coeff[2]*fSkin[2]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 

  }

  return 0.;
}

