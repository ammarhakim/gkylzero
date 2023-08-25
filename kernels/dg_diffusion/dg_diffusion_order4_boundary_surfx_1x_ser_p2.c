#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[3] = {0.0}; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[2]-6.708203932499369*coeff[0]*fEdge[2]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.16059535957086*coeff[0]*fSkin[2]-9.077304717673634*coeff[0]*fEdge[2]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 20.34375*coeff[0]*fSkin[2]-0.65625*coeff[0]*fEdge[2]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[2]; 
  boundSurf_incr[2] = 19.6875*coeff[0]*fSkin[2]-10.89276566120836*coeff[0]*fSkin[1]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[2]-6.708203932499369*coeff[0]*fEdge[2]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-14.16059535957086*coeff[0]*fSkin[2])+9.077304717673634*coeff[0]*fEdge[2]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 20.34375*coeff[0]*fSkin[2]-0.65625*coeff[0]*fEdge[2]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[2]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 19.6875*coeff[0]*fSkin[2]+10.89276566120836*coeff[0]*fSkin[1]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order4_boundary_surfx_1x_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[3] = {0.0}; 
  vol_incr[2] = -31.81980515339464*fSkin[0]*coeff[2]; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-5.303300858899105*coeff[2]*fSkin[2])-1.797214641813825*coeff[1]*fSkin[2]+4.743416490252569*coeff[0]*fSkin[2]+5.303300858899105*coeff[2]*fEdge[2]-1.797214641813825*coeff[1]*fEdge[2]-4.743416490252569*coeff[0]*fEdge[2]-6.418623720763661*fSkin[1]*coeff[2]-6.418623720763661*fEdge[1]*coeff[2]-3.705794133009818*fSkin[0]*coeff[2]+3.705794133009818*fEdge[0]*coeff[2]-0.9943689110435817*coeff[1]*fSkin[1]+5.74099158464807*coeff[0]*fSkin[1]+0.9943689110435817*coeff[1]*fEdge[1]+5.74099158464807*coeff[0]*fEdge[1]+3.31456303681194*coeff[0]*fSkin[0]-3.31456303681194*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-11.19493359006374*coeff[2]*fSkin[2])-3.112867071728247*coeff[1]*fSkin[2]+10.01305300439132*coeff[0]*fSkin[2]+7.176239480810091*coeff[2]*fEdge[2]-3.112867071728247*coeff[1]*fEdge[2]-6.418623720763665*coeff[0]*fEdge[2]-12.2291206389324*fSkin[1]*coeff[2]-10.00564415912651*fEdge[1]*coeff[2]-6.418623720763661*fSkin[0]*coeff[2]+6.418623720763661*fEdge[0]*coeff[2]-1.72229747539442*coeff[1]*fSkin[1]+10.9380580214794*coeff[0]*fSkin[1]+1.72229747539442*coeff[1]*fEdge[1]+8.949320199392238*coeff[0]*fEdge[1]+5.74099158464807*coeff[0]*fSkin[0]-5.74099158464807*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-34.53800131965151*coeff[2]*fSkin[2])-11.53939308514262*coeff[1]*fSkin[2]+14.38520357976382*coeff[0]*fSkin[2]+3.409330602369041*coeff[2]*fEdge[2]-0.516689242618324*coeff[1]*fEdge[2]-0.4640388251536773*coeff[0]*fEdge[2]-42.4833377263957*fSkin[1]*coeff[2]-11.48198316929614*fEdge[1]*coeff[2]-26.18504799081432*fSkin[0]*coeff[2]+10.27514541411701*fEdge[0]*coeff[2]-14.89729241469947*coeff[1]*fSkin[1]+11.0400327997135*coeff[0]*fSkin[1]-4.669300607592371*coeff[1]*fEdge[1]+3.337684334797107*coeff[0]*fEdge[1]-9.756308055560766*fSkin[0]*coeff[1]+5.648388874272021*fEdge[0]*coeff[1]+2.964635306407856*coeff[0]*fSkin[0]-2.964635306407856*coeff[0]*fEdge[0]; 

  boundSurf_incr[0] = 3.59442928362765*coeff[1]*fSkin[2]-1.988737822087164*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 4.018694109253648*coeff[2]*fSkin[2]-6.225734143456494*coeff[1]*fSkin[2]-3.59442928362765*coeff[0]*fSkin[2]-2.223476479805891*fSkin[1]*coeff[2]+3.444594950788841*coeff[1]*fSkin[1]+1.988737822087164*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-31.12867071728247*coeff[2]*fSkin[2])+12.05608232776095*coeff[1]*fSkin[2]+13.92116475461015*coeff[0]*fSkin[2]+31.00135455709956*fSkin[1]*coeff[2]-15.90990257669731*fSkin[0]*coeff[2]-10.2279918071071*coeff[1]*fSkin[1]-7.702348464916393*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]; 

  } else { 

  edgeSurf_incr[0] = (-5.303300858899105*coeff[2]*fSkin[2])+1.797214641813825*coeff[1]*fSkin[2]+4.743416490252569*coeff[0]*fSkin[2]+5.303300858899105*coeff[2]*fEdge[2]+1.797214641813825*coeff[1]*fEdge[2]-4.743416490252569*coeff[0]*fEdge[2]+6.418623720763661*fSkin[1]*coeff[2]+6.418623720763661*fEdge[1]*coeff[2]-3.705794133009818*fSkin[0]*coeff[2]+3.705794133009818*fEdge[0]*coeff[2]-0.9943689110435817*coeff[1]*fSkin[1]-5.74099158464807*coeff[0]*fSkin[1]+0.9943689110435817*coeff[1]*fEdge[1]-5.74099158464807*coeff[0]*fEdge[1]+3.31456303681194*coeff[0]*fSkin[0]-3.31456303681194*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 11.19493359006374*coeff[2]*fSkin[2]-3.112867071728247*coeff[1]*fSkin[2]-10.01305300439132*coeff[0]*fSkin[2]-7.176239480810091*coeff[2]*fEdge[2]-3.112867071728247*coeff[1]*fEdge[2]+6.418623720763665*coeff[0]*fEdge[2]-12.2291206389324*fSkin[1]*coeff[2]-10.00564415912651*fEdge[1]*coeff[2]+6.418623720763661*fSkin[0]*coeff[2]-6.418623720763661*fEdge[0]*coeff[2]+1.72229747539442*coeff[1]*fSkin[1]+10.9380580214794*coeff[0]*fSkin[1]-1.72229747539442*coeff[1]*fEdge[1]+8.949320199392238*coeff[0]*fEdge[1]-5.74099158464807*coeff[0]*fSkin[0]+5.74099158464807*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-34.53800131965151*coeff[2]*fSkin[2])+11.53939308514262*coeff[1]*fSkin[2]+14.38520357976382*coeff[0]*fSkin[2]+3.409330602369041*coeff[2]*fEdge[2]+0.516689242618324*coeff[1]*fEdge[2]-0.4640388251536773*coeff[0]*fEdge[2]+42.4833377263957*fSkin[1]*coeff[2]+11.48198316929614*fEdge[1]*coeff[2]-26.18504799081432*fSkin[0]*coeff[2]+10.27514541411701*fEdge[0]*coeff[2]-14.89729241469947*coeff[1]*fSkin[1]-11.0400327997135*coeff[0]*fSkin[1]-4.669300607592371*coeff[1]*fEdge[1]-3.337684334797107*coeff[0]*fEdge[1]+9.756308055560766*fSkin[0]*coeff[1]-5.648388874272021*fEdge[0]*coeff[1]+2.964635306407856*coeff[0]*fSkin[0]-2.964635306407856*coeff[0]*fEdge[0]; 

  boundSurf_incr[0] = (-3.59442928362765*coeff[1]*fSkin[2])-1.988737822087164*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = (-4.018694109253648*coeff[2]*fSkin[2])-6.225734143456494*coeff[1]*fSkin[2]+3.59442928362765*coeff[0]*fSkin[2]-2.223476479805891*fSkin[1]*coeff[2]-3.444594950788841*coeff[1]*fSkin[1]+1.988737822087164*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-31.12867071728247*coeff[2]*fSkin[2])-12.05608232776095*coeff[1]*fSkin[2]+13.92116475461015*coeff[0]*fSkin[2]-31.00135455709956*fSkin[1]*coeff[2]-15.90990257669731*fSkin[0]*coeff[2]-10.2279918071071*coeff[1]*fSkin[1]+7.702348464916393*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  }

  return 0.;
}

