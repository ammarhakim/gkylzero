#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_tensor_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[3] = {0.0}; 
  vol_incr[2] = 6.708203932499369*coeff[0]*fSkin[0]; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.6708203932499369*coeff[0]*fSkin[2])+0.6708203932499369*coeff[0]*fEdge[2]-1.190784930203603*coeff[0]*fSkin[1]-1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.585502557353661*coeff[0]*fSkin[2])+0.7382874503707886*coeff[0]*fEdge[2]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]-2.056810333988042*coeff[0]*fSkin[0]+1.190784930203603*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-3.140625*coeff[0]*fSkin[2])-0.140625*coeff[0]*fEdge[2]-5.022775277112744*coeff[0]*fSkin[1]-0.3025768239224549*coeff[0]*fEdge[1]-3.773364712030896*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*fSkin[2]-1.25*coeff[0]*fSkin[1]+0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = (-3.75*coeff[0]*fSkin[2])+4.841229182759272*coeff[0]*fSkin[1]-3.354101966249685*coeff[0]*fSkin[0]; 

  } else { 

  edgeSurf_incr[0] = (-0.6708203932499369*coeff[0]*fSkin[2])+0.6708203932499369*coeff[0]*fEdge[2]+1.190784930203603*coeff[0]*fSkin[1]+1.190784930203603*coeff[0]*fEdge[1]-0.9375*coeff[0]*fSkin[0]+0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 1.585502557353661*coeff[0]*fSkin[2]-0.7382874503707886*coeff[0]*fEdge[2]-2.671875*coeff[0]*fSkin[1]-1.453125*coeff[0]*fEdge[1]+2.056810333988042*coeff[0]*fSkin[0]-1.190784930203603*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-3.140625*coeff[0]*fSkin[2])-0.140625*coeff[0]*fEdge[2]+5.022775277112744*coeff[0]*fSkin[1]+0.3025768239224549*coeff[0]*fEdge[1]-3.773364712030896*coeff[0]*fSkin[0]+0.4192627457812108*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = (-0.9682458365518543*coeff[0]*fSkin[2])-1.25*coeff[0]*fSkin[1]-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = (-3.75*coeff[0]*fSkin[2])-4.841229182759272*coeff[0]*fSkin[1]-3.354101966249685*coeff[0]*fSkin[0]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_tensor_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[3] = {0.0}; 
  vol_incr[1] = 4.743416490252569*fSkin[1]*coeff[2]+2.121320343559642*fSkin[0]*coeff[1]; 
  vol_incr[2] = 14.23024947075771*coeff[2]*fSkin[2]+10.60660171779821*fSkin[0]*coeff[2]+9.48683298050514*coeff[1]*fSkin[1]+4.743416490252569*coeff[0]*fSkin[0]; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[2]-0.4743416490252568*coeff[0]*fSkin[2]-0.5303300858899105*coeff[2]*fEdge[2]+0.4743416490252568*coeff[0]*fEdge[2]+0.9413981457120035*fSkin[1]*coeff[2]+0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]-0.8420120990817169*coeff[0]*fSkin[1]-0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.2487763020014165*coeff[2]*fSkin[2]-0.5188111786213743*coeff[1]*fSkin[2]-1.121119609893386*coeff[0]*fSkin[2]-1.588341005085966*coeff[2]*fEdge[2]-0.5188111786213743*coeff[1]*fEdge[2]+0.5220480626221115*coeff[0]*fEdge[2]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]+0.5990715472712751*fSkin[0]*coeff[2]-1.96837794103419*fEdge[0]*coeff[2]-0.7463289060042488*coeff[1]*fSkin[1]-1.889300930982805*coeff[0]*fSkin[1]+0.7463289060042488*coeff[1]*fEdge[1]-1.027514541411701*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-1.454384534777511*coeff[0]*fSkin[0]+0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.40820177054373*coeff[2]*fSkin[2])-2.009347054626824*coeff[1]*fSkin[2]-2.220757234663999*coeff[0]*fSkin[2]-3.779910015670014*coeff[2]*fEdge[2]-2.009347054626824*coeff[1]*fEdge[2]-0.0994368911043575*coeff[0]*fEdge[2]-1.626614282316952*fSkin[1]*coeff[2]+5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.308931947855521*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]-3.551638458822559*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]-0.2139541240254559*coeff[0]*fEdge[1]-2.053959590644372*fSkin[0]*coeff[1]-2.053959590644372*fEdge[0]*coeff[1]-2.668171775767069*coeff[0]*fSkin[0]+0.2964635306407852*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 1.530931089239486*coeff[2]*fSkin[2]-1.185854122563142*coeff[1]*fSkin[2]+0.6846531968814573*coeff[0]*fSkin[2]-1.976423537605237*fSkin[1]*coeff[2]+1.369306393762915*fSkin[0]*coeff[2]+1.530931089239486*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = (-5.929270612815711*coeff[2]*fSkin[2])+4.592793267718456*coeff[1]*fSkin[2]-2.651650429449552*coeff[0]*fSkin[2]+7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]+3.423265984407287*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 

  } else { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[2]-0.4743416490252568*coeff[0]*fSkin[2]-0.5303300858899105*coeff[2]*fEdge[2]+0.4743416490252568*coeff[0]*fEdge[2]-0.9413981457120035*fSkin[1]*coeff[2]-0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]+0.8420120990817169*coeff[0]*fSkin[1]+0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.2487763020014165*coeff[2]*fSkin[2])-0.5188111786213743*coeff[1]*fSkin[2]+1.121119609893386*coeff[0]*fSkin[2]+1.588341005085966*coeff[2]*fEdge[2]-0.5188111786213743*coeff[1]*fEdge[2]-0.5220480626221115*coeff[0]*fEdge[2]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]-0.5990715472712751*fSkin[0]*coeff[2]+1.96837794103419*fEdge[0]*coeff[2]+0.7463289060042488*coeff[1]*fSkin[1]-1.889300930982805*coeff[0]*fSkin[1]-0.7463289060042488*coeff[1]*fEdge[1]-1.027514541411701*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+1.454384534777511*coeff[0]*fSkin[0]-0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.40820177054373*coeff[2]*fSkin[2])+2.009347054626824*coeff[1]*fSkin[2]-2.220757234663999*coeff[0]*fSkin[2]-3.779910015670014*coeff[2]*fEdge[2]+2.009347054626824*coeff[1]*fEdge[2]-0.0994368911043575*coeff[0]*fEdge[2]+1.626614282316952*fSkin[1]*coeff[2]-5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.308931947855521*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]+3.551638458822559*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]+0.2139541240254559*coeff[0]*fEdge[1]+2.053959590644372*fSkin[0]*coeff[1]+2.053959590644372*fEdge[0]*coeff[1]-2.668171775767069*coeff[0]*fSkin[0]+0.2964635306407852*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = (-1.530931089239486*coeff[2]*fSkin[2])-1.185854122563142*coeff[1]*fSkin[2]-0.6846531968814573*coeff[0]*fSkin[2]-1.976423537605237*fSkin[1]*coeff[2]-1.369306393762915*fSkin[0]*coeff[2]-1.530931089239486*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = (-5.929270612815711*coeff[2]*fSkin[2])-4.592793267718456*coeff[1]*fSkin[2]-2.651650429449552*coeff[0]*fSkin[2]-7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]-3.423265984407287*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  }

  return 0.;
}

