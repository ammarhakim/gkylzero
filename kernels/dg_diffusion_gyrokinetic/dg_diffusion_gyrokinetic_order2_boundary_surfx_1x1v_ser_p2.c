#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double vol_incr[8] = {0.0}; 
  vol_incr[4] = 6.708203932499369*coeff[0]*qSkin[0]; 
  vol_incr[6] = 6.708203932499369*coeff[0]*qSkin[2]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.6708203932499369*coeff[0]*qSkin[4])+0.6708203932499369*coeff[0]*qEdge[4]-1.190784930203603*coeff[0]*qSkin[1]-1.190784930203603*coeff[0]*qEdge[1]-0.9375*coeff[0]*qSkin[0]+0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = (-1.585502557353661*coeff[0]*qSkin[4])+0.7382874503707886*coeff[0]*qEdge[4]-2.671875*coeff[0]*qSkin[1]-1.453125*coeff[0]*qEdge[1]-2.056810333988042*coeff[0]*qSkin[0]+1.190784930203603*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = (-0.6708203932499369*coeff[0]*qSkin[6])+0.6708203932499369*coeff[0]*qEdge[6]-1.190784930203603*coeff[0]*qSkin[3]-1.190784930203603*coeff[0]*qEdge[3]-0.9375*coeff[0]*qSkin[2]+0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = (-1.585502557353661*coeff[0]*qSkin[6])+0.7382874503707888*coeff[0]*qEdge[6]-2.671875*coeff[0]*qSkin[3]-1.453125*coeff[0]*qEdge[3]-2.056810333988042*coeff[0]*qSkin[2]+1.190784930203603*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = (-3.140625*coeff[0]*qSkin[4])-0.140625*coeff[0]*qEdge[4]-5.022775277112744*coeff[0]*qSkin[1]-0.3025768239224549*coeff[0]*qEdge[1]-3.773364712030896*coeff[0]*qSkin[0]+0.4192627457812108*coeff[0]*qEdge[0]; 
  edgeSurf_incr[5] = (-1.190784930203603*coeff[0]*qSkin[7])-1.190784930203603*coeff[0]*qEdge[7]-0.9375*coeff[0]*qSkin[5]+0.9375*coeff[0]*qEdge[5]; 
  edgeSurf_incr[6] = (-3.140625*coeff[0]*qSkin[6])-0.140625*coeff[0]*qEdge[6]-5.022775277112744*coeff[0]*qSkin[3]-0.3025768239224544*coeff[0]*qEdge[3]-3.773364712030894*coeff[0]*qSkin[2]+0.4192627457812105*coeff[0]*qEdge[2]; 
  edgeSurf_incr[7] = (-2.671875*coeff[0]*qSkin[7])-1.453125*coeff[0]*qEdge[7]-2.056810333988042*coeff[0]*qSkin[5]+1.190784930203603*coeff[0]*qEdge[5]; 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*qSkin[4]-1.25*coeff[0]*qSkin[1]+0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = 0.9682458365518543*coeff[0]*qSkin[6]-1.25*coeff[0]*qSkin[3]+0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[4] = (-3.75*coeff[0]*qSkin[4])+4.841229182759272*coeff[0]*qSkin[1]-3.354101966249685*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = (-3.75*coeff[0]*qSkin[6])+4.841229182759271*coeff[0]*qSkin[3]-3.354101966249684*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = 0.8660254037844387*coeff[0]*qSkin[5]-1.25*coeff[0]*qSkin[7]; 

  } else { 

  edgeSurf_incr[0] = (-0.6708203932499369*coeff[0]*qSkin[4])+0.6708203932499369*coeff[0]*qEdge[4]+1.190784930203603*coeff[0]*qSkin[1]+1.190784930203603*coeff[0]*qEdge[1]-0.9375*coeff[0]*qSkin[0]+0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 1.585502557353661*coeff[0]*qSkin[4]-0.7382874503707886*coeff[0]*qEdge[4]-2.671875*coeff[0]*qSkin[1]-1.453125*coeff[0]*qEdge[1]+2.056810333988042*coeff[0]*qSkin[0]-1.190784930203603*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = (-0.6708203932499369*coeff[0]*qSkin[6])+0.6708203932499369*coeff[0]*qEdge[6]+1.190784930203603*coeff[0]*qSkin[3]+1.190784930203603*coeff[0]*qEdge[3]-0.9375*coeff[0]*qSkin[2]+0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 1.585502557353661*coeff[0]*qSkin[6]-0.7382874503707888*coeff[0]*qEdge[6]-2.671875*coeff[0]*qSkin[3]-1.453125*coeff[0]*qEdge[3]+2.056810333988042*coeff[0]*qSkin[2]-1.190784930203603*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = (-3.140625*coeff[0]*qSkin[4])-0.140625*coeff[0]*qEdge[4]+5.022775277112744*coeff[0]*qSkin[1]+0.3025768239224549*coeff[0]*qEdge[1]-3.773364712030896*coeff[0]*qSkin[0]+0.4192627457812108*coeff[0]*qEdge[0]; 
  edgeSurf_incr[5] = 1.190784930203603*coeff[0]*qSkin[7]+1.190784930203603*coeff[0]*qEdge[7]-0.9375*coeff[0]*qSkin[5]+0.9375*coeff[0]*qEdge[5]; 
  edgeSurf_incr[6] = (-3.140625*coeff[0]*qSkin[6])-0.140625*coeff[0]*qEdge[6]+5.022775277112744*coeff[0]*qSkin[3]+0.3025768239224544*coeff[0]*qEdge[3]-3.773364712030894*coeff[0]*qSkin[2]+0.4192627457812105*coeff[0]*qEdge[2]; 
  edgeSurf_incr[7] = (-2.671875*coeff[0]*qSkin[7])-1.453125*coeff[0]*qEdge[7]+2.056810333988042*coeff[0]*qSkin[5]-1.190784930203603*coeff[0]*qEdge[5]; 

  boundSurf_incr[1] = (-0.9682458365518543*coeff[0]*qSkin[4])-1.25*coeff[0]*qSkin[1]-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = (-0.9682458365518543*coeff[0]*qSkin[6])-1.25*coeff[0]*qSkin[3]-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[4] = (-3.75*coeff[0]*qSkin[4])-4.841229182759272*coeff[0]*qSkin[1]-3.354101966249685*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = (-3.75*coeff[0]*qSkin[6])-4.841229182759271*coeff[0]*qSkin[3]-3.354101966249684*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = (-1.25*coeff[0]*qSkin[7])-0.8660254037844387*coeff[0]*qSkin[5]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fSkin[8];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[2]*qSkin[4]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qSkin[4]+qSkin[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1])); 
  fSkin[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qSkin[6]+21.21320343559643*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2])); 
  fSkin[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qSkin[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qSkin[3]+21.21320343559643*jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qSkin[4]+49.49747468305833*qSkin[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qSkin[1]); 
  fSkin[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[7]+21.21320343559643*jacobgeo_inv[0]*qSkin[5]); 
  fSkin[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*qSkin[6]+132.815661727072*jacobgeo_inv[1]*qSkin[3]+148.492424049175*jacobgeo_inv[2]*qSkin[2]); 
  fSkin[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qSkin[7]+21.21320343559643*jacobgeo_inv[1]*qSkin[5]); 

  double fEdge[8];
  fEdge[0] = 0.7071067811865476*(jacobgeo_inv[2]*qEdge[4]+jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qEdge[4]+qEdge[1]*jacobgeo_inv[2])+7.071067811865476*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1])); 
  fEdge[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qEdge[6]+21.21320343559643*(jacobgeo_inv[1]*qEdge[3]+jacobgeo_inv[0]*qEdge[2])); 
  fEdge[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qEdge[6]+(18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qEdge[3]+21.21320343559643*jacobgeo_inv[1]*qEdge[2]); 
  fEdge[4] = 0.01428571428571429*((31.62277660168381*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qEdge[4]+49.49747468305833*qEdge[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qEdge[1]); 
  fEdge[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qEdge[7]+21.21320343559643*jacobgeo_inv[0]*qEdge[5]); 
  fEdge[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.492424049175*jacobgeo_inv[0])*qEdge[6]+132.815661727072*jacobgeo_inv[1]*qEdge[3]+148.492424049175*jacobgeo_inv[2]*qEdge[2]); 
  fEdge[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.21320343559643*jacobgeo_inv[0])*qEdge[7]+21.21320343559643*jacobgeo_inv[1]*qEdge[5]); 

  double vol_incr[8] = {0.0}; 
  vol_incr[1] = 4.743416490252569*fSkin[1]*coeff[2]+2.121320343559642*fSkin[0]*coeff[1]; 
  vol_incr[3] = 4.743416490252569*coeff[2]*fSkin[3]+2.121320343559642*coeff[1]*fSkin[2]; 
  vol_incr[4] = 14.23024947075771*coeff[2]*fSkin[4]+10.60660171779821*fSkin[0]*coeff[2]+9.48683298050514*coeff[1]*fSkin[1]+4.743416490252569*coeff[0]*fSkin[0]; 
  vol_incr[6] = 14.23024947075771*coeff[2]*fSkin[6]+9.48683298050514*coeff[1]*fSkin[3]+10.60660171779821*coeff[2]*fSkin[2]+4.743416490252569*coeff[0]*fSkin[2]; 
  vol_incr[7] = 4.743416490252569*coeff[2]*fSkin[7]+2.121320343559642*coeff[1]*fSkin[5]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[4]-0.4743416490252568*coeff[0]*fSkin[4]-0.5303300858899105*coeff[2]*fEdge[4]+0.4743416490252568*coeff[0]*fEdge[4]+0.9413981457120035*fSkin[1]*coeff[2]+0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]-0.8420120990817169*coeff[0]*fSkin[1]-0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.2487763020014165*coeff[2]*fSkin[4]-0.5188111786213743*coeff[1]*fSkin[4]-1.121119609893386*coeff[0]*fSkin[4]-1.588341005085966*coeff[2]*fEdge[4]-0.5188111786213743*coeff[1]*fEdge[4]+0.5220480626221115*coeff[0]*fEdge[4]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]+0.5990715472712751*fSkin[0]*coeff[2]-1.96837794103419*fEdge[0]*coeff[2]-0.7463289060042488*coeff[1]*fSkin[1]-1.889300930982805*coeff[0]*fSkin[1]+0.7463289060042488*coeff[1]*fEdge[1]-1.027514541411701*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-1.454384534777511*coeff[0]*fSkin[0]+0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5303300858899104*coeff[2]*fSkin[6]-0.4743416490252568*coeff[0]*fSkin[6]-0.5303300858899104*coeff[2]*fEdge[6]+0.4743416490252568*coeff[0]*fEdge[6]+0.9413981457120035*coeff[2]*fSkin[3]-0.8420120990817169*coeff[0]*fSkin[3]+0.9413981457120035*coeff[2]*fEdge[3]-0.8420120990817169*coeff[0]*fEdge[3]+0.7411588266019635*coeff[2]*fSkin[2]-0.6629126073623879*coeff[0]*fSkin[2]-0.7411588266019635*coeff[2]*fEdge[2]+0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.2487763020014169*coeff[2]*fSkin[6]-0.5188111786213743*coeff[1]*fSkin[6]-1.121119609893386*coeff[0]*fSkin[6]-1.588341005085966*coeff[2]*fEdge[6]-0.5188111786213743*coeff[1]*fEdge[6]+0.5220480626221116*coeff[0]*fEdge[6]+0.6670429439417671*coeff[2]*fSkin[3]-0.7463289060042488*coeff[1]*fSkin[3]-1.889300930982805*coeff[0]*fSkin[3]+2.594055893106872*coeff[2]*fEdge[3]+0.7463289060042488*coeff[1]*fEdge[3]-1.027514541411701*coeff[0]*fEdge[3]+0.5990715472712751*coeff[2]*fSkin[2]-0.5303300858899105*coeff[1]*fSkin[2]-1.454384534777511*coeff[0]*fSkin[2]-1.96837794103419*coeff[2]*fEdge[2]-0.5303300858899105*coeff[1]*fEdge[2]+0.8420120990817168*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = (-1.40820177054373*coeff[2]*fSkin[4])-2.009347054626824*coeff[1]*fSkin[4]-2.220757234663999*coeff[0]*fSkin[4]-3.779910015670014*coeff[2]*fEdge[4]-2.009347054626824*coeff[1]*fEdge[4]-0.0994368911043575*coeff[0]*fEdge[4]-1.626614282316952*fSkin[1]*coeff[2]+5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.308931947855521*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]-3.551638458822559*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]-0.2139541240254559*coeff[0]*fEdge[1]-2.053959590644372*fSkin[0]*coeff[1]-2.053959590644372*fEdge[0]*coeff[1]-2.668171775767069*coeff[0]*fSkin[0]+0.2964635306407852*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 0.9413981457120036*coeff[2]*fSkin[7]-0.842012099081717*coeff[0]*fSkin[7]+0.9413981457120036*coeff[2]*fEdge[7]-0.842012099081717*coeff[0]*fEdge[7]+0.7411588266019635*coeff[2]*fSkin[5]-0.6629126073623879*coeff[0]*fSkin[5]-0.7411588266019635*coeff[2]*fEdge[5]+0.6629126073623879*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = (-1.40820177054373*coeff[2]*fSkin[6])-2.009347054626824*coeff[1]*fSkin[6]-2.220757234663999*coeff[0]*fSkin[6]-3.779910015670014*coeff[2]*fEdge[6]-2.009347054626824*coeff[1]*fEdge[6]-0.0994368911043575*coeff[0]*fEdge[6]-1.626614282316953*coeff[2]*fSkin[3]-2.890519423747656*coeff[1]*fSkin[3]-3.551638458822559*coeff[0]*fSkin[3]+5.836674777725538*coeff[2]*fEdge[3]+2.890519423747656*coeff[1]*fEdge[3]-0.2139541240254559*coeff[0]*fEdge[3]-0.9943689110435823*coeff[2]*fSkin[2]-2.053959590644372*coeff[1]*fSkin[2]-2.668171775767069*coeff[0]*fSkin[2]-4.308931947855521*coeff[2]*fEdge[2]-2.053959590644372*coeff[1]*fEdge[2]+0.2964635306407852*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 0.6670429439417671*coeff[2]*fSkin[7]-0.7463289060042488*coeff[1]*fSkin[7]-1.889300930982805*coeff[0]*fSkin[7]+2.594055893106872*coeff[2]*fEdge[7]+0.7463289060042488*coeff[1]*fEdge[7]-1.027514541411701*coeff[0]*fEdge[7]+0.599071547271275*coeff[2]*fSkin[5]-0.5303300858899104*coeff[1]*fSkin[5]-1.454384534777511*coeff[0]*fSkin[5]-1.96837794103419*coeff[2]*fEdge[5]-0.5303300858899104*coeff[1]*fEdge[5]+0.842012099081717*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 1.530931089239486*coeff[2]*fSkin[4]-1.185854122563142*coeff[1]*fSkin[4]+0.6846531968814573*coeff[0]*fSkin[4]-1.976423537605237*fSkin[1]*coeff[2]+1.369306393762915*fSkin[0]*coeff[2]+1.530931089239486*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = 1.530931089239486*coeff[2]*fSkin[6]-1.185854122563142*coeff[1]*fSkin[6]+0.6846531968814574*coeff[0]*fSkin[6]-1.976423537605237*coeff[2]*fSkin[3]+1.530931089239486*coeff[1]*fSkin[3]-0.883883476483184*coeff[0]*fSkin[3]+1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = (-5.929270612815711*coeff[2]*fSkin[4])+4.592793267718456*coeff[1]*fSkin[4]-2.651650429449552*coeff[0]*fSkin[4]+7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]+3.423265984407287*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = (-5.929270612815711*coeff[2]*fSkin[6])+4.592793267718456*coeff[1]*fSkin[6]-2.651650429449552*coeff[0]*fSkin[6]+7.65465544619743*coeff[2]*fSkin[3]-5.929270612815709*coeff[1]*fSkin[3]+3.423265984407287*coeff[0]*fSkin[3]-5.303300858899106*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = (-1.976423537605237*coeff[2]*fSkin[7])+1.530931089239486*coeff[1]*fSkin[7]-0.883883476483184*coeff[0]*fSkin[7]+1.369306393762915*coeff[2]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[5]+0.6123724356957944*coeff[0]*fSkin[5]; 

  } else { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[4]-0.4743416490252568*coeff[0]*fSkin[4]-0.5303300858899105*coeff[2]*fEdge[4]+0.4743416490252568*coeff[0]*fEdge[4]-0.9413981457120035*fSkin[1]*coeff[2]-0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]+0.8420120990817169*coeff[0]*fSkin[1]+0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.2487763020014165*coeff[2]*fSkin[4])-0.5188111786213743*coeff[1]*fSkin[4]+1.121119609893386*coeff[0]*fSkin[4]+1.588341005085966*coeff[2]*fEdge[4]-0.5188111786213743*coeff[1]*fEdge[4]-0.5220480626221115*coeff[0]*fEdge[4]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]-0.5990715472712751*fSkin[0]*coeff[2]+1.96837794103419*fEdge[0]*coeff[2]+0.7463289060042488*coeff[1]*fSkin[1]-1.889300930982805*coeff[0]*fSkin[1]-0.7463289060042488*coeff[1]*fEdge[1]-1.027514541411701*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+1.454384534777511*coeff[0]*fSkin[0]-0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5303300858899104*coeff[2]*fSkin[6]-0.4743416490252568*coeff[0]*fSkin[6]-0.5303300858899104*coeff[2]*fEdge[6]+0.4743416490252568*coeff[0]*fEdge[6]-0.9413981457120035*coeff[2]*fSkin[3]+0.8420120990817169*coeff[0]*fSkin[3]-0.9413981457120035*coeff[2]*fEdge[3]+0.8420120990817169*coeff[0]*fEdge[3]+0.7411588266019635*coeff[2]*fSkin[2]-0.6629126073623879*coeff[0]*fSkin[2]-0.7411588266019635*coeff[2]*fEdge[2]+0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.2487763020014169*coeff[2]*fSkin[6])-0.5188111786213743*coeff[1]*fSkin[6]+1.121119609893386*coeff[0]*fSkin[6]+1.588341005085966*coeff[2]*fEdge[6]-0.5188111786213743*coeff[1]*fEdge[6]-0.5220480626221116*coeff[0]*fEdge[6]+0.6670429439417671*coeff[2]*fSkin[3]+0.7463289060042488*coeff[1]*fSkin[3]-1.889300930982805*coeff[0]*fSkin[3]+2.594055893106872*coeff[2]*fEdge[3]-0.7463289060042488*coeff[1]*fEdge[3]-1.027514541411701*coeff[0]*fEdge[3]-0.5990715472712751*coeff[2]*fSkin[2]-0.5303300858899105*coeff[1]*fSkin[2]+1.454384534777511*coeff[0]*fSkin[2]+1.96837794103419*coeff[2]*fEdge[2]-0.5303300858899105*coeff[1]*fEdge[2]-0.8420120990817168*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = (-1.40820177054373*coeff[2]*fSkin[4])+2.009347054626824*coeff[1]*fSkin[4]-2.220757234663999*coeff[0]*fSkin[4]-3.779910015670014*coeff[2]*fEdge[4]+2.009347054626824*coeff[1]*fEdge[4]-0.0994368911043575*coeff[0]*fEdge[4]+1.626614282316952*fSkin[1]*coeff[2]-5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.308931947855521*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]+3.551638458822559*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]+0.2139541240254559*coeff[0]*fEdge[1]+2.053959590644372*fSkin[0]*coeff[1]+2.053959590644372*fEdge[0]*coeff[1]-2.668171775767069*coeff[0]*fSkin[0]+0.2964635306407852*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = (-0.9413981457120036*coeff[2]*fSkin[7])+0.842012099081717*coeff[0]*fSkin[7]-0.9413981457120036*coeff[2]*fEdge[7]+0.842012099081717*coeff[0]*fEdge[7]+0.7411588266019635*coeff[2]*fSkin[5]-0.6629126073623879*coeff[0]*fSkin[5]-0.7411588266019635*coeff[2]*fEdge[5]+0.6629126073623879*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = (-1.40820177054373*coeff[2]*fSkin[6])+2.009347054626824*coeff[1]*fSkin[6]-2.220757234663999*coeff[0]*fSkin[6]-3.779910015670014*coeff[2]*fEdge[6]+2.009347054626824*coeff[1]*fEdge[6]-0.0994368911043575*coeff[0]*fEdge[6]+1.626614282316953*coeff[2]*fSkin[3]-2.890519423747656*coeff[1]*fSkin[3]+3.551638458822559*coeff[0]*fSkin[3]-5.836674777725538*coeff[2]*fEdge[3]+2.890519423747656*coeff[1]*fEdge[3]+0.2139541240254559*coeff[0]*fEdge[3]-0.9943689110435823*coeff[2]*fSkin[2]+2.053959590644372*coeff[1]*fSkin[2]-2.668171775767069*coeff[0]*fSkin[2]-4.308931947855521*coeff[2]*fEdge[2]+2.053959590644372*coeff[1]*fEdge[2]+0.2964635306407852*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 0.6670429439417671*coeff[2]*fSkin[7]+0.7463289060042488*coeff[1]*fSkin[7]-1.889300930982805*coeff[0]*fSkin[7]+2.594055893106872*coeff[2]*fEdge[7]-0.7463289060042488*coeff[1]*fEdge[7]-1.027514541411701*coeff[0]*fEdge[7]-0.599071547271275*coeff[2]*fSkin[5]-0.5303300858899104*coeff[1]*fSkin[5]+1.454384534777511*coeff[0]*fSkin[5]+1.96837794103419*coeff[2]*fEdge[5]-0.5303300858899104*coeff[1]*fEdge[5]-0.842012099081717*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = (-1.530931089239486*coeff[2]*fSkin[4])-1.185854122563142*coeff[1]*fSkin[4]-0.6846531968814573*coeff[0]*fSkin[4]-1.976423537605237*fSkin[1]*coeff[2]-1.369306393762915*fSkin[0]*coeff[2]-1.530931089239486*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = (-1.530931089239486*coeff[2]*fSkin[6])-1.185854122563142*coeff[1]*fSkin[6]-0.6846531968814574*coeff[0]*fSkin[6]-1.976423537605237*coeff[2]*fSkin[3]-1.530931089239486*coeff[1]*fSkin[3]-0.883883476483184*coeff[0]*fSkin[3]-1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = (-5.929270612815711*coeff[2]*fSkin[4])-4.592793267718456*coeff[1]*fSkin[4]-2.651650429449552*coeff[0]*fSkin[4]-7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]-3.423265984407287*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = (-5.929270612815711*coeff[2]*fSkin[6])-4.592793267718456*coeff[1]*fSkin[6]-2.651650429449552*coeff[0]*fSkin[6]-7.65465544619743*coeff[2]*fSkin[3]-5.929270612815709*coeff[1]*fSkin[3]-3.423265984407287*coeff[0]*fSkin[3]-5.303300858899106*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = (-1.976423537605237*coeff[2]*fSkin[7])-1.530931089239486*coeff[1]*fSkin[7]-0.883883476483184*coeff[0]*fSkin[7]-1.369306393762915*coeff[2]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[5]-0.6123724356957944*coeff[0]*fSkin[5]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 

  }

  return 0.;
}

