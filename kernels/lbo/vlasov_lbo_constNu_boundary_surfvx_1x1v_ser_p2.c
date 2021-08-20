#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double Gdiff2_l[3]; 
  double Gdiff2_r[3]; 


  Gdiff2_l[0] = 0.015625*(15.65247584249852*nuVtSqSum_l[1]*fSkin[7]+15.65247584249852*nuVtSqSum_l[1]*fEdge[7]-22.5166604983954*nuVtSqSum_l[2]*fSkin[6]+22.5166604983954*nuVtSqSum_l[2]*fEdge[6]+15.65247584249852*nuVtSqSum_l[0]*fSkin[5]+15.65247584249852*nuVtSqSum_l[0]*fEdge[5]+16.0*nuVtSqSum_l[2]*fSkin[4]+16.0*nuVtSqSum_l[2]*fEdge[4]-22.5166604983954*nuVtSqSum_l[1]*fSkin[3]+22.5166604983954*nuVtSqSum_l[1]*fEdge[3]-22.5166604983954*nuVtSqSum_l[0]*fSkin[2]+22.5166604983954*nuVtSqSum_l[0]*fEdge[2]+(16.0*fSkin[1]+16.0*fEdge[1])*nuVtSqSum_l[1]+(16.0*fSkin[0]+16.0*fEdge[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = 0.003125*((70.0*nuVtSqSum_l[2]+78.26237921249266*nuVtSqSum_l[0])*fSkin[7]+(70.0*nuVtSqSum_l[2]+78.26237921249266*nuVtSqSum_l[0])*fEdge[7]-100.6975670013928*nuVtSqSum_l[1]*fSkin[6]+100.6975670013928*nuVtSqSum_l[1]*fEdge[6]+78.26237921249266*nuVtSqSum_l[1]*fSkin[5]+78.26237921249266*nuVtSqSum_l[1]*fEdge[5]+71.55417527999328*nuVtSqSum_l[1]*fSkin[4]+71.55417527999328*nuVtSqSum_l[1]*fEdge[4]+((-100.6975670013928*nuVtSqSum_l[2])-112.583302491977*nuVtSqSum_l[0])*fSkin[3]+(100.6975670013928*nuVtSqSum_l[2]+112.583302491977*nuVtSqSum_l[0])*fEdge[3]+(71.55417527999328*fSkin[1]+71.55417527999328*fEdge[1])*nuVtSqSum_l[2]-112.583302491977*nuVtSqSum_l[1]*fSkin[2]+112.583302491977*nuVtSqSum_l[1]*fEdge[2]+(80.0*fSkin[0]+80.0*fEdge[0])*nuVtSqSum_l[1]+80.0*nuVtSqSum_l[0]*fSkin[1]+80.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff2_l[2] = 4.464285714285714E-4*(490.0*nuVtSqSum_l[1]*fSkin[7]+490.0*nuVtSqSum_l[1]*fEdge[7]+((-503.4878350069642*nuVtSqSum_l[2])-788.0831174438391*nuVtSqSum_l[0])*fSkin[6]+(503.4878350069642*nuVtSqSum_l[2]+788.0831174438391*nuVtSqSum_l[0])*fEdge[6]+547.8366544874486*nuVtSqSum_l[2]*fSkin[5]+547.8366544874486*nuVtSqSum_l[2]*fEdge[5]+(357.7708763999664*nuVtSqSum_l[2]+560.0*nuVtSqSum_l[0])*fSkin[4]+(357.7708763999664*nuVtSqSum_l[2]+560.0*nuVtSqSum_l[0])*fEdge[4]-704.8829690097499*nuVtSqSum_l[1]*fSkin[3]+704.8829690097499*nuVtSqSum_l[1]*fEdge[3]+((-788.0831174438391*fSkin[2])+788.0831174438391*fEdge[2]+560.0*fSkin[0]+560.0*fEdge[0])*nuVtSqSum_l[2]+(500.8792269599529*fSkin[1]+500.8792269599529*fEdge[1])*nuVtSqSum_l[1]); 


  Gdiff2_r[0] = 0.015625*(15.65247584249852*nuVtSqSum_r[1]*fSkin[7]+15.65247584249852*nuVtSqSum_r[1]*fEdge[7]+22.5166604983954*nuVtSqSum_r[2]*fSkin[6]-22.5166604983954*nuVtSqSum_r[2]*fEdge[6]+15.65247584249852*nuVtSqSum_r[0]*fSkin[5]+15.65247584249852*nuVtSqSum_r[0]*fEdge[5]+16.0*nuVtSqSum_r[2]*fSkin[4]+16.0*nuVtSqSum_r[2]*fEdge[4]+22.5166604983954*nuVtSqSum_r[1]*fSkin[3]-22.5166604983954*nuVtSqSum_r[1]*fEdge[3]+22.5166604983954*nuVtSqSum_r[0]*fSkin[2]-22.5166604983954*nuVtSqSum_r[0]*fEdge[2]+(16.0*fSkin[1]+16.0*fEdge[1])*nuVtSqSum_r[1]+(16.0*fSkin[0]+16.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = 0.003125*((70.0*nuVtSqSum_r[2]+78.26237921249266*nuVtSqSum_r[0])*fSkin[7]+(70.0*nuVtSqSum_r[2]+78.26237921249266*nuVtSqSum_r[0])*fEdge[7]+100.6975670013928*nuVtSqSum_r[1]*fSkin[6]-100.6975670013928*nuVtSqSum_r[1]*fEdge[6]+78.26237921249266*nuVtSqSum_r[1]*fSkin[5]+78.26237921249266*nuVtSqSum_r[1]*fEdge[5]+71.55417527999328*nuVtSqSum_r[1]*fSkin[4]+71.55417527999328*nuVtSqSum_r[1]*fEdge[4]+(100.6975670013928*nuVtSqSum_r[2]+112.583302491977*nuVtSqSum_r[0])*fSkin[3]+((-100.6975670013928*nuVtSqSum_r[2])-112.583302491977*nuVtSqSum_r[0])*fEdge[3]+(71.55417527999328*fSkin[1]+71.55417527999328*fEdge[1])*nuVtSqSum_r[2]+112.583302491977*nuVtSqSum_r[1]*fSkin[2]-112.583302491977*nuVtSqSum_r[1]*fEdge[2]+(80.0*fSkin[0]+80.0*fEdge[0])*nuVtSqSum_r[1]+80.0*nuVtSqSum_r[0]*fSkin[1]+80.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff2_r[2] = 4.464285714285714E-4*(490.0*nuVtSqSum_r[1]*fSkin[7]+490.0*nuVtSqSum_r[1]*fEdge[7]+(503.4878350069642*nuVtSqSum_r[2]+788.0831174438391*nuVtSqSum_r[0])*fSkin[6]+((-503.4878350069642*nuVtSqSum_r[2])-788.0831174438391*nuVtSqSum_r[0])*fEdge[6]+547.8366544874486*nuVtSqSum_r[2]*fSkin[5]+547.8366544874486*nuVtSqSum_r[2]*fEdge[5]+(357.7708763999664*nuVtSqSum_r[2]+560.0*nuVtSqSum_r[0])*fSkin[4]+(357.7708763999664*nuVtSqSum_r[2]+560.0*nuVtSqSum_r[0])*fEdge[4]+704.8829690097499*nuVtSqSum_r[1]*fSkin[3]-704.8829690097499*nuVtSqSum_r[1]*fEdge[3]+(788.0831174438391*fSkin[2]-788.0831174438391*fEdge[2]+560.0*fSkin[0]+560.0*fEdge[0])*nuVtSqSum_r[2]+(500.8792269599529*fSkin[1]+500.8792269599529*fEdge[1])*nuVtSqSum_r[1]); 

  if (edge == -1) { 

  const double *sumNuUx_r = &nuUSum[0]; 

  double alphaDrSurf_r[3]; 
  alphaDrSurf_r[0] = -0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUx_r[0]); 
  alphaDrSurf_r[1] = sumNuUx_r[1]; 
  alphaDrSurf_r[2] = sumNuUx_r[2]; 

  double fUpwindQuad_r[3];
  if (0.6324555320336759*alphaDrSurf_r[2]-0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = (-1.5*fSkin[7])+0.7745966692414833*fSkin[6]+1.118033988749895*fSkin[5]+0.4472135954999579*fSkin[4]-1.161895003862224*fSkin[3]+0.8660254037844386*fSkin[2]-0.6708203932499369*fSkin[1]+0.5*fSkin[0]; 
  } else { 
    fUpwindQuad_r[0] = (-1.5*fEdge[7])-0.7745966692414833*fEdge[6]+1.118033988749895*fEdge[5]+0.4472135954999579*fEdge[4]+1.161895003862224*fEdge[3]-0.8660254037844386*fEdge[2]-0.6708203932499369*fEdge[1]+0.5*fEdge[0]; 
  } 
  if (0.7071067811865475*alphaDrSurf_r[0]-0.7905694150420947*alphaDrSurf_r[2] > 0) { 
    fUpwindQuad_r[1] = (-0.9682458365518543*fSkin[6])+1.118033988749895*fSkin[5]-0.5590169943749475*fSkin[4]+0.8660254037844386*fSkin[2]+0.5*fSkin[0]; 
  } else { 
    fUpwindQuad_r[1] = 0.9682458365518543*fEdge[6]+1.118033988749895*fEdge[5]-0.5590169943749475*fEdge[4]-0.8660254037844386*fEdge[2]+0.5*fEdge[0]; 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]+0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = 1.5*fSkin[7]+0.7745966692414833*fSkin[6]+1.118033988749895*fSkin[5]+0.4472135954999579*fSkin[4]+1.161895003862224*fSkin[3]+0.8660254037844386*fSkin[2]+0.6708203932499369*fSkin[1]+0.5*fSkin[0]; 
  } else { 
    fUpwindQuad_r[2] = 1.5*fEdge[7]-0.7745966692414833*fEdge[6]+1.118033988749895*fEdge[5]+0.4472135954999579*fEdge[4]-1.161895003862224*fEdge[3]-0.8660254037844386*fEdge[2]+0.6708203932499369*fEdge[1]+0.5*fEdge[0]; 
  } 

  double fUpwind_r[3];
  fUpwind_r[0] = 0.07856742013183861*(5.0*fUpwindQuad_r[2]+8.0*fUpwindQuad_r[1]+5.0*fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5270462766947298*(fUpwindQuad_r[2]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3513641844631532*(fUpwindQuad_r[2]-2.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 

  double Gdiff_r[3]; 
  double Ghat_r[3]; 


  Gdiff_r[0] = -0.0125*(53.66563145999495*nuVtSqSum_r[1]*fSkin[7]-53.66563145999495*nuVtSqSum_r[1]*fEdge[7]+95.26279441628826*nuVtSqSum_r[2]*fSkin[6]+95.26279441628826*nuVtSqSum_r[2]*fEdge[6]+53.66563145999495*nuVtSqSum_r[0]*fSkin[5]-53.66563145999495*nuVtSqSum_r[0]*fEdge[5]+75.0*nuVtSqSum_r[2]*fSkin[4]-75.0*nuVtSqSum_r[2]*fEdge[4]+95.26279441628825*nuVtSqSum_r[1]*fSkin[3]+95.26279441628825*nuVtSqSum_r[1]*fEdge[3]+95.26279441628825*nuVtSqSum_r[0]*fSkin[2]+95.26279441628825*nuVtSqSum_r[0]*fEdge[2]+(75.0*fSkin[1]-75.0*fEdge[1])*nuVtSqSum_r[1]+(75.0*fSkin[0]-75.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = -0.0025*((240.0*nuVtSqSum_r[2]+268.3281572999747*nuVtSqSum_r[0])*fSkin[7]+((-240.0*nuVtSqSum_r[2])-268.3281572999747*nuVtSqSum_r[0])*fEdge[7]+426.0281680828158*nuVtSqSum_r[1]*fSkin[6]+426.0281680828158*nuVtSqSum_r[1]*fEdge[6]+268.3281572999748*nuVtSqSum_r[1]*fSkin[5]-268.3281572999748*nuVtSqSum_r[1]*fEdge[5]+335.4101966249686*nuVtSqSum_r[1]*fSkin[4]-335.4101966249686*nuVtSqSum_r[1]*fEdge[4]+(426.0281680828159*nuVtSqSum_r[2]+476.3139720814412*nuVtSqSum_r[0])*fSkin[3]+(426.0281680828159*nuVtSqSum_r[2]+476.3139720814412*nuVtSqSum_r[0])*fEdge[3]+(335.4101966249686*fSkin[1]-335.4101966249686*fEdge[1])*nuVtSqSum_r[2]+476.3139720814412*nuVtSqSum_r[1]*fSkin[2]+476.3139720814412*nuVtSqSum_r[1]*fEdge[2]+(375.0*fSkin[0]-375.0*fEdge[0])*nuVtSqSum_r[1]+375.0*nuVtSqSum_r[0]*fSkin[1]-375.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff_r[2] = -3.571428571428571E-4*(1680.0*nuVtSqSum_r[1]*fSkin[7]-1680.0*nuVtSqSum_r[1]*fEdge[7]+(2130.140840414079*nuVtSqSum_r[2]+3334.19780457009*nuVtSqSum_r[0])*fSkin[6]+(2130.140840414079*nuVtSqSum_r[2]+3334.19780457009*nuVtSqSum_r[0])*fEdge[6]+1878.297101099823*nuVtSqSum_r[2]*fSkin[5]-1878.297101099823*nuVtSqSum_r[2]*fEdge[5]+(1677.050983124843*nuVtSqSum_r[2]+2625.0*nuVtSqSum_r[0])*fSkin[4]+((-1677.050983124843*nuVtSqSum_r[2])-2625.0*nuVtSqSum_r[0])*fEdge[4]+2982.197176579711*nuVtSqSum_r[1]*fSkin[3]+2982.197176579711*nuVtSqSum_r[1]*fEdge[3]+(3334.197804570088*fSkin[2]+3334.197804570088*fEdge[2]+2625.0*fSkin[0]-2625.0*fEdge[0])*nuVtSqSum_r[2]+(2347.87137637478*fSkin[1]-2347.87137637478*fEdge[1])*nuVtSqSum_r[1]); 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.7071067811865475*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaDrSurf_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.4517539514526256*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[2]+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[1]; 

  out[0] += -0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += (-1.369306393762915*nuVtSqSum_l[1]*fSkin[7]*rdvSq4)+1.060660171779821*nuVtSqSum_l[2]*fSkin[6]*rdvSq4-1.369306393762915*nuVtSqSum_l[0]*fSkin[5]*rdvSq4-0.6123724356957944*nuVtSqSum_l[2]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[1]*nuVtSqSum_l[1]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[0]*rdvSq4+1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2; 
  out[3] += (-1.224744871391589*nuVtSqSum_l[2]*fSkin[7]*rdvSq4)-1.369306393762915*nuVtSqSum_l[0]*fSkin[7]*rdvSq4+0.9486832980505138*nuVtSqSum_l[1]*fSkin[6]*rdvSq4-1.369306393762915*nuVtSqSum_l[1]*fSkin[5]*rdvSq4-0.5477225575051661*nuVtSqSum_l[1]*fSkin[4]*rdvSq4+0.9486832980505137*nuVtSqSum_l[2]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[3]*rdvSq4-0.5477225575051661*fSkin[1]*nuVtSqSum_l[2]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[1]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2; 
  out[4] += -0.7071067811865475*Ghat_r[2]*rdv2; 
  out[5] += 5.303300858899106*nuVtSqSum_l[1]*fSkin[7]*rdvSq4-4.107919181288745*nuVtSqSum_l[2]*fSkin[6]*rdvSq4+5.303300858899105*nuVtSqSum_l[0]*fSkin[5]*rdvSq4+2.371708245126284*nuVtSqSum_l[2]*fSkin[4]*rdvSq4-4.107919181288745*nuVtSqSum_l[1]*fSkin[3]*rdvSq4-4.107919181288745*nuVtSqSum_l[0]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[1]*nuVtSqSum_l[1]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum_l[0]*rdvSq4+4.743416490252569*Gdiff2_r[0]*rdvSq4-1.581138830084189*Ghat_r[0]*rdv2; 
  out[6] += (-1.224744871391589*nuVtSqSum_l[1]*fSkin[7]*rdvSq4)+0.6776309271789384*nuVtSqSum_l[2]*fSkin[6]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[6]*rdvSq4-1.369306393762915*nuVtSqSum_l[2]*fSkin[5]*rdvSq4-0.3912303982179757*nuVtSqSum_l[2]*fSkin[4]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[4]*rdvSq4+0.9486832980505138*nuVtSqSum_l[1]*fSkin[3]*rdvSq4+1.060660171779821*fSkin[2]*nuVtSqSum_l[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[2]*rdvSq4+1.224744871391589*Gdiff2_r[2]*rdvSq4-0.5477225575051661*fSkin[1]*nuVtSqSum_l[1]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2; 
  out[7] += 4.743416490252569*nuVtSqSum_l[2]*fSkin[7]*rdvSq4+5.303300858899105*nuVtSqSum_l[0]*fSkin[7]*rdvSq4-3.674234614174766*nuVtSqSum_l[1]*fSkin[6]*rdvSq4+5.303300858899106*nuVtSqSum_l[1]*fSkin[5]*rdvSq4+2.121320343559642*nuVtSqSum_l[1]*fSkin[4]*rdvSq4-3.674234614174767*nuVtSqSum_l[2]*fSkin[3]*rdvSq4-4.107919181288746*nuVtSqSum_l[0]*fSkin[3]*rdvSq4+2.121320343559642*fSkin[1]*nuVtSqSum_l[2]*rdvSq4-4.107919181288746*nuVtSqSum_l[1]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum_l[1]*rdvSq4+2.371708245126284*nuVtSqSum_l[0]*fSkin[1]*rdvSq4+4.743416490252569*Gdiff2_r[1]*rdvSq4-1.581138830084189*Ghat_r[1]*rdv2; 

  } else {

  const double *sumNuUx_l = &nuUSum[0]; 

  double alphaDrSurf_l[3]; 
  alphaDrSurf_l[0] = 0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUx_l[0]); 
  alphaDrSurf_l[1] = -1.0*sumNuUx_l[1]; 
  alphaDrSurf_l[2] = -1.0*sumNuUx_l[2]; 

  double fUpwindQuad_l[3];
  if (0.6324555320336759*alphaDrSurf_l[2]-0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = (-1.5*fEdge[7])+0.7745966692414833*fEdge[6]+1.118033988749895*fEdge[5]+0.4472135954999579*fEdge[4]-1.161895003862224*fEdge[3]+0.8660254037844386*fEdge[2]-0.6708203932499369*fEdge[1]+0.5*fEdge[0]; 
  } else { 
    fUpwindQuad_l[0] = (-1.5*fSkin[7])-0.7745966692414833*fSkin[6]+1.118033988749895*fSkin[5]+0.4472135954999579*fSkin[4]+1.161895003862224*fSkin[3]-0.8660254037844386*fSkin[2]-0.6708203932499369*fSkin[1]+0.5*fSkin[0]; 
  } 
  if (0.7071067811865475*alphaDrSurf_l[0]-0.7905694150420947*alphaDrSurf_l[2] > 0) { 
    fUpwindQuad_l[1] = (-0.9682458365518543*fEdge[6])+1.118033988749895*fEdge[5]-0.5590169943749475*fEdge[4]+0.8660254037844386*fEdge[2]+0.5*fEdge[0]; 
  } else { 
    fUpwindQuad_l[1] = 0.9682458365518543*fSkin[6]+1.118033988749895*fSkin[5]-0.5590169943749475*fSkin[4]-0.8660254037844386*fSkin[2]+0.5*fSkin[0]; 
  } 
  if (0.6324555320336759*alphaDrSurf_l[2]+0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = 1.5*fEdge[7]+0.7745966692414833*fEdge[6]+1.118033988749895*fEdge[5]+0.4472135954999579*fEdge[4]+1.161895003862224*fEdge[3]+0.8660254037844386*fEdge[2]+0.6708203932499369*fEdge[1]+0.5*fEdge[0]; 
  } else { 
    fUpwindQuad_l[2] = 1.5*fSkin[7]-0.7745966692414833*fSkin[6]+1.118033988749895*fSkin[5]+0.4472135954999579*fSkin[4]-1.161895003862224*fSkin[3]-0.8660254037844386*fSkin[2]+0.6708203932499369*fSkin[1]+0.5*fSkin[0]; 
  } 

  double fUpwind_l[3];
  fUpwind_l[0] = 0.07856742013183861*(5.0*fUpwindQuad_l[2]+8.0*fUpwindQuad_l[1]+5.0*fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5270462766947298*(fUpwindQuad_l[2]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3513641844631532*(fUpwindQuad_l[2]-2.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 

  double Gdiff_l[3]; 
  double Ghat_l[3]; 


  Gdiff_l[0] = 0.0125*(53.66563145999495*nuVtSqSum_l[1]*fSkin[7]-53.66563145999495*nuVtSqSum_l[1]*fEdge[7]-95.26279441628826*nuVtSqSum_l[2]*fSkin[6]-95.26279441628826*nuVtSqSum_l[2]*fEdge[6]+53.66563145999495*nuVtSqSum_l[0]*fSkin[5]-53.66563145999495*nuVtSqSum_l[0]*fEdge[5]+75.0*nuVtSqSum_l[2]*fSkin[4]-75.0*nuVtSqSum_l[2]*fEdge[4]-95.26279441628825*nuVtSqSum_l[1]*fSkin[3]-95.26279441628825*nuVtSqSum_l[1]*fEdge[3]-95.26279441628825*nuVtSqSum_l[0]*fSkin[2]-95.26279441628825*nuVtSqSum_l[0]*fEdge[2]+(75.0*fSkin[1]-75.0*fEdge[1])*nuVtSqSum_l[1]+(75.0*fSkin[0]-75.0*fEdge[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = 0.0025*((240.0*nuVtSqSum_l[2]+268.3281572999747*nuVtSqSum_l[0])*fSkin[7]+((-240.0*nuVtSqSum_l[2])-268.3281572999747*nuVtSqSum_l[0])*fEdge[7]-426.0281680828158*nuVtSqSum_l[1]*fSkin[6]-426.0281680828158*nuVtSqSum_l[1]*fEdge[6]+268.3281572999748*nuVtSqSum_l[1]*fSkin[5]-268.3281572999748*nuVtSqSum_l[1]*fEdge[5]+335.4101966249686*nuVtSqSum_l[1]*fSkin[4]-335.4101966249686*nuVtSqSum_l[1]*fEdge[4]+((-426.0281680828159*nuVtSqSum_l[2])-476.3139720814412*nuVtSqSum_l[0])*fSkin[3]+((-426.0281680828159*nuVtSqSum_l[2])-476.3139720814412*nuVtSqSum_l[0])*fEdge[3]+(335.4101966249686*fSkin[1]-335.4101966249686*fEdge[1])*nuVtSqSum_l[2]-476.3139720814412*nuVtSqSum_l[1]*fSkin[2]-476.3139720814412*nuVtSqSum_l[1]*fEdge[2]+(375.0*fSkin[0]-375.0*fEdge[0])*nuVtSqSum_l[1]+375.0*nuVtSqSum_l[0]*fSkin[1]-375.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff_l[2] = 3.571428571428571E-4*(1680.0*nuVtSqSum_l[1]*fSkin[7]-1680.0*nuVtSqSum_l[1]*fEdge[7]+((-2130.140840414079*nuVtSqSum_l[2])-3334.19780457009*nuVtSqSum_l[0])*fSkin[6]+((-2130.140840414079*nuVtSqSum_l[2])-3334.19780457009*nuVtSqSum_l[0])*fEdge[6]+1878.297101099823*nuVtSqSum_l[2]*fSkin[5]-1878.297101099823*nuVtSqSum_l[2]*fEdge[5]+(1677.050983124843*nuVtSqSum_l[2]+2625.0*nuVtSqSum_l[0])*fSkin[4]+((-1677.050983124843*nuVtSqSum_l[2])-2625.0*nuVtSqSum_l[0])*fEdge[4]-2982.197176579711*nuVtSqSum_l[1]*fSkin[3]-2982.197176579711*nuVtSqSum_l[1]*fEdge[3]+((-3334.197804570088*fSkin[2])-3334.197804570088*fEdge[2]+2625.0*fSkin[0]-2625.0*fEdge[0])*nuVtSqSum_l[2]+(2347.87137637478*fSkin[1]-2347.87137637478*fEdge[1])*nuVtSqSum_l[1]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.7071067811865475*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaDrSurf_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.4517539514526256*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[2]+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[1]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.369306393762915*nuVtSqSum_r[1]*fSkin[7]*rdvSq4+1.060660171779821*nuVtSqSum_r[2]*fSkin[6]*rdvSq4+1.369306393762915*nuVtSqSum_r[0]*fSkin[5]*rdvSq4+0.6123724356957944*nuVtSqSum_r[2]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[1]*nuVtSqSum_r[1]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 1.224744871391589*nuVtSqSum_r[2]*fSkin[7]*rdvSq4+1.369306393762915*nuVtSqSum_r[0]*fSkin[7]*rdvSq4+0.9486832980505138*nuVtSqSum_r[1]*fSkin[6]*rdvSq4+1.369306393762915*nuVtSqSum_r[1]*fSkin[5]*rdvSq4+0.5477225575051661*nuVtSqSum_r[1]*fSkin[4]*rdvSq4+0.9486832980505137*nuVtSqSum_r[2]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[3]*rdvSq4+0.5477225575051661*fSkin[1]*nuVtSqSum_r[2]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[1]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_l[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[2]*rdv2; 
  out[5] += 5.303300858899106*nuVtSqSum_r[1]*fSkin[7]*rdvSq4+4.107919181288745*nuVtSqSum_r[2]*fSkin[6]*rdvSq4+5.303300858899105*nuVtSqSum_r[0]*fSkin[5]*rdvSq4+2.371708245126284*nuVtSqSum_r[2]*fSkin[4]*rdvSq4+4.107919181288745*nuVtSqSum_r[1]*fSkin[3]*rdvSq4+4.107919181288745*nuVtSqSum_r[0]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[1]*nuVtSqSum_r[1]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum_r[0]*rdvSq4+4.743416490252569*Gdiff2_l[0]*rdvSq4+1.581138830084189*Ghat_l[0]*rdv2; 
  out[6] += 1.224744871391589*nuVtSqSum_r[1]*fSkin[7]*rdvSq4+0.6776309271789384*nuVtSqSum_r[2]*fSkin[6]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[6]*rdvSq4+1.369306393762915*nuVtSqSum_r[2]*fSkin[5]*rdvSq4+0.3912303982179757*nuVtSqSum_r[2]*fSkin[4]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[4]*rdvSq4+0.9486832980505138*nuVtSqSum_r[1]*fSkin[3]*rdvSq4+1.060660171779821*fSkin[2]*nuVtSqSum_r[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4+0.5477225575051661*fSkin[1]*nuVtSqSum_r[1]*rdvSq4-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 4.743416490252569*nuVtSqSum_r[2]*fSkin[7]*rdvSq4+5.303300858899105*nuVtSqSum_r[0]*fSkin[7]*rdvSq4+3.674234614174766*nuVtSqSum_r[1]*fSkin[6]*rdvSq4+5.303300858899106*nuVtSqSum_r[1]*fSkin[5]*rdvSq4+2.121320343559642*nuVtSqSum_r[1]*fSkin[4]*rdvSq4+3.674234614174767*nuVtSqSum_r[2]*fSkin[3]*rdvSq4+4.107919181288746*nuVtSqSum_r[0]*fSkin[3]*rdvSq4+2.121320343559642*fSkin[1]*nuVtSqSum_r[2]*rdvSq4+4.107919181288746*nuVtSqSum_r[1]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum_r[1]*rdvSq4+2.371708245126284*nuVtSqSum_r[0]*fSkin[1]*rdvSq4+4.743416490252569*Gdiff2_l[1]*rdvSq4+1.581138830084189*Ghat_l[1]*rdv2; 
  } 
} 
