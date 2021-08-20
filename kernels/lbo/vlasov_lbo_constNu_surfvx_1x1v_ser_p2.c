#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx_l = &nuUSum_l[0]; 
  const double *sumNuUx_r = &nuUSum_r[0]; 

  double alphaDrSurf_l[3]; 
  alphaDrSurf_l[0] = 0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUx_l[0]); 
  alphaDrSurf_l[1] = -1.0*sumNuUx_l[1]; 
  alphaDrSurf_l[2] = -1.0*sumNuUx_l[2]; 

  double alphaDrSurf_r[3]; 
  alphaDrSurf_r[0] = -0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUx_r[0]); 
  alphaDrSurf_r[1] = sumNuUx_r[1]; 
  alphaDrSurf_r[2] = sumNuUx_r[2]; 

  double fUpwindQuad_l[3];
  if (0.6324555320336759*alphaDrSurf_l[2]-0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = (-1.5*fl[7])+0.7745966692414833*fl[6]+1.118033988749895*fl[5]+0.4472135954999579*fl[4]-1.161895003862224*fl[3]+0.8660254037844386*fl[2]-0.6708203932499369*fl[1]+0.5*fl[0]; 
  } else { 
    fUpwindQuad_l[0] = (-1.5*fc[7])-0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]+1.161895003862224*fc[3]-0.8660254037844386*fc[2]-0.6708203932499369*fc[1]+0.5*fc[0]; 
  } 
  if (0.7071067811865475*alphaDrSurf_l[0]-0.7905694150420947*alphaDrSurf_l[2] > 0) { 
    fUpwindQuad_l[1] = (-0.9682458365518543*fl[6])+1.118033988749895*fl[5]-0.5590169943749475*fl[4]+0.8660254037844386*fl[2]+0.5*fl[0]; 
  } else { 
    fUpwindQuad_l[1] = 0.9682458365518543*fc[6]+1.118033988749895*fc[5]-0.5590169943749475*fc[4]-0.8660254037844386*fc[2]+0.5*fc[0]; 
  } 
  if (0.6324555320336759*alphaDrSurf_l[2]+0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = 1.5*fl[7]+0.7745966692414833*fl[6]+1.118033988749895*fl[5]+0.4472135954999579*fl[4]+1.161895003862224*fl[3]+0.8660254037844386*fl[2]+0.6708203932499369*fl[1]+0.5*fl[0]; 
  } else { 
    fUpwindQuad_l[2] = 1.5*fc[7]-0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]-1.161895003862224*fc[3]-0.8660254037844386*fc[2]+0.6708203932499369*fc[1]+0.5*fc[0]; 
  } 

  double fUpwind_l[3];
  fUpwind_l[0] = 0.07856742013183861*(5.0*fUpwindQuad_l[2]+8.0*fUpwindQuad_l[1]+5.0*fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5270462766947298*(fUpwindQuad_l[2]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3513641844631532*(fUpwindQuad_l[2]-2.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 

  double fUpwindQuad_r[3];
  if (0.6324555320336759*alphaDrSurf_r[2]-0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = (-1.5*fc[7])+0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]-1.161895003862224*fc[3]+0.8660254037844386*fc[2]-0.6708203932499369*fc[1]+0.5*fc[0]; 
  } else { 
    fUpwindQuad_r[0] = (-1.5*fr[7])-0.7745966692414833*fr[6]+1.118033988749895*fr[5]+0.4472135954999579*fr[4]+1.161895003862224*fr[3]-0.8660254037844386*fr[2]-0.6708203932499369*fr[1]+0.5*fr[0]; 
  } 
  if (0.7071067811865475*alphaDrSurf_r[0]-0.7905694150420947*alphaDrSurf_r[2] > 0) { 
    fUpwindQuad_r[1] = (-0.9682458365518543*fc[6])+1.118033988749895*fc[5]-0.5590169943749475*fc[4]+0.8660254037844386*fc[2]+0.5*fc[0]; 
  } else { 
    fUpwindQuad_r[1] = 0.9682458365518543*fr[6]+1.118033988749895*fr[5]-0.5590169943749475*fr[4]-0.8660254037844386*fr[2]+0.5*fr[0]; 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]+0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = 1.5*fc[7]+0.7745966692414833*fc[6]+1.118033988749895*fc[5]+0.4472135954999579*fc[4]+1.161895003862224*fc[3]+0.8660254037844386*fc[2]+0.6708203932499369*fc[1]+0.5*fc[0]; 
  } else { 
    fUpwindQuad_r[2] = 1.5*fr[7]-0.7745966692414833*fr[6]+1.118033988749895*fr[5]+0.4472135954999579*fr[4]-1.161895003862224*fr[3]-0.8660254037844386*fr[2]+0.6708203932499369*fr[1]+0.5*fr[0]; 
  } 

  double fUpwind_r[3];
  fUpwind_r[0] = 0.07856742013183861*(5.0*fUpwindQuad_r[2]+8.0*fUpwindQuad_r[1]+5.0*fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5270462766947298*(fUpwindQuad_r[2]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3513641844631532*(fUpwindQuad_r[2]-2.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 

  double Gdiff_l[3]; 
  double Gdiff_r[3]; 
  double Ghat_l[3]; 
  double Ghat_r[3]; 
  double Gdiff2_l[3]; 
  double Gdiff2_r[3]; 


  Gdiff2_l[0] = 0.015625*(15.65247584249852*nuVtSqSum_l[1]*fl[7]+15.65247584249852*nuVtSqSum_l[1]*fc[7]+22.5166604983954*nuVtSqSum_l[2]*fl[6]-22.5166604983954*nuVtSqSum_l[2]*fc[6]+15.65247584249852*nuVtSqSum_l[0]*fl[5]+15.65247584249852*nuVtSqSum_l[0]*fc[5]+16.0*nuVtSqSum_l[2]*fl[4]+16.0*nuVtSqSum_l[2]*fc[4]+22.5166604983954*nuVtSqSum_l[1]*fl[3]-22.5166604983954*nuVtSqSum_l[1]*fc[3]+22.5166604983954*nuVtSqSum_l[0]*fl[2]-22.5166604983954*nuVtSqSum_l[0]*fc[2]+(16.0*fl[1]+16.0*fc[1])*nuVtSqSum_l[1]+(16.0*fl[0]+16.0*fc[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = 0.003125*((70.0*nuVtSqSum_l[2]+78.26237921249266*nuVtSqSum_l[0])*fl[7]+(70.0*nuVtSqSum_l[2]+78.26237921249266*nuVtSqSum_l[0])*fc[7]+100.6975670013928*nuVtSqSum_l[1]*fl[6]-100.6975670013928*nuVtSqSum_l[1]*fc[6]+78.26237921249266*nuVtSqSum_l[1]*fl[5]+78.26237921249266*nuVtSqSum_l[1]*fc[5]+71.55417527999328*nuVtSqSum_l[1]*fl[4]+71.55417527999328*nuVtSqSum_l[1]*fc[4]+(100.6975670013928*nuVtSqSum_l[2]+112.583302491977*nuVtSqSum_l[0])*fl[3]+((-100.6975670013928*nuVtSqSum_l[2])-112.583302491977*nuVtSqSum_l[0])*fc[3]+(71.55417527999328*fl[1]+71.55417527999328*fc[1])*nuVtSqSum_l[2]+112.583302491977*nuVtSqSum_l[1]*fl[2]-112.583302491977*nuVtSqSum_l[1]*fc[2]+(80.0*fl[0]+80.0*fc[0])*nuVtSqSum_l[1]+80.0*nuVtSqSum_l[0]*fl[1]+80.0*nuVtSqSum_l[0]*fc[1]); 
  Gdiff2_l[2] = 4.464285714285714E-4*(490.0*nuVtSqSum_l[1]*fl[7]+490.0*nuVtSqSum_l[1]*fc[7]+(503.4878350069642*nuVtSqSum_l[2]+788.0831174438391*nuVtSqSum_l[0])*fl[6]+((-503.4878350069642*nuVtSqSum_l[2])-788.0831174438391*nuVtSqSum_l[0])*fc[6]+547.8366544874486*nuVtSqSum_l[2]*fl[5]+547.8366544874486*nuVtSqSum_l[2]*fc[5]+(357.7708763999664*nuVtSqSum_l[2]+560.0*nuVtSqSum_l[0])*fl[4]+(357.7708763999664*nuVtSqSum_l[2]+560.0*nuVtSqSum_l[0])*fc[4]+704.8829690097499*nuVtSqSum_l[1]*fl[3]-704.8829690097499*nuVtSqSum_l[1]*fc[3]+(788.0831174438391*fl[2]-788.0831174438391*fc[2]+560.0*fl[0]+560.0*fc[0])*nuVtSqSum_l[2]+(500.8792269599529*fl[1]+500.8792269599529*fc[1])*nuVtSqSum_l[1]); 


  Gdiff2_r[0] = 0.015625*(15.65247584249852*nuVtSqSum_r[1]*fr[7]+15.65247584249852*nuVtSqSum_r[1]*fc[7]-22.5166604983954*nuVtSqSum_r[2]*fr[6]+22.5166604983954*nuVtSqSum_r[2]*fc[6]+15.65247584249852*nuVtSqSum_r[0]*fr[5]+15.65247584249852*nuVtSqSum_r[0]*fc[5]+16.0*nuVtSqSum_r[2]*fr[4]+16.0*nuVtSqSum_r[2]*fc[4]-22.5166604983954*nuVtSqSum_r[1]*fr[3]+22.5166604983954*nuVtSqSum_r[1]*fc[3]-22.5166604983954*nuVtSqSum_r[0]*fr[2]+22.5166604983954*nuVtSqSum_r[0]*fc[2]+(16.0*fr[1]+16.0*fc[1])*nuVtSqSum_r[1]+(16.0*fr[0]+16.0*fc[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = 0.003125*((70.0*nuVtSqSum_r[2]+78.26237921249266*nuVtSqSum_r[0])*fr[7]+(70.0*nuVtSqSum_r[2]+78.26237921249266*nuVtSqSum_r[0])*fc[7]-100.6975670013928*nuVtSqSum_r[1]*fr[6]+100.6975670013928*nuVtSqSum_r[1]*fc[6]+78.26237921249266*nuVtSqSum_r[1]*fr[5]+78.26237921249266*nuVtSqSum_r[1]*fc[5]+71.55417527999328*nuVtSqSum_r[1]*fr[4]+71.55417527999328*nuVtSqSum_r[1]*fc[4]+((-100.6975670013928*nuVtSqSum_r[2])-112.583302491977*nuVtSqSum_r[0])*fr[3]+(100.6975670013928*nuVtSqSum_r[2]+112.583302491977*nuVtSqSum_r[0])*fc[3]+(71.55417527999328*fr[1]+71.55417527999328*fc[1])*nuVtSqSum_r[2]-112.583302491977*nuVtSqSum_r[1]*fr[2]+112.583302491977*nuVtSqSum_r[1]*fc[2]+(80.0*fr[0]+80.0*fc[0])*nuVtSqSum_r[1]+80.0*nuVtSqSum_r[0]*fr[1]+80.0*nuVtSqSum_r[0]*fc[1]); 
  Gdiff2_r[2] = 4.464285714285714E-4*(490.0*nuVtSqSum_r[1]*fr[7]+490.0*nuVtSqSum_r[1]*fc[7]+((-503.4878350069642*nuVtSqSum_r[2])-788.0831174438391*nuVtSqSum_r[0])*fr[6]+(503.4878350069642*nuVtSqSum_r[2]+788.0831174438391*nuVtSqSum_r[0])*fc[6]+547.8366544874486*nuVtSqSum_r[2]*fr[5]+547.8366544874486*nuVtSqSum_r[2]*fc[5]+(357.7708763999664*nuVtSqSum_r[2]+560.0*nuVtSqSum_r[0])*fr[4]+(357.7708763999664*nuVtSqSum_r[2]+560.0*nuVtSqSum_r[0])*fc[4]-704.8829690097499*nuVtSqSum_r[1]*fr[3]+704.8829690097499*nuVtSqSum_r[1]*fc[3]+((-788.0831174438391*fr[2])+788.0831174438391*fc[2]+560.0*fr[0]+560.0*fc[0])*nuVtSqSum_r[2]+(500.8792269599529*fr[1]+500.8792269599529*fc[1])*nuVtSqSum_r[1]); 


  Gdiff_l[0] = -0.0125*(53.66563145999495*nuVtSqSum_l[1]*fl[7]-53.66563145999495*nuVtSqSum_l[1]*fc[7]+95.26279441628826*nuVtSqSum_l[2]*fl[6]+95.26279441628826*nuVtSqSum_l[2]*fc[6]+53.66563145999495*nuVtSqSum_l[0]*fl[5]-53.66563145999495*nuVtSqSum_l[0]*fc[5]+75.0*nuVtSqSum_l[2]*fl[4]-75.0*nuVtSqSum_l[2]*fc[4]+95.26279441628825*nuVtSqSum_l[1]*fl[3]+95.26279441628825*nuVtSqSum_l[1]*fc[3]+95.26279441628825*nuVtSqSum_l[0]*fl[2]+95.26279441628825*nuVtSqSum_l[0]*fc[2]+(75.0*fl[1]-75.0*fc[1])*nuVtSqSum_l[1]+(75.0*fl[0]-75.0*fc[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = -0.0025*((240.0*nuVtSqSum_l[2]+268.3281572999747*nuVtSqSum_l[0])*fl[7]+((-240.0*nuVtSqSum_l[2])-268.3281572999747*nuVtSqSum_l[0])*fc[7]+426.0281680828158*nuVtSqSum_l[1]*fl[6]+426.0281680828158*nuVtSqSum_l[1]*fc[6]+268.3281572999748*nuVtSqSum_l[1]*fl[5]-268.3281572999748*nuVtSqSum_l[1]*fc[5]+335.4101966249686*nuVtSqSum_l[1]*fl[4]-335.4101966249686*nuVtSqSum_l[1]*fc[4]+(426.0281680828159*nuVtSqSum_l[2]+476.3139720814412*nuVtSqSum_l[0])*fl[3]+(426.0281680828159*nuVtSqSum_l[2]+476.3139720814412*nuVtSqSum_l[0])*fc[3]+(335.4101966249686*fl[1]-335.4101966249686*fc[1])*nuVtSqSum_l[2]+476.3139720814412*nuVtSqSum_l[1]*fl[2]+476.3139720814412*nuVtSqSum_l[1]*fc[2]+(375.0*fl[0]-375.0*fc[0])*nuVtSqSum_l[1]+375.0*nuVtSqSum_l[0]*fl[1]-375.0*nuVtSqSum_l[0]*fc[1]); 
  Gdiff_l[2] = -3.571428571428571E-4*(1680.0*nuVtSqSum_l[1]*fl[7]-1680.0*nuVtSqSum_l[1]*fc[7]+(2130.140840414079*nuVtSqSum_l[2]+3334.19780457009*nuVtSqSum_l[0])*fl[6]+(2130.140840414079*nuVtSqSum_l[2]+3334.19780457009*nuVtSqSum_l[0])*fc[6]+1878.297101099823*nuVtSqSum_l[2]*fl[5]-1878.297101099823*nuVtSqSum_l[2]*fc[5]+(1677.050983124843*nuVtSqSum_l[2]+2625.0*nuVtSqSum_l[0])*fl[4]+((-1677.050983124843*nuVtSqSum_l[2])-2625.0*nuVtSqSum_l[0])*fc[4]+2982.197176579711*nuVtSqSum_l[1]*fl[3]+2982.197176579711*nuVtSqSum_l[1]*fc[3]+(3334.197804570088*fl[2]+3334.197804570088*fc[2]+2625.0*fl[0]-2625.0*fc[0])*nuVtSqSum_l[2]+(2347.87137637478*fl[1]-2347.87137637478*fc[1])*nuVtSqSum_l[1]); 


  Gdiff_r[0] = 0.0125*(53.66563145999495*nuVtSqSum_r[1]*fr[7]-53.66563145999495*nuVtSqSum_r[1]*fc[7]-95.26279441628826*nuVtSqSum_r[2]*fr[6]-95.26279441628826*nuVtSqSum_r[2]*fc[6]+53.66563145999495*nuVtSqSum_r[0]*fr[5]-53.66563145999495*nuVtSqSum_r[0]*fc[5]+75.0*nuVtSqSum_r[2]*fr[4]-75.0*nuVtSqSum_r[2]*fc[4]-95.26279441628825*nuVtSqSum_r[1]*fr[3]-95.26279441628825*nuVtSqSum_r[1]*fc[3]-95.26279441628825*nuVtSqSum_r[0]*fr[2]-95.26279441628825*nuVtSqSum_r[0]*fc[2]+(75.0*fr[1]-75.0*fc[1])*nuVtSqSum_r[1]+(75.0*fr[0]-75.0*fc[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = 0.0025*((240.0*nuVtSqSum_r[2]+268.3281572999747*nuVtSqSum_r[0])*fr[7]+((-240.0*nuVtSqSum_r[2])-268.3281572999747*nuVtSqSum_r[0])*fc[7]-426.0281680828158*nuVtSqSum_r[1]*fr[6]-426.0281680828158*nuVtSqSum_r[1]*fc[6]+268.3281572999748*nuVtSqSum_r[1]*fr[5]-268.3281572999748*nuVtSqSum_r[1]*fc[5]+335.4101966249686*nuVtSqSum_r[1]*fr[4]-335.4101966249686*nuVtSqSum_r[1]*fc[4]+((-426.0281680828159*nuVtSqSum_r[2])-476.3139720814412*nuVtSqSum_r[0])*fr[3]+((-426.0281680828159*nuVtSqSum_r[2])-476.3139720814412*nuVtSqSum_r[0])*fc[3]+(335.4101966249686*fr[1]-335.4101966249686*fc[1])*nuVtSqSum_r[2]-476.3139720814412*nuVtSqSum_r[1]*fr[2]-476.3139720814412*nuVtSqSum_r[1]*fc[2]+(375.0*fr[0]-375.0*fc[0])*nuVtSqSum_r[1]+375.0*nuVtSqSum_r[0]*fr[1]-375.0*nuVtSqSum_r[0]*fc[1]); 
  Gdiff_r[2] = 3.571428571428571E-4*(1680.0*nuVtSqSum_r[1]*fr[7]-1680.0*nuVtSqSum_r[1]*fc[7]+((-2130.140840414079*nuVtSqSum_r[2])-3334.19780457009*nuVtSqSum_r[0])*fr[6]+((-2130.140840414079*nuVtSqSum_r[2])-3334.19780457009*nuVtSqSum_r[0])*fc[6]+1878.297101099823*nuVtSqSum_r[2]*fr[5]-1878.297101099823*nuVtSqSum_r[2]*fc[5]+(1677.050983124843*nuVtSqSum_r[2]+2625.0*nuVtSqSum_r[0])*fr[4]+((-1677.050983124843*nuVtSqSum_r[2])-2625.0*nuVtSqSum_r[0])*fc[4]-2982.197176579711*nuVtSqSum_r[1]*fr[3]-2982.197176579711*nuVtSqSum_r[1]*fc[3]+((-3334.197804570088*fr[2])-3334.197804570088*fc[2]+2625.0*fr[0]-2625.0*fc[0])*nuVtSqSum_r[2]+(2347.87137637478*fr[1]-2347.87137637478*fc[1])*nuVtSqSum_r[1]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.7071067811865475*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaDrSurf_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.4517539514526256*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[2]+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[1]; 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.7071067811865475*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaDrSurf_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.4517539514526256*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[2]+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[1]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[5] += 4.743416490252569*Gdiff2_r[0]*rdvSq4+4.743416490252569*Gdiff2_l[0]*rdvSq4-1.581138830084189*Ghat_r[0]*rdv2+1.581138830084189*Ghat_l[0]*rdv2; 
  out[6] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 4.743416490252569*Gdiff2_r[1]*rdvSq4+4.743416490252569*Gdiff2_l[1]*rdvSq4-1.581138830084189*Ghat_r[1]*rdv2+1.581138830084189*Ghat_l[1]*rdv2; 
} 
