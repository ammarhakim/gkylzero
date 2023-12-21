#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nI, const double *vnu, const double *vsqnu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: cell-center coordinates. 
  // dxv[2]: cell spacing. 
  // nI: atomic density 
  // vnu: 2/pi*v*nu(v) dg field representation (v'(v||,mu) in notes) 
  // vsqnu: sqrt(mu*me/2B)*v^2*nu(v) dg field representation (v''(v||,mu) in notes) 
  // f: input distribution function.
  // out: incremented output 

  double rdv2[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  printf("\n Volume: w[0]=%e, w[1]= %e, rdv= %e, f[0]=%e,",w[0],w[1],rdv2[0],f[0]);
  double alphaDrag[6] = {0,0}; 
  alphaDrag[0] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[1]+0.7071067811865475*nI[0]*rdv2[0]*vnu[0]; 
  alphaDrag[1] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[1]+0.7071067811865475*rdv2[0]*vnu[0]*nI[1]; 
  alphaDrag[2] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[3]+0.7071067811865475*nI[0]*rdv2[0]*vnu[2]; 
  alphaDrag[3] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[3]+0.7071067811865475*rdv2[0]*nI[1]*vnu[2]; 
  alphaDrag[4] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[5]+0.7071067811865475*nI[0]*rdv2[0]*vnu[4]; 
  alphaDrag[5] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[5]+0.7071067811865475*rdv2[0]*nI[1]*vnu[4]; 

  out[2] += 0.8660254037844386*(alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[4]*f[5]+f[4]*alphaDrag[5]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[4] += 1.732050807568877*(alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4])+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[5] += 1.732050807568877*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  /*  printf("\n Volume: w[0]=%e, w[1]= %e, rdv= %e, f[0]=%e,",w[0],w[1],rdv2[0],f[0]);
  for (int i=0;i<6;i++){
    printf("out[%i]=%e ",i,out[i]);
  }
  printf("\n");
  printf("     ");
  for (int i=0;i<6;i++){
    printf("nufield[%i]=%e ",i,vnu[i]);
  }
  printf("\n\n");*/
  return fabs(1.25*alphaDrag[0]-1.397542485937369*alphaDrag[4]);
} 
