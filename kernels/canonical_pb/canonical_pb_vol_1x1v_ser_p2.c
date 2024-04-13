#include <gkyl_canonical_pb_kernels.h> 
GKYL_CU_DH double canonical_pb_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // fin: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvx = w[1];
  double rdvx2 = 2.0/dxv[1];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvxSq = w[1]*w[1];
  double rdvx2Sq = rdvx2*rdvx2;

  double alphax[8] = {0.}; 
  alphax[0] = 1.732050807568877*hamil[2]*rdvx2*rdx2; 
  alphax[1] = 1.732050807568877*hamil[3]*rdvx2*rdx2; 
  alphax[2] = 3.872983346207417*hamil[5]*rdvx2*rdx2; 
  alphax[3] = 3.872983346207417*hamil[7]*rdvx2*rdx2; 
  alphax[4] = 1.732050807568877*hamil[6]*rdvx2*rdx2; 

  double alphavx[8] = {0.}; 
  alphavx[0] = -1.732050807568877*hamil[1]*rdvx2*rdx2; 
  alphavx[1] = -3.872983346207417*hamil[4]*rdvx2*rdx2; 
  alphavx[2] = -1.732050807568877*hamil[3]*rdvx2*rdx2; 
  alphavx[3] = -3.872983346207417*hamil[6]*rdvx2*rdx2; 
  alphavx[5] = -1.732050807568877*hamil[7]*rdvx2*rdx2; 

  out[1] += 0.8660254037844386*(alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.8660254037844386*(alphavx[5]*fin[5]+alphavx[3]*fin[3]+alphavx[2]*fin[2]+alphavx[1]*fin[1]+alphavx[0]*fin[0]); 
  out[3] += 0.1*((8.660254037844387*alphavx[5]+7.745966692414834*alphax[3])*fin[7]+8.660254037844387*alphax[4]*fin[6]+7.745966692414834*(alphavx[3]*fin[6]+alphax[2]*fin[5]+alphavx[1]*fin[4])+8.660254037844386*((alphavx[2]+alphax[1])*fin[3]+fin[1]*alphax[3]+fin[2]*(alphavx[3]+alphax[0])+fin[0]*alphax[2]+alphavx[0]*fin[1]+fin[0]*alphavx[1])); 
  out[4] += 0.1*(17.32050807568877*alphax[3]*fin[6]+17.32050807568877*(alphax[1]*fin[4]+fin[1]*alphax[4])+19.36491673103709*(alphax[2]*fin[3]+fin[2]*alphax[3]+alphax[0]*fin[1]+fin[0]*alphax[1])); 
  out[5] += 0.1*(17.32050807568877*alphavx[3]*fin[7]+17.32050807568877*(alphavx[2]*fin[5]+fin[2]*alphavx[5])+19.36491673103709*(alphavx[1]*fin[3]+fin[1]*alphavx[3]+alphavx[0]*fin[2]+fin[0]*alphavx[2])); 
  out[6] += 0.1*(17.32050807568877*alphax[2]*fin[7]+(8.660254037844386*alphavx[2]+17.32050807568877*alphax[1])*fin[6]+17.32050807568877*alphax[3]*fin[5]+(17.32050807568877*alphax[3]+8.660254037844387*alphavx[0])*fin[4]+fin[3]*(17.32050807568877*alphax[4]+7.745966692414834*alphavx[3])+19.36491673103708*(alphax[0]*fin[3]+fin[0]*alphax[3]+alphax[1]*fin[2])+fin[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphavx[1])); 
  out[7] += 0.1*((17.32050807568877*alphavx[2]+8.660254037844386*alphax[1])*fin[7]+17.32050807568877*alphavx[1]*fin[6]+(17.32050807568877*alphavx[3]+8.660254037844387*alphax[0])*fin[5]+17.32050807568877*(fin[3]*alphavx[5]+alphavx[3]*fin[4])+7.745966692414834*alphax[3]*fin[3]+19.36491673103708*(alphavx[0]*fin[3]+fin[0]*alphavx[3])+7.745966692414834*alphax[2]*fin[2]+19.36491673103708*(alphavx[1]*fin[2]+fin[1]*alphavx[2])); 

  return 0.; 
} 
