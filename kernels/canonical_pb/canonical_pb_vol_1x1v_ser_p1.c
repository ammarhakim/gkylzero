#include <gkyl_canonical_pb_kernels.h> 
GKYL_CU_DH double canonical_pb_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, const double *fin, double* GKYL_RESTRICT out) 
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

  double alphax[6] = {0.}; 
  alphax[0] = 1.732050807568877*hamil[2]*rdvx2*rdx2Sq; 
  alphax[1] = 1.732050807568877*hamil[3]*rdvx2*rdx2Sq; 
  alphax[2] = 3.872983346207417*hamil[4]*rdvx2*rdx2Sq; 
  alphax[3] = 3.872983346207417*hamil[5]*rdvx2*rdx2Sq; 

  double alphavx[6] = {0.}; 
  alphavx[0] = -1.732050807568877*hamil[1]*rdvx2Sq*rdx2; 
  alphavx[2] = -1.732050807568877*hamil[3]*rdvx2Sq*rdx2; 
  alphavx[4] = -1.732050807568877*hamil[5]*rdvx2Sq*rdx2; 

  out[1] += 0.8660254037844386*(alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.8660254037844386*(alphavx[4]*fin[4]+alphavx[2]*fin[2]+alphavx[0]*fin[0]); 
  out[3] += 0.1*(8.660254037844387*alphavx[4]*fin[5]+7.745966692414834*(alphax[3]*fin[5]+alphax[2]*fin[4])+8.660254037844386*((alphavx[2]+alphax[1])*fin[3]+fin[1]*alphax[3]+alphax[0]*fin[2]+fin[0]*alphax[2]+alphavx[0]*fin[1])); 
  out[4] += 0.8660254037844386*(2.0*(alphavx[2]*fin[4]+fin[2]*alphavx[4])+2.23606797749979*(alphavx[0]*fin[2]+fin[0]*alphavx[2])); 
  out[5] += 0.1*((17.32050807568877*alphavx[2]+8.660254037844386*alphax[1])*fin[5]+8.660254037844387*alphax[0]*fin[4]+fin[3]*(17.32050807568877*alphavx[4]+7.745966692414834*alphax[3]+19.36491673103708*alphavx[0])+7.745966692414834*alphax[2]*fin[2]+19.36491673103708*fin[1]*alphavx[2]); 

  return 0.; 
} 
