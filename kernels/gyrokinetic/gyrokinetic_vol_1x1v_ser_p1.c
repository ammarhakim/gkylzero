#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_i: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[4]; 
  hamil[0] = (0.3333333333333333*(3.0*rdvpar2Sq*(m_*wvparSq+1.414213562373095*phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 

  double BstarZdBmag[4]; 
  BstarZdBmag[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmag[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmag[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[4]; 
  alphax[0] = (0.8660254037844386*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.8660254037844386*BstarZdBmag[1]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.8660254037844386*BstarZdBmag[2]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.8660254037844386*hamil[2]*BstarZdBmag[3]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(0.8660254037844386*alphax[3]-0.4999999999999999*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.25*((-0.8660254037844386*alphax[3])+0.4999999999999999*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*((-0.8660254037844386*alphax[3])-0.4999999999999999*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.25*(0.8660254037844386*alphax[3]+0.4999999999999999*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  double alphavpar[4]; 
  alphavpar[0] = -(0.8660254037844386*BstarZdBmag[0]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(0.8660254037844386*BstarZdBmag[1]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.8660254037844386*hamil[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.8660254037844386*hamil[1]*BstarZdBmag[3]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*(0.8660254037844386*alphavpar[3]-0.8660254037844386*alphavpar[2]-0.4999999999999999*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.25*((-0.8660254037844386*(alphavpar[3]+alphavpar[2]))+0.4999999999999999*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*((-0.8660254037844386*alphavpar[3])+0.8660254037844386*alphavpar[2]-0.4999999999999999*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.25*(0.8660254037844386*(alphavpar[3]+alphavpar[2])+0.4999999999999999*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.8660254037844386*(alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[3] += 0.8660254037844386*((alphavpar[2]+alphax[1])*f[3]+f[1]*alphax[3]+f[2]*(alphavpar[3]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 

  return cflFreq; 
} 
GKYL_CU_DH double gyrokinetic_vol_step2_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[3] += -(1.224744871391589*(dApardt[0]*f[1]+f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.25*(0.7071067811865475*dApardt[0]-0.7071067811865474*dApardt[1])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.25*(0.7071067811865474*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.25*(0.7071067811865475*dApardt[0]-0.7071067811865474*dApardt[1])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.25*(0.7071067811865474*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
