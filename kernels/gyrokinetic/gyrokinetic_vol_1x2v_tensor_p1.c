#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out) 
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
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.6123724356957944*BstarZdBmag[1]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.6123724356957944*BstarZdBmag[2]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.6123724356957944*hamil[2]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932736*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932736*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932736*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932736*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932736*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932736*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932736*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932736*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  double alphavpar[8]; 
  alphavpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(0.6123724356957944*BstarZdBmag[1]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.6123724356957944*hamil[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.6123724356957944*BstarZdBmag[0]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[4] = -(0.6123724356957944*hamil[1]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.6123724356957944*BstarZdBmag[1]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.6123724356957944*BstarZdBmag[2]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[7] = -(0.6123724356957944*BstarZdBmag[4]*hamil[5]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.6123724356957942*alphavpar[7])+0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932736*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6123724356957942*alphavpar[7]-0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932736*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6123724356957942*alphavpar[7]+0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932736*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957942*alphavpar[7])-0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932736*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.6123724356957942*alphavpar[7]-0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932736*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957942*alphavpar[7])+0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932736*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957942*alphavpar[7])-0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932736*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.6123724356957942*alphavpar[7]+0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932736*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932736*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.6123724356957944*(alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphavpar[6]*f[7]+f[6]*alphavpar[7]+alphavpar[3]*f[5]+f[3]*alphavpar[5]+(alphavpar[2]+alphax[1])*f[4]+f[1]*alphax[4]+f[2]*(alphavpar[4]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[5] += 0.6123724356957944*(alphax[4]*f[7]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphavpar[4]*f[7]+f[4]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*((alphavpar[2]+alphax[1])*f[7]+f[2]*alphavpar[7]+(alphavpar[4]+alphax[0])*f[6]+f[4]*alphavpar[6]+(alphax[4]+alphavpar[0])*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*f[3]+f[1]*alphavpar[3]); 

  return cflFreq; 
} 
GKYL_CU_DH double gyrokinetic_vol_step2_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[4] += -(1.224744871391589*(dApardt[0]*f[1]+f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[6] += -(1.224744871391589*(dApardt[1]*f[5]+dApardt[0]*f[3])*q_*rdvpar2)/m_; 
  out[7] += -(1.224744871391589*(dApardt[0]*f[5]+dApardt[1]*f[3])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865474*dApardt[1])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865474*dApardt[1])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865474*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865474*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865474*dApardt[1])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865474*dApardt[1])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865474*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865474*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
