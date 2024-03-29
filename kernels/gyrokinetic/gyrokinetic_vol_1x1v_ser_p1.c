#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
  // fin: Distribution function.
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

  double hamil[6] = {0.}; 
  hamil[0] = (0.3333333333333333*(3.0*rdvpar2Sq*(m_*wvparSq+1.414213562373095*phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double BstarZdBmag[6] = {0.}; 
  BstarZdBmag[0] = (1.732050807568877*jacobtot_inv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_)/q_; 
  BstarZdBmag[1] = (1.732050807568877*b_y[1]*jacobtot_inv[1]*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_)/q_; 
  BstarZdBmag[2] = (jacobtot_inv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = (b_y[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double alphax[6] = {0.}; 
  alphax[0] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[4]+BstarZdBmag[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[4]+BstarZdBmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[0]*hamil[4]+BstarZdBmag[2]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[4]+hamil[2]*BstarZdBmag[3])*rdvpar2*rdx2)/m_; 
  alphax[4] = (1.732050807568877*BstarZdBmag[2]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[5] = (1.732050807568877*BstarZdBmag[3]*hamil[4]*rdvpar2*rdx2)/m_; 


  double alphavpar[6] = {0.}; 
  alphavpar[0] = -(0.8660254037844386*BstarZdBmag[0]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(0.8660254037844386*BstarZdBmag[1]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.8660254037844386*hamil[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.8660254037844386*hamil[1]*BstarZdBmag[3]*rdvpar2*rdx2)/m_; 


  out[1] += 0.8660254037844386*(alphax[5]*fin[5]+alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.8660254037844386*(alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[3]*fin[5]+fin[3]*alphax[5]+alphax[2]*fin[4]+fin[2]*alphax[4])+8.660254037844386*((alphavpar[2]+alphax[1])*fin[3]+fin[1]*alphax[3]+fin[2]*(alphavpar[3]+alphax[0])+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1])); 
  out[4] += 0.1*(17.32050807568877*alphavpar[3]*fin[5]+17.32050807568877*alphavpar[2]*fin[4]+19.36491673103709*(alphavpar[1]*fin[3]+fin[1]*alphavpar[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2])); 
  out[5] += 0.01428571428571429*((38.72983346207417*alphax[5]+121.2435565298214*alphavpar[2])*fin[5]+60.6217782649107*(alphax[1]*fin[5]+fin[1]*alphax[5])+(38.72983346207417*alphax[4]+121.2435565298214*alphavpar[3])*fin[4]+60.62177826491071*(alphax[0]*fin[4]+fin[0]*alphax[4])+54.22176684690384*alphax[3]*fin[3]+135.5544171172596*(alphavpar[0]*fin[3]+fin[0]*alphavpar[3])+54.22176684690384*alphax[2]*fin[2]+135.5544171172596*(alphavpar[1]*fin[2]+fin[1]*alphavpar[2])); 

  return 0.; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[3] += -(1.224744871391589*(apardot[0]*fin[1]+fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[4] += -(2.738612787525831*(apardot[1]*fin[3]+apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[5] += -(2.738612787525831*(apardot[0]*fin[3]+apardot[1]*fin[2])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.25*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.1767766952966369*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.25*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.1767766952966369*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
