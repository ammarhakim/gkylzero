#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *Bstar_Bmag, const double *fin, double* GKYL_RESTRICT out) 
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
  // Bstar_Bmag: Bstar/Bmag volume expansion, pre-computed time-independent part.
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
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  const double *BstarZdBmag = &Bstar_Bmag[0]; 

  double hamil[8] = {0.}; 
  hamil[0] = m_*wvparSq+(0.3333333333333333*m_)/rdvpar2Sq+1.414213562373095*phi[0]*q_; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = 1.414213562373095*phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double alphax[8] = {0.}; 
  alphax[0] = (1.224744871391589*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (1.224744871391589*BstarZdBmag[1]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (2.738612787525831*BstarZdBmag[0]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[3] = (2.738612787525831*BstarZdBmag[1]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[4] = (1.224744871391589*BstarZdBmag[2]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[6] = (2.738612787525831*BstarZdBmag[2]*hamil[5]*rdvpar2*rdx2)/m_; 

  double alphavpar[8] = {0.}; 
  alphavpar[0] = (((-2.738612787525831*BstarZdBmag[1]*hamil[4])-1.224744871391589*BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = ((((-2.449489742783178*BstarZdBmag[2])-2.738612787525831*BstarZdBmag[0])*hamil[4]-1.224744871391589*BstarZdBmag[1]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (((-2.449489742783178*BstarZdBmag[1]*hamil[4])-1.224744871391589*hamil[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 

  out[1] += 0.8660254037844386*(alphax[6]*fin[6]+alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.8660254037844386*(alphavpar[4]*fin[4]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.7745966692414834*alphax[3]*fin[7]+0.8660254037844387*(alphax[4]*fin[6]+fin[4]*alphax[6])+0.7745966692414833*(alphax[2]*fin[5]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4])+0.8660254037844386*(alphax[1]*fin[3]+fin[1]*alphax[3]+alphax[0]*fin[2]+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[4] += 1.732050807568877*(alphax[3]*fin[6]+fin[3]*alphax[6])+1.732050807568877*(alphax[1]*fin[4]+fin[1]*alphax[4])+1.936491673103709*(alphax[2]*fin[3]+fin[2]*alphax[3]+alphax[0]*fin[1]+fin[0]*alphax[1]); 
  out[5] += 1.936491673103709*(alphavpar[4]*fin[6]+alphavpar[1]*fin[3]+alphavpar[0]*fin[2]); 
  out[6] += 1.549193338482967*alphax[6]*fin[7]+1.732050807568877*(alphax[2]*fin[7]+alphax[1]*fin[6]+fin[1]*alphax[6])+1.732050807568877*alphax[3]*fin[5]+(0.5532833351724881*alphavpar[4]+1.732050807568877*alphax[3]+0.8660254037844387*alphavpar[0])*fin[4]+1.732050807568877*fin[3]*alphax[4]+0.8660254037844387*fin[0]*alphavpar[4]+1.936491673103709*(alphax[0]*fin[3]+fin[0]*alphax[3]+alphax[1]*fin[2])+fin[1]*(1.936491673103709*alphax[2]+0.7745966692414834*alphavpar[1]); 
  out[7] += 0.8660254037844386*alphax[1]*fin[7]+(0.7745966692414834*alphax[6]+1.732050807568877*alphavpar[1])*fin[6]+0.8660254037844387*alphax[0]*fin[5]+fin[3]*(1.732050807568877*alphavpar[4]+0.7745966692414834*alphax[3]+1.936491673103709*alphavpar[0])+(0.7745966692414834*alphax[2]+1.936491673103709*alphavpar[1])*fin[2]; 

  return 0.; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_1x1v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(apardot[2]*fin[4]+apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*apardot[1]*fin[4]+4.47213595499958*fin[1]*apardot[2]+5.0*apardot[0]*fin[1]+5.0*fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[5] += -(0.7071067811865475*(3.872983346207417*apardot[2]*fin[6]+3.872983346207417*apardot[1]*fin[3]+3.872983346207417*apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[6] += -(0.07824607964359516*(10.0*apardot[2]*fin[4]+15.65247584249853*apardot[0]*fin[4]+15.65247584249853*fin[0]*apardot[2]+14.0*apardot[1]*fin[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.1414213562373095*(17.32050807568877*apardot[1]*fin[6]+17.32050807568877*apardot[2]*fin[3]+19.36491673103708*apardot[0]*fin[3]+19.36491673103708*apardot[1]*fin[2])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.25*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.25*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.25*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.25*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.25*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.25*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
