#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out) 
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
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[8] = {0.}; 
  hamil[0] = (0.3333333333333333*(3.0*rdvpar2Sq*(m_*wvparSq+1.414213562373095*phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = 1.414213562373095*phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double BstarZdBmag[8] = {0.}; 
  BstarZdBmag[0] = (1.732050807568877*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobtot_inv[2]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_)/q_; 
  BstarZdBmag[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobtot_inv[2]+11.18033988749895*jacobtot_inv[0])+5.0*b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobtot_inv[2]+jacobtot_inv[1]*cmag[2])+5.0*(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = ((2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = ((b_y[2]*(2.0*jacobtot_inv[2]+2.23606797749979*jacobtot_inv[0])+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobtot_inv[2]+7.0*(5.0*jacobtot_inv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobtot_inv[1]))*q_))/q_; 
  BstarZdBmag[6] = (1.0*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 

  double alphax[8] = {0.}; 
  alphax[0] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[5]+BstarZdBmag[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[5]+BstarZdBmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[0]*hamil[5]+BstarZdBmag[2]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[5]+hamil[2]*BstarZdBmag[3])*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.5*(3.872983346207417*hamil[5]*BstarZdBmag[6]+1.732050807568877*hamil[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphax[5] = (1.732050807568877*BstarZdBmag[2]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.5*(1.732050807568877*hamil[2]*BstarZdBmag[6]+3.872983346207417*BstarZdBmag[4]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphax[7] = (1.732050807568877*BstarZdBmag[3]*hamil[5]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-0.7745966692414833*alphax[7])-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*alphax[7]-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.25*((-0.7745966692414833*alphax[7])+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.7745966692414833*alphax[7]-1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]-1.161895003862225*alphax[3]-0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*alphax[7])-0.5590169943749475*alphax[5]+1.118033988749895*alphax[4]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.25*(0.7745966692414833*alphax[7]+1.5*alphax[6]+0.4472135954999579*alphax[5]+1.118033988749895*alphax[4]+1.161895003862225*alphax[3]+0.6708203932499369*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  double alphavpar[8] = {0.}; 
  alphavpar[0] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[1]*hamil[4]+BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(1.732050807568877*(0.5*(2.23606797749979*BstarZdBmag[0]*hamil[4]+BstarZdBmag[1]*hamil[1])+BstarZdBmag[4]*hamil[4])*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.8660254037844386*(2.23606797749979*BstarZdBmag[3]*hamil[4]+hamil[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(1.0*(1.732050807568877*hamil[4]*BstarZdBmag[6]+0.8660254037844386*(2.23606797749979*BstarZdBmag[2]*hamil[4]+hamil[1]*BstarZdBmag[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (1.732050807568877*((-1.0*BstarZdBmag[1]*hamil[4])-0.5*hamil[1]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphavpar[6] = (((-0.8660254037844386*hamil[1]*BstarZdBmag[6])-1.732050807568877*BstarZdBmag[3]*hamil[4])*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.25*((-0.7745966692414833*alphavpar[6])+0.4472135954999579*alphavpar[4]+1.161895003862225*alphavpar[3]-0.8660254037844386*alphavpar[2]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*alphavpar[6]-0.5590169943749475*alphavpar[4]-0.8660254037844386*alphavpar[2]+0.5*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.25*((-0.7745966692414833*alphavpar[6])+0.4472135954999579*alphavpar[4]-1.161895003862225*alphavpar[3]-0.8660254037844386*alphavpar[2]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.25*(0.7745966692414833*alphavpar[6]+0.4472135954999579*alphavpar[4]-1.161895003862225*alphavpar[3]+0.8660254037844386*alphavpar[2]-0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*alphavpar[6])-0.5590169943749475*alphavpar[4]+0.8660254037844386*alphavpar[2]+0.5*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.25*(0.7745966692414833*alphavpar[6]+0.4472135954999579*alphavpar[4]+1.161895003862225*alphavpar[3]+0.8660254037844386*alphavpar[2]+0.6708203932499369*alphavpar[1]+0.5*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.8660254037844386*(alphax[7]*fin[7]+alphax[6]*fin[6]+alphax[5]*fin[5]+alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.8660254037844386*(alphavpar[6]*fin[6]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[3] += 0.1*(7.745966692414834*(alphax[3]*fin[7]+fin[3]*alphax[7])+(8.660254037844387*alphax[4]+7.745966692414834*alphavpar[3])*fin[6]+8.660254037844387*fin[4]*alphax[6]+7.745966692414834*(fin[3]*alphavpar[6]+alphax[2]*fin[5]+fin[2]*alphax[5]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4])+8.660254037844386*((alphavpar[2]+alphax[1])*fin[3]+fin[1]*alphax[3]+fin[2]*(alphavpar[3]+alphax[0])+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1])); 
  out[4] += 0.1*(19.36491673103708*(alphax[5]*fin[7]+fin[5]*alphax[7])+17.32050807568877*(alphax[3]*fin[6]+fin[3]*alphax[6])+17.32050807568877*(alphax[1]*fin[4]+fin[1]*alphax[4])+19.36491673103709*(alphax[2]*fin[3]+fin[2]*alphax[3]+alphax[0]*fin[1]+fin[0]*alphax[1])); 
  out[5] += 0.1*(17.32050807568877*alphavpar[3]*fin[7]+19.36491673103708*(alphavpar[4]*fin[6]+fin[4]*alphavpar[6])+17.32050807568877*alphavpar[2]*fin[5]+19.36491673103709*(alphavpar[1]*fin[3]+fin[1]*alphavpar[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2])); 
  out[6] += 0.01428571428571429*((108.4435336938077*alphax[6]+121.2435565298214*alphax[2])*fin[7]+(108.4435336938077*fin[6]+121.2435565298214*fin[2])*alphax[7]+(38.72983346207417*alphavpar[6]+60.6217782649107*alphavpar[2])*fin[6]+121.2435565298214*(alphax[1]*fin[6]+fin[1]*alphax[6])+60.6217782649107*fin[2]*alphavpar[6]+121.2435565298214*(alphax[3]*fin[5]+fin[3]*alphax[5])+(38.72983346207417*alphavpar[4]+121.2435565298214*alphax[3]+60.62177826491071*alphavpar[0])*fin[4]+121.2435565298214*fin[3]*alphax[4]+60.62177826491071*fin[0]*alphavpar[4]+54.22176684690384*alphavpar[3]*fin[3]+135.5544171172596*(alphax[0]*fin[3]+fin[0]*alphax[3]+alphax[1]*fin[2])+fin[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphavpar[1])); 
  out[7] += 0.01428571428571429*((38.72983346207417*alphax[7]+108.4435336938077*alphavpar[6]+121.2435565298214*alphavpar[2])*fin[7]+60.6217782649107*(alphax[1]*fin[7]+fin[1]*alphax[7])+54.22176684690384*alphax[6]*fin[6]+121.2435565298214*(alphavpar[1]*fin[6]+fin[1]*alphavpar[6])+(38.72983346207417*alphax[5]+121.2435565298214*alphavpar[3])*fin[5]+60.62177826491071*(alphax[0]*fin[5]+fin[0]*alphax[5])+121.2435565298214*alphavpar[3]*fin[4]+fin[3]*(121.2435565298214*alphavpar[4]+54.22176684690384*alphax[3])+135.5544171172596*(alphavpar[0]*fin[3]+fin[0]*alphavpar[3])+54.22176684690384*alphax[2]*fin[2]+135.5544171172596*(alphavpar[1]*fin[2]+fin[1]*alphavpar[2])); 

  return cflFreq; 
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
