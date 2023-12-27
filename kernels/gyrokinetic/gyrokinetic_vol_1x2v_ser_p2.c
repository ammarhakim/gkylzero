#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out) 
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
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wvparSq = wvpar*wvpar;
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[20] = {0.}; 
  hamil[0] = 1.414213562373095*m_*wvparSq+2.0*bmag[0]*wmu+(0.4714045207910317*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double alphax[20] = {0.}; 
  alphax[0] = ((0.8660254037844386*cmag[2]*hamil[2]*jacobtot_inv[2]+0.8660254037844386*cmag[1]*jacobtot_inv[1]*hamil[2]+0.8660254037844386*cmag[0]*jacobtot_inv[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = ((0.7745966692414833*cmag[1]*hamil[2]*jacobtot_inv[2]+0.7745966692414833*jacobtot_inv[1]*cmag[2]*hamil[2]+0.8660254037844386*cmag[0]*jacobtot_inv[1]*hamil[2]+0.8660254037844386*jacobtot_inv[0]*cmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = ((1.936491673103709*cmag[2]*jacobtot_inv[2]*hamil[8]+1.936491673103709*cmag[1]*jacobtot_inv[1]*hamil[8]+1.936491673103709*cmag[0]*jacobtot_inv[0]*hamil[8])*rdvpar2*rdx2)/m_; 
  alphax[4] = ((1.732050807568877*cmag[1]*jacobtot_inv[2]*hamil[8]+1.732050807568877*jacobtot_inv[1]*cmag[2]*hamil[8]+1.936491673103709*cmag[0]*jacobtot_inv[1]*hamil[8]+1.936491673103709*jacobtot_inv[0]*cmag[1]*hamil[8])*rdvpar2*rdx2)/m_; 
  alphax[7] = ((0.5532833351724881*cmag[2]*hamil[2]*jacobtot_inv[2]+0.8660254037844386*cmag[0]*hamil[2]*jacobtot_inv[2]+0.8660254037844386*jacobtot_inv[0]*cmag[2]*hamil[2]+0.7745966692414833*cmag[1]*jacobtot_inv[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[11] = ((1.237179148263484*cmag[2]*jacobtot_inv[2]*hamil[8]+1.936491673103709*cmag[0]*jacobtot_inv[2]*hamil[8]+1.936491673103709*jacobtot_inv[0]*cmag[2]*hamil[8]+1.732050807568877*cmag[1]*jacobtot_inv[1]*hamil[8])*rdvpar2*rdx2)/m_; 

  double alphavpar[20] = {0.}; 
  alphavpar[0] = (((-1.732050807568877*cmag[1]*jacobtot_inv[2]*hamil[7])-1.732050807568877*jacobtot_inv[1]*cmag[2]*hamil[7]-1.936491673103709*cmag[0]*jacobtot_inv[1]*hamil[7]-1.936491673103709*jacobtot_inv[0]*cmag[1]*hamil[7]-0.8660254037844386*hamil[1]*cmag[2]*jacobtot_inv[2]-0.8660254037844386*cmag[1]*hamil[1]*jacobtot_inv[1]-0.8660254037844386*cmag[0]*jacobtot_inv[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = (((-3.043058343448684*cmag[2]*jacobtot_inv[2]*hamil[7])-1.732050807568877*cmag[0]*jacobtot_inv[2]*hamil[7]-1.732050807568877*jacobtot_inv[0]*cmag[2]*hamil[7]-3.485685011586674*cmag[1]*jacobtot_inv[1]*hamil[7]-1.936491673103709*cmag[0]*jacobtot_inv[0]*hamil[7]-0.7745966692414833*cmag[1]*hamil[1]*jacobtot_inv[2]-0.7745966692414833*hamil[1]*jacobtot_inv[1]*cmag[2]-0.8660254037844386*cmag[0]*hamil[1]*jacobtot_inv[1]-0.8660254037844386*jacobtot_inv[0]*cmag[1]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[3] = (((-1.732050807568877*cmag[1]*jacobtot_inv[2]*hamil[13])-1.732050807568877*jacobtot_inv[1]*cmag[2]*hamil[13]-1.936491673103709*cmag[0]*jacobtot_inv[1]*hamil[13]-1.936491673103709*jacobtot_inv[0]*cmag[1]*hamil[13]-0.8660254037844386*cmag[2]*jacobtot_inv[2]*hamil[5]-0.8660254037844386*cmag[1]*jacobtot_inv[1]*hamil[5]-0.8660254037844386*cmag[0]*jacobtot_inv[0]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[5] = (((-3.043058343448685*cmag[2]*jacobtot_inv[2]*hamil[13])-1.732050807568877*cmag[0]*jacobtot_inv[2]*hamil[13]-1.732050807568877*jacobtot_inv[0]*cmag[2]*hamil[13]-3.485685011586675*cmag[1]*jacobtot_inv[1]*hamil[13]-1.936491673103709*cmag[0]*jacobtot_inv[0]*hamil[13]-0.7745966692414833*cmag[1]*jacobtot_inv[2]*hamil[5]-0.7745966692414833*jacobtot_inv[1]*cmag[2]*hamil[5]-0.8660254037844386*cmag[0]*jacobtot_inv[1]*hamil[5]-0.8660254037844386*jacobtot_inv[0]*cmag[1]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[7] = (((-3.043058343448684*cmag[1]*jacobtot_inv[2]*hamil[7])-3.043058343448684*jacobtot_inv[1]*cmag[2]*hamil[7]-1.732050807568877*cmag[0]*jacobtot_inv[1]*hamil[7]-1.732050807568877*jacobtot_inv[0]*cmag[1]*hamil[7]-0.5532833351724881*hamil[1]*cmag[2]*jacobtot_inv[2]-0.8660254037844386*cmag[0]*hamil[1]*jacobtot_inv[2]-0.8660254037844386*jacobtot_inv[0]*hamil[1]*cmag[2]-0.7745966692414833*cmag[1]*hamil[1]*jacobtot_inv[1])*rdvpar2*rdx2)/m_; 
  alphavpar[13] = (((-3.043058343448684*cmag[1]*jacobtot_inv[2]*hamil[13])-3.043058343448684*jacobtot_inv[1]*cmag[2]*hamil[13]-1.732050807568877*cmag[0]*jacobtot_inv[1]*hamil[13]-1.732050807568877*jacobtot_inv[0]*cmag[1]*hamil[13]-0.5532833351724881*cmag[2]*jacobtot_inv[2]*hamil[5]-0.8660254037844387*cmag[0]*jacobtot_inv[2]*hamil[5]-0.8660254037844387*jacobtot_inv[0]*cmag[2]*hamil[5]-0.7745966692414834*cmag[1]*jacobtot_inv[1]*hamil[5])*rdvpar2*rdx2)/m_; 

  out[1] += 0.6123724356957944*(alphax[11]*fin[11]+alphax[7]*fin[7]+alphax[4]*fin[4]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.6123724356957944*(alphavpar[13]*fin[13]+alphavpar[7]*fin[7]+alphavpar[5]*fin[5]+alphavpar[3]*fin[3]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.5477225575051661*(alphavpar[5]*fin[13]+fin[5]*alphavpar[13]+alphax[4]*fin[12])+0.6123724356957944*(alphax[7]*fin[11]+fin[7]*alphax[11])+0.5477225575051661*(alphax[2]*fin[8]+alphavpar[1]*fin[7]+fin[1]*alphavpar[7])+0.6123724356957944*(alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+alphax[1]*fin[4]+fin[1]*alphax[4]+alphax[0]*fin[2]+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[5] += 0.6123724356957944*(alphax[11]*fin[17]+alphax[7]*fin[13]+alphax[4]*fin[10]+alphax[2]*fin[6]+alphax[1]*fin[5]+alphax[0]*fin[3]); 
  out[6] += 0.5477225575051661*alphavpar[5]*fin[15]+0.6123724356957944*(alphavpar[7]*fin[13]+fin[7]*alphavpar[13])+0.5477225575051661*alphavpar[3]*fin[9]+0.6123724356957944*(alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 1.224744871391589*(alphax[4]*fin[11]+fin[4]*alphax[11]+alphax[1]*fin[7]+fin[1]*alphax[7])+1.369306393762915*(alphax[2]*fin[4]+fin[2]*alphax[4]+alphax[0]*fin[1]+fin[0]*alphax[1]); 
  out[8] += 1.369306393762915*(alphavpar[13]*fin[17]+alphavpar[7]*fin[11]+alphavpar[5]*fin[10]+alphavpar[3]*fin[6]+alphavpar[1]*fin[4]+alphavpar[0]*fin[2]); 
  out[10] += 0.5477225575051661*alphax[4]*fin[18]+0.6123724356957944*alphax[7]*fin[17]+0.4898979485566357*alphavpar[13]*fin[15]+0.5477225575051661*(alphavpar[3]*fin[15]+alphax[2]*fin[14])+0.6123724356957944*alphax[11]*fin[13]+0.5477225575051661*(alphavpar[1]*fin[13]+fin[1]*alphavpar[13])+0.6123724356957944*alphax[1]*fin[10]+0.5477225575051661*(alphavpar[5]*(fin[9]+fin[7])+fin[5]*alphavpar[7])+0.6123724356957944*(alphax[0]*fin[6]+(alphax[4]+alphavpar[0])*fin[5]+fin[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*fin[3]+fin[1]*alphavpar[3]); 
  out[11] += 0.3912303982179757*alphavpar[13]*fin[13]+0.6123724356957944*(alphavpar[3]*fin[13]+fin[3]*alphavpar[13])+1.095445115010332*alphax[11]*fin[12]+1.224744871391589*(alphax[2]*fin[12]+alphax[1]*fin[11]+fin[1]*alphax[11]+alphax[4]*fin[8])+(0.3912303982179757*alphavpar[7]+1.224744871391589*alphax[4]+0.6123724356957944*alphavpar[0])*fin[7]+1.224744871391589*fin[4]*alphax[7]+0.6123724356957944*fin[0]*alphavpar[7]+0.5477225575051661*alphavpar[5]*fin[5]+1.369306393762915*(alphax[0]*fin[4]+fin[0]*alphax[4]+alphax[1]*fin[2])+fin[1]*(1.369306393762915*alphax[2]+0.5477225575051661*alphavpar[1]); 
  out[12] += 1.224744871391589*(alphavpar[5]*fin[17]+fin[10]*alphavpar[13])+0.6123724356957944*alphax[1]*fin[12]+(0.5477225575051661*alphax[11]+1.224744871391589*alphavpar[1])*fin[11]+1.369306393762915*alphavpar[3]*fin[10]+0.6123724356957944*alphax[0]*fin[8]+1.224744871391589*fin[4]*alphavpar[7]+1.369306393762915*alphavpar[5]*fin[6]+(0.5477225575051661*alphax[4]+1.369306393762915*alphavpar[0])*fin[4]+(0.5477225575051661*alphax[2]+1.369306393762915*alphavpar[1])*fin[2]; 
  out[13] += 1.224744871391589*(alphax[4]*fin[17]+alphax[1]*fin[13])+fin[10]*(1.224744871391589*alphax[11]+1.369306393762915*alphax[2])+1.224744871391589*fin[5]*alphax[7]+1.369306393762915*(alphax[4]*fin[6]+alphax[0]*fin[5]+alphax[1]*fin[3]); 
  out[14] += 1.224744871391589*alphavpar[5]*fin[19]+1.369306393762915*alphavpar[7]*fin[17]+1.224744871391589*alphavpar[3]*fin[16]+1.369306393762915*(fin[11]*alphavpar[13]+alphavpar[1]*fin[10]+alphavpar[0]*fin[6]+fin[4]*alphavpar[5]+fin[2]*alphavpar[3]); 
  out[15] += 0.6123724356957944*(alphax[4]*fin[19]+alphax[2]*fin[16]+alphax[1]*fin[15]+alphax[0]*fin[9]); 
  out[16] += 0.6123724356957944*alphavpar[1]*fin[15]+0.5477225575051661*alphavpar[13]*fin[13]+0.6123724356957944*alphavpar[0]*fin[9]+0.5477225575051661*(alphavpar[5]*fin[5]+alphavpar[3]*fin[3]); 
  out[17] += 1.095445115010332*alphax[11]*fin[18]+1.224744871391589*(alphax[2]*fin[18]+alphax[1]*fin[17])+0.4898979485566356*alphavpar[5]*fin[15]+1.224744871391589*alphax[4]*fin[14]+(0.3912303982179757*alphavpar[7]+1.224744871391589*alphax[4]+0.6123724356957944*alphavpar[0])*fin[13]+(0.5477225575051661*fin[9]+0.3912303982179757*fin[7]+0.6123724356957944*fin[0])*alphavpar[13]+1.224744871391589*fin[5]*alphax[11]+(1.224744871391589*alphax[7]+1.369306393762915*alphax[0])*fin[10]+0.6123724356957944*(alphavpar[3]*fin[7]+fin[3]*alphavpar[7])+1.369306393762915*(alphax[1]*fin[6]+alphax[2]*fin[5])+0.5477225575051661*(alphavpar[1]*fin[5]+fin[1]*alphavpar[5])+1.369306393762915*fin[3]*alphax[4]; 
  out[18] += (1.095445115010332*alphavpar[13]+1.224744871391589*alphavpar[3])*fin[19]+0.6123724356957944*alphax[1]*fin[18]+0.5477225575051661*alphax[11]*fin[17]+1.224744871391589*(alphavpar[1]*fin[17]+alphavpar[5]*fin[16])+0.6123724356957944*alphax[0]*fin[14]+1.224744871391589*(fin[4]*alphavpar[13]+alphavpar[5]*fin[11])+(1.224744871391589*alphavpar[7]+0.5477225575051661*alphax[4]+1.369306393762915*alphavpar[0])*fin[10]+0.5477225575051661*alphax[2]*fin[6]+1.369306393762915*(alphavpar[1]*fin[6]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]); 
  out[19] += 0.6123724356957944*(alphax[1]*fin[19]+alphax[0]*fin[16])+(0.5477225575051661*alphavpar[7]+0.6123724356957944*(alphax[4]+alphavpar[0]))*fin[15]+0.4898979485566356*(alphavpar[5]*fin[13]+fin[5]*alphavpar[13])+0.6123724356957944*(alphax[2]+alphavpar[1])*fin[9]+0.5477225575051661*(alphavpar[3]*fin[5]+fin[3]*alphavpar[5]); 

  return 0.; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(apardot[2]*fin[7]+apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[4] += -(0.2449489742783178*(4.47213595499958*apardot[1]*fin[7]+4.47213595499958*fin[1]*apardot[2]+5.0*apardot[0]*fin[1]+5.0*fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[6] += -(0.1414213562373095*(8.660254037844387*apardot[2]*fin[13]+8.660254037844386*apardot[1]*fin[5]+8.660254037844386*apardot[0]*fin[3])*q_*rdvpar2)/m_; 
  out[8] += -(0.7071067811865475*(3.872983346207417*apardot[2]*fin[11]+3.872983346207417*apardot[1]*fin[4]+3.872983346207417*apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.1414213562373095*(7.745966692414834*apardot[1]*fin[13]+7.745966692414834*apardot[2]*fin[5]+8.660254037844386*apardot[0]*fin[5]+8.660254037844386*apardot[1]*fin[3])*q_*rdvpar2)/m_; 
  out[11] += -(0.07824607964359516*(10.0*apardot[2]*fin[7]+15.65247584249853*apardot[0]*fin[7]+15.65247584249853*fin[0]*apardot[2]+14.0*apardot[1]*fin[1])*q_*rdvpar2)/m_; 
  out[12] += -(0.1414213562373095*(17.32050807568877*apardot[1]*fin[11]+17.32050807568877*apardot[2]*fin[4]+19.36491673103708*apardot[0]*fin[4]+19.36491673103708*apardot[1]*fin[2])*q_*rdvpar2)/m_; 
  out[14] += -(2.738612787525831*(apardot[2]*fin[17]+apardot[1]*fin[10]+apardot[0]*fin[6])*q_*rdvpar2)/m_; 
  out[16] += -(0.1414213562373095*(8.660254037844386*apardot[1]*fin[15]+8.660254037844387*apardot[0]*fin[9])*q_*rdvpar2)/m_; 
  out[17] += -(0.02020305089104421*(38.72983346207417*apardot[2]*fin[13]+60.62177826491071*apardot[0]*fin[13]+54.22176684690384*apardot[1]*fin[5]+60.6217782649107*apardot[2]*fin[3])*q_*rdvpar2)/m_; 
  out[18] += -(1.224744871391589*(2.0*apardot[1]*fin[17]+2.0*apardot[2]*fin[10]+2.23606797749979*apardot[0]*fin[10]+2.23606797749979*apardot[1]*fin[6])*q_*rdvpar2)/m_; 
  out[19] += -(0.1414213562373095*(7.745966692414834*apardot[2]*fin[15]+8.660254037844387*apardot[0]*fin[15]+8.660254037844386*apardot[1]*fin[9])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.125*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.125*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*apardot[2]-0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7905694150420947*apardot[2])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*apardot[2]+0.9486832980505137*apardot[1]+0.7071067811865475*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
