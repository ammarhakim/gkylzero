#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out) 
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

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[12] = {0.}; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 

  double BstarZdBmag[12] = {0.}; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobtot_inv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobtot_inv[1]*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobtot_inv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 

  double alphax[12] = {0.}; 
  alphax[0] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[2]*hamil[8]+BstarZdBmag[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[4]*hamil[8]+BstarZdBmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[0]*hamil[8]+BstarZdBmag[2]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil[8]+hamil[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphax[8] = (1.224744871391589*BstarZdBmag[2]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphax[9] = (1.224744871391589*BstarZdBmag[4]*hamil[8]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.5477225575051661*alphax[9])+0.3162277660168379*alphax[8]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6846531968814574*alphax[9]-0.3952847075210473*alphax[8]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alphax[9])+0.3162277660168379*alphax[8]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alphax[9])+0.3162277660168379*alphax[8]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6846531968814574*alphax[9]-0.3952847075210473*alphax[8]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alphax[9])+0.3162277660168379*alphax[8]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.5477225575051661*alphax[9]+0.3162277660168379*alphax[8]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6846531968814574*alphax[9])-0.3952847075210473*alphax[8]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alphax[9]+0.3162277660168379*alphax[8]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alphax[9]+0.3162277660168379*alphax[8]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6846531968814574*alphax[9])-0.3952847075210473*alphax[8]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alphax[9]+0.3162277660168379*alphax[8]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  double alphavpar[12] = {0.}; 
  alphavpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(0.6123724356957944*BstarZdBmag[1]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.6123724356957944*hamil[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.6123724356957944*BstarZdBmag[0]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[4] = -(0.6123724356957944*hamil[1]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.6123724356957944*BstarZdBmag[1]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.6123724356957944*BstarZdBmag[2]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[7] = -(0.6123724356957944*BstarZdBmag[4]*hamil[5]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.6123724356957944*alphavpar[7])+0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alphavpar[7]-0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*(alphavpar[7]+alphavpar[6])-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*(alphavpar[7]+alphavpar[6]))+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.6123724356957944*alphavpar[7]-0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alphavpar[7])+0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*(alphavpar[7]+alphavpar[6]))-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*(alphavpar[7]+alphavpar[6])+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.6123724356957944*(alphax[9]*fin[9]+alphax[8]*fin[8]+alphax[4]*fin[4]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.6123724356957944*(alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[4] += 0.07071067811865474*(7.745966692414834*(alphax[4]*fin[9]+fin[4]*alphax[9]+alphax[2]*fin[8]+fin[2]*alphax[8])+8.660254037844386*(alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]+(alphavpar[2]+alphax[1])*fin[4]+fin[1]*alphax[4]+fin[2]*(alphavpar[4]+alphax[0])+fin[0]*alphax[2]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1])); 
  out[5] += 0.07071067811865474*(8.660254037844387*(alphax[9]*fin[11]+alphax[8]*fin[10])+8.660254037844386*(alphax[4]*fin[7]+alphax[2]*fin[6]+alphax[1]*fin[5]+alphax[0]*fin[3])); 
  out[6] += 0.6123724356957944*(alphavpar[4]*fin[7]+fin[4]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[7] += 0.07071067811865474*(7.745966692414834*(alphax[4]*fin[11]+alphax[2]*fin[10]+fin[7]*alphax[9]+fin[6]*alphax[8])+8.660254037844386*((alphavpar[2]+alphax[1])*fin[7]+fin[2]*alphavpar[7]+(alphavpar[4]+alphax[0])*fin[6]+fin[4]*alphavpar[6]+(alphax[4]+alphavpar[0])*fin[5]+fin[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*fin[3]+fin[1]*alphavpar[3])); 
  out[8] += 0.07071067811865474*(17.32050807568877*alphavpar[7]*fin[11]+17.32050807568877*(alphavpar[6]*fin[10]+alphavpar[4]*fin[9])+17.32050807568877*alphavpar[2]*fin[8]+19.36491673103709*(alphavpar[5]*fin[7]+fin[5]*alphavpar[7]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2])); 
  out[9] += 0.01010152544552211*(121.2435565298214*alphavpar[6]*fin[11]+121.2435565298214*alphavpar[7]*fin[10]+(38.72983346207417*alphax[9]+121.2435565298214*alphavpar[2])*fin[9]+60.6217782649107*(alphax[1]*fin[9]+fin[1]*alphax[9])+(38.72983346207417*alphax[8]+121.2435565298214*alphavpar[4])*fin[8]+60.62177826491071*(alphax[0]*fin[8]+fin[0]*alphax[8])+135.5544171172596*(alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6])+54.22176684690384*alphax[4]*fin[4]+135.5544171172596*(alphavpar[0]*fin[4]+fin[0]*alphavpar[4])+54.22176684690384*alphax[2]*fin[2]+135.5544171172596*(alphavpar[1]*fin[2]+fin[1]*alphavpar[2])); 
  out[10] += 0.07071067811865474*(17.32050807568877*alphavpar[4]*fin[11]+17.32050807568877*(alphavpar[2]*fin[10]+alphavpar[7]*fin[9])+17.32050807568877*alphavpar[6]*fin[8]+19.36491673103708*(alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3])); 
  out[11] += 0.01010152544552211*((38.72983346207417*alphax[9]+121.2435565298214*alphavpar[2]+60.6217782649107*alphax[1])*fin[11]+(38.72983346207417*alphax[8]+121.2435565298214*alphavpar[4]+60.62177826491071*alphax[0])*fin[10]+121.2435565298214*alphavpar[6]*fin[9]+60.62177826491071*fin[5]*alphax[9]+121.2435565298214*alphavpar[7]*fin[8]+60.6217782649107*fin[3]*alphax[8]+54.22176684690384*alphax[4]*fin[7]+135.5544171172596*(alphavpar[0]*fin[7]+fin[0]*alphavpar[7])+54.22176684690384*alphax[2]*fin[6]+135.5544171172596*(alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4])); 

  return cflFreq; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_1x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[4] += -(1.224744871391589*(apardot[0]*fin[1]+fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[6] += -(1.224744871391589*(apardot[1]*fin[5]+apardot[0]*fin[3])*q_*rdvpar2)/m_; 
  out[7] += -(1.224744871391589*(apardot[0]*fin[5]+apardot[1]*fin[3])*q_*rdvpar2)/m_; 
  out[8] += -(2.738612787525831*(apardot[1]*fin[4]+apardot[0]*fin[2])*q_*rdvpar2)/m_; 
  out[9] += -(2.738612787525831*(apardot[0]*fin[4]+apardot[1]*fin[2])*q_*rdvpar2)/m_; 
  out[10] += -(2.738612787525831*(apardot[1]*fin[7]+apardot[0]*fin[6])*q_*rdvpar2)/m_; 
  out[11] += -(2.738612787525831*(apardot[0]*fin[7]+apardot[1]*fin[6])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*apardot[0]-0.7071067811865475*apardot[1])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.08838834764831843*(apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
