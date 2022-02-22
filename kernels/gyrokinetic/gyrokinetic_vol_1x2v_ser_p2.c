#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out) 
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
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[20]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double BstarZdBmag[20]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmag[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[20]; 
  alphax[0] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[2]*hamil[8]+BstarZdBmag[0]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[4]*hamil[8]+BstarZdBmag[1]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[0]*hamil[8]+BstarZdBmag[2]*hamil[2])*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil[8]+hamil[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alphax[7] = (0.3535533905932737*(3.872983346207417*hamil[8]*BstarZdBmag[11]+1.732050807568877*hamil[2]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alphax[8] = (1.224744871391589*BstarZdBmag[2]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphax[11] = (0.3535533905932737*(1.732050807568877*hamil[2]*BstarZdBmag[11]+3.872983346207417*BstarZdBmag[7]*hamil[8])*rdvpar2*rdx2)/m_; 
  alphax[12] = (1.224744871391589*BstarZdBmag[4]*hamil[8]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6846531968814574*alphax[12]-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6846531968814574*alphax[12]-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.6846531968814574*alphax[12]-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphax[12])-1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]+0.821583836257749*alphax[4]-0.4743416490252567*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphax[12])-1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]+0.821583836257749*alphax[4]-0.4743416490252567*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphax[12])-1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]+0.821583836257749*alphax[4]-0.4743416490252567*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphax[12])+1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]-0.821583836257749*alphax[4]+0.4743416490252567*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphax[12])+1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]-0.821583836257749*alphax[4]+0.4743416490252567*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphax[12])+1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]-0.821583836257749*alphax[4]+0.4743416490252567*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6846531968814574*alphax[12])-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6846531968814574*alphax[12])-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.6846531968814574*alphax[12])-0.3952847075210473*alphax[8]+0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphax[12]-1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]-0.821583836257749*alphax[4]-0.4743416490252567*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphax[12]-1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]-0.821583836257749*alphax[4]-0.4743416490252567*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphax[12]-1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]-0.821583836257749*alphax[4]-0.4743416490252567*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphax[12]+1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]+0.821583836257749*alphax[4]+0.4743416490252567*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphax[12]+1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]+0.821583836257749*alphax[4]+0.4743416490252567*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphax[12]+1.060660171779821*alphax[11]+0.3162277660168378*alphax[8]+0.7905694150420947*alphax[7]+0.821583836257749*alphax[4]+0.4743416490252567*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  double alphavpar[20]; 
  alphavpar[0] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(0.6123724356957944*((2.0*BstarZdBmag[7]+2.23606797749979*BstarZdBmag[0])*hamil[7]+BstarZdBmag[1]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[2] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[4]*hamil[7]+hamil[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.3535533905932737*(3.872983346207417*BstarZdBmag[1]*hamil[13]+1.732050807568877*BstarZdBmag[0]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[4] = -(0.3535533905932737*(3.464101615137755*hamil[7]*BstarZdBmag[11]+1.732050807568877*(2.23606797749979*BstarZdBmag[2]*hamil[7]+hamil[1]*BstarZdBmag[4]))*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.3535533905932737*(3.872983346207417*(0.8944271909999159*BstarZdBmag[7]+BstarZdBmag[0])*hamil[13]+1.732050807568877*BstarZdBmag[1]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.3535533905932737*(3.872983346207417*BstarZdBmag[4]*hamil[13]+1.732050807568877*BstarZdBmag[2]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[7] = -(0.6123724356957944*(2.0*BstarZdBmag[1]*hamil[7]+hamil[1]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alphavpar[10] = -(0.3535533905932737*((3.464101615137754*BstarZdBmag[11]+3.872983346207417*BstarZdBmag[2])*hamil[13]+1.732050807568877*BstarZdBmag[4]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[11] = -(0.3535533905932737*(1.732050807568877*hamil[1]*BstarZdBmag[11]+3.464101615137755*BstarZdBmag[4]*hamil[7])*rdvpar2*rdx2)/m_; 
  alphavpar[13] = -(0.3535533905932737*(3.464101615137754*BstarZdBmag[1]*hamil[13]+1.732050807568877*hamil[5]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alphavpar[17] = -(0.6123724356957944*(2.0*BstarZdBmag[4]*hamil[13]+hamil[5]*BstarZdBmag[11])*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6846531968814574*alphavpar[11]-0.3952847075210473*alphavpar[7]-0.6123724356957944*alphavpar[2]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.9185586535436915*alphavpar[17])+0.5303300858899105*alphavpar[13]+0.6846531968814574*alphavpar[11]-0.3952847075210473*alphavpar[7]+0.821583836257749*alphavpar[6]-0.4743416490252567*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.9185586535436915*alphavpar[17]-0.5303300858899105*alphavpar[13]+0.6846531968814574*alphavpar[11]-0.3952847075210473*alphavpar[7]-0.821583836257749*alphavpar[6]+0.4743416490252567*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphavpar[11])+0.3162277660168378*alphavpar[7]+0.821583836257749*alphavpar[4]-0.6123724356957944*alphavpar[2]-0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.7348469228349532*alphavpar[17]-0.4242640687119284*alphavpar[13]-0.547722557505166*alphavpar[11]-1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]+0.821583836257749*alphavpar[6]+0.6363961030678926*alphavpar[5]+0.821583836257749*alphavpar[4]-0.4743416490252567*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.7348469228349532*alphavpar[17])+0.4242640687119284*alphavpar[13]-0.547722557505166*alphavpar[11]+1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]-0.821583836257749*alphavpar[6]-0.6363961030678926*alphavpar[5]+0.821583836257749*alphavpar[4]+0.4743416490252567*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.547722557505166*alphavpar[11])+0.3162277660168378*alphavpar[7]-0.821583836257749*alphavpar[4]-0.6123724356957944*alphavpar[2]+0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*(0.7348469228349532*alphavpar[17]-0.4242640687119284*alphavpar[13]-0.547722557505166*alphavpar[11]+1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]+0.821583836257749*alphavpar[6]-0.6363961030678926*alphavpar[5]-0.821583836257749*alphavpar[4]-0.4743416490252567*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.125*((-0.7348469228349532*alphavpar[17])+0.4242640687119284*alphavpar[13]-0.547722557505166*alphavpar[11]-1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]-0.821583836257749*alphavpar[6]+0.6363961030678926*alphavpar[5]-0.821583836257749*alphavpar[4]+0.4743416490252567*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6846531968814574*alphavpar[11])-0.3952847075210473*alphavpar[7]+0.6123724356957944*alphavpar[2]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.9185586535436915*alphavpar[17]+0.5303300858899105*alphavpar[13]-0.6846531968814574*alphavpar[11]-0.3952847075210473*alphavpar[7]-0.821583836257749*alphavpar[6]-0.4743416490252567*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.9185586535436915*alphavpar[17])-0.5303300858899105*alphavpar[13]-0.6846531968814574*alphavpar[11]-0.3952847075210473*alphavpar[7]+0.821583836257749*alphavpar[6]+0.4743416490252567*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphavpar[11]+0.3162277660168378*alphavpar[7]-0.821583836257749*alphavpar[4]+0.6123724356957944*alphavpar[2]-0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.7348469228349532*alphavpar[17])-0.4242640687119284*alphavpar[13]+0.547722557505166*alphavpar[11]+1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]-0.821583836257749*alphavpar[6]+0.6363961030678926*alphavpar[5]-0.821583836257749*alphavpar[4]-0.4743416490252567*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.7348469228349532*alphavpar[17]+0.4242640687119284*alphavpar[13]+0.547722557505166*alphavpar[11]-1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]+0.821583836257749*alphavpar[6]-0.6363961030678926*alphavpar[5]-0.821583836257749*alphavpar[4]+0.4743416490252567*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.547722557505166*alphavpar[11]+0.3162277660168378*alphavpar[7]+0.821583836257749*alphavpar[4]+0.6123724356957944*alphavpar[2]+0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*((-0.7348469228349532*alphavpar[17])-0.4242640687119284*alphavpar[13]+0.547722557505166*alphavpar[11]-1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]-0.821583836257749*alphavpar[6]-0.6363961030678926*alphavpar[5]+0.821583836257749*alphavpar[4]-0.4743416490252567*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.125*(0.7348469228349532*alphavpar[17]+0.4242640687119284*alphavpar[13]+0.547722557505166*alphavpar[11]+1.10227038425243*alphavpar[10]+0.3162277660168378*alphavpar[7]+0.821583836257749*alphavpar[6]+0.6363961030678926*alphavpar[5]+0.821583836257749*alphavpar[4]+0.4743416490252567*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.4743416490252567*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.6123724356957944*(alphax[12]*f[12]+alphax[11]*f[11]+alphax[8]*f[8]+alphax[7]*f[7]+alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphavpar[17]*f[17]+alphavpar[13]*f[13]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[4] += 0.07071067811865474*(7.745966692414834*(alphavpar[10]*f[17]+f[10]*alphavpar[17]+alphavpar[5]*f[13]+f[5]*alphavpar[13]+alphax[4]*f[12]+f[4]*alphax[12])+(8.660254037844387*alphax[7]+7.745966692414834*alphavpar[4])*f[11]+8.660254037844387*f[7]*alphax[11]+7.745966692414834*f[4]*alphavpar[11]+8.660254037844386*(alphavpar[6]*f[10]+f[6]*alphavpar[10])+7.745966692414834*(alphax[2]*f[8]+f[2]*alphax[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7])+8.660254037844386*(alphavpar[3]*f[5]+f[3]*alphavpar[5]+(alphavpar[2]+alphax[1])*f[4]+f[1]*alphax[4]+f[2]*(alphavpar[4]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[5] += 0.07071067811865474*(8.660254037844387*(alphax[12]*f[18]+alphax[11]*f[17]+alphax[8]*f[14]+alphax[7]*f[13])+8.660254037844386*(alphax[4]*f[10]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3])); 
  out[6] += 0.07071067811865474*(7.745966692414834*alphavpar[10]*f[19]+8.660254037844387*(alphavpar[11]*f[17]+f[11]*alphavpar[17])+7.745966692414834*(alphavpar[6]*f[16]+alphavpar[5]*f[15])+8.660254037844387*(alphavpar[7]*f[13]+f[7]*alphavpar[13])+8.660254037844386*(alphavpar[4]*f[10]+f[4]*alphavpar[10])+7.745966692414834*alphavpar[3]*f[9]+8.660254037844386*(alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphavpar[0]*f[3]+f[0]*alphavpar[3])); 
  out[7] += 0.07071067811865474*(19.36491673103708*(alphax[8]*f[12]+f[8]*alphax[12])+17.32050807568877*(alphax[4]*f[11]+f[4]*alphax[11])+17.32050807568877*(alphax[1]*f[7]+f[1]*alphax[7])+19.36491673103709*(alphax[2]*f[4]+f[2]*alphax[4]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[8] += 0.07071067811865474*(17.32050807568877*alphavpar[10]*f[18]+19.36491673103708*(alphavpar[13]*f[17]+f[13]*alphavpar[17])+17.32050807568877*(alphavpar[6]*f[14]+alphavpar[4]*f[12])+19.36491673103708*(alphavpar[7]*f[11]+f[7]*alphavpar[11])+19.36491673103709*(alphavpar[5]*f[10]+f[5]*alphavpar[10])+17.32050807568877*alphavpar[2]*f[8]+19.36491673103709*(alphavpar[3]*f[6]+f[3]*alphavpar[6]+alphavpar[1]*f[4]+f[1]*alphavpar[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2])); 
  out[10] += 0.07071067811865474*(6.928203230275509*alphavpar[17]*f[19]+7.745966692414834*(alphavpar[6]*f[19]+alphax[4]*f[18])+8.660254037844386*alphax[7]*f[17]+7.745966692414834*(alphavpar[4]*f[17]+f[4]*alphavpar[17]+alphavpar[10]*f[16])+6.928203230275509*alphavpar[13]*f[15]+7.745966692414834*(alphavpar[3]*f[15]+alphax[2]*f[14])+8.660254037844386*alphax[11]*f[13]+7.745966692414834*(alphavpar[1]*f[13]+f[1]*alphavpar[13]+f[10]*alphax[12]+alphavpar[10]*f[11]+f[10]*alphavpar[11])+8.660254037844386*((alphavpar[2]+alphax[1])*f[10]+f[2]*alphavpar[10])+7.745966692414834*(alphavpar[5]*f[9]+f[6]*alphax[8]+alphavpar[5]*f[7]+f[5]*alphavpar[7])+8.660254037844386*((alphavpar[4]+alphax[0])*f[6]+f[4]*alphavpar[6]+(alphax[4]+alphavpar[0])*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*f[3]+f[1]*alphavpar[3])); 
  out[11] += 0.01010152544552211*(38.72983346207417*alphavpar[17]*f[17]+60.62177826491071*(alphavpar[6]*f[17]+f[6]*alphavpar[17])+38.72983346207417*alphavpar[13]*f[13]+60.6217782649107*(alphavpar[3]*f[13]+f[3]*alphavpar[13])+(108.4435336938077*alphax[11]+121.2435565298214*alphax[2])*f[12]+(108.4435336938077*f[11]+121.2435565298214*f[2])*alphax[12]+(38.72983346207417*alphavpar[11]+60.6217782649107*alphavpar[2])*f[11]+121.2435565298214*(alphax[1]*f[11]+f[1]*alphax[11])+60.6217782649107*f[2]*alphavpar[11]+54.22176684690384*alphavpar[10]*f[10]+121.2435565298214*(alphax[4]*f[8]+f[4]*alphax[8])+(38.72983346207417*alphavpar[7]+121.2435565298214*alphax[4]+60.62177826491071*alphavpar[0])*f[7]+121.2435565298214*f[4]*alphax[7]+60.62177826491071*f[0]*alphavpar[7]+54.22176684690384*(alphavpar[5]*f[5]+alphavpar[4]*f[4])+135.5544171172596*(alphax[0]*f[4]+f[0]*alphax[4]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphavpar[1])); 
  out[12] += 0.01010152544552211*(108.4435336938077*alphavpar[17]*f[18]+121.2435565298214*(alphavpar[6]*f[18]+alphavpar[5]*f[17]+f[5]*alphavpar[17])+121.2435565298214*(alphavpar[10]*(f[14]+f[13])+f[10]*alphavpar[13])+(38.72983346207417*alphax[12]+108.4435336938077*alphavpar[11]+121.2435565298214*alphavpar[2])*f[12]+60.6217782649107*(alphax[1]*f[12]+f[1]*alphax[12])+54.22176684690384*alphax[11]*f[11]+121.2435565298214*(alphavpar[1]*f[11]+f[1]*alphavpar[11])+135.5544171172596*(alphavpar[3]*f[10]+f[3]*alphavpar[10])+(38.72983346207417*alphax[8]+121.2435565298214*alphavpar[4])*f[8]+60.62177826491071*(alphax[0]*f[8]+f[0]*alphax[8])+121.2435565298214*(alphavpar[4]*f[7]+f[4]*alphavpar[7])+135.5544171172596*(alphavpar[5]*f[6]+f[5]*alphavpar[6])+54.22176684690384*alphax[4]*f[4]+135.5544171172596*(alphavpar[0]*f[4]+f[0]*alphavpar[4])+54.22176684690384*alphax[2]*f[2]+135.5544171172596*(alphavpar[1]*f[2]+f[1]*alphavpar[2])); 
  out[13] += 0.07071067811865474*(19.36491673103708*alphax[8]*f[18]+17.32050807568877*alphax[4]*f[17]+19.36491673103708*alphax[12]*f[14]+17.32050807568877*alphax[1]*f[13]+f[10]*(17.32050807568877*alphax[11]+19.36491673103708*alphax[2])+17.32050807568877*f[5]*alphax[7]+19.36491673103708*(alphax[4]*f[6]+alphax[0]*f[5]+alphax[1]*f[3])); 
  out[14] += 0.07071067811865474*(17.32050807568877*(alphavpar[5]*f[19]+alphavpar[4]*f[18])+19.36491673103708*(alphavpar[7]*f[17]+f[7]*alphavpar[17])+17.32050807568877*(alphavpar[3]*f[16]+alphavpar[10]*f[15]+alphavpar[2]*f[14])+19.36491673103708*(alphavpar[11]*f[13]+f[11]*alphavpar[13])+17.32050807568877*alphavpar[10]*f[12]+19.36491673103708*(alphavpar[1]*f[10]+f[1]*alphavpar[10])+17.32050807568877*alphavpar[6]*(f[9]+f[8])+19.36491673103708*(alphavpar[0]*f[6]+f[0]*alphavpar[6]+alphavpar[4]*f[5]+f[4]*alphavpar[5]+alphavpar[2]*f[3]+f[2]*alphavpar[3])); 
  out[15] += 0.07071067811865474*(8.660254037844387*alphax[4]*f[19]+8.660254037844386*(alphax[2]*f[16]+alphax[1]*f[15])+8.660254037844387*alphax[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844387*alphavpar[4]*f[19]+7.745966692414834*alphavpar[17]*f[17]+8.660254037844386*(alphavpar[2]*f[16]+alphavpar[1]*f[15])+7.745966692414834*(alphavpar[13]*f[13]+alphavpar[10]*f[10])+8.660254037844387*alphavpar[0]*f[9]+7.745966692414834*(alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[3]*f[3])); 
  out[17] += 0.002020305089104421*(242.4871130596428*alphavpar[10]*f[19]+(542.2176684690384*alphax[11]+606.217782649107*alphax[2])*f[18]+(542.2176684690384*alphax[12]+193.6491673103708*alphavpar[11]+303.1088913245535*alphavpar[2]+606.217782649107*alphax[1])*f[17]+(271.1088342345192*f[16]+193.6491673103708*f[11]+303.1088913245535*f[2])*alphavpar[17]+242.4871130596428*alphavpar[5]*f[15]+606.2177826491072*alphax[4]*f[14]+(193.6491673103708*alphavpar[7]+606.2177826491072*alphax[4]+303.1088913245536*alphavpar[0])*f[13]+(271.1088342345192*f[9]+193.6491673103708*f[7]+303.1088913245536*f[0])*alphavpar[13]+606.2177826491072*f[6]*alphax[12]+303.1088913245536*alphavpar[6]*f[11]+606.2177826491072*f[5]*alphax[11]+303.1088913245536*f[6]*alphavpar[11]+(606.217782649107*(alphax[8]+alphax[7])+271.1088342345192*alphavpar[4]+677.7720855862981*alphax[0])*f[10]+271.1088342345192*f[4]*alphavpar[10]+303.1088913245535*(alphavpar[3]*f[7]+f[3]*alphavpar[7])+677.7720855862981*(alphax[1]*f[6]+alphax[2]*f[5])+271.1088342345192*(alphavpar[1]*f[5]+f[1]*alphavpar[5])+677.7720855862981*f[3]*alphax[4]); 
  out[18] += 0.01010152544552211*((108.4435336938077*alphavpar[13]+121.2435565298214*alphavpar[3])*f[19]+(38.72983346207417*alphax[12]+108.4435336938077*alphavpar[11]+121.2435565298214*alphavpar[2]+60.6217782649107*alphax[1])*f[18]+(54.22176684690384*alphax[11]+121.2435565298214*alphavpar[1])*f[17]+(108.4435336938077*(f[15]+f[12])+121.2435565298214*f[1])*alphavpar[17]+121.2435565298214*(alphavpar[5]*f[16]+alphavpar[6]*f[15])+(38.72983346207417*alphax[8]+121.2435565298214*alphavpar[4]+60.62177826491071*alphax[0])*f[14]+121.2435565298214*(alphavpar[4]*f[13]+f[4]*alphavpar[13]+alphavpar[6]*f[12])+60.62177826491071*f[5]*alphax[12]+121.2435565298214*(alphavpar[5]*f[11]+f[5]*alphavpar[11])+(121.2435565298214*alphavpar[7]+54.22176684690384*alphax[4]+135.5544171172596*alphavpar[0])*f[10]+(121.2435565298214*(f[9]+f[8]+f[7])+135.5544171172596*f[0])*alphavpar[10]+60.6217782649107*f[3]*alphax[8]+54.22176684690384*alphax[2]*f[6]+135.5544171172596*(alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphavpar[3]*f[4]+f[3]*alphavpar[4])); 
  out[19] += 0.01414213562373095*((38.72983346207417*(alphax[12]+alphavpar[11])+43.30127018922193*(alphavpar[2]+alphax[1]))*f[19]+34.64101615137754*(alphavpar[10]*f[17]+f[10]*alphavpar[17])+(38.72983346207417*alphax[8]+43.30127018922195*(alphavpar[4]+alphax[0]))*f[16]+(38.72983346207417*alphavpar[7]+43.30127018922195*(alphax[4]+alphavpar[0]))*f[15]+34.64101615137755*(alphavpar[5]*f[13]+f[5]*alphavpar[13])+38.72983346207418*(alphavpar[6]*f[10]+f[6]*alphavpar[10])+43.30127018922193*(alphax[2]+alphavpar[1])*f[9]+38.72983346207418*(alphavpar[3]*f[5]+f[3]*alphavpar[5])); 

  return cflFreq; 
} 
GKYL_CU_DH double gyrokinetic_vol_step2_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[7]+dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[4] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[7]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[6] += -(0.1414213562373095*(8.660254037844387*dApardt[2]*f[13]+8.660254037844386*dApardt[1]*f[5]+8.660254037844386*dApardt[0]*f[3])*q_*rdvpar2)/m_; 
  out[8] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[11]+3.872983346207417*dApardt[1]*f[4]+3.872983346207417*dApardt[0]*f[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.1414213562373095*(7.745966692414834*dApardt[1]*f[13]+7.745966692414834*dApardt[2]*f[5]+8.660254037844386*dApardt[0]*f[5]+8.660254037844386*dApardt[1]*f[3])*q_*rdvpar2)/m_; 
  out[11] += -(0.07824607964359516*(10.0*dApardt[2]*f[7]+15.65247584249853*dApardt[0]*f[7]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*q_*rdvpar2)/m_; 
  out[12] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[11]+17.32050807568877*dApardt[2]*f[4]+19.36491673103708*dApardt[0]*f[4]+19.36491673103708*dApardt[1]*f[2])*q_*rdvpar2)/m_; 
  out[14] += -(2.738612787525831*(dApardt[2]*f[17]+dApardt[1]*f[10]+dApardt[0]*f[6])*q_*rdvpar2)/m_; 
  out[16] += -(0.1414213562373095*(8.660254037844386*dApardt[1]*f[15]+8.660254037844387*dApardt[0]*f[9])*q_*rdvpar2)/m_; 
  out[17] += -(0.02020305089104421*(38.72983346207417*dApardt[2]*f[13]+60.62177826491071*dApardt[0]*f[13]+54.22176684690384*dApardt[1]*f[5]+60.6217782649107*dApardt[2]*f[3])*q_*rdvpar2)/m_; 
  out[18] += -(1.224744871391589*(2.0*dApardt[1]*f[17]+2.0*dApardt[2]*f[10]+2.23606797749979*dApardt[0]*f[10]+2.23606797749979*dApardt[1]*f[6])*q_*rdvpar2)/m_; 
  out[19] += -(0.1414213562373095*(7.745966692414834*dApardt[2]*f[15]+8.660254037844387*dApardt[0]*f[15]+8.660254037844386*dApardt[1]*f[9])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336758*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336758*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336758*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336758*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336758*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336758*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336758*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336758*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336758*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336758*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336758*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336758*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
