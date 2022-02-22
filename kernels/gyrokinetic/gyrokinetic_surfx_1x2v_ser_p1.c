#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out) 
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
  // fL,fC,fR: distribution function in left, center and right cells.
  // out: output increment in center cell.

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

  double *b_x = &b_i[0];
  double *b_y = &b_i[2];
  double *b_z = &b_i[4];

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

  double alphaL[4]; 
  alphaL[0] = -(0.25*(3.0*BstarZdBmag[1]-1.732050807568877*BstarZdBmag[0])*hamil[2]*rdvpar2)/m_; 
  alphaL[1] = -(0.25*hamil[2]*(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double alphaR[4]; 
  alphaR[0] = (0.25*(3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0])*hamil[2]*rdvpar2)/m_; 
  alphaR[1] = (0.25*hamil[2]*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double fUpOrdL[4];
  if (alphaL[0]-1.0*alphaL[1] > 0.) {
    fUpOrdL[0] = 0.6123724356957942*fL[7]+0.3535533905932737*fL[6]-0.6123724356957944*(fL[5]+fL[4])-0.3535533905932736*(fL[3]+fL[2])+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[0] = (-0.6123724356957942*fC[7])+0.3535533905932737*fC[6]+0.6123724356957944*(fC[5]+fC[4])-0.3535533905932736*(fC[3]+fC[2])-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[0]-1.0*alphaL[1] > 0.) {
    fUpOrdL[1] = (-0.6123724356957942*fL[7])-0.3535533905932737*fL[6]+0.6123724356957944*fL[5]-0.6123724356957944*fL[4]+0.3535533905932736*fL[3]-0.3535533905932736*fL[2]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[1] = 0.6123724356957942*fC[7]-0.3535533905932737*fC[6]-0.6123724356957944*fC[5]+0.6123724356957944*fC[4]+0.3535533905932736*fC[3]-0.3535533905932736*fC[2]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[2] = (-0.6123724356957942*fL[7])-0.3535533905932737*fL[6]-0.6123724356957944*fL[5]+0.6123724356957944*fL[4]-0.3535533905932736*fL[3]+0.3535533905932736*fL[2]+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[2] = 0.6123724356957942*fC[7]-0.3535533905932737*fC[6]+0.6123724356957944*fC[5]-0.6123724356957944*fC[4]-0.3535533905932736*fC[3]+0.3535533905932736*fC[2]-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }
  if (alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[3] = 0.6123724356957942*fL[7]+0.3535533905932737*fL[6]+0.6123724356957944*(fL[5]+fL[4])+0.3535533905932736*(fL[3]+fL[2])+0.6123724356957944*fL[1]+0.3535533905932737*fL[0]; 
  } else {
    fUpOrdL[3] = (-0.6123724356957942*fC[7])+0.3535533905932737*fC[6]-0.6123724356957944*(fC[5]+fC[4])+0.3535533905932736*(fC[3]+fC[2])-0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  }

  double fUpL[4] = {0.};
  fUpL[0] = 0.5*fUpOrdL[3]+0.5*fUpOrdL[2]+0.5*fUpOrdL[1]+0.5*fUpOrdL[0]; 
  fUpL[1] = 0.5000000000000001*fUpOrdL[3]+0.5000000000000001*fUpOrdL[2]-0.5000000000000001*fUpOrdL[1]-0.5000000000000001*fUpOrdL[0]; 
  fUpL[2] = 0.5000000000000001*fUpOrdL[3]-0.5000000000000001*fUpOrdL[2]+0.5000000000000001*fUpOrdL[1]-0.5000000000000001*fUpOrdL[0]; 
  fUpL[3] = 0.5*fUpOrdL[3]-0.5*fUpOrdL[2]-0.5*fUpOrdL[1]+0.5*fUpOrdL[0]; 

  double GhatL[8] = {0.}; 
  GhatL[0] += 0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] += 0.5*alphaL[1]*fUpL[3]+0.5*alphaL[0]*fUpL[2]; 
  GhatL[3] += 0.5*alphaL[0]*fUpL[3]+0.5*alphaL[1]*fUpL[2]; 

  double fUpOrdR[4];
  if (alphaR[0]-1.0*alphaR[1] > 0.) {
    fUpOrdR[0] = 0.6123724356957942*fC[7]+0.3535533905932737*fC[6]-0.6123724356957944*(fC[5]+fC[4])-0.3535533905932736*(fC[3]+fC[2])+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[0] = (-0.6123724356957942*fR[7])+0.3535533905932737*fR[6]+0.6123724356957944*(fR[5]+fR[4])-0.3535533905932736*(fR[3]+fR[2])-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[0]-1.0*alphaR[1] > 0.) {
    fUpOrdR[1] = (-0.6123724356957942*fC[7])-0.3535533905932737*fC[6]+0.6123724356957944*fC[5]-0.6123724356957944*fC[4]+0.3535533905932736*fC[3]-0.3535533905932736*fC[2]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[1] = 0.6123724356957942*fR[7]-0.3535533905932737*fR[6]-0.6123724356957944*fR[5]+0.6123724356957944*fR[4]+0.3535533905932736*fR[3]-0.3535533905932736*fR[2]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[2] = (-0.6123724356957942*fC[7])-0.3535533905932737*fC[6]-0.6123724356957944*fC[5]+0.6123724356957944*fC[4]-0.3535533905932736*fC[3]+0.3535533905932736*fC[2]+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[2] = 0.6123724356957942*fR[7]-0.3535533905932737*fR[6]+0.6123724356957944*fR[5]-0.6123724356957944*fR[4]-0.3535533905932736*fR[3]+0.3535533905932736*fR[2]-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }
  if (alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[3] = 0.6123724356957942*fC[7]+0.3535533905932737*fC[6]+0.6123724356957944*(fC[5]+fC[4])+0.3535533905932736*(fC[3]+fC[2])+0.6123724356957944*fC[1]+0.3535533905932737*fC[0]; 
  } else {
    fUpOrdR[3] = (-0.6123724356957942*fR[7])+0.3535533905932737*fR[6]-0.6123724356957944*(fR[5]+fR[4])+0.3535533905932736*(fR[3]+fR[2])-0.6123724356957944*fR[1]+0.3535533905932737*fR[0]; 
  }

  double fUpR[4] = {0.};
  fUpR[0] = 0.5*fUpOrdR[3]+0.5*fUpOrdR[2]+0.5*fUpOrdR[1]+0.5*fUpOrdR[0]; 
  fUpR[1] = 0.5000000000000001*fUpOrdR[3]+0.5000000000000001*fUpOrdR[2]-0.5000000000000001*fUpOrdR[1]-0.5000000000000001*fUpOrdR[0]; 
  fUpR[2] = 0.5000000000000001*fUpOrdR[3]-0.5000000000000001*fUpOrdR[2]+0.5000000000000001*fUpOrdR[1]-0.5000000000000001*fUpOrdR[0]; 
  fUpR[3] = 0.5*fUpOrdR[3]-0.5*fUpOrdR[2]-0.5*fUpOrdR[1]+0.5*fUpOrdR[0]; 

  double GhatR[8] = {0.}; 
  GhatR[0] += 0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] += 0.5*alphaR[1]*fUpR[3]+0.5*alphaR[0]*fUpR[2]; 
  GhatR[3] += 0.5*alphaR[0]*fUpR[3]+0.5*alphaR[1]*fUpR[2]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[5] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[6] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[7] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 

} 
