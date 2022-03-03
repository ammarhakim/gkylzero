#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
  // fl,fc,fr: distribution function in left, center and right cells.
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
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobtot_inv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobtot_inv[1]*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobtot_inv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[4]; 
  alphaL[0] = (0.25*(hamil[1]*(3.0*BstarZdBmag[2]-1.732050807568877*BstarZdBmag[0])*rdx2-5.656854249492382*apardot[0]*q_))/m_; 
  alphaL[1] = (0.25*(hamil[1]*(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[1])*rdx2-5.656854249492382*apardot[1]*q_))/m_; 
  alphaL[2] = (0.25*(3.0*BstarZdBmag[2]-1.732050807568877*BstarZdBmag[0])*hamil[5]*rdx2)/m_; 
  alphaL[3] = (0.25*(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[1])*hamil[5]*rdx2)/m_; 

  double alphaR[4]; 
  alphaR[0] = -(0.25*(hamil[1]*(3.0*BstarZdBmag[2]+1.732050807568877*BstarZdBmag[0])*rdx2+5.656854249492382*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.25*(hamil[1]*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[1])*rdx2+5.656854249492382*apardot[1]*q_))/m_; 
  alphaR[2] = -(0.25*(3.0*BstarZdBmag[2]+1.732050807568877*BstarZdBmag[0])*hamil[5]*rdx2)/m_; 
  alphaR[3] = -(0.25*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[1])*hamil[5]*rdx2)/m_; 

  double fUpOrdL[4];
  if (alphaL[3]-1.0*alphaL[2]-1.0*alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[0] = 0.6123724356957942*fl[7]-0.6123724356957944*fl[6]+0.3535533905932737*fl[5]-0.6123724356957944*fl[4]-0.3535533905932736*fl[3]+0.6123724356957944*fl[2]-0.3535533905932736*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[0] = (-0.6123724356957942*fc[7])+0.6123724356957944*fc[6]+0.3535533905932737*fc[5]+0.6123724356957944*fc[4]-0.3535533905932736*fc[3]-0.6123724356957944*fc[2]-0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  }
  if ((-1.0*alphaL[3])+alphaL[2]-1.0*alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[1] = (-0.6123724356957942*fl[7])+0.6123724356957944*fl[6]-0.3535533905932737*fl[5]-0.6123724356957944*fl[4]+0.3535533905932736*fl[3]+0.6123724356957944*fl[2]-0.3535533905932736*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[1] = 0.6123724356957942*fc[7]-0.6123724356957944*fc[6]-0.3535533905932737*fc[5]+0.6123724356957944*fc[4]+0.3535533905932736*fc[3]-0.6123724356957944*fc[2]-0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  }
  if ((-1.0*alphaL[3])-1.0*alphaL[2]+alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[2] = (-0.6123724356957942*fl[7])-0.6123724356957944*fl[6]-0.3535533905932737*fl[5]+0.6123724356957944*fl[4]-0.3535533905932736*fl[3]+0.6123724356957944*fl[2]+0.3535533905932736*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[2] = 0.6123724356957942*fc[7]+0.6123724356957944*fc[6]-0.3535533905932737*fc[5]-0.6123724356957944*fc[4]-0.3535533905932736*fc[3]-0.6123724356957944*fc[2]+0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  }
  if (alphaL[3]+alphaL[2]+alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[3] = 0.6123724356957942*fl[7]+0.6123724356957944*fl[6]+0.3535533905932737*fl[5]+0.6123724356957944*fl[4]+0.3535533905932736*fl[3]+0.6123724356957944*fl[2]+0.3535533905932736*fl[1]+0.3535533905932737*fl[0]; 
  } else {
    fUpOrdL[3] = (-0.6123724356957942*fc[7])-0.6123724356957944*fc[6]+0.3535533905932737*fc[5]-0.6123724356957944*fc[4]+0.3535533905932736*fc[3]-0.6123724356957944*fc[2]+0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  }

  double fUpL[4] = {0.};
  fUpL[0] = 0.5*fUpOrdL[3]+0.5*fUpOrdL[2]+0.5*fUpOrdL[1]+0.5*fUpOrdL[0]; 
  fUpL[1] = 0.5000000000000001*fUpOrdL[3]+0.5000000000000001*fUpOrdL[2]-0.5000000000000001*fUpOrdL[1]-0.5000000000000001*fUpOrdL[0]; 
  fUpL[2] = 0.5000000000000001*fUpOrdL[3]-0.5000000000000001*fUpOrdL[2]+0.5000000000000001*fUpOrdL[1]-0.5000000000000001*fUpOrdL[0]; 
  fUpL[3] = 0.5*fUpOrdL[3]-0.5*fUpOrdL[2]-0.5*fUpOrdL[1]+0.5*fUpOrdL[0]; 

  double GhatL[8] = {0.}; 
  GhatL[0] += 0.5*alphaL[3]*fUpL[3]+0.5*alphaL[2]*fUpL[2]+0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.5*alphaL[2]*fUpL[3]+0.5*fUpL[2]*alphaL[3]+0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] += 0.5*alphaL[1]*fUpL[3]+0.5*fUpL[1]*alphaL[3]+0.5*alphaL[0]*fUpL[2]+0.5*fUpL[0]*alphaL[2]; 
  GhatL[3] += 0.5*alphaL[0]*fUpL[3]+0.5*fUpL[0]*alphaL[3]+0.5*alphaL[1]*fUpL[2]+0.5*fUpL[1]*alphaL[2]; 

  double fUpOrdR[4];
  if (alphaR[3]-1.0*alphaR[2]-1.0*alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[0] = 0.6123724356957942*fc[7]-0.6123724356957944*fc[6]+0.3535533905932737*fc[5]-0.6123724356957944*fc[4]-0.3535533905932736*fc[3]+0.6123724356957944*fc[2]-0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[0] = (-0.6123724356957942*fr[7])+0.6123724356957944*fr[6]+0.3535533905932737*fr[5]+0.6123724356957944*fr[4]-0.3535533905932736*fr[3]-0.6123724356957944*fr[2]-0.3535533905932736*fr[1]+0.3535533905932737*fr[0]; 
  }
  if ((-1.0*alphaR[3])+alphaR[2]-1.0*alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[1] = (-0.6123724356957942*fc[7])+0.6123724356957944*fc[6]-0.3535533905932737*fc[5]-0.6123724356957944*fc[4]+0.3535533905932736*fc[3]+0.6123724356957944*fc[2]-0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[1] = 0.6123724356957942*fr[7]-0.6123724356957944*fr[6]-0.3535533905932737*fr[5]+0.6123724356957944*fr[4]+0.3535533905932736*fr[3]-0.6123724356957944*fr[2]-0.3535533905932736*fr[1]+0.3535533905932737*fr[0]; 
  }
  if ((-1.0*alphaR[3])-1.0*alphaR[2]+alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[2] = (-0.6123724356957942*fc[7])-0.6123724356957944*fc[6]-0.3535533905932737*fc[5]+0.6123724356957944*fc[4]-0.3535533905932736*fc[3]+0.6123724356957944*fc[2]+0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[2] = 0.6123724356957942*fr[7]+0.6123724356957944*fr[6]-0.3535533905932737*fr[5]-0.6123724356957944*fr[4]-0.3535533905932736*fr[3]-0.6123724356957944*fr[2]+0.3535533905932736*fr[1]+0.3535533905932737*fr[0]; 
  }
  if (alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[3] = 0.6123724356957942*fc[7]+0.6123724356957944*fc[6]+0.3535533905932737*fc[5]+0.6123724356957944*fc[4]+0.3535533905932736*fc[3]+0.6123724356957944*fc[2]+0.3535533905932736*fc[1]+0.3535533905932737*fc[0]; 
  } else {
    fUpOrdR[3] = (-0.6123724356957942*fr[7])-0.6123724356957944*fr[6]+0.3535533905932737*fr[5]-0.6123724356957944*fr[4]+0.3535533905932736*fr[3]-0.6123724356957944*fr[2]+0.3535533905932736*fr[1]+0.3535533905932737*fr[0]; 
  }

  double fUpR[4] = {0.};
  fUpR[0] = 0.5*fUpOrdR[3]+0.5*fUpOrdR[2]+0.5*fUpOrdR[1]+0.5*fUpOrdR[0]; 
  fUpR[1] = 0.5000000000000001*fUpOrdR[3]+0.5000000000000001*fUpOrdR[2]-0.5000000000000001*fUpOrdR[1]-0.5000000000000001*fUpOrdR[0]; 
  fUpR[2] = 0.5000000000000001*fUpOrdR[3]-0.5000000000000001*fUpOrdR[2]+0.5000000000000001*fUpOrdR[1]-0.5000000000000001*fUpOrdR[0]; 
  fUpR[3] = 0.5*fUpOrdR[3]-0.5*fUpOrdR[2]-0.5*fUpOrdR[1]+0.5*fUpOrdR[0]; 

  double GhatR[8] = {0.}; 
  GhatR[0] += 0.5*alphaR[3]*fUpR[3]+0.5*alphaR[2]*fUpR[2]+0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.5*alphaR[2]*fUpR[3]+0.5*fUpR[2]*alphaR[3]+0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] += 0.5*alphaR[1]*fUpR[3]+0.5*fUpR[1]*alphaR[3]+0.5*alphaR[0]*fUpR[2]+0.5*fUpR[0]*alphaR[2]; 
  GhatR[3] += 0.5*alphaR[0]*fUpR[3]+0.5*fUpR[0]*alphaR[3]+0.5*alphaR[1]*fUpR[2]+0.5*fUpR[1]*alphaR[2]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdvpar2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[5] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdvpar2; 
  out[6] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdvpar2; 
  out[7] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdvpar2; 

} 
