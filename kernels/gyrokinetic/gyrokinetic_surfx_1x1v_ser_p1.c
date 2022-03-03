#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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
  BstarZdBmag[0] = (1.732050807568877*jacobtot_inv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_)/q_; 
  BstarZdBmag[1] = (1.732050807568877*b_y[1]*jacobtot_inv[1]*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_)/q_; 
  BstarZdBmag[2] = (jacobtot_inv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = (b_y[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[2]; 
  alphaL[0] = -(0.3535533905932737*(3.0*BstarZdBmag[1]-1.732050807568877*BstarZdBmag[0])*hamil[2]*rdvpar2)/m_; 
  alphaL[1] = -(0.3535533905932737*hamil[2]*(3.0*BstarZdBmag[3]-1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double alphaR[2]; 
  alphaR[0] = (0.3535533905932737*(3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0])*hamil[2]*rdvpar2)/m_; 
  alphaR[1] = (0.3535533905932737*hamil[2]*(3.0*BstarZdBmag[3]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double fUpOrdL[2];
  if (alphaL[0]-1.0*alphaL[1] > 0.) {
    fUpOrdL[0] = (-0.8660254037844386*fl[3])-0.4999999999999999*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  } else {
    fUpOrdL[0] = 0.8660254037844386*fc[3]-0.4999999999999999*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  }
  if (alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[1] = 0.8660254037844386*fl[3]+0.4999999999999999*fl[2]+0.8660254037844386*fl[1]+0.5*fl[0]; 
  } else {
    fUpOrdL[1] = (-0.8660254037844386*fc[3])+0.4999999999999999*fc[2]-0.8660254037844386*fc[1]+0.5*fc[0]; 
  }

  double fUpL[2] = {0.};
  fUpL[0] = 0.7071067811865475*fUpOrdL[1]+0.7071067811865475*fUpOrdL[0]; 
  fUpL[1] = 0.7071067811865471*fUpOrdL[1]-0.7071067811865471*fUpOrdL[0]; 

  double GhatL[4] = {0.}; 
  GhatL[0] += 0.7071067811865475*alphaL[1]*fUpL[1]+0.7071067811865475*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.7071067811865475*alphaL[0]*fUpL[1]+0.7071067811865475*fUpL[0]*alphaL[1]; 

  double fUpOrdR[2];
  if (alphaR[0]-1.0*alphaR[1] > 0.) {
    fUpOrdR[0] = (-0.8660254037844386*fc[3])-0.4999999999999999*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  } else {
    fUpOrdR[0] = 0.8660254037844386*fr[3]-0.4999999999999999*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  }
  if (alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[1] = 0.8660254037844386*fc[3]+0.4999999999999999*fc[2]+0.8660254037844386*fc[1]+0.5*fc[0]; 
  } else {
    fUpOrdR[1] = (-0.8660254037844386*fr[3])+0.4999999999999999*fr[2]-0.8660254037844386*fr[1]+0.5*fr[0]; 
  }

  double fUpR[2] = {0.};
  fUpR[0] = 0.7071067811865475*fUpOrdR[1]+0.7071067811865475*fUpOrdR[0]; 
  fUpR[1] = 0.7071067811865471*fUpOrdR[1]-0.7071067811865471*fUpOrdR[0]; 

  double GhatR[4] = {0.}; 
  GhatR[0] += 0.7071067811865475*alphaR[1]*fUpR[1]+0.7071067811865475*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.7071067811865475*alphaR[0]*fUpR[1]+0.7071067811865475*fUpR[0]*alphaR[1]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 

} 
