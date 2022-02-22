#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_boundary_surfvpar_1x1v_tensor_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
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

  if (edge == -1) { 

  double alphaR[2]; 
  alphaR[0] = -(0.3535533905932737*(hamil[1]*(3.0*BstarZdBmag[2]+1.732050807568877*BstarZdBmag[0])*rdx2+2.828427124746191*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.3535533905932737*(hamil[1]*(3.0*BstarZdBmag[3]+1.732050807568877*BstarZdBmag[1])*rdx2+2.828427124746191*apardot[1]*q_))/m_; 

  double fUpOrdR[2];
  if (alphaR[0]-1.0*alphaR[1] > 0.) {
    fUpOrdR[0] = (-0.8660254037844386*fskin[3])+0.8660254037844386*fskin[2]-0.4999999999999999*fskin[1]+0.5*fskin[0]; 
  } else {
    fUpOrdR[0] = 0.8660254037844386*fedge[3]-0.8660254037844386*fedge[2]-0.4999999999999999*fedge[1]+0.5*fedge[0]; 
  }
  if (alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[1] = 0.8660254037844386*(fskin[3]+fskin[2])+0.4999999999999999*fskin[1]+0.5*fskin[0]; 
  } else {
    fUpOrdR[1] = (-0.8660254037844386*(fedge[3]+fedge[2]))+0.4999999999999999*fedge[1]+0.5*fedge[0]; 
  }

  double fUpR[2] = {0.};
  fUpR[0] = 0.7071067811865475*fUpOrdR[1]+0.7071067811865475*fUpOrdR[0]; 
  fUpR[1] = 0.7071067811865471*fUpOrdR[1]-0.7071067811865471*fUpOrdR[0]; 

  double GhatR[4] = {0.}; 
  GhatR[0] += 0.7071067811865475*alphaR[1]*fUpR[1]+0.7071067811865475*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.7071067811865475*alphaR[0]*fUpR[1]+0.7071067811865475*fUpR[0]*alphaR[1]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[3] += -1.224744871391589*GhatR[1]*rdvpar2; 

  } else { 

  double alphaL[2]; 
  alphaL[0] = -(0.3535533905932737*(hamil[1]*(1.732050807568877*BstarZdBmag[0]-3.0*BstarZdBmag[2])*rdx2+2.828427124746191*apardot[0]*q_))/m_; 
  alphaL[1] = -(0.3535533905932737*(hamil[1]*(1.732050807568877*BstarZdBmag[1]-3.0*BstarZdBmag[3])*rdx2+2.828427124746191*apardot[1]*q_))/m_; 

  double fUpOrdL[2];
  if (alphaL[0]-1.0*alphaL[1] > 0.) {
    fUpOrdL[0] = (-0.8660254037844386*fedge[3])+0.8660254037844386*fedge[2]-0.4999999999999999*fedge[1]+0.5*fedge[0]; 
  } else {
    fUpOrdL[0] = 0.8660254037844386*fskin[3]-0.8660254037844386*fskin[2]-0.4999999999999999*fskin[1]+0.5*fskin[0]; 
  }
  if (alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[1] = 0.8660254037844386*(fedge[3]+fedge[2])+0.4999999999999999*fedge[1]+0.5*fedge[0]; 
  } else {
    fUpOrdL[1] = (-0.8660254037844386*(fskin[3]+fskin[2]))+0.4999999999999999*fskin[1]+0.5*fskin[0]; 
  }

  double fUpL[2] = {0.};
  fUpL[0] = 0.7071067811865475*fUpOrdL[1]+0.7071067811865475*fUpOrdL[0]; 
  fUpL[1] = 0.7071067811865471*fUpOrdL[1]-0.7071067811865471*fUpOrdL[0]; 

  double GhatL[4] = {0.}; 
  GhatL[0] += 0.7071067811865475*alphaL[1]*fUpL[1]+0.7071067811865475*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.7071067811865475*alphaL[0]*fUpL[1]+0.7071067811865475*fUpL[0]*alphaL[1]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[3] += -1.224744871391589*GhatL[1]*rdvpar2; 

  } 

} 
