#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[4];
  const double *b_z = &b_i[8];

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = 0.8660254037844386*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[2])*rdy2; 
  BstarXdBmag[1] = 0.1732050807568877*((9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[3]+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[2])*rdy2; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -(0.8660254037844386*rdx2*(2.0*jacobtot_inv[0]*b_z[1]*m_*wvpar+(2.0*apar[1]*b_z[1]*jacobtot_inv[1]+jacobtot_inv[0]*(apar[0]*b_z[1]+b_z[0]*apar[1]))*q_))/q_; 
  BstarYdBmag[1] = -(0.8660254037844386*rdx2*(2.0*b_z[1]*jacobtot_inv[1]*m_*wvpar+((apar[0]*b_z[1]+b_z[0]*apar[1])*jacobtot_inv[1]+2.0*jacobtot_inv[0]*apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmag[2] = -0.8660254037844386*((2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[3]+jacobtot_inv[0]*b_z[1]*apar[2])*rdx2; 
  BstarYdBmag[3] = -(1.0*jacobtot_inv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[5] = -0.8660254037844386*((b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1])*apar[3]+b_z[1]*jacobtot_inv[1]*apar[2])*rdx2; 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  if (edge == -1) { 

  double alphaR[8]; 
  alphaR[0] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[1]*rdx2)+3.0*(hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdy2+8.0*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+BstarXdBmag[1]*hamil[1]*rdx2)+3.0*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5])*rdy2+8.0*apardot[1]*q_))/m_; 
  alphaR[2] = -(0.1020620726159657*(3.0*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[5]*rdx2)+13.85640646055102*apardot[2]*q_))/m_; 
  alphaR[3] = -(0.3061862178478971*BstarXdBmag[0]*hamil[8]*rdx2)/m_; 
  alphaR[4] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+BstarXdBmag[1]*hamil[5]*rdx2)+8.0*apardot[3]*q_))/m_; 
  alphaR[5] = -(0.3061862178478971*BstarXdBmag[1]*hamil[8]*rdx2)/m_; 

  double fUpOrdR[8];
  if (alphaR[5]+alphaR[4]-1.0*alphaR[3]-1.0*alphaR[2]-1.0*alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[0] = (-0.4330127018922193*fskin[15])+0.4330127018922192*(fskin[14]+fskin[13])-0.25*fskin[12]+0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]+0.25*(fskin[9]+fskin[8])-0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]-0.25*fskin[4]+0.4330127018922193*fskin[3]-0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  } else {
    fUpOrdR[0] = 0.4330127018922193*fedge[15]-0.4330127018922192*(fedge[14]+fedge[13])-0.25*fedge[12]-0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]+0.25*(fedge[9]+fedge[8])+0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]-0.25*fedge[4]-0.4330127018922193*fedge[3]-0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  }
  if ((-1.0*alphaR[5])+alphaR[4]+alphaR[3]-1.0*alphaR[2]-1.0*alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[1] = 0.4330127018922193*fskin[15]-0.4330127018922192*(fskin[14]+fskin[13])+0.25*fskin[12]+0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]-0.25*(fskin[9]+fskin[8])-0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]+0.25*fskin[4]+0.4330127018922193*fskin[3]-0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  } else {
    fUpOrdR[1] = (-0.4330127018922193*fedge[15])+0.4330127018922192*(fedge[14]+fedge[13])+0.25*fedge[12]-0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]-0.25*(fedge[9]+fedge[8])+0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]+0.25*fedge[4]-0.4330127018922193*fedge[3]-0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  }
  if (alphaR[5]-1.0*alphaR[4]-1.0*alphaR[3]+alphaR[2]-1.0*alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[2] = 0.4330127018922193*fskin[15]-0.4330127018922192*fskin[14]+0.4330127018922192*fskin[13]+0.25*fskin[12]-0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]-0.25*fskin[9]+0.25*fskin[8]+0.4330127018922193*fskin[7]-0.4330127018922193*fskin[6]-0.25*fskin[5]-0.25*fskin[4]+0.4330127018922193*fskin[3]+0.25*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  } else {
    fUpOrdR[2] = (-0.4330127018922193*fedge[15])+0.4330127018922192*fedge[14]-0.4330127018922192*fedge[13]+0.25*fedge[12]+0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]-0.25*fedge[9]+0.25*fedge[8]-0.4330127018922193*fedge[7]+0.4330127018922193*fedge[6]-0.25*fedge[5]-0.25*fedge[4]-0.4330127018922193*fedge[3]+0.25*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  }
  if ((-1.0*alphaR[5])-1.0*alphaR[4]+alphaR[3]+alphaR[2]-1.0*alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[3] = (-0.4330127018922193*fskin[15])+0.4330127018922192*fskin[14]-0.4330127018922192*fskin[13]-0.25*fskin[12]-0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]+0.25*fskin[9]-0.25*fskin[8]+0.4330127018922193*fskin[7]-0.4330127018922193*fskin[6]-0.25*fskin[5]+0.25*fskin[4]+0.4330127018922193*fskin[3]+0.25*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  } else {
    fUpOrdR[3] = 0.4330127018922193*fedge[15]-0.4330127018922192*fedge[14]+0.4330127018922192*fedge[13]-0.25*fedge[12]+0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]+0.25*fedge[9]-0.25*fedge[8]-0.4330127018922193*fedge[7]+0.4330127018922193*fedge[6]-0.25*fedge[5]+0.25*fedge[4]-0.4330127018922193*fedge[3]+0.25*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  }
  if ((-1.0*alphaR[5])-1.0*alphaR[4]-1.0*alphaR[3]-1.0*alphaR[2]+alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[4] = 0.4330127018922193*fskin[15]+0.4330127018922192*fskin[14]-0.4330127018922192*fskin[13]+0.25*fskin[12]-0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]+0.25*fskin[9]-0.25*fskin[8]-0.4330127018922193*fskin[7]+0.4330127018922193*fskin[6]-0.25*fskin[5]-0.25*fskin[4]+0.4330127018922193*fskin[3]-0.25*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  } else {
    fUpOrdR[4] = (-0.4330127018922193*fedge[15])-0.4330127018922192*fedge[14]+0.4330127018922192*fedge[13]+0.25*fedge[12]+0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]+0.25*fedge[9]-0.25*fedge[8]+0.4330127018922193*fedge[7]-0.4330127018922193*fedge[6]-0.25*fedge[5]-0.25*fedge[4]-0.4330127018922193*fedge[3]-0.25*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  }
  if (alphaR[5]-1.0*alphaR[4]+alphaR[3]-1.0*alphaR[2]+alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[5] = (-0.4330127018922193*fskin[15])-0.4330127018922192*fskin[14]+0.4330127018922192*fskin[13]-0.25*fskin[12]-0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]-0.25*fskin[9]+0.25*fskin[8]-0.4330127018922193*fskin[7]+0.4330127018922193*fskin[6]-0.25*fskin[5]+0.25*fskin[4]+0.4330127018922193*fskin[3]-0.25*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  } else {
    fUpOrdR[5] = 0.4330127018922193*fedge[15]+0.4330127018922192*fedge[14]-0.4330127018922192*fedge[13]-0.25*fedge[12]+0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]-0.25*fedge[9]+0.25*fedge[8]+0.4330127018922193*fedge[7]-0.4330127018922193*fedge[6]-0.25*fedge[5]+0.25*fedge[4]-0.4330127018922193*fedge[3]-0.25*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  }
  if ((-1.0*alphaR[5])+alphaR[4]-1.0*alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[6] = (-0.4330127018922193*fskin[15])-0.4330127018922192*(fskin[14]+fskin[13])-0.25*fskin[12]+0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]-0.25*(fskin[9]+fskin[8])+0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]-0.25*fskin[4]+0.4330127018922193*fskin[3]+0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  } else {
    fUpOrdR[6] = 0.4330127018922193*fedge[15]+0.4330127018922192*(fedge[14]+fedge[13])-0.25*fedge[12]-0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]-0.25*(fedge[9]+fedge[8])-0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]-0.25*fedge[4]-0.4330127018922193*fedge[3]+0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  }
  if (alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0] > 0.) {
    fUpOrdR[7] = 0.4330127018922193*fskin[15]+0.4330127018922192*(fskin[14]+fskin[13])+0.25*fskin[12]+0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]+0.25*(fskin[9]+fskin[8])+0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]+0.25*fskin[4]+0.4330127018922193*fskin[3]+0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  } else {
    fUpOrdR[7] = (-0.4330127018922193*fedge[15])-0.4330127018922192*(fedge[14]+fedge[13])+0.25*fedge[12]-0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]+0.25*(fedge[9]+fedge[8])-0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]+0.25*fedge[4]-0.4330127018922193*fedge[3]+0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  }

  double fUpR[8] = {0.};
  fUpR[0] = 0.3535533905932737*fUpOrdR[7]+0.3535533905932737*fUpOrdR[6]+0.3535533905932737*fUpOrdR[5]+0.3535533905932737*fUpOrdR[4]+0.3535533905932737*fUpOrdR[3]+0.3535533905932737*fUpOrdR[2]+0.3535533905932737*fUpOrdR[1]+0.3535533905932737*fUpOrdR[0]; 
  fUpR[1] = 0.3535533905932736*fUpOrdR[7]+0.3535533905932736*fUpOrdR[6]+0.3535533905932736*fUpOrdR[5]+0.3535533905932736*fUpOrdR[4]-0.3535533905932736*fUpOrdR[3]-0.3535533905932736*fUpOrdR[2]-0.3535533905932736*fUpOrdR[1]-0.3535533905932736*fUpOrdR[0]; 
  fUpR[2] = 0.3535533905932736*fUpOrdR[7]+0.3535533905932736*fUpOrdR[6]-0.3535533905932736*fUpOrdR[5]-0.3535533905932736*fUpOrdR[4]+0.3535533905932736*fUpOrdR[3]+0.3535533905932736*fUpOrdR[2]-0.3535533905932736*fUpOrdR[1]-0.3535533905932736*fUpOrdR[0]; 
  fUpR[3] = 0.3535533905932736*fUpOrdR[7]-0.3535533905932736*fUpOrdR[6]+0.3535533905932736*fUpOrdR[5]-0.3535533905932736*fUpOrdR[4]+0.3535533905932736*fUpOrdR[3]-0.3535533905932736*fUpOrdR[2]+0.3535533905932736*fUpOrdR[1]-0.3535533905932736*fUpOrdR[0]; 
  fUpR[4] = 0.3535533905932737*fUpOrdR[7]+0.3535533905932737*fUpOrdR[6]-0.3535533905932737*fUpOrdR[5]-0.3535533905932737*fUpOrdR[4]-0.3535533905932737*fUpOrdR[3]-0.3535533905932737*fUpOrdR[2]+0.3535533905932737*fUpOrdR[1]+0.3535533905932737*fUpOrdR[0]; 
  fUpR[5] = 0.3535533905932737*fUpOrdR[7]-0.3535533905932737*fUpOrdR[6]+0.3535533905932737*fUpOrdR[5]-0.3535533905932737*fUpOrdR[4]-0.3535533905932737*fUpOrdR[3]+0.3535533905932737*fUpOrdR[2]-0.3535533905932737*fUpOrdR[1]+0.3535533905932737*fUpOrdR[0]; 
  fUpR[6] = 0.3535533905932737*fUpOrdR[7]-0.3535533905932737*fUpOrdR[6]-0.3535533905932737*fUpOrdR[5]+0.3535533905932737*fUpOrdR[4]+0.3535533905932737*fUpOrdR[3]-0.3535533905932737*fUpOrdR[2]-0.3535533905932737*fUpOrdR[1]+0.3535533905932737*fUpOrdR[0]; 
  fUpR[7] = 0.3535533905932736*fUpOrdR[7]-0.3535533905932736*fUpOrdR[6]-0.3535533905932736*fUpOrdR[5]+0.3535533905932736*fUpOrdR[4]-0.3535533905932736*fUpOrdR[3]+0.3535533905932736*fUpOrdR[2]+0.3535533905932736*fUpOrdR[1]-0.3535533905932736*fUpOrdR[0]; 

  double GhatR[16] = {0.}; 
  GhatR[0] += 0.3535533905932737*alphaR[5]*fUpR[5]+0.3535533905932737*alphaR[4]*fUpR[4]+0.3535533905932737*alphaR[3]*fUpR[3]+0.3535533905932737*alphaR[2]*fUpR[2]+0.3535533905932737*alphaR[1]*fUpR[1]+0.3535533905932737*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.3535533905932737*alphaR[3]*fUpR[5]+0.3535533905932737*fUpR[3]*alphaR[5]+0.3535533905932737*alphaR[2]*fUpR[4]+0.3535533905932737*fUpR[2]*alphaR[4]+0.3535533905932737*alphaR[0]*fUpR[1]+0.3535533905932737*fUpR[0]*alphaR[1]; 
  GhatR[2] += 0.3535533905932737*alphaR[5]*fUpR[7]+0.3535533905932737*alphaR[3]*fUpR[6]+0.3535533905932737*alphaR[1]*fUpR[4]+0.3535533905932737*fUpR[1]*alphaR[4]+0.3535533905932737*alphaR[0]*fUpR[2]+0.3535533905932737*fUpR[0]*alphaR[2]; 
  GhatR[3] += 0.3535533905932737*alphaR[4]*fUpR[7]+0.3535533905932737*alphaR[2]*fUpR[6]+0.3535533905932737*alphaR[1]*fUpR[5]+0.3535533905932737*fUpR[1]*alphaR[5]+0.3535533905932737*alphaR[0]*fUpR[3]+0.3535533905932737*fUpR[0]*alphaR[3]; 
  GhatR[4] += 0.3535533905932737*alphaR[3]*fUpR[7]+0.3535533905932737*alphaR[5]*fUpR[6]+0.3535533905932737*alphaR[0]*fUpR[4]+0.3535533905932737*fUpR[0]*alphaR[4]+0.3535533905932737*alphaR[1]*fUpR[2]+0.3535533905932737*fUpR[1]*alphaR[2]; 
  GhatR[5] += 0.3535533905932737*alphaR[2]*fUpR[7]+0.3535533905932737*alphaR[4]*fUpR[6]+0.3535533905932737*alphaR[0]*fUpR[5]+0.3535533905932737*fUpR[0]*alphaR[5]+0.3535533905932737*alphaR[1]*fUpR[3]+0.3535533905932737*fUpR[1]*alphaR[3]; 
  GhatR[6] += 0.3535533905932737*alphaR[1]*fUpR[7]+0.3535533905932737*alphaR[0]*fUpR[6]+0.3535533905932737*alphaR[4]*fUpR[5]+0.3535533905932737*fUpR[4]*alphaR[5]+0.3535533905932737*alphaR[2]*fUpR[3]+0.3535533905932737*fUpR[2]*alphaR[3]; 
  GhatR[7] += 0.3535533905932737*alphaR[0]*fUpR[7]+0.3535533905932737*alphaR[1]*fUpR[6]+0.3535533905932737*alphaR[2]*fUpR[5]+0.3535533905932737*fUpR[2]*alphaR[5]+0.3535533905932737*alphaR[3]*fUpR[4]+0.3535533905932737*fUpR[3]*alphaR[4]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[3] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[6] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[7] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[8] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[9] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[12] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[13] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[14] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[15] += -1.224744871391589*GhatR[7]*rdvpar2; 

  } else { 

  double alphaL[8]; 
  alphaL[0] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[1]*rdx2)-3.0*(hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdy2+8.0*apardot[0]*q_))/m_; 
  alphaL[1] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+BstarXdBmag[1]*hamil[1]*rdx2)-3.0*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5])*rdy2+8.0*apardot[1]*q_))/m_; 
  alphaL[2] = -(0.1020620726159657*(3.0*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[5]*rdx2)+13.85640646055102*apardot[2]*q_))/m_; 
  alphaL[3] = -(0.3061862178478971*BstarXdBmag[0]*hamil[8]*rdx2)/m_; 
  alphaL[4] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+BstarXdBmag[1]*hamil[5]*rdx2)+8.0*apardot[3]*q_))/m_; 
  alphaL[5] = -(0.3061862178478971*BstarXdBmag[1]*hamil[8]*rdx2)/m_; 

  double fUpOrdL[8];
  if (alphaL[5]+alphaL[4]-1.0*alphaL[3]-1.0*alphaL[2]-1.0*alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[0] = (-0.4330127018922193*fedge[15])+0.4330127018922192*(fedge[14]+fedge[13])-0.25*fedge[12]+0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]+0.25*(fedge[9]+fedge[8])-0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]-0.25*fedge[4]+0.4330127018922193*fedge[3]-0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  } else {
    fUpOrdL[0] = 0.4330127018922193*fskin[15]-0.4330127018922192*(fskin[14]+fskin[13])-0.25*fskin[12]-0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]+0.25*(fskin[9]+fskin[8])+0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]-0.25*fskin[4]-0.4330127018922193*fskin[3]-0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  }
  if ((-1.0*alphaL[5])+alphaL[4]+alphaL[3]-1.0*alphaL[2]-1.0*alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[1] = 0.4330127018922193*fedge[15]-0.4330127018922192*(fedge[14]+fedge[13])+0.25*fedge[12]+0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]-0.25*(fedge[9]+fedge[8])-0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]+0.25*fedge[4]+0.4330127018922193*fedge[3]-0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  } else {
    fUpOrdL[1] = (-0.4330127018922193*fskin[15])+0.4330127018922192*(fskin[14]+fskin[13])+0.25*fskin[12]-0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]-0.25*(fskin[9]+fskin[8])+0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]+0.25*fskin[4]-0.4330127018922193*fskin[3]-0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  }
  if (alphaL[5]-1.0*alphaL[4]-1.0*alphaL[3]+alphaL[2]-1.0*alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[2] = 0.4330127018922193*fedge[15]-0.4330127018922192*fedge[14]+0.4330127018922192*fedge[13]+0.25*fedge[12]-0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]-0.25*fedge[9]+0.25*fedge[8]+0.4330127018922193*fedge[7]-0.4330127018922193*fedge[6]-0.25*fedge[5]-0.25*fedge[4]+0.4330127018922193*fedge[3]+0.25*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  } else {
    fUpOrdL[2] = (-0.4330127018922193*fskin[15])+0.4330127018922192*fskin[14]-0.4330127018922192*fskin[13]+0.25*fskin[12]+0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]-0.25*fskin[9]+0.25*fskin[8]-0.4330127018922193*fskin[7]+0.4330127018922193*fskin[6]-0.25*fskin[5]-0.25*fskin[4]-0.4330127018922193*fskin[3]+0.25*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  }
  if ((-1.0*alphaL[5])-1.0*alphaL[4]+alphaL[3]+alphaL[2]-1.0*alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[3] = (-0.4330127018922193*fedge[15])+0.4330127018922192*fedge[14]-0.4330127018922192*fedge[13]-0.25*fedge[12]-0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]+0.25*fedge[9]-0.25*fedge[8]+0.4330127018922193*fedge[7]-0.4330127018922193*fedge[6]-0.25*fedge[5]+0.25*fedge[4]+0.4330127018922193*fedge[3]+0.25*fedge[2]-0.25*fedge[1]+0.25*fedge[0]; 
  } else {
    fUpOrdL[3] = 0.4330127018922193*fskin[15]-0.4330127018922192*fskin[14]+0.4330127018922192*fskin[13]-0.25*fskin[12]+0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]+0.25*fskin[9]-0.25*fskin[8]-0.4330127018922193*fskin[7]+0.4330127018922193*fskin[6]-0.25*fskin[5]+0.25*fskin[4]-0.4330127018922193*fskin[3]+0.25*fskin[2]-0.25*fskin[1]+0.25*fskin[0]; 
  }
  if ((-1.0*alphaL[5])-1.0*alphaL[4]-1.0*alphaL[3]-1.0*alphaL[2]+alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[4] = 0.4330127018922193*fedge[15]+0.4330127018922192*fedge[14]-0.4330127018922192*fedge[13]+0.25*fedge[12]-0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]+0.25*fedge[9]-0.25*fedge[8]-0.4330127018922193*fedge[7]+0.4330127018922193*fedge[6]-0.25*fedge[5]-0.25*fedge[4]+0.4330127018922193*fedge[3]-0.25*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  } else {
    fUpOrdL[4] = (-0.4330127018922193*fskin[15])-0.4330127018922192*fskin[14]+0.4330127018922192*fskin[13]+0.25*fskin[12]+0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]+0.25*fskin[9]-0.25*fskin[8]+0.4330127018922193*fskin[7]-0.4330127018922193*fskin[6]-0.25*fskin[5]-0.25*fskin[4]-0.4330127018922193*fskin[3]-0.25*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  }
  if (alphaL[5]-1.0*alphaL[4]+alphaL[3]-1.0*alphaL[2]+alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[5] = (-0.4330127018922193*fedge[15])-0.4330127018922192*fedge[14]+0.4330127018922192*fedge[13]-0.25*fedge[12]-0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]-0.25*fedge[9]+0.25*fedge[8]-0.4330127018922193*fedge[7]+0.4330127018922193*fedge[6]-0.25*fedge[5]+0.25*fedge[4]+0.4330127018922193*fedge[3]-0.25*fedge[2]+0.25*fedge[1]+0.25*fedge[0]; 
  } else {
    fUpOrdL[5] = 0.4330127018922193*fskin[15]+0.4330127018922192*fskin[14]-0.4330127018922192*fskin[13]-0.25*fskin[12]+0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]-0.25*fskin[9]+0.25*fskin[8]+0.4330127018922193*fskin[7]-0.4330127018922193*fskin[6]-0.25*fskin[5]+0.25*fskin[4]-0.4330127018922193*fskin[3]-0.25*fskin[2]+0.25*fskin[1]+0.25*fskin[0]; 
  }
  if ((-1.0*alphaL[5])+alphaL[4]-1.0*alphaL[3]+alphaL[2]+alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[6] = (-0.4330127018922193*fedge[15])-0.4330127018922192*(fedge[14]+fedge[13])-0.25*fedge[12]+0.4330127018922192*fedge[11]-0.4330127018922193*fedge[10]-0.25*(fedge[9]+fedge[8])+0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]-0.25*fedge[4]+0.4330127018922193*fedge[3]+0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  } else {
    fUpOrdL[6] = 0.4330127018922193*fskin[15]+0.4330127018922192*(fskin[14]+fskin[13])-0.25*fskin[12]-0.4330127018922192*fskin[11]+0.4330127018922193*fskin[10]-0.25*(fskin[9]+fskin[8])-0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]-0.25*fskin[4]-0.4330127018922193*fskin[3]+0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  }
  if (alphaL[5]+alphaL[4]+alphaL[3]+alphaL[2]+alphaL[1]+alphaL[0] > 0.) {
    fUpOrdL[7] = 0.4330127018922193*fedge[15]+0.4330127018922192*(fedge[14]+fedge[13])+0.25*fedge[12]+0.4330127018922192*fedge[11]+0.4330127018922193*fedge[10]+0.25*(fedge[9]+fedge[8])+0.4330127018922193*(fedge[7]+fedge[6])+0.25*fedge[5]+0.25*fedge[4]+0.4330127018922193*fedge[3]+0.25*(fedge[2]+fedge[1])+0.25*fedge[0]; 
  } else {
    fUpOrdL[7] = (-0.4330127018922193*fskin[15])-0.4330127018922192*(fskin[14]+fskin[13])+0.25*fskin[12]-0.4330127018922192*fskin[11]-0.4330127018922193*fskin[10]+0.25*(fskin[9]+fskin[8])-0.4330127018922193*(fskin[7]+fskin[6])+0.25*fskin[5]+0.25*fskin[4]-0.4330127018922193*fskin[3]+0.25*(fskin[2]+fskin[1])+0.25*fskin[0]; 
  }

  double fUpL[8] = {0.};
  fUpL[0] = 0.3535533905932737*fUpOrdL[7]+0.3535533905932737*fUpOrdL[6]+0.3535533905932737*fUpOrdL[5]+0.3535533905932737*fUpOrdL[4]+0.3535533905932737*fUpOrdL[3]+0.3535533905932737*fUpOrdL[2]+0.3535533905932737*fUpOrdL[1]+0.3535533905932737*fUpOrdL[0]; 
  fUpL[1] = 0.3535533905932736*fUpOrdL[7]+0.3535533905932736*fUpOrdL[6]+0.3535533905932736*fUpOrdL[5]+0.3535533905932736*fUpOrdL[4]-0.3535533905932736*fUpOrdL[3]-0.3535533905932736*fUpOrdL[2]-0.3535533905932736*fUpOrdL[1]-0.3535533905932736*fUpOrdL[0]; 
  fUpL[2] = 0.3535533905932736*fUpOrdL[7]+0.3535533905932736*fUpOrdL[6]-0.3535533905932736*fUpOrdL[5]-0.3535533905932736*fUpOrdL[4]+0.3535533905932736*fUpOrdL[3]+0.3535533905932736*fUpOrdL[2]-0.3535533905932736*fUpOrdL[1]-0.3535533905932736*fUpOrdL[0]; 
  fUpL[3] = 0.3535533905932736*fUpOrdL[7]-0.3535533905932736*fUpOrdL[6]+0.3535533905932736*fUpOrdL[5]-0.3535533905932736*fUpOrdL[4]+0.3535533905932736*fUpOrdL[3]-0.3535533905932736*fUpOrdL[2]+0.3535533905932736*fUpOrdL[1]-0.3535533905932736*fUpOrdL[0]; 
  fUpL[4] = 0.3535533905932737*fUpOrdL[7]+0.3535533905932737*fUpOrdL[6]-0.3535533905932737*fUpOrdL[5]-0.3535533905932737*fUpOrdL[4]-0.3535533905932737*fUpOrdL[3]-0.3535533905932737*fUpOrdL[2]+0.3535533905932737*fUpOrdL[1]+0.3535533905932737*fUpOrdL[0]; 
  fUpL[5] = 0.3535533905932737*fUpOrdL[7]-0.3535533905932737*fUpOrdL[6]+0.3535533905932737*fUpOrdL[5]-0.3535533905932737*fUpOrdL[4]-0.3535533905932737*fUpOrdL[3]+0.3535533905932737*fUpOrdL[2]-0.3535533905932737*fUpOrdL[1]+0.3535533905932737*fUpOrdL[0]; 
  fUpL[6] = 0.3535533905932737*fUpOrdL[7]-0.3535533905932737*fUpOrdL[6]-0.3535533905932737*fUpOrdL[5]+0.3535533905932737*fUpOrdL[4]+0.3535533905932737*fUpOrdL[3]-0.3535533905932737*fUpOrdL[2]-0.3535533905932737*fUpOrdL[1]+0.3535533905932737*fUpOrdL[0]; 
  fUpL[7] = 0.3535533905932736*fUpOrdL[7]-0.3535533905932736*fUpOrdL[6]-0.3535533905932736*fUpOrdL[5]+0.3535533905932736*fUpOrdL[4]-0.3535533905932736*fUpOrdL[3]+0.3535533905932736*fUpOrdL[2]+0.3535533905932736*fUpOrdL[1]-0.3535533905932736*fUpOrdL[0]; 

  double GhatL[16] = {0.}; 
  GhatL[0] += 0.3535533905932737*alphaL[5]*fUpL[5]+0.3535533905932737*alphaL[4]*fUpL[4]+0.3535533905932737*alphaL[3]*fUpL[3]+0.3535533905932737*alphaL[2]*fUpL[2]+0.3535533905932737*alphaL[1]*fUpL[1]+0.3535533905932737*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.3535533905932737*alphaL[3]*fUpL[5]+0.3535533905932737*fUpL[3]*alphaL[5]+0.3535533905932737*alphaL[2]*fUpL[4]+0.3535533905932737*fUpL[2]*alphaL[4]+0.3535533905932737*alphaL[0]*fUpL[1]+0.3535533905932737*fUpL[0]*alphaL[1]; 
  GhatL[2] += 0.3535533905932737*alphaL[5]*fUpL[7]+0.3535533905932737*alphaL[3]*fUpL[6]+0.3535533905932737*alphaL[1]*fUpL[4]+0.3535533905932737*fUpL[1]*alphaL[4]+0.3535533905932737*alphaL[0]*fUpL[2]+0.3535533905932737*fUpL[0]*alphaL[2]; 
  GhatL[3] += 0.3535533905932737*alphaL[4]*fUpL[7]+0.3535533905932737*alphaL[2]*fUpL[6]+0.3535533905932737*alphaL[1]*fUpL[5]+0.3535533905932737*fUpL[1]*alphaL[5]+0.3535533905932737*alphaL[0]*fUpL[3]+0.3535533905932737*fUpL[0]*alphaL[3]; 
  GhatL[4] += 0.3535533905932737*alphaL[3]*fUpL[7]+0.3535533905932737*alphaL[5]*fUpL[6]+0.3535533905932737*alphaL[0]*fUpL[4]+0.3535533905932737*fUpL[0]*alphaL[4]+0.3535533905932737*alphaL[1]*fUpL[2]+0.3535533905932737*fUpL[1]*alphaL[2]; 
  GhatL[5] += 0.3535533905932737*alphaL[2]*fUpL[7]+0.3535533905932737*alphaL[4]*fUpL[6]+0.3535533905932737*alphaL[0]*fUpL[5]+0.3535533905932737*fUpL[0]*alphaL[5]+0.3535533905932737*alphaL[1]*fUpL[3]+0.3535533905932737*fUpL[1]*alphaL[3]; 
  GhatL[6] += 0.3535533905932737*alphaL[1]*fUpL[7]+0.3535533905932737*alphaL[0]*fUpL[6]+0.3535533905932737*alphaL[4]*fUpL[5]+0.3535533905932737*fUpL[4]*alphaL[5]+0.3535533905932737*alphaL[2]*fUpL[3]+0.3535533905932737*fUpL[2]*alphaL[3]; 
  GhatL[7] += 0.3535533905932737*alphaL[0]*fUpL[7]+0.3535533905932737*alphaL[1]*fUpL[6]+0.3535533905932737*alphaL[2]*fUpL[5]+0.3535533905932737*fUpL[2]*alphaL[5]+0.3535533905932737*alphaL[3]*fUpL[4]+0.3535533905932737*fUpL[3]*alphaL[4]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[3] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[6] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[7] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[8] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[9] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[12] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[13] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[14] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[15] += -1.224744871391589*GhatL[7]*rdvpar2; 

  } 

} 
