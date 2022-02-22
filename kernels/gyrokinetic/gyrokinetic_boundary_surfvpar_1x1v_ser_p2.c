#include <gkyl_gyrokinetic_kernels.h>
GKYL_CU_DH void gyrokinetic_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fskin, const double *fedge, double* GKYL_RESTRICT out) 
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
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[8]; 
  hamil[0] = (0.3333333333333333*(3.0*rdvpar2Sq*(m_*wvparSq+1.414213562373095*phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = 1.414213562373095*phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.732050807568877*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobtot_inv[2]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_)/q_; 
  BstarZdBmag[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobtot_inv[2]+11.18033988749895*jacobtot_inv[0])+5.0*b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobtot_inv[2]+jacobtot_inv[1]*cmag[2])+5.0*(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = ((2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[3] = ((b_y[2]*(2.0*jacobtot_inv[2]+2.23606797749979*jacobtot_inv[0])+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobtot_inv[2]+7.0*(5.0*jacobtot_inv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobtot_inv[1]))*q_))/q_; 
  BstarZdBmag[6] = (1.0*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  if (edge == -1) { 

  double alphaR[3]; 
  alphaR[0] = -(0.3535533905932737*(((6.708203932499369*BstarZdBmag[3]+3.872983346207417*BstarZdBmag[1])*hamil[4]+hamil[1]*(3.0*BstarZdBmag[2]+1.732050807568877*BstarZdBmag[0]))*rdx2+2.828427124746191*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.07071067811865474*((hamil[4]*(30.0*BstarZdBmag[6]+17.32050807568877*BstarZdBmag[4]+33.54101966249685*BstarZdBmag[2]+19.36491673103709*BstarZdBmag[0])+hamil[1]*(15.0*BstarZdBmag[3]+8.660254037844386*BstarZdBmag[1]))*rdx2+14.14213562373095*apardot[1]*q_))/m_; 
  alphaR[2] = -(0.07071067811865474*((15.0*hamil[1]*BstarZdBmag[6]+(30.0*BstarZdBmag[3]+17.32050807568877*BstarZdBmag[1])*hamil[4]+8.660254037844386*hamil[1]*BstarZdBmag[4])*rdx2+14.14213562373095*apardot[2]*q_))/m_; 

  double fUpOrdR[3];
  if (alphaR[0]-1.118033988749896*alphaR[2] > 0.) {
    fUpOrdR[0] = (-0.9682458365518543*fskin[6])+1.118033988749895*fskin[5]-0.5590169943749475*fskin[4]+0.8660254037844386*fskin[2]+0.5*fskin[0]; 
  } else {
    fUpOrdR[0] = 0.9682458365518543*fedge[6]+1.118033988749895*fedge[5]-0.5590169943749475*fedge[4]-0.8660254037844386*fedge[2]+0.5*fedge[0]; 
  }
  if (alphaR[2]-1.499999999999997*alphaR[1]+1.118033988749892*alphaR[0] > 0.) {
    fUpOrdR[1] = (-1.5*fskin[7])+0.7745966692414834*fskin[6]+1.118033988749895*fskin[5]+0.447213595499958*fskin[4]-1.161895003862225*fskin[3]+0.8660254037844386*fskin[2]-0.6708203932499369*fskin[1]+0.5*fskin[0]; 
  } else {
    fUpOrdR[1] = (-1.5*fedge[7])-0.7745966692414834*fedge[6]+1.118033988749895*fedge[5]+0.447213595499958*fedge[4]+1.161895003862225*fedge[3]-0.8660254037844386*fedge[2]-0.6708203932499369*fedge[1]+0.5*fedge[0]; 
  }
  if (alphaR[2]+1.499999999999997*alphaR[1]+1.118033988749892*alphaR[0] > 0.) {
    fUpOrdR[2] = 1.5*fskin[7]+0.7745966692414834*fskin[6]+1.118033988749895*fskin[5]+0.447213595499958*fskin[4]+1.161895003862225*fskin[3]+0.8660254037844386*fskin[2]+0.6708203932499369*fskin[1]+0.5*fskin[0]; 
  } else {
    fUpOrdR[2] = 1.5*fedge[7]-0.7745966692414834*fedge[6]+1.118033988749895*fedge[5]+0.447213595499958*fedge[4]-1.161895003862225*fedge[3]-0.8660254037844386*fedge[2]+0.6708203932499369*fedge[1]+0.5*fedge[0]; 
  }

  double fUpR[3] = {0.};
  fUpR[0] = 0.392837100659193*fUpOrdR[2]+0.392837100659193*fUpOrdR[1]+0.6285393610547091*fUpOrdR[0]; 
  fUpR[1] = 0.5270462766947293*fUpOrdR[2]-0.5270462766947293*fUpOrdR[1]; 
  fUpR[2] = 0.3513641844631533*fUpOrdR[2]+0.3513641844631533*fUpOrdR[1]-0.7027283689263066*fUpOrdR[0]; 

  double GhatR[8] = {0.}; 
  GhatR[0] += 0.7071067811865475*alphaR[2]*fUpR[2]+0.7071067811865475*alphaR[1]*fUpR[1]+0.7071067811865475*alphaR[0]*fUpR[0]; 
  GhatR[1] += 0.6324555320336759*alphaR[1]*fUpR[2]+0.6324555320336759*fUpR[1]*alphaR[2]+0.7071067811865475*alphaR[0]*fUpR[1]+0.7071067811865475*fUpR[0]*alphaR[1]; 
  GhatR[2] += 0.4517539514526256*alphaR[2]*fUpR[2]+0.7071067811865475*alphaR[0]*fUpR[2]+0.7071067811865475*fUpR[0]*alphaR[2]+0.6324555320336759*alphaR[1]*fUpR[1]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[3] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[4] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[5] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[6] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[7] += -1.58113883008419*GhatR[1]*rdvpar2; 

  } else { 

  double alphaL[3]; 
  alphaL[0] = (0.3535533905932737*(((6.708203932499369*BstarZdBmag[3]-3.872983346207417*BstarZdBmag[1])*hamil[4]+hamil[1]*(3.0*BstarZdBmag[2]-1.732050807568877*BstarZdBmag[0]))*rdx2-2.828427124746191*apardot[0]*q_))/m_; 
  alphaL[1] = (0.07071067811865474*((hamil[4]*(30.0*BstarZdBmag[6]-17.32050807568877*BstarZdBmag[4]+33.54101966249685*BstarZdBmag[2]-19.36491673103709*BstarZdBmag[0])+hamil[1]*(15.0*BstarZdBmag[3]-8.660254037844386*BstarZdBmag[1]))*rdx2-14.14213562373095*apardot[1]*q_))/m_; 
  alphaL[2] = (0.07071067811865474*((15.0*hamil[1]*BstarZdBmag[6]+(30.0*BstarZdBmag[3]-17.32050807568877*BstarZdBmag[1])*hamil[4]-8.660254037844386*hamil[1]*BstarZdBmag[4])*rdx2-14.14213562373095*apardot[2]*q_))/m_; 

  double fUpOrdL[3];
  if (alphaL[0]-1.118033988749896*alphaL[2] > 0.) {
    fUpOrdL[0] = (-0.9682458365518543*fedge[6])+1.118033988749895*fedge[5]-0.5590169943749475*fedge[4]+0.8660254037844386*fedge[2]+0.5*fedge[0]; 
  } else {
    fUpOrdL[0] = 0.9682458365518543*fskin[6]+1.118033988749895*fskin[5]-0.5590169943749475*fskin[4]-0.8660254037844386*fskin[2]+0.5*fskin[0]; 
  }
  if (alphaL[2]-1.499999999999997*alphaL[1]+1.118033988749892*alphaL[0] > 0.) {
    fUpOrdL[1] = (-1.5*fedge[7])+0.7745966692414834*fedge[6]+1.118033988749895*fedge[5]+0.447213595499958*fedge[4]-1.161895003862225*fedge[3]+0.8660254037844386*fedge[2]-0.6708203932499369*fedge[1]+0.5*fedge[0]; 
  } else {
    fUpOrdL[1] = (-1.5*fskin[7])-0.7745966692414834*fskin[6]+1.118033988749895*fskin[5]+0.447213595499958*fskin[4]+1.161895003862225*fskin[3]-0.8660254037844386*fskin[2]-0.6708203932499369*fskin[1]+0.5*fskin[0]; 
  }
  if (alphaL[2]+1.499999999999997*alphaL[1]+1.118033988749892*alphaL[0] > 0.) {
    fUpOrdL[2] = 1.5*fedge[7]+0.7745966692414834*fedge[6]+1.118033988749895*fedge[5]+0.447213595499958*fedge[4]+1.161895003862225*fedge[3]+0.8660254037844386*fedge[2]+0.6708203932499369*fedge[1]+0.5*fedge[0]; 
  } else {
    fUpOrdL[2] = 1.5*fskin[7]-0.7745966692414834*fskin[6]+1.118033988749895*fskin[5]+0.447213595499958*fskin[4]-1.161895003862225*fskin[3]-0.8660254037844386*fskin[2]+0.6708203932499369*fskin[1]+0.5*fskin[0]; 
  }

  double fUpL[3] = {0.};
  fUpL[0] = 0.392837100659193*fUpOrdL[2]+0.392837100659193*fUpOrdL[1]+0.6285393610547091*fUpOrdL[0]; 
  fUpL[1] = 0.5270462766947293*fUpOrdL[2]-0.5270462766947293*fUpOrdL[1]; 
  fUpL[2] = 0.3513641844631533*fUpOrdL[2]+0.3513641844631533*fUpOrdL[1]-0.7027283689263066*fUpOrdL[0]; 

  double GhatL[8] = {0.}; 
  GhatL[0] += 0.7071067811865475*alphaL[2]*fUpL[2]+0.7071067811865475*alphaL[1]*fUpL[1]+0.7071067811865475*alphaL[0]*fUpL[0]; 
  GhatL[1] += 0.6324555320336759*alphaL[1]*fUpL[2]+0.6324555320336759*fUpL[1]*alphaL[2]+0.7071067811865475*alphaL[0]*fUpL[1]+0.7071067811865475*fUpL[0]*alphaL[1]; 
  GhatL[2] += 0.4517539514526256*alphaL[2]*fUpL[2]+0.7071067811865475*alphaL[0]*fUpL[2]+0.7071067811865475*fUpL[0]*alphaL[2]+0.6324555320336759*alphaL[1]*fUpL[1]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[3] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[4] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[5] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[6] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[7] += 1.58113883008419*GhatL[1]*rdvpar2; 

  } 

} 
