#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_4x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void gyrokinetic_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[11] = 2.0*(bmag[4]*wmu+phi[4]*q_); 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 

  double BstarXdBmag[48] = {0.}; 
  BstarXdBmag[0] = 0.01428571428571429*(3.872983346207417*((10.0*b_z[4]+15.65247584249853*b_z[0])*jacobtot_inv[4]+7.0*(2.23606797749979*jacobtot_inv[0]*b_z[4]+2.0*b_z[1]*jacobtot_inv[1]))*apar[6]+12.12435565298214*((5.0*apar[2]*b_z[4]+4.47213595499958*b_z[1]*apar[3])*jacobtot_inv[4]+4.47213595499958*jacobtot_inv[1]*apar[3]*b_z[4]+5.0*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[2])))*rdy2; 
  BstarXdBmag[1] = 0.01428571428571429*(3.872983346207417*(24.59674775249769*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+14.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[6]+1.732050807568877*((55.0*apar[3]*b_z[4]+31.30495168499706*(b_z[0]*apar[3]+b_z[1]*apar[2]))*jacobtot_inv[4]+7.0*(4.47213595499958*(jacobtot_inv[0]*apar[3]+jacobtot_inv[1]*apar[2])*b_z[4]+(9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[3]+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[2])))*rdy2; 
  BstarXdBmag[2] = 0.1*(3.872983346207417*(4.47213595499958*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[7]+19.36491673103709*(b_z[4]*jacobtot_inv[4]+b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[5])*rdy2; 
  BstarXdBmag[5] = 0.01428571428571429*(3.872983346207417*((55.0*b_z[4]+31.30495168499706*b_z[0])*jacobtot_inv[4]+7.0*(4.47213595499958*jacobtot_inv[0]*b_z[4]+9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0]))*apar[7]+12.12435565298214*(10.0*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+11.18033988749895*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[5])*rdy2; 
  BstarXdBmag[11] = 0.01428571428571429*(3.872983346207417*((33.54101966249685*b_z[4]+10.0*b_z[0])*jacobtot_inv[4]+10.0*jacobtot_inv[0]*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobtot_inv[1]+7.0*b_z[0]*jacobtot_inv[0]))*apar[6]+1.732050807568877*((22.3606797749979*apar[2]*b_z[4]+5.0*(11.0*b_z[1]*apar[3]+7.0*b_z[0]*apar[2]))*jacobtot_inv[4]+5.0*(11.0*jacobtot_inv[1]*apar[3]+7.0*jacobtot_inv[0]*apar[2])*b_z[4]+31.30495168499706*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+b_z[1]*jacobtot_inv[1]*apar[2])))*rdy2; 
  BstarXdBmag[19] = 0.01428571428571429*(1.732050807568877*(122.9837387624885*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+70.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[7]+3.872983346207417*((22.3606797749979*b_z[4]+35.0*b_z[0])*jacobtot_inv[4]+7.0*(5.0*jacobtot_inv[0]*b_z[4]+4.47213595499958*b_z[1]*jacobtot_inv[1]))*apar[5])*rdy2; 

  double BstarYdBmag[48] = {0.}; 
  BstarYdBmag[0] = -(0.8660254037844386*rdx2*(2.0*(2.23606797749979*jacobtot_inv[1]*b_z[4]+jacobtot_inv[0]*b_z[1])*m_*wvpar+(3.0*(apar[1]*b_z[4]+b_z[1]*apar[4])*jacobtot_inv[4]+(4.0*jacobtot_inv[1]*apar[4]+2.23606797749979*(apar[0]*jacobtot_inv[1]+jacobtot_inv[0]*apar[1]))*b_z[4]+2.23606797749979*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[4]+2.0*apar[1]*b_z[1]*jacobtot_inv[1]+jacobtot_inv[0]*(apar[0]*b_z[1]+b_z[0]*apar[1]))*q_))/q_; 
  BstarYdBmag[1] = -(0.02474358296526967*rdx2*(14.0*(b_z[4]*(10.0*jacobtot_inv[4]+11.18033988749895*jacobtot_inv[0])+5.0*b_z[1]*jacobtot_inv[1])*m_*wvpar+(2.0*((122.9837387624885*apar[4]+35.0*apar[0])*b_z[4]+7.0*(5.0*b_z[0]*apar[4]+4.47213595499958*apar[1]*b_z[1]))*jacobtot_inv[4]+7.0*((20.0*jacobtot_inv[0]*apar[4]+2.23606797749979*(11.0*apar[1]*jacobtot_inv[1]+5.0*apar[0]*jacobtot_inv[0]))*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[4]+5.0*((apar[0]*b_z[1]+b_z[0]*apar[1])*jacobtot_inv[1]+2.0*jacobtot_inv[0]*apar[1]*b_z[1])))*q_))/q_; 
  BstarYdBmag[2] = -0.05773502691896259*(6.708203932499369*(2.23606797749979*(3.0*b_z[1]*jacobtot_inv[4]+4.0*jacobtot_inv[1]*b_z[4])+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[6]+3.0*(b_z[4]*(15.0*apar[3]*jacobtot_inv[4]+11.18033988749895*(jacobtot_inv[0]*apar[3]+jacobtot_inv[1]*apar[2]))+5.0*((2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[3]+jacobtot_inv[0]*b_z[1]*apar[2])))*rdx2; 
  BstarYdBmag[3] = -(1.0*(2.23606797749979*jacobtot_inv[1]*b_z[4]+jacobtot_inv[0]*b_z[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[5] = -0.01428571428571429*(3.872983346207417*(2.0*(55.0*b_z[4]+15.65247584249853*b_z[0])*jacobtot_inv[4]+7.0*(8.94427190999916*jacobtot_inv[0]*b_z[4]+11.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0]))*apar[6]+12.12435565298214*(2.0*(5.0*apar[2]*b_z[4]+4.47213595499958*b_z[1]*apar[3])*jacobtot_inv[4]+2.23606797749979*(11.0*jacobtot_inv[1]*apar[3]+5.0*jacobtot_inv[0]*apar[2])*b_z[4]+5.0*((b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1])*apar[3]+b_z[1]*jacobtot_inv[1]*apar[2])))*rdx2; 
  BstarYdBmag[6] = -(1.0*(b_z[4]*(2.0*jacobtot_inv[4]+2.23606797749979*jacobtot_inv[0])+b_z[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[11] = -(0.02474358296526967*rdx2*(70.0*(b_z[1]*jacobtot_inv[4]+2.0*jacobtot_inv[1]*b_z[4])*m_*wvpar+((145.3444185374863*(apar[1]*b_z[4]+b_z[1]*apar[4])+35.0*(apar[0]*b_z[1]+b_z[0]*apar[1]))*jacobtot_inv[4]+(245.9674775249769*jacobtot_inv[1]*apar[4]+35.0*(2.0*apar[0]*jacobtot_inv[1]+3.0*jacobtot_inv[0]*apar[1]))*b_z[4]+7.0*(5.0*(2.0*b_z[0]*jacobtot_inv[1]+3.0*jacobtot_inv[0]*b_z[1])*apar[4]+8.94427190999916*apar[1]*b_z[1]*jacobtot_inv[1]))*q_))/q_; 
  BstarYdBmag[12] = -0.1*(3.872983346207417*(b_z[4]*(6.708203932499369*jacobtot_inv[4]+5.0*jacobtot_inv[0])+2.23606797749979*(2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0]))*apar[7]+1.732050807568877*(11.18033988749895*jacobtot_inv[1]*b_z[4]+5.0*jacobtot_inv[0]*b_z[1])*apar[5])*rdx2; 
  BstarYdBmag[19] = -0.01844277783908294*(6.708203932499369*(2.23606797749979*(13.0*b_z[1]*jacobtot_inv[4]+22.0*jacobtot_inv[1]*b_z[4])+7.0*(2.0*b_z[0]*jacobtot_inv[1]+3.0*jacobtot_inv[0]*b_z[1]))*apar[6]+3.0*((65.0*apar[3]*b_z[4]+15.65247584249853*(b_z[0]*apar[3]+b_z[1]*apar[2]))*jacobtot_inv[4]+7.0*(2.23606797749979*(3.0*jacobtot_inv[0]*apar[3]+2.0*jacobtot_inv[1]*apar[2])*b_z[4]+4.0*b_z[1]*jacobtot_inv[1]*apar[3])))*rdx2; 
  BstarYdBmag[20] = -0.1*(1.732050807568877*(2.23606797749979*(4.0*b_z[1]*jacobtot_inv[4]+11.0*jacobtot_inv[1]*b_z[4])+5.0*(b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1]))*apar[7]+3.872983346207417*(b_z[4]*(4.47213595499958*jacobtot_inv[4]+5.0*jacobtot_inv[0])+2.23606797749979*b_z[1]*jacobtot_inv[1])*apar[5])*rdx2; 
  BstarYdBmag[21] = -(1.0*(b_z[1]*jacobtot_inv[4]+2.0*jacobtot_inv[1]*b_z[4])*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[20] = {0.}; 
  alphaL[0] = (0.03535533905932736*((15.0*hamil[19]*BstarYdBmag[21]-19.36491673103708*BstarYdBmag[5]*hamil[20]-8.660254037844387*BstarYdBmag[11]*hamil[19]-19.36491673103709*BstarYdBmag[2]*hamil[12]+hamil[5]*(15.0*BstarYdBmag[6]-8.660254037844386*BstarYdBmag[1])+hamil[2]*(15.0*BstarYdBmag[3]-8.660254037844386*BstarYdBmag[0]))*rdy2+((-19.36491673103708*BstarXdBmag[5]*hamil[19])-19.36491673103709*BstarXdBmag[1]*hamil[11]-8.660254037844386*(BstarXdBmag[2]*hamil[5]+BstarXdBmag[0]*hamil[1]))*rdx2-40.0*apardot[0]*q_))/m_; 
  alphaL[1] = (0.03535533905932736*((13.41640786499874*hamil[5]*BstarYdBmag[21]+((-17.32050807568877*BstarYdBmag[19])-19.36491673103708*BstarYdBmag[2])*hamil[20]+(13.41640786499874*BstarYdBmag[6]-7.745966692414834*BstarYdBmag[1])*hamil[19]-19.36491673103709*BstarYdBmag[5]*hamil[12]-7.745966692414834*hamil[5]*BstarYdBmag[11]+15.0*hamil[2]*BstarYdBmag[6]+(15.0*BstarYdBmag[3]-8.660254037844386*BstarYdBmag[0])*hamil[5]-8.660254037844386*BstarYdBmag[1]*hamil[2])*rdy2+(((-17.32050807568877*BstarXdBmag[19])-19.36491673103708*BstarXdBmag[2])*hamil[19]+((-17.32050807568877*BstarXdBmag[11])-19.36491673103709*BstarXdBmag[0])*hamil[11]-8.660254037844386*(BstarXdBmag[5]*hamil[5]+BstarXdBmag[1]*hamil[1]))*rdx2-40.0*apardot[1]*q_))/m_; 
  alphaL[2] = -(0.03535533905932736*(((17.32050807568877*BstarYdBmag[20]-33.54101966249684*BstarYdBmag[6]+19.36491673103708*BstarYdBmag[1])*hamil[20]+8.660254037844386*BstarYdBmag[19]*hamil[19]+(17.32050807568877*BstarYdBmag[12]-33.54101966249685*BstarYdBmag[3]+19.36491673103709*BstarYdBmag[0])*hamil[12]+8.660254037844386*(BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2]))*rdy2+(7.745966692414834*BstarXdBmag[2]*hamil[20]+19.36491673103708*BstarXdBmag[1]*hamil[19]+19.36491673103709*BstarXdBmag[5]*hamil[11]+8.660254037844386*(BstarXdBmag[0]*hamil[5]+hamil[1]*BstarXdBmag[2]))*rdx2+40.0*apardot[2]*q_))/m_; 
  alphaL[3] = -(0.1767766952966368*(3.872983346207417*BstarXdBmag[1]*hamil[25]+1.732050807568877*BstarXdBmag[0]*hamil[8])*rdx2)/m_; 
  alphaL[4] = (0.03535533905932736*((hamil[20]*(30.0*BstarYdBmag[21]-17.32050807568877*(BstarYdBmag[12]+BstarYdBmag[11])+33.54101966249684*BstarYdBmag[3]-19.36491673103708*BstarYdBmag[0])-17.32050807568877*hamil[12]*BstarYdBmag[20]-7.745966692414834*(BstarYdBmag[5]*hamil[19]+hamil[5]*BstarYdBmag[19])+(33.54101966249685*BstarYdBmag[6]-19.36491673103709*BstarYdBmag[1])*hamil[12]-8.660254037844386*(BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5]))*rdy2+((-7.745966692414834*BstarXdBmag[5]*hamil[20])+((-17.32050807568877*BstarXdBmag[11])-19.36491673103708*BstarXdBmag[0])*hamil[19]+hamil[11]*((-17.32050807568877*BstarXdBmag[19])-19.36491673103709*BstarXdBmag[2])-8.660254037844386*(BstarXdBmag[1]*hamil[5]+hamil[1]*BstarXdBmag[5]))*rdx2-40.0*apardot[3]*q_))/m_; 
  alphaL[5] = -(0.03535533905932736*((17.32050807568877*BstarXdBmag[11]+19.36491673103708*BstarXdBmag[0])*hamil[25]+8.660254037844386*BstarXdBmag[1]*hamil[8])*rdx2)/m_; 
  alphaL[6] = -(0.1767766952966368*(3.872983346207417*BstarXdBmag[5]*hamil[25]+1.732050807568877*BstarXdBmag[2]*hamil[8])*rdx2)/m_; 
  alphaL[7] = (0.005050762722761052*(((67.0820393249937*hamil[19]+105.0*hamil[2])*BstarYdBmag[21]-121.2435565298214*BstarYdBmag[5]*hamil[20]+((-38.72983346207417*BstarYdBmag[11])+105.0*BstarYdBmag[3]-60.62177826491071*BstarYdBmag[0])*hamil[19]-135.5544171172596*hamil[12]*BstarYdBmag[19]-60.6217782649107*hamil[2]*BstarYdBmag[11]+hamil[5]*(93.91485505499116*BstarYdBmag[6]-54.22176684690384*BstarYdBmag[1]))*rdy2+((-121.2435565298214*BstarXdBmag[5]*hamil[19])-60.62177826491071*hamil[5]*BstarXdBmag[19]-121.2435565298214*BstarXdBmag[1]*hamil[11]-60.6217782649107*hamil[1]*BstarXdBmag[11])*rdx2-280.0*apardot[4]*q_))/m_; 
  alphaL[8] = -(0.07905694150420944*((7.745966692414834*BstarYdBmag[5]*hamil[20]+3.872983346207417*hamil[5]*BstarYdBmag[20]+7.745966692414834*BstarYdBmag[2]*hamil[12]+3.872983346207417*hamil[2]*BstarYdBmag[12])*rdy2+(3.872983346207417*BstarXdBmag[0]*hamil[20]+7.745966692414834*BstarXdBmag[5]*hamil[19]+3.464101615137754*BstarXdBmag[2]*hamil[5])*rdx2+17.88854381999832*apardot[5]*q_))/m_; 
  alphaL[10] = -(0.1767766952966368*((3.464101615137754*BstarXdBmag[19]+3.872983346207417*BstarXdBmag[2])*hamil[25]+1.732050807568877*BstarXdBmag[5]*hamil[8])*rdx2)/m_; 
  alphaL[11] = (0.005050762722761052*((234.787137637478*hamil[12]*BstarYdBmag[21]+((-108.4435336938077*BstarYdBmag[20])+210.0*BstarYdBmag[6]-121.2435565298214*BstarYdBmag[1])*hamil[20]+((-38.72983346207417*BstarYdBmag[19])-60.6217782649107*BstarYdBmag[2])*hamil[19]-60.6217782649107*hamil[2]*BstarYdBmag[19]-135.5544171172596*BstarYdBmag[11]*hamil[12]-54.22176684690384*BstarYdBmag[5]*hamil[5])*rdy2+((-54.22176684690384*BstarXdBmag[19]*hamil[20])-121.2435565298214*BstarXdBmag[1]*hamil[19]-60.6217782649107*hamil[1]*BstarXdBmag[19]-121.2435565298214*BstarXdBmag[5]*hamil[11]-60.62177826491071*hamil[5]*BstarXdBmag[11])*rdx2-280.0*apardot[6]*q_))/m_; 
  alphaL[12] = -(0.04564354645876382*(((12.0*BstarYdBmag[19]+13.41640786499874*BstarYdBmag[2])*hamil[20]+(6.0*hamil[19]+6.708203932499369*hamil[2])*BstarYdBmag[20]+13.41640786499874*BstarYdBmag[5]*hamil[12]+6.708203932499369*hamil[5]*BstarYdBmag[12])*rdy2+(6.708203932499369*BstarXdBmag[1]*hamil[20]+(12.0*BstarXdBmag[19]+13.41640786499874*BstarXdBmag[2])*hamil[19]+6.0*BstarXdBmag[5]*hamil[5])*rdx2+30.98386676965934*apardot[7]*q_))/m_; 
  alphaL[13] = -(0.03535533905932736*(17.32050807568877*BstarXdBmag[1]*hamil[25]+8.660254037844387*hamil[8]*BstarXdBmag[11])*rdx2)/m_; 
  alphaL[17] = -(0.03535533905932736*(17.32050807568877*BstarXdBmag[5]*hamil[25]+8.660254037844387*hamil[8]*BstarXdBmag[19])*rdx2)/m_; 

  double alphaR[20] = {0.}; 
  alphaR[0] = -(0.03535533905932736*((15.0*hamil[19]*BstarYdBmag[21]+19.36491673103708*BstarYdBmag[5]*hamil[20]+8.660254037844387*BstarYdBmag[11]*hamil[19]+19.36491673103709*BstarYdBmag[2]*hamil[12]+hamil[5]*(15.0*BstarYdBmag[6]+8.660254037844386*BstarYdBmag[1])+hamil[2]*(15.0*BstarYdBmag[3]+8.660254037844386*BstarYdBmag[0]))*rdy2+(19.36491673103708*BstarXdBmag[5]*hamil[19]+19.36491673103709*BstarXdBmag[1]*hamil[11]+8.660254037844386*(BstarXdBmag[2]*hamil[5]+BstarXdBmag[0]*hamil[1]))*rdx2+40.0*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.03535533905932736*((13.41640786499874*hamil[5]*BstarYdBmag[21]+(17.32050807568877*BstarYdBmag[19]+19.36491673103708*BstarYdBmag[2])*hamil[20]+(13.41640786499874*BstarYdBmag[6]+7.745966692414834*BstarYdBmag[1])*hamil[19]+19.36491673103709*BstarYdBmag[5]*hamil[12]+7.745966692414834*hamil[5]*BstarYdBmag[11]+15.0*hamil[2]*BstarYdBmag[6]+(15.0*BstarYdBmag[3]+8.660254037844386*BstarYdBmag[0])*hamil[5]+8.660254037844386*BstarYdBmag[1]*hamil[2])*rdy2+((17.32050807568877*BstarXdBmag[19]+19.36491673103708*BstarXdBmag[2])*hamil[19]+(17.32050807568877*BstarXdBmag[11]+19.36491673103709*BstarXdBmag[0])*hamil[11]+8.660254037844386*(BstarXdBmag[5]*hamil[5]+BstarXdBmag[1]*hamil[1]))*rdx2+40.0*apardot[1]*q_))/m_; 
  alphaR[2] = -(0.03535533905932736*(((17.32050807568877*BstarYdBmag[20]+33.54101966249684*BstarYdBmag[6]+19.36491673103708*BstarYdBmag[1])*hamil[20]+8.660254037844386*BstarYdBmag[19]*hamil[19]+(17.32050807568877*BstarYdBmag[12]+33.54101966249685*BstarYdBmag[3]+19.36491673103709*BstarYdBmag[0])*hamil[12]+8.660254037844386*(BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2]))*rdy2+(7.745966692414834*BstarXdBmag[2]*hamil[20]+19.36491673103708*BstarXdBmag[1]*hamil[19]+19.36491673103709*BstarXdBmag[5]*hamil[11]+8.660254037844386*(BstarXdBmag[0]*hamil[5]+hamil[1]*BstarXdBmag[2]))*rdx2+40.0*apardot[2]*q_))/m_; 
  alphaR[3] = -(0.1767766952966368*(3.872983346207417*BstarXdBmag[1]*hamil[25]+1.732050807568877*BstarXdBmag[0]*hamil[8])*rdx2)/m_; 
  alphaR[4] = -(0.03535533905932736*((hamil[20]*(30.0*BstarYdBmag[21]+17.32050807568877*(BstarYdBmag[12]+BstarYdBmag[11])+33.54101966249684*BstarYdBmag[3]+19.36491673103708*BstarYdBmag[0])+17.32050807568877*hamil[12]*BstarYdBmag[20]+7.745966692414834*(BstarYdBmag[5]*hamil[19]+hamil[5]*BstarYdBmag[19])+(33.54101966249685*BstarYdBmag[6]+19.36491673103709*BstarYdBmag[1])*hamil[12]+8.660254037844386*(BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5]))*rdy2+(7.745966692414834*BstarXdBmag[5]*hamil[20]+(17.32050807568877*BstarXdBmag[11]+19.36491673103708*BstarXdBmag[0])*hamil[19]+hamil[11]*(17.32050807568877*BstarXdBmag[19]+19.36491673103709*BstarXdBmag[2])+8.660254037844386*(BstarXdBmag[1]*hamil[5]+hamil[1]*BstarXdBmag[5]))*rdx2+40.0*apardot[3]*q_))/m_; 
  alphaR[5] = -(0.03535533905932736*((17.32050807568877*BstarXdBmag[11]+19.36491673103708*BstarXdBmag[0])*hamil[25]+8.660254037844386*BstarXdBmag[1]*hamil[8])*rdx2)/m_; 
  alphaR[6] = -(0.1767766952966368*(3.872983346207417*BstarXdBmag[5]*hamil[25]+1.732050807568877*BstarXdBmag[2]*hamil[8])*rdx2)/m_; 
  alphaR[7] = -(0.005050762722761052*(((67.0820393249937*hamil[19]+105.0*hamil[2])*BstarYdBmag[21]+121.2435565298214*BstarYdBmag[5]*hamil[20]+(38.72983346207417*BstarYdBmag[11]+105.0*BstarYdBmag[3]+60.62177826491071*BstarYdBmag[0])*hamil[19]+135.5544171172596*hamil[12]*BstarYdBmag[19]+60.6217782649107*hamil[2]*BstarYdBmag[11]+hamil[5]*(93.91485505499116*BstarYdBmag[6]+54.22176684690384*BstarYdBmag[1]))*rdy2+(121.2435565298214*BstarXdBmag[5]*hamil[19]+60.62177826491071*hamil[5]*BstarXdBmag[19]+121.2435565298214*BstarXdBmag[1]*hamil[11]+60.6217782649107*hamil[1]*BstarXdBmag[11])*rdx2+280.0*apardot[4]*q_))/m_; 
  alphaR[8] = -(0.07905694150420944*((7.745966692414834*BstarYdBmag[5]*hamil[20]+3.872983346207417*hamil[5]*BstarYdBmag[20]+7.745966692414834*BstarYdBmag[2]*hamil[12]+3.872983346207417*hamil[2]*BstarYdBmag[12])*rdy2+(3.872983346207417*BstarXdBmag[0]*hamil[20]+7.745966692414834*BstarXdBmag[5]*hamil[19]+3.464101615137754*BstarXdBmag[2]*hamil[5])*rdx2+17.88854381999832*apardot[5]*q_))/m_; 
  alphaR[10] = -(0.1767766952966368*((3.464101615137754*BstarXdBmag[19]+3.872983346207417*BstarXdBmag[2])*hamil[25]+1.732050807568877*BstarXdBmag[5]*hamil[8])*rdx2)/m_; 
  alphaR[11] = -(0.005050762722761052*((234.787137637478*hamil[12]*BstarYdBmag[21]+(108.4435336938077*BstarYdBmag[20]+210.0*BstarYdBmag[6]+121.2435565298214*BstarYdBmag[1])*hamil[20]+(38.72983346207417*BstarYdBmag[19]+60.6217782649107*BstarYdBmag[2])*hamil[19]+60.6217782649107*hamil[2]*BstarYdBmag[19]+135.5544171172596*BstarYdBmag[11]*hamil[12]+54.22176684690384*BstarYdBmag[5]*hamil[5])*rdy2+(54.22176684690384*BstarXdBmag[19]*hamil[20]+121.2435565298214*BstarXdBmag[1]*hamil[19]+60.6217782649107*hamil[1]*BstarXdBmag[19]+121.2435565298214*BstarXdBmag[5]*hamil[11]+60.62177826491071*hamil[5]*BstarXdBmag[11])*rdx2+280.0*apardot[6]*q_))/m_; 
  alphaR[12] = -(0.04564354645876382*(((12.0*BstarYdBmag[19]+13.41640786499874*BstarYdBmag[2])*hamil[20]+(6.0*hamil[19]+6.708203932499369*hamil[2])*BstarYdBmag[20]+13.41640786499874*BstarYdBmag[5]*hamil[12]+6.708203932499369*hamil[5]*BstarYdBmag[12])*rdy2+(6.708203932499369*BstarXdBmag[1]*hamil[20]+(12.0*BstarXdBmag[19]+13.41640786499874*BstarXdBmag[2])*hamil[19]+6.0*BstarXdBmag[5]*hamil[5])*rdx2+30.98386676965934*apardot[7]*q_))/m_; 
  alphaR[13] = -(0.03535533905932736*(17.32050807568877*BstarXdBmag[1]*hamil[25]+8.660254037844387*hamil[8]*BstarXdBmag[11])*rdx2)/m_; 
  alphaR[17] = -(0.03535533905932736*(17.32050807568877*BstarXdBmag[5]*hamil[25]+8.660254037844387*hamil[8]*BstarXdBmag[19])*rdx2)/m_; 

  double fUpOrdL[27] = {0.};
  if (alphaL[17]-0.7453559924999302*alphaL[13]-0.7453559924999297*alphaL[12]-0.7453559924999302*alphaL[11]-1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]+1.118033988749895*alphaL[6]+1.118033988749895*alphaL[5]+1.118033988749895*alphaL[4]-0.8333333333333334*alphaL[3]-0.8333333333333334*alphaL[2]-0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fc); 
  } 
  if ((-1.0*alphaL[12])-1.000000000000001*alphaL[11]+0.7453559924999301*alphaL[8]+0.7453559924999301*alphaL[7]+1.500000000000001*alphaL[4]-1.118033988749895*alphaL[2]-1.118033988749895*alphaL[1]+0.833333333333333*alphaL[0] > 0.) {
    fUpOrdL[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fc); 
  } 
  if ((-1.0*alphaL[17])+0.7453559924999302*alphaL[13]-0.7453559924999297*alphaL[12]-0.7453559924999302*alphaL[11]+1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]-1.118033988749895*alphaL[6]-1.118033988749895*alphaL[5]+1.118033988749895*alphaL[4]+0.8333333333333334*alphaL[3]-0.8333333333333334*alphaL[2]-0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fc); 
  } 
  if ((-1.0*alphaL[13])+1.249999999999999*alphaL[12]-0.9316949906249115*alphaL[8]+0.7453559924999296*alphaL[7]+1.5*alphaL[5]-1.118033988749894*alphaL[3]-1.118033988749894*alphaL[1]+0.8333333333333324*alphaL[0] > 0.) {
    fUpOrdL[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fc); 
  } 
  if (alphaL[12]-0.74535599249993*alphaL[8]+0.5962847939999443*alphaL[7]-0.8944271909999164*alphaL[1]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fc); 
  } 
  if (alphaL[13]+1.249999999999999*alphaL[12]-0.9316949906249115*alphaL[8]+0.7453559924999296*alphaL[7]-1.5*alphaL[5]+1.118033988749894*alphaL[3]-1.118033988749894*alphaL[1]+0.8333333333333324*alphaL[0] > 0.) {
    fUpOrdL[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fc); 
  } 
  if ((-1.0*alphaL[17])-0.7453559924999302*alphaL[13]-0.7453559924999297*alphaL[12]+0.7453559924999302*alphaL[11]+1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]-1.118033988749895*alphaL[6]+1.118033988749895*alphaL[5]-1.118033988749895*alphaL[4]-0.8333333333333334*alphaL[3]+0.8333333333333334*alphaL[2]-0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fc); 
  } 
  if ((-1.0*alphaL[12])+1.000000000000001*alphaL[11]+0.7453559924999301*alphaL[8]+0.7453559924999301*alphaL[7]-1.500000000000001*alphaL[4]+1.118033988749895*alphaL[2]-1.118033988749895*alphaL[1]+0.833333333333333*alphaL[0] > 0.) {
    fUpOrdL[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fc); 
  } 
  if (alphaL[17]+0.7453559924999302*alphaL[13]-0.7453559924999297*alphaL[12]+0.7453559924999302*alphaL[11]-1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]+1.118033988749895*alphaL[6]-1.118033988749895*alphaL[5]-1.118033988749895*alphaL[4]+0.8333333333333334*alphaL[3]+0.8333333333333334*alphaL[2]-0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fc); 
  } 
  if ((-1.0*alphaL[17])+0.7453559924999298*alphaL[13]+0.7453559924999298*alphaL[11]+0.4444444444444446*alphaL[8]-0.5555555555555556*alphaL[7]+0.8944271909999165*alphaL[6]-0.666666666666667*alphaL[3]-0.666666666666667*alphaL[2]+0.4969039949999532*alphaL[0] > 0.) {
    fUpOrdL[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fl); 
  } else { 
    fUpOrdL[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fc); 
  } 
  if (alphaL[11]+0.5962847939999443*alphaL[8]-0.74535599249993*alphaL[7]-0.8944271909999164*alphaL[2]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fl); 
  } else { 
    fUpOrdL[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fc); 
  } 
  if (alphaL[17]-0.7453559924999298*alphaL[13]+0.7453559924999298*alphaL[11]+0.4444444444444446*alphaL[8]-0.5555555555555556*alphaL[7]-0.8944271909999165*alphaL[6]+0.666666666666667*alphaL[3]-0.666666666666667*alphaL[2]+0.4969039949999532*alphaL[0] > 0.) {
    fUpOrdL[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fl); 
  } else { 
    fUpOrdL[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fc); 
  } 
  if (alphaL[13]-0.74535599249993*alphaL[8]-0.74535599249993*alphaL[7]-0.8944271909999164*alphaL[3]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fl); 
  } else { 
    fUpOrdL[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fc); 
  } 
  if ((-1.118033988749895*alphaL[8])-1.118033988749895*alphaL[7]+alphaL[0] > 0.) {
    fUpOrdL[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fl); 
  } else { 
    fUpOrdL[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fc); 
  } 
  if ((-1.0*alphaL[13])-0.74535599249993*alphaL[8]-0.74535599249993*alphaL[7]+0.8944271909999164*alphaL[3]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fl); 
  } else { 
    fUpOrdL[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fc); 
  } 
  if (alphaL[17]+0.7453559924999298*alphaL[13]-0.7453559924999298*alphaL[11]+0.4444444444444446*alphaL[8]-0.5555555555555556*alphaL[7]-0.8944271909999165*alphaL[6]-0.666666666666667*alphaL[3]+0.666666666666667*alphaL[2]+0.4969039949999532*alphaL[0] > 0.) {
    fUpOrdL[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fl); 
  } else { 
    fUpOrdL[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fc); 
  } 
  if ((-1.0*alphaL[11])+0.5962847939999443*alphaL[8]-0.74535599249993*alphaL[7]+0.8944271909999164*alphaL[2]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fl); 
  } else { 
    fUpOrdL[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fc); 
  } 
  if ((-1.0*alphaL[17])-0.7453559924999298*alphaL[13]-0.7453559924999298*alphaL[11]+0.4444444444444446*alphaL[8]-0.5555555555555556*alphaL[7]+0.8944271909999165*alphaL[6]+0.666666666666667*alphaL[3]+0.666666666666667*alphaL[2]+0.4969039949999532*alphaL[0] > 0.) {
    fUpOrdL[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fl); 
  } else { 
    fUpOrdL[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fc); 
  } 
  if (alphaL[17]-0.7453559924999302*alphaL[13]+0.7453559924999297*alphaL[12]-0.7453559924999302*alphaL[11]+1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]+1.118033988749895*alphaL[6]-1.118033988749895*alphaL[5]-1.118033988749895*alphaL[4]-0.8333333333333334*alphaL[3]-0.8333333333333334*alphaL[2]+0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fl); 
  } else { 
    fUpOrdL[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fc); 
  } 
  if (alphaL[12]-1.000000000000001*alphaL[11]+0.7453559924999301*alphaL[8]+0.7453559924999301*alphaL[7]-1.500000000000001*alphaL[4]-1.118033988749895*alphaL[2]+1.118033988749895*alphaL[1]+0.833333333333333*alphaL[0] > 0.) {
    fUpOrdL[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fl); 
  } else { 
    fUpOrdL[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fc); 
  } 
  if ((-1.0*alphaL[17])+0.7453559924999302*alphaL[13]+0.7453559924999297*alphaL[12]-0.7453559924999302*alphaL[11]-1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]-1.118033988749895*alphaL[6]+1.118033988749895*alphaL[5]-1.118033988749895*alphaL[4]+0.8333333333333334*alphaL[3]-0.8333333333333334*alphaL[2]+0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fl); 
  } else { 
    fUpOrdL[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fc); 
  } 
  if ((-1.0*alphaL[13])-1.249999999999999*alphaL[12]-0.9316949906249115*alphaL[8]+0.7453559924999296*alphaL[7]-1.5*alphaL[5]-1.118033988749894*alphaL[3]+1.118033988749894*alphaL[1]+0.8333333333333324*alphaL[0] > 0.) {
    fUpOrdL[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fl); 
  } else { 
    fUpOrdL[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fc); 
  } 
  if ((-1.0*alphaL[12])-0.74535599249993*alphaL[8]+0.5962847939999443*alphaL[7]+0.8944271909999164*alphaL[1]+0.6666666666666666*alphaL[0] > 0.) {
    fUpOrdL[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fl); 
  } else { 
    fUpOrdL[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fc); 
  } 
  if (alphaL[13]-1.249999999999999*alphaL[12]-0.9316949906249115*alphaL[8]+0.7453559924999296*alphaL[7]+1.5*alphaL[5]+1.118033988749894*alphaL[3]+1.118033988749894*alphaL[1]+0.8333333333333324*alphaL[0] > 0.) {
    fUpOrdL[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fl); 
  } else { 
    fUpOrdL[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fc); 
  } 
  if ((-1.0*alphaL[17])-0.7453559924999302*alphaL[13]+0.7453559924999297*alphaL[12]+0.7453559924999302*alphaL[11]-1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]-1.118033988749895*alphaL[6]-1.118033988749895*alphaL[5]+1.118033988749895*alphaL[4]-0.8333333333333334*alphaL[3]+0.8333333333333334*alphaL[2]+0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fl); 
  } else { 
    fUpOrdL[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fc); 
  } 
  if (alphaL[12]+1.000000000000001*alphaL[11]+0.7453559924999301*alphaL[8]+0.7453559924999301*alphaL[7]+1.500000000000001*alphaL[4]+1.118033988749895*alphaL[2]+1.118033988749895*alphaL[1]+0.833333333333333*alphaL[0] > 0.) {
    fUpOrdL[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fl); 
  } else { 
    fUpOrdL[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fc); 
  } 
  if (alphaL[17]+0.7453559924999302*alphaL[13]+0.7453559924999297*alphaL[12]+0.7453559924999302*alphaL[11]+1.499999999999999*alphaL[10]+0.5555555555555556*alphaL[8]+0.5555555555555556*alphaL[7]+1.118033988749895*alphaL[6]+1.118033988749895*alphaL[5]+1.118033988749895*alphaL[4]+0.8333333333333334*alphaL[3]+0.8333333333333334*alphaL[2]+0.8333333333333334*alphaL[1]+0.6211299937499412*alphaL[0] > 0.) {
    fUpOrdL[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fl); 
  } else { 
    fUpOrdL[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fc); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[48] = {0.}; 
  GhatL[0] = 0.3535533905932737*alphaL[17]*fUpL[17]+0.3535533905932737*alphaL[13]*fUpL[13]+0.3535533905932737*alphaL[12]*fUpL[12]+0.3535533905932737*alphaL[11]*fUpL[11]+0.3535533905932737*alphaL[10]*fUpL[10]+0.3535533905932737*alphaL[8]*fUpL[8]+0.3535533905932737*alphaL[7]*fUpL[7]+0.3535533905932737*alphaL[6]*fUpL[6]+0.3535533905932737*alphaL[5]*fUpL[5]+0.3535533905932737*alphaL[4]*fUpL[4]+0.3535533905932737*alphaL[3]*fUpL[3]+0.3535533905932737*alphaL[2]*fUpL[2]+0.3535533905932737*alphaL[1]*fUpL[1]+0.3535533905932737*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.3162277660168379*alphaL[10]*fUpL[17]+0.3162277660168379*fUpL[10]*alphaL[17]+0.3162277660168379*alphaL[5]*fUpL[13]+0.3162277660168379*fUpL[5]*alphaL[13]+0.3535533905932737*alphaL[8]*fUpL[12]+0.3535533905932737*fUpL[8]*alphaL[12]+0.3162277660168379*alphaL[4]*fUpL[11]+0.3162277660168379*fUpL[4]*alphaL[11]+0.3535533905932737*alphaL[6]*fUpL[10]+0.3535533905932737*fUpL[6]*alphaL[10]+0.3162277660168379*alphaL[1]*fUpL[7]+0.3162277660168379*fUpL[1]*alphaL[7]+0.3535533905932737*alphaL[3]*fUpL[5]+0.3535533905932737*fUpL[3]*alphaL[5]+0.3535533905932737*alphaL[2]*fUpL[4]+0.3535533905932737*fUpL[2]*alphaL[4]+0.3535533905932737*alphaL[0]*fUpL[1]+0.3535533905932737*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.3162277660168379*alphaL[10]*fUpL[18]+0.3535533905932737*alphaL[13]*fUpL[17]+0.3535533905932737*fUpL[13]*alphaL[17]+0.3162277660168379*alphaL[6]*fUpL[14]+0.3162277660168379*alphaL[4]*fUpL[12]+0.3162277660168379*fUpL[4]*alphaL[12]+0.3535533905932737*alphaL[7]*fUpL[11]+0.3535533905932737*fUpL[7]*alphaL[11]+0.3535533905932737*alphaL[5]*fUpL[10]+0.3535533905932737*fUpL[5]*alphaL[10]+0.3162277660168379*alphaL[2]*fUpL[8]+0.3162277660168379*fUpL[2]*alphaL[8]+0.3535533905932737*alphaL[3]*fUpL[6]+0.3535533905932737*fUpL[3]*alphaL[6]+0.3535533905932737*alphaL[1]*fUpL[4]+0.3535533905932737*fUpL[1]*alphaL[4]+0.3535533905932737*alphaL[0]*fUpL[2]+0.3535533905932737*fUpL[0]*alphaL[2]; 
  GhatL[3] = 0.3162277660168379*alphaL[10]*fUpL[19]+0.3535533905932737*alphaL[12]*fUpL[18]+0.3535533905932737*alphaL[11]*fUpL[17]+0.3535533905932737*fUpL[11]*alphaL[17]+0.3162277660168379*alphaL[6]*fUpL[16]+0.3162277660168379*alphaL[5]*fUpL[15]+0.3535533905932737*alphaL[8]*fUpL[14]+0.3535533905932737*alphaL[7]*fUpL[13]+0.3535533905932737*fUpL[7]*alphaL[13]+0.3535533905932737*alphaL[4]*fUpL[10]+0.3535533905932737*fUpL[4]*alphaL[10]+0.3162277660168379*alphaL[3]*fUpL[9]+0.3535533905932737*alphaL[2]*fUpL[6]+0.3535533905932737*fUpL[2]*alphaL[6]+0.3535533905932737*alphaL[1]*fUpL[5]+0.3535533905932737*fUpL[1]*alphaL[5]+0.3535533905932737*alphaL[0]*fUpL[3]+0.3535533905932737*fUpL[0]*alphaL[3]; 
  GhatL[4] = 0.2828427124746191*alphaL[17]*fUpL[18]+0.3162277660168379*alphaL[6]*fUpL[18]+0.3162277660168379*alphaL[5]*fUpL[17]+0.3162277660168379*fUpL[5]*alphaL[17]+0.3162277660168379*alphaL[10]*fUpL[14]+0.3162277660168379*alphaL[10]*fUpL[13]+0.3162277660168379*fUpL[10]*alphaL[13]+0.2828427124746191*alphaL[11]*fUpL[12]+0.3162277660168379*alphaL[2]*fUpL[12]+0.2828427124746191*fUpL[11]*alphaL[12]+0.3162277660168379*fUpL[2]*alphaL[12]+0.3162277660168379*alphaL[1]*fUpL[11]+0.3162277660168379*fUpL[1]*alphaL[11]+0.3535533905932737*alphaL[3]*fUpL[10]+0.3535533905932737*fUpL[3]*alphaL[10]+0.3162277660168379*alphaL[4]*fUpL[8]+0.3162277660168379*fUpL[4]*alphaL[8]+0.3162277660168379*alphaL[4]*fUpL[7]+0.3162277660168379*fUpL[4]*alphaL[7]+0.3535533905932737*alphaL[5]*fUpL[6]+0.3535533905932737*fUpL[5]*alphaL[6]+0.3535533905932737*alphaL[0]*fUpL[4]+0.3535533905932737*fUpL[0]*alphaL[4]+0.3535533905932737*alphaL[1]*fUpL[2]+0.3535533905932737*fUpL[1]*alphaL[2]; 
  GhatL[5] = 0.2828427124746191*alphaL[17]*fUpL[19]+0.3162277660168379*alphaL[6]*fUpL[19]+0.3535533905932737*alphaL[8]*fUpL[18]+0.3162277660168379*alphaL[4]*fUpL[17]+0.3162277660168379*fUpL[4]*alphaL[17]+0.3162277660168379*alphaL[10]*fUpL[16]+0.2828427124746191*alphaL[13]*fUpL[15]+0.3162277660168379*alphaL[3]*fUpL[15]+0.3535533905932737*alphaL[12]*fUpL[14]+0.3162277660168379*alphaL[1]*fUpL[13]+0.3162277660168379*fUpL[1]*alphaL[13]+0.3162277660168379*alphaL[10]*fUpL[11]+0.3162277660168379*fUpL[10]*alphaL[11]+0.3535533905932737*alphaL[2]*fUpL[10]+0.3535533905932737*fUpL[2]*alphaL[10]+0.3162277660168379*alphaL[5]*fUpL[9]+0.3162277660168379*alphaL[5]*fUpL[7]+0.3162277660168379*fUpL[5]*alphaL[7]+0.3535533905932737*alphaL[4]*fUpL[6]+0.3535533905932737*fUpL[4]*alphaL[6]+0.3535533905932737*alphaL[0]*fUpL[5]+0.3535533905932737*fUpL[0]*alphaL[5]+0.3535533905932737*alphaL[1]*fUpL[3]+0.3535533905932737*fUpL[1]*alphaL[3]; 
  GhatL[6] = 0.3162277660168379*alphaL[5]*fUpL[19]+0.3162277660168379*alphaL[4]*fUpL[18]+0.3535533905932737*alphaL[7]*fUpL[17]+0.3535533905932737*fUpL[7]*alphaL[17]+0.3162277660168379*alphaL[3]*fUpL[16]+0.3162277660168379*alphaL[10]*fUpL[15]+0.3162277660168379*alphaL[2]*fUpL[14]+0.3535533905932737*alphaL[11]*fUpL[13]+0.3535533905932737*fUpL[11]*alphaL[13]+0.3162277660168379*alphaL[10]*fUpL[12]+0.3162277660168379*fUpL[10]*alphaL[12]+0.3535533905932737*alphaL[1]*fUpL[10]+0.3535533905932737*fUpL[1]*alphaL[10]+0.3162277660168379*alphaL[6]*fUpL[9]+0.3162277660168379*alphaL[6]*fUpL[8]+0.3162277660168379*fUpL[6]*alphaL[8]+0.3535533905932737*alphaL[0]*fUpL[6]+0.3535533905932737*fUpL[0]*alphaL[6]+0.3535533905932737*alphaL[4]*fUpL[5]+0.3535533905932737*fUpL[4]*alphaL[5]+0.3535533905932737*alphaL[2]*fUpL[3]+0.3535533905932737*fUpL[2]*alphaL[3]; 
  GhatL[7] = 0.2258769757263128*alphaL[17]*fUpL[17]+0.3535533905932737*alphaL[6]*fUpL[17]+0.3535533905932737*fUpL[6]*alphaL[17]+0.2258769757263128*alphaL[13]*fUpL[13]+0.3535533905932737*alphaL[3]*fUpL[13]+0.3535533905932737*fUpL[3]*alphaL[13]+0.3162277660168379*alphaL[12]*fUpL[12]+0.2258769757263128*alphaL[11]*fUpL[11]+0.3535533905932737*alphaL[2]*fUpL[11]+0.3535533905932737*fUpL[2]*alphaL[11]+0.3162277660168379*alphaL[10]*fUpL[10]+0.2258769757263128*alphaL[7]*fUpL[7]+0.3535533905932737*alphaL[0]*fUpL[7]+0.3535533905932737*fUpL[0]*alphaL[7]+0.3162277660168379*alphaL[5]*fUpL[5]+0.3162277660168379*alphaL[4]*fUpL[4]+0.3162277660168379*alphaL[1]*fUpL[1]; 
  GhatL[8] = 0.3535533905932737*alphaL[5]*fUpL[18]+0.3162277660168379*alphaL[17]*fUpL[17]+0.3535533905932737*alphaL[3]*fUpL[14]+0.2258769757263128*alphaL[12]*fUpL[12]+0.3535533905932737*alphaL[1]*fUpL[12]+0.3535533905932737*fUpL[1]*alphaL[12]+0.3162277660168379*alphaL[11]*fUpL[11]+0.3162277660168379*alphaL[10]*fUpL[10]+0.2258769757263128*alphaL[8]*fUpL[8]+0.3535533905932737*alphaL[0]*fUpL[8]+0.3535533905932737*fUpL[0]*alphaL[8]+0.3162277660168379*alphaL[6]*fUpL[6]+0.3162277660168379*alphaL[4]*fUpL[4]+0.3162277660168379*alphaL[2]*fUpL[2]; 
  GhatL[9] = 0.3535533905932737*alphaL[4]*fUpL[19]+0.3162277660168379*alphaL[17]*fUpL[17]+0.3535533905932737*alphaL[2]*fUpL[16]+0.3535533905932737*alphaL[1]*fUpL[15]+0.3162277660168379*alphaL[13]*fUpL[13]+0.3162277660168379*alphaL[10]*fUpL[10]+0.3535533905932737*alphaL[0]*fUpL[9]+0.3162277660168379*alphaL[6]*fUpL[6]+0.3162277660168379*alphaL[5]*fUpL[5]+0.3162277660168379*alphaL[3]*fUpL[3]; 
  GhatL[10] = 0.282842712474619*alphaL[13]*fUpL[19]+0.3162277660168379*alphaL[3]*fUpL[19]+0.282842712474619*alphaL[11]*fUpL[18]+0.3162277660168379*alphaL[2]*fUpL[18]+0.282842712474619*alphaL[12]*fUpL[17]+0.3162277660168379*alphaL[1]*fUpL[17]+0.282842712474619*fUpL[15]*alphaL[17]+0.282842712474619*fUpL[12]*alphaL[17]+0.3162277660168379*fUpL[1]*alphaL[17]+0.3162277660168379*alphaL[5]*fUpL[16]+0.3162277660168379*alphaL[6]*fUpL[15]+0.3162277660168379*alphaL[4]*fUpL[14]+0.3162277660168379*alphaL[4]*fUpL[13]+0.3162277660168379*fUpL[4]*alphaL[13]+0.3162277660168379*alphaL[6]*fUpL[12]+0.3162277660168379*fUpL[6]*alphaL[12]+0.3162277660168379*alphaL[5]*fUpL[11]+0.3162277660168379*fUpL[5]*alphaL[11]+0.3162277660168379*alphaL[8]*fUpL[10]+0.3162277660168379*alphaL[7]*fUpL[10]+0.3535533905932737*alphaL[0]*fUpL[10]+0.3162277660168379*fUpL[9]*alphaL[10]+0.3162277660168379*fUpL[8]*alphaL[10]+0.3162277660168379*fUpL[7]*alphaL[10]+0.3535533905932737*fUpL[0]*alphaL[10]+0.3535533905932737*alphaL[1]*fUpL[6]+0.3535533905932737*fUpL[1]*alphaL[6]+0.3535533905932737*alphaL[2]*fUpL[5]+0.3535533905932737*fUpL[2]*alphaL[5]+0.3535533905932737*alphaL[3]*fUpL[4]+0.3535533905932737*fUpL[3]*alphaL[4]; 
  GhatL[11] = 0.282842712474619*alphaL[10]*fUpL[18]+0.2258769757263128*alphaL[13]*fUpL[17]+0.3535533905932737*alphaL[3]*fUpL[17]+0.3162277660168379*fUpL[14]*alphaL[17]+0.2258769757263128*fUpL[13]*alphaL[17]+0.3535533905932737*fUpL[3]*alphaL[17]+0.3535533905932737*alphaL[6]*fUpL[13]+0.3535533905932737*fUpL[6]*alphaL[13]+0.2828427124746191*alphaL[4]*fUpL[12]+0.2828427124746191*fUpL[4]*alphaL[12]+0.3162277660168379*alphaL[8]*fUpL[11]+0.2258769757263128*alphaL[7]*fUpL[11]+0.3535533905932737*alphaL[0]*fUpL[11]+0.3162277660168379*fUpL[8]*alphaL[11]+0.2258769757263128*fUpL[7]*alphaL[11]+0.3535533905932737*fUpL[0]*alphaL[11]+0.3162277660168379*alphaL[5]*fUpL[10]+0.3162277660168379*fUpL[5]*alphaL[10]+0.3535533905932737*alphaL[2]*fUpL[7]+0.3535533905932737*fUpL[2]*alphaL[7]+0.3162277660168379*alphaL[1]*fUpL[4]+0.3162277660168379*fUpL[1]*alphaL[4]; 
  GhatL[12] = 0.3162277660168379*alphaL[13]*fUpL[18]+0.3535533905932737*alphaL[3]*fUpL[18]+0.282842712474619*alphaL[10]*fUpL[17]+0.282842712474619*fUpL[10]*alphaL[17]+0.3535533905932737*alphaL[5]*fUpL[14]+0.2258769757263128*alphaL[8]*fUpL[12]+0.3162277660168379*alphaL[7]*fUpL[12]+0.3535533905932737*alphaL[0]*fUpL[12]+0.2258769757263128*fUpL[8]*alphaL[12]+0.3162277660168379*fUpL[7]*alphaL[12]+0.3535533905932737*fUpL[0]*alphaL[12]+0.2828427124746191*alphaL[4]*fUpL[11]+0.2828427124746191*fUpL[4]*alphaL[11]+0.3162277660168379*alphaL[6]*fUpL[10]+0.3162277660168379*fUpL[6]*alphaL[10]+0.3535533905932737*alphaL[1]*fUpL[8]+0.3535533905932737*fUpL[1]*alphaL[8]+0.3162277660168379*alphaL[2]*fUpL[4]+0.3162277660168379*fUpL[2]*alphaL[4]; 
  GhatL[13] = 0.282842712474619*alphaL[10]*fUpL[19]+0.3162277660168379*alphaL[12]*fUpL[18]+0.2258769757263128*alphaL[11]*fUpL[17]+0.3535533905932737*alphaL[2]*fUpL[17]+0.3162277660168379*fUpL[16]*alphaL[17]+0.2258769757263128*fUpL[11]*alphaL[17]+0.3535533905932737*fUpL[2]*alphaL[17]+0.2828427124746191*alphaL[5]*fUpL[15]+0.2258769757263128*alphaL[7]*fUpL[13]+0.3535533905932737*alphaL[0]*fUpL[13]+0.3162277660168379*fUpL[9]*alphaL[13]+0.2258769757263128*fUpL[7]*alphaL[13]+0.3535533905932737*fUpL[0]*alphaL[13]+0.3535533905932737*alphaL[6]*fUpL[11]+0.3535533905932737*fUpL[6]*alphaL[11]+0.3162277660168379*alphaL[4]*fUpL[10]+0.3162277660168379*fUpL[4]*alphaL[10]+0.3535533905932737*alphaL[3]*fUpL[7]+0.3535533905932737*fUpL[3]*alphaL[7]+0.3162277660168379*alphaL[1]*fUpL[5]+0.3162277660168379*fUpL[1]*alphaL[5]; 
  GhatL[14] = 0.282842712474619*alphaL[10]*fUpL[19]+0.2258769757263128*alphaL[12]*fUpL[18]+0.3535533905932737*alphaL[1]*fUpL[18]+0.3162277660168379*alphaL[11]*fUpL[17]+0.3162277660168379*fUpL[11]*alphaL[17]+0.2828427124746191*alphaL[6]*fUpL[16]+0.2258769757263128*alphaL[8]*fUpL[14]+0.3535533905932737*alphaL[0]*fUpL[14]+0.3535533905932737*alphaL[5]*fUpL[12]+0.3535533905932737*fUpL[5]*alphaL[12]+0.3162277660168379*alphaL[4]*fUpL[10]+0.3162277660168379*fUpL[4]*alphaL[10]+0.3535533905932737*alphaL[3]*fUpL[8]+0.3535533905932737*fUpL[3]*alphaL[8]+0.3162277660168379*alphaL[2]*fUpL[6]+0.3162277660168379*fUpL[2]*alphaL[6]; 
  GhatL[15] = 0.3162277660168379*alphaL[11]*fUpL[19]+0.3535533905932737*alphaL[2]*fUpL[19]+0.282842712474619*alphaL[10]*fUpL[17]+0.282842712474619*fUpL[10]*alphaL[17]+0.3535533905932737*alphaL[4]*fUpL[16]+0.3162277660168379*alphaL[7]*fUpL[15]+0.3535533905932737*alphaL[0]*fUpL[15]+0.2828427124746191*alphaL[5]*fUpL[13]+0.2828427124746191*fUpL[5]*alphaL[13]+0.3162277660168379*alphaL[6]*fUpL[10]+0.3162277660168379*fUpL[6]*alphaL[10]+0.3535533905932737*alphaL[1]*fUpL[9]+0.3162277660168379*alphaL[3]*fUpL[5]+0.3162277660168379*fUpL[3]*alphaL[5]; 
  GhatL[16] = 0.3162277660168379*alphaL[12]*fUpL[19]+0.3535533905932737*alphaL[1]*fUpL[19]+0.282842712474619*alphaL[10]*fUpL[18]+0.3162277660168379*alphaL[13]*fUpL[17]+0.3162277660168379*fUpL[13]*alphaL[17]+0.3162277660168379*alphaL[8]*fUpL[16]+0.3535533905932737*alphaL[0]*fUpL[16]+0.3535533905932737*alphaL[4]*fUpL[15]+0.2828427124746191*alphaL[6]*fUpL[14]+0.3162277660168379*alphaL[5]*fUpL[10]+0.3162277660168379*fUpL[5]*alphaL[10]+0.3535533905932737*alphaL[2]*fUpL[9]+0.3162277660168379*alphaL[3]*fUpL[6]+0.3162277660168379*fUpL[3]*alphaL[6]; 
  GhatL[17] = 0.2828427124746191*alphaL[5]*fUpL[19]+0.2828427124746191*alphaL[4]*fUpL[18]+0.3162277660168379*alphaL[8]*fUpL[17]+0.2258769757263128*alphaL[7]*fUpL[17]+0.3535533905932737*alphaL[0]*fUpL[17]+0.3162277660168379*fUpL[9]*alphaL[17]+0.3162277660168379*fUpL[8]*alphaL[17]+0.2258769757263128*fUpL[7]*alphaL[17]+0.3535533905932737*fUpL[0]*alphaL[17]+0.3162277660168379*alphaL[13]*fUpL[16]+0.282842712474619*alphaL[10]*fUpL[15]+0.3162277660168379*alphaL[11]*fUpL[14]+0.2258769757263128*alphaL[11]*fUpL[13]+0.3535533905932737*alphaL[2]*fUpL[13]+0.2258769757263128*fUpL[11]*alphaL[13]+0.3535533905932737*fUpL[2]*alphaL[13]+0.282842712474619*alphaL[10]*fUpL[12]+0.282842712474619*fUpL[10]*alphaL[12]+0.3535533905932737*alphaL[3]*fUpL[11]+0.3535533905932737*fUpL[3]*alphaL[11]+0.3162277660168379*alphaL[1]*fUpL[10]+0.3162277660168379*fUpL[1]*alphaL[10]+0.3535533905932737*alphaL[6]*fUpL[7]+0.3535533905932737*fUpL[6]*alphaL[7]+0.3162277660168379*alphaL[4]*fUpL[5]+0.3162277660168379*fUpL[4]*alphaL[5]; 
  GhatL[18] = 0.2529822128134704*alphaL[17]*fUpL[19]+0.2828427124746191*alphaL[6]*fUpL[19]+0.2258769757263128*alphaL[8]*fUpL[18]+0.3162277660168379*alphaL[7]*fUpL[18]+0.3535533905932737*alphaL[0]*fUpL[18]+0.2828427124746191*alphaL[4]*fUpL[17]+0.2828427124746191*fUpL[4]*alphaL[17]+0.282842712474619*alphaL[10]*fUpL[16]+0.2258769757263128*alphaL[12]*fUpL[14]+0.3535533905932737*alphaL[1]*fUpL[14]+0.3162277660168379*alphaL[12]*fUpL[13]+0.3162277660168379*fUpL[12]*alphaL[13]+0.3535533905932737*alphaL[3]*fUpL[12]+0.3535533905932737*fUpL[3]*alphaL[12]+0.282842712474619*alphaL[10]*fUpL[11]+0.282842712474619*fUpL[10]*alphaL[11]+0.3162277660168379*alphaL[2]*fUpL[10]+0.3162277660168379*fUpL[2]*alphaL[10]+0.3535533905932737*alphaL[5]*fUpL[8]+0.3535533905932737*fUpL[5]*alphaL[8]+0.3162277660168379*alphaL[4]*fUpL[6]+0.3162277660168379*fUpL[4]*alphaL[6]; 
  GhatL[19] = 0.3162277660168379*alphaL[8]*fUpL[19]+0.3162277660168379*alphaL[7]*fUpL[19]+0.3535533905932737*alphaL[0]*fUpL[19]+0.2529822128134704*alphaL[17]*fUpL[18]+0.2828427124746191*alphaL[6]*fUpL[18]+0.2828427124746191*alphaL[5]*fUpL[17]+0.2828427124746191*fUpL[5]*alphaL[17]+0.3162277660168379*alphaL[12]*fUpL[16]+0.3535533905932737*alphaL[1]*fUpL[16]+0.3162277660168379*alphaL[11]*fUpL[15]+0.3535533905932737*alphaL[2]*fUpL[15]+0.282842712474619*alphaL[10]*fUpL[14]+0.282842712474619*alphaL[10]*fUpL[13]+0.282842712474619*fUpL[10]*alphaL[13]+0.3162277660168379*alphaL[3]*fUpL[10]+0.3162277660168379*fUpL[3]*alphaL[10]+0.3535533905932737*alphaL[4]*fUpL[9]+0.3162277660168379*alphaL[5]*fUpL[6]+0.3162277660168379*fUpL[5]*alphaL[6]; 

  double fUpOrdR[27] = {0.};
  if (alphaR[17]-0.7453559924999302*alphaR[13]-0.7453559924999297*alphaR[12]-0.7453559924999302*alphaR[11]-1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]+1.118033988749895*alphaR[6]+1.118033988749895*alphaR[5]+1.118033988749895*alphaR[4]-0.8333333333333334*alphaR[3]-0.8333333333333334*alphaR[2]-0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fr); 
  } 
  if ((-1.0*alphaR[12])-1.000000000000001*alphaR[11]+0.7453559924999301*alphaR[8]+0.7453559924999301*alphaR[7]+1.500000000000001*alphaR[4]-1.118033988749895*alphaR[2]-1.118033988749895*alphaR[1]+0.833333333333333*alphaR[0] > 0.) {
    fUpOrdR[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-1.0*alphaR[17])+0.7453559924999302*alphaR[13]-0.7453559924999297*alphaR[12]-0.7453559924999302*alphaR[11]+1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]-1.118033988749895*alphaR[6]-1.118033988749895*alphaR[5]+1.118033988749895*alphaR[4]+0.8333333333333334*alphaR[3]-0.8333333333333334*alphaR[2]-0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if ((-1.0*alphaR[13])+1.249999999999999*alphaR[12]-0.9316949906249115*alphaR[8]+0.7453559924999296*alphaR[7]+1.5*alphaR[5]-1.118033988749894*alphaR[3]-1.118033988749894*alphaR[1]+0.8333333333333324*alphaR[0] > 0.) {
    fUpOrdR[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (alphaR[12]-0.74535599249993*alphaR[8]+0.5962847939999443*alphaR[7]-0.8944271909999164*alphaR[1]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (alphaR[13]+1.249999999999999*alphaR[12]-0.9316949906249115*alphaR[8]+0.7453559924999296*alphaR[7]-1.5*alphaR[5]+1.118033988749894*alphaR[3]-1.118033988749894*alphaR[1]+0.8333333333333324*alphaR[0] > 0.) {
    fUpOrdR[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if ((-1.0*alphaR[17])-0.7453559924999302*alphaR[13]-0.7453559924999297*alphaR[12]+0.7453559924999302*alphaR[11]+1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]-1.118033988749895*alphaR[6]+1.118033988749895*alphaR[5]-1.118033988749895*alphaR[4]-0.8333333333333334*alphaR[3]+0.8333333333333334*alphaR[2]-0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fr); 
  } 
  if ((-1.0*alphaR[12])+1.000000000000001*alphaR[11]+0.7453559924999301*alphaR[8]+0.7453559924999301*alphaR[7]-1.500000000000001*alphaR[4]+1.118033988749895*alphaR[2]-1.118033988749895*alphaR[1]+0.833333333333333*alphaR[0] > 0.) {
    fUpOrdR[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (alphaR[17]+0.7453559924999302*alphaR[13]-0.7453559924999297*alphaR[12]+0.7453559924999302*alphaR[11]-1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]+1.118033988749895*alphaR[6]-1.118033988749895*alphaR[5]-1.118033988749895*alphaR[4]+0.8333333333333334*alphaR[3]+0.8333333333333334*alphaR[2]-0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 
  if ((-1.0*alphaR[17])+0.7453559924999298*alphaR[13]+0.7453559924999298*alphaR[11]+0.4444444444444446*alphaR[8]-0.5555555555555556*alphaR[7]+0.8944271909999165*alphaR[6]-0.666666666666667*alphaR[3]-0.666666666666667*alphaR[2]+0.4969039949999532*alphaR[0] > 0.) {
    fUpOrdR[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpOrdR[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fr); 
  } 
  if (alphaR[11]+0.5962847939999443*alphaR[8]-0.74535599249993*alphaR[7]-0.8944271909999164*alphaR[2]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpOrdR[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fr); 
  } 
  if (alphaR[17]-0.7453559924999298*alphaR[13]+0.7453559924999298*alphaR[11]+0.4444444444444446*alphaR[8]-0.5555555555555556*alphaR[7]-0.8944271909999165*alphaR[6]+0.666666666666667*alphaR[3]-0.666666666666667*alphaR[2]+0.4969039949999532*alphaR[0] > 0.) {
    fUpOrdR[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpOrdR[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (alphaR[13]-0.74535599249993*alphaR[8]-0.74535599249993*alphaR[7]-0.8944271909999164*alphaR[3]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpOrdR[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fr); 
  } 
  if ((-1.118033988749895*alphaR[8])-1.118033988749895*alphaR[7]+alphaR[0] > 0.) {
    fUpOrdR[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpOrdR[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fr); 
  } 
  if ((-1.0*alphaR[13])-0.74535599249993*alphaR[8]-0.74535599249993*alphaR[7]+0.8944271909999164*alphaR[3]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpOrdR[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (alphaR[17]+0.7453559924999298*alphaR[13]-0.7453559924999298*alphaR[11]+0.4444444444444446*alphaR[8]-0.5555555555555556*alphaR[7]-0.8944271909999165*alphaR[6]-0.666666666666667*alphaR[3]+0.666666666666667*alphaR[2]+0.4969039949999532*alphaR[0] > 0.) {
    fUpOrdR[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpOrdR[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fr); 
  } 
  if ((-1.0*alphaR[11])+0.5962847939999443*alphaR[8]-0.74535599249993*alphaR[7]+0.8944271909999164*alphaR[2]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fc); 
  } else { 
    fUpOrdR[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fr); 
  } 
  if ((-1.0*alphaR[17])-0.7453559924999298*alphaR[13]-0.7453559924999298*alphaR[11]+0.4444444444444446*alphaR[8]-0.5555555555555556*alphaR[7]+0.8944271909999165*alphaR[6]+0.666666666666667*alphaR[3]+0.666666666666667*alphaR[2]+0.4969039949999532*alphaR[0] > 0.) {
    fUpOrdR[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpOrdR[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fr); 
  } 
  if (alphaR[17]-0.7453559924999302*alphaR[13]+0.7453559924999297*alphaR[12]-0.7453559924999302*alphaR[11]+1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]+1.118033988749895*alphaR[6]-1.118033988749895*alphaR[5]-1.118033988749895*alphaR[4]-0.8333333333333334*alphaR[3]-0.8333333333333334*alphaR[2]+0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fc); 
  } else { 
    fUpOrdR[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fr); 
  } 
  if (alphaR[12]-1.000000000000001*alphaR[11]+0.7453559924999301*alphaR[8]+0.7453559924999301*alphaR[7]-1.500000000000001*alphaR[4]-1.118033988749895*alphaR[2]+1.118033988749895*alphaR[1]+0.833333333333333*alphaR[0] > 0.) {
    fUpOrdR[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fc); 
  } else { 
    fUpOrdR[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fr); 
  } 
  if ((-1.0*alphaR[17])+0.7453559924999302*alphaR[13]+0.7453559924999297*alphaR[12]-0.7453559924999302*alphaR[11]-1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]-1.118033988749895*alphaR[6]+1.118033988749895*alphaR[5]-1.118033988749895*alphaR[4]+0.8333333333333334*alphaR[3]-0.8333333333333334*alphaR[2]+0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fc); 
  } else { 
    fUpOrdR[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fr); 
  } 
  if ((-1.0*alphaR[13])-1.249999999999999*alphaR[12]-0.9316949906249115*alphaR[8]+0.7453559924999296*alphaR[7]-1.5*alphaR[5]-1.118033988749894*alphaR[3]+1.118033988749894*alphaR[1]+0.8333333333333324*alphaR[0] > 0.) {
    fUpOrdR[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fc); 
  } else { 
    fUpOrdR[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fr); 
  } 
  if ((-1.0*alphaR[12])-0.74535599249993*alphaR[8]+0.5962847939999443*alphaR[7]+0.8944271909999164*alphaR[1]+0.6666666666666666*alphaR[0] > 0.) {
    fUpOrdR[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fc); 
  } else { 
    fUpOrdR[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fr); 
  } 
  if (alphaR[13]-1.249999999999999*alphaR[12]-0.9316949906249115*alphaR[8]+0.7453559924999296*alphaR[7]+1.5*alphaR[5]+1.118033988749894*alphaR[3]+1.118033988749894*alphaR[1]+0.8333333333333324*alphaR[0] > 0.) {
    fUpOrdR[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fc); 
  } else { 
    fUpOrdR[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fr); 
  } 
  if ((-1.0*alphaR[17])-0.7453559924999302*alphaR[13]+0.7453559924999297*alphaR[12]+0.7453559924999302*alphaR[11]-1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]-1.118033988749895*alphaR[6]-1.118033988749895*alphaR[5]+1.118033988749895*alphaR[4]-0.8333333333333334*alphaR[3]+0.8333333333333334*alphaR[2]+0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fc); 
  } else { 
    fUpOrdR[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fr); 
  } 
  if (alphaR[12]+1.000000000000001*alphaR[11]+0.7453559924999301*alphaR[8]+0.7453559924999301*alphaR[7]+1.500000000000001*alphaR[4]+1.118033988749895*alphaR[2]+1.118033988749895*alphaR[1]+0.833333333333333*alphaR[0] > 0.) {
    fUpOrdR[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fc); 
  } else { 
    fUpOrdR[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fr); 
  } 
  if (alphaR[17]+0.7453559924999302*alphaR[13]+0.7453559924999297*alphaR[12]+0.7453559924999302*alphaR[11]+1.499999999999999*alphaR[10]+0.5555555555555556*alphaR[8]+0.5555555555555556*alphaR[7]+1.118033988749895*alphaR[6]+1.118033988749895*alphaR[5]+1.118033988749895*alphaR[4]+0.8333333333333334*alphaR[3]+0.8333333333333334*alphaR[2]+0.8333333333333334*alphaR[1]+0.6211299937499412*alphaR[0] > 0.) {
    fUpOrdR[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fc); 
  } else { 
    fUpOrdR[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[48] = {0.}; 
  GhatR[0] = 0.3535533905932737*alphaR[17]*fUpR[17]+0.3535533905932737*alphaR[13]*fUpR[13]+0.3535533905932737*alphaR[12]*fUpR[12]+0.3535533905932737*alphaR[11]*fUpR[11]+0.3535533905932737*alphaR[10]*fUpR[10]+0.3535533905932737*alphaR[8]*fUpR[8]+0.3535533905932737*alphaR[7]*fUpR[7]+0.3535533905932737*alphaR[6]*fUpR[6]+0.3535533905932737*alphaR[5]*fUpR[5]+0.3535533905932737*alphaR[4]*fUpR[4]+0.3535533905932737*alphaR[3]*fUpR[3]+0.3535533905932737*alphaR[2]*fUpR[2]+0.3535533905932737*alphaR[1]*fUpR[1]+0.3535533905932737*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.3162277660168379*alphaR[10]*fUpR[17]+0.3162277660168379*fUpR[10]*alphaR[17]+0.3162277660168379*alphaR[5]*fUpR[13]+0.3162277660168379*fUpR[5]*alphaR[13]+0.3535533905932737*alphaR[8]*fUpR[12]+0.3535533905932737*fUpR[8]*alphaR[12]+0.3162277660168379*alphaR[4]*fUpR[11]+0.3162277660168379*fUpR[4]*alphaR[11]+0.3535533905932737*alphaR[6]*fUpR[10]+0.3535533905932737*fUpR[6]*alphaR[10]+0.3162277660168379*alphaR[1]*fUpR[7]+0.3162277660168379*fUpR[1]*alphaR[7]+0.3535533905932737*alphaR[3]*fUpR[5]+0.3535533905932737*fUpR[3]*alphaR[5]+0.3535533905932737*alphaR[2]*fUpR[4]+0.3535533905932737*fUpR[2]*alphaR[4]+0.3535533905932737*alphaR[0]*fUpR[1]+0.3535533905932737*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.3162277660168379*alphaR[10]*fUpR[18]+0.3535533905932737*alphaR[13]*fUpR[17]+0.3535533905932737*fUpR[13]*alphaR[17]+0.3162277660168379*alphaR[6]*fUpR[14]+0.3162277660168379*alphaR[4]*fUpR[12]+0.3162277660168379*fUpR[4]*alphaR[12]+0.3535533905932737*alphaR[7]*fUpR[11]+0.3535533905932737*fUpR[7]*alphaR[11]+0.3535533905932737*alphaR[5]*fUpR[10]+0.3535533905932737*fUpR[5]*alphaR[10]+0.3162277660168379*alphaR[2]*fUpR[8]+0.3162277660168379*fUpR[2]*alphaR[8]+0.3535533905932737*alphaR[3]*fUpR[6]+0.3535533905932737*fUpR[3]*alphaR[6]+0.3535533905932737*alphaR[1]*fUpR[4]+0.3535533905932737*fUpR[1]*alphaR[4]+0.3535533905932737*alphaR[0]*fUpR[2]+0.3535533905932737*fUpR[0]*alphaR[2]; 
  GhatR[3] = 0.3162277660168379*alphaR[10]*fUpR[19]+0.3535533905932737*alphaR[12]*fUpR[18]+0.3535533905932737*alphaR[11]*fUpR[17]+0.3535533905932737*fUpR[11]*alphaR[17]+0.3162277660168379*alphaR[6]*fUpR[16]+0.3162277660168379*alphaR[5]*fUpR[15]+0.3535533905932737*alphaR[8]*fUpR[14]+0.3535533905932737*alphaR[7]*fUpR[13]+0.3535533905932737*fUpR[7]*alphaR[13]+0.3535533905932737*alphaR[4]*fUpR[10]+0.3535533905932737*fUpR[4]*alphaR[10]+0.3162277660168379*alphaR[3]*fUpR[9]+0.3535533905932737*alphaR[2]*fUpR[6]+0.3535533905932737*fUpR[2]*alphaR[6]+0.3535533905932737*alphaR[1]*fUpR[5]+0.3535533905932737*fUpR[1]*alphaR[5]+0.3535533905932737*alphaR[0]*fUpR[3]+0.3535533905932737*fUpR[0]*alphaR[3]; 
  GhatR[4] = 0.2828427124746191*alphaR[17]*fUpR[18]+0.3162277660168379*alphaR[6]*fUpR[18]+0.3162277660168379*alphaR[5]*fUpR[17]+0.3162277660168379*fUpR[5]*alphaR[17]+0.3162277660168379*alphaR[10]*fUpR[14]+0.3162277660168379*alphaR[10]*fUpR[13]+0.3162277660168379*fUpR[10]*alphaR[13]+0.2828427124746191*alphaR[11]*fUpR[12]+0.3162277660168379*alphaR[2]*fUpR[12]+0.2828427124746191*fUpR[11]*alphaR[12]+0.3162277660168379*fUpR[2]*alphaR[12]+0.3162277660168379*alphaR[1]*fUpR[11]+0.3162277660168379*fUpR[1]*alphaR[11]+0.3535533905932737*alphaR[3]*fUpR[10]+0.3535533905932737*fUpR[3]*alphaR[10]+0.3162277660168379*alphaR[4]*fUpR[8]+0.3162277660168379*fUpR[4]*alphaR[8]+0.3162277660168379*alphaR[4]*fUpR[7]+0.3162277660168379*fUpR[4]*alphaR[7]+0.3535533905932737*alphaR[5]*fUpR[6]+0.3535533905932737*fUpR[5]*alphaR[6]+0.3535533905932737*alphaR[0]*fUpR[4]+0.3535533905932737*fUpR[0]*alphaR[4]+0.3535533905932737*alphaR[1]*fUpR[2]+0.3535533905932737*fUpR[1]*alphaR[2]; 
  GhatR[5] = 0.2828427124746191*alphaR[17]*fUpR[19]+0.3162277660168379*alphaR[6]*fUpR[19]+0.3535533905932737*alphaR[8]*fUpR[18]+0.3162277660168379*alphaR[4]*fUpR[17]+0.3162277660168379*fUpR[4]*alphaR[17]+0.3162277660168379*alphaR[10]*fUpR[16]+0.2828427124746191*alphaR[13]*fUpR[15]+0.3162277660168379*alphaR[3]*fUpR[15]+0.3535533905932737*alphaR[12]*fUpR[14]+0.3162277660168379*alphaR[1]*fUpR[13]+0.3162277660168379*fUpR[1]*alphaR[13]+0.3162277660168379*alphaR[10]*fUpR[11]+0.3162277660168379*fUpR[10]*alphaR[11]+0.3535533905932737*alphaR[2]*fUpR[10]+0.3535533905932737*fUpR[2]*alphaR[10]+0.3162277660168379*alphaR[5]*fUpR[9]+0.3162277660168379*alphaR[5]*fUpR[7]+0.3162277660168379*fUpR[5]*alphaR[7]+0.3535533905932737*alphaR[4]*fUpR[6]+0.3535533905932737*fUpR[4]*alphaR[6]+0.3535533905932737*alphaR[0]*fUpR[5]+0.3535533905932737*fUpR[0]*alphaR[5]+0.3535533905932737*alphaR[1]*fUpR[3]+0.3535533905932737*fUpR[1]*alphaR[3]; 
  GhatR[6] = 0.3162277660168379*alphaR[5]*fUpR[19]+0.3162277660168379*alphaR[4]*fUpR[18]+0.3535533905932737*alphaR[7]*fUpR[17]+0.3535533905932737*fUpR[7]*alphaR[17]+0.3162277660168379*alphaR[3]*fUpR[16]+0.3162277660168379*alphaR[10]*fUpR[15]+0.3162277660168379*alphaR[2]*fUpR[14]+0.3535533905932737*alphaR[11]*fUpR[13]+0.3535533905932737*fUpR[11]*alphaR[13]+0.3162277660168379*alphaR[10]*fUpR[12]+0.3162277660168379*fUpR[10]*alphaR[12]+0.3535533905932737*alphaR[1]*fUpR[10]+0.3535533905932737*fUpR[1]*alphaR[10]+0.3162277660168379*alphaR[6]*fUpR[9]+0.3162277660168379*alphaR[6]*fUpR[8]+0.3162277660168379*fUpR[6]*alphaR[8]+0.3535533905932737*alphaR[0]*fUpR[6]+0.3535533905932737*fUpR[0]*alphaR[6]+0.3535533905932737*alphaR[4]*fUpR[5]+0.3535533905932737*fUpR[4]*alphaR[5]+0.3535533905932737*alphaR[2]*fUpR[3]+0.3535533905932737*fUpR[2]*alphaR[3]; 
  GhatR[7] = 0.2258769757263128*alphaR[17]*fUpR[17]+0.3535533905932737*alphaR[6]*fUpR[17]+0.3535533905932737*fUpR[6]*alphaR[17]+0.2258769757263128*alphaR[13]*fUpR[13]+0.3535533905932737*alphaR[3]*fUpR[13]+0.3535533905932737*fUpR[3]*alphaR[13]+0.3162277660168379*alphaR[12]*fUpR[12]+0.2258769757263128*alphaR[11]*fUpR[11]+0.3535533905932737*alphaR[2]*fUpR[11]+0.3535533905932737*fUpR[2]*alphaR[11]+0.3162277660168379*alphaR[10]*fUpR[10]+0.2258769757263128*alphaR[7]*fUpR[7]+0.3535533905932737*alphaR[0]*fUpR[7]+0.3535533905932737*fUpR[0]*alphaR[7]+0.3162277660168379*alphaR[5]*fUpR[5]+0.3162277660168379*alphaR[4]*fUpR[4]+0.3162277660168379*alphaR[1]*fUpR[1]; 
  GhatR[8] = 0.3535533905932737*alphaR[5]*fUpR[18]+0.3162277660168379*alphaR[17]*fUpR[17]+0.3535533905932737*alphaR[3]*fUpR[14]+0.2258769757263128*alphaR[12]*fUpR[12]+0.3535533905932737*alphaR[1]*fUpR[12]+0.3535533905932737*fUpR[1]*alphaR[12]+0.3162277660168379*alphaR[11]*fUpR[11]+0.3162277660168379*alphaR[10]*fUpR[10]+0.2258769757263128*alphaR[8]*fUpR[8]+0.3535533905932737*alphaR[0]*fUpR[8]+0.3535533905932737*fUpR[0]*alphaR[8]+0.3162277660168379*alphaR[6]*fUpR[6]+0.3162277660168379*alphaR[4]*fUpR[4]+0.3162277660168379*alphaR[2]*fUpR[2]; 
  GhatR[9] = 0.3535533905932737*alphaR[4]*fUpR[19]+0.3162277660168379*alphaR[17]*fUpR[17]+0.3535533905932737*alphaR[2]*fUpR[16]+0.3535533905932737*alphaR[1]*fUpR[15]+0.3162277660168379*alphaR[13]*fUpR[13]+0.3162277660168379*alphaR[10]*fUpR[10]+0.3535533905932737*alphaR[0]*fUpR[9]+0.3162277660168379*alphaR[6]*fUpR[6]+0.3162277660168379*alphaR[5]*fUpR[5]+0.3162277660168379*alphaR[3]*fUpR[3]; 
  GhatR[10] = 0.282842712474619*alphaR[13]*fUpR[19]+0.3162277660168379*alphaR[3]*fUpR[19]+0.282842712474619*alphaR[11]*fUpR[18]+0.3162277660168379*alphaR[2]*fUpR[18]+0.282842712474619*alphaR[12]*fUpR[17]+0.3162277660168379*alphaR[1]*fUpR[17]+0.282842712474619*fUpR[15]*alphaR[17]+0.282842712474619*fUpR[12]*alphaR[17]+0.3162277660168379*fUpR[1]*alphaR[17]+0.3162277660168379*alphaR[5]*fUpR[16]+0.3162277660168379*alphaR[6]*fUpR[15]+0.3162277660168379*alphaR[4]*fUpR[14]+0.3162277660168379*alphaR[4]*fUpR[13]+0.3162277660168379*fUpR[4]*alphaR[13]+0.3162277660168379*alphaR[6]*fUpR[12]+0.3162277660168379*fUpR[6]*alphaR[12]+0.3162277660168379*alphaR[5]*fUpR[11]+0.3162277660168379*fUpR[5]*alphaR[11]+0.3162277660168379*alphaR[8]*fUpR[10]+0.3162277660168379*alphaR[7]*fUpR[10]+0.3535533905932737*alphaR[0]*fUpR[10]+0.3162277660168379*fUpR[9]*alphaR[10]+0.3162277660168379*fUpR[8]*alphaR[10]+0.3162277660168379*fUpR[7]*alphaR[10]+0.3535533905932737*fUpR[0]*alphaR[10]+0.3535533905932737*alphaR[1]*fUpR[6]+0.3535533905932737*fUpR[1]*alphaR[6]+0.3535533905932737*alphaR[2]*fUpR[5]+0.3535533905932737*fUpR[2]*alphaR[5]+0.3535533905932737*alphaR[3]*fUpR[4]+0.3535533905932737*fUpR[3]*alphaR[4]; 
  GhatR[11] = 0.282842712474619*alphaR[10]*fUpR[18]+0.2258769757263128*alphaR[13]*fUpR[17]+0.3535533905932737*alphaR[3]*fUpR[17]+0.3162277660168379*fUpR[14]*alphaR[17]+0.2258769757263128*fUpR[13]*alphaR[17]+0.3535533905932737*fUpR[3]*alphaR[17]+0.3535533905932737*alphaR[6]*fUpR[13]+0.3535533905932737*fUpR[6]*alphaR[13]+0.2828427124746191*alphaR[4]*fUpR[12]+0.2828427124746191*fUpR[4]*alphaR[12]+0.3162277660168379*alphaR[8]*fUpR[11]+0.2258769757263128*alphaR[7]*fUpR[11]+0.3535533905932737*alphaR[0]*fUpR[11]+0.3162277660168379*fUpR[8]*alphaR[11]+0.2258769757263128*fUpR[7]*alphaR[11]+0.3535533905932737*fUpR[0]*alphaR[11]+0.3162277660168379*alphaR[5]*fUpR[10]+0.3162277660168379*fUpR[5]*alphaR[10]+0.3535533905932737*alphaR[2]*fUpR[7]+0.3535533905932737*fUpR[2]*alphaR[7]+0.3162277660168379*alphaR[1]*fUpR[4]+0.3162277660168379*fUpR[1]*alphaR[4]; 
  GhatR[12] = 0.3162277660168379*alphaR[13]*fUpR[18]+0.3535533905932737*alphaR[3]*fUpR[18]+0.282842712474619*alphaR[10]*fUpR[17]+0.282842712474619*fUpR[10]*alphaR[17]+0.3535533905932737*alphaR[5]*fUpR[14]+0.2258769757263128*alphaR[8]*fUpR[12]+0.3162277660168379*alphaR[7]*fUpR[12]+0.3535533905932737*alphaR[0]*fUpR[12]+0.2258769757263128*fUpR[8]*alphaR[12]+0.3162277660168379*fUpR[7]*alphaR[12]+0.3535533905932737*fUpR[0]*alphaR[12]+0.2828427124746191*alphaR[4]*fUpR[11]+0.2828427124746191*fUpR[4]*alphaR[11]+0.3162277660168379*alphaR[6]*fUpR[10]+0.3162277660168379*fUpR[6]*alphaR[10]+0.3535533905932737*alphaR[1]*fUpR[8]+0.3535533905932737*fUpR[1]*alphaR[8]+0.3162277660168379*alphaR[2]*fUpR[4]+0.3162277660168379*fUpR[2]*alphaR[4]; 
  GhatR[13] = 0.282842712474619*alphaR[10]*fUpR[19]+0.3162277660168379*alphaR[12]*fUpR[18]+0.2258769757263128*alphaR[11]*fUpR[17]+0.3535533905932737*alphaR[2]*fUpR[17]+0.3162277660168379*fUpR[16]*alphaR[17]+0.2258769757263128*fUpR[11]*alphaR[17]+0.3535533905932737*fUpR[2]*alphaR[17]+0.2828427124746191*alphaR[5]*fUpR[15]+0.2258769757263128*alphaR[7]*fUpR[13]+0.3535533905932737*alphaR[0]*fUpR[13]+0.3162277660168379*fUpR[9]*alphaR[13]+0.2258769757263128*fUpR[7]*alphaR[13]+0.3535533905932737*fUpR[0]*alphaR[13]+0.3535533905932737*alphaR[6]*fUpR[11]+0.3535533905932737*fUpR[6]*alphaR[11]+0.3162277660168379*alphaR[4]*fUpR[10]+0.3162277660168379*fUpR[4]*alphaR[10]+0.3535533905932737*alphaR[3]*fUpR[7]+0.3535533905932737*fUpR[3]*alphaR[7]+0.3162277660168379*alphaR[1]*fUpR[5]+0.3162277660168379*fUpR[1]*alphaR[5]; 
  GhatR[14] = 0.282842712474619*alphaR[10]*fUpR[19]+0.2258769757263128*alphaR[12]*fUpR[18]+0.3535533905932737*alphaR[1]*fUpR[18]+0.3162277660168379*alphaR[11]*fUpR[17]+0.3162277660168379*fUpR[11]*alphaR[17]+0.2828427124746191*alphaR[6]*fUpR[16]+0.2258769757263128*alphaR[8]*fUpR[14]+0.3535533905932737*alphaR[0]*fUpR[14]+0.3535533905932737*alphaR[5]*fUpR[12]+0.3535533905932737*fUpR[5]*alphaR[12]+0.3162277660168379*alphaR[4]*fUpR[10]+0.3162277660168379*fUpR[4]*alphaR[10]+0.3535533905932737*alphaR[3]*fUpR[8]+0.3535533905932737*fUpR[3]*alphaR[8]+0.3162277660168379*alphaR[2]*fUpR[6]+0.3162277660168379*fUpR[2]*alphaR[6]; 
  GhatR[15] = 0.3162277660168379*alphaR[11]*fUpR[19]+0.3535533905932737*alphaR[2]*fUpR[19]+0.282842712474619*alphaR[10]*fUpR[17]+0.282842712474619*fUpR[10]*alphaR[17]+0.3535533905932737*alphaR[4]*fUpR[16]+0.3162277660168379*alphaR[7]*fUpR[15]+0.3535533905932737*alphaR[0]*fUpR[15]+0.2828427124746191*alphaR[5]*fUpR[13]+0.2828427124746191*fUpR[5]*alphaR[13]+0.3162277660168379*alphaR[6]*fUpR[10]+0.3162277660168379*fUpR[6]*alphaR[10]+0.3535533905932737*alphaR[1]*fUpR[9]+0.3162277660168379*alphaR[3]*fUpR[5]+0.3162277660168379*fUpR[3]*alphaR[5]; 
  GhatR[16] = 0.3162277660168379*alphaR[12]*fUpR[19]+0.3535533905932737*alphaR[1]*fUpR[19]+0.282842712474619*alphaR[10]*fUpR[18]+0.3162277660168379*alphaR[13]*fUpR[17]+0.3162277660168379*fUpR[13]*alphaR[17]+0.3162277660168379*alphaR[8]*fUpR[16]+0.3535533905932737*alphaR[0]*fUpR[16]+0.3535533905932737*alphaR[4]*fUpR[15]+0.2828427124746191*alphaR[6]*fUpR[14]+0.3162277660168379*alphaR[5]*fUpR[10]+0.3162277660168379*fUpR[5]*alphaR[10]+0.3535533905932737*alphaR[2]*fUpR[9]+0.3162277660168379*alphaR[3]*fUpR[6]+0.3162277660168379*fUpR[3]*alphaR[6]; 
  GhatR[17] = 0.2828427124746191*alphaR[5]*fUpR[19]+0.2828427124746191*alphaR[4]*fUpR[18]+0.3162277660168379*alphaR[8]*fUpR[17]+0.2258769757263128*alphaR[7]*fUpR[17]+0.3535533905932737*alphaR[0]*fUpR[17]+0.3162277660168379*fUpR[9]*alphaR[17]+0.3162277660168379*fUpR[8]*alphaR[17]+0.2258769757263128*fUpR[7]*alphaR[17]+0.3535533905932737*fUpR[0]*alphaR[17]+0.3162277660168379*alphaR[13]*fUpR[16]+0.282842712474619*alphaR[10]*fUpR[15]+0.3162277660168379*alphaR[11]*fUpR[14]+0.2258769757263128*alphaR[11]*fUpR[13]+0.3535533905932737*alphaR[2]*fUpR[13]+0.2258769757263128*fUpR[11]*alphaR[13]+0.3535533905932737*fUpR[2]*alphaR[13]+0.282842712474619*alphaR[10]*fUpR[12]+0.282842712474619*fUpR[10]*alphaR[12]+0.3535533905932737*alphaR[3]*fUpR[11]+0.3535533905932737*fUpR[3]*alphaR[11]+0.3162277660168379*alphaR[1]*fUpR[10]+0.3162277660168379*fUpR[1]*alphaR[10]+0.3535533905932737*alphaR[6]*fUpR[7]+0.3535533905932737*fUpR[6]*alphaR[7]+0.3162277660168379*alphaR[4]*fUpR[5]+0.3162277660168379*fUpR[4]*alphaR[5]; 
  GhatR[18] = 0.2529822128134704*alphaR[17]*fUpR[19]+0.2828427124746191*alphaR[6]*fUpR[19]+0.2258769757263128*alphaR[8]*fUpR[18]+0.3162277660168379*alphaR[7]*fUpR[18]+0.3535533905932737*alphaR[0]*fUpR[18]+0.2828427124746191*alphaR[4]*fUpR[17]+0.2828427124746191*fUpR[4]*alphaR[17]+0.282842712474619*alphaR[10]*fUpR[16]+0.2258769757263128*alphaR[12]*fUpR[14]+0.3535533905932737*alphaR[1]*fUpR[14]+0.3162277660168379*alphaR[12]*fUpR[13]+0.3162277660168379*fUpR[12]*alphaR[13]+0.3535533905932737*alphaR[3]*fUpR[12]+0.3535533905932737*fUpR[3]*alphaR[12]+0.282842712474619*alphaR[10]*fUpR[11]+0.282842712474619*fUpR[10]*alphaR[11]+0.3162277660168379*alphaR[2]*fUpR[10]+0.3162277660168379*fUpR[2]*alphaR[10]+0.3535533905932737*alphaR[5]*fUpR[8]+0.3535533905932737*fUpR[5]*alphaR[8]+0.3162277660168379*alphaR[4]*fUpR[6]+0.3162277660168379*fUpR[4]*alphaR[6]; 
  GhatR[19] = 0.3162277660168379*alphaR[8]*fUpR[19]+0.3162277660168379*alphaR[7]*fUpR[19]+0.3535533905932737*alphaR[0]*fUpR[19]+0.2529822128134704*alphaR[17]*fUpR[18]+0.2828427124746191*alphaR[6]*fUpR[18]+0.2828427124746191*alphaR[5]*fUpR[17]+0.2828427124746191*fUpR[5]*alphaR[17]+0.3162277660168379*alphaR[12]*fUpR[16]+0.3535533905932737*alphaR[1]*fUpR[16]+0.3162277660168379*alphaR[11]*fUpR[15]+0.3535533905932737*alphaR[2]*fUpR[15]+0.282842712474619*alphaR[10]*fUpR[14]+0.282842712474619*alphaR[10]*fUpR[13]+0.282842712474619*fUpR[10]*alphaR[13]+0.3162277660168379*alphaR[3]*fUpR[10]+0.3162277660168379*fUpR[3]*alphaR[10]+0.3535533905932737*alphaR[4]*fUpR[9]+0.3162277660168379*alphaR[5]*fUpR[6]+0.3162277660168379*fUpR[5]*alphaR[6]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdvpar2; 
  out[3] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdvpar2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdvpar2; 
  out[6] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[7] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdvpar2; 
  out[8] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdvpar2; 
  out[9] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdvpar2; 
  out[10] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdvpar2; 
  out[11] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdvpar2; 
  out[12] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdvpar2; 
  out[13] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[14] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdvpar2; 
  out[15] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdvpar2; 
  out[16] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdvpar2; 
  out[17] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdvpar2; 
  out[18] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdvpar2; 
  out[19] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdvpar2; 
  out[20] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdvpar2; 
  out[21] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdvpar2; 
  out[22] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdvpar2; 
  out[23] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 
  out[24] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdvpar2; 
  out[25] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdvpar2; 
  out[26] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdvpar2; 
  out[27] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdvpar2; 
  out[28] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdvpar2; 
  out[29] += (0.7071067811865475*GhatL[16]-0.7071067811865475*GhatR[16])*rdvpar2; 
  out[30] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdvpar2; 
  out[31] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdvpar2; 
  out[32] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdvpar2; 
  out[33] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdvpar2; 
  out[34] += (1.58113883008419*GhatL[4]-1.58113883008419*GhatR[4])*rdvpar2; 
  out[35] += (0.7071067811865475*GhatL[17]-0.7071067811865475*GhatR[17])*rdvpar2; 
  out[36] += (0.7071067811865475*GhatL[18]-0.7071067811865475*GhatR[18])*rdvpar2; 
  out[37] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdvpar2; 
  out[38] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdvpar2; 
  out[39] += (1.58113883008419*GhatL[5]-1.58113883008419*GhatR[5])*rdvpar2; 
  out[40] += (1.58113883008419*GhatL[6]-1.58113883008419*GhatR[6])*rdvpar2; 
  out[41] += (0.7071067811865475*GhatL[19]-0.7071067811865475*GhatR[19])*rdvpar2; 
  out[42] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdvpar2; 
  out[43] += ((-1.224744871391589*GhatR[16])-1.224744871391589*GhatL[16])*rdvpar2; 
  out[44] += ((-1.224744871391589*GhatR[17])-1.224744871391589*GhatL[17])*rdvpar2; 
  out[45] += ((-1.224744871391589*GhatR[18])-1.224744871391589*GhatL[18])*rdvpar2; 
  out[46] += (1.58113883008419*GhatL[10]-1.58113883008419*GhatR[10])*rdvpar2; 
  out[47] += ((-1.224744871391589*GhatR[19])-1.224744871391589*GhatL[19])*rdvpar2; 

} 
