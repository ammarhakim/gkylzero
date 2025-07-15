#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *phi,
    const double *dualcurlbhatoverB, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // fin: Distribution function.
  // out: output increment.

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];
  double rdmu2 = 2.0/dxv[3];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[2]*dxv[2];

  const double *bioverJB_x = &bioverJB[0];
  const double *bioverJB_y = &bioverJB[4];
  const double *bioverJB_z = &bioverJB[8];

  const double *dualcurlbhatoverB_x = &dualcurlbhatoverB[0];
  const double *dualcurlbhatoverB_y = &dualcurlbhatoverB[4];
  const double *dualcurlbhatoverB_z = &dualcurlbhatoverB[8];

  double hamil[24] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+1.414213562373095*(vmapSq[0]*m_+bmag[0]*vmap[2]); 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_+1.414213562373095*bmag[2]*vmap[2]; 
  hamil[3] = 1.414213562373095*vmapSq[1]*m_; 
  hamil[4] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 2.0*phi[3]*q_+1.414213562373095*vmap[2]*bmag[3]; 
  hamil[8] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[9] = 1.414213562373095*bmag[2]*vmap[3]; 
  hamil[12] = 1.414213562373095*bmag[3]*vmap[3]; 
  hamil[16] = 1.414213562373095*vmapSq[2]*m_; 

  double vmap2 = vmap[1]*vmap[1]; 

  double hamil2[2] = {0.}; 
  hamil2[0] = hamil[3]*hamil[3]; 
  hamil2[1] = hamil[16]*hamil[16]; 

  double alphax[24] = {0.}; 
  alphax[0] = (rdx2*((1.25*dualcurlbhatoverB_x[0]*hamil2[1])/vmap2+(0.25*dualcurlbhatoverB_x[0]*hamil2[0])/vmap2))/(m_*q_)+((0.8660254037844386*bioverJB_y[1]*hamil[5]+0.8660254037844386*bioverJB_y[0]*hamil[2])*rdx2*rdz2)/q_; 
  alphax[1] = (rdx2*((1.25*dualcurlbhatoverB_x[1]*hamil2[1])/vmap2+(0.25*hamil2[0]*dualcurlbhatoverB_x[1])/vmap2))/(m_*q_)+((0.8660254037844386*bioverJB_y[0]*hamil[5]+0.8660254037844386*bioverJB_y[1]*hamil[2])*rdx2*rdz2)/q_; 
  alphax[2] = (rdx2*((1.25*hamil2[1]*dualcurlbhatoverB_x[2])/vmap2+(0.25*hamil2[0]*dualcurlbhatoverB_x[2])/vmap2))/(m_*q_)+((0.8660254037844386*bioverJB_y[3]*hamil[5]+0.8660254037844386*bioverJB_y[2]*hamil[2])*rdx2*rdz2)/q_; 
  alphax[3] = (1.118033988749895*dualcurlbhatoverB_x[0]*hamil[3]*hamil[16]*rdx2)/(m_*q_*vmap2); 
  alphax[4] = ((0.8660254037844386*bioverJB_y[1]*hamil[12]+0.8660254037844386*bioverJB_y[0]*hamil[9])*rdx2*rdz2)/q_; 
  alphax[5] = (rdx2*((1.25*hamil2[1]*dualcurlbhatoverB_x[3])/vmap2+(0.25*hamil2[0]*dualcurlbhatoverB_x[3])/vmap2))/(m_*q_)+((0.8660254037844386*bioverJB_y[2]*hamil[5]+0.8660254037844386*hamil[2]*bioverJB_y[3])*rdx2*rdz2)/q_; 
  alphax[6] = (1.118033988749895*dualcurlbhatoverB_x[1]*hamil[3]*hamil[16]*rdx2)/(m_*q_*vmap2); 
  alphax[7] = (1.118033988749895*dualcurlbhatoverB_x[2]*hamil[3]*hamil[16]*rdx2)/(m_*q_*vmap2); 
  alphax[8] = ((0.8660254037844386*bioverJB_y[0]*hamil[12]+0.8660254037844386*bioverJB_y[1]*hamil[9])*rdx2*rdz2)/q_; 
  alphax[9] = ((0.8660254037844386*bioverJB_y[3]*hamil[12]+0.8660254037844386*bioverJB_y[2]*hamil[9])*rdx2*rdz2)/q_; 
  alphax[11] = (1.118033988749895*dualcurlbhatoverB_x[3]*hamil[3]*hamil[16]*rdx2)/(m_*q_*vmap2); 
  alphax[12] = ((0.8660254037844386*bioverJB_y[2]*hamil[12]+0.8660254037844386*bioverJB_y[3]*hamil[9])*rdx2*rdz2)/q_; 
  alphax[16] = (1.118033988749895*dualcurlbhatoverB_x[0]*hamil2[1]*rdx2)/(m_*q_*vmap2); 
  alphax[17] = (1.118033988749895*dualcurlbhatoverB_x[1]*hamil2[1]*rdx2)/(m_*q_*vmap2); 
  alphax[18] = (1.118033988749895*hamil2[1]*dualcurlbhatoverB_x[2]*rdx2)/(m_*q_*vmap2); 
  alphax[20] = (1.118033988749895*hamil2[1]*dualcurlbhatoverB_x[3]*rdx2)/(m_*q_*vmap2); 


  out[1] += 0.4330127018922193*(alphax[20]*fin[20]+alphax[18]*fin[18]+alphax[17]*fin[17]+alphax[16]*fin[16]+alphax[12]*fin[12]+alphax[11]*fin[11]+alphax[9]*fin[9]+alphax[8]*fin[8]+alphax[7]*fin[7]+alphax[6]*fin[6]+alphax[5]*fin[5]+alphax[4]*fin[4]+alphax[3]*fin[3]+alphax[2]*fin[2]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[5] += 0.4330127018922194*(alphax[17]*fin[20]+fin[17]*alphax[20]+alphax[16]*fin[18]+fin[16]*alphax[18])+0.4330127018922193*(alphax[8]*fin[12]+fin[8]*alphax[12]+alphax[6]*fin[11]+fin[6]*alphax[11]+alphax[4]*fin[9]+fin[4]*alphax[9]+alphax[3]*fin[7]+fin[3]*alphax[7]+alphax[1]*fin[5]+fin[1]*alphax[5]+alphax[0]*fin[2]+fin[0]*alphax[2]); 
  out[6] += 0.3872983346207416*(alphax[11]*fin[20]+fin[11]*alphax[20])+0.3872983346207417*(alphax[7]*fin[18]+fin[7]*alphax[18]+alphax[6]*fin[17]+fin[6]*alphax[17])+0.3872983346207416*(alphax[3]*fin[16]+fin[3]*alphax[16])+0.4330127018922193*(alphax[12]*fin[15]+alphax[9]*fin[14]+alphax[8]*fin[13]+alphax[5]*fin[11]+fin[5]*alphax[11]+alphax[4]*fin[10]+alphax[2]*fin[7]+fin[2]*alphax[7]+alphax[1]*fin[6]+fin[1]*alphax[6]+alphax[0]*fin[3]+fin[0]*alphax[3]); 
  out[8] += 0.4330127018922194*(alphax[20]*fin[23]+alphax[18]*fin[22]+alphax[17]*fin[21]+alphax[16]*fin[19])+0.4330127018922193*(alphax[11]*fin[15]+alphax[7]*fin[14]+alphax[6]*fin[13]+alphax[5]*fin[12]+fin[5]*alphax[12]+alphax[3]*fin[10]+alphax[2]*fin[9]+fin[2]*alphax[9]+alphax[1]*fin[8]+fin[1]*alphax[8]+alphax[0]*fin[4]+fin[0]*alphax[4]); 
  out[11] += 0.3872983346207416*(alphax[6]*fin[20]+fin[6]*alphax[20])+0.3872983346207417*(alphax[3]*fin[18]+fin[3]*alphax[18]+alphax[11]*fin[17]+fin[11]*alphax[17])+0.3872983346207416*(alphax[7]*fin[16]+fin[7]*alphax[16])+0.4330127018922193*(alphax[8]*fin[15]+alphax[4]*fin[14]+alphax[12]*fin[13]+alphax[1]*fin[11]+fin[1]*alphax[11]+alphax[9]*fin[10]+alphax[0]*fin[7]+fin[0]*alphax[7]+alphax[5]*fin[6]+fin[5]*alphax[6]+alphax[2]*fin[3]+fin[2]*alphax[3]); 
  out[12] += 0.4330127018922193*(alphax[17]*fin[23]+alphax[16]*fin[22]+alphax[20]*fin[21]+alphax[18]*fin[19]+alphax[6]*fin[15]+alphax[3]*fin[14]+alphax[11]*fin[13]+alphax[1]*fin[12]+fin[1]*alphax[12]+alphax[7]*fin[10]+alphax[0]*fin[9]+fin[0]*alphax[9]+alphax[5]*fin[8]+fin[5]*alphax[8]+alphax[2]*fin[4]+fin[2]*alphax[4]); 
  out[13] += 0.3872983346207417*alphax[11]*fin[23]+0.3872983346207416*(alphax[7]*fin[22]+alphax[6]*fin[21]+fin[15]*alphax[20])+0.3872983346207417*(alphax[3]*fin[19]+fin[14]*alphax[18]+fin[13]*alphax[17])+0.3872983346207416*fin[10]*alphax[16]+0.4330127018922193*(alphax[5]*fin[15]+alphax[2]*fin[14]+alphax[1]*fin[13]+alphax[11]*fin[12]+fin[11]*alphax[12]+alphax[0]*fin[10]+alphax[7]*fin[9]+fin[7]*alphax[9]+alphax[6]*fin[8]+fin[6]*alphax[8]+alphax[3]*fin[4]+fin[3]*alphax[4]); 
  out[15] += 0.3872983346207417*alphax[6]*fin[23]+0.3872983346207416*(alphax[3]*fin[22]+alphax[11]*fin[21]+fin[13]*alphax[20])+0.3872983346207417*(alphax[7]*fin[19]+fin[10]*alphax[18]+fin[15]*alphax[17])+0.3872983346207416*fin[14]*alphax[16]+0.4330127018922193*(alphax[1]*fin[15]+alphax[0]*fin[14]+alphax[5]*fin[13]+alphax[6]*fin[12]+fin[6]*alphax[12]+alphax[8]*fin[11]+fin[8]*alphax[11]+alphax[2]*fin[10]+alphax[3]*fin[9]+fin[3]*alphax[9]+alphax[4]*fin[7]+fin[4]*alphax[7]); 
  out[17] += 0.4330127018922193*alphax[12]*fin[23]+0.4330127018922194*(alphax[9]*fin[22]+alphax[8]*fin[21])+0.276641667586244*alphax[20]*fin[20]+0.4330127018922194*(alphax[5]*fin[20]+fin[5]*alphax[20])+0.4330127018922193*alphax[4]*fin[19]+0.276641667586244*alphax[18]*fin[18]+0.4330127018922193*(alphax[2]*fin[18]+fin[2]*alphax[18])+0.276641667586244*alphax[17]*fin[17]+0.4330127018922193*(alphax[1]*fin[17]+fin[1]*alphax[17])+0.276641667586244*alphax[16]*fin[16]+0.4330127018922194*(alphax[0]*fin[16]+fin[0]*alphax[16])+0.3872983346207417*(alphax[11]*fin[11]+alphax[7]*fin[7]+alphax[6]*fin[6]+alphax[3]*fin[3]); 
  out[20] += 0.4330127018922194*alphax[8]*fin[23]+0.4330127018922193*(alphax[4]*fin[22]+alphax[12]*fin[21])+(0.276641667586244*alphax[17]+0.4330127018922193*alphax[1])*fin[20]+(0.276641667586244*fin[17]+0.4330127018922193*fin[1])*alphax[20]+0.4330127018922194*alphax[9]*fin[19]+(0.276641667586244*alphax[16]+0.4330127018922194*alphax[0])*fin[18]+0.276641667586244*fin[16]*alphax[18]+0.4330127018922194*(fin[0]*alphax[18]+alphax[5]*fin[17]+fin[5]*alphax[17])+0.4330127018922193*(alphax[2]*fin[16]+fin[2]*alphax[16])+0.3872983346207416*(alphax[6]*fin[11]+fin[6]*alphax[11]+alphax[3]*fin[7]+fin[3]*alphax[7]); 
  out[21] += (0.276641667586244*alphax[20]+0.4330127018922194*alphax[5])*fin[23]+(0.276641667586244*alphax[18]+0.4330127018922193*alphax[2])*fin[22]+0.276641667586244*alphax[17]*fin[21]+0.4330127018922193*(alphax[1]*fin[21]+alphax[12]*fin[20]+fin[12]*alphax[20])+0.276641667586244*alphax[16]*fin[19]+0.4330127018922194*(alphax[0]*fin[19]+alphax[9]*fin[18]+fin[9]*alphax[18]+alphax[8]*fin[17]+fin[8]*alphax[17])+0.4330127018922193*(alphax[4]*fin[16]+fin[4]*alphax[16])+0.3872983346207416*(alphax[11]*fin[15]+alphax[7]*fin[14]+alphax[6]*fin[13]+alphax[3]*fin[10]); 
  out[23] += (0.276641667586244*alphax[17]+0.4330127018922193*alphax[1])*fin[23]+(0.276641667586244*alphax[16]+0.4330127018922194*alphax[0])*fin[22]+0.276641667586244*alphax[20]*fin[21]+0.4330127018922194*(alphax[5]*fin[21]+alphax[8]*fin[20]+fin[8]*alphax[20])+0.276641667586244*alphax[18]*fin[19]+0.4330127018922193*(alphax[2]*fin[19]+alphax[4]*fin[18]+fin[4]*alphax[18]+alphax[12]*fin[17]+fin[12]*alphax[17])+0.4330127018922194*(alphax[9]*fin[16]+fin[9]*alphax[16])+0.3872983346207417*(alphax[6]*fin[15]+alphax[3]*fin[14]+alphax[11]*fin[13]+alphax[7]*fin[10]); 

  double alphaz[24] = {0.}; 
  alphaz[0] = (rdz2*((1.25*dualcurlbhatoverB_y[0]*hamil2[1])/vmap2+(0.25*dualcurlbhatoverB_y[0]*hamil2[0])/vmap2))/(m_*q_)+(0.7071067811865475*rtg33inv[0]*vmap[1]*hamil[3]*rdz2)/(m_*vmap2)+(((-0.8660254037844386*bioverJB_y[2]*hamil[5])-0.8660254037844386*bioverJB_y[0]*hamil[1])*rdx2*rdz2)/q_; 
  alphaz[1] = (rdz2*((1.25*dualcurlbhatoverB_y[1]*hamil2[1])/vmap2+(0.25*hamil2[0]*dualcurlbhatoverB_y[1])/vmap2))/(m_*q_)+(0.7071067811865475*rtg33inv[1]*vmap[1]*hamil[3]*rdz2)/(m_*vmap2)+(((-0.8660254037844386*bioverJB_y[3]*hamil[5])-0.8660254037844386*bioverJB_y[1]*hamil[1])*rdx2*rdz2)/q_; 
  alphaz[2] = (rdz2*((1.25*hamil2[1]*dualcurlbhatoverB_y[2])/vmap2+(0.25*hamil2[0]*dualcurlbhatoverB_y[2])/vmap2))/(m_*q_)+(0.7071067811865475*vmap[1]*rtg33inv[2]*hamil[3]*rdz2)/(m_*vmap2)+(((-0.8660254037844386*bioverJB_y[0]*hamil[5])-0.8660254037844386*hamil[1]*bioverJB_y[2])*rdx2*rdz2)/q_; 
  alphaz[3] = (1.118033988749895*dualcurlbhatoverB_y[0]*hamil[3]*hamil[16]*rdz2)/(m_*q_*vmap2)+(1.58113883008419*rtg33inv[0]*vmap[1]*hamil[16]*rdz2)/(m_*vmap2); 
  alphaz[4] = (((-0.8660254037844386*bioverJB_y[2]*hamil[12])-0.8660254037844386*bioverJB_y[0]*hamil[8])*rdx2*rdz2)/q_; 
  alphaz[5] = (rdz2*((1.25*hamil2[1]*dualcurlbhatoverB_y[3])/vmap2+(0.25*hamil2[0]*dualcurlbhatoverB_y[3])/vmap2))/(m_*q_)+(0.7071067811865475*vmap[1]*hamil[3]*rtg33inv[3]*rdz2)/(m_*vmap2)+(((-0.8660254037844386*bioverJB_y[1]*hamil[5])-0.8660254037844386*hamil[1]*bioverJB_y[3])*rdx2*rdz2)/q_; 
  alphaz[6] = (1.118033988749895*dualcurlbhatoverB_y[1]*hamil[3]*hamil[16]*rdz2)/(m_*q_*vmap2)+(1.58113883008419*rtg33inv[1]*vmap[1]*hamil[16]*rdz2)/(m_*vmap2); 
  alphaz[7] = (1.118033988749895*dualcurlbhatoverB_y[2]*hamil[3]*hamil[16]*rdz2)/(m_*q_*vmap2)+(1.58113883008419*vmap[1]*rtg33inv[2]*hamil[16]*rdz2)/(m_*vmap2); 
  alphaz[8] = (((-0.8660254037844386*bioverJB_y[3]*hamil[12])-0.8660254037844386*bioverJB_y[1]*hamil[8])*rdx2*rdz2)/q_; 
  alphaz[9] = (((-0.8660254037844386*bioverJB_y[0]*hamil[12])-0.8660254037844386*bioverJB_y[2]*hamil[8])*rdx2*rdz2)/q_; 
  alphaz[11] = (1.118033988749895*dualcurlbhatoverB_y[3]*hamil[3]*hamil[16]*rdz2)/(m_*q_*vmap2)+(1.58113883008419*vmap[1]*rtg33inv[3]*hamil[16]*rdz2)/(m_*vmap2); 
  alphaz[12] = (((-0.8660254037844386*bioverJB_y[1]*hamil[12])-0.8660254037844386*bioverJB_y[3]*hamil[8])*rdx2*rdz2)/q_; 
  alphaz[16] = (1.118033988749895*dualcurlbhatoverB_y[0]*hamil2[1]*rdz2)/(m_*q_*vmap2); 
  alphaz[17] = (1.118033988749895*dualcurlbhatoverB_y[1]*hamil2[1]*rdz2)/(m_*q_*vmap2); 
  alphaz[18] = (1.118033988749895*hamil2[1]*dualcurlbhatoverB_y[2]*rdz2)/(m_*q_*vmap2); 
  alphaz[20] = (1.118033988749895*hamil2[1]*dualcurlbhatoverB_y[3]*rdz2)/(m_*q_*vmap2); 


  out[2] += 0.4330127018922193*(alphaz[20]*fin[20]+alphaz[18]*fin[18]+alphaz[17]*fin[17]+alphaz[16]*fin[16]+alphaz[12]*fin[12]+alphaz[11]*fin[11]+alphaz[9]*fin[9]+alphaz[8]*fin[8]+alphaz[7]*fin[7]+alphaz[6]*fin[6]+alphaz[5]*fin[5]+alphaz[4]*fin[4]+alphaz[3]*fin[3]+alphaz[2]*fin[2]+alphaz[1]*fin[1]+alphaz[0]*fin[0]); 
  out[5] += 0.4330127018922194*(alphaz[18]*fin[20]+fin[18]*alphaz[20]+alphaz[16]*fin[17]+fin[16]*alphaz[17])+0.4330127018922193*(alphaz[9]*fin[12]+fin[9]*alphaz[12]+alphaz[7]*fin[11]+fin[7]*alphaz[11]+alphaz[4]*fin[8]+fin[4]*alphaz[8]+alphaz[3]*fin[6]+fin[3]*alphaz[6]+alphaz[2]*fin[5]+fin[2]*alphaz[5]+alphaz[0]*fin[1]+fin[0]*alphaz[1]); 
  out[7] += 0.3872983346207416*(alphaz[11]*fin[20]+fin[11]*alphaz[20])+0.3872983346207417*(alphaz[7]*fin[18]+fin[7]*alphaz[18]+alphaz[6]*fin[17]+fin[6]*alphaz[17])+0.3872983346207416*(alphaz[3]*fin[16]+fin[3]*alphaz[16])+0.4330127018922193*(alphaz[12]*fin[15]+alphaz[9]*fin[14]+alphaz[8]*fin[13]+alphaz[5]*fin[11]+fin[5]*alphaz[11]+alphaz[4]*fin[10]+alphaz[2]*fin[7]+fin[2]*alphaz[7]+alphaz[1]*fin[6]+fin[1]*alphaz[6]+alphaz[0]*fin[3]+fin[0]*alphaz[3]); 
  out[9] += 0.4330127018922194*(alphaz[20]*fin[23]+alphaz[18]*fin[22]+alphaz[17]*fin[21]+alphaz[16]*fin[19])+0.4330127018922193*(alphaz[11]*fin[15]+alphaz[7]*fin[14]+alphaz[6]*fin[13]+alphaz[5]*fin[12]+fin[5]*alphaz[12]+alphaz[3]*fin[10]+alphaz[2]*fin[9]+fin[2]*alphaz[9]+alphaz[1]*fin[8]+fin[1]*alphaz[8]+alphaz[0]*fin[4]+fin[0]*alphaz[4]); 
  out[11] += 0.3872983346207416*(alphaz[7]*fin[20]+fin[7]*alphaz[20])+0.3872983346207417*(alphaz[11]*fin[18]+fin[11]*alphaz[18]+alphaz[3]*fin[17]+fin[3]*alphaz[17])+0.3872983346207416*(alphaz[6]*fin[16]+fin[6]*alphaz[16])+0.4330127018922193*(alphaz[9]*fin[15]+alphaz[12]*fin[14]+alphaz[4]*fin[13]+alphaz[2]*fin[11]+fin[2]*alphaz[11]+alphaz[8]*fin[10]+alphaz[5]*fin[7]+fin[5]*alphaz[7]+alphaz[0]*fin[6]+fin[0]*alphaz[6]+alphaz[1]*fin[3]+fin[1]*alphaz[3]); 
  out[12] += 0.4330127018922193*(alphaz[18]*fin[23]+alphaz[20]*fin[22]+alphaz[16]*fin[21]+alphaz[17]*fin[19]+alphaz[7]*fin[15]+alphaz[11]*fin[14]+alphaz[3]*fin[13]+alphaz[2]*fin[12]+fin[2]*alphaz[12]+alphaz[6]*fin[10]+alphaz[5]*fin[9]+fin[5]*alphaz[9]+alphaz[0]*fin[8]+fin[0]*alphaz[8]+alphaz[1]*fin[4]+fin[1]*alphaz[4]); 
  out[14] += 0.3872983346207417*alphaz[11]*fin[23]+0.3872983346207416*(alphaz[7]*fin[22]+alphaz[6]*fin[21]+fin[15]*alphaz[20])+0.3872983346207417*(alphaz[3]*fin[19]+fin[14]*alphaz[18]+fin[13]*alphaz[17])+0.3872983346207416*fin[10]*alphaz[16]+0.4330127018922193*(alphaz[5]*fin[15]+alphaz[2]*fin[14]+alphaz[1]*fin[13]+alphaz[11]*fin[12]+fin[11]*alphaz[12]+alphaz[0]*fin[10]+alphaz[7]*fin[9]+fin[7]*alphaz[9]+alphaz[6]*fin[8]+fin[6]*alphaz[8]+alphaz[3]*fin[4]+fin[3]*alphaz[4]); 
  out[15] += 0.3872983346207417*alphaz[7]*fin[23]+0.3872983346207416*(alphaz[11]*fin[22]+alphaz[3]*fin[21]+fin[14]*alphaz[20])+0.3872983346207417*(alphaz[6]*fin[19]+fin[15]*alphaz[18]+fin[10]*alphaz[17])+0.3872983346207416*fin[13]*alphaz[16]+0.4330127018922193*(alphaz[2]*fin[15]+alphaz[5]*fin[14]+alphaz[0]*fin[13]+alphaz[7]*fin[12]+fin[7]*alphaz[12]+alphaz[9]*fin[11]+fin[9]*alphaz[11]+alphaz[1]*fin[10]+alphaz[3]*fin[8]+fin[3]*alphaz[8]+alphaz[4]*fin[6]+fin[4]*alphaz[6]); 
  out[18] += 0.4330127018922193*alphaz[12]*fin[23]+0.4330127018922194*(alphaz[9]*fin[22]+alphaz[8]*fin[21])+0.276641667586244*alphaz[20]*fin[20]+0.4330127018922194*(alphaz[5]*fin[20]+fin[5]*alphaz[20])+0.4330127018922193*alphaz[4]*fin[19]+0.276641667586244*alphaz[18]*fin[18]+0.4330127018922193*(alphaz[2]*fin[18]+fin[2]*alphaz[18])+0.276641667586244*alphaz[17]*fin[17]+0.4330127018922193*(alphaz[1]*fin[17]+fin[1]*alphaz[17])+0.276641667586244*alphaz[16]*fin[16]+0.4330127018922194*(alphaz[0]*fin[16]+fin[0]*alphaz[16])+0.3872983346207417*(alphaz[11]*fin[11]+alphaz[7]*fin[7]+alphaz[6]*fin[6]+alphaz[3]*fin[3]); 
  out[20] += 0.4330127018922194*alphaz[9]*fin[23]+0.4330127018922193*(alphaz[12]*fin[22]+alphaz[4]*fin[21])+(0.276641667586244*alphaz[18]+0.4330127018922193*alphaz[2])*fin[20]+(0.276641667586244*fin[18]+0.4330127018922193*fin[2])*alphaz[20]+0.4330127018922194*(alphaz[8]*fin[19]+alphaz[5]*fin[18]+fin[5]*alphaz[18])+(0.276641667586244*alphaz[16]+0.4330127018922194*alphaz[0])*fin[17]+(0.276641667586244*fin[16]+0.4330127018922194*fin[0])*alphaz[17]+0.4330127018922193*(alphaz[1]*fin[16]+fin[1]*alphaz[16])+0.3872983346207416*(alphaz[7]*fin[11]+fin[7]*alphaz[11]+alphaz[3]*fin[6]+fin[3]*alphaz[6]); 
  out[22] += (0.276641667586244*alphaz[20]+0.4330127018922194*alphaz[5])*fin[23]+(0.276641667586244*alphaz[18]+0.4330127018922193*alphaz[2])*fin[22]+0.276641667586244*alphaz[17]*fin[21]+0.4330127018922193*(alphaz[1]*fin[21]+alphaz[12]*fin[20]+fin[12]*alphaz[20])+0.276641667586244*alphaz[16]*fin[19]+0.4330127018922194*(alphaz[0]*fin[19]+alphaz[9]*fin[18]+fin[9]*alphaz[18]+alphaz[8]*fin[17]+fin[8]*alphaz[17])+0.4330127018922193*(alphaz[4]*fin[16]+fin[4]*alphaz[16])+0.3872983346207416*(alphaz[11]*fin[15]+alphaz[7]*fin[14]+alphaz[6]*fin[13]+alphaz[3]*fin[10]); 
  out[23] += (0.276641667586244*alphaz[18]+0.4330127018922193*alphaz[2])*fin[23]+(0.276641667586244*alphaz[20]+0.4330127018922194*alphaz[5])*fin[22]+0.276641667586244*alphaz[16]*fin[21]+0.4330127018922194*(alphaz[0]*fin[21]+alphaz[9]*fin[20]+fin[9]*alphaz[20])+0.276641667586244*alphaz[17]*fin[19]+0.4330127018922193*(alphaz[1]*fin[19]+alphaz[12]*fin[18]+fin[12]*alphaz[18]+alphaz[4]*fin[17]+fin[4]*alphaz[17])+0.4330127018922194*(alphaz[8]*fin[16]+fin[8]*alphaz[16])+0.3872983346207417*(alphaz[7]*fin[15]+alphaz[11]*fin[14]+alphaz[3]*fin[13]+alphaz[6]*fin[10]); 

  double alphavpar[24] = {0.}; 
  alphavpar[0] = (rdx2*((0.25*dualcurlbhatoverB_x[2]*hamil[3]*hamil[5])/vmap2+(0.25*dualcurlbhatoverB_x[0]*hamil[1]*hamil[3])/vmap2)+rdz2*((0.25*dualcurlbhatoverB_y[1]*hamil[3]*hamil[5])/vmap2+(0.25*dualcurlbhatoverB_y[0]*hamil[2]*hamil[3])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[1]*hamil[5])/vmap2)-(0.7071067811865475*rtg33inv[0]*hamil[2])/vmap2))/m_; 
  alphavpar[1] = (rdx2*((0.25*dualcurlbhatoverB_x[3]*hamil[3]*hamil[5])/vmap2+(0.25*dualcurlbhatoverB_x[1]*hamil[1]*hamil[3])/vmap2)+rdz2*((0.25*dualcurlbhatoverB_y[0]*hamil[3]*hamil[5])/vmap2+(0.25*dualcurlbhatoverB_y[1]*hamil[2]*hamil[3])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[0]*hamil[5])/vmap2)-(0.7071067811865475*rtg33inv[1]*hamil[2])/vmap2))/m_; 
  alphavpar[2] = (rdz2*((0.25*dualcurlbhatoverB_y[3]*hamil[3]*hamil[5])/vmap2+(0.25*dualcurlbhatoverB_y[2]*hamil[2]*hamil[3])/vmap2)+rdx2*((0.25*dualcurlbhatoverB_x[0]*hamil[3]*hamil[5])/vmap2+(0.25*hamil[1]*dualcurlbhatoverB_x[2]*hamil[3])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[3]*hamil[5])/vmap2)-(0.7071067811865475*hamil[2]*rtg33inv[2])/vmap2))/m_; 
  alphavpar[3] = (rdx2*((0.5590169943749475*dualcurlbhatoverB_x[2]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_x[0]*hamil[1]*hamil[16])/vmap2)+rdz2*((0.5590169943749475*dualcurlbhatoverB_y[1]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[0]*hamil[2]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[4] = (rdx2*((0.25*dualcurlbhatoverB_x[2]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_x[0]*hamil[3]*hamil[8])/vmap2)+rdz2*((0.25*dualcurlbhatoverB_y[1]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_y[0]*hamil[3]*hamil[9])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[1]*hamil[12])/vmap2)-(0.7071067811865475*rtg33inv[0]*hamil[9])/vmap2))/m_; 
  alphavpar[5] = (rdz2*((0.25*dualcurlbhatoverB_y[2]*hamil[3]*hamil[5])/vmap2+(0.25*hamil[2]*dualcurlbhatoverB_y[3]*hamil[3])/vmap2)+rdx2*((0.25*dualcurlbhatoverB_x[1]*hamil[3]*hamil[5])/vmap2+(0.25*hamil[1]*dualcurlbhatoverB_x[3]*hamil[3])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[2]*hamil[5])/vmap2)-(0.7071067811865475*hamil[2]*rtg33inv[3])/vmap2))/m_; 
  alphavpar[6] = (rdx2*((0.5590169943749475*dualcurlbhatoverB_x[3]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_x[1]*hamil[1]*hamil[16])/vmap2)+rdz2*((0.5590169943749475*dualcurlbhatoverB_y[0]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[1]*hamil[2]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[7] = (rdz2*((0.5590169943749475*dualcurlbhatoverB_y[3]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[2]*hamil[2]*hamil[16])/vmap2)+rdx2*((0.5590169943749475*dualcurlbhatoverB_x[0]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*hamil[1]*dualcurlbhatoverB_x[2]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[8] = (rdx2*((0.25*dualcurlbhatoverB_x[3]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_x[1]*hamil[3]*hamil[8])/vmap2)+rdz2*((0.25*dualcurlbhatoverB_y[0]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_y[1]*hamil[3]*hamil[9])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[0]*hamil[12])/vmap2)-(0.7071067811865475*rtg33inv[1]*hamil[9])/vmap2))/m_; 
  alphavpar[9] = (rdz2*((0.25*dualcurlbhatoverB_y[3]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_y[2]*hamil[3]*hamil[9])/vmap2)+rdx2*((0.25*dualcurlbhatoverB_x[0]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_x[2]*hamil[3]*hamil[8])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[3]*hamil[12])/vmap2)-(0.7071067811865475*rtg33inv[2]*hamil[9])/vmap2))/m_; 
  alphavpar[10] = (rdx2*((0.5590169943749475*dualcurlbhatoverB_x[2]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_x[0]*hamil[8]*hamil[16])/vmap2)+rdz2*((0.5590169943749475*dualcurlbhatoverB_y[1]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[0]*hamil[9]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[11] = (rdz2*((0.5590169943749475*dualcurlbhatoverB_y[2]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*hamil[2]*dualcurlbhatoverB_y[3]*hamil[16])/vmap2)+rdx2*((0.5590169943749475*dualcurlbhatoverB_x[1]*hamil[5]*hamil[16])/vmap2+(0.5590169943749475*hamil[1]*dualcurlbhatoverB_x[3]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[12] = (rdz2*((0.25*dualcurlbhatoverB_y[2]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_y[3]*hamil[3]*hamil[9])/vmap2)+rdx2*((0.25*dualcurlbhatoverB_x[1]*hamil[3]*hamil[12])/vmap2+(0.25*dualcurlbhatoverB_x[3]*hamil[3]*hamil[8])/vmap2))/(m_*q_)+(vmap[1]*rdz2*((-(0.7071067811865475*rtg33inv[2]*hamil[12])/vmap2)-(0.7071067811865475*rtg33inv[3]*hamil[9])/vmap2))/m_; 
  alphavpar[13] = (rdx2*((0.5590169943749475*dualcurlbhatoverB_x[3]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_x[1]*hamil[8]*hamil[16])/vmap2)+rdz2*((0.5590169943749475*dualcurlbhatoverB_y[0]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[1]*hamil[9]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[14] = (rdz2*((0.5590169943749475*dualcurlbhatoverB_y[3]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[2]*hamil[9]*hamil[16])/vmap2)+rdx2*((0.5590169943749475*dualcurlbhatoverB_x[0]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_x[2]*hamil[8]*hamil[16])/vmap2))/(m_*q_); 
  alphavpar[15] = (rdz2*((0.5590169943749475*dualcurlbhatoverB_y[2]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_y[3]*hamil[9]*hamil[16])/vmap2)+rdx2*((0.5590169943749475*dualcurlbhatoverB_x[1]*hamil[12]*hamil[16])/vmap2+(0.5590169943749475*dualcurlbhatoverB_x[3]*hamil[8]*hamil[16])/vmap2))/(m_*q_); 


  out[3] += 0.4330127018922193*(alphavpar[15]*fin[15]+alphavpar[14]*fin[14]+alphavpar[13]*fin[13]+alphavpar[12]*fin[12]+alphavpar[11]*fin[11]+alphavpar[10]*fin[10]+alphavpar[9]*fin[9]+alphavpar[8]*fin[8]+alphavpar[7]*fin[7]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[6] += 0.4330127018922193*(alphavpar[14]*fin[15]+fin[14]*alphavpar[15]+alphavpar[10]*fin[13]+fin[10]*alphavpar[13]+alphavpar[9]*fin[12]+fin[9]*alphavpar[12]+alphavpar[7]*fin[11]+fin[7]*alphavpar[11]+alphavpar[4]*fin[8]+fin[4]*alphavpar[8]+alphavpar[3]*fin[6]+fin[3]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphavpar[13]*fin[15]+fin[13]*alphavpar[15]+alphavpar[10]*fin[14]+fin[10]*alphavpar[14]+alphavpar[8]*fin[12]+fin[8]*alphavpar[12]+alphavpar[6]*fin[11]+fin[6]*alphavpar[11]+alphavpar[4]*fin[9]+fin[4]*alphavpar[9]+alphavpar[3]*fin[7]+fin[3]*alphavpar[7]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[10] += 0.4330127018922193*(alphavpar[11]*fin[15]+fin[11]*alphavpar[15]+alphavpar[7]*fin[14]+fin[7]*alphavpar[14]+alphavpar[6]*fin[13]+fin[6]*alphavpar[13]+alphavpar[5]*fin[12]+fin[5]*alphavpar[12]+alphavpar[3]*fin[10]+fin[3]*alphavpar[10]+alphavpar[2]*fin[9]+fin[2]*alphavpar[9]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]); 
  out[11] += 0.4330127018922193*(alphavpar[10]*fin[15]+fin[10]*alphavpar[15]+alphavpar[13]*fin[14]+fin[13]*alphavpar[14]+alphavpar[4]*fin[12]+fin[4]*alphavpar[12]+alphavpar[3]*fin[11]+fin[3]*alphavpar[11]+alphavpar[8]*fin[9]+fin[8]*alphavpar[9]+alphavpar[6]*fin[7]+fin[6]*alphavpar[7]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[13] += 0.4330127018922193*(alphavpar[7]*fin[15]+fin[7]*alphavpar[15]+alphavpar[11]*fin[14]+fin[11]*alphavpar[14]+alphavpar[3]*fin[13]+fin[3]*alphavpar[13]+alphavpar[2]*fin[12]+fin[2]*alphavpar[12]+alphavpar[6]*fin[10]+fin[6]*alphavpar[10]+alphavpar[5]*fin[9]+fin[5]*alphavpar[9]+alphavpar[0]*fin[8]+fin[0]*alphavpar[8]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]); 
  out[14] += 0.4330127018922193*(alphavpar[6]*fin[15]+fin[6]*alphavpar[15]+alphavpar[3]*fin[14]+fin[3]*alphavpar[14]+alphavpar[11]*fin[13]+fin[11]*alphavpar[13]+alphavpar[1]*fin[12]+fin[1]*alphavpar[12]+alphavpar[7]*fin[10]+fin[7]*alphavpar[10]+alphavpar[0]*fin[9]+fin[0]*alphavpar[9]+alphavpar[5]*fin[8]+fin[5]*alphavpar[8]+alphavpar[2]*fin[4]+fin[2]*alphavpar[4]); 
  out[15] += 0.4330127018922193*(alphavpar[3]*fin[15]+fin[3]*alphavpar[15]+alphavpar[6]*fin[14]+fin[6]*alphavpar[14]+alphavpar[7]*fin[13]+fin[7]*alphavpar[13]+alphavpar[0]*fin[12]+fin[0]*alphavpar[12]+alphavpar[10]*fin[11]+fin[10]*alphavpar[11]+alphavpar[1]*fin[9]+fin[1]*alphavpar[9]+alphavpar[2]*fin[8]+fin[2]*alphavpar[8]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]); 
  out[16] += 0.8660254037844387*alphavpar[15]*fin[23]+0.8660254037844386*(alphavpar[14]*fin[22]+alphavpar[13]*fin[21]+alphavpar[11]*fin[20])+0.8660254037844387*(alphavpar[10]*fin[19]+alphavpar[7]*fin[18]+alphavpar[6]*fin[17])+0.8660254037844386*alphavpar[3]*fin[16]+0.9682458365518543*(alphavpar[12]*fin[15]+fin[12]*alphavpar[15]+alphavpar[9]*fin[14]+fin[9]*alphavpar[14]+alphavpar[8]*fin[13]+fin[8]*alphavpar[13]+alphavpar[5]*fin[11]+fin[5]*alphavpar[11]+alphavpar[4]*fin[10]+fin[4]*alphavpar[10]+alphavpar[2]*fin[7]+fin[2]*alphavpar[7]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[17] += 0.8660254037844386*alphavpar[14]*fin[23]+0.8660254037844387*(alphavpar[15]*fin[22]+alphavpar[10]*fin[21]+alphavpar[7]*fin[20])+0.8660254037844386*(alphavpar[13]*fin[19]+alphavpar[11]*fin[18]+alphavpar[3]*fin[17])+0.8660254037844387*alphavpar[6]*fin[16]+0.9682458365518543*(alphavpar[9]*fin[15]+fin[9]*alphavpar[15]+alphavpar[12]*fin[14]+fin[12]*alphavpar[14]+alphavpar[4]*fin[13]+fin[4]*alphavpar[13]+alphavpar[2]*fin[11]+fin[2]*alphavpar[11]+alphavpar[8]*fin[10]+fin[8]*alphavpar[10]+alphavpar[5]*fin[7]+fin[5]*alphavpar[7]+alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[18] += 0.8660254037844386*alphavpar[13]*fin[23]+0.8660254037844387*(alphavpar[10]*fin[22]+alphavpar[15]*fin[21]+alphavpar[6]*fin[20])+0.8660254037844386*(alphavpar[14]*fin[19]+alphavpar[3]*fin[18]+alphavpar[11]*fin[17])+0.8660254037844387*alphavpar[7]*fin[16]+0.9682458365518543*(alphavpar[8]*fin[15]+fin[8]*alphavpar[15]+alphavpar[4]*fin[14]+fin[4]*alphavpar[14]+alphavpar[12]*fin[13]+fin[12]*alphavpar[13]+alphavpar[1]*fin[11]+fin[1]*alphavpar[11]+alphavpar[9]*fin[10]+fin[9]*alphavpar[10]+alphavpar[0]*fin[7]+fin[0]*alphavpar[7]+alphavpar[5]*fin[6]+fin[5]*alphavpar[6]+alphavpar[2]*fin[3]+fin[2]*alphavpar[3]); 
  out[19] += 0.8660254037844386*alphavpar[11]*fin[23]+0.8660254037844387*(alphavpar[7]*fin[22]+alphavpar[6]*fin[21]+alphavpar[15]*fin[20])+0.8660254037844386*(alphavpar[3]*fin[19]+alphavpar[14]*fin[18]+alphavpar[13]*fin[17])+0.8660254037844387*alphavpar[10]*fin[16]+0.9682458365518543*(alphavpar[5]*fin[15]+fin[5]*alphavpar[15]+alphavpar[2]*fin[14]+fin[2]*alphavpar[14]+alphavpar[1]*fin[13]+fin[1]*alphavpar[13]+alphavpar[11]*fin[12]+fin[11]*alphavpar[12]+alphavpar[0]*fin[10]+fin[0]*alphavpar[10]+alphavpar[7]*fin[9]+fin[7]*alphavpar[9]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4]); 
  out[20] += 0.8660254037844387*alphavpar[10]*fin[23]+0.8660254037844386*(alphavpar[13]*fin[22]+alphavpar[14]*fin[21]+alphavpar[3]*fin[20])+0.8660254037844387*(alphavpar[15]*fin[19]+alphavpar[6]*fin[18]+alphavpar[7]*fin[17])+0.8660254037844386*alphavpar[11]*fin[16]+0.9682458365518543*(alphavpar[4]*fin[15]+fin[4]*alphavpar[15]+alphavpar[8]*fin[14]+fin[8]*alphavpar[14]+alphavpar[9]*fin[13]+fin[9]*alphavpar[13]+alphavpar[10]*fin[12]+fin[10]*alphavpar[12]+alphavpar[0]*fin[11]+fin[0]*alphavpar[11]+alphavpar[1]*fin[7]+fin[1]*alphavpar[7]+alphavpar[2]*fin[6]+fin[2]*alphavpar[6]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5]); 
  out[21] += 0.8660254037844387*alphavpar[7]*fin[23]+0.8660254037844386*(alphavpar[11]*fin[22]+alphavpar[3]*fin[21]+alphavpar[14]*fin[20])+0.8660254037844387*(alphavpar[6]*fin[19]+alphavpar[15]*fin[18]+alphavpar[10]*fin[17])+0.8660254037844386*alphavpar[13]*fin[16]+0.9682458365518543*(alphavpar[2]*fin[15]+fin[2]*alphavpar[15]+alphavpar[5]*fin[14]+fin[5]*alphavpar[14]+alphavpar[0]*fin[13]+fin[0]*alphavpar[13]+alphavpar[7]*fin[12]+fin[7]*alphavpar[12]+alphavpar[9]*fin[11]+fin[9]*alphavpar[11]+alphavpar[1]*fin[10]+fin[1]*alphavpar[10]+alphavpar[3]*fin[8]+fin[3]*alphavpar[8]+alphavpar[4]*fin[6]+fin[4]*alphavpar[6]); 
  out[22] += 0.8660254037844387*alphavpar[6]*fin[23]+0.8660254037844386*(alphavpar[3]*fin[22]+alphavpar[11]*fin[21]+alphavpar[13]*fin[20])+0.8660254037844387*(alphavpar[7]*fin[19]+alphavpar[10]*fin[18]+alphavpar[15]*fin[17])+0.8660254037844386*alphavpar[14]*fin[16]+0.9682458365518543*(alphavpar[1]*fin[15]+fin[1]*alphavpar[15]+alphavpar[0]*fin[14]+fin[0]*alphavpar[14]+alphavpar[5]*fin[13]+fin[5]*alphavpar[13]+alphavpar[6]*fin[12]+fin[6]*alphavpar[12]+alphavpar[8]*fin[11]+fin[8]*alphavpar[11]+alphavpar[2]*fin[10]+fin[2]*alphavpar[10]+alphavpar[3]*fin[9]+fin[3]*alphavpar[9]+alphavpar[4]*fin[7]+fin[4]*alphavpar[7]); 
  out[23] += 0.8660254037844386*alphavpar[3]*fin[23]+0.8660254037844387*(alphavpar[6]*fin[22]+alphavpar[7]*fin[21]+alphavpar[10]*fin[20])+0.8660254037844386*(alphavpar[11]*fin[19]+alphavpar[13]*fin[18]+alphavpar[14]*fin[17])+0.8660254037844387*alphavpar[15]*fin[16]+0.9682458365518543*(alphavpar[0]*fin[15]+fin[0]*alphavpar[15]+alphavpar[1]*fin[14]+fin[1]*alphavpar[14]+alphavpar[2]*fin[13]+fin[2]*alphavpar[13]+alphavpar[3]*fin[12]+fin[3]*alphavpar[12]+alphavpar[4]*fin[11]+fin[4]*alphavpar[11]+alphavpar[5]*fin[10]+fin[5]*alphavpar[10]+alphavpar[6]*fin[9]+fin[6]*alphavpar[9]+alphavpar[7]*fin[8]+fin[7]*alphavpar[8]); 

  return 0.; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[2]; 
  out[3] += -(0.8660254037844386*(apardot[3]*fin[5]+apardot[2]*fin[2]+apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[6] += -(0.8660254037844386*(apardot[2]*fin[5]+fin[2]*apardot[3]+apardot[0]*fin[1]+fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.8660254037844386*(apardot[1]*fin[5]+fin[1]*apardot[3]+apardot[0]*fin[2]+fin[0]*apardot[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.8660254037844386*(apardot[3]*fin[12]+apardot[2]*fin[9]+apardot[1]*fin[8]+apardot[0]*fin[4])*q_*rdvpar2)/m_; 
  out[11] += -(0.8660254037844386*(apardot[0]*fin[5]+fin[0]*apardot[3]+apardot[1]*fin[2]+fin[1]*apardot[2])*q_*rdvpar2)/m_; 
  out[13] += -(0.8660254037844386*(apardot[2]*fin[12]+apardot[3]*fin[9]+apardot[0]*fin[8]+apardot[1]*fin[4])*q_*rdvpar2)/m_; 
  out[14] += -(0.8660254037844386*(apardot[1]*fin[12]+apardot[0]*fin[9]+apardot[3]*fin[8]+apardot[2]*fin[4])*q_*rdvpar2)/m_; 
  out[15] += -(0.8660254037844386*(apardot[0]*fin[12]+apardot[1]*fin[9]+apardot[2]*fin[8]+apardot[3]*fin[4])*q_*rdvpar2)/m_; 
  out[16] += -(1.936491673103709*(apardot[3]*fin[11]+apardot[2]*fin[7]+apardot[1]*fin[6]+apardot[0]*fin[3])*q_*rdvpar2)/m_; 
  out[17] += -(1.936491673103709*(apardot[2]*fin[11]+apardot[3]*fin[7]+apardot[0]*fin[6]+apardot[1]*fin[3])*q_*rdvpar2)/m_; 
  out[18] += -(1.936491673103709*(apardot[1]*fin[11]+apardot[0]*fin[7]+apardot[3]*fin[6]+apardot[2]*fin[3])*q_*rdvpar2)/m_; 
  out[19] += -(1.936491673103709*(apardot[3]*fin[15]+apardot[2]*fin[14]+apardot[1]*fin[13]+apardot[0]*fin[10])*q_*rdvpar2)/m_; 
  out[20] += -(1.936491673103709*(apardot[0]*fin[11]+apardot[1]*fin[7]+apardot[2]*fin[6]+apardot[3]*fin[3])*q_*rdvpar2)/m_; 
  out[21] += -(1.936491673103709*(apardot[2]*fin[15]+apardot[3]*fin[14]+apardot[0]*fin[13]+apardot[1]*fin[10])*q_*rdvpar2)/m_; 
  out[22] += -(1.936491673103709*(apardot[1]*fin[15]+apardot[0]*fin[14]+apardot[3]*fin[13]+apardot[2]*fin[10])*q_*rdvpar2)/m_; 
  out[23] += -(1.936491673103709*(apardot[0]*fin[15]+apardot[1]*fin[14]+apardot[2]*fin[13]+apardot[3]*fin[10])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
