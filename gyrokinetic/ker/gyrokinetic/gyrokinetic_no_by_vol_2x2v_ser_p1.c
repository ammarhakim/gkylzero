#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_no_by_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv,
    const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
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

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];
  double rdmu2 = 2.0/dxv[3];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[2]*dxv[2];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[4];
  const double *b_z = &b_i[8];

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

  double alphax[24] = {0.}; 

  double alphaz[24] = {0.}; 
  alphaz[0] = ((0.3535533905932737*cmag[3]*hamil[3]*jacobtot_inv[3]+0.3535533905932737*cmag[2]*jacobtot_inv[2]*hamil[3]+0.3535533905932737*cmag[1]*jacobtot_inv[1]*hamil[3]+0.3535533905932737*cmag[0]*jacobtot_inv[0]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphaz[1] = ((0.3535533905932737*cmag[2]*hamil[3]*jacobtot_inv[3]+0.3535533905932737*jacobtot_inv[2]*cmag[3]*hamil[3]+0.3535533905932737*cmag[0]*jacobtot_inv[1]*hamil[3]+0.3535533905932737*jacobtot_inv[0]*cmag[1]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphaz[2] = ((0.3535533905932737*cmag[1]*hamil[3]*jacobtot_inv[3]+0.3535533905932737*jacobtot_inv[1]*cmag[3]*hamil[3]+0.3535533905932737*cmag[0]*jacobtot_inv[2]*hamil[3]+0.3535533905932737*jacobtot_inv[0]*cmag[2]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphaz[3] = ((0.7905694150420947*cmag[3]*jacobtot_inv[3]*hamil[16]+0.7905694150420947*cmag[2]*jacobtot_inv[2]*hamil[16]+0.7905694150420947*cmag[1]*jacobtot_inv[1]*hamil[16]+0.7905694150420947*cmag[0]*jacobtot_inv[0]*hamil[16])*rdz2)/(vmap[1]*m_); 
  alphaz[5] = ((0.3535533905932737*cmag[0]*hamil[3]*jacobtot_inv[3]+0.3535533905932737*jacobtot_inv[0]*cmag[3]*hamil[3]+0.3535533905932737*cmag[1]*jacobtot_inv[2]*hamil[3]+0.3535533905932737*jacobtot_inv[1]*cmag[2]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphaz[6] = ((0.7905694150420947*cmag[2]*jacobtot_inv[3]*hamil[16]+0.7905694150420947*jacobtot_inv[2]*cmag[3]*hamil[16]+0.7905694150420947*cmag[0]*jacobtot_inv[1]*hamil[16]+0.7905694150420947*jacobtot_inv[0]*cmag[1]*hamil[16])*rdz2)/(vmap[1]*m_); 
  alphaz[7] = ((0.7905694150420947*cmag[1]*jacobtot_inv[3]*hamil[16]+0.7905694150420947*jacobtot_inv[1]*cmag[3]*hamil[16]+0.7905694150420947*cmag[0]*jacobtot_inv[2]*hamil[16]+0.7905694150420947*jacobtot_inv[0]*cmag[2]*hamil[16])*rdz2)/(vmap[1]*m_); 
  alphaz[11] = ((0.7905694150420947*cmag[0]*jacobtot_inv[3]*hamil[16]+0.7905694150420947*jacobtot_inv[0]*cmag[3]*hamil[16]+0.7905694150420947*cmag[1]*jacobtot_inv[2]*hamil[16]+0.7905694150420947*jacobtot_inv[1]*cmag[2]*hamil[16])*rdz2)/(vmap[1]*m_); 

  double alphavpar[24] = {0.}; 
  alphavpar[0] = (((-0.3535533905932737*cmag[2]*jacobtot_inv[3]*hamil[5])-0.3535533905932737*jacobtot_inv[2]*cmag[3]*hamil[5]-0.3535533905932737*cmag[0]*jacobtot_inv[1]*hamil[5]-0.3535533905932737*jacobtot_inv[0]*cmag[1]*hamil[5]-0.3535533905932737*hamil[2]*cmag[3]*jacobtot_inv[3]-0.3535533905932737*cmag[2]*hamil[2]*jacobtot_inv[2]-0.3535533905932737*cmag[1]*jacobtot_inv[1]*hamil[2]-0.3535533905932737*cmag[0]*jacobtot_inv[0]*hamil[2])*rdz2)/(vmap[1]*m_); 
  alphavpar[1] = (((-0.6363961030678926*cmag[3]*jacobtot_inv[3]*hamil[5])-0.3535533905932737*cmag[2]*jacobtot_inv[2]*hamil[5]-0.6363961030678926*cmag[1]*jacobtot_inv[1]*hamil[5]-0.3535533905932737*cmag[0]*jacobtot_inv[0]*hamil[5]-0.3535533905932737*cmag[2]*hamil[2]*jacobtot_inv[3]-0.3535533905932737*hamil[2]*jacobtot_inv[2]*cmag[3]-0.3535533905932737*cmag[0]*jacobtot_inv[1]*hamil[2]-0.3535533905932737*jacobtot_inv[0]*cmag[1]*hamil[2])*rdz2)/(vmap[1]*m_); 
  alphavpar[2] = (((-0.3535533905932737*cmag[0]*jacobtot_inv[3]*hamil[5])-0.3535533905932737*jacobtot_inv[0]*cmag[3]*hamil[5]-0.3535533905932737*cmag[1]*jacobtot_inv[2]*hamil[5]-0.3535533905932737*jacobtot_inv[1]*cmag[2]*hamil[5]-0.3535533905932737*cmag[1]*hamil[2]*jacobtot_inv[3]-0.3535533905932737*jacobtot_inv[1]*hamil[2]*cmag[3]-0.3535533905932737*cmag[0]*hamil[2]*jacobtot_inv[2]-0.3535533905932737*jacobtot_inv[0]*cmag[2]*hamil[2])*rdz2)/(vmap[1]*m_); 
  alphavpar[4] = (((-0.3535533905932737*cmag[2]*jacobtot_inv[3]*hamil[12])-0.3535533905932737*jacobtot_inv[2]*cmag[3]*hamil[12]-0.3535533905932737*cmag[0]*jacobtot_inv[1]*hamil[12]-0.3535533905932737*jacobtot_inv[0]*cmag[1]*hamil[12]-0.3535533905932737*cmag[3]*jacobtot_inv[3]*hamil[9]-0.3535533905932737*cmag[2]*jacobtot_inv[2]*hamil[9]-0.3535533905932737*cmag[1]*jacobtot_inv[1]*hamil[9]-0.3535533905932737*cmag[0]*jacobtot_inv[0]*hamil[9])*rdz2)/(vmap[1]*m_); 
  alphavpar[5] = (((-0.6363961030678926*cmag[1]*jacobtot_inv[3]*hamil[5])-0.6363961030678926*jacobtot_inv[1]*cmag[3]*hamil[5]-0.3535533905932737*cmag[0]*jacobtot_inv[2]*hamil[5]-0.3535533905932737*jacobtot_inv[0]*cmag[2]*hamil[5]-0.3535533905932737*cmag[0]*hamil[2]*jacobtot_inv[3]-0.3535533905932737*jacobtot_inv[0]*hamil[2]*cmag[3]-0.3535533905932737*cmag[1]*hamil[2]*jacobtot_inv[2]-0.3535533905932737*jacobtot_inv[1]*cmag[2]*hamil[2])*rdz2)/(vmap[1]*m_); 
  alphavpar[8] = (((-0.6363961030678926*cmag[3]*jacobtot_inv[3]*hamil[12])-0.3535533905932737*cmag[2]*jacobtot_inv[2]*hamil[12]-0.6363961030678926*cmag[1]*jacobtot_inv[1]*hamil[12]-0.3535533905932737*cmag[0]*jacobtot_inv[0]*hamil[12]-0.3535533905932737*cmag[2]*jacobtot_inv[3]*hamil[9]-0.3535533905932737*jacobtot_inv[2]*cmag[3]*hamil[9]-0.3535533905932737*cmag[0]*jacobtot_inv[1]*hamil[9]-0.3535533905932737*jacobtot_inv[0]*cmag[1]*hamil[9])*rdz2)/(vmap[1]*m_); 
  alphavpar[9] = (((-0.3535533905932737*cmag[0]*jacobtot_inv[3]*hamil[12])-0.3535533905932737*jacobtot_inv[0]*cmag[3]*hamil[12]-0.3535533905932737*cmag[1]*jacobtot_inv[2]*hamil[12]-0.3535533905932737*jacobtot_inv[1]*cmag[2]*hamil[12]-0.3535533905932737*cmag[1]*jacobtot_inv[3]*hamil[9]-0.3535533905932737*jacobtot_inv[1]*cmag[3]*hamil[9]-0.3535533905932737*cmag[0]*jacobtot_inv[2]*hamil[9]-0.3535533905932737*jacobtot_inv[0]*cmag[2]*hamil[9])*rdz2)/(vmap[1]*m_); 
  alphavpar[12] = (((-0.6363961030678926*cmag[1]*jacobtot_inv[3]*hamil[12])-0.6363961030678926*jacobtot_inv[1]*cmag[3]*hamil[12]-0.3535533905932737*cmag[0]*jacobtot_inv[2]*hamil[12]-0.3535533905932737*jacobtot_inv[0]*cmag[2]*hamil[12]-0.3535533905932737*cmag[0]*jacobtot_inv[3]*hamil[9]-0.3535533905932737*jacobtot_inv[0]*cmag[3]*hamil[9]-0.3535533905932737*cmag[1]*jacobtot_inv[2]*hamil[9]-0.3535533905932737*jacobtot_inv[1]*cmag[2]*hamil[9])*rdz2)/(vmap[1]*m_); 

  out[2] += 0.4330127018922193*(alphaz[11]*fin[11]+alphaz[7]*fin[7]+alphaz[6]*fin[6]+alphaz[5]*fin[5]+alphaz[3]*fin[3]+alphaz[2]*fin[2]+alphaz[1]*fin[1]+alphaz[0]*fin[0]); 
  out[3] += 0.4330127018922193*(alphavpar[12]*fin[12]+alphavpar[9]*fin[9]+alphavpar[8]*fin[8]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[5] += 0.4330127018922193*(alphaz[7]*fin[11]+fin[7]*alphaz[11]+alphaz[3]*fin[6]+fin[3]*alphaz[6]+alphaz[2]*fin[5]+fin[2]*alphaz[5]+alphaz[0]*fin[1]+fin[0]*alphaz[1]); 
  out[6] += 0.4330127018922193*(alphavpar[9]*fin[12]+fin[9]*alphavpar[12]+alphavpar[4]*fin[8]+fin[4]*alphavpar[8]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[7] += 0.3872983346207416*alphaz[11]*fin[20]+0.3872983346207417*(alphaz[7]*fin[18]+alphaz[6]*fin[17])+0.3872983346207416*alphaz[3]*fin[16]+0.4330127018922193*(alphavpar[8]*fin[12]+fin[8]*alphavpar[12]+alphaz[5]*fin[11]+fin[5]*alphaz[11]+alphavpar[4]*fin[9]+fin[4]*alphavpar[9]+alphaz[2]*fin[7]+fin[2]*alphaz[7]+alphaz[1]*fin[6]+fin[1]*alphaz[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphaz[0]*fin[3]+fin[0]*alphaz[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[9] += 0.4330127018922193*(alphaz[11]*fin[15]+alphaz[7]*fin[14]+alphaz[6]*fin[13]+alphaz[5]*fin[12]+alphaz[3]*fin[10]+alphaz[2]*fin[9]+alphaz[1]*fin[8]+alphaz[0]*fin[4]); 
  out[10] += 0.4330127018922193*(alphavpar[5]*fin[12]+fin[5]*alphavpar[12]+alphavpar[2]*fin[9]+fin[2]*alphavpar[9]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]); 
  out[11] += 0.3872983346207416*alphaz[7]*fin[20]+0.3872983346207417*(alphaz[11]*fin[18]+alphaz[3]*fin[17])+0.3872983346207416*alphaz[6]*fin[16]+0.4330127018922193*(alphavpar[4]*fin[12]+fin[4]*alphavpar[12]+alphaz[2]*fin[11]+fin[2]*alphaz[11]+alphavpar[8]*fin[9]+fin[8]*alphavpar[9]+alphaz[5]*fin[7]+fin[5]*alphaz[7]+alphaz[0]*fin[6]+fin[0]*alphaz[6]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphaz[1]*fin[3]+fin[1]*alphaz[3]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*(alphaz[7]*fin[15]+alphaz[11]*fin[14]+alphaz[3]*fin[13]+alphaz[2]*fin[12]+alphaz[6]*fin[10]+alphaz[5]*fin[9]+alphaz[0]*fin[8]+alphaz[1]*fin[4]); 
  out[13] += 0.4330127018922193*(alphavpar[2]*fin[12]+fin[2]*alphavpar[12]+alphavpar[5]*fin[9]+fin[5]*alphavpar[9]+alphavpar[0]*fin[8]+fin[0]*alphavpar[8]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]); 
  out[14] += 0.3872983346207417*alphaz[11]*fin[23]+0.3872983346207416*(alphaz[7]*fin[22]+alphaz[6]*fin[21])+0.3872983346207417*alphaz[3]*fin[19]+0.4330127018922193*(alphaz[5]*fin[15]+alphaz[2]*fin[14]+alphaz[1]*fin[13]+(alphaz[11]+alphavpar[1])*fin[12]+fin[1]*alphavpar[12]+alphaz[0]*fin[10]+(alphaz[7]+alphavpar[0])*fin[9]+fin[0]*alphavpar[9]+(alphaz[6]+alphavpar[5])*fin[8]+fin[5]*alphavpar[8]+(alphaz[3]+alphavpar[2])*fin[4]+fin[2]*alphavpar[4]); 
  out[15] += 0.3872983346207417*alphaz[7]*fin[23]+0.3872983346207416*(alphaz[11]*fin[22]+alphaz[3]*fin[21])+0.3872983346207417*alphaz[6]*fin[19]+0.4330127018922193*(alphaz[2]*fin[15]+alphaz[5]*fin[14]+alphaz[0]*fin[13]+(alphaz[7]+alphavpar[0])*fin[12]+fin[0]*alphavpar[12]+fin[9]*alphaz[11]+alphaz[1]*fin[10]+alphavpar[1]*fin[9]+fin[1]*alphavpar[9]+(alphaz[3]+alphavpar[2])*fin[8]+fin[2]*alphavpar[8]+fin[4]*alphaz[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]); 
  out[16] += 0.9682458365518543*(alphavpar[12]*fin[15]+alphavpar[9]*fin[14]+alphavpar[8]*fin[13]+alphavpar[5]*fin[11]+alphavpar[4]*fin[10]+alphavpar[2]*fin[7]+alphavpar[1]*fin[6]+alphavpar[0]*fin[3]); 
  out[17] += 0.9682458365518543*(alphavpar[9]*fin[15]+alphavpar[12]*fin[14]+alphavpar[4]*fin[13]+alphavpar[2]*fin[11]+alphavpar[8]*fin[10]+alphavpar[5]*fin[7]+alphavpar[0]*fin[6]+alphavpar[1]*fin[3]); 
  out[18] += 0.4330127018922194*alphaz[5]*fin[20]+0.4330127018922193*(alphaz[2]*fin[18]+alphaz[1]*fin[17])+0.4330127018922194*alphaz[0]*fin[16]+0.9682458365518543*(alphavpar[8]*fin[15]+alphavpar[4]*fin[14]+alphavpar[12]*fin[13])+0.3872983346207417*alphaz[11]*fin[11]+0.9682458365518543*(alphavpar[1]*fin[11]+alphavpar[9]*fin[10])+(0.3872983346207417*alphaz[7]+0.9682458365518543*alphavpar[0])*fin[7]+(0.3872983346207417*alphaz[6]+0.9682458365518543*alphavpar[5])*fin[6]+(0.3872983346207417*alphaz[3]+0.9682458365518543*alphavpar[2])*fin[3]; 
  out[19] += 0.9682458365518543*(alphavpar[5]*fin[15]+alphavpar[2]*fin[14]+alphavpar[1]*fin[13]+fin[11]*alphavpar[12]+alphavpar[0]*fin[10]+fin[7]*alphavpar[9]+fin[6]*alphavpar[8]+fin[3]*alphavpar[4]); 
  out[20] += 0.4330127018922193*alphaz[2]*fin[20]+0.4330127018922194*(alphaz[5]*fin[18]+alphaz[0]*fin[17])+0.4330127018922193*alphaz[1]*fin[16]+0.9682458365518543*(alphavpar[4]*fin[15]+alphavpar[8]*fin[14]+alphavpar[9]*fin[13]+fin[10]*alphavpar[12])+(0.3872983346207416*alphaz[7]+0.9682458365518543*alphavpar[0])*fin[11]+fin[7]*(0.3872983346207416*alphaz[11]+0.9682458365518543*alphavpar[1])+(0.3872983346207416*alphaz[3]+0.9682458365518543*alphavpar[2])*fin[6]+fin[3]*(0.3872983346207416*alphaz[6]+0.9682458365518543*alphavpar[5]); 
  out[21] += 0.9682458365518543*(alphavpar[2]*fin[15]+alphavpar[5]*fin[14]+alphavpar[0]*fin[13]+fin[7]*alphavpar[12]+alphavpar[9]*fin[11]+alphavpar[1]*fin[10]+fin[3]*alphavpar[8]+alphavpar[4]*fin[6]); 
  out[22] += 0.4330127018922194*alphaz[5]*fin[23]+0.4330127018922193*(alphaz[2]*fin[22]+alphaz[1]*fin[21])+0.4330127018922194*alphaz[0]*fin[19]+(0.3872983346207416*alphaz[11]+0.9682458365518543*alphavpar[1])*fin[15]+(0.3872983346207416*alphaz[7]+0.9682458365518543*alphavpar[0])*fin[14]+0.3872983346207416*alphaz[6]*fin[13]+0.9682458365518543*(alphavpar[5]*fin[13]+fin[6]*alphavpar[12]+alphavpar[8]*fin[11])+0.3872983346207416*alphaz[3]*fin[10]+0.9682458365518543*(alphavpar[2]*fin[10]+fin[3]*alphavpar[9]+alphavpar[4]*fin[7]); 
  out[23] += 0.4330127018922193*alphaz[2]*fin[23]+0.4330127018922194*(alphaz[5]*fin[22]+alphaz[0]*fin[21])+0.4330127018922193*alphaz[1]*fin[19]+(0.3872983346207417*alphaz[7]+0.9682458365518543*alphavpar[0])*fin[15]+(0.3872983346207417*alphaz[11]+0.9682458365518543*alphavpar[1])*fin[14]+0.3872983346207417*alphaz[3]*fin[13]+0.9682458365518543*(alphavpar[2]*fin[13]+fin[3]*alphavpar[12]+alphavpar[4]*fin[11])+0.3872983346207417*alphaz[6]*fin[10]+0.9682458365518543*(alphavpar[5]*fin[10]+fin[6]*alphavpar[9]+fin[7]*alphavpar[8]); 

  return 0.; 
} 
