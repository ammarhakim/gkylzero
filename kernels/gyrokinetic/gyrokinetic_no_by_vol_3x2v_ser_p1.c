#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_no_by_vol_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  double rdy2 = 2.0/dxv[1];
  double rdz2 = 2.0/dxv[2];
  double rdvpar2 = 2.0/dxv[3];
  double rdmu2 = 2.0/dxv[4];

  double rdvpar2Sq = rdvpar2*rdvpar2;
  double dvparSq = dxv[3]*dxv[3];

  const double *bioverJB_x = &bioverJB[0];
  const double *bioverJB_y = &bioverJB[8];
  const double *bioverJB_z = &bioverJB[16];

  const double *dualcurlbhatoverB_x = &dualcurlbhatoverB[0];
  const double *dualcurlbhatoverB_y = &dualcurlbhatoverB[8];
  const double *dualcurlbhatoverB_z = &dualcurlbhatoverB[16];

  double hamil[48] = {0.}; 
  hamil[0] = 2.0*(phi[0]*q_+vmapSq[0]*m_)+1.414213562373095*bmag[0]*vmap[2]; 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*phi[3]*q_+1.414213562373095*vmap[2]*bmag[3]; 
  hamil[4] = 2.0*vmapSq[1]*m_; 
  hamil[5] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*phi[5]*q_+1.414213562373095*vmap[2]*bmag[5]; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[14] = 1.414213562373095*bmag[3]*vmap[3]; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = 1.414213562373095*vmap[3]*bmag[5]; 
  hamil[32] = 2.0*vmapSq[2]*m_; 

  double hamil2[48] = {0.}; 
  hamil2[0] = hamil[0]*hamil[0]; 
  hamil2[1] = hamil[1]*hamil[1]; 
  hamil2[2] = hamil[2]*hamil[2]; 
  hamil2[3] = hamil[3]*hamil[3]; 
  hamil2[4] = hamil[4]*hamil[4]; 
  hamil2[5] = hamil[5]*hamil[5]; 
  hamil2[6] = hamil[6]*hamil[6]; 
  hamil2[7] = hamil[7]*hamil[7]; 
  hamil2[8] = hamil[8]*hamil[8]; 
  hamil2[9] = hamil[9]*hamil[9]; 
  hamil2[10] = hamil[10]*hamil[10]; 
  hamil2[11] = hamil[11]*hamil[11]; 
  hamil2[12] = hamil[12]*hamil[12]; 
  hamil2[13] = hamil[13]*hamil[13]; 
  hamil2[14] = hamil[14]*hamil[14]; 
  hamil2[15] = hamil[15]*hamil[15]; 
  hamil2[16] = hamil[16]*hamil[16]; 
  hamil2[17] = hamil[17]*hamil[17]; 
  hamil2[18] = hamil[18]*hamil[18]; 
  hamil2[19] = hamil[19]*hamil[19]; 
  hamil2[20] = hamil[20]*hamil[20]; 
  hamil2[21] = hamil[21]*hamil[21]; 
  hamil2[22] = hamil[22]*hamil[22]; 
  hamil2[23] = hamil[23]*hamil[23]; 
  hamil2[24] = hamil[24]*hamil[24]; 
  hamil2[25] = hamil[25]*hamil[25]; 
  hamil2[26] = hamil[26]*hamil[26]; 
  hamil2[27] = hamil[27]*hamil[27]; 
  hamil2[28] = hamil[28]*hamil[28]; 
  hamil2[29] = hamil[29]*hamil[29]; 
  hamil2[30] = hamil[30]*hamil[30]; 
  hamil2[31] = hamil[31]*hamil[31]; 
  hamil2[32] = hamil[32]*hamil[32]; 
  hamil2[33] = hamil[33]*hamil[33]; 
  hamil2[34] = hamil[34]*hamil[34]; 
  hamil2[35] = hamil[35]*hamil[35]; 
  hamil2[36] = hamil[36]*hamil[36]; 
  hamil2[37] = hamil[37]*hamil[37]; 
  hamil2[38] = hamil[38]*hamil[38]; 
  hamil2[39] = hamil[39]*hamil[39]; 
  hamil2[40] = hamil[40]*hamil[40]; 
  hamil2[41] = hamil[41]*hamil[41]; 
  hamil2[42] = hamil[42]*hamil[42]; 
  hamil2[43] = hamil[43]*hamil[43]; 
  hamil2[44] = hamil[44]*hamil[44]; 
  hamil2[45] = hamil[45]*hamil[45]; 
  hamil2[46] = hamil[46]*hamil[46]; 
  hamil2[47] = hamil[47]*hamil[47]; 

  double vmap2 = vmap[1]*vmap[1]; 

  double alphax[48] = {0.}; 



  double alphay[48] = {0.}; 



  double alphaz[48] = {0.}; 
  alphaz[0] = (0.5*rtg33inv[0]*hamil[4]*rdz2)/(vmap[1]*m_); 
  alphaz[1] = (0.5*rtg33inv[1]*hamil[4]*rdz2)/(vmap[1]*m_); 
  alphaz[3] = (0.5*rtg33inv[3]*hamil[4]*rdz2)/(vmap[1]*m_); 
  alphaz[4] = (1.118033988749895*rtg33inv[0]*hamil[32]*rdz2)/(vmap[1]*m_); 
  alphaz[7] = (0.5*hamil[4]*rtg33inv[5]*rdz2)/(vmap[1]*m_); 
  alphaz[9] = (1.118033988749895*rtg33inv[1]*hamil[32]*rdz2)/(vmap[1]*m_); 
  alphaz[11] = (1.118033988749895*rtg33inv[3]*hamil[32]*rdz2)/(vmap[1]*m_); 
  alphaz[18] = (1.118033988749895*rtg33inv[5]*hamil[32]*rdz2)/(vmap[1]*m_); 


  out[3] += 0.3061862178478971*alphaz[18]*fin[18]+0.3061862178478971*alphaz[11]*fin[11]+0.3061862178478971*alphaz[9]*fin[9]+0.3061862178478971*alphaz[7]*fin[7]+0.3061862178478971*alphaz[4]*fin[4]+0.3061862178478971*alphaz[3]*fin[3]+0.3061862178478971*alphaz[1]*fin[1]+0.3061862178478971*alphaz[0]*fin[0]; 
  out[7] += 0.3061862178478971*alphaz[11]*fin[18]+0.3061862178478971*fin[11]*alphaz[18]+0.3061862178478971*alphaz[4]*fin[9]+0.3061862178478971*fin[4]*alphaz[9]+0.3061862178478971*alphaz[3]*fin[7]+0.3061862178478971*fin[3]*alphaz[7]+0.3061862178478971*alphaz[0]*fin[1]+0.3061862178478971*fin[0]*alphaz[1]; 
  out[8] += 0.3061862178478971*alphaz[18]*fin[26]+0.3061862178478971*alphaz[11]*fin[19]+0.3061862178478971*alphaz[9]*fin[17]+0.3061862178478971*alphaz[7]*fin[16]+0.3061862178478971*alphaz[4]*fin[10]+0.3061862178478971*alphaz[3]*fin[8]+0.3061862178478971*alphaz[1]*fin[6]+0.3061862178478971*alphaz[0]*fin[2]; 
  out[11] += 0.273861278752583*alphaz[18]*fin[38]+0.2738612787525831*alphaz[11]*fin[35]+0.2738612787525831*alphaz[9]*fin[33]+0.273861278752583*alphaz[4]*fin[32]+0.3061862178478971*alphaz[7]*fin[18]+0.3061862178478971*fin[7]*alphaz[18]+0.3061862178478971*alphaz[3]*fin[11]+0.3061862178478971*fin[3]*alphaz[11]+0.3061862178478971*alphaz[1]*fin[9]+0.3061862178478971*fin[1]*alphaz[9]+0.3061862178478971*alphaz[0]*fin[4]+0.3061862178478971*fin[0]*alphaz[4]; 
  out[14] += 0.3061862178478971*alphaz[18]*fin[29]+0.3061862178478971*alphaz[11]*fin[25]+0.3061862178478971*alphaz[9]*fin[23]+0.3061862178478971*alphaz[7]*fin[21]+0.3061862178478971*alphaz[4]*fin[15]+0.3061862178478971*alphaz[3]*fin[14]+0.3061862178478971*alphaz[1]*fin[12]+0.3061862178478971*alphaz[0]*fin[5]; 
  out[16] += 0.3061862178478971*alphaz[11]*fin[26]+0.3061862178478971*alphaz[18]*fin[19]+0.3061862178478971*alphaz[4]*fin[17]+0.3061862178478971*alphaz[3]*fin[16]+0.3061862178478971*alphaz[9]*fin[10]+0.3061862178478971*alphaz[7]*fin[8]+0.3061862178478971*alphaz[0]*fin[6]+0.3061862178478971*alphaz[1]*fin[2]; 
  out[18] += 0.273861278752583*alphaz[11]*fin[38]+0.2738612787525831*alphaz[18]*fin[35]+0.2738612787525831*alphaz[4]*fin[33]+0.273861278752583*alphaz[9]*fin[32]+0.3061862178478971*alphaz[3]*fin[18]+0.3061862178478971*fin[3]*alphaz[18]+0.3061862178478971*alphaz[7]*fin[11]+0.3061862178478971*fin[7]*alphaz[11]+0.3061862178478971*alphaz[0]*fin[9]+0.3061862178478971*fin[0]*alphaz[9]+0.3061862178478971*alphaz[1]*fin[4]+0.3061862178478971*fin[1]*alphaz[4]; 
  out[19] += 0.2738612787525831*alphaz[18]*fin[43]+0.273861278752583*alphaz[11]*fin[39]+0.273861278752583*alphaz[9]*fin[37]+0.2738612787525831*alphaz[4]*fin[34]+0.3061862178478971*alphaz[7]*fin[26]+0.3061862178478971*alphaz[3]*fin[19]+0.3061862178478971*fin[16]*alphaz[18]+0.3061862178478971*alphaz[1]*fin[17]+0.3061862178478971*fin[8]*alphaz[11]+0.3061862178478971*alphaz[0]*fin[10]+0.3061862178478971*fin[6]*alphaz[9]+0.3061862178478971*fin[2]*alphaz[4]; 
  out[21] += 0.3061862178478971*alphaz[11]*fin[29]+0.3061862178478971*alphaz[18]*fin[25]+0.3061862178478971*alphaz[4]*fin[23]+0.3061862178478971*alphaz[3]*fin[21]+0.3061862178478971*alphaz[9]*fin[15]+0.3061862178478971*alphaz[7]*fin[14]+0.3061862178478971*alphaz[0]*fin[12]+0.3061862178478971*alphaz[1]*fin[5]; 
  out[22] += 0.3061862178478971*alphaz[18]*fin[31]+0.3061862178478971*alphaz[11]*fin[30]+0.3061862178478971*alphaz[9]*fin[28]+0.3061862178478971*alphaz[7]*fin[27]+0.3061862178478971*alphaz[4]*fin[24]+0.3061862178478971*alphaz[3]*fin[22]+0.3061862178478971*alphaz[1]*fin[20]+0.3061862178478971*alphaz[0]*fin[13]; 
  out[25] += 0.2738612787525831*alphaz[18]*fin[45]+0.273861278752583*alphaz[11]*fin[42]+0.273861278752583*alphaz[9]*fin[40]+0.2738612787525831*alphaz[4]*fin[36]+0.3061862178478971*alphaz[7]*fin[29]+0.3061862178478971*alphaz[3]*fin[25]+0.3061862178478971*alphaz[1]*fin[23]+0.3061862178478971*alphaz[18]*fin[21]+0.3061862178478971*alphaz[0]*fin[15]+0.3061862178478971*alphaz[11]*fin[14]+0.3061862178478971*alphaz[9]*fin[12]+0.3061862178478971*alphaz[4]*fin[5]; 
  out[26] += 0.2738612787525831*alphaz[11]*fin[43]+0.273861278752583*alphaz[18]*fin[39]+0.273861278752583*alphaz[4]*fin[37]+0.2738612787525831*alphaz[9]*fin[34]+0.3061862178478971*alphaz[3]*fin[26]+0.3061862178478971*alphaz[7]*fin[19]+0.3061862178478971*fin[8]*alphaz[18]+0.3061862178478971*alphaz[0]*fin[17]+0.3061862178478971*alphaz[11]*fin[16]+0.3061862178478971*alphaz[1]*fin[10]+0.3061862178478971*fin[2]*alphaz[9]+0.3061862178478971*alphaz[4]*fin[6]; 
  out[27] += 0.3061862178478971*alphaz[11]*fin[31]+0.3061862178478971*alphaz[18]*fin[30]+0.3061862178478971*alphaz[4]*fin[28]+0.3061862178478971*alphaz[3]*fin[27]+0.3061862178478971*alphaz[9]*fin[24]+0.3061862178478971*alphaz[7]*fin[22]+0.3061862178478971*alphaz[0]*fin[20]+0.3061862178478971*alphaz[1]*fin[13]; 
  out[29] += 0.2738612787525831*alphaz[11]*fin[45]+0.273861278752583*alphaz[18]*fin[42]+0.273861278752583*alphaz[4]*fin[40]+0.2738612787525831*alphaz[9]*fin[36]+0.3061862178478971*alphaz[3]*fin[29]+0.3061862178478971*alphaz[7]*fin[25]+0.3061862178478971*alphaz[0]*fin[23]+0.3061862178478971*alphaz[11]*fin[21]+0.3061862178478971*fin[14]*alphaz[18]+0.3061862178478971*alphaz[1]*fin[15]+0.3061862178478971*alphaz[4]*fin[12]+0.3061862178478971*fin[5]*alphaz[9]; 
  out[30] += 0.273861278752583*alphaz[18]*fin[47]+0.2738612787525831*alphaz[11]*fin[46]+0.2738612787525831*alphaz[9]*fin[44]+0.273861278752583*alphaz[4]*fin[41]+0.3061862178478971*alphaz[7]*fin[31]+0.3061862178478971*alphaz[3]*fin[30]+0.3061862178478971*alphaz[1]*fin[28]+0.3061862178478971*alphaz[18]*fin[27]+0.3061862178478971*alphaz[0]*fin[24]+0.3061862178478971*alphaz[11]*fin[22]+0.3061862178478971*alphaz[9]*fin[20]+0.3061862178478971*alphaz[4]*fin[13]; 
  out[31] += 0.273861278752583*alphaz[11]*fin[47]+0.2738612787525831*alphaz[18]*fin[46]+0.2738612787525831*alphaz[4]*fin[44]+0.273861278752583*alphaz[9]*fin[41]+0.3061862178478971*alphaz[3]*fin[31]+0.3061862178478971*alphaz[7]*fin[30]+0.3061862178478971*alphaz[0]*fin[28]+0.3061862178478971*alphaz[11]*fin[27]+0.3061862178478971*alphaz[1]*fin[24]+0.3061862178478971*alphaz[18]*fin[22]+0.3061862178478971*alphaz[4]*fin[20]+0.3061862178478971*alphaz[9]*fin[13]; 
  out[35] += 0.3061862178478971*alphaz[7]*fin[38]+0.3061862178478971*alphaz[3]*fin[35]+0.3061862178478971*alphaz[1]*fin[33]+0.3061862178478971*alphaz[0]*fin[32]+0.2738612787525831*alphaz[18]*fin[18]+0.2738612787525831*alphaz[11]*fin[11]+0.2738612787525831*alphaz[9]*fin[9]+0.2738612787525831*alphaz[4]*fin[4]; 
  out[38] += 0.3061862178478971*alphaz[3]*fin[38]+0.3061862178478971*alphaz[7]*fin[35]+0.3061862178478971*alphaz[0]*fin[33]+0.3061862178478971*alphaz[1]*fin[32]+0.273861278752583*alphaz[11]*fin[18]+0.273861278752583*fin[11]*alphaz[18]+0.273861278752583*alphaz[4]*fin[9]+0.273861278752583*fin[4]*alphaz[9]; 
  out[39] += 0.3061862178478971*alphaz[7]*fin[43]+0.3061862178478971*alphaz[3]*fin[39]+0.3061862178478971*alphaz[1]*fin[37]+0.3061862178478971*alphaz[0]*fin[34]+0.273861278752583*alphaz[18]*fin[26]+0.273861278752583*alphaz[11]*fin[19]+0.273861278752583*alphaz[9]*fin[17]+0.273861278752583*alphaz[4]*fin[10]; 
  out[42] += 0.3061862178478971*alphaz[7]*fin[45]+0.3061862178478971*alphaz[3]*fin[42]+0.3061862178478971*alphaz[1]*fin[40]+0.3061862178478971*alphaz[0]*fin[36]+0.273861278752583*alphaz[18]*fin[29]+0.273861278752583*alphaz[11]*fin[25]+0.273861278752583*alphaz[9]*fin[23]+0.273861278752583*alphaz[4]*fin[15]; 
  out[43] += 0.3061862178478971*alphaz[3]*fin[43]+0.3061862178478971*alphaz[7]*fin[39]+0.3061862178478971*alphaz[0]*fin[37]+0.3061862178478971*alphaz[1]*fin[34]+0.2738612787525831*alphaz[11]*fin[26]+0.2738612787525831*alphaz[18]*fin[19]+0.2738612787525831*alphaz[4]*fin[17]+0.2738612787525831*alphaz[9]*fin[10]; 
  out[45] += 0.3061862178478971*alphaz[3]*fin[45]+0.3061862178478971*alphaz[7]*fin[42]+0.3061862178478971*alphaz[0]*fin[40]+0.3061862178478971*alphaz[1]*fin[36]+0.2738612787525831*alphaz[11]*fin[29]+0.2738612787525831*alphaz[18]*fin[25]+0.2738612787525831*alphaz[4]*fin[23]+0.2738612787525831*alphaz[9]*fin[15]; 
  out[46] += 0.3061862178478971*alphaz[7]*fin[47]+0.3061862178478971*alphaz[3]*fin[46]+0.3061862178478971*alphaz[1]*fin[44]+0.3061862178478971*alphaz[0]*fin[41]+0.2738612787525831*alphaz[18]*fin[31]+0.2738612787525831*alphaz[11]*fin[30]+0.2738612787525831*alphaz[9]*fin[28]+0.2738612787525831*alphaz[4]*fin[24]; 
  out[47] += 0.3061862178478971*alphaz[3]*fin[47]+0.3061862178478971*alphaz[7]*fin[46]+0.3061862178478971*alphaz[0]*fin[44]+0.3061862178478971*alphaz[1]*fin[41]+0.273861278752583*alphaz[11]*fin[31]+0.273861278752583*alphaz[18]*fin[30]+0.273861278752583*alphaz[4]*fin[28]+0.273861278752583*alphaz[9]*fin[24]; 

  double alphavpar[48] = {0.}; 
  alphavpar[0] = (((-0.5*rtg33inv[1]*hamil[7])-0.5*rtg33inv[0]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphavpar[1] = (((-0.5*rtg33inv[0]*hamil[7])-0.5*rtg33inv[1]*hamil[3])*rdz2)/(vmap[1]*m_); 
  alphavpar[2] = (((-0.5*rtg33inv[1]*hamil[16])-0.5*rtg33inv[0]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[3] = (((-0.5*rtg33inv[5]*hamil[7])-0.5*hamil[3]*rtg33inv[3])*rdz2)/(vmap[1]*m_); 
  alphavpar[5] = (((-0.5*rtg33inv[1]*hamil[21])-0.5*rtg33inv[0]*hamil[14])*rdz2)/(vmap[1]*m_); 
  alphavpar[6] = (((-0.5*rtg33inv[0]*hamil[16])-0.5*rtg33inv[1]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[7] = (((-0.5*rtg33inv[3]*hamil[7])-0.5*hamil[3]*rtg33inv[5])*rdz2)/(vmap[1]*m_); 
  alphavpar[8] = (((-0.5*rtg33inv[5]*hamil[16])-0.5*rtg33inv[3]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[12] = (((-0.5*rtg33inv[0]*hamil[21])-0.5*rtg33inv[1]*hamil[14])*rdz2)/(vmap[1]*m_); 
  alphavpar[14] = (((-0.5*rtg33inv[5]*hamil[21])-0.5*rtg33inv[3]*hamil[14])*rdz2)/(vmap[1]*m_); 
  alphavpar[16] = (((-0.5*rtg33inv[3]*hamil[16])-0.5*rtg33inv[5]*hamil[8])*rdz2)/(vmap[1]*m_); 
  alphavpar[21] = (((-0.5*rtg33inv[3]*hamil[21])-0.5*rtg33inv[5]*hamil[14])*rdz2)/(vmap[1]*m_); 


  out[4] += 0.3061862178478971*alphavpar[21]*fin[21]+0.3061862178478971*alphavpar[16]*fin[16]+0.3061862178478971*alphavpar[14]*fin[14]+0.3061862178478971*alphavpar[12]*fin[12]+0.3061862178478971*alphavpar[8]*fin[8]+0.3061862178478971*alphavpar[7]*fin[7]+0.3061862178478971*alphavpar[6]*fin[6]+0.3061862178478971*alphavpar[5]*fin[5]+0.3061862178478971*alphavpar[3]*fin[3]+0.3061862178478971*alphavpar[2]*fin[2]+0.3061862178478971*alphavpar[1]*fin[1]+0.3061862178478971*alphavpar[0]*fin[0]; 
  out[9] += 0.3061862178478971*alphavpar[14]*fin[21]+0.3061862178478971*fin[14]*alphavpar[21]+0.3061862178478971*alphavpar[8]*fin[16]+0.3061862178478971*fin[8]*alphavpar[16]+0.3061862178478971*alphavpar[5]*fin[12]+0.3061862178478971*fin[5]*alphavpar[12]+0.3061862178478971*alphavpar[3]*fin[7]+0.3061862178478971*fin[3]*alphavpar[7]+0.3061862178478971*alphavpar[2]*fin[6]+0.3061862178478971*fin[2]*alphavpar[6]+0.3061862178478971*alphavpar[0]*fin[1]+0.3061862178478971*fin[0]*alphavpar[1]; 
  out[10] += 0.3061862178478971*alphavpar[21]*fin[27]+0.3061862178478971*alphavpar[14]*fin[22]+0.3061862178478971*alphavpar[12]*fin[20]+0.3061862178478971*alphavpar[7]*fin[16]+0.3061862178478971*fin[7]*alphavpar[16]+0.3061862178478971*alphavpar[5]*fin[13]+0.3061862178478971*alphavpar[3]*fin[8]+0.3061862178478971*fin[3]*alphavpar[8]+0.3061862178478971*alphavpar[1]*fin[6]+0.3061862178478971*fin[1]*alphavpar[6]+0.3061862178478971*alphavpar[0]*fin[2]+0.3061862178478971*fin[0]*alphavpar[2]; 
  out[11] += 0.3061862178478971*alphavpar[12]*fin[21]+0.3061862178478971*fin[12]*alphavpar[21]+0.3061862178478971*alphavpar[6]*fin[16]+0.3061862178478971*fin[6]*alphavpar[16]+0.3061862178478971*alphavpar[5]*fin[14]+0.3061862178478971*fin[5]*alphavpar[14]+0.3061862178478971*alphavpar[2]*fin[8]+0.3061862178478971*fin[2]*alphavpar[8]+0.3061862178478971*alphavpar[1]*fin[7]+0.3061862178478971*fin[1]*alphavpar[7]+0.3061862178478971*alphavpar[0]*fin[3]+0.3061862178478971*fin[0]*alphavpar[3]; 
  out[15] += 0.3061862178478971*alphavpar[16]*fin[27]+0.3061862178478971*alphavpar[8]*fin[22]+0.3061862178478971*alphavpar[7]*fin[21]+0.3061862178478971*fin[7]*alphavpar[21]+0.3061862178478971*alphavpar[6]*fin[20]+0.3061862178478971*alphavpar[3]*fin[14]+0.3061862178478971*fin[3]*alphavpar[14]+0.3061862178478971*alphavpar[2]*fin[13]+0.3061862178478971*alphavpar[1]*fin[12]+0.3061862178478971*fin[1]*alphavpar[12]+0.3061862178478971*alphavpar[0]*fin[5]+0.3061862178478971*fin[0]*alphavpar[5]; 
  out[17] += 0.3061862178478971*alphavpar[14]*fin[27]+0.3061862178478971*alphavpar[21]*fin[22]+0.3061862178478971*alphavpar[5]*fin[20]+0.3061862178478971*alphavpar[3]*fin[16]+0.3061862178478971*fin[3]*alphavpar[16]+0.3061862178478971*alphavpar[12]*fin[13]+0.3061862178478971*alphavpar[7]*fin[8]+0.3061862178478971*fin[7]*alphavpar[8]+0.3061862178478971*alphavpar[0]*fin[6]+0.3061862178478971*fin[0]*alphavpar[6]+0.3061862178478971*alphavpar[1]*fin[2]+0.3061862178478971*fin[1]*alphavpar[2]; 
  out[18] += 0.3061862178478971*alphavpar[5]*fin[21]+0.3061862178478971*fin[5]*alphavpar[21]+0.3061862178478971*alphavpar[2]*fin[16]+0.3061862178478971*fin[2]*alphavpar[16]+0.3061862178478971*alphavpar[12]*fin[14]+0.3061862178478971*fin[12]*alphavpar[14]+0.3061862178478971*alphavpar[6]*fin[8]+0.3061862178478971*fin[6]*alphavpar[8]+0.3061862178478971*alphavpar[0]*fin[7]+0.3061862178478971*fin[0]*alphavpar[7]+0.3061862178478971*alphavpar[1]*fin[3]+0.3061862178478971*fin[1]*alphavpar[3]; 
  out[19] += 0.3061862178478971*alphavpar[12]*fin[27]+0.3061862178478971*alphavpar[5]*fin[22]+0.3061862178478971*fin[20]*alphavpar[21]+0.3061862178478971*alphavpar[1]*fin[16]+0.3061862178478971*fin[1]*alphavpar[16]+0.3061862178478971*fin[13]*alphavpar[14]+0.3061862178478971*alphavpar[0]*fin[8]+0.3061862178478971*fin[0]*alphavpar[8]+0.3061862178478971*alphavpar[6]*fin[7]+0.3061862178478971*fin[6]*alphavpar[7]+0.3061862178478971*alphavpar[2]*fin[3]+0.3061862178478971*fin[2]*alphavpar[3]; 
  out[23] += 0.3061862178478971*alphavpar[8]*fin[27]+0.3061862178478971*alphavpar[16]*fin[22]+0.3061862178478971*alphavpar[3]*fin[21]+0.3061862178478971*fin[3]*alphavpar[21]+0.3061862178478971*alphavpar[2]*fin[20]+0.3061862178478971*alphavpar[7]*fin[14]+0.3061862178478971*fin[7]*alphavpar[14]+0.3061862178478971*alphavpar[6]*fin[13]+0.3061862178478971*alphavpar[0]*fin[12]+0.3061862178478971*fin[0]*alphavpar[12]+0.3061862178478971*alphavpar[1]*fin[5]+0.3061862178478971*fin[1]*alphavpar[5]; 
  out[24] += 0.3061862178478971*alphavpar[7]*fin[27]+0.3061862178478971*alphavpar[3]*fin[22]+0.3061862178478971*alphavpar[16]*fin[21]+0.3061862178478971*fin[16]*alphavpar[21]+0.3061862178478971*alphavpar[1]*fin[20]+0.3061862178478971*alphavpar[8]*fin[14]+0.3061862178478971*fin[8]*alphavpar[14]+0.3061862178478971*alphavpar[0]*fin[13]+0.3061862178478971*alphavpar[6]*fin[12]+0.3061862178478971*fin[6]*alphavpar[12]+0.3061862178478971*alphavpar[2]*fin[5]+0.3061862178478971*fin[2]*alphavpar[5]; 
  out[25] += 0.3061862178478971*alphavpar[6]*fin[27]+0.3061862178478971*alphavpar[2]*fin[22]+0.3061862178478971*alphavpar[1]*fin[21]+0.3061862178478971*fin[1]*alphavpar[21]+0.3061862178478971*alphavpar[16]*fin[20]+0.3061862178478971*alphavpar[0]*fin[14]+0.3061862178478971*fin[0]*alphavpar[14]+0.3061862178478971*alphavpar[8]*fin[13]+0.3061862178478971*alphavpar[7]*fin[12]+0.3061862178478971*fin[7]*alphavpar[12]+0.3061862178478971*alphavpar[3]*fin[5]+0.3061862178478971*fin[3]*alphavpar[5]; 
  out[26] += 0.3061862178478971*alphavpar[5]*fin[27]+0.3061862178478971*alphavpar[12]*fin[22]+0.3061862178478971*fin[13]*alphavpar[21]+0.3061862178478971*alphavpar[14]*fin[20]+0.3061862178478971*alphavpar[0]*fin[16]+0.3061862178478971*fin[0]*alphavpar[16]+0.3061862178478971*alphavpar[1]*fin[8]+0.3061862178478971*fin[1]*alphavpar[8]+0.3061862178478971*alphavpar[2]*fin[7]+0.3061862178478971*fin[2]*alphavpar[7]+0.3061862178478971*alphavpar[3]*fin[6]+0.3061862178478971*fin[3]*alphavpar[6]; 
  out[28] += 0.3061862178478971*alphavpar[3]*fin[27]+0.3061862178478971*alphavpar[7]*fin[22]+0.3061862178478971*alphavpar[8]*fin[21]+0.3061862178478971*fin[8]*alphavpar[21]+0.3061862178478971*alphavpar[0]*fin[20]+0.3061862178478971*alphavpar[14]*fin[16]+0.3061862178478971*fin[14]*alphavpar[16]+0.3061862178478971*alphavpar[1]*fin[13]+0.3061862178478971*alphavpar[2]*fin[12]+0.3061862178478971*fin[2]*alphavpar[12]+0.3061862178478971*alphavpar[5]*fin[6]+0.3061862178478971*fin[5]*alphavpar[6]; 
  out[29] += 0.3061862178478971*alphavpar[2]*fin[27]+0.3061862178478971*alphavpar[6]*fin[22]+0.3061862178478971*alphavpar[0]*fin[21]+0.3061862178478971*fin[0]*alphavpar[21]+0.3061862178478971*alphavpar[8]*fin[20]+0.3061862178478971*fin[13]*alphavpar[16]+0.3061862178478971*alphavpar[1]*fin[14]+0.3061862178478971*fin[1]*alphavpar[14]+0.3061862178478971*alphavpar[3]*fin[12]+0.3061862178478971*fin[3]*alphavpar[12]+0.3061862178478971*alphavpar[5]*fin[7]+0.3061862178478971*fin[5]*alphavpar[7]; 
  out[30] += 0.3061862178478971*alphavpar[1]*fin[27]+0.3061862178478971*alphavpar[0]*fin[22]+0.3061862178478971*alphavpar[6]*fin[21]+0.3061862178478971*fin[6]*alphavpar[21]+0.3061862178478971*alphavpar[7]*fin[20]+0.3061862178478971*alphavpar[12]*fin[16]+0.3061862178478971*fin[12]*alphavpar[16]+0.3061862178478971*alphavpar[2]*fin[14]+0.3061862178478971*fin[2]*alphavpar[14]+0.3061862178478971*alphavpar[3]*fin[13]+0.3061862178478971*alphavpar[5]*fin[8]+0.3061862178478971*fin[5]*alphavpar[8]; 
  out[31] += 0.3061862178478971*alphavpar[0]*fin[27]+0.3061862178478971*alphavpar[1]*fin[22]+0.3061862178478971*alphavpar[2]*fin[21]+0.3061862178478971*fin[2]*alphavpar[21]+0.3061862178478971*alphavpar[3]*fin[20]+0.3061862178478971*alphavpar[5]*fin[16]+0.3061862178478971*fin[5]*alphavpar[16]+0.3061862178478971*alphavpar[6]*fin[14]+0.3061862178478971*fin[6]*alphavpar[14]+0.3061862178478971*alphavpar[7]*fin[13]+0.3061862178478971*alphavpar[8]*fin[12]+0.3061862178478971*fin[8]*alphavpar[12]; 
  out[32] += 0.6846531968814572*alphavpar[21]*fin[29]+0.6846531968814572*alphavpar[16]*fin[26]+0.6846531968814572*alphavpar[14]*fin[25]+0.6846531968814572*alphavpar[12]*fin[23]+0.6846531968814572*alphavpar[8]*fin[19]+0.6846531968814572*alphavpar[7]*fin[18]+0.6846531968814572*alphavpar[6]*fin[17]+0.6846531968814572*alphavpar[5]*fin[15]+0.6846531968814572*alphavpar[3]*fin[11]+0.6846531968814572*alphavpar[2]*fin[10]+0.6846531968814572*alphavpar[1]*fin[9]+0.6846531968814572*alphavpar[0]*fin[4]; 
  out[33] += 0.6846531968814573*alphavpar[14]*fin[29]+0.6846531968814573*alphavpar[8]*fin[26]+0.6846531968814573*alphavpar[21]*fin[25]+0.6846531968814573*alphavpar[5]*fin[23]+0.6846531968814573*alphavpar[16]*fin[19]+0.6846531968814573*alphavpar[3]*fin[18]+0.6846531968814573*alphavpar[2]*fin[17]+0.6846531968814573*alphavpar[12]*fin[15]+0.6846531968814573*alphavpar[7]*fin[11]+0.6846531968814573*alphavpar[6]*fin[10]+0.6846531968814573*alphavpar[0]*fin[9]+0.6846531968814573*alphavpar[1]*fin[4]; 
  out[34] += 0.6846531968814573*alphavpar[21]*fin[31]+0.6846531968814573*alphavpar[14]*fin[30]+0.6846531968814573*alphavpar[12]*fin[28]+0.6846531968814573*alphavpar[7]*fin[26]+0.6846531968814573*alphavpar[5]*fin[24]+0.6846531968814573*alphavpar[3]*fin[19]+0.6846531968814573*alphavpar[16]*fin[18]+0.6846531968814573*alphavpar[1]*fin[17]+0.6846531968814573*alphavpar[8]*fin[11]+0.6846531968814573*alphavpar[0]*fin[10]+0.6846531968814573*alphavpar[6]*fin[9]+0.6846531968814573*alphavpar[2]*fin[4]; 
  out[35] += 0.6846531968814573*alphavpar[12]*fin[29]+0.6846531968814573*alphavpar[6]*fin[26]+0.6846531968814573*alphavpar[5]*fin[25]+0.6846531968814573*alphavpar[21]*fin[23]+0.6846531968814573*alphavpar[2]*fin[19]+0.6846531968814573*alphavpar[1]*fin[18]+0.6846531968814573*alphavpar[16]*fin[17]+0.6846531968814573*alphavpar[14]*fin[15]+0.6846531968814573*alphavpar[0]*fin[11]+0.6846531968814573*alphavpar[8]*fin[10]+0.6846531968814573*alphavpar[7]*fin[9]+0.6846531968814573*alphavpar[3]*fin[4]; 
  out[36] += 0.6846531968814573*alphavpar[16]*fin[31]+0.6846531968814573*alphavpar[8]*fin[30]+0.6846531968814573*alphavpar[7]*fin[29]+0.6846531968814573*alphavpar[6]*fin[28]+0.6846531968814573*alphavpar[3]*fin[25]+0.6846531968814573*alphavpar[2]*fin[24]+0.6846531968814573*alphavpar[1]*fin[23]+0.6846531968814573*fin[18]*alphavpar[21]+0.6846531968814573*alphavpar[0]*fin[15]+0.6846531968814573*fin[11]*alphavpar[14]+0.6846531968814573*fin[9]*alphavpar[12]+0.6846531968814573*fin[4]*alphavpar[5]; 
  out[37] += 0.6846531968814572*alphavpar[14]*fin[31]+0.6846531968814572*alphavpar[21]*fin[30]+0.6846531968814572*alphavpar[5]*fin[28]+0.6846531968814572*alphavpar[3]*fin[26]+0.6846531968814572*alphavpar[12]*fin[24]+0.6846531968814572*alphavpar[7]*fin[19]+0.6846531968814572*alphavpar[8]*fin[18]+0.6846531968814572*alphavpar[0]*fin[17]+0.6846531968814572*fin[11]*alphavpar[16]+0.6846531968814572*alphavpar[1]*fin[10]+0.6846531968814572*alphavpar[2]*fin[9]+0.6846531968814572*fin[4]*alphavpar[6]; 
  out[38] += 0.6846531968814572*alphavpar[5]*fin[29]+0.6846531968814572*alphavpar[2]*fin[26]+0.6846531968814572*alphavpar[12]*fin[25]+0.6846531968814572*alphavpar[14]*fin[23]+0.6846531968814572*fin[15]*alphavpar[21]+0.6846531968814572*alphavpar[6]*fin[19]+0.6846531968814572*alphavpar[0]*fin[18]+0.6846531968814572*alphavpar[8]*fin[17]+0.6846531968814572*fin[10]*alphavpar[16]+0.6846531968814572*alphavpar[1]*fin[11]+0.6846531968814572*alphavpar[3]*fin[9]+0.6846531968814572*fin[4]*alphavpar[7]; 
  out[39] += 0.6846531968814572*alphavpar[12]*fin[31]+0.6846531968814572*alphavpar[5]*fin[30]+0.6846531968814572*alphavpar[21]*fin[28]+0.6846531968814572*alphavpar[1]*fin[26]+0.6846531968814572*alphavpar[14]*fin[24]+0.6846531968814572*alphavpar[0]*fin[19]+0.6846531968814572*alphavpar[6]*fin[18]+0.6846531968814572*alphavpar[7]*fin[17]+0.6846531968814572*fin[9]*alphavpar[16]+0.6846531968814572*alphavpar[2]*fin[11]+0.6846531968814572*alphavpar[3]*fin[10]+0.6846531968814572*fin[4]*alphavpar[8]; 
  out[40] += 0.6846531968814572*alphavpar[8]*fin[31]+0.6846531968814572*alphavpar[16]*fin[30]+0.6846531968814572*alphavpar[3]*fin[29]+0.6846531968814572*alphavpar[2]*fin[28]+0.6846531968814572*alphavpar[7]*fin[25]+0.6846531968814572*alphavpar[6]*fin[24]+0.6846531968814572*alphavpar[0]*fin[23]+0.6846531968814572*fin[11]*alphavpar[21]+0.6846531968814572*alphavpar[14]*fin[18]+0.6846531968814572*alphavpar[1]*fin[15]+0.6846531968814572*fin[4]*alphavpar[12]+0.6846531968814572*alphavpar[5]*fin[9]; 
  out[41] += 0.6846531968814572*alphavpar[7]*fin[31]+0.6846531968814572*alphavpar[3]*fin[30]+0.6846531968814572*alphavpar[16]*fin[29]+0.6846531968814572*alphavpar[1]*fin[28]+0.6846531968814572*alphavpar[21]*fin[26]+0.6846531968814572*alphavpar[8]*fin[25]+0.6846531968814572*alphavpar[0]*fin[24]+0.6846531968814572*alphavpar[6]*fin[23]+0.6846531968814572*alphavpar[14]*fin[19]+0.6846531968814572*alphavpar[12]*fin[17]+0.6846531968814572*alphavpar[2]*fin[15]+0.6846531968814572*alphavpar[5]*fin[10]; 
  out[42] += 0.6846531968814572*alphavpar[6]*fin[31]+0.6846531968814572*alphavpar[2]*fin[30]+0.6846531968814572*alphavpar[1]*fin[29]+0.6846531968814572*alphavpar[16]*fin[28]+0.6846531968814572*alphavpar[0]*fin[25]+0.6846531968814572*alphavpar[8]*fin[24]+0.6846531968814572*alphavpar[7]*fin[23]+0.6846531968814572*fin[9]*alphavpar[21]+0.6846531968814572*alphavpar[12]*fin[18]+0.6846531968814572*alphavpar[3]*fin[15]+0.6846531968814572*fin[4]*alphavpar[14]+0.6846531968814572*alphavpar[5]*fin[11]; 
  out[43] += 0.6846531968814573*alphavpar[5]*fin[31]+0.6846531968814573*alphavpar[12]*fin[30]+0.6846531968814573*alphavpar[14]*fin[28]+0.6846531968814573*alphavpar[0]*fin[26]+0.6846531968814573*alphavpar[21]*fin[24]+0.6846531968814573*alphavpar[1]*fin[19]+0.6846531968814573*alphavpar[2]*fin[18]+0.6846531968814573*alphavpar[3]*fin[17]+0.6846531968814573*fin[4]*alphavpar[16]+0.6846531968814573*alphavpar[6]*fin[11]+0.6846531968814573*alphavpar[7]*fin[10]+0.6846531968814573*alphavpar[8]*fin[9]; 
  out[44] += 0.6846531968814573*alphavpar[3]*fin[31]+0.6846531968814573*alphavpar[7]*fin[30]+0.6846531968814573*alphavpar[8]*fin[29]+0.6846531968814573*alphavpar[0]*fin[28]+0.6846531968814573*alphavpar[14]*fin[26]+0.6846531968814573*alphavpar[16]*fin[25]+0.6846531968814573*alphavpar[1]*fin[24]+0.6846531968814573*alphavpar[2]*fin[23]+0.6846531968814573*fin[19]*alphavpar[21]+0.6846531968814573*alphavpar[5]*fin[17]+0.6846531968814573*alphavpar[6]*fin[15]+0.6846531968814573*fin[10]*alphavpar[12]; 
  out[45] += 0.6846531968814573*alphavpar[2]*fin[31]+0.6846531968814573*alphavpar[6]*fin[30]+0.6846531968814573*alphavpar[0]*fin[29]+0.6846531968814573*alphavpar[8]*fin[28]+0.6846531968814573*alphavpar[1]*fin[25]+0.6846531968814573*alphavpar[16]*fin[24]+0.6846531968814573*alphavpar[3]*fin[23]+0.6846531968814573*fin[4]*alphavpar[21]+0.6846531968814573*alphavpar[5]*fin[18]+0.6846531968814573*alphavpar[7]*fin[15]+0.6846531968814573*fin[9]*alphavpar[14]+0.6846531968814573*fin[11]*alphavpar[12]; 
  out[46] += 0.6846531968814573*alphavpar[1]*fin[31]+0.6846531968814573*alphavpar[0]*fin[30]+0.6846531968814573*alphavpar[6]*fin[29]+0.6846531968814573*alphavpar[7]*fin[28]+0.6846531968814573*alphavpar[12]*fin[26]+0.6846531968814573*alphavpar[2]*fin[25]+0.6846531968814573*alphavpar[3]*fin[24]+0.6846531968814573*alphavpar[16]*fin[23]+0.6846531968814573*fin[17]*alphavpar[21]+0.6846531968814573*alphavpar[5]*fin[19]+0.6846531968814573*alphavpar[8]*fin[15]+0.6846531968814573*fin[10]*alphavpar[14]; 
  out[47] += 0.6846531968814572*alphavpar[0]*fin[31]+0.6846531968814572*alphavpar[1]*fin[30]+0.6846531968814572*alphavpar[2]*fin[29]+0.6846531968814572*alphavpar[3]*fin[28]+0.6846531968814572*alphavpar[5]*fin[26]+0.6846531968814572*alphavpar[6]*fin[25]+0.6846531968814572*alphavpar[7]*fin[24]+0.6846531968814572*alphavpar[8]*fin[23]+0.6846531968814572*fin[10]*alphavpar[21]+0.6846531968814572*alphavpar[12]*fin[19]+0.6846531968814572*alphavpar[14]*fin[17]+0.6846531968814572*fin[15]*alphavpar[16]; 

  return 0.; 
} 
