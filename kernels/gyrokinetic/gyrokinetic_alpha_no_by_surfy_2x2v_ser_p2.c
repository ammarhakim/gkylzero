#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv,
    const double *cmag, const double *b_i, const double *phi, double* GKYL_RESTRICT alpha_surf,
    double* GKYL_RESTRICT sgn_alpha_surf) 
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
  // phi: electrostatic potential.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+1.414213562373095*(vmapSq[0]*m_+bmag[0]*vmap[2]); 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_+1.414213562373095*bmag[2]*vmap[2]; 
  hamil[3] = 1.414213562373095*vmapSq[1]*m_; 
  hamil[4] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 2.0*phi[3]*q_+1.414213562373095*vmap[2]*bmag[3]; 
  hamil[8] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[9] = 1.414213562373095*bmag[2]*vmap[3]; 
  hamil[11] = 2.0*phi[4]*q_+1.414213562373095*vmap[2]*bmag[4]; 
  hamil[12] = 2.0*phi[5]*q_+1.414213562373095*vmap[2]*bmag[5]; 
  hamil[13] = 1.414213562373095*vmapSq[2]*m_; 
  hamil[16] = 1.414213562373095*bmag[3]*vmap[3]; 
  hamil[19] = 2.0*phi[6]*q_+1.414213562373095*vmap[2]*bmag[6]; 
  hamil[20] = 2.0*phi[7]*q_+1.414213562373095*vmap[2]*bmag[7]; 
  hamil[25] = 1.414213562373095*vmap[3]*bmag[4]; 
  hamil[26] = 1.414213562373095*vmap[3]*bmag[5]; 
  hamil[35] = 1.414213562373095*vmap[3]*bmag[6]; 
  hamil[36] = 1.414213562373095*vmap[3]*bmag[7]; 

  double *alphaL = &alpha_surf[20];
  double *sgn_alpha_surfL = &sgn_alpha_surf[27];
  alphaL[0] = (1.25*hamil[3]*cmag[7]*jacobtot_inv[7]-0.9682458365518543*cmag[3]*hamil[3]*jacobtot_inv[7]+0.5590169943749476*cmag[1]*hamil[3]*jacobtot_inv[7]-0.9682458365518543*hamil[3]*jacobtot_inv[3]*cmag[7]+0.5590169943749476*jacobtot_inv[1]*hamil[3]*cmag[7]+0.75*hamil[3]*cmag[6]*jacobtot_inv[6]-0.4330127018922194*hamil[3]*cmag[4]*jacobtot_inv[6]-0.4330127018922194*hamil[3]*jacobtot_inv[4]*cmag[6]+1.25*hamil[3]*cmag[5]*jacobtot_inv[5]-0.9682458365518543*cmag[2]*hamil[3]*jacobtot_inv[5]+0.5590169943749475*cmag[0]*hamil[3]*jacobtot_inv[5]-0.9682458365518543*jacobtot_inv[2]*hamil[3]*cmag[5]+0.5590169943749475*jacobtot_inv[0]*hamil[3]*cmag[5]+0.25*hamil[3]*cmag[4]*jacobtot_inv[4]+0.75*cmag[3]*hamil[3]*jacobtot_inv[3]-0.4330127018922193*cmag[1]*hamil[3]*jacobtot_inv[3]-0.4330127018922193*jacobtot_inv[1]*cmag[3]*hamil[3]+0.75*cmag[2]*jacobtot_inv[2]*hamil[3]-0.4330127018922193*cmag[0]*jacobtot_inv[2]*hamil[3]-0.4330127018922193*jacobtot_inv[0]*cmag[2]*hamil[3]+0.25*cmag[1]*jacobtot_inv[1]*hamil[3]+0.25*cmag[0]*jacobtot_inv[0]*hamil[3])/(vmap[1]*m_); 
  alphaL[1] = ((-0.8660254037844386*hamil[3]*cmag[6]*jacobtot_inv[7])+1.25*hamil[3]*cmag[5]*jacobtot_inv[7]+0.5000000000000001*hamil[3]*cmag[4]*jacobtot_inv[7]-0.9682458365518543*cmag[2]*hamil[3]*jacobtot_inv[7]+0.5590169943749476*cmag[0]*hamil[3]*jacobtot_inv[7]-0.8660254037844386*hamil[3]*jacobtot_inv[6]*cmag[7]+1.25*hamil[3]*jacobtot_inv[5]*cmag[7]+0.5000000000000001*hamil[3]*jacobtot_inv[4]*cmag[7]-0.9682458365518543*jacobtot_inv[2]*hamil[3]*cmag[7]+0.5590169943749476*jacobtot_inv[0]*hamil[3]*cmag[7]+0.6708203932499369*cmag[3]*hamil[3]*jacobtot_inv[6]-0.3872983346207417*cmag[1]*hamil[3]*jacobtot_inv[6]+0.6708203932499369*hamil[3]*jacobtot_inv[3]*cmag[6]-0.3872983346207417*jacobtot_inv[1]*hamil[3]*cmag[6]-0.9682458365518543*cmag[3]*hamil[3]*jacobtot_inv[5]+0.5590169943749475*cmag[1]*hamil[3]*jacobtot_inv[5]-0.9682458365518543*hamil[3]*jacobtot_inv[3]*cmag[5]+0.5590169943749475*jacobtot_inv[1]*hamil[3]*cmag[5]-0.3872983346207416*cmag[3]*hamil[3]*jacobtot_inv[4]+0.223606797749979*cmag[1]*hamil[3]*jacobtot_inv[4]-0.3872983346207416*hamil[3]*jacobtot_inv[3]*cmag[4]+0.223606797749979*jacobtot_inv[1]*hamil[3]*cmag[4]+0.75*cmag[2]*hamil[3]*jacobtot_inv[3]-0.4330127018922193*cmag[0]*hamil[3]*jacobtot_inv[3]+0.75*jacobtot_inv[2]*cmag[3]*hamil[3]-0.4330127018922193*jacobtot_inv[0]*cmag[3]*hamil[3]-0.4330127018922193*cmag[1]*jacobtot_inv[2]*hamil[3]-0.4330127018922193*jacobtot_inv[1]*cmag[2]*hamil[3]+0.25*cmag[0]*jacobtot_inv[1]*hamil[3]+0.25*jacobtot_inv[0]*cmag[1]*hamil[3])/(vmap[1]*m_); 
  alphaL[2] = (2.795084971874738*cmag[7]*jacobtot_inv[7]*hamil[13]-2.165063509461097*cmag[3]*jacobtot_inv[7]*hamil[13]+1.25*cmag[1]*jacobtot_inv[7]*hamil[13]-2.165063509461097*jacobtot_inv[3]*cmag[7]*hamil[13]+1.25*jacobtot_inv[1]*cmag[7]*hamil[13]+1.677050983124842*cmag[6]*jacobtot_inv[6]*hamil[13]-0.9682458365518543*cmag[4]*jacobtot_inv[6]*hamil[13]-0.9682458365518543*jacobtot_inv[4]*cmag[6]*hamil[13]+2.795084971874738*cmag[5]*jacobtot_inv[5]*hamil[13]-2.165063509461096*cmag[2]*jacobtot_inv[5]*hamil[13]+1.25*cmag[0]*jacobtot_inv[5]*hamil[13]-2.165063509461096*jacobtot_inv[2]*cmag[5]*hamil[13]+1.25*jacobtot_inv[0]*cmag[5]*hamil[13]+0.5590169943749475*cmag[4]*jacobtot_inv[4]*hamil[13]+1.677050983124842*cmag[3]*jacobtot_inv[3]*hamil[13]-0.9682458365518543*cmag[1]*jacobtot_inv[3]*hamil[13]-0.9682458365518543*jacobtot_inv[1]*cmag[3]*hamil[13]+1.677050983124842*cmag[2]*jacobtot_inv[2]*hamil[13]-0.9682458365518543*cmag[0]*jacobtot_inv[2]*hamil[13]-0.9682458365518543*jacobtot_inv[0]*cmag[2]*hamil[13]+0.5590169943749475*cmag[1]*jacobtot_inv[1]*hamil[13]+0.5590169943749475*cmag[0]*jacobtot_inv[0]*hamil[13])/(vmap[1]*m_); 
  alphaL[4] = ((-1.936491673103709*cmag[6]*jacobtot_inv[7]*hamil[13])+2.795084971874738*cmag[5]*jacobtot_inv[7]*hamil[13]+1.118033988749895*cmag[4]*jacobtot_inv[7]*hamil[13]-2.165063509461097*cmag[2]*jacobtot_inv[7]*hamil[13]+1.25*cmag[0]*jacobtot_inv[7]*hamil[13]-1.936491673103709*jacobtot_inv[6]*cmag[7]*hamil[13]+2.795084971874738*jacobtot_inv[5]*cmag[7]*hamil[13]+1.118033988749895*jacobtot_inv[4]*cmag[7]*hamil[13]-2.165063509461097*jacobtot_inv[2]*cmag[7]*hamil[13]+1.25*jacobtot_inv[0]*cmag[7]*hamil[13]+1.5*cmag[3]*jacobtot_inv[6]*hamil[13]-0.8660254037844387*cmag[1]*jacobtot_inv[6]*hamil[13]+1.5*jacobtot_inv[3]*cmag[6]*hamil[13]-0.8660254037844387*jacobtot_inv[1]*cmag[6]*hamil[13]-2.165063509461096*cmag[3]*jacobtot_inv[5]*hamil[13]+1.25*cmag[1]*jacobtot_inv[5]*hamil[13]-2.165063509461096*jacobtot_inv[3]*cmag[5]*hamil[13]+1.25*jacobtot_inv[1]*cmag[5]*hamil[13]-0.8660254037844386*cmag[3]*jacobtot_inv[4]*hamil[13]+0.5*cmag[1]*jacobtot_inv[4]*hamil[13]-0.8660254037844386*jacobtot_inv[3]*cmag[4]*hamil[13]+0.5*jacobtot_inv[1]*cmag[4]*hamil[13]+1.677050983124842*cmag[2]*jacobtot_inv[3]*hamil[13]-0.9682458365518543*cmag[0]*jacobtot_inv[3]*hamil[13]+1.677050983124842*jacobtot_inv[2]*cmag[3]*hamil[13]-0.9682458365518543*jacobtot_inv[0]*cmag[3]*hamil[13]-0.9682458365518543*cmag[1]*jacobtot_inv[2]*hamil[13]-0.9682458365518543*jacobtot_inv[1]*cmag[2]*hamil[13]+0.5590169943749475*cmag[0]*jacobtot_inv[1]*hamil[13]+0.5590169943749475*jacobtot_inv[0]*cmag[1]*hamil[13])/(vmap[1]*m_); 
  alphaL[7] = (1.118033988749895*hamil[3]*cmag[7]*jacobtot_inv[7]-0.8660254037844387*cmag[3]*hamil[3]*jacobtot_inv[7]+0.5000000000000001*cmag[1]*hamil[3]*jacobtot_inv[7]-0.8660254037844387*hamil[3]*jacobtot_inv[3]*cmag[7]+0.5000000000000001*jacobtot_inv[1]*hamil[3]*cmag[7]+0.479157423749955*hamil[3]*cmag[6]*jacobtot_inv[6]-0.9682458365518543*hamil[3]*cmag[5]*jacobtot_inv[6]-0.276641667586244*hamil[3]*cmag[4]*jacobtot_inv[6]+0.75*cmag[2]*hamil[3]*jacobtot_inv[6]-0.4330127018922194*cmag[0]*hamil[3]*jacobtot_inv[6]-0.9682458365518543*hamil[3]*jacobtot_inv[5]*cmag[6]-0.276641667586244*hamil[3]*jacobtot_inv[4]*cmag[6]+0.75*jacobtot_inv[2]*hamil[3]*cmag[6]-0.4330127018922194*jacobtot_inv[0]*hamil[3]*cmag[6]+0.5590169943749475*hamil[3]*cmag[4]*jacobtot_inv[5]+0.5590169943749475*hamil[3]*jacobtot_inv[4]*cmag[5]+0.159719141249985*hamil[3]*cmag[4]*jacobtot_inv[4]-0.4330127018922193*cmag[2]*hamil[3]*jacobtot_inv[4]+0.25*cmag[0]*hamil[3]*jacobtot_inv[4]-0.4330127018922193*jacobtot_inv[2]*hamil[3]*cmag[4]+0.25*jacobtot_inv[0]*hamil[3]*cmag[4]+0.6708203932499369*cmag[3]*hamil[3]*jacobtot_inv[3]-0.3872983346207416*cmag[1]*hamil[3]*jacobtot_inv[3]-0.3872983346207416*jacobtot_inv[1]*cmag[3]*hamil[3]+0.223606797749979*cmag[1]*jacobtot_inv[1]*hamil[3])/(vmap[1]*m_); 
  alphaL[11] = (2.5*cmag[7]*jacobtot_inv[7]*hamil[13]-1.936491673103709*cmag[3]*jacobtot_inv[7]*hamil[13]+1.118033988749895*cmag[1]*jacobtot_inv[7]*hamil[13]-1.936491673103709*jacobtot_inv[3]*cmag[7]*hamil[13]+1.118033988749895*jacobtot_inv[1]*cmag[7]*hamil[13]+1.071428571428571*cmag[6]*jacobtot_inv[6]*hamil[13]-2.165063509461096*cmag[5]*jacobtot_inv[6]*hamil[13]-0.6185895741317419*cmag[4]*jacobtot_inv[6]*hamil[13]+1.677050983124842*cmag[2]*jacobtot_inv[6]*hamil[13]-0.9682458365518543*cmag[0]*jacobtot_inv[6]*hamil[13]-2.165063509461096*jacobtot_inv[5]*cmag[6]*hamil[13]-0.6185895741317419*jacobtot_inv[4]*cmag[6]*hamil[13]+1.677050983124842*jacobtot_inv[2]*cmag[6]*hamil[13]-0.9682458365518543*jacobtot_inv[0]*cmag[6]*hamil[13]+1.25*cmag[4]*jacobtot_inv[5]*hamil[13]+1.25*jacobtot_inv[4]*cmag[5]*hamil[13]+0.3571428571428572*cmag[4]*jacobtot_inv[4]*hamil[13]-0.9682458365518543*cmag[2]*jacobtot_inv[4]*hamil[13]+0.5590169943749476*cmag[0]*jacobtot_inv[4]*hamil[13]-0.9682458365518543*jacobtot_inv[2]*cmag[4]*hamil[13]+0.5590169943749476*jacobtot_inv[0]*cmag[4]*hamil[13]+1.5*cmag[3]*jacobtot_inv[3]*hamil[13]-0.8660254037844387*cmag[1]*jacobtot_inv[3]*hamil[13]-0.8660254037844387*jacobtot_inv[1]*cmag[3]*hamil[13]+0.5000000000000001*cmag[1]*jacobtot_inv[1]*hamil[13])/(vmap[1]*m_); 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.4242640687119286*alphaL[11])+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if ((-0.4242640687119286*alphaL[11])+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaL[11])+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]-0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaL[7]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaL[7]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaL[7]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[11]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[11]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaL[11]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaL[0]-0.3952847075210471*alphaL[7] > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaL[0]-0.3952847075210471*alphaL[7] > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaL[0]-0.3952847075210471*alphaL[7] > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[11])-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[11])-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaL[11])-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaL[11])+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaL[11])+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaL[11])+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaL[7]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaL[7]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaL[7]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[24] = 1.0; 
  else  
    sgn_alpha_surfL[24] = -1.0; 
  
  if (sgn_alpha_surfL[24] == sgn_alpha_surfL[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[25] = 1.0; 
  else  
    sgn_alpha_surfL[25] = -1.0; 
  
  if (sgn_alpha_surfL[25] == sgn_alpha_surfL[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]+0.4743416490252568*(alphaL[2]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[26] = 1.0; 
  else  
    sgn_alpha_surfL[26] = -1.0; 
  
  if (sgn_alpha_surfL[26] == sgn_alpha_surfL[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
