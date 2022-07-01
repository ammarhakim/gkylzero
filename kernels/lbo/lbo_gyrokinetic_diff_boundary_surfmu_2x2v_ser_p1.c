#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_diff_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // m_:           species mass.
  // bmag_inv:     1/(magnetic field magnitude). 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fskin/edge:   Distribution function in cells 
  // out:          Incremented distribution function in cell 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double facDiff[4]; 
  // Expand diffusion coefficient in conf basis.
  facDiff[0] = (bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0])*m_; 
  facDiff[1] = (bmag_inv[2]*nuVtSqSum[3]+nuVtSqSum[2]*bmag_inv[3]+bmag_inv[0]*nuVtSqSum[1]+nuVtSqSum[0]*bmag_inv[1])*m_; 
  facDiff[2] = (bmag_inv[1]*nuVtSqSum[3]+nuVtSqSum[1]*bmag_inv[3]+bmag_inv[0]*nuVtSqSum[2]+nuVtSqSum[0]*bmag_inv[2])*m_; 
  facDiff[3] = (bmag_inv[0]*nuVtSqSum[3]+nuVtSqSum[0]*bmag_inv[3]+bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2])*m_; 

  double vol_incr[16] = {0.0};
  vol_incr[4] = 0.4330127018922193*dxv[3]*facDiff[3]*fskin[5]*rdvSq4+0.4330127018922193*facDiff[2]*fskin[2]*dxv[3]*rdvSq4+0.4330127018922193*facDiff[1]*fskin[1]*dxv[3]*rdvSq4+0.4330127018922193*facDiff[0]*fskin[0]*dxv[3]*rdvSq4; 
  vol_incr[8] = 0.4330127018922193*facDiff[2]*dxv[3]*fskin[5]*rdvSq4+0.4330127018922193*fskin[2]*dxv[3]*facDiff[3]*rdvSq4+0.4330127018922193*facDiff[0]*fskin[1]*dxv[3]*rdvSq4+0.4330127018922193*fskin[0]*facDiff[1]*dxv[3]*rdvSq4; 
  vol_incr[9] = 0.4330127018922193*facDiff[1]*dxv[3]*fskin[5]*rdvSq4+0.4330127018922193*fskin[1]*dxv[3]*facDiff[3]*rdvSq4+0.4330127018922193*facDiff[0]*fskin[2]*dxv[3]*rdvSq4+0.4330127018922193*fskin[0]*facDiff[2]*dxv[3]*rdvSq4; 
  vol_incr[10] = 0.4330127018922193*dxv[3]*facDiff[3]*fskin[11]*rdvSq4+0.4330127018922193*facDiff[2]*dxv[3]*fskin[7]*rdvSq4+0.4330127018922193*facDiff[1]*dxv[3]*fskin[6]*rdvSq4+0.4330127018922193*facDiff[0]*dxv[3]*fskin[3]*rdvSq4; 
  vol_incr[12] = 0.4330127018922193*facDiff[0]*dxv[3]*fskin[5]*rdvSq4+0.4330127018922193*fskin[0]*dxv[3]*facDiff[3]*rdvSq4+0.4330127018922193*facDiff[1]*fskin[2]*dxv[3]*rdvSq4+0.4330127018922193*fskin[1]*facDiff[2]*dxv[3]*rdvSq4; 
  vol_incr[13] = 0.4330127018922193*facDiff[2]*dxv[3]*fskin[11]*rdvSq4+0.4330127018922193*dxv[3]*facDiff[3]*fskin[7]*rdvSq4+0.4330127018922193*facDiff[0]*dxv[3]*fskin[6]*rdvSq4+0.4330127018922193*facDiff[1]*dxv[3]*fskin[3]*rdvSq4; 
  vol_incr[14] = 0.4330127018922193*facDiff[1]*dxv[3]*fskin[11]*rdvSq4+0.4330127018922193*facDiff[0]*dxv[3]*fskin[7]*rdvSq4+0.4330127018922193*dxv[3]*facDiff[3]*fskin[6]*rdvSq4+0.4330127018922193*facDiff[2]*dxv[3]*fskin[3]*rdvSq4; 
  vol_incr[15] = 0.4330127018922193*facDiff[0]*dxv[3]*fskin[11]*rdvSq4+0.4330127018922193*facDiff[1]*dxv[3]*fskin[7]*rdvSq4+0.4330127018922193*facDiff[2]*dxv[3]*fskin[6]*rdvSq4+0.4330127018922193*dxv[3]*facDiff[3]*fskin[3]*rdvSq4; 

  double temp_diff[16] = {0.0}; 
  double temp_edge[16] = {0.0}; 
  double diff_incr[16] = {0.0}; 
  double edge_incr[16] = {0.0}; 

  if (edge == -1) { 

  temp_diff[0] = (-0.5412658773652741*w[3]*fskin[4])-0.270632938682637*dxv[3]*fskin[4]-0.5412658773652741*w[3]*fedge[4]-0.270632938682637*dxv[3]*fedge[4]-0.5625*fskin[0]*w[3]+0.5625*fedge[0]*w[3]-0.28125*fskin[0]*dxv[3]+0.28125*fedge[0]*dxv[3]; 
  temp_diff[1] = (-0.5412658773652741*w[3]*fskin[8])-0.270632938682637*dxv[3]*fskin[8]-0.5412658773652741*w[3]*fedge[8]-0.270632938682637*dxv[3]*fedge[8]-0.5625*fskin[1]*w[3]+0.5625*fedge[1]*w[3]-0.28125*fskin[1]*dxv[3]+0.28125*fedge[1]*dxv[3]; 
  temp_diff[2] = (-0.5412658773652741*w[3]*fskin[9])-0.270632938682637*dxv[3]*fskin[9]-0.5412658773652741*w[3]*fedge[9]-0.270632938682637*dxv[3]*fedge[9]-0.5625*fskin[2]*w[3]+0.5625*fedge[2]*w[3]-0.28125*fskin[2]*dxv[3]+0.28125*fedge[2]*dxv[3]; 
  temp_diff[3] = (-0.5412658773652741*w[3]*fskin[10])-0.270632938682637*dxv[3]*fskin[10]-0.5412658773652741*w[3]*fedge[10]-0.270632938682637*dxv[3]*fedge[10]-0.5625*fskin[3]*w[3]+0.5625*fedge[3]*w[3]-0.28125*dxv[3]*fskin[3]+0.28125*dxv[3]*fedge[3]; 
  temp_diff[4] = (-1.4375*w[3]*fskin[4])-0.71875*dxv[3]*fskin[4]-0.4375*w[3]*fedge[4]-0.21875*dxv[3]*fedge[4]-1.407291281149713*fskin[0]*w[3]+0.5412658773652741*fedge[0]*w[3]-0.7036456405748563*fskin[0]*dxv[3]+0.270632938682637*fedge[0]*dxv[3]; 
  temp_diff[5] = (-0.5412658773652741*w[3]*fskin[12])-0.270632938682637*dxv[3]*fskin[12]-0.5412658773652741*w[3]*fedge[12]-0.270632938682637*dxv[3]*fedge[12]-0.5625*w[3]*fskin[5]-0.28125*dxv[3]*fskin[5]+0.5625*w[3]*fedge[5]+0.28125*dxv[3]*fedge[5]; 
  temp_diff[6] = (-0.5412658773652741*w[3]*fskin[13])-0.270632938682637*dxv[3]*fskin[13]-0.5412658773652741*w[3]*fedge[13]-0.270632938682637*dxv[3]*fedge[13]-0.5625*w[3]*fskin[6]-0.28125*dxv[3]*fskin[6]+0.5625*w[3]*fedge[6]+0.28125*dxv[3]*fedge[6]; 
  temp_diff[7] = (-0.5412658773652741*w[3]*fskin[14])-0.270632938682637*dxv[3]*fskin[14]-0.5412658773652741*w[3]*fedge[14]-0.270632938682637*dxv[3]*fedge[14]-0.5625*w[3]*fskin[7]-0.28125*dxv[3]*fskin[7]+0.5625*w[3]*fedge[7]+0.28125*dxv[3]*fedge[7]; 
  temp_diff[8] = (-1.4375*w[3]*fskin[8])-0.71875*dxv[3]*fskin[8]-0.4375*w[3]*fedge[8]-0.21875*dxv[3]*fedge[8]-1.407291281149713*fskin[1]*w[3]+0.5412658773652741*fedge[1]*w[3]-0.7036456405748563*fskin[1]*dxv[3]+0.270632938682637*fedge[1]*dxv[3]; 
  temp_diff[9] = (-1.4375*w[3]*fskin[9])-0.71875*dxv[3]*fskin[9]-0.4375*w[3]*fedge[9]-0.21875*dxv[3]*fedge[9]-1.407291281149713*fskin[2]*w[3]+0.5412658773652741*fedge[2]*w[3]-0.7036456405748563*fskin[2]*dxv[3]+0.270632938682637*fedge[2]*dxv[3]; 
  temp_diff[10] = (-1.4375*w[3]*fskin[10])-0.71875*dxv[3]*fskin[10]-0.4375*w[3]*fedge[10]-0.21875*dxv[3]*fedge[10]-1.407291281149713*fskin[3]*w[3]+0.5412658773652741*fedge[3]*w[3]-0.7036456405748563*dxv[3]*fskin[3]+0.270632938682637*dxv[3]*fedge[3]; 
  temp_diff[11] = (-0.5412658773652741*w[3]*fskin[15])-0.270632938682637*dxv[3]*fskin[15]-0.5412658773652741*w[3]*fedge[15]-0.270632938682637*dxv[3]*fedge[15]-0.5625*w[3]*fskin[11]-0.28125*dxv[3]*fskin[11]+0.5625*w[3]*fedge[11]+0.28125*dxv[3]*fedge[11]; 
  temp_diff[12] = (-1.4375*w[3]*fskin[12])-0.71875*dxv[3]*fskin[12]-0.4375*w[3]*fedge[12]-0.21875*dxv[3]*fedge[12]-1.407291281149713*w[3]*fskin[5]-0.7036456405748563*dxv[3]*fskin[5]+0.5412658773652741*w[3]*fedge[5]+0.270632938682637*dxv[3]*fedge[5]; 
  temp_diff[13] = (-1.4375*w[3]*fskin[13])-0.71875*dxv[3]*fskin[13]-0.4375*w[3]*fedge[13]-0.21875*dxv[3]*fedge[13]-1.407291281149713*w[3]*fskin[6]-0.7036456405748563*dxv[3]*fskin[6]+0.5412658773652741*w[3]*fedge[6]+0.270632938682637*dxv[3]*fedge[6]; 
  temp_diff[14] = (-1.4375*w[3]*fskin[14])-0.71875*dxv[3]*fskin[14]-0.4375*w[3]*fedge[14]-0.21875*dxv[3]*fedge[14]-1.407291281149713*w[3]*fskin[7]-0.7036456405748563*dxv[3]*fskin[7]+0.5412658773652741*w[3]*fedge[7]+0.270632938682637*dxv[3]*fedge[7]; 
  temp_diff[15] = (-1.4375*w[3]*fskin[15])-0.71875*dxv[3]*fskin[15]-0.4375*w[3]*fedge[15]-0.21875*dxv[3]*fedge[15]-1.407291281149713*w[3]*fskin[11]-0.7036456405748563*dxv[3]*fskin[11]+0.5412658773652741*w[3]*fedge[11]+0.270632938682637*dxv[3]*fedge[11]; 

  temp_edge[4] = (-1.5*w[3]*fskin[4])+0.75*dxv[3]*fskin[4]+0.8660254037844386*fskin[0]*w[3]-0.4330127018922193*fskin[0]*dxv[3]; 
  temp_edge[8] = (-1.5*w[3]*fskin[8])+0.75*dxv[3]*fskin[8]+0.8660254037844386*fskin[1]*w[3]-0.4330127018922193*fskin[1]*dxv[3]; 
  temp_edge[9] = (-1.5*w[3]*fskin[9])+0.75*dxv[3]*fskin[9]+0.8660254037844386*fskin[2]*w[3]-0.4330127018922193*fskin[2]*dxv[3]; 
  temp_edge[10] = (-1.5*w[3]*fskin[10])+0.75*dxv[3]*fskin[10]+0.8660254037844386*fskin[3]*w[3]-0.4330127018922193*dxv[3]*fskin[3]; 
  temp_edge[12] = (-1.5*w[3]*fskin[12])+0.75*dxv[3]*fskin[12]+0.8660254037844386*w[3]*fskin[5]-0.4330127018922193*dxv[3]*fskin[5]; 
  temp_edge[13] = (-1.5*w[3]*fskin[13])+0.75*dxv[3]*fskin[13]+0.8660254037844386*w[3]*fskin[6]-0.4330127018922193*dxv[3]*fskin[6]; 
  temp_edge[14] = (-1.5*w[3]*fskin[14])+0.75*dxv[3]*fskin[14]+0.8660254037844386*w[3]*fskin[7]-0.4330127018922193*dxv[3]*fskin[7]; 
  temp_edge[15] = (-1.5*w[3]*fskin[15])+0.75*dxv[3]*fskin[15]+0.8660254037844386*w[3]*fskin[11]-0.4330127018922193*dxv[3]*fskin[11]; 

  diff_incr[0] = 0.5*facDiff[3]*temp_diff[5]+0.5*facDiff[2]*temp_diff[2]+0.5*facDiff[1]*temp_diff[1]+0.5*facDiff[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*facDiff[2]*temp_diff[5]+0.5*temp_diff[2]*facDiff[3]+0.5*facDiff[0]*temp_diff[1]+0.5*temp_diff[0]*facDiff[1]; 
  diff_incr[2] = 0.5*facDiff[1]*temp_diff[5]+0.5*temp_diff[1]*facDiff[3]+0.5*facDiff[0]*temp_diff[2]+0.5*temp_diff[0]*facDiff[2]; 
  diff_incr[3] = 0.5*facDiff[3]*temp_diff[11]+0.5*facDiff[2]*temp_diff[7]+0.5*facDiff[1]*temp_diff[6]+0.5*facDiff[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*facDiff[3]*temp_diff[12]+0.5*facDiff[2]*temp_diff[9]+0.5*facDiff[1]*temp_diff[8]+0.5*facDiff[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*facDiff[0]*temp_diff[5]+0.5*temp_diff[0]*facDiff[3]+0.5*facDiff[1]*temp_diff[2]+0.5*temp_diff[1]*facDiff[2]; 
  diff_incr[6] = 0.5*facDiff[2]*temp_diff[11]+0.5*facDiff[3]*temp_diff[7]+0.5*facDiff[0]*temp_diff[6]+0.5*facDiff[1]*temp_diff[3]; 
  diff_incr[7] = 0.5*facDiff[1]*temp_diff[11]+0.5*facDiff[0]*temp_diff[7]+0.5*facDiff[3]*temp_diff[6]+0.5*facDiff[2]*temp_diff[3]; 
  diff_incr[8] = 0.5*facDiff[2]*temp_diff[12]+0.5*facDiff[3]*temp_diff[9]+0.5*facDiff[0]*temp_diff[8]+0.5*facDiff[1]*temp_diff[4]; 
  diff_incr[9] = 0.5*facDiff[1]*temp_diff[12]+0.5*facDiff[0]*temp_diff[9]+0.5*facDiff[3]*temp_diff[8]+0.5*facDiff[2]*temp_diff[4]; 
  diff_incr[10] = 0.5*facDiff[3]*temp_diff[15]+0.5*facDiff[2]*temp_diff[14]+0.5*facDiff[1]*temp_diff[13]+0.5*facDiff[0]*temp_diff[10]; 
  diff_incr[11] = 0.5*facDiff[0]*temp_diff[11]+0.5*facDiff[1]*temp_diff[7]+0.5*facDiff[2]*temp_diff[6]+0.5*facDiff[3]*temp_diff[3]; 
  diff_incr[12] = 0.5*facDiff[0]*temp_diff[12]+0.5*facDiff[1]*temp_diff[9]+0.5*facDiff[2]*temp_diff[8]+0.5*facDiff[3]*temp_diff[4]; 
  diff_incr[13] = 0.5*facDiff[2]*temp_diff[15]+0.5*facDiff[3]*temp_diff[14]+0.5*facDiff[0]*temp_diff[13]+0.5*facDiff[1]*temp_diff[10]; 
  diff_incr[14] = 0.5*facDiff[1]*temp_diff[15]+0.5*facDiff[0]*temp_diff[14]+0.5*facDiff[3]*temp_diff[13]+0.5*facDiff[2]*temp_diff[10]; 
  diff_incr[15] = 0.5*facDiff[0]*temp_diff[15]+0.5*facDiff[1]*temp_diff[14]+0.5*facDiff[2]*temp_diff[13]+0.5*facDiff[3]*temp_diff[10]; 

  edge_incr[0] = 0.5*facDiff[3]*temp_edge[5]+0.5*facDiff[2]*temp_edge[2]+0.5*facDiff[1]*temp_edge[1]+0.5*facDiff[0]*temp_edge[0]; 
  edge_incr[1] = 0.5*facDiff[2]*temp_edge[5]+0.5*temp_edge[2]*facDiff[3]+0.5*facDiff[0]*temp_edge[1]+0.5*temp_edge[0]*facDiff[1]; 
  edge_incr[2] = 0.5*facDiff[1]*temp_edge[5]+0.5*temp_edge[1]*facDiff[3]+0.5*facDiff[0]*temp_edge[2]+0.5*temp_edge[0]*facDiff[2]; 
  edge_incr[3] = 0.5*facDiff[3]*temp_edge[11]+0.5*facDiff[2]*temp_edge[7]+0.5*facDiff[1]*temp_edge[6]+0.5*facDiff[0]*temp_edge[3]; 
  edge_incr[4] = 0.5*facDiff[3]*temp_edge[12]+0.5*facDiff[2]*temp_edge[9]+0.5*facDiff[1]*temp_edge[8]+0.5*facDiff[0]*temp_edge[4]; 
  edge_incr[5] = 0.5*facDiff[0]*temp_edge[5]+0.5*temp_edge[0]*facDiff[3]+0.5*facDiff[1]*temp_edge[2]+0.5*temp_edge[1]*facDiff[2]; 
  edge_incr[6] = 0.5*facDiff[2]*temp_edge[11]+0.5*facDiff[3]*temp_edge[7]+0.5*facDiff[0]*temp_edge[6]+0.5*facDiff[1]*temp_edge[3]; 
  edge_incr[7] = 0.5*facDiff[1]*temp_edge[11]+0.5*facDiff[0]*temp_edge[7]+0.5*facDiff[3]*temp_edge[6]+0.5*facDiff[2]*temp_edge[3]; 
  edge_incr[8] = 0.5*facDiff[2]*temp_edge[12]+0.5*facDiff[3]*temp_edge[9]+0.5*facDiff[0]*temp_edge[8]+0.5*facDiff[1]*temp_edge[4]; 
  edge_incr[9] = 0.5*facDiff[1]*temp_edge[12]+0.5*facDiff[0]*temp_edge[9]+0.5*facDiff[3]*temp_edge[8]+0.5*facDiff[2]*temp_edge[4]; 
  edge_incr[10] = 0.5*facDiff[3]*temp_edge[15]+0.5*facDiff[2]*temp_edge[14]+0.5*facDiff[1]*temp_edge[13]+0.5*facDiff[0]*temp_edge[10]; 
  edge_incr[11] = 0.5*facDiff[0]*temp_edge[11]+0.5*facDiff[1]*temp_edge[7]+0.5*facDiff[2]*temp_edge[6]+0.5*facDiff[3]*temp_edge[3]; 
  edge_incr[12] = 0.5*facDiff[0]*temp_edge[12]+0.5*facDiff[1]*temp_edge[9]+0.5*facDiff[2]*temp_edge[8]+0.5*facDiff[3]*temp_edge[4]; 
  edge_incr[13] = 0.5*facDiff[2]*temp_edge[15]+0.5*facDiff[3]*temp_edge[14]+0.5*facDiff[0]*temp_edge[13]+0.5*facDiff[1]*temp_edge[10]; 
  edge_incr[14] = 0.5*facDiff[1]*temp_edge[15]+0.5*facDiff[0]*temp_edge[14]+0.5*facDiff[3]*temp_edge[13]+0.5*facDiff[2]*temp_edge[10]; 
  edge_incr[15] = 0.5*facDiff[0]*temp_edge[15]+0.5*facDiff[1]*temp_edge[14]+0.5*facDiff[2]*temp_edge[13]+0.5*facDiff[3]*temp_edge[10]; 


  } else { 

  temp_diff[0] = 0.5412658773652741*w[3]*fskin[4]-0.270632938682637*dxv[3]*fskin[4]+0.5412658773652741*w[3]*fedge[4]-0.270632938682637*dxv[3]*fedge[4]-0.5625*fskin[0]*w[3]+0.5625*fedge[0]*w[3]+0.28125*fskin[0]*dxv[3]-0.28125*fedge[0]*dxv[3]; 
  temp_diff[1] = 0.5412658773652741*w[3]*fskin[8]-0.270632938682637*dxv[3]*fskin[8]+0.5412658773652741*w[3]*fedge[8]-0.270632938682637*dxv[3]*fedge[8]-0.5625*fskin[1]*w[3]+0.5625*fedge[1]*w[3]+0.28125*fskin[1]*dxv[3]-0.28125*fedge[1]*dxv[3]; 
  temp_diff[2] = 0.5412658773652741*w[3]*fskin[9]-0.270632938682637*dxv[3]*fskin[9]+0.5412658773652741*w[3]*fedge[9]-0.270632938682637*dxv[3]*fedge[9]-0.5625*fskin[2]*w[3]+0.5625*fedge[2]*w[3]+0.28125*fskin[2]*dxv[3]-0.28125*fedge[2]*dxv[3]; 
  temp_diff[3] = 0.5412658773652741*w[3]*fskin[10]-0.270632938682637*dxv[3]*fskin[10]+0.5412658773652741*w[3]*fedge[10]-0.270632938682637*dxv[3]*fedge[10]-0.5625*fskin[3]*w[3]+0.5625*fedge[3]*w[3]+0.28125*dxv[3]*fskin[3]-0.28125*dxv[3]*fedge[3]; 
  temp_diff[4] = (-1.4375*w[3]*fskin[4])+0.71875*dxv[3]*fskin[4]-0.4375*w[3]*fedge[4]+0.21875*dxv[3]*fedge[4]+1.407291281149713*fskin[0]*w[3]-0.5412658773652741*fedge[0]*w[3]-0.7036456405748563*fskin[0]*dxv[3]+0.270632938682637*fedge[0]*dxv[3]; 
  temp_diff[5] = 0.5412658773652741*w[3]*fskin[12]-0.270632938682637*dxv[3]*fskin[12]+0.5412658773652741*w[3]*fedge[12]-0.270632938682637*dxv[3]*fedge[12]-0.5625*w[3]*fskin[5]+0.28125*dxv[3]*fskin[5]+0.5625*w[3]*fedge[5]-0.28125*dxv[3]*fedge[5]; 
  temp_diff[6] = 0.5412658773652741*w[3]*fskin[13]-0.270632938682637*dxv[3]*fskin[13]+0.5412658773652741*w[3]*fedge[13]-0.270632938682637*dxv[3]*fedge[13]-0.5625*w[3]*fskin[6]+0.28125*dxv[3]*fskin[6]+0.5625*w[3]*fedge[6]-0.28125*dxv[3]*fedge[6]; 
  temp_diff[7] = 0.5412658773652741*w[3]*fskin[14]-0.270632938682637*dxv[3]*fskin[14]+0.5412658773652741*w[3]*fedge[14]-0.270632938682637*dxv[3]*fedge[14]-0.5625*w[3]*fskin[7]+0.28125*dxv[3]*fskin[7]+0.5625*w[3]*fedge[7]-0.28125*dxv[3]*fedge[7]; 
  temp_diff[8] = (-1.4375*w[3]*fskin[8])+0.71875*dxv[3]*fskin[8]-0.4375*w[3]*fedge[8]+0.21875*dxv[3]*fedge[8]+1.407291281149713*fskin[1]*w[3]-0.5412658773652741*fedge[1]*w[3]-0.7036456405748563*fskin[1]*dxv[3]+0.270632938682637*fedge[1]*dxv[3]; 
  temp_diff[9] = (-1.4375*w[3]*fskin[9])+0.71875*dxv[3]*fskin[9]-0.4375*w[3]*fedge[9]+0.21875*dxv[3]*fedge[9]+1.407291281149713*fskin[2]*w[3]-0.5412658773652741*fedge[2]*w[3]-0.7036456405748563*fskin[2]*dxv[3]+0.270632938682637*fedge[2]*dxv[3]; 
  temp_diff[10] = (-1.4375*w[3]*fskin[10])+0.71875*dxv[3]*fskin[10]-0.4375*w[3]*fedge[10]+0.21875*dxv[3]*fedge[10]+1.407291281149713*fskin[3]*w[3]-0.5412658773652741*fedge[3]*w[3]-0.7036456405748563*dxv[3]*fskin[3]+0.270632938682637*dxv[3]*fedge[3]; 
  temp_diff[11] = 0.5412658773652741*w[3]*fskin[15]-0.270632938682637*dxv[3]*fskin[15]+0.5412658773652741*w[3]*fedge[15]-0.270632938682637*dxv[3]*fedge[15]-0.5625*w[3]*fskin[11]+0.28125*dxv[3]*fskin[11]+0.5625*w[3]*fedge[11]-0.28125*dxv[3]*fedge[11]; 
  temp_diff[12] = (-1.4375*w[3]*fskin[12])+0.71875*dxv[3]*fskin[12]-0.4375*w[3]*fedge[12]+0.21875*dxv[3]*fedge[12]+1.407291281149713*w[3]*fskin[5]-0.7036456405748563*dxv[3]*fskin[5]-0.5412658773652741*w[3]*fedge[5]+0.270632938682637*dxv[3]*fedge[5]; 
  temp_diff[13] = (-1.4375*w[3]*fskin[13])+0.71875*dxv[3]*fskin[13]-0.4375*w[3]*fedge[13]+0.21875*dxv[3]*fedge[13]+1.407291281149713*w[3]*fskin[6]-0.7036456405748563*dxv[3]*fskin[6]-0.5412658773652741*w[3]*fedge[6]+0.270632938682637*dxv[3]*fedge[6]; 
  temp_diff[14] = (-1.4375*w[3]*fskin[14])+0.71875*dxv[3]*fskin[14]-0.4375*w[3]*fedge[14]+0.21875*dxv[3]*fedge[14]+1.407291281149713*w[3]*fskin[7]-0.7036456405748563*dxv[3]*fskin[7]-0.5412658773652741*w[3]*fedge[7]+0.270632938682637*dxv[3]*fedge[7]; 
  temp_diff[15] = (-1.4375*w[3]*fskin[15])+0.71875*dxv[3]*fskin[15]-0.4375*w[3]*fedge[15]+0.21875*dxv[3]*fedge[15]+1.407291281149713*w[3]*fskin[11]-0.7036456405748563*dxv[3]*fskin[11]-0.5412658773652741*w[3]*fedge[11]+0.270632938682637*dxv[3]*fedge[11]; 

  temp_edge[4] = (-1.5*w[3]*fskin[4])-0.75*dxv[3]*fskin[4]-0.8660254037844386*fskin[0]*w[3]-0.4330127018922193*fskin[0]*dxv[3]; 
  temp_edge[8] = (-1.5*w[3]*fskin[8])-0.75*dxv[3]*fskin[8]-0.8660254037844386*fskin[1]*w[3]-0.4330127018922193*fskin[1]*dxv[3]; 
  temp_edge[9] = (-1.5*w[3]*fskin[9])-0.75*dxv[3]*fskin[9]-0.8660254037844386*fskin[2]*w[3]-0.4330127018922193*fskin[2]*dxv[3]; 
  temp_edge[10] = (-1.5*w[3]*fskin[10])-0.75*dxv[3]*fskin[10]-0.8660254037844386*fskin[3]*w[3]-0.4330127018922193*dxv[3]*fskin[3]; 
  temp_edge[12] = (-1.5*w[3]*fskin[12])-0.75*dxv[3]*fskin[12]-0.8660254037844386*w[3]*fskin[5]-0.4330127018922193*dxv[3]*fskin[5]; 
  temp_edge[13] = (-1.5*w[3]*fskin[13])-0.75*dxv[3]*fskin[13]-0.8660254037844386*w[3]*fskin[6]-0.4330127018922193*dxv[3]*fskin[6]; 
  temp_edge[14] = (-1.5*w[3]*fskin[14])-0.75*dxv[3]*fskin[14]-0.8660254037844386*w[3]*fskin[7]-0.4330127018922193*dxv[3]*fskin[7]; 
  temp_edge[15] = (-1.5*w[3]*fskin[15])-0.75*dxv[3]*fskin[15]-0.8660254037844386*w[3]*fskin[11]-0.4330127018922193*dxv[3]*fskin[11]; 

  diff_incr[0] = 0.5*facDiff[3]*temp_diff[5]+0.5*facDiff[2]*temp_diff[2]+0.5*facDiff[1]*temp_diff[1]+0.5*facDiff[0]*temp_diff[0]; 
  diff_incr[1] = 0.5*facDiff[2]*temp_diff[5]+0.5*temp_diff[2]*facDiff[3]+0.5*facDiff[0]*temp_diff[1]+0.5*temp_diff[0]*facDiff[1]; 
  diff_incr[2] = 0.5*facDiff[1]*temp_diff[5]+0.5*temp_diff[1]*facDiff[3]+0.5*facDiff[0]*temp_diff[2]+0.5*temp_diff[0]*facDiff[2]; 
  diff_incr[3] = 0.5*facDiff[3]*temp_diff[11]+0.5*facDiff[2]*temp_diff[7]+0.5*facDiff[1]*temp_diff[6]+0.5*facDiff[0]*temp_diff[3]; 
  diff_incr[4] = 0.5*facDiff[3]*temp_diff[12]+0.5*facDiff[2]*temp_diff[9]+0.5*facDiff[1]*temp_diff[8]+0.5*facDiff[0]*temp_diff[4]; 
  diff_incr[5] = 0.5*facDiff[0]*temp_diff[5]+0.5*temp_diff[0]*facDiff[3]+0.5*facDiff[1]*temp_diff[2]+0.5*temp_diff[1]*facDiff[2]; 
  diff_incr[6] = 0.5*facDiff[2]*temp_diff[11]+0.5*facDiff[3]*temp_diff[7]+0.5*facDiff[0]*temp_diff[6]+0.5*facDiff[1]*temp_diff[3]; 
  diff_incr[7] = 0.5*facDiff[1]*temp_diff[11]+0.5*facDiff[0]*temp_diff[7]+0.5*facDiff[3]*temp_diff[6]+0.5*facDiff[2]*temp_diff[3]; 
  diff_incr[8] = 0.5*facDiff[2]*temp_diff[12]+0.5*facDiff[3]*temp_diff[9]+0.5*facDiff[0]*temp_diff[8]+0.5*facDiff[1]*temp_diff[4]; 
  diff_incr[9] = 0.5*facDiff[1]*temp_diff[12]+0.5*facDiff[0]*temp_diff[9]+0.5*facDiff[3]*temp_diff[8]+0.5*facDiff[2]*temp_diff[4]; 
  diff_incr[10] = 0.5*facDiff[3]*temp_diff[15]+0.5*facDiff[2]*temp_diff[14]+0.5*facDiff[1]*temp_diff[13]+0.5*facDiff[0]*temp_diff[10]; 
  diff_incr[11] = 0.5*facDiff[0]*temp_diff[11]+0.5*facDiff[1]*temp_diff[7]+0.5*facDiff[2]*temp_diff[6]+0.5*facDiff[3]*temp_diff[3]; 
  diff_incr[12] = 0.5*facDiff[0]*temp_diff[12]+0.5*facDiff[1]*temp_diff[9]+0.5*facDiff[2]*temp_diff[8]+0.5*facDiff[3]*temp_diff[4]; 
  diff_incr[13] = 0.5*facDiff[2]*temp_diff[15]+0.5*facDiff[3]*temp_diff[14]+0.5*facDiff[0]*temp_diff[13]+0.5*facDiff[1]*temp_diff[10]; 
  diff_incr[14] = 0.5*facDiff[1]*temp_diff[15]+0.5*facDiff[0]*temp_diff[14]+0.5*facDiff[3]*temp_diff[13]+0.5*facDiff[2]*temp_diff[10]; 
  diff_incr[15] = 0.5*facDiff[0]*temp_diff[15]+0.5*facDiff[1]*temp_diff[14]+0.5*facDiff[2]*temp_diff[13]+0.5*facDiff[3]*temp_diff[10]; 

  edge_incr[0] = 0.5*facDiff[3]*temp_edge[5]+0.5*facDiff[2]*temp_edge[2]+0.5*facDiff[1]*temp_edge[1]+0.5*facDiff[0]*temp_edge[0]; 
  edge_incr[1] = 0.5*facDiff[2]*temp_edge[5]+0.5*temp_edge[2]*facDiff[3]+0.5*facDiff[0]*temp_edge[1]+0.5*temp_edge[0]*facDiff[1]; 
  edge_incr[2] = 0.5*facDiff[1]*temp_edge[5]+0.5*temp_edge[1]*facDiff[3]+0.5*facDiff[0]*temp_edge[2]+0.5*temp_edge[0]*facDiff[2]; 
  edge_incr[3] = 0.5*facDiff[3]*temp_edge[11]+0.5*facDiff[2]*temp_edge[7]+0.5*facDiff[1]*temp_edge[6]+0.5*facDiff[0]*temp_edge[3]; 
  edge_incr[4] = 0.5*facDiff[3]*temp_edge[12]+0.5*facDiff[2]*temp_edge[9]+0.5*facDiff[1]*temp_edge[8]+0.5*facDiff[0]*temp_edge[4]; 
  edge_incr[5] = 0.5*facDiff[0]*temp_edge[5]+0.5*temp_edge[0]*facDiff[3]+0.5*facDiff[1]*temp_edge[2]+0.5*temp_edge[1]*facDiff[2]; 
  edge_incr[6] = 0.5*facDiff[2]*temp_edge[11]+0.5*facDiff[3]*temp_edge[7]+0.5*facDiff[0]*temp_edge[6]+0.5*facDiff[1]*temp_edge[3]; 
  edge_incr[7] = 0.5*facDiff[1]*temp_edge[11]+0.5*facDiff[0]*temp_edge[7]+0.5*facDiff[3]*temp_edge[6]+0.5*facDiff[2]*temp_edge[3]; 
  edge_incr[8] = 0.5*facDiff[2]*temp_edge[12]+0.5*facDiff[3]*temp_edge[9]+0.5*facDiff[0]*temp_edge[8]+0.5*facDiff[1]*temp_edge[4]; 
  edge_incr[9] = 0.5*facDiff[1]*temp_edge[12]+0.5*facDiff[0]*temp_edge[9]+0.5*facDiff[3]*temp_edge[8]+0.5*facDiff[2]*temp_edge[4]; 
  edge_incr[10] = 0.5*facDiff[3]*temp_edge[15]+0.5*facDiff[2]*temp_edge[14]+0.5*facDiff[1]*temp_edge[13]+0.5*facDiff[0]*temp_edge[10]; 
  edge_incr[11] = 0.5*facDiff[0]*temp_edge[11]+0.5*facDiff[1]*temp_edge[7]+0.5*facDiff[2]*temp_edge[6]+0.5*facDiff[3]*temp_edge[3]; 
  edge_incr[12] = 0.5*facDiff[0]*temp_edge[12]+0.5*facDiff[1]*temp_edge[9]+0.5*facDiff[2]*temp_edge[8]+0.5*facDiff[3]*temp_edge[4]; 
  edge_incr[13] = 0.5*facDiff[2]*temp_edge[15]+0.5*facDiff[3]*temp_edge[14]+0.5*facDiff[0]*temp_edge[13]+0.5*facDiff[1]*temp_edge[10]; 
  edge_incr[14] = 0.5*facDiff[1]*temp_edge[15]+0.5*facDiff[0]*temp_edge[14]+0.5*facDiff[3]*temp_edge[13]+0.5*facDiff[2]*temp_edge[10]; 
  edge_incr[15] = 0.5*facDiff[0]*temp_edge[15]+0.5*facDiff[1]*temp_edge[14]+0.5*facDiff[2]*temp_edge[13]+0.5*facDiff[3]*temp_edge[10]; 

  } 

  out[0] += edge_incr[0]*rdvSq4+diff_incr[0]*rdvSq4+vol_incr[0]; 
  out[1] += edge_incr[1]*rdvSq4+diff_incr[1]*rdvSq4+vol_incr[1]; 
  out[2] += edge_incr[2]*rdvSq4+diff_incr[2]*rdvSq4+vol_incr[2]; 
  out[3] += edge_incr[3]*rdvSq4+diff_incr[3]*rdvSq4+vol_incr[3]; 
  out[4] += edge_incr[4]*rdvSq4+diff_incr[4]*rdvSq4+vol_incr[4]; 
  out[5] += edge_incr[5]*rdvSq4+diff_incr[5]*rdvSq4+vol_incr[5]; 
  out[6] += edge_incr[6]*rdvSq4+diff_incr[6]*rdvSq4+vol_incr[6]; 
  out[7] += edge_incr[7]*rdvSq4+diff_incr[7]*rdvSq4+vol_incr[7]; 
  out[8] += edge_incr[8]*rdvSq4+diff_incr[8]*rdvSq4+vol_incr[8]; 
  out[9] += edge_incr[9]*rdvSq4+diff_incr[9]*rdvSq4+vol_incr[9]; 
  out[10] += edge_incr[10]*rdvSq4+diff_incr[10]*rdvSq4+vol_incr[10]; 
  out[11] += edge_incr[11]*rdvSq4+diff_incr[11]*rdvSq4+vol_incr[11]; 
  out[12] += edge_incr[12]*rdvSq4+diff_incr[12]*rdvSq4+vol_incr[12]; 
  out[13] += edge_incr[13]*rdvSq4+diff_incr[13]*rdvSq4+vol_incr[13]; 
  out[14] += edge_incr[14]*rdvSq4+diff_incr[14]*rdvSq4+vol_incr[14]; 
  out[15] += edge_incr[15]*rdvSq4+diff_incr[15]*rdvSq4+vol_incr[15]; 
} 
