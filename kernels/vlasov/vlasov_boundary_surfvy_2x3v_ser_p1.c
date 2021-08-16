#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E1 = &qmem[4]; 
  const double *B0 = &qmem[12]; 
  const double *B1 = &qmem[16]; 
  const double *B2 = &qmem[20]; 

  double Ghat[16]; 
  double alpha[16]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = 0.5773502691896258*B0[0]*dv3; 
  alpha[5] = 2.0*B0[3]*wv3-2.0*B2[3]*wv1+2.0*E1[3]; 
  alpha[6] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = -0.5773502691896258*B2[2]*dv1; 
  alpha[8] = 0.5773502691896258*B0[1]*dv3; 
  alpha[9] = 0.5773502691896258*B0[2]*dv3; 
  alpha[11] = -0.5773502691896258*B2[3]*dv1; 
  alpha[12] = 0.5773502691896258*B0[3]*dv3; 

  double fUpwindQuad[16];
  double fUpwind[16];

  if (edge == -1) { 

  if ((-alpha[12])-alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = 0.3061862178478971*fSkin[31]-0.3061862178478971*(fSkin[30]+fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*(fSkin[25]+fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])+0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])-0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6])-0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[0] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*(fEdge[30]+fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*(fEdge[25]+fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])-0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])+0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6])-0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } 
  if (alpha[12]+alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = (-0.3061862178478971*(fSkin[31]+fSkin[30]))+0.3061862178478971*(fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24])-0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*(fSkin[21]+fSkin[20])+0.3061862178478971*fSkin[19]-0.3061862178478971*(fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13])-0.1767766952966368*fSkin[12]-0.3061862178478971*(fSkin[11]+fSkin[10])+0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2])+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[1] = 0.3061862178478971*(fEdge[31]+fEdge[30])-0.3061862178478971*(fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24])+0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*(fEdge[21]+fEdge[20])-0.3061862178478971*fEdge[19]+0.3061862178478971*(fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13])-0.1767766952966368*fEdge[12]+0.3061862178478971*(fEdge[11]+fEdge[10])-0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2])+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*fSkin[30]-0.3061862178478971*fSkin[29]+0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25])-0.3061862178478971*fSkin[24]+0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*fSkin[21]+0.1767766952966368*fSkin[20]-0.3061862178478971*fSkin[19]+0.3061862178478971*fSkin[18]-0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*fSkin[13]+0.1767766952966368*fSkin[12]-0.3061862178478971*fSkin[11]+0.3061862178478971*fSkin[10]-0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*(fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*fSkin[2]-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[2] = 0.3061862178478971*fEdge[31]-0.3061862178478971*fEdge[30]+0.3061862178478971*fEdge[29]-0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25])+0.3061862178478971*fEdge[24]-0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*fEdge[21]+0.1767766952966368*fEdge[20]+0.3061862178478971*fEdge[19]-0.3061862178478971*fEdge[18]+0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*fEdge[13]+0.1767766952966368*fEdge[12]+0.3061862178478971*fEdge[11]-0.3061862178478971*fEdge[10]+0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*(fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*fEdge[2]-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = 0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29])-0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*fSkin[25]-0.3061862178478971*(fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21])-0.1767766952966368*fSkin[20]-0.3061862178478971*(fSkin[19]+fSkin[18])+0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*(fSkin[13]+fSkin[12])-0.3061862178478971*fSkin[11]+0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*fSkin[6]-0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[3] = (-0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29]))+0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*fEdge[25]+0.3061862178478971*(fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21])-0.1767766952966368*fEdge[20]+0.3061862178478971*(fEdge[19]+fEdge[18])-0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*(fEdge[13]+fEdge[12])+0.3061862178478971*fEdge[11]-0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*fEdge[6]-0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*(fSkin[30]+fSkin[29])-0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*fSkin[25]+0.3061862178478971*(fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21])-0.1767766952966368*fSkin[20]-0.3061862178478971*(fSkin[19]+fSkin[18])+0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*(fSkin[13]+fSkin[12])+0.3061862178478971*fSkin[11]-0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*fSkin[6]-0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*(fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[4] = 0.3061862178478971*fEdge[31]-0.3061862178478971*(fEdge[30]+fEdge[29])+0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*fEdge[25]-0.3061862178478971*(fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21])-0.1767766952966368*fEdge[20]+0.3061862178478971*(fEdge[19]+fEdge[18])-0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*(fEdge[13]+fEdge[12])-0.3061862178478971*fEdge[11]+0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*fEdge[6]-0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*(fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } 
  if (alpha[12]-alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = 0.3061862178478971*(fSkin[31]+fSkin[30])-0.3061862178478971*fSkin[29]+0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25])+0.3061862178478971*fSkin[24]-0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*fSkin[21]+0.1767766952966368*fSkin[20]-0.3061862178478971*fSkin[19]+0.3061862178478971*fSkin[18]-0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*fSkin[13]-0.1767766952966368*fSkin[12]+0.3061862178478971*fSkin[11]-0.3061862178478971*fSkin[10]+0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*(fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*fSkin[2]+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[5] = (-0.3061862178478971*(fEdge[31]+fEdge[30]))+0.3061862178478971*fEdge[29]-0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25])-0.3061862178478971*fEdge[24]+0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*fEdge[21]+0.1767766952966368*fEdge[20]+0.3061862178478971*fEdge[19]-0.3061862178478971*fEdge[18]+0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*fEdge[13]-0.1767766952966368*fEdge[12]-0.3061862178478971*fEdge[11]+0.3061862178478971*fEdge[10]-0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*(fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*fEdge[2]+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = 0.3061862178478971*fSkin[31]-0.3061862178478971*fSkin[30]+0.3061862178478971*(fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24])+0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*(fSkin[21]+fSkin[20])+0.3061862178478971*fSkin[19]-0.3061862178478971*(fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13])+0.1767766952966368*fSkin[12]+0.3061862178478971*(fSkin[11]+fSkin[10])-0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2])-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[6] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*fEdge[30]-0.3061862178478971*(fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24])-0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*(fEdge[21]+fEdge[20])-0.3061862178478971*fEdge[19]+0.3061862178478971*(fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13])+0.1767766952966368*fEdge[12]-0.3061862178478971*(fEdge[11]+fEdge[10])+0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2])-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = (-0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29]+fSkin[28]))-0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*(fSkin[25]+fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])+0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])+0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6])-0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[7] = 0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*(fEdge[25]+fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])-0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])-0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6])-0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1]+fEdge[0]); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[8] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*(fSkin[30]+fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])+0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])-0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[8] = 0.3061862178478971*fEdge[31]-0.3061862178478971*(fEdge[30]+fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])-0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])+0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[9] = 0.3061862178478971*(fSkin[31]+fSkin[30])-0.3061862178478971*(fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*(fSkin[25]+fSkin[24])+0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*(fSkin[21]+fSkin[20])+0.3061862178478971*fSkin[19]-0.3061862178478971*(fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13])+0.1767766952966368*fSkin[12]-0.3061862178478971*(fSkin[11]+fSkin[10])+0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6])+0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2])+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[9] = (-0.3061862178478971*(fEdge[31]+fEdge[30]))+0.3061862178478971*(fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*(fEdge[25]+fEdge[24])-0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*(fEdge[21]+fEdge[20])-0.3061862178478971*fEdge[19]+0.3061862178478971*(fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13])+0.1767766952966368*fEdge[12]+0.3061862178478971*(fEdge[11]+fEdge[10])-0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6])+0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2])+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[10] = 0.3061862178478971*fSkin[31]-0.3061862178478971*fSkin[30]+0.3061862178478971*fSkin[29]-0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*fSkin[25]+0.3061862178478971*fSkin[24]-0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*fSkin[21]-0.1767766952966368*fSkin[20]-0.3061862178478971*fSkin[19]+0.3061862178478971*fSkin[18]-0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*fSkin[13]-0.1767766952966368*fSkin[12]-0.3061862178478971*fSkin[11]+0.3061862178478971*fSkin[10]-0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*fSkin[6]+0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*fSkin[2]-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[10] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*fEdge[30]-0.3061862178478971*fEdge[29]+0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*fEdge[25]-0.3061862178478971*fEdge[24]+0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*fEdge[21]-0.1767766952966368*fEdge[20]+0.3061862178478971*fEdge[19]-0.3061862178478971*fEdge[18]+0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*fEdge[13]-0.1767766952966368*fEdge[12]+0.3061862178478971*fEdge[11]-0.3061862178478971*fEdge[10]+0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*fEdge[6]+0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*fEdge[2]-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } 
  if (alpha[12]-alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[11] = (-0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29]))+0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25])+0.3061862178478971*(fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21])+0.1767766952966368*fSkin[20]-0.3061862178478971*(fSkin[19]+fSkin[18])+0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*(fSkin[13]+fSkin[12])-0.3061862178478971*fSkin[11]+0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*(fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[11] = 0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29])-0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25])-0.3061862178478971*(fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21])+0.1767766952966368*fEdge[20]+0.3061862178478971*(fEdge[19]+fEdge[18])-0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*(fEdge[13]+fEdge[12])+0.3061862178478971*fEdge[11]-0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*(fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[12] = 0.3061862178478971*fSkin[31]-0.3061862178478971*(fSkin[30]+fSkin[29])+0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25])-0.3061862178478971*(fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21])+0.1767766952966368*fSkin[20]-0.3061862178478971*(fSkin[19]+fSkin[18])+0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*(fSkin[13]+fSkin[12])+0.3061862178478971*fSkin[11]-0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*(fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*(fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[12] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*(fEdge[30]+fEdge[29])-0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25])+0.3061862178478971*(fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21])+0.1767766952966368*fEdge[20]+0.3061862178478971*(fEdge[19]+fEdge[18])-0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*(fEdge[13]+fEdge[12])-0.3061862178478971*fEdge[11]+0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*(fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*(fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[13] = (-0.3061862178478971*(fSkin[31]+fSkin[30]))+0.3061862178478971*fSkin[29]-0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*fSkin[25]-0.3061862178478971*fSkin[24]+0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*fSkin[21]-0.1767766952966368*fSkin[20]-0.3061862178478971*fSkin[19]+0.3061862178478971*fSkin[18]-0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*fSkin[13]+0.1767766952966368*fSkin[12]+0.3061862178478971*fSkin[11]-0.3061862178478971*fSkin[10]+0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*fSkin[6]+0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*fSkin[2]+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[13] = 0.3061862178478971*(fEdge[31]+fEdge[30])-0.3061862178478971*fEdge[29]+0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*fEdge[25]+0.3061862178478971*fEdge[24]-0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*fEdge[21]-0.1767766952966368*fEdge[20]+0.3061862178478971*fEdge[19]-0.3061862178478971*fEdge[18]+0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*fEdge[13]+0.1767766952966368*fEdge[12]-0.3061862178478971*fEdge[11]+0.3061862178478971*fEdge[10]-0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*fEdge[6]+0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*fEdge[2]+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } 
  if ((-alpha[12])-alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[14] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*fSkin[30]-0.3061862178478971*(fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*(fSkin[25]+fSkin[24])-0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*(fSkin[21]+fSkin[20])+0.3061862178478971*fSkin[19]-0.3061862178478971*(fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13])-0.1767766952966368*fSkin[12]+0.3061862178478971*(fSkin[11]+fSkin[10])-0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6])+0.1767766952966368*fSkin[5]+0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2])-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } else { 

    fUpwindQuad[14] = 0.3061862178478971*fEdge[31]-0.3061862178478971*fEdge[30]+0.3061862178478971*(fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*(fEdge[25]+fEdge[24])+0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*(fEdge[21]+fEdge[20])-0.3061862178478971*fEdge[19]+0.3061862178478971*(fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13])-0.1767766952966368*fEdge[12]-0.3061862178478971*(fEdge[11]+fEdge[10])+0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6])+0.1767766952966368*fEdge[5]-0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2])-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } 
  if (alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[15] = 0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])+0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])+0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6]+fSkin[5])+0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[15] = (-0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29]+fEdge[28]))+0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])-0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])-0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6]+fEdge[5])-0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1]+fEdge[0]); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] = 0.25*(alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.25*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.25*(alpha[12]*fUpwind[15]+alpha[9]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.25*(alpha[11]*fUpwind[15]+alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] = 0.25*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] = 0.25*(alpha[9]*fUpwind[15]+alpha[12]*fUpwind[14]+alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[8]*fUpwind[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] = 0.25*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[12]*fUpwind[13]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] = 0.25*(alpha[7]*fUpwind[15]+alpha[11]*fUpwind[14]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] = 0.25*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[11]*fUpwind[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] = 0.25*(alpha[5]*fUpwind[15]+alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] = 0.25*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[9]*fUpwind[13]+fUpwind[10]*alpha[12]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[12] = 0.25*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+fUpwind[10]*alpha[11]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[13] = 0.25*(alpha[2]*fUpwind[15]+alpha[5]*fUpwind[14]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[14] = 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[5]*fUpwind[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] = 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += -0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -1.224744871391589*Ghat[1]*dv11; 
  out[10] += -1.224744871391589*Ghat[2]*dv11; 
  out[11] += -1.224744871391589*Ghat[3]*dv11; 
  out[12] += -0.7071067811865475*Ghat[8]*dv11; 
  out[13] += -0.7071067811865475*Ghat[9]*dv11; 
  out[14] += -0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -1.224744871391589*Ghat[4]*dv11; 
  out[16] += -0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -1.224744871391589*Ghat[7]*dv11; 
  out[20] += -0.7071067811865475*Ghat[12]*dv11; 
  out[21] += -0.7071067811865475*Ghat[13]*dv11; 
  out[22] += -0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.224744871391589*Ghat[9]*dv11; 
  out[25] += -1.224744871391589*Ghat[10]*dv11; 
  out[26] += -1.224744871391589*Ghat[11]*dv11; 
  out[27] += -0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -1.224744871391589*Ghat[12]*dv11; 
  out[29] += -1.224744871391589*Ghat[13]*dv11; 
  out[30] += -1.224744871391589*Ghat[14]*dv11; 
  out[31] += -1.224744871391589*Ghat[15]*dv11; 

  } else { 

  if ((-alpha[12])-alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = 0.3061862178478971*fEdge[31]-0.3061862178478971*(fEdge[30]+fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*(fEdge[25]+fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])+0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])-0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6])-0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[0] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*(fSkin[30]+fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*(fSkin[25]+fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])-0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])+0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6])-0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } 
  if (alpha[12]+alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = (-0.3061862178478971*(fEdge[31]+fEdge[30]))+0.3061862178478971*(fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24])-0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*(fEdge[21]+fEdge[20])+0.3061862178478971*fEdge[19]-0.3061862178478971*(fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13])-0.1767766952966368*fEdge[12]-0.3061862178478971*(fEdge[11]+fEdge[10])+0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2])+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[1] = 0.3061862178478971*(fSkin[31]+fSkin[30])-0.3061862178478971*(fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24])+0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*(fSkin[21]+fSkin[20])-0.3061862178478971*fSkin[19]+0.3061862178478971*(fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13])-0.1767766952966368*fSkin[12]+0.3061862178478971*(fSkin[11]+fSkin[10])-0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2])+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*fEdge[30]-0.3061862178478971*fEdge[29]+0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25])-0.3061862178478971*fEdge[24]+0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*fEdge[21]+0.1767766952966368*fEdge[20]-0.3061862178478971*fEdge[19]+0.3061862178478971*fEdge[18]-0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*fEdge[13]+0.1767766952966368*fEdge[12]-0.3061862178478971*fEdge[11]+0.3061862178478971*fEdge[10]-0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*(fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*fEdge[2]-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[2] = 0.3061862178478971*fSkin[31]-0.3061862178478971*fSkin[30]+0.3061862178478971*fSkin[29]-0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25])+0.3061862178478971*fSkin[24]-0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*fSkin[21]+0.1767766952966368*fSkin[20]+0.3061862178478971*fSkin[19]-0.3061862178478971*fSkin[18]+0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*fSkin[13]+0.1767766952966368*fSkin[12]+0.3061862178478971*fSkin[11]-0.3061862178478971*fSkin[10]+0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*(fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*fSkin[2]-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = 0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29])-0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*fEdge[25]-0.3061862178478971*(fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21])-0.1767766952966368*fEdge[20]-0.3061862178478971*(fEdge[19]+fEdge[18])+0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*(fEdge[13]+fEdge[12])-0.3061862178478971*fEdge[11]+0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*fEdge[6]-0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[3] = (-0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29]))+0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*fSkin[25]+0.3061862178478971*(fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21])-0.1767766952966368*fSkin[20]+0.3061862178478971*(fSkin[19]+fSkin[18])-0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*(fSkin[13]+fSkin[12])+0.3061862178478971*fSkin[11]-0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*fSkin[6]-0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*(fEdge[30]+fEdge[29])-0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*fEdge[25]+0.3061862178478971*(fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21])-0.1767766952966368*fEdge[20]-0.3061862178478971*(fEdge[19]+fEdge[18])+0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*(fEdge[13]+fEdge[12])+0.3061862178478971*fEdge[11]-0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*fEdge[6]-0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*(fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[4] = 0.3061862178478971*fSkin[31]-0.3061862178478971*(fSkin[30]+fSkin[29])+0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*fSkin[25]-0.3061862178478971*(fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21])-0.1767766952966368*fSkin[20]+0.3061862178478971*(fSkin[19]+fSkin[18])-0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*(fSkin[13]+fSkin[12])-0.3061862178478971*fSkin[11]+0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*fSkin[6]-0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*(fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } 
  if (alpha[12]-alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = 0.3061862178478971*(fEdge[31]+fEdge[30])-0.3061862178478971*fEdge[29]+0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25])+0.3061862178478971*fEdge[24]-0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*fEdge[21]+0.1767766952966368*fEdge[20]-0.3061862178478971*fEdge[19]+0.3061862178478971*fEdge[18]-0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*fEdge[13]-0.1767766952966368*fEdge[12]+0.3061862178478971*fEdge[11]-0.3061862178478971*fEdge[10]+0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*(fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*fEdge[2]+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[5] = (-0.3061862178478971*(fSkin[31]+fSkin[30]))+0.3061862178478971*fSkin[29]-0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25])-0.3061862178478971*fSkin[24]+0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*fSkin[21]+0.1767766952966368*fSkin[20]+0.3061862178478971*fSkin[19]-0.3061862178478971*fSkin[18]+0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*fSkin[13]-0.1767766952966368*fSkin[12]-0.3061862178478971*fSkin[11]+0.3061862178478971*fSkin[10]-0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*(fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*fSkin[2]+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = 0.3061862178478971*fEdge[31]-0.3061862178478971*fEdge[30]+0.3061862178478971*(fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24])+0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*(fEdge[21]+fEdge[20])+0.3061862178478971*fEdge[19]-0.3061862178478971*(fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13])+0.1767766952966368*fEdge[12]+0.3061862178478971*(fEdge[11]+fEdge[10])-0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2])-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[6] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*fSkin[30]-0.3061862178478971*(fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24])-0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*(fSkin[21]+fSkin[20])-0.3061862178478971*fSkin[19]+0.3061862178478971*(fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13])+0.1767766952966368*fSkin[12]-0.3061862178478971*(fSkin[11]+fSkin[10])+0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2])-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = (-0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29]+fEdge[28]))-0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*(fEdge[25]+fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])+0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]-0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])+0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6])-0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[7] = 0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*(fSkin[25]+fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])-0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]+0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])-0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6])-0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1]+fSkin[0]); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[8] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*(fEdge[30]+fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])+0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])-0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[8] = 0.3061862178478971*fSkin[31]-0.3061862178478971*(fSkin[30]+fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])-0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])+0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[9] = 0.3061862178478971*(fEdge[31]+fEdge[30])-0.3061862178478971*(fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*(fEdge[25]+fEdge[24])+0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*(fEdge[21]+fEdge[20])+0.3061862178478971*fEdge[19]-0.3061862178478971*(fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*(fEdge[14]+fEdge[13])+0.1767766952966368*fEdge[12]-0.3061862178478971*(fEdge[11]+fEdge[10])+0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6])+0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]-0.1767766952966368*(fEdge[3]+fEdge[2])+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[9] = (-0.3061862178478971*(fSkin[31]+fSkin[30]))+0.3061862178478971*(fSkin[29]+fSkin[28])+0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*(fSkin[25]+fSkin[24])-0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*(fSkin[21]+fSkin[20])-0.3061862178478971*fSkin[19]+0.3061862178478971*(fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*(fSkin[14]+fSkin[13])+0.1767766952966368*fSkin[12]+0.3061862178478971*(fSkin[11]+fSkin[10])-0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6])+0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]-0.1767766952966368*(fSkin[3]+fSkin[2])+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[10] = 0.3061862178478971*fEdge[31]-0.3061862178478971*fEdge[30]+0.3061862178478971*fEdge[29]-0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]+0.3061862178478971*fEdge[26]-0.3061862178478971*fEdge[25]+0.3061862178478971*fEdge[24]-0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*fEdge[21]-0.1767766952966368*fEdge[20]-0.3061862178478971*fEdge[19]+0.3061862178478971*fEdge[18]-0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*fEdge[13]-0.1767766952966368*fEdge[12]-0.3061862178478971*fEdge[11]+0.3061862178478971*fEdge[10]-0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*fEdge[6]+0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*fEdge[2]-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[10] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*fSkin[30]-0.3061862178478971*fSkin[29]+0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]-0.3061862178478971*fSkin[26]+0.3061862178478971*fSkin[25]-0.3061862178478971*fSkin[24]+0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*fSkin[21]-0.1767766952966368*fSkin[20]+0.3061862178478971*fSkin[19]-0.3061862178478971*fSkin[18]+0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*fSkin[13]-0.1767766952966368*fSkin[12]+0.3061862178478971*fSkin[11]-0.3061862178478971*fSkin[10]+0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*fSkin[6]+0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*fSkin[2]-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } 
  if (alpha[12]-alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[11] = (-0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29]))+0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]-0.3061862178478971*(fEdge[26]+fEdge[25])+0.3061862178478971*(fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21])+0.1767766952966368*fEdge[20]-0.3061862178478971*(fEdge[19]+fEdge[18])+0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]-0.1767766952966368*fEdge[14]+0.1767766952966368*(fEdge[13]+fEdge[12])-0.3061862178478971*fEdge[11]+0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*(fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]-0.1767766952966368*fEdge[3]+0.1767766952966368*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[11] = 0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29])-0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]+0.3061862178478971*(fSkin[26]+fSkin[25])-0.3061862178478971*(fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21])+0.1767766952966368*fSkin[20]+0.3061862178478971*(fSkin[19]+fSkin[18])-0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]-0.1767766952966368*fSkin[14]+0.1767766952966368*(fSkin[13]+fSkin[12])+0.3061862178478971*fSkin[11]-0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*(fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]-0.1767766952966368*fSkin[3]+0.1767766952966368*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[12] = 0.3061862178478971*fEdge[31]-0.3061862178478971*(fEdge[30]+fEdge[29])+0.3061862178478971*fEdge[28]+0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25])-0.3061862178478971*(fEdge[24]+fEdge[23])-0.1767766952966368*(fEdge[22]+fEdge[21])+0.1767766952966368*fEdge[20]-0.3061862178478971*(fEdge[19]+fEdge[18])+0.3061862178478971*fEdge[17]+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*(fEdge[13]+fEdge[12])+0.3061862178478971*fEdge[11]-0.3061862178478971*(fEdge[10]+fEdge[9])-0.1767766952966368*(fEdge[8]+fEdge[7])+0.1767766952966368*(fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*(fEdge[2]+fEdge[1])+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[12] = (-0.3061862178478971*fSkin[31])+0.3061862178478971*(fSkin[30]+fSkin[29])-0.3061862178478971*fSkin[28]+0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25])+0.3061862178478971*(fSkin[24]+fSkin[23])-0.1767766952966368*(fSkin[22]+fSkin[21])+0.1767766952966368*fSkin[20]+0.3061862178478971*(fSkin[19]+fSkin[18])-0.3061862178478971*fSkin[17]+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*(fSkin[13]+fSkin[12])-0.3061862178478971*fSkin[11]+0.3061862178478971*(fSkin[10]+fSkin[9])-0.1767766952966368*(fSkin[8]+fSkin[7])+0.1767766952966368*(fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*(fSkin[2]+fSkin[1])+0.1767766952966368*fSkin[0]; 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[13] = (-0.3061862178478971*(fEdge[31]+fEdge[30]))+0.3061862178478971*fEdge[29]-0.3061862178478971*fEdge[28]-0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*fEdge[25]-0.3061862178478971*fEdge[24]+0.3061862178478971*fEdge[23]-0.1767766952966368*fEdge[22]+0.1767766952966368*fEdge[21]-0.1767766952966368*fEdge[20]-0.3061862178478971*fEdge[19]+0.3061862178478971*fEdge[18]-0.3061862178478971*fEdge[17]-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*fEdge[14]-0.1767766952966368*fEdge[13]+0.1767766952966368*fEdge[12]+0.3061862178478971*fEdge[11]-0.3061862178478971*fEdge[10]+0.3061862178478971*fEdge[9]-0.1767766952966368*fEdge[8]+0.1767766952966368*fEdge[7]-0.1767766952966368*fEdge[6]+0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]+0.1767766952966368*fEdge[3]-0.1767766952966368*fEdge[2]+0.1767766952966368*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[13] = 0.3061862178478971*(fSkin[31]+fSkin[30])-0.3061862178478971*fSkin[29]+0.3061862178478971*fSkin[28]-0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*fSkin[25]+0.3061862178478971*fSkin[24]-0.3061862178478971*fSkin[23]-0.1767766952966368*fSkin[22]+0.1767766952966368*fSkin[21]-0.1767766952966368*fSkin[20]+0.3061862178478971*fSkin[19]-0.3061862178478971*fSkin[18]+0.3061862178478971*fSkin[17]-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*fSkin[14]-0.1767766952966368*fSkin[13]+0.1767766952966368*fSkin[12]-0.3061862178478971*fSkin[11]+0.3061862178478971*fSkin[10]-0.3061862178478971*fSkin[9]-0.1767766952966368*fSkin[8]+0.1767766952966368*fSkin[7]-0.1767766952966368*fSkin[6]+0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]+0.1767766952966368*fSkin[3]-0.1767766952966368*fSkin[2]+0.1767766952966368*(fSkin[1]+fSkin[0]); 
  } 
  if ((-alpha[12])-alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[14] = (-0.3061862178478971*fEdge[31])+0.3061862178478971*fEdge[30]-0.3061862178478971*(fEdge[29]+fEdge[28])-0.1767766952966368*fEdge[27]-0.3061862178478971*fEdge[26]+0.3061862178478971*(fEdge[25]+fEdge[24])-0.3061862178478971*fEdge[23]+0.1767766952966368*fEdge[22]-0.1767766952966368*(fEdge[21]+fEdge[20])+0.3061862178478971*fEdge[19]-0.3061862178478971*(fEdge[18]+fEdge[17])-0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13])-0.1767766952966368*fEdge[12]+0.3061862178478971*(fEdge[11]+fEdge[10])-0.3061862178478971*fEdge[9]+0.1767766952966368*fEdge[8]-0.1767766952966368*(fEdge[7]+fEdge[6])+0.1767766952966368*fEdge[5]+0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2])-0.1767766952966368*fEdge[1]+0.1767766952966368*fEdge[0]; 
  } else { 

    fUpwindQuad[14] = 0.3061862178478971*fSkin[31]-0.3061862178478971*fSkin[30]+0.3061862178478971*(fSkin[29]+fSkin[28])-0.1767766952966368*fSkin[27]+0.3061862178478971*fSkin[26]-0.3061862178478971*(fSkin[25]+fSkin[24])+0.3061862178478971*fSkin[23]+0.1767766952966368*fSkin[22]-0.1767766952966368*(fSkin[21]+fSkin[20])-0.3061862178478971*fSkin[19]+0.3061862178478971*(fSkin[18]+fSkin[17])-0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13])-0.1767766952966368*fSkin[12]-0.3061862178478971*(fSkin[11]+fSkin[10])+0.3061862178478971*fSkin[9]+0.1767766952966368*fSkin[8]-0.1767766952966368*(fSkin[7]+fSkin[6])+0.1767766952966368*fSkin[5]-0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2])-0.1767766952966368*fSkin[1]+0.1767766952966368*fSkin[0]; 
  } 
  if (alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[15] = 0.3061862178478971*(fEdge[31]+fEdge[30]+fEdge[29]+fEdge[28])+0.1767766952966368*fEdge[27]+0.3061862178478971*(fEdge[26]+fEdge[25]+fEdge[24]+fEdge[23])+0.1767766952966368*(fEdge[22]+fEdge[21]+fEdge[20])+0.3061862178478971*(fEdge[19]+fEdge[18]+fEdge[17])+0.1767766952966368*fEdge[16]+0.3061862178478971*fEdge[15]+0.1767766952966368*(fEdge[14]+fEdge[13]+fEdge[12])+0.3061862178478971*(fEdge[11]+fEdge[10]+fEdge[9])+0.1767766952966368*(fEdge[8]+fEdge[7]+fEdge[6]+fEdge[5])+0.3061862178478971*fEdge[4]+0.1767766952966368*(fEdge[3]+fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[15] = (-0.3061862178478971*(fSkin[31]+fSkin[30]+fSkin[29]+fSkin[28]))+0.1767766952966368*fSkin[27]-0.3061862178478971*(fSkin[26]+fSkin[25]+fSkin[24]+fSkin[23])+0.1767766952966368*(fSkin[22]+fSkin[21]+fSkin[20])-0.3061862178478971*(fSkin[19]+fSkin[18]+fSkin[17])+0.1767766952966368*fSkin[16]-0.3061862178478971*fSkin[15]+0.1767766952966368*(fSkin[14]+fSkin[13]+fSkin[12])-0.3061862178478971*(fSkin[11]+fSkin[10]+fSkin[9])+0.1767766952966368*(fSkin[8]+fSkin[7]+fSkin[6]+fSkin[5])-0.3061862178478971*fSkin[4]+0.1767766952966368*(fSkin[3]+fSkin[2]+fSkin[1]+fSkin[0]); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] = 0.25*(alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.25*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.25*(alpha[12]*fUpwind[15]+alpha[9]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.25*(alpha[11]*fUpwind[15]+alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] = 0.25*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] = 0.25*(alpha[9]*fUpwind[15]+alpha[12]*fUpwind[14]+alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[8]*fUpwind[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] = 0.25*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[12]*fUpwind[13]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] = 0.25*(alpha[7]*fUpwind[15]+alpha[11]*fUpwind[14]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] = 0.25*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[11]*fUpwind[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] = 0.25*(alpha[5]*fUpwind[15]+alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] = 0.25*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[9]*fUpwind[13]+fUpwind[10]*alpha[12]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[12] = 0.25*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+fUpwind[10]*alpha[11]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[13] = 0.25*(alpha[2]*fUpwind[15]+alpha[5]*fUpwind[14]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[14] = 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[5]*fUpwind[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] = 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += 0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -1.224744871391589*Ghat[1]*dv11; 
  out[10] += -1.224744871391589*Ghat[2]*dv11; 
  out[11] += -1.224744871391589*Ghat[3]*dv11; 
  out[12] += 0.7071067811865475*Ghat[8]*dv11; 
  out[13] += 0.7071067811865475*Ghat[9]*dv11; 
  out[14] += 0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -1.224744871391589*Ghat[4]*dv11; 
  out[16] += 0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -1.224744871391589*Ghat[7]*dv11; 
  out[20] += 0.7071067811865475*Ghat[12]*dv11; 
  out[21] += 0.7071067811865475*Ghat[13]*dv11; 
  out[22] += 0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.224744871391589*Ghat[9]*dv11; 
  out[25] += -1.224744871391589*Ghat[10]*dv11; 
  out[26] += -1.224744871391589*Ghat[11]*dv11; 
  out[27] += 0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -1.224744871391589*Ghat[12]*dv11; 
  out[29] += -1.224744871391589*Ghat[13]*dv11; 
  out[30] += -1.224744871391589*Ghat[14]*dv11; 
  out[31] += -1.224744871391589*Ghat[15]*dv11; 

  } 
} 
