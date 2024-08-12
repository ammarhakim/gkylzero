#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_lorentz_vmap_2v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT gamma, double* GKYL_RESTRICT gamma_inv) 
{ 
  // w:   Cell-center coordinates.
  // dxv: Cell spacing.
  // vmap:  Momentum-space nonuniform mapping.
  // gamma:  Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv: Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
 
  const double *p0 = &vmap[0]; 
  const double *p1 = &vmap[4]; 
  double gamma_nodal[8] = {0.0};
  double gamma_inv_nodal[8] = {0.0};

  gamma_nodal[0] = sqrt(1.0 + pow((-1.870828693386971*p0[3])+1.58113883008419*p0[2]-1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0) + pow((-1.870828693386971*p1[3])+1.58113883008419*p1[2]-1.224744871391589*p1[1]+0.7071067811865475*p1[0], 2.0));
  gamma_inv_nodal[0] = 1.0/gamma_nodal[0];
  gamma_nodal[1] = sqrt(1.0 + pow(0.7071067811865475*p0[0]-0.7905694150420947*p0[2], 2.0) + pow((-1.870828693386971*p1[3])+1.58113883008419*p1[2]-1.224744871391589*p1[1]+0.7071067811865475*p1[0], 2.0));
  gamma_inv_nodal[1] = 1.0/gamma_nodal[1];
  gamma_nodal[2] = sqrt(1.0 + pow(1.870828693386971*p0[3]+1.58113883008419*p0[2]+1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0) + pow((-1.870828693386971*p1[3])+1.58113883008419*p1[2]-1.224744871391589*p1[1]+0.7071067811865475*p1[0], 2.0));
  gamma_inv_nodal[2] = 1.0/gamma_nodal[2];
  gamma_nodal[3] = sqrt(1.0 + pow((-1.870828693386971*p0[3])+1.58113883008419*p0[2]-1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0) + pow(0.7071067811865475*p1[0]-0.7905694150420947*p1[2], 2.0));
  gamma_inv_nodal[3] = 1.0/gamma_nodal[3];
  gamma_nodal[4] = sqrt(1.0 + pow(1.870828693386971*p0[3]+1.58113883008419*p0[2]+1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0) + pow(0.7071067811865475*p1[0]-0.7905694150420947*p1[2], 2.0));
  gamma_inv_nodal[4] = 1.0/gamma_nodal[4];
  gamma_nodal[5] = sqrt(1.0 + pow((-1.870828693386971*p0[3])+1.58113883008419*p0[2]-1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0) + pow(1.870828693386971*p1[3]+1.58113883008419*p1[2]+1.224744871391589*p1[1]+0.7071067811865475*p1[0], 2.0));
  gamma_inv_nodal[5] = 1.0/gamma_nodal[5];
  gamma_nodal[6] = sqrt(1.0 + pow(0.7071067811865475*p0[0]-0.7905694150420947*p0[2], 2.0) + pow(1.870828693386971*p1[3]+1.58113883008419*p1[2]+1.224744871391589*p1[1]+0.7071067811865475*p1[0], 2.0));
  gamma_inv_nodal[6] = 1.0/gamma_nodal[6];
  gamma_nodal[7] = sqrt(1.0 + pow(1.870828693386971*p0[3]+1.58113883008419*p0[2]+1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0) + pow(1.870828693386971*p1[3]+1.58113883008419*p1[2]+1.224744871391589*p1[1]+0.7071067811865475*p1[0], 2.0));
  gamma_inv_nodal[7] = 1.0/gamma_nodal[7];

  gamma[0] = (-0.1666666666666667*gamma_nodal[7])+0.6666666666666666*gamma_nodal[6]-0.1666666666666667*gamma_nodal[5]+0.6666666666666666*gamma_nodal[4]+0.6666666666666666*gamma_nodal[3]-0.1666666666666667*gamma_nodal[2]+0.6666666666666666*gamma_nodal[1]-0.1666666666666667*gamma_nodal[0]; 
  gamma[1] = 0.09622504486493764*gamma_nodal[7]-0.09622504486493764*gamma_nodal[5]+0.3849001794597506*gamma_nodal[4]-0.3849001794597506*gamma_nodal[3]+0.09622504486493764*gamma_nodal[2]-0.09622504486493764*gamma_nodal[0]; 
  gamma[2] = 0.09622504486493764*gamma_nodal[7]+0.3849001794597506*gamma_nodal[6]+0.09622504486493764*gamma_nodal[5]-0.09622504486493764*gamma_nodal[2]-0.3849001794597506*gamma_nodal[1]-0.09622504486493764*gamma_nodal[0]; 
  gamma[3] = 0.1666666666666667*gamma_nodal[7]-0.1666666666666667*gamma_nodal[5]-0.1666666666666667*gamma_nodal[2]+0.1666666666666667*gamma_nodal[0]; 
  gamma[4] = 0.149071198499986*gamma_nodal[7]-0.2981423969999719*gamma_nodal[6]+0.149071198499986*gamma_nodal[5]+0.149071198499986*gamma_nodal[2]-0.2981423969999719*gamma_nodal[1]+0.149071198499986*gamma_nodal[0]; 
  gamma[5] = 0.149071198499986*gamma_nodal[7]+0.149071198499986*gamma_nodal[5]-0.2981423969999719*gamma_nodal[4]-0.2981423969999719*gamma_nodal[3]+0.149071198499986*gamma_nodal[2]+0.149071198499986*gamma_nodal[0]; 
  gamma[6] = 0.08606629658238703*gamma_nodal[7]-0.1721325931647741*gamma_nodal[6]+0.08606629658238703*gamma_nodal[5]-0.08606629658238703*gamma_nodal[2]+0.1721325931647741*gamma_nodal[1]-0.08606629658238703*gamma_nodal[0]; 
  gamma[7] = 0.08606629658238703*gamma_nodal[7]-0.08606629658238703*gamma_nodal[5]-0.1721325931647741*gamma_nodal[4]+0.1721325931647741*gamma_nodal[3]+0.08606629658238703*gamma_nodal[2]-0.08606629658238703*gamma_nodal[0]; 

  gamma_inv[0] = (-0.1666666666666667*gamma_inv_nodal[7])+0.6666666666666666*gamma_inv_nodal[6]-0.1666666666666667*gamma_inv_nodal[5]+0.6666666666666666*gamma_inv_nodal[4]+0.6666666666666666*gamma_inv_nodal[3]-0.1666666666666667*gamma_inv_nodal[2]+0.6666666666666666*gamma_inv_nodal[1]-0.1666666666666667*gamma_inv_nodal[0]; 
  gamma_inv[1] = 0.09622504486493764*gamma_inv_nodal[7]-0.09622504486493764*gamma_inv_nodal[5]+0.3849001794597506*gamma_inv_nodal[4]-0.3849001794597506*gamma_inv_nodal[3]+0.09622504486493764*gamma_inv_nodal[2]-0.09622504486493764*gamma_inv_nodal[0]; 
  gamma_inv[2] = 0.09622504486493764*gamma_inv_nodal[7]+0.3849001794597506*gamma_inv_nodal[6]+0.09622504486493764*gamma_inv_nodal[5]-0.09622504486493764*gamma_inv_nodal[2]-0.3849001794597506*gamma_inv_nodal[1]-0.09622504486493764*gamma_inv_nodal[0]; 
  gamma_inv[3] = 0.1666666666666667*gamma_inv_nodal[7]-0.1666666666666667*gamma_inv_nodal[5]-0.1666666666666667*gamma_inv_nodal[2]+0.1666666666666667*gamma_inv_nodal[0]; 
  gamma_inv[4] = 0.149071198499986*gamma_inv_nodal[7]-0.2981423969999719*gamma_inv_nodal[6]+0.149071198499986*gamma_inv_nodal[5]+0.149071198499986*gamma_inv_nodal[2]-0.2981423969999719*gamma_inv_nodal[1]+0.149071198499986*gamma_inv_nodal[0]; 
  gamma_inv[5] = 0.149071198499986*gamma_inv_nodal[7]+0.149071198499986*gamma_inv_nodal[5]-0.2981423969999719*gamma_inv_nodal[4]-0.2981423969999719*gamma_inv_nodal[3]+0.149071198499986*gamma_inv_nodal[2]+0.149071198499986*gamma_inv_nodal[0]; 
  gamma_inv[6] = 0.08606629658238703*gamma_inv_nodal[7]-0.1721325931647741*gamma_inv_nodal[6]+0.08606629658238703*gamma_inv_nodal[5]-0.08606629658238703*gamma_inv_nodal[2]+0.1721325931647741*gamma_inv_nodal[1]-0.08606629658238703*gamma_inv_nodal[0]; 
  gamma_inv[7] = 0.08606629658238703*gamma_inv_nodal[7]-0.08606629658238703*gamma_inv_nodal[5]-0.1721325931647741*gamma_inv_nodal[4]+0.1721325931647741*gamma_inv_nodal[3]+0.08606629658238703*gamma_inv_nodal[2]-0.08606629658238703*gamma_inv_nodal[0]; 
} 
