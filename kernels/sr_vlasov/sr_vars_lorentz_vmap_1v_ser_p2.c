#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_lorentz_vmap_1v_ser_p2(const double *w, const double *dxv, const double *vmap, double* GKYL_RESTRICT gamma, double* GKYL_RESTRICT gamma_inv) 
{ 
  // w:   Cell-center coordinates.
  // dxv: Cell spacing.
  // vmap:  Momentum-space nonuniform mapping.
  // gamma:  Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv: Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
 
  const double *p0 = &vmap[0]; 
  double gamma_nodal[3] = {0.0};
  double gamma_inv_nodal[3] = {0.0};

  gamma_nodal[0] = sqrt(1.0 + pow((-1.870828693386971*p0[3])+1.58113883008419*p0[2]-1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0));
  gamma_inv_nodal[0] = 1.0/gamma_nodal[0];
  gamma_nodal[1] = sqrt(1.0 + pow(0.7071067811865475*p0[0]-0.7905694150420947*p0[2], 2.0));
  gamma_inv_nodal[1] = 1.0/gamma_nodal[1];
  gamma_nodal[2] = sqrt(1.0 + pow(1.870828693386971*p0[3]+1.58113883008419*p0[2]+1.224744871391589*p0[1]+0.7071067811865475*p0[0], 2.0));
  gamma_inv_nodal[2] = 1.0/gamma_nodal[2];

  gamma[0] = 0.2357022603955158*gamma_nodal[2]+0.9428090415820636*gamma_nodal[1]+0.2357022603955158*gamma_nodal[0]; 
  gamma[1] = 0.408248290463863*gamma_nodal[2]-0.408248290463863*gamma_nodal[0]; 
  gamma[2] = 0.210818510677892*gamma_nodal[2]-0.421637021355784*gamma_nodal[1]+0.210818510677892*gamma_nodal[0]; 

  gamma_inv[0] = 0.2357022603955158*gamma_inv_nodal[2]+0.9428090415820636*gamma_inv_nodal[1]+0.2357022603955158*gamma_inv_nodal[0]; 
  gamma_inv[1] = 0.408248290463863*gamma_inv_nodal[2]-0.408248290463863*gamma_inv_nodal[0]; 
  gamma_inv[2] = 0.210818510677892*gamma_inv_nodal[2]-0.421637021355784*gamma_inv_nodal[1]+0.210818510677892*gamma_inv_nodal[0]; 
} 
