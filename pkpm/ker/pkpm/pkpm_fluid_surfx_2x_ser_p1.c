#include <gkyl_pkpm_fluid_kernels.h> 
GKYL_CU_DH double pkpm_fluid_surfx_2x_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  // flux_l/flux_r: Stored momentum fluxes [F_x, F_y, F_z] on the lower/upper face of the cell. 
  // out:           Incremented momentum increment [rho ux, rho uy, rho uz]. 
  const double dx1 = 2.0/dxv[0]; 
  const double *gl0 = &flux_l[0]; 
  const double *gr0 = &flux_r[0]; 
  double *out0 = &out[0]; 
  const double *gl1 = &flux_l[2]; 
  const double *gr1 = &flux_r[2]; 
  double *out1 = &out[4]; 
  const double *gl2 = &flux_l[4]; 
  const double *gr2 = &flux_r[4]; 
  double *out2 = &out[8]; 
  out0[0] += (0.7071067811865475*gl0[0]-0.7071067811865475*gr0[0])*dx1; 
  out0[1] += (-(1.224744871391589*gr0[0])-1.224744871391589*gl0[0])*dx1; 
  out0[2] += (0.7071067811865475*gl0[1]-0.7071067811865475*gr0[1])*dx1; 
  out0[3] += (-(1.224744871391589*gr0[1])-1.224744871391589*gl0[1])*dx1; 
  out1[0] += (0.7071067811865475*gl1[0]-0.7071067811865475*gr1[0])*dx1; 
  out1[1] += (-(1.224744871391589*gr1[0])-1.224744871391589*gl1[0])*dx1; 
  out1[2] += (0.7071067811865475*gl1[1]-0.7071067811865475*gr1[1])*dx1; 
  out1[3] += (-(1.224744871391589*gr1[1])-1.224744871391589*gl1[1])*dx1; 
  out2[0] += (0.7071067811865475*gl2[0]-0.7071067811865475*gr2[0])*dx1; 
  out2[1] += (-(1.224744871391589*gr2[0])-1.224744871391589*gl2[0])*dx1; 
  out2[2] += (0.7071067811865475*gl2[1]-0.7071067811865475*gr2[1])*dx1; 
  out2[3] += (-(1.224744871391589*gr2[1])-1.224744871391589*gl2[1])*dx1; 
  return 0.0; 
} 
