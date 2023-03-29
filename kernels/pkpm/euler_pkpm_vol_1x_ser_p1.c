#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:       bulk flow velocity (ux, uy, uz).
  // div_p:     divergence of the pressure tensor.
  // statevec: [rho ux, rho uy, rho uz], Fluid input state vector.
  // out: Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[2]; 
  const double *rhouz = &statevec[4]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[2]; 
  const double *uz = &u_i[4]; 

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[2]; 
  const double *div_p_z = &div_p[4]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[2]; 
  double *outrhouz = &out[4]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.7071067811865475*ux[0])); 

  outrhoux[0] += -1.0*div_p_x[0]; 
  outrhoux[1] += 1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10-1.0*div_p_x[1]; 

  outrhouy[0] += -1.0*div_p_y[0]; 
  outrhouy[1] += 1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10-1.0*div_p_y[1]; 

  outrhouz[0] += -1.0*div_p_z[0]; 
  outrhouz[1] += 1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10-1.0*div_p_z[1]; 

  return cflFreq_mid; 
} 
