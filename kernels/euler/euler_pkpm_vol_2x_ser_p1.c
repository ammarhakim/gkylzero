#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      Input bulk flow velocity [ux, uy, uz].
  // div_p:    Input divergence of the pressure tensor.
  // u_perp_i: Input perpendicular bulk velocity [u_perp_x, u_perp_y, u_perp_z].
  // p_perp:   Input perpendicular pressure.
  // statevec: [rho ux, rho uy, rho uz, p_perp], Fluid input state vector.
  // out: Incremented output.

  double dx10 = 2./dxv[0]; 
  double dx11 = 2./dxv[1]; 

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[4]; 
  const double *rhouz = &statevec[8]; 
  const double *E_perp = &statevec[12]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[4]; 
  const double *div_p_z = &div_p[8]; 

  const double *u_perp_x = &u_perp_i[0]; 
  const double *u_perp_y = &u_perp_i[4]; 
  const double *u_perp_z = &u_perp_i[8]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[8]; 
  double *outE_perp = &out[12]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.5*ux[0])); 
  cflFreq_mid += 0.5*3.0*dx11*(fabs(0.5*uy[0])); 

  outrhoux[0] += -1.0*div_p_x[0]; 
  outrhoux[1] += 0.8660254037844386*rhoux[3]*ux[3]*dx10+0.8660254037844386*rhoux[2]*ux[2]*dx10+0.8660254037844386*rhoux[1]*ux[1]*dx10+0.8660254037844386*rhoux[0]*ux[0]*dx10-1.0*div_p_x[1]; 
  outrhoux[2] += 0.8660254037844386*rhoux[3]*uy[3]*dx11+0.8660254037844386*rhoux[2]*uy[2]*dx11+0.8660254037844386*rhoux[1]*uy[1]*dx11+0.8660254037844386*rhoux[0]*uy[0]*dx11-1.0*div_p_x[2]; 
  outrhoux[3] += 0.8660254037844386*rhoux[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhoux[3]*dx11+0.8660254037844386*rhoux[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhoux[1]*dx11+0.8660254037844386*rhoux[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhoux[3]*dx10+0.8660254037844386*rhoux[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhoux[2]*dx10-1.0*div_p_x[3]; 

  outrhouy[0] += -1.0*div_p_y[0]; 
  outrhouy[1] += 0.8660254037844386*rhouy[3]*ux[3]*dx10+0.8660254037844386*rhouy[2]*ux[2]*dx10+0.8660254037844386*rhouy[1]*ux[1]*dx10+0.8660254037844386*rhouy[0]*ux[0]*dx10-1.0*div_p_y[1]; 
  outrhouy[2] += 0.8660254037844386*rhouy[3]*uy[3]*dx11+0.8660254037844386*rhouy[2]*uy[2]*dx11+0.8660254037844386*rhouy[1]*uy[1]*dx11+0.8660254037844386*rhouy[0]*uy[0]*dx11-1.0*div_p_y[2]; 
  outrhouy[3] += 0.8660254037844386*rhouy[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouy[3]*dx11+0.8660254037844386*rhouy[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouy[1]*dx11+0.8660254037844386*rhouy[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouy[3]*dx10+0.8660254037844386*rhouy[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouy[2]*dx10-1.0*div_p_y[3]; 

  outrhouz[0] += -1.0*div_p_z[0]; 
  outrhouz[1] += 0.8660254037844386*rhouz[3]*ux[3]*dx10+0.8660254037844386*rhouz[2]*ux[2]*dx10+0.8660254037844386*rhouz[1]*ux[1]*dx10+0.8660254037844386*rhouz[0]*ux[0]*dx10-1.0*div_p_z[1]; 
  outrhouz[2] += 0.8660254037844386*rhouz[3]*uy[3]*dx11+0.8660254037844386*rhouz[2]*uy[2]*dx11+0.8660254037844386*rhouz[1]*uy[1]*dx11+0.8660254037844386*rhouz[0]*uy[0]*dx11-1.0*div_p_z[2]; 
  outrhouz[3] += 0.8660254037844386*rhouz[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouz[3]*dx11+0.8660254037844386*rhouz[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouz[1]*dx11+0.8660254037844386*rhouz[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouz[3]*dx10+0.8660254037844386*rhouz[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouz[2]*dx10-1.0*div_p_z[3]; 

  outE_perp[1] += 0.8660254037844386*E_perp[3]*ux[3]*dx10+0.8660254037844386*p_perp[3]*u_perp_x[3]*dx10+0.8660254037844386*E_perp[2]*ux[2]*dx10+0.8660254037844386*p_perp[2]*u_perp_x[2]*dx10+0.8660254037844386*E_perp[1]*ux[1]*dx10+0.8660254037844386*p_perp[1]*u_perp_x[1]*dx10+0.8660254037844386*E_perp[0]*ux[0]*dx10+0.8660254037844386*p_perp[0]*u_perp_x[0]*dx10; 
  outE_perp[2] += 0.8660254037844386*E_perp[3]*uy[3]*dx11+0.8660254037844386*p_perp[3]*u_perp_y[3]*dx11+0.8660254037844386*E_perp[2]*uy[2]*dx11+0.8660254037844386*p_perp[2]*u_perp_y[2]*dx11+0.8660254037844386*E_perp[1]*uy[1]*dx11+0.8660254037844386*p_perp[1]*u_perp_y[1]*dx11+0.8660254037844386*E_perp[0]*uy[0]*dx11+0.8660254037844386*p_perp[0]*u_perp_y[0]*dx11; 
  outE_perp[3] += 0.8660254037844386*E_perp[2]*uy[3]*dx11+0.8660254037844386*p_perp[2]*u_perp_y[3]*dx11+0.8660254037844386*u_perp_y[2]*p_perp[3]*dx11+0.8660254037844386*uy[2]*E_perp[3]*dx11+0.8660254037844386*E_perp[0]*uy[1]*dx11+0.8660254037844386*p_perp[0]*u_perp_y[1]*dx11+0.8660254037844386*u_perp_y[0]*p_perp[1]*dx11+0.8660254037844386*uy[0]*E_perp[1]*dx11+0.8660254037844386*E_perp[1]*ux[3]*dx10+0.8660254037844386*p_perp[1]*u_perp_x[3]*dx10+0.8660254037844386*u_perp_x[1]*p_perp[3]*dx10+0.8660254037844386*ux[1]*E_perp[3]*dx10+0.8660254037844386*E_perp[0]*ux[2]*dx10+0.8660254037844386*p_perp[0]*u_perp_x[2]*dx10+0.8660254037844386*u_perp_x[0]*p_perp[2]*dx10+0.8660254037844386*ux[0]*E_perp[2]*dx10; 

  return cflFreq_mid; 
} 
