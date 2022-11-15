#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out) 
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

  const double *rhoux = &statevec[0]; 
  const double *rhouy = &statevec[3]; 
  const double *rhouz = &statevec[6]; 
  const double *E_perp = &statevec[9]; 

  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[3]; 
  const double *div_p_z = &div_p[6]; 

  const double *u_perp_x = &u_perp_i[0]; 
  const double *u_perp_y = &u_perp_i[3]; 
  const double *u_perp_z = &u_perp_i[6]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[3]; 
  double *outrhouz = &out[6]; 
  double *outE_perp = &out[9]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*5.0*dx10*(fabs(0.7071067811865475*ux[0]-0.7905694150420947*ux[2])); 

  outrhoux[0] += -1.0*div_p_x[0]; 
  outrhoux[1] += 1.224744871391589*rhoux[2]*ux[2]*dx10+1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10-1.0*div_p_x[1]; 
  outrhoux[2] += 2.449489742783178*rhoux[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhoux[2]*dx10+2.738612787525831*rhoux[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhoux[1]*dx10-1.0*div_p_x[2]; 

  outrhouy[0] += -1.0*div_p_y[0]; 
  outrhouy[1] += 1.224744871391589*rhouy[2]*ux[2]*dx10+1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10-1.0*div_p_y[1]; 
  outrhouy[2] += 2.449489742783178*rhouy[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhouy[2]*dx10+2.738612787525831*rhouy[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhouy[1]*dx10-1.0*div_p_y[2]; 

  outrhouz[0] += -1.0*div_p_z[0]; 
  outrhouz[1] += 1.224744871391589*rhouz[2]*ux[2]*dx10+1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10-1.0*div_p_z[1]; 
  outrhouz[2] += 2.449489742783178*rhouz[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhouz[2]*dx10+2.738612787525831*rhouz[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhouz[1]*dx10-1.0*div_p_z[2]; 

  outE_perp[1] += 1.224744871391589*E_perp[2]*ux[2]*dx10+1.224744871391589*p_perp[2]*u_perp_x[2]*dx10+1.224744871391589*E_perp[1]*ux[1]*dx10+1.224744871391589*p_perp[1]*u_perp_x[1]*dx10+1.224744871391589*E_perp[0]*ux[0]*dx10+1.224744871391589*p_perp[0]*u_perp_x[0]*dx10; 
  outE_perp[2] += 2.449489742783178*E_perp[1]*ux[2]*dx10+2.449489742783178*p_perp[1]*u_perp_x[2]*dx10+2.449489742783178*u_perp_x[1]*p_perp[2]*dx10+2.449489742783178*ux[1]*E_perp[2]*dx10+2.738612787525831*E_perp[0]*ux[1]*dx10+2.738612787525831*p_perp[0]*u_perp_x[1]*dx10+2.738612787525831*u_perp_x[0]*p_perp[1]*dx10+2.738612787525831*ux[0]*E_perp[1]*dx10; 

  return cflFreq_mid; 
} 
