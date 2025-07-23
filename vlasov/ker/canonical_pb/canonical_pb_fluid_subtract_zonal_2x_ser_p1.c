#include <gkyl_canonical_pb_kernels.h> 
void canonical_pb_fluid_subtract_zonal_2x_ser_p1(const double *phi_zonal, const double *n_zonal, double* GKYL_RESTRICT adiabatic_coupling_phi_n) 
{ 
  // phi_zonal: 1/Ly int phi dy.
  // n_zonal: 1/Ly int n dy.
  // adiabatic_coupling_phi_n: (phi, n) array for adiabatic coupling, with zonal component subtracted out.
  double *out_phi = &adiabatic_coupling_phi_n[0]; 
  double *out_n = &adiabatic_coupling_phi_n[4]; 
  out_phi[0] -= 1.4142135623730951*phi_zonal[0]; 
  out_n[0] -= 1.4142135623730951*n_zonal[0]; 
  out_phi[1] -= 1.4142135623730951*phi_zonal[1]; 
  out_n[1] -= 1.4142135623730951*n_zonal[1]; 
} 
