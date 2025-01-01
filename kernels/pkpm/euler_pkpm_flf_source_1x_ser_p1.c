#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_flf_source_1x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *div_b, double* out) 
{ 
  // qmem:             q/m*EM fields.
  // div_b:            Input geometry factor div(b) = -1/B*dB/dz.
  // vlasov_pkpm_moms: [Jrho, Jp_parallel, Jp_perp], Moments computed from kinetic equation in pkpm model.
  //                   Includes Jacobian factor 1/B.
  // out:              Output increment
  const double *Jrho = &vlasov_pkpm_moms[0]; 
  const double *Jp_perp = &vlasov_pkpm_moms[4]; 
  const double *Epar = &qmem[0]; 

  double *outrhoux = &out[0]; 

  outrhoux[0] += 0.7071067811865475*Jp_perp[1]*div_b[1]+0.7071067811865475*Epar[1]*Jrho[1]+0.7071067811865475*Jp_perp[0]*div_b[0]+0.7071067811865475*Epar[0]*Jrho[0]; 
  outrhoux[1] += 0.7071067811865475*Jp_perp[0]*div_b[1]+0.7071067811865475*Epar[0]*Jrho[1]+0.7071067811865475*div_b[0]*Jp_perp[1]+0.7071067811865475*Jrho[0]*Epar[1]; 

} 
