#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_1x_tensor_p1(const double *w, const double *dxv, 
  const double *vlasov_pkpm_moms, const double *pkpm_u, const double *p_ij, 
  const double *euler_pkpm, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing.
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // pkpm_u:     Input flow velocity [ux, uy, uz].
  // p_ij:       Input pressure tensor.
  // euler_pkpm: [rho ux, rho uy, rho uz], Fluid input state vector.
  // out:        Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[2]; 
  const double *rhouz = &euler_pkpm[4]; 

  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[2]; 
  const double *uz = &pkpm_u[4]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[3]; 
  const double *Pxz = &p_ij[6]; 
  const double *Pyy = &p_ij[9]; 
  const double *Pyz = &p_ij[12]; 
  const double *Pzz = &p_ij[15]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[2]; 
  double *outrhouz = &out[4]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.7071067811865475*ux[0])); 

  double rhoux_2p[3] = {0.0}; 
  double rhouy_2p[3] = {0.0}; 
  double rhouz_2p[3] = {0.0}; 
  rhoux_2p[0] = 0.7071067811865475*rho[1]*ux[1]+0.7071067811865475*rho[0]*ux[0]; 
  rhoux_2p[1] = 0.6324555320336759*ux[1]*rho[2]+0.7071067811865475*rho[0]*ux[1]+0.7071067811865475*ux[0]*rho[1]; 
  rhoux_2p[2] = 0.7071067811865475*ux[0]*rho[2]+0.6324555320336759*rho[1]*ux[1]; 

  rhouy_2p[0] = 0.7071067811865475*rho[1]*uy[1]+0.7071067811865475*rho[0]*uy[0]; 
  rhouy_2p[1] = 0.6324555320336759*uy[1]*rho[2]+0.7071067811865475*rho[0]*uy[1]+0.7071067811865475*uy[0]*rho[1]; 
  rhouy_2p[2] = 0.7071067811865475*uy[0]*rho[2]+0.6324555320336759*rho[1]*uy[1]; 

  rhouz_2p[0] = 0.7071067811865475*rho[1]*uz[1]+0.7071067811865475*rho[0]*uz[0]; 
  rhouz_2p[1] = 0.6324555320336759*uz[1]*rho[2]+0.7071067811865475*rho[0]*uz[1]+0.7071067811865475*uz[0]*rho[1]; 
  rhouz_2p[2] = 0.7071067811865475*uz[0]*rho[2]+0.6324555320336759*rho[1]*uz[1]; 

  outrhoux[1] += 1.224744871391589*rhoux_2p[1]*ux[1]*dx10+1.224744871391589*rhoux_2p[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy_2p[1]*ux[1]*dx10+1.224744871391589*rhouy_2p[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz_2p[1]*ux[1]*dx10+1.224744871391589*rhouz_2p[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 

  return cflFreq_mid; 
} 
