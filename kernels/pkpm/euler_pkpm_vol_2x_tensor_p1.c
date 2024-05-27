#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_2x_tensor_p1(const double *w, const double *dxv, 
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
  double dx11 = 2./dxv[1]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 

  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *ux = &pkpm_u[0]; 
  const double *uy = &pkpm_u[4]; 
  const double *uz = &pkpm_u[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[9]; 
  const double *Pxz = &p_ij[18]; 
  const double *Pyy = &p_ij[27]; 
  const double *Pyz = &p_ij[36]; 
  const double *Pzz = &p_ij[45]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[8]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.5*ux[0])); 
  cflFreq_mid += 0.5*3.0*dx11*(fabs(0.5*uy[0])); 

  double rhoux_2p[9] = {0.0}; 
  double rhouy_2p[9] = {0.0}; 
  double rhouz_2p[9] = {0.0}; 
  rhoux_2p[0] = 0.5*rho[3]*ux[3]+0.5*rho[2]*ux[2]+0.5*rho[1]*ux[1]+0.5*rho[0]*ux[0]; 
  rhoux_2p[1] = 0.447213595499958*ux[3]*rho[6]+0.4472135954999579*ux[1]*rho[4]+0.5*rho[2]*ux[3]+0.5*ux[2]*rho[3]+0.5*rho[0]*ux[1]+0.5*ux[0]*rho[1]; 
  rhoux_2p[2] = 0.447213595499958*ux[3]*rho[7]+0.4472135954999579*ux[2]*rho[5]+0.5*rho[1]*ux[3]+0.5*ux[1]*rho[3]+0.5*rho[0]*ux[2]+0.5*ux[0]*rho[2]; 
  rhoux_2p[3] = 0.4*ux[3]*rho[8]+0.447213595499958*ux[2]*rho[7]+0.447213595499958*ux[1]*rho[6]+0.4472135954999579*ux[3]*rho[5]+0.4472135954999579*ux[3]*rho[4]+0.5*rho[0]*ux[3]+0.5*ux[0]*rho[3]+0.5*rho[1]*ux[2]+0.5*ux[1]*rho[2]; 
  rhoux_2p[4] = 0.5000000000000001*ux[2]*rho[6]+0.5*ux[0]*rho[4]+0.4472135954999579*rho[3]*ux[3]+0.4472135954999579*rho[1]*ux[1]; 
  rhoux_2p[5] = 0.5000000000000001*ux[1]*rho[7]+0.5*ux[0]*rho[5]+0.4472135954999579*rho[3]*ux[3]+0.4472135954999579*rho[2]*ux[2]; 
  rhoux_2p[6] = 0.447213595499958*ux[2]*rho[8]+0.4*ux[3]*rho[7]+0.5*ux[0]*rho[6]+0.5000000000000001*ux[2]*rho[4]+0.447213595499958*rho[1]*ux[3]+0.447213595499958*ux[1]*rho[3]; 
  rhoux_2p[7] = 0.447213595499958*ux[1]*rho[8]+0.5*ux[0]*rho[7]+0.4*ux[3]*rho[6]+0.5000000000000001*ux[1]*rho[5]+0.447213595499958*rho[2]*ux[3]+0.447213595499958*ux[2]*rho[3]; 
  rhoux_2p[8] = 0.5*ux[0]*rho[8]+0.447213595499958*ux[1]*rho[7]+0.447213595499958*ux[2]*rho[6]+0.4*rho[3]*ux[3]; 

  rhouy_2p[0] = 0.5*rho[3]*uy[3]+0.5*rho[2]*uy[2]+0.5*rho[1]*uy[1]+0.5*rho[0]*uy[0]; 
  rhouy_2p[1] = 0.447213595499958*uy[3]*rho[6]+0.4472135954999579*uy[1]*rho[4]+0.5*rho[2]*uy[3]+0.5*uy[2]*rho[3]+0.5*rho[0]*uy[1]+0.5*uy[0]*rho[1]; 
  rhouy_2p[2] = 0.447213595499958*uy[3]*rho[7]+0.4472135954999579*uy[2]*rho[5]+0.5*rho[1]*uy[3]+0.5*uy[1]*rho[3]+0.5*rho[0]*uy[2]+0.5*uy[0]*rho[2]; 
  rhouy_2p[3] = 0.4*uy[3]*rho[8]+0.447213595499958*uy[2]*rho[7]+0.447213595499958*uy[1]*rho[6]+0.4472135954999579*uy[3]*rho[5]+0.4472135954999579*uy[3]*rho[4]+0.5*rho[0]*uy[3]+0.5*uy[0]*rho[3]+0.5*rho[1]*uy[2]+0.5*uy[1]*rho[2]; 
  rhouy_2p[4] = 0.5000000000000001*uy[2]*rho[6]+0.5*uy[0]*rho[4]+0.4472135954999579*rho[3]*uy[3]+0.4472135954999579*rho[1]*uy[1]; 
  rhouy_2p[5] = 0.5000000000000001*uy[1]*rho[7]+0.5*uy[0]*rho[5]+0.4472135954999579*rho[3]*uy[3]+0.4472135954999579*rho[2]*uy[2]; 
  rhouy_2p[6] = 0.447213595499958*uy[2]*rho[8]+0.4*uy[3]*rho[7]+0.5*uy[0]*rho[6]+0.5000000000000001*uy[2]*rho[4]+0.447213595499958*rho[1]*uy[3]+0.447213595499958*uy[1]*rho[3]; 
  rhouy_2p[7] = 0.447213595499958*uy[1]*rho[8]+0.5*uy[0]*rho[7]+0.4*uy[3]*rho[6]+0.5000000000000001*uy[1]*rho[5]+0.447213595499958*rho[2]*uy[3]+0.447213595499958*uy[2]*rho[3]; 
  rhouy_2p[8] = 0.5*uy[0]*rho[8]+0.447213595499958*uy[1]*rho[7]+0.447213595499958*uy[2]*rho[6]+0.4*rho[3]*uy[3]; 

  rhouz_2p[0] = 0.5*rho[3]*uz[3]+0.5*rho[2]*uz[2]+0.5*rho[1]*uz[1]+0.5*rho[0]*uz[0]; 
  rhouz_2p[1] = 0.447213595499958*uz[3]*rho[6]+0.4472135954999579*uz[1]*rho[4]+0.5*rho[2]*uz[3]+0.5*uz[2]*rho[3]+0.5*rho[0]*uz[1]+0.5*uz[0]*rho[1]; 
  rhouz_2p[2] = 0.447213595499958*uz[3]*rho[7]+0.4472135954999579*uz[2]*rho[5]+0.5*rho[1]*uz[3]+0.5*uz[1]*rho[3]+0.5*rho[0]*uz[2]+0.5*uz[0]*rho[2]; 
  rhouz_2p[3] = 0.4*uz[3]*rho[8]+0.447213595499958*uz[2]*rho[7]+0.447213595499958*uz[1]*rho[6]+0.4472135954999579*uz[3]*rho[5]+0.4472135954999579*uz[3]*rho[4]+0.5*rho[0]*uz[3]+0.5*uz[0]*rho[3]+0.5*rho[1]*uz[2]+0.5*uz[1]*rho[2]; 
  rhouz_2p[4] = 0.5000000000000001*uz[2]*rho[6]+0.5*uz[0]*rho[4]+0.4472135954999579*rho[3]*uz[3]+0.4472135954999579*rho[1]*uz[1]; 
  rhouz_2p[5] = 0.5000000000000001*uz[1]*rho[7]+0.5*uz[0]*rho[5]+0.4472135954999579*rho[3]*uz[3]+0.4472135954999579*rho[2]*uz[2]; 
  rhouz_2p[6] = 0.447213595499958*uz[2]*rho[8]+0.4*uz[3]*rho[7]+0.5*uz[0]*rho[6]+0.5000000000000001*uz[2]*rho[4]+0.447213595499958*rho[1]*uz[3]+0.447213595499958*uz[1]*rho[3]; 
  rhouz_2p[7] = 0.447213595499958*uz[1]*rho[8]+0.5*uz[0]*rho[7]+0.4*uz[3]*rho[6]+0.5000000000000001*uz[1]*rho[5]+0.447213595499958*rho[2]*uz[3]+0.447213595499958*uz[2]*rho[3]; 
  rhouz_2p[8] = 0.5*uz[0]*rho[8]+0.447213595499958*uz[1]*rho[7]+0.447213595499958*uz[2]*rho[6]+0.4*rho[3]*uz[3]; 

  outrhoux[1] += 0.8660254037844386*rhoux_2p[3]*ux[3]*dx10+0.8660254037844386*rhoux_2p[2]*ux[2]*dx10+0.8660254037844386*rhoux_2p[1]*ux[1]*dx10+0.8660254037844386*rhoux_2p[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 
  outrhoux[2] += 0.8660254037844386*rhoux_2p[3]*uy[3]*dx11+0.8660254037844386*rhoux_2p[2]*uy[2]*dx11+0.8660254037844386*rhoux_2p[1]*uy[1]*dx11+0.8660254037844386*rhoux_2p[0]*uy[0]*dx11+1.732050807568877*Pxy[0]*dx11; 
  outrhoux[3] += 0.7745966692414834*uy[3]*rhoux_2p[6]*dx11+0.7745966692414833*uy[1]*rhoux_2p[4]*dx11+0.8660254037844386*rhoux_2p[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhoux_2p[3]*dx11+0.8660254037844386*rhoux_2p[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhoux_2p[1]*dx11+1.732050807568877*Pxy[1]*dx11+0.7745966692414834*ux[3]*rhoux_2p[7]*dx10+0.7745966692414833*ux[2]*rhoux_2p[5]*dx10+0.8660254037844386*rhoux_2p[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhoux_2p[3]*dx10+0.8660254037844386*rhoux_2p[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhoux_2p[2]*dx10+1.732050807568877*Pxx[2]*dx10; 

  outrhouy[1] += 0.8660254037844386*rhouy_2p[3]*ux[3]*dx10+0.8660254037844386*rhouy_2p[2]*ux[2]*dx10+0.8660254037844386*rhouy_2p[1]*ux[1]*dx10+0.8660254037844386*rhouy_2p[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 
  outrhouy[2] += 0.8660254037844386*rhouy_2p[3]*uy[3]*dx11+0.8660254037844386*rhouy_2p[2]*uy[2]*dx11+0.8660254037844386*rhouy_2p[1]*uy[1]*dx11+0.8660254037844386*rhouy_2p[0]*uy[0]*dx11+1.732050807568877*Pyy[0]*dx11; 
  outrhouy[3] += 0.7745966692414834*uy[3]*rhouy_2p[6]*dx11+0.7745966692414833*uy[1]*rhouy_2p[4]*dx11+0.8660254037844386*rhouy_2p[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouy_2p[3]*dx11+0.8660254037844386*rhouy_2p[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouy_2p[1]*dx11+1.732050807568877*Pyy[1]*dx11+0.7745966692414834*ux[3]*rhouy_2p[7]*dx10+0.7745966692414833*ux[2]*rhouy_2p[5]*dx10+0.8660254037844386*rhouy_2p[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouy_2p[3]*dx10+0.8660254037844386*rhouy_2p[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouy_2p[2]*dx10+1.732050807568877*Pxy[2]*dx10; 

  outrhouz[1] += 0.8660254037844386*rhouz_2p[3]*ux[3]*dx10+0.8660254037844386*rhouz_2p[2]*ux[2]*dx10+0.8660254037844386*rhouz_2p[1]*ux[1]*dx10+0.8660254037844386*rhouz_2p[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 
  outrhouz[2] += 0.8660254037844386*rhouz_2p[3]*uy[3]*dx11+0.8660254037844386*rhouz_2p[2]*uy[2]*dx11+0.8660254037844386*rhouz_2p[1]*uy[1]*dx11+0.8660254037844386*rhouz_2p[0]*uy[0]*dx11+1.732050807568877*Pyz[0]*dx11; 
  outrhouz[3] += 0.7745966692414834*uy[3]*rhouz_2p[6]*dx11+0.7745966692414833*uy[1]*rhouz_2p[4]*dx11+0.8660254037844386*rhouz_2p[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouz_2p[3]*dx11+0.8660254037844386*rhouz_2p[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouz_2p[1]*dx11+1.732050807568877*Pyz[1]*dx11+0.7745966692414834*ux[3]*rhouz_2p[7]*dx10+0.7745966692414833*ux[2]*rhouz_2p[5]*dx10+0.8660254037844386*rhouz_2p[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouz_2p[3]*dx10+0.8660254037844386*rhouz_2p[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouz_2p[2]*dx10+1.732050807568877*Pxz[2]*dx10; 

  return cflFreq_mid; 
} 
