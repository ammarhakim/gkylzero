#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p3_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p3(const double *w, const double *dxv, 
     const double *u_i, const double *p_ij, const double *bvar, const double *rho_inv_b, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // p_ij:     pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[4]; 
  const double *rho_inv_bz = &rho_inv_b[8]; 

  double alphaSurf_l[4] = {0.0}; 
  alphaSurf_l[0] = ((-4.183300132670378*bxbz[2]*uz[3])-1.870828693386971*bxbz[0]*uz[3]-4.183300132670378*bxby[2]*uy[3]-1.870828693386971*bxby[0]*uy[3]-4.183300132670378*bxbx[2]*ux[3]-1.870828693386971*bxbx[0]*ux[3]-2.738612787525831*bxbz[1]*uz[2]-2.738612787525831*bxby[1]*uy[2]-2.738612787525831*bxbx[1]*ux[2]-1.224744871391589*bxbz[0]*uz[1]-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+(2.091650066335188*bxbz[2]*uz[3]+0.9354143466934851*bxbz[0]*uz[3]+2.091650066335188*bxby[2]*uy[3]+0.9354143466934851*bxby[0]*uy[3]+2.091650066335188*bxbx[2]*ux[3]+0.9354143466934851*bxbx[0]*ux[3]+1.369306393762915*bxbz[1]*uz[2]+1.369306393762915*bxby[1]*uy[2]+1.369306393762915*bxbx[1]*ux[2]+0.6123724356957944*bxbz[0]*uz[1]+0.6123724356957944*bxby[0]*uy[1]+0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+(4.183300132670378*rho_inv_bz[2]*Pxz[3]+1.870828693386971*rho_inv_bz[0]*Pxz[3]+4.183300132670378*rho_inv_by[2]*Pxy[3]+1.870828693386971*rho_inv_by[0]*Pxy[3]+4.183300132670378*rho_inv_bx[2]*Pxx[3]+1.870828693386971*rho_inv_bx[0]*Pxx[3]+2.738612787525831*rho_inv_bz[1]*Pxz[2]+2.738612787525831*rho_inv_by[1]*Pxy[2]+2.738612787525831*rho_inv_bx[1]*Pxx[2]+1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf_l[1] = ((-3.674234614174766*bxbz[3]*uz[3])-5.612486080160912*bxbz[1]*uz[3]-3.674234614174766*bxby[3]*uy[3]-5.612486080160912*bxby[1]*uy[3]-3.674234614174766*bxbx[3]*ux[3]-5.612486080160912*bxbx[1]*ux[3]-2.449489742783178*bxbz[2]*uz[2]-2.738612787525831*bxbz[0]*uz[2]-2.449489742783178*bxby[2]*uy[2]-2.738612787525831*bxby[0]*uy[2]-2.449489742783178*bxbx[2]*ux[2]-2.738612787525831*bxbx[0]*ux[2]-1.224744871391589*bxbz[1]*uz[1]-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+(1.837117307087383*bxbz[3]*uz[3]+2.806243040080455*bxbz[1]*uz[3]+1.837117307087383*bxby[3]*uy[3]+2.806243040080455*bxby[1]*uy[3]+1.837117307087383*bxbx[3]*ux[3]+2.806243040080455*bxbx[1]*ux[3]+1.224744871391589*bxbz[2]*uz[2]+1.369306393762915*bxbz[0]*uz[2]+1.224744871391589*bxby[2]*uy[2]+1.369306393762915*bxby[0]*uy[2]+1.224744871391589*bxbx[2]*ux[2]+1.369306393762915*bxbx[0]*ux[2]+0.6123724356957944*bxbz[1]*uz[1]+0.6123724356957944*bxby[1]*uy[1]+0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+(3.674234614174766*Pxz[3]*rho_inv_bz[3]+3.674234614174766*Pxy[3]*rho_inv_by[3]+3.674234614174766*Pxx[3]*rho_inv_bx[3]+5.612486080160912*rho_inv_bz[1]*Pxz[3]+5.612486080160912*rho_inv_by[1]*Pxy[3]+5.612486080160912*rho_inv_bx[1]*Pxx[3]+2.449489742783178*Pxz[2]*rho_inv_bz[2]+2.449489742783178*Pxy[2]*rho_inv_by[2]+2.449489742783178*Pxx[2]*rho_inv_bx[2]+2.738612787525831*rho_inv_bz[0]*Pxz[2]+2.738612787525831*rho_inv_by[0]*Pxy[2]+2.738612787525831*rho_inv_bx[0]*Pxx[2]+1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dx0; 
  alphaSurf_l[2] = ((-4.543441112511214*bxbz[2]*uz[3])-4.183300132670378*bxbz[0]*uz[3]-4.543441112511214*bxby[2]*uy[3]-4.183300132670378*bxby[0]*uy[3]-4.543441112511214*bxbx[2]*ux[3]-4.183300132670378*bxbx[0]*ux[3]-2.405351177211819*uz[2]*bxbz[3]-2.405351177211819*uy[2]*bxby[3]-2.405351177211819*ux[2]*bxbx[3]-2.449489742783178*bxbz[1]*uz[2]-2.449489742783178*bxby[1]*uy[2]-2.449489742783178*bxbx[1]*ux[2]-1.224744871391589*uz[1]*bxbz[2]-1.224744871391589*uy[1]*bxby[2]-1.224744871391589*ux[1]*bxbx[2])*dx0*wvpar+(2.271720556255607*bxbz[2]*uz[3]+2.091650066335188*bxbz[0]*uz[3]+2.271720556255607*bxby[2]*uy[3]+2.091650066335188*bxby[0]*uy[3]+2.271720556255607*bxbx[2]*ux[3]+2.091650066335188*bxbx[0]*ux[3]+1.202675588605909*uz[2]*bxbz[3]+1.202675588605909*uy[2]*bxby[3]+1.202675588605909*ux[2]*bxbx[3]+1.224744871391589*bxbz[1]*uz[2]+1.224744871391589*bxby[1]*uy[2]+1.224744871391589*bxbx[1]*ux[2]+0.6123724356957944*uz[1]*bxbz[2]+0.6123724356957944*uy[1]*bxby[2]+0.6123724356957944*ux[1]*bxbx[2])*dvpar*dx0+(2.405351177211819*Pxz[2]*rho_inv_bz[3]+2.405351177211819*Pxy[2]*rho_inv_by[3]+2.405351177211819*Pxx[2]*rho_inv_bx[3]+4.543441112511214*rho_inv_bz[2]*Pxz[3]+4.183300132670378*rho_inv_bz[0]*Pxz[3]+4.543441112511214*rho_inv_by[2]*Pxy[3]+4.183300132670378*rho_inv_by[0]*Pxy[3]+4.543441112511214*rho_inv_bx[2]*Pxx[3]+4.183300132670378*rho_inv_bx[0]*Pxx[3]+1.224744871391589*Pxz[1]*rho_inv_bz[2]+1.224744871391589*Pxy[1]*rho_inv_by[2]+1.224744871391589*Pxx[1]*rho_inv_bx[2]+2.449489742783178*rho_inv_bz[1]*Pxz[2]+2.449489742783178*rho_inv_by[1]*Pxy[2]+2.449489742783178*rho_inv_bx[1]*Pxx[2])*dx0; 
  alphaSurf_l[3] = ((-4.365266951236265*bxbz[3]*uz[3])-3.674234614174766*bxbz[1]*uz[3]-4.365266951236265*bxby[3]*uy[3]-3.674234614174766*bxby[1]*uy[3]-4.365266951236265*bxbx[3]*ux[3]-3.674234614174766*bxbx[1]*ux[3]-1.224744871391589*uz[1]*bxbz[3]-1.224744871391589*uy[1]*bxby[3]-1.224744871391589*ux[1]*bxbx[3]-2.405351177211819*bxbz[2]*uz[2]-2.405351177211819*bxby[2]*uy[2]-2.405351177211819*bxbx[2]*ux[2])*dx0*wvpar+(2.182633475618132*bxbz[3]*uz[3]+1.837117307087383*bxbz[1]*uz[3]+2.182633475618132*bxby[3]*uy[3]+1.837117307087383*bxby[1]*uy[3]+2.182633475618132*bxbx[3]*ux[3]+1.837117307087383*bxbx[1]*ux[3]+0.6123724356957944*uz[1]*bxbz[3]+0.6123724356957944*uy[1]*bxby[3]+0.6123724356957944*ux[1]*bxbx[3]+1.202675588605909*bxbz[2]*uz[2]+1.202675588605909*bxby[2]*uy[2]+1.202675588605909*bxbx[2]*ux[2])*dvpar*dx0+(4.365266951236265*Pxz[3]*rho_inv_bz[3]+1.224744871391589*Pxz[1]*rho_inv_bz[3]+4.365266951236265*Pxy[3]*rho_inv_by[3]+1.224744871391589*Pxy[1]*rho_inv_by[3]+4.365266951236265*Pxx[3]*rho_inv_bx[3]+1.224744871391589*Pxx[1]*rho_inv_bx[3]+3.674234614174766*rho_inv_bz[1]*Pxz[3]+3.674234614174766*rho_inv_by[1]*Pxy[3]+3.674234614174766*rho_inv_bx[1]*Pxx[3]+2.405351177211819*Pxz[2]*rho_inv_bz[2]+2.405351177211819*Pxy[2]*rho_inv_by[2]+2.405351177211819*Pxx[2]*rho_inv_bx[2])*dx0; 

  double alphaSurf_r[4] = {0.0}; 
  alphaSurf_r[0] = ((-4.183300132670378*bxbz[2]*uz[3])-1.870828693386971*bxbz[0]*uz[3]-4.183300132670378*bxby[2]*uy[3]-1.870828693386971*bxby[0]*uy[3]-4.183300132670378*bxbx[2]*ux[3]-1.870828693386971*bxbx[0]*ux[3]-2.738612787525831*bxbz[1]*uz[2]-2.738612787525831*bxby[1]*uy[2]-2.738612787525831*bxbx[1]*ux[2]-1.224744871391589*bxbz[0]*uz[1]-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+((-2.091650066335188*bxbz[2]*uz[3])-0.9354143466934851*bxbz[0]*uz[3]-2.091650066335188*bxby[2]*uy[3]-0.9354143466934851*bxby[0]*uy[3]-2.091650066335188*bxbx[2]*ux[3]-0.9354143466934851*bxbx[0]*ux[3]-1.369306393762915*bxbz[1]*uz[2]-1.369306393762915*bxby[1]*uy[2]-1.369306393762915*bxbx[1]*ux[2]-0.6123724356957944*bxbz[0]*uz[1]-0.6123724356957944*bxby[0]*uy[1]-0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+(4.183300132670378*rho_inv_bz[2]*Pxz[3]+1.870828693386971*rho_inv_bz[0]*Pxz[3]+4.183300132670378*rho_inv_by[2]*Pxy[3]+1.870828693386971*rho_inv_by[0]*Pxy[3]+4.183300132670378*rho_inv_bx[2]*Pxx[3]+1.870828693386971*rho_inv_bx[0]*Pxx[3]+2.738612787525831*rho_inv_bz[1]*Pxz[2]+2.738612787525831*rho_inv_by[1]*Pxy[2]+2.738612787525831*rho_inv_bx[1]*Pxx[2]+1.224744871391589*rho_inv_bz[0]*Pxz[1]+1.224744871391589*rho_inv_by[0]*Pxy[1]+1.224744871391589*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf_r[1] = ((-3.674234614174766*bxbz[3]*uz[3])-5.612486080160912*bxbz[1]*uz[3]-3.674234614174766*bxby[3]*uy[3]-5.612486080160912*bxby[1]*uy[3]-3.674234614174766*bxbx[3]*ux[3]-5.612486080160912*bxbx[1]*ux[3]-2.449489742783178*bxbz[2]*uz[2]-2.738612787525831*bxbz[0]*uz[2]-2.449489742783178*bxby[2]*uy[2]-2.738612787525831*bxby[0]*uy[2]-2.449489742783178*bxbx[2]*ux[2]-2.738612787525831*bxbx[0]*ux[2]-1.224744871391589*bxbz[1]*uz[1]-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+((-1.837117307087383*bxbz[3]*uz[3])-2.806243040080455*bxbz[1]*uz[3]-1.837117307087383*bxby[3]*uy[3]-2.806243040080455*bxby[1]*uy[3]-1.837117307087383*bxbx[3]*ux[3]-2.806243040080455*bxbx[1]*ux[3]-1.224744871391589*bxbz[2]*uz[2]-1.369306393762915*bxbz[0]*uz[2]-1.224744871391589*bxby[2]*uy[2]-1.369306393762915*bxby[0]*uy[2]-1.224744871391589*bxbx[2]*ux[2]-1.369306393762915*bxbx[0]*ux[2]-0.6123724356957944*bxbz[1]*uz[1]-0.6123724356957944*bxby[1]*uy[1]-0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+(3.674234614174766*Pxz[3]*rho_inv_bz[3]+3.674234614174766*Pxy[3]*rho_inv_by[3]+3.674234614174766*Pxx[3]*rho_inv_bx[3]+5.612486080160912*rho_inv_bz[1]*Pxz[3]+5.612486080160912*rho_inv_by[1]*Pxy[3]+5.612486080160912*rho_inv_bx[1]*Pxx[3]+2.449489742783178*Pxz[2]*rho_inv_bz[2]+2.449489742783178*Pxy[2]*rho_inv_by[2]+2.449489742783178*Pxx[2]*rho_inv_bx[2]+2.738612787525831*rho_inv_bz[0]*Pxz[2]+2.738612787525831*rho_inv_by[0]*Pxy[2]+2.738612787525831*rho_inv_bx[0]*Pxx[2]+1.224744871391589*Pxz[1]*rho_inv_bz[1]+1.224744871391589*Pxy[1]*rho_inv_by[1]+1.224744871391589*Pxx[1]*rho_inv_bx[1])*dx0; 
  alphaSurf_r[2] = ((-4.543441112511214*bxbz[2]*uz[3])-4.183300132670378*bxbz[0]*uz[3]-4.543441112511214*bxby[2]*uy[3]-4.183300132670378*bxby[0]*uy[3]-4.543441112511214*bxbx[2]*ux[3]-4.183300132670378*bxbx[0]*ux[3]-2.405351177211819*uz[2]*bxbz[3]-2.405351177211819*uy[2]*bxby[3]-2.405351177211819*ux[2]*bxbx[3]-2.449489742783178*bxbz[1]*uz[2]-2.449489742783178*bxby[1]*uy[2]-2.449489742783178*bxbx[1]*ux[2]-1.224744871391589*uz[1]*bxbz[2]-1.224744871391589*uy[1]*bxby[2]-1.224744871391589*ux[1]*bxbx[2])*dx0*wvpar+((-2.271720556255607*bxbz[2]*uz[3])-2.091650066335188*bxbz[0]*uz[3]-2.271720556255607*bxby[2]*uy[3]-2.091650066335188*bxby[0]*uy[3]-2.271720556255607*bxbx[2]*ux[3]-2.091650066335188*bxbx[0]*ux[3]-1.202675588605909*uz[2]*bxbz[3]-1.202675588605909*uy[2]*bxby[3]-1.202675588605909*ux[2]*bxbx[3]-1.224744871391589*bxbz[1]*uz[2]-1.224744871391589*bxby[1]*uy[2]-1.224744871391589*bxbx[1]*ux[2]-0.6123724356957944*uz[1]*bxbz[2]-0.6123724356957944*uy[1]*bxby[2]-0.6123724356957944*ux[1]*bxbx[2])*dvpar*dx0+(2.405351177211819*Pxz[2]*rho_inv_bz[3]+2.405351177211819*Pxy[2]*rho_inv_by[3]+2.405351177211819*Pxx[2]*rho_inv_bx[3]+4.543441112511214*rho_inv_bz[2]*Pxz[3]+4.183300132670378*rho_inv_bz[0]*Pxz[3]+4.543441112511214*rho_inv_by[2]*Pxy[3]+4.183300132670378*rho_inv_by[0]*Pxy[3]+4.543441112511214*rho_inv_bx[2]*Pxx[3]+4.183300132670378*rho_inv_bx[0]*Pxx[3]+1.224744871391589*Pxz[1]*rho_inv_bz[2]+1.224744871391589*Pxy[1]*rho_inv_by[2]+1.224744871391589*Pxx[1]*rho_inv_bx[2]+2.449489742783178*rho_inv_bz[1]*Pxz[2]+2.449489742783178*rho_inv_by[1]*Pxy[2]+2.449489742783178*rho_inv_bx[1]*Pxx[2])*dx0; 
  alphaSurf_r[3] = ((-4.365266951236265*bxbz[3]*uz[3])-3.674234614174766*bxbz[1]*uz[3]-4.365266951236265*bxby[3]*uy[3]-3.674234614174766*bxby[1]*uy[3]-4.365266951236265*bxbx[3]*ux[3]-3.674234614174766*bxbx[1]*ux[3]-1.224744871391589*uz[1]*bxbz[3]-1.224744871391589*uy[1]*bxby[3]-1.224744871391589*ux[1]*bxbx[3]-2.405351177211819*bxbz[2]*uz[2]-2.405351177211819*bxby[2]*uy[2]-2.405351177211819*bxbx[2]*ux[2])*dx0*wvpar+((-2.182633475618132*bxbz[3]*uz[3])-1.837117307087383*bxbz[1]*uz[3]-2.182633475618132*bxby[3]*uy[3]-1.837117307087383*bxby[1]*uy[3]-2.182633475618132*bxbx[3]*ux[3]-1.837117307087383*bxbx[1]*ux[3]-0.6123724356957944*uz[1]*bxbz[3]-0.6123724356957944*uy[1]*bxby[3]-0.6123724356957944*ux[1]*bxbx[3]-1.202675588605909*bxbz[2]*uz[2]-1.202675588605909*bxby[2]*uy[2]-1.202675588605909*bxbx[2]*ux[2])*dvpar*dx0+(4.365266951236265*Pxz[3]*rho_inv_bz[3]+1.224744871391589*Pxz[1]*rho_inv_bz[3]+4.365266951236265*Pxy[3]*rho_inv_by[3]+1.224744871391589*Pxy[1]*rho_inv_by[3]+4.365266951236265*Pxx[3]*rho_inv_bx[3]+1.224744871391589*Pxx[1]*rho_inv_bx[3]+3.674234614174766*rho_inv_bz[1]*Pxz[3]+3.674234614174766*rho_inv_by[1]*Pxy[3]+3.674234614174766*rho_inv_bx[1]*Pxx[3]+2.405351177211819*Pxz[2]*rho_inv_bz[2]+2.405351177211819*Pxy[2]*rho_inv_by[2]+2.405351177211819*Pxx[2]*rho_inv_bx[2])*dx0; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if ((-0.5701294036773671*alphaSurf_l[3])+0.9681844646844029*alphaSurf_l[2]-1.054672281193885*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fc); 
  } 
  if ((-0.5701294036773671*alphaSurf_r[3])+0.9681844646844029*alphaSurf_r[2]-1.054672281193885*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.7702725556588816*alphaSurf_l[3]-0.5164305132317774*alphaSurf_l[2]-0.4163900395009129*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.7702725556588816*alphaSurf_r[3]-0.5164305132317774*alphaSurf_r[2]-0.4163900395009129*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fr); 
  } 
  if ((-0.7702725556588816*alphaSurf_l[3])-0.5164305132317774*alphaSurf_l[2]+0.4163900395009129*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fc); 
  } 
  if ((-0.7702725556588816*alphaSurf_r[3])-0.5164305132317774*alphaSurf_r[2]+0.4163900395009129*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fr); 
  } 
  if (0.5701294036773671*alphaSurf_l[3]+0.9681844646844029*alphaSurf_l[2]+1.054672281193885*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fc); 
  } 
  if (0.5701294036773671*alphaSurf_r[3]+0.9681844646844029*alphaSurf_r[2]+1.054672281193885*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaSurf_l[3]*fUpwind_l[3]+0.7071067811865475*alphaSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.6210590034081186*alphaSurf_l[2]*fUpwind_l[3]+0.6210590034081186*fUpwind_l[2]*alphaSurf_l[3]+0.6324555320336759*alphaSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaSurf_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.421637021355784*alphaSurf_l[3]*fUpwind_l[3]+0.6210590034081186*alphaSurf_l[1]*fUpwind_l[3]+0.6210590034081186*fUpwind_l[1]*alphaSurf_l[3]+0.4517539514526256*alphaSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[2]+0.6324555320336759*alphaSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[3] = 0.421637021355784*alphaSurf_l[2]*fUpwind_l[3]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[3]+0.421637021355784*fUpwind_l[2]*alphaSurf_l[3]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[3]+0.6210590034081186*alphaSurf_l[1]*fUpwind_l[2]+0.6210590034081186*fUpwind_l[1]*alphaSurf_l[2]; 

  Ghat_r[0] = 0.7071067811865475*alphaSurf_r[3]*fUpwind_r[3]+0.7071067811865475*alphaSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.6210590034081186*alphaSurf_r[2]*fUpwind_r[3]+0.6210590034081186*fUpwind_r[2]*alphaSurf_r[3]+0.6324555320336759*alphaSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaSurf_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.421637021355784*alphaSurf_r[3]*fUpwind_r[3]+0.6210590034081186*alphaSurf_r[1]*fUpwind_r[3]+0.6210590034081186*fUpwind_r[1]*alphaSurf_r[3]+0.4517539514526256*alphaSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[2]+0.6324555320336759*alphaSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[3] = 0.421637021355784*alphaSurf_r[2]*fUpwind_r[3]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[3]+0.421637021355784*fUpwind_r[2]*alphaSurf_r[3]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[3]+0.6210590034081186*alphaSurf_r[1]*fUpwind_r[2]+0.6210590034081186*fUpwind_r[1]*alphaSurf_r[2]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv1par; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv1par; 
  out[5] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv1par; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv1par; 
  out[7] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv1par; 
  out[8] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv1par; 
  out[9] += -1.870828693386971*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv1par; 
  out[11] += -1.870828693386971*(Ghat_r[1]+Ghat_l[1])*dv1par; 

} 
