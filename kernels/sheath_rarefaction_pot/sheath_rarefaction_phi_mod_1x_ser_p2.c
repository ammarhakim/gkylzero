#include <gkyl_sheath_rarefaction_pot_kernels.h> 


GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_1x_ser_p2(double elem_q, double mElc, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1, m2par (and m2perp, but not used) moments of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1, m2par (and m2perp, but not used) moments of ions.
  // phi: electrostatic potential.

  double phiNodal[3];
  phiNodal[0] = 1.5811388300841895*phi[2]-1.224744871391589*phi[1]+0.7071067811865475*phi[0]; 
  phiNodal[1] = 0.7071067811865475*phi[0]-0.7905694150420947*phi[2]; 
  phiNodal[2] = 1.5811388300841895*phi[2]+1.224744871391589*phi[1]+0.7071067811865475*phi[0]; 

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[3];
  const double *m1Ion = &momsIon[3];
  const double *m2parElc = &momsElc[6];
  const double *m2parIon = &momsIon[6];

  double m0eSurfNodal[1];
  double m0iSurfNodal[1];
  double m1eSurfNodal[1];
  double m1iSurfNodal[1];
  double m2pareSurfNodal[1];
  double m2pariSurfNodal[1];
  double phiWallNodal[1];
  m0eSurfNodal[0] = 1.5811388300841895*m0Elc[2]-1.224744871391589*m0Elc[1]+0.7071067811865475*m0Elc[0]; 
  m0iSurfNodal[0] = 1.5811388300841895*m0Ion[2]-1.224744871391589*m0Ion[1]+0.7071067811865475*m0Ion[0]; 
  m1eSurfNodal[0] = 1.5811388300841895*m1Elc[2]-1.224744871391589*m1Elc[1]+0.7071067811865475*m1Elc[0]; 
  m1iSurfNodal[0] = 1.5811388300841895*m1Ion[2]-1.224744871391589*m1Ion[1]+0.7071067811865475*m1Ion[0]; 
  m2pareSurfNodal[0] = 1.5811388300841895*m2parElc[2]-1.224744871391589*m2parElc[1]+0.7071067811865475*m2parElc[0]; 
  m2pariSurfNodal[0] = 1.5811388300841895*m2parIon[2]-1.224744871391589*m2parIon[1]+0.7071067811865475*m2parIon[0]; 
  phiWallNodal[0] = phiWall[0]; 

  double upar_e;
  double upar_i;
  double Tpar_e;
  double Tpar_i;
  double c_s;
  double Delta_phi_r;

  upar_e = m1eSurfNodal[0]/m0eSurfNodal[0];
  upar_i = m1iSurfNodal[0]/m0iSurfNodal[0];

  Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[0]*mElc*upar_e-1.0*m2pareSurfNodal[0]*mElc))/m0eSurfNodal[0]));
  Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[0]*mIon*upar_i-1.0*m2pariSurfNodal[0]*mIon))/m0iSurfNodal[0]));

  c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  if (c_s > 0.0)
    Delta_phi_r = ((fmin(1.0,fabs(upar_i)/c_s)-1.0)*Tpar_e)/elem_q;
  else
    Delta_phi_r = 0.0;

  phiNodal[0] = fmax(phiNodal[0] + Delta_phi_r, phiWallNodal[0]);

  phi[0] = 0.2357022603955158*phiNodal[2]+0.9428090415820636*phiNodal[1]+0.2357022603955158*phiNodal[0]; 
  phi[1] = 0.408248290463863*phiNodal[2]-0.408248290463863*phiNodal[0]; 
  phi[2] = 0.21081851067789195*phiNodal[2]-0.421637021355784*phiNodal[1]+0.21081851067789195*phiNodal[0]; 
 
}

GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_1x_ser_p2(double elem_q, double mElc, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1, m2par (and m2perp, but not used) moments of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1, m2par (and m2perp, but not used) moments of ions.
  // phi: electrostatic potential.

  double phiNodal[3];
  phiNodal[0] = 1.5811388300841895*phi[2]-1.224744871391589*phi[1]+0.7071067811865475*phi[0]; 
  phiNodal[1] = 0.7071067811865475*phi[0]-0.7905694150420947*phi[2]; 
  phiNodal[2] = 1.5811388300841895*phi[2]+1.224744871391589*phi[1]+0.7071067811865475*phi[0]; 

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[3];
  const double *m1Ion = &momsIon[3];
  const double *m2parElc = &momsElc[6];
  const double *m2parIon = &momsIon[6];

  double m0eSurfNodal[1];
  double m0iSurfNodal[1];
  double m1eSurfNodal[1];
  double m1iSurfNodal[1];
  double m2pareSurfNodal[1];
  double m2pariSurfNodal[1];
  double phiWallNodal[1];
  m0eSurfNodal[0] = 1.5811388300841895*m0Elc[2]+1.224744871391589*m0Elc[1]+0.7071067811865475*m0Elc[0]; 
  m0iSurfNodal[0] = 1.5811388300841895*m0Ion[2]+1.224744871391589*m0Ion[1]+0.7071067811865475*m0Ion[0]; 
  m1eSurfNodal[0] = 1.5811388300841895*m1Elc[2]+1.224744871391589*m1Elc[1]+0.7071067811865475*m1Elc[0]; 
  m1iSurfNodal[0] = 1.5811388300841895*m1Ion[2]+1.224744871391589*m1Ion[1]+0.7071067811865475*m1Ion[0]; 
  m2pareSurfNodal[0] = 1.5811388300841895*m2parElc[2]+1.224744871391589*m2parElc[1]+0.7071067811865475*m2parElc[0]; 
  m2pariSurfNodal[0] = 1.5811388300841895*m2parIon[2]+1.224744871391589*m2parIon[1]+0.7071067811865475*m2parIon[0]; 
  phiWallNodal[0] = phiWall[0]; 

  double upar_e;
  double upar_i;
  double Tpar_e;
  double Tpar_i;
  double c_s;
  double Delta_phi_r;

  upar_e = m1eSurfNodal[0]/m0eSurfNodal[0];
  upar_i = m1iSurfNodal[0]/m0iSurfNodal[0];

  Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[0]*mElc*upar_e-1.0*m2pareSurfNodal[0]*mElc))/m0eSurfNodal[0]));
  Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[0]*mIon*upar_i-1.0*m2pariSurfNodal[0]*mIon))/m0iSurfNodal[0]));

  c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  if (c_s > 0.0)
    Delta_phi_r = ((fmin(1.0,fabs(upar_i)/c_s)-1.0)*Tpar_e)/elem_q;
  else
    Delta_phi_r = 0.0;

  phiNodal[2] = fmax(phiNodal[2] + Delta_phi_r, phiWallNodal[0]);

  phi[0] = 0.2357022603955158*phiNodal[2]+0.9428090415820636*phiNodal[1]+0.2357022603955158*phiNodal[0]; 
  phi[1] = 0.408248290463863*phiNodal[2]-0.408248290463863*phiNodal[0]; 
  phi[2] = 0.21081851067789195*phiNodal[2]-0.421637021355784*phiNodal[1]+0.21081851067789195*phiNodal[0]; 
 
}
