#include <gkyl_sheath_rarefaction_pot_kernels.h> 


GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_1x_tensor_p2(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1 and m2 moments of electrons.
  // m2parElc: parallel m2 moment of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1 and m2 moments of ions.
  // m2parIon: parallel m2 moment of ions.
  // phi: electrostatic potential.

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[3];
  const double *m1Ion = &momsIon[3];

  double m0_e = 1.58113883008419*m0Elc[2]-1.224744871391589*m0Elc[1]+0.7071067811865475*m0Elc[0];
  double m0_i = 1.58113883008419*m0Ion[2]-1.224744871391589*m0Ion[1]+0.7071067811865475*m0Ion[0];
  double m1_e = 1.58113883008419*m1Elc[2]-1.224744871391589*m1Elc[1]+0.7071067811865475*m1Elc[0];
  double m1_i = 1.58113883008419*m1Ion[2]-1.224744871391589*m1Ion[1]+0.7071067811865475*m1Ion[0];
  double m2par_e = 1.58113883008419*m2parElc[2]-1.224744871391589*m2parElc[1]+0.7071067811865475*m2parElc[0];
  double m2par_i = 1.58113883008419*m2parIon[2]-1.224744871391589*m2parIon[1]+0.7071067811865475*m2parIon[0];

  double upar_e = m1_e/m0_e;
  double upar_i = m1_i/m0_i;

  double Tpar_e = -(1.0*(m1_e*mElc*upar_e-1.0*m2par_e*mElc))/m0_e;
  double Tpar_i = -(1.0*(m1_i*mIon*upar_i-1.0*m2par_i*mIon))/m0_i;

  if (Tpar_e>0. && Tpar_i>0.) {

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);

  double Delta_phi_r = (Tpar_e*fabs(upar_i)-1.0*Tpar_e*c_s)/(c_s*elem_q);

  double phiNod[3] = {0.};
  phiNod[0] = fmax(Delta_phi_r+1.58113883008419*phi[2]-1.224744871391589*phi[1]+0.7071067811865475*phi[0],1.58113883008419*phiWall[2]-1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]); 
  phiNod[1] = 0.7071067811865475*phi[0]-0.7905694150420947*phi[2]; 
  phiNod[2] = 1.58113883008419*phi[2]+1.224744871391589*phi[1]+0.7071067811865475*phi[0]; 

  phi[0] = 0.2357022603955158*phiNod[2]+0.9428090415820636*phiNod[1]+0.2357022603955158*phiNod[0]; 
  phi[1] = 0.408248290463863*phiNod[2]-0.408248290463863*phiNod[0]; 
  phi[2] = 0.210818510677892*phiNod[2]-0.421637021355784*phiNod[1]+0.210818510677892*phiNod[0]; 

  }
 
}

GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_1x_tensor_p2(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1 and m2 moments of electrons.
  // m2parElc: parallel m2 moment of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1 and m2 moments of ions.
  // m2parIon: parallel m2 moment of ions.
  // phi: electrostatic potential.

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[3];
  const double *m1Ion = &momsIon[3];

  double m0_e = 1.58113883008419*m0Elc[2]+1.224744871391589*m0Elc[1]+0.7071067811865475*m0Elc[0];
  double m0_i = 1.58113883008419*m0Ion[2]+1.224744871391589*m0Ion[1]+0.7071067811865475*m0Ion[0];
  double m1_e = 1.58113883008419*m1Elc[2]+1.224744871391589*m1Elc[1]+0.7071067811865475*m1Elc[0];
  double m1_i = 1.58113883008419*m1Ion[2]+1.224744871391589*m1Ion[1]+0.7071067811865475*m1Ion[0];
  double m2par_e = 1.58113883008419*m2parElc[2]+1.224744871391589*m2parElc[1]+0.7071067811865475*m2parElc[0];
  double m2par_i = 1.58113883008419*m2parIon[2]+1.224744871391589*m2parIon[1]+0.7071067811865475*m2parIon[0];

  double upar_e = m1_e/m0_e;
  double upar_i = m1_i/m0_i;

  double Tpar_e = -(1.0*(m1_e*mElc*upar_e-1.0*m2par_e*mElc))/m0_e;
  double Tpar_i = -(1.0*(m1_i*mIon*upar_i-1.0*m2par_i*mIon))/m0_i;

  if (Tpar_e>0. && Tpar_i>0.) {

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);

  double Delta_phi_r = (Tpar_e*fabs(upar_i)-1.0*Tpar_e*c_s)/(c_s*elem_q);

  double phiNod[3] = {0.};
  phiNod[0] = 1.58113883008419*phi[2]-1.224744871391589*phi[1]+0.7071067811865475*phi[0]; 
  phiNod[1] = 0.7071067811865475*phi[0]-0.7905694150420947*phi[2]; 
  phiNod[2] = fmax(Delta_phi_r+1.58113883008419*phi[2]+1.224744871391589*phi[1]+0.7071067811865475*phi[0],1.58113883008419*phiWall[2]+1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]); 

  phi[0] = 0.2357022603955158*phiNod[2]+0.9428090415820636*phiNod[1]+0.2357022603955158*phiNod[0]; 
  phi[1] = 0.408248290463863*phiNod[2]-0.408248290463863*phiNod[0]; 
  phi[2] = 0.210818510677892*phiNod[2]-0.421637021355784*phiNod[1]+0.210818510677892*phiNod[0]; 

  }
 
}
