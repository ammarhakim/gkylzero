#include <gkyl_sheath_rarefaction_pot_kernels.h> 


GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_2x_ser_p2(double elem_q, double mElc, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1, m2par (and m2perp, but not used) moments of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1, m2par (and m2perp, but not used) moments of ions.
  // phi: electrostatic potential.

  double phiNodal[8];
  phiNodal[0] = -(1.9364916731037085*phi[7])-1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]+1.5*phi[3]-0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[1] = 0.9682458365518543*phi[6]+1.118033988749895*phi[5]-0.5590169943749475*phi[4]-0.8660254037844386*phi[2]+0.5*phi[0]; 
  phiNodal[2] = 1.9364916731037085*phi[7]-1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]-1.5*phi[3]-0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[3] = 0.9682458365518543*phi[7]-0.5590169943749475*phi[5]+1.118033988749895*phi[4]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[4] = -(0.9682458365518543*phi[7])-0.5590169943749475*phi[5]+1.118033988749895*phi[4]+0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[5] = -(1.9364916731037085*phi[7])+1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]-1.5*phi[3]+0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[6] = -(0.9682458365518543*phi[6])+1.118033988749895*phi[5]-0.5590169943749475*phi[4]+0.8660254037844386*phi[2]+0.5*phi[0]; 
  phiNodal[7] = 1.9364916731037085*phi[7]+1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]+1.5*phi[3]+0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[8];
  const double *m1Ion = &momsIon[8];
  const double *m2parElc = &momsElc[16];
  const double *m2parIon = &momsIon[16];

  double m0eSurfNodal[3];
  double m0iSurfNodal[3];
  double m1eSurfNodal[3];
  double m1iSurfNodal[3];
  double m2pareSurfNodal[3];
  double m2pariSurfNodal[3];
  double phiWallNodal[3];
  m0eSurfNodal[0] = -(1.9364916731037085*m0Elc[7])-1.9364916731037085*m0Elc[6]+1.118033988749895*m0Elc[5]+1.118033988749895*m0Elc[4]+1.5*m0Elc[3]-0.8660254037844386*m0Elc[2]-0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0eSurfNodal[1] = 0.9682458365518543*m0Elc[6]+1.118033988749895*m0Elc[5]-0.5590169943749475*m0Elc[4]-0.8660254037844386*m0Elc[2]+0.5*m0Elc[0]; 
  m0eSurfNodal[2] = 1.9364916731037085*m0Elc[7]-1.9364916731037085*m0Elc[6]+1.118033988749895*m0Elc[5]+1.118033988749895*m0Elc[4]-1.5*m0Elc[3]-0.8660254037844386*m0Elc[2]+0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0iSurfNodal[0] = -(1.9364916731037085*m0Ion[7])-1.9364916731037085*m0Ion[6]+1.118033988749895*m0Ion[5]+1.118033988749895*m0Ion[4]+1.5*m0Ion[3]-0.8660254037844386*m0Ion[2]-0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m0iSurfNodal[1] = 0.9682458365518543*m0Ion[6]+1.118033988749895*m0Ion[5]-0.5590169943749475*m0Ion[4]-0.8660254037844386*m0Ion[2]+0.5*m0Ion[0]; 
  m0iSurfNodal[2] = 1.9364916731037085*m0Ion[7]-1.9364916731037085*m0Ion[6]+1.118033988749895*m0Ion[5]+1.118033988749895*m0Ion[4]-1.5*m0Ion[3]-0.8660254037844386*m0Ion[2]+0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m1eSurfNodal[0] = -(1.9364916731037085*m1Elc[7])-1.9364916731037085*m1Elc[6]+1.118033988749895*m1Elc[5]+1.118033988749895*m1Elc[4]+1.5*m1Elc[3]-0.8660254037844386*m1Elc[2]-0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1eSurfNodal[1] = 0.9682458365518543*m1Elc[6]+1.118033988749895*m1Elc[5]-0.5590169943749475*m1Elc[4]-0.8660254037844386*m1Elc[2]+0.5*m1Elc[0]; 
  m1eSurfNodal[2] = 1.9364916731037085*m1Elc[7]-1.9364916731037085*m1Elc[6]+1.118033988749895*m1Elc[5]+1.118033988749895*m1Elc[4]-1.5*m1Elc[3]-0.8660254037844386*m1Elc[2]+0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1iSurfNodal[0] = -(1.9364916731037085*m1Ion[7])-1.9364916731037085*m1Ion[6]+1.118033988749895*m1Ion[5]+1.118033988749895*m1Ion[4]+1.5*m1Ion[3]-0.8660254037844386*m1Ion[2]-0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m1iSurfNodal[1] = 0.9682458365518543*m1Ion[6]+1.118033988749895*m1Ion[5]-0.5590169943749475*m1Ion[4]-0.8660254037844386*m1Ion[2]+0.5*m1Ion[0]; 
  m1iSurfNodal[2] = 1.9364916731037085*m1Ion[7]-1.9364916731037085*m1Ion[6]+1.118033988749895*m1Ion[5]+1.118033988749895*m1Ion[4]-1.5*m1Ion[3]-0.8660254037844386*m1Ion[2]+0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m2pareSurfNodal[0] = -(1.9364916731037085*m2parElc[7])-1.9364916731037085*m2parElc[6]+1.118033988749895*m2parElc[5]+1.118033988749895*m2parElc[4]+1.5*m2parElc[3]-0.8660254037844386*m2parElc[2]-0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pareSurfNodal[1] = 0.9682458365518543*m2parElc[6]+1.118033988749895*m2parElc[5]-0.5590169943749475*m2parElc[4]-0.8660254037844386*m2parElc[2]+0.5*m2parElc[0]; 
  m2pareSurfNodal[2] = 1.9364916731037085*m2parElc[7]-1.9364916731037085*m2parElc[6]+1.118033988749895*m2parElc[5]+1.118033988749895*m2parElc[4]-1.5*m2parElc[3]-0.8660254037844386*m2parElc[2]+0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pariSurfNodal[0] = -(1.9364916731037085*m2parIon[7])-1.9364916731037085*m2parIon[6]+1.118033988749895*m2parIon[5]+1.118033988749895*m2parIon[4]+1.5*m2parIon[3]-0.8660254037844386*m2parIon[2]-0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  m2pariSurfNodal[1] = 0.9682458365518543*m2parIon[6]+1.118033988749895*m2parIon[5]-0.5590169943749475*m2parIon[4]-0.8660254037844386*m2parIon[2]+0.5*m2parIon[0]; 
  m2pariSurfNodal[2] = 1.9364916731037085*m2parIon[7]-1.9364916731037085*m2parIon[6]+1.118033988749895*m2parIon[5]+1.118033988749895*m2parIon[4]-1.5*m2parIon[3]-0.8660254037844386*m2parIon[2]+0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  phiWallNodal[0] = 1.5811388300841895*phiWall[2]-1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]; 
  phiWallNodal[1] = 0.7071067811865475*phiWall[0]-0.7905694150420947*phiWall[2]; 
  phiWallNodal[2] = 1.5811388300841895*phiWall[2]+1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]; 

  double upar_e = m1eSurfNodal[0]/m0eSurfNodal[0];
  double upar_i = m1iSurfNodal[0]/m0iSurfNodal[0];

  double Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[0]*mElc*upar_e-1.0*m2pareSurfNodal[0]*mElc))/m0eSurfNodal[0]));
  double Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[0]*mIon*upar_i-1.0*m2pariSurfNodal[0]*mIon))/m0iSurfNodal[0]));

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  double Delta_phi_r;
  if (c_s > 0.0)
    Delta_phi_r = 0;
  else
    Delta_phi_r = 0.0;

  phiNodal[0] = fmax(phiNodal[0] + Delta_phi_r, phi_wallNodal[0];

  double upar_e = m1eSurfNodal[1]/m0eSurfNodal[1];
  double upar_i = m1iSurfNodal[1]/m0iSurfNodal[1];

  double Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[1]*mElc*upar_e-1.0*m2pareSurfNodal[1]*mElc))/m0eSurfNodal[1]));
  double Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[1]*mIon*upar_i-1.0*m2pariSurfNodal[1]*mIon))/m0iSurfNodal[1]));

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  double Delta_phi_r;
  if (c_s > 0.0)
    Delta_phi_r = 1;
  else
    Delta_phi_r = 0.0;

  phiNodal[1] = fmax(phiNodal[1] + Delta_phi_r, phi_wallNodal[1];

  double upar_e = m1eSurfNodal[2]/m0eSurfNodal[2];
  double upar_i = m1iSurfNodal[2]/m0iSurfNodal[2];

  double Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[2]*mElc*upar_e-1.0*m2pareSurfNodal[2]*mElc))/m0eSurfNodal[2]));
  double Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[2]*mIon*upar_i-1.0*m2pariSurfNodal[2]*mIon))/m0iSurfNodal[2]));

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  double Delta_phi_r;
  if (c_s > 0.0)
    Delta_phi_r = 2;
  else
    Delta_phi_r = 0.0;

  phiNodal[2] = fmax(phiNodal[2] + Delta_phi_r, phi_wallNodal[2];

  phi[0] = -(0.16666666666666666*phiNodal[7])+0.6666666666666666*phiNodal[6]-0.16666666666666666*phiNodal[5]+0.6666666666666666*phiNodal[4]+0.6666666666666666*phiNodal[3]-0.16666666666666666*phiNodal[2]+0.6666666666666666*phiNodal[1]-0.16666666666666666*phiNodal[0]; 
  phi[1] = 0.09622504486493764*phiNodal[7]-0.09622504486493764*phiNodal[5]+0.3849001794597506*phiNodal[4]-0.3849001794597506*phiNodal[3]+0.09622504486493764*phiNodal[2]-0.09622504486493764*phiNodal[0]; 
  phi[2] = 0.09622504486493764*phiNodal[7]+0.3849001794597506*phiNodal[6]+0.09622504486493764*phiNodal[5]-0.09622504486493764*phiNodal[2]-0.3849001794597506*phiNodal[1]-0.09622504486493764*phiNodal[0]; 
  phi[3] = 0.16666666666666666*phiNodal[7]-0.16666666666666666*phiNodal[5]-0.16666666666666666*phiNodal[2]+0.16666666666666666*phiNodal[0]; 
  phi[4] = 0.14907119849998596*phiNodal[7]-0.2981423969999719*phiNodal[6]+0.14907119849998596*phiNodal[5]+0.14907119849998596*phiNodal[2]-0.2981423969999719*phiNodal[1]+0.14907119849998596*phiNodal[0]; 
  phi[5] = 0.14907119849998596*phiNodal[7]+0.14907119849998596*phiNodal[5]-0.2981423969999719*phiNodal[4]-0.2981423969999719*phiNodal[3]+0.14907119849998596*phiNodal[2]+0.14907119849998596*phiNodal[0]; 
  phi[6] = 0.08606629658238703*phiNodal[7]-0.17213259316477406*phiNodal[6]+0.08606629658238703*phiNodal[5]-0.08606629658238703*phiNodal[2]+0.17213259316477406*phiNodal[1]-0.08606629658238703*phiNodal[0]; 
  phi[7] = 0.08606629658238703*phiNodal[7]-0.08606629658238703*phiNodal[5]-0.17213259316477406*phiNodal[4]+0.17213259316477406*phiNodal[3]+0.08606629658238703*phiNodal[2]-0.08606629658238703*phiNodal[0]; 
 
}

GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_2x_ser_p2(double elem_q, double mElc, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1, m2par (and m2perp, but not used) moments of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1, m2par (and m2perp, but not used) moments of ions.
  // phi: electrostatic potential.

  double phiNodal[8];
  phiNodal[0] = -(1.9364916731037085*phi[7])-1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]+1.5*phi[3]-0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[1] = 0.9682458365518543*phi[6]+1.118033988749895*phi[5]-0.5590169943749475*phi[4]-0.8660254037844386*phi[2]+0.5*phi[0]; 
  phiNodal[2] = 1.9364916731037085*phi[7]-1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]-1.5*phi[3]-0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[3] = 0.9682458365518543*phi[7]-0.5590169943749475*phi[5]+1.118033988749895*phi[4]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[4] = -(0.9682458365518543*phi[7])-0.5590169943749475*phi[5]+1.118033988749895*phi[4]+0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[5] = -(1.9364916731037085*phi[7])+1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]-1.5*phi[3]+0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[6] = -(0.9682458365518543*phi[6])+1.118033988749895*phi[5]-0.5590169943749475*phi[4]+0.8660254037844386*phi[2]+0.5*phi[0]; 
  phiNodal[7] = 1.9364916731037085*phi[7]+1.9364916731037085*phi[6]+1.118033988749895*phi[5]+1.118033988749895*phi[4]+1.5*phi[3]+0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[8];
  const double *m1Ion = &momsIon[8];
  const double *m2parElc = &momsElc[16];
  const double *m2parIon = &momsIon[16];

  double m0eSurfNodal[3];
  double m0iSurfNodal[3];
  double m1eSurfNodal[3];
  double m1iSurfNodal[3];
  double m2pareSurfNodal[3];
  double m2pariSurfNodal[3];
  double phiWallNodal[3];
  m0eSurfNodal[0] = -(1.9364916731037085*m0Elc[7])+1.9364916731037085*m0Elc[6]+1.118033988749895*m0Elc[5]+1.118033988749895*m0Elc[4]-1.5*m0Elc[3]+0.8660254037844386*m0Elc[2]-0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0eSurfNodal[1] = -(0.9682458365518543*m0Elc[6])+1.118033988749895*m0Elc[5]-0.5590169943749475*m0Elc[4]+0.8660254037844386*m0Elc[2]+0.5*m0Elc[0]; 
  m0eSurfNodal[2] = 1.9364916731037085*m0Elc[7]+1.9364916731037085*m0Elc[6]+1.118033988749895*m0Elc[5]+1.118033988749895*m0Elc[4]+1.5*m0Elc[3]+0.8660254037844386*m0Elc[2]+0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0iSurfNodal[0] = -(1.9364916731037085*m0Ion[7])+1.9364916731037085*m0Ion[6]+1.118033988749895*m0Ion[5]+1.118033988749895*m0Ion[4]-1.5*m0Ion[3]+0.8660254037844386*m0Ion[2]-0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m0iSurfNodal[1] = -(0.9682458365518543*m0Ion[6])+1.118033988749895*m0Ion[5]-0.5590169943749475*m0Ion[4]+0.8660254037844386*m0Ion[2]+0.5*m0Ion[0]; 
  m0iSurfNodal[2] = 1.9364916731037085*m0Ion[7]+1.9364916731037085*m0Ion[6]+1.118033988749895*m0Ion[5]+1.118033988749895*m0Ion[4]+1.5*m0Ion[3]+0.8660254037844386*m0Ion[2]+0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m1eSurfNodal[0] = -(1.9364916731037085*m1Elc[7])+1.9364916731037085*m1Elc[6]+1.118033988749895*m1Elc[5]+1.118033988749895*m1Elc[4]-1.5*m1Elc[3]+0.8660254037844386*m1Elc[2]-0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1eSurfNodal[1] = -(0.9682458365518543*m1Elc[6])+1.118033988749895*m1Elc[5]-0.5590169943749475*m1Elc[4]+0.8660254037844386*m1Elc[2]+0.5*m1Elc[0]; 
  m1eSurfNodal[2] = 1.9364916731037085*m1Elc[7]+1.9364916731037085*m1Elc[6]+1.118033988749895*m1Elc[5]+1.118033988749895*m1Elc[4]+1.5*m1Elc[3]+0.8660254037844386*m1Elc[2]+0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1iSurfNodal[0] = -(1.9364916731037085*m1Ion[7])+1.9364916731037085*m1Ion[6]+1.118033988749895*m1Ion[5]+1.118033988749895*m1Ion[4]-1.5*m1Ion[3]+0.8660254037844386*m1Ion[2]-0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m1iSurfNodal[1] = -(0.9682458365518543*m1Ion[6])+1.118033988749895*m1Ion[5]-0.5590169943749475*m1Ion[4]+0.8660254037844386*m1Ion[2]+0.5*m1Ion[0]; 
  m1iSurfNodal[2] = 1.9364916731037085*m1Ion[7]+1.9364916731037085*m1Ion[6]+1.118033988749895*m1Ion[5]+1.118033988749895*m1Ion[4]+1.5*m1Ion[3]+0.8660254037844386*m1Ion[2]+0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m2pareSurfNodal[0] = -(1.9364916731037085*m2parElc[7])+1.9364916731037085*m2parElc[6]+1.118033988749895*m2parElc[5]+1.118033988749895*m2parElc[4]-1.5*m2parElc[3]+0.8660254037844386*m2parElc[2]-0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pareSurfNodal[1] = -(0.9682458365518543*m2parElc[6])+1.118033988749895*m2parElc[5]-0.5590169943749475*m2parElc[4]+0.8660254037844386*m2parElc[2]+0.5*m2parElc[0]; 
  m2pareSurfNodal[2] = 1.9364916731037085*m2parElc[7]+1.9364916731037085*m2parElc[6]+1.118033988749895*m2parElc[5]+1.118033988749895*m2parElc[4]+1.5*m2parElc[3]+0.8660254037844386*m2parElc[2]+0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pariSurfNodal[0] = -(1.9364916731037085*m2parIon[7])+1.9364916731037085*m2parIon[6]+1.118033988749895*m2parIon[5]+1.118033988749895*m2parIon[4]-1.5*m2parIon[3]+0.8660254037844386*m2parIon[2]-0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  m2pariSurfNodal[1] = -(0.9682458365518543*m2parIon[6])+1.118033988749895*m2parIon[5]-0.5590169943749475*m2parIon[4]+0.8660254037844386*m2parIon[2]+0.5*m2parIon[0]; 
  m2pariSurfNodal[2] = 1.9364916731037085*m2parIon[7]+1.9364916731037085*m2parIon[6]+1.118033988749895*m2parIon[5]+1.118033988749895*m2parIon[4]+1.5*m2parIon[3]+0.8660254037844386*m2parIon[2]+0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  phiWallNodal[0] = 1.5811388300841895*phiWall[2]-1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]; 
  phiWallNodal[1] = 0.7071067811865475*phiWall[0]-0.7905694150420947*phiWall[2]; 
  phiWallNodal[2] = 1.5811388300841895*phiWall[2]+1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]; 

  double upar_e = m1eSurfNodal[0]/m0eSurfNodal[0];
  double upar_i = m1iSurfNodal[0]/m0iSurfNodal[0];

  double Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[0]*mElc*upar_e-1.0*m2pareSurfNodal[0]*mElc))/m0eSurfNodal[0]));
  double Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[0]*mIon*upar_i-1.0*m2pariSurfNodal[0]*mIon))/m0iSurfNodal[0]));

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  double Delta_phi_r;
  if (c_s > 0.0)
    Delta_phi_r = 0;
  else
    Delta_phi_r = 0.0;

  phiNodal[5] = fmax(phiNodal[5] + Delta_phi_r, phi_wallNodal[0];

  double upar_e = m1eSurfNodal[1]/m0eSurfNodal[1];
  double upar_i = m1iSurfNodal[1]/m0iSurfNodal[1];

  double Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[1]*mElc*upar_e-1.0*m2pareSurfNodal[1]*mElc))/m0eSurfNodal[1]));
  double Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[1]*mIon*upar_i-1.0*m2pariSurfNodal[1]*mIon))/m0iSurfNodal[1]));

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  double Delta_phi_r;
  if (c_s > 0.0)
    Delta_phi_r = 1;
  else
    Delta_phi_r = 0.0;

  phiNodal[6] = fmax(phiNodal[6] + Delta_phi_r, phi_wallNodal[1];

  double upar_e = m1eSurfNodal[2]/m0eSurfNodal[2];
  double upar_i = m1iSurfNodal[2]/m0iSurfNodal[2];

  double Tpar_e = fmax(0.0, -((1.0*(m1eSurfNodal[2]*mElc*upar_e-1.0*m2pareSurfNodal[2]*mElc))/m0eSurfNodal[2]));
  double Tpar_i = fmax(0.0, -((1.0*(m1iSurfNodal[2]*mIon*upar_i-1.0*m2pariSurfNodal[2]*mIon))/m0iSurfNodal[2]));

  double c_s = sqrt((3.0*Tpar_i+Tpar_e)/mIon);
  double Delta_phi_r;
  if (c_s > 0.0)
    Delta_phi_r = 2;
  else
    Delta_phi_r = 0.0;

  phiNodal[7] = fmax(phiNodal[7] + Delta_phi_r, phi_wallNodal[2];

  phi[0] = -(0.16666666666666666*phiNodal[7])+0.6666666666666666*phiNodal[6]-0.16666666666666666*phiNodal[5]+0.6666666666666666*phiNodal[4]+0.6666666666666666*phiNodal[3]-0.16666666666666666*phiNodal[2]+0.6666666666666666*phiNodal[1]-0.16666666666666666*phiNodal[0]; 
  phi[1] = 0.09622504486493764*phiNodal[7]-0.09622504486493764*phiNodal[5]+0.3849001794597506*phiNodal[4]-0.3849001794597506*phiNodal[3]+0.09622504486493764*phiNodal[2]-0.09622504486493764*phiNodal[0]; 
  phi[2] = 0.09622504486493764*phiNodal[7]+0.3849001794597506*phiNodal[6]+0.09622504486493764*phiNodal[5]-0.09622504486493764*phiNodal[2]-0.3849001794597506*phiNodal[1]-0.09622504486493764*phiNodal[0]; 
  phi[3] = 0.16666666666666666*phiNodal[7]-0.16666666666666666*phiNodal[5]-0.16666666666666666*phiNodal[2]+0.16666666666666666*phiNodal[0]; 
  phi[4] = 0.14907119849998596*phiNodal[7]-0.2981423969999719*phiNodal[6]+0.14907119849998596*phiNodal[5]+0.14907119849998596*phiNodal[2]-0.2981423969999719*phiNodal[1]+0.14907119849998596*phiNodal[0]; 
  phi[5] = 0.14907119849998596*phiNodal[7]+0.14907119849998596*phiNodal[5]-0.2981423969999719*phiNodal[4]-0.2981423969999719*phiNodal[3]+0.14907119849998596*phiNodal[2]+0.14907119849998596*phiNodal[0]; 
  phi[6] = 0.08606629658238703*phiNodal[7]-0.17213259316477406*phiNodal[6]+0.08606629658238703*phiNodal[5]-0.08606629658238703*phiNodal[2]+0.17213259316477406*phiNodal[1]-0.08606629658238703*phiNodal[0]; 
  phi[7] = 0.08606629658238703*phiNodal[7]-0.08606629658238703*phiNodal[5]-0.17213259316477406*phiNodal[4]+0.17213259316477406*phiNodal[3]+0.08606629658238703*phiNodal[2]-0.08606629658238703*phiNodal[0]; 
 
}
