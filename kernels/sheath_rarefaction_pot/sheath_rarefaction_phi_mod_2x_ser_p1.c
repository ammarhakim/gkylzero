#include <gkyl_sheath_rarefaction_pot_kernels.h> 


GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_2x_ser_p1(double elem_q, double mElc, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1, m2par (and m2perp, but not used) moments of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1, m2par (and m2perp, but not used) moments of ions.
  // phi: electrostatic potential.

  double phiNodal[4];
  phiNodal[0] = 1.5*phi[3]-0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[1] = -(1.5*phi[3])-0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[2] = -(1.5*phi[3])+0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[3] = 1.5*phi[3]+0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[4];
  const double *m1Ion = &momsIon[4];
  const double *m2parElc = &momsElc[8];
  const double *m2parIon = &momsIon[8];

  double m0eSurfNodal[2];
  double m0iSurfNodal[2];
  double m1eSurfNodal[2];
  double m1iSurfNodal[2];
  double m2pareSurfNodal[2];
  double m2pariSurfNodal[2];
  double phiWallNodal[2];
  m0eSurfNodal[0] = 1.5*m0Elc[3]-0.8660254037844386*m0Elc[2]-0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0eSurfNodal[1] = -(1.5*m0Elc[3])-0.8660254037844386*m0Elc[2]+0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0iSurfNodal[0] = 1.5*m0Ion[3]-0.8660254037844386*m0Ion[2]-0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m0iSurfNodal[1] = -(1.5*m0Ion[3])-0.8660254037844386*m0Ion[2]+0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m1eSurfNodal[0] = 1.5*m1Elc[3]-0.8660254037844386*m1Elc[2]-0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1eSurfNodal[1] = -(1.5*m1Elc[3])-0.8660254037844386*m1Elc[2]+0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1iSurfNodal[0] = 1.5*m1Ion[3]-0.8660254037844386*m1Ion[2]-0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m1iSurfNodal[1] = -(1.5*m1Ion[3])-0.8660254037844386*m1Ion[2]+0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m2pareSurfNodal[0] = 1.5*m2parElc[3]-0.8660254037844386*m2parElc[2]-0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pareSurfNodal[1] = -(1.5*m2parElc[3])-0.8660254037844386*m2parElc[2]+0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pariSurfNodal[0] = 1.5*m2parIon[3]-0.8660254037844386*m2parIon[2]-0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  m2pariSurfNodal[1] = -(1.5*m2parIon[3])-0.8660254037844386*m2parIon[2]+0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  phiWallNodal[0] = 0.7071067811865475*phiWall[0]-1.224744871391589*phiWall[1]; 
  phiWallNodal[1] = 1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]; 

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

  phi[0] = 0.5*phiNodal[3]+0.5*phiNodal[2]+0.5*phiNodal[1]+0.5*phiNodal[0]; 
  phi[1] = 0.2886751345948129*phiNodal[3]-0.2886751345948129*phiNodal[2]+0.2886751345948129*phiNodal[1]-0.2886751345948129*phiNodal[0]; 
  phi[2] = 0.2886751345948129*phiNodal[3]+0.2886751345948129*phiNodal[2]-0.2886751345948129*phiNodal[1]-0.2886751345948129*phiNodal[0]; 
  phi[3] = 0.16666666666666666*phiNodal[3]-0.16666666666666666*phiNodal[2]-0.16666666666666666*phiNodal[1]+0.16666666666666666*phiNodal[0]; 
 
}

GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_2x_ser_p1(double elem_q, double mElc, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi) 
{ 
  // elem_q: elementary charge.
  // mElc: electron mass.
  // momsElc: m0, m1, m2par (and m2perp, but not used) moments of electrons.
  // mIon: ion mass.
  // momsIon: m0, m1, m2par (and m2perp, but not used) moments of ions.
  // phi: electrostatic potential.

  double phiNodal[4];
  phiNodal[0] = 1.5*phi[3]-0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[1] = -(1.5*phi[3])-0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[2] = -(1.5*phi[3])+0.8660254037844386*phi[2]-0.8660254037844386*phi[1]+0.5*phi[0]; 
  phiNodal[3] = 1.5*phi[3]+0.8660254037844386*phi[2]+0.8660254037844386*phi[1]+0.5*phi[0]; 

  const double *m0Elc = &momsElc[0];
  const double *m0Ion = &momsIon[0];
  const double *m1Elc = &momsElc[4];
  const double *m1Ion = &momsIon[4];
  const double *m2parElc = &momsElc[8];
  const double *m2parIon = &momsIon[8];

  double m0eSurfNodal[2];
  double m0iSurfNodal[2];
  double m1eSurfNodal[2];
  double m1iSurfNodal[2];
  double m2pareSurfNodal[2];
  double m2pariSurfNodal[2];
  double phiWallNodal[2];
  m0eSurfNodal[0] = -(1.5*m0Elc[3])+0.8660254037844386*m0Elc[2]-0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0eSurfNodal[1] = 1.5*m0Elc[3]+0.8660254037844386*m0Elc[2]+0.8660254037844386*m0Elc[1]+0.5*m0Elc[0]; 
  m0iSurfNodal[0] = -(1.5*m0Ion[3])+0.8660254037844386*m0Ion[2]-0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m0iSurfNodal[1] = 1.5*m0Ion[3]+0.8660254037844386*m0Ion[2]+0.8660254037844386*m0Ion[1]+0.5*m0Ion[0]; 
  m1eSurfNodal[0] = -(1.5*m1Elc[3])+0.8660254037844386*m1Elc[2]-0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1eSurfNodal[1] = 1.5*m1Elc[3]+0.8660254037844386*m1Elc[2]+0.8660254037844386*m1Elc[1]+0.5*m1Elc[0]; 
  m1iSurfNodal[0] = -(1.5*m1Ion[3])+0.8660254037844386*m1Ion[2]-0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m1iSurfNodal[1] = 1.5*m1Ion[3]+0.8660254037844386*m1Ion[2]+0.8660254037844386*m1Ion[1]+0.5*m1Ion[0]; 
  m2pareSurfNodal[0] = -(1.5*m2parElc[3])+0.8660254037844386*m2parElc[2]-0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pareSurfNodal[1] = 1.5*m2parElc[3]+0.8660254037844386*m2parElc[2]+0.8660254037844386*m2parElc[1]+0.5*m2parElc[0]; 
  m2pariSurfNodal[0] = -(1.5*m2parIon[3])+0.8660254037844386*m2parIon[2]-0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  m2pariSurfNodal[1] = 1.5*m2parIon[3]+0.8660254037844386*m2parIon[2]+0.8660254037844386*m2parIon[1]+0.5*m2parIon[0]; 
  phiWallNodal[0] = 0.7071067811865475*phiWall[0]-1.224744871391589*phiWall[1]; 
  phiWallNodal[1] = 1.224744871391589*phiWall[1]+0.7071067811865475*phiWall[0]; 

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

  phiNodal[2] = fmax(phiNodal[2] + Delta_phi_r, phi_wallNodal[0];

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

  phiNodal[3] = fmax(phiNodal[3] + Delta_phi_r, phi_wallNodal[1];

  phi[0] = 0.5*phiNodal[3]+0.5*phiNodal[2]+0.5*phiNodal[1]+0.5*phiNodal[0]; 
  phi[1] = 0.2886751345948129*phiNodal[3]-0.2886751345948129*phiNodal[2]+0.2886751345948129*phiNodal[1]-0.2886751345948129*phiNodal[0]; 
  phi[2] = 0.2886751345948129*phiNodal[3]+0.2886751345948129*phiNodal[2]-0.2886751345948129*phiNodal[1]-0.2886751345948129*phiNodal[0]; 
  phi[3] = 0.16666666666666666*phiNodal[3]-0.16666666666666666*phiNodal[2]-0.16666666666666666*phiNodal[1]+0.16666666666666666*phiNodal[0]; 
 
}
