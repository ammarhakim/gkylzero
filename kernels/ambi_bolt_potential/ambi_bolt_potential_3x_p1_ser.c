#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_3x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // cmag:        Clebsch function in definition of magnetic field.
  // jacobtotInv: reciprocal of the phase-space and conf-space Jacobians (1/(J*B)).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  // Particle number density evaluate at the sheath entrance
  out[0] = 0.3535533905932737*((jacInv[5]-1.7320508075688772*jacInv[1])*m0JacIon[5]-1.7320508075688772*m0JacIon[1]*jacInv[5]+(jacInv[3]-1.7320508075688772*jacInv[0])*m0JacIon[3]-1.7320508075688772*m0JacIon[0]*jacInv[3]+jacInv[1]*m0JacIon[1]+jacInv[0]*m0JacIon[0]); 
  out[1] = 0.3535533905932737*((jacInv[3]-1.7320508075688772*jacInv[0])*m0JacIon[5]+(m0JacIon[3]-1.7320508075688772*m0JacIon[0])*jacInv[5]-1.7320508075688772*jacInv[1]*m0JacIon[3]+m0JacIon[1]*(jacInv[0]-1.7320508075688772*jacInv[3])+m0JacIon[0]*jacInv[1]); 
  out[2] = 0.3535533905932737*((jacInv[5]-1.7320508075688772*jacInv[1])*m0JacIon[7]+(jacInv[3]-1.7320508075688772*jacInv[0])*m0JacIon[6]+m0JacIon[4]*(jacInv[1]-1.7320508075688772*jacInv[5])+m0JacIon[2]*(jacInv[0]-1.7320508075688772*jacInv[3])); 
  out[4] = 0.3535533905932737*((jacInv[3]-1.7320508075688772*jacInv[0])*m0JacIon[7]+(jacInv[5]-1.7320508075688772*jacInv[1])*m0JacIon[6]-1.7320508075688772*m0JacIon[2]*jacInv[5]+(jacInv[0]-1.7320508075688772*jacInv[3])*m0JacIon[4]+jacInv[1]*m0JacIon[2]); 

  double GammaJacIonB[4];
  GammaJacIonB[0] = GammaJacIonB[0]; 
  GammaJacIonB[1] = GammaJacIonB[1]; 
  GammaJacIonB[2] = GammaJacIonB[2]; 
  GammaJacIonB[3] = GammaJacIonB[3]; 

  double x3HalfMomJacElcB[4];
  x3HalfMomJacElcB[0] = x3HalfMomJacElcB[0]; 
  x3HalfMomJacElcB[1] = x3HalfMomJacElcB[1]; 
  x3HalfMomJacElcB[2] = x3HalfMomJacElcB[2]; 
  x3HalfMomJacElcB[3] = x3HalfMomJacElcB[3]; 

  double m0JacIonB[4];
  m0JacIonB[0] = m0JacIonB[0]; 
  m0JacIonB[1] = m0JacIonB[1]; 
  m0JacIonB[2] = m0JacIonB[2]; 
  m0JacIonB[3] = m0JacIonB[3]; 

  double phiS_qp[4];
  if ((isfinite(0.5*GammaJacIonB[3]-0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (0.5*GammaJacIonB[3]-0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (0.5*m0JacIonB[3]-0.5*m0JacIonB[2]-0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[3]/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])-(1.0*GammaJacIonB[2])/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])-(1.0*GammaJacIonB[1])/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }
  if ((isfinite(-(0.5*GammaJacIonB[3])+0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (-(0.5*GammaJacIonB[3])+0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (-(0.5*m0JacIonB[3])+0.5*m0JacIonB[2]-0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[1] = (log(-((1.0*GammaJacIonB[3])/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))+GammaJacIonB[2]/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])-(1.0*GammaJacIonB[1])/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[1] = 0.0;
  }
  if ((isfinite(-(0.5*GammaJacIonB[3])-0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (-(0.5*GammaJacIonB[3])-0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (-(0.5*m0JacIonB[3])-0.5*m0JacIonB[2]+0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[2] = (log(-((1.0*GammaJacIonB[3])/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))-(1.0*GammaJacIonB[2])/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[1]/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[2] = 0.0;
  }
  if ((isfinite(0.5*GammaJacIonB[3]+0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (0.5*GammaJacIonB[3]+0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (0.5*m0JacIonB[3]+0.5*m0JacIonB[2]+0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[3] = (log(GammaJacIonB[3]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[2]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[1]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[3] = 0.0;
  }

  // Sheath potential
  out[8] = 0.7071067811865475*(phiS_qp[3]+phiS_qp[2]+phiS_qp[1]+phiS_qp[0]); 
  out[9] = 0.7071067811865475*(phiS_qp[3]+phiS_qp[2]-1.0*(phiS_qp[1]+phiS_qp[0])); 
  out[10] = 0.7071067811865475*(phiS_qp[3]-1.0*phiS_qp[2]+phiS_qp[1]-1.0*phiS_qp[0]); 
  out[12] = 0.7071067811865475*(phiS_qp[3]-1.0*(phiS_qp[2]+phiS_qp[1])+phiS_qp[0]); 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_3x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // cmag:        Clebsch function in definition of magnetic field.
  // jacobtotInv: reciprocal of the phase-space and conf-space Jacobians (1/(J*B)).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  // Particle number density evaluate at the sheath entrance
  out[0] = 0.3535533905932737*((jacInv[5]+1.7320508075688772*jacInv[1])*m0JacIon[5]+1.7320508075688772*m0JacIon[1]*jacInv[5]+(jacInv[3]+1.7320508075688772*jacInv[0])*m0JacIon[3]+1.7320508075688772*m0JacIon[0]*jacInv[3]+jacInv[1]*m0JacIon[1]+jacInv[0]*m0JacIon[0]); 
  out[1] = 0.3535533905932737*((jacInv[3]+1.7320508075688772*jacInv[0])*m0JacIon[5]+(m0JacIon[3]+1.7320508075688772*m0JacIon[0])*jacInv[5]+1.7320508075688772*jacInv[1]*m0JacIon[3]+m0JacIon[1]*(1.7320508075688772*jacInv[3]+jacInv[0])+m0JacIon[0]*jacInv[1]); 
  out[2] = 0.3535533905932737*((jacInv[5]+1.7320508075688772*jacInv[1])*m0JacIon[7]+(jacInv[3]+1.7320508075688772*jacInv[0])*m0JacIon[6]+m0JacIon[4]*(1.7320508075688772*jacInv[5]+jacInv[1])+m0JacIon[2]*(1.7320508075688772*jacInv[3]+jacInv[0])); 
  out[4] = 0.3535533905932737*((jacInv[3]+1.7320508075688772*jacInv[0])*m0JacIon[7]+(jacInv[5]+1.7320508075688772*jacInv[1])*m0JacIon[6]+1.7320508075688772*m0JacIon[2]*jacInv[5]+(1.7320508075688772*jacInv[3]+jacInv[0])*m0JacIon[4]+jacInv[1]*m0JacIon[2]); 

  double GammaJacIonB[4];
  GammaJacIonB[0] = GammaJacIonB[0]; 
  GammaJacIonB[1] = GammaJacIonB[1]; 
  GammaJacIonB[2] = GammaJacIonB[2]; 
  GammaJacIonB[3] = GammaJacIonB[3]; 

  double x3HalfMomJacElcB[4];
  x3HalfMomJacElcB[0] = x3HalfMomJacElcB[0]; 
  x3HalfMomJacElcB[1] = x3HalfMomJacElcB[1]; 
  x3HalfMomJacElcB[2] = x3HalfMomJacElcB[2]; 
  x3HalfMomJacElcB[3] = x3HalfMomJacElcB[3]; 

  double m0JacIonB[4];
  m0JacIonB[0] = m0JacIonB[0]; 
  m0JacIonB[1] = m0JacIonB[1]; 
  m0JacIonB[2] = m0JacIonB[2]; 
  m0JacIonB[3] = m0JacIonB[3]; 

  double phiS_qp[4];
  if ((isfinite(0.5*GammaJacIonB[3]-0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (0.5*GammaJacIonB[3]-0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (0.5*m0JacIonB[3]-0.5*m0JacIonB[2]-0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[3]/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])-(1.0*GammaJacIonB[2])/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])-(1.0*GammaJacIonB[1])/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(x3HalfMomJacElcB[3]-1.0*x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }
  if ((isfinite(-(0.5*GammaJacIonB[3])+0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (-(0.5*GammaJacIonB[3])+0.5*GammaJacIonB[2]-0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (-(0.5*m0JacIonB[3])+0.5*m0JacIonB[2]-0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[1] = (log(-((1.0*GammaJacIonB[3])/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))+GammaJacIonB[2]/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])-(1.0*GammaJacIonB[1])/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(-(1.0*x3HalfMomJacElcB[3])+x3HalfMomJacElcB[2]-1.0*x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[1] = 0.0;
  }
  if ((isfinite(-(0.5*GammaJacIonB[3])-0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (-(0.5*GammaJacIonB[3])-0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (-(0.5*m0JacIonB[3])-0.5*m0JacIonB[2]+0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[2] = (log(-((1.0*GammaJacIonB[3])/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))-(1.0*GammaJacIonB[2])/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[1]/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(-(1.0*x3HalfMomJacElcB[3])-1.0*x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[2] = 0.0;
  }
  if ((isfinite(0.5*GammaJacIonB[3]+0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0])) && (0.5*GammaJacIonB[3]+0.5*GammaJacIonB[2]+0.5*GammaJacIonB[1]+0.5*GammaJacIonB[0]>0.) && (0.5*m0JacIonB[3]+0.5*m0JacIonB[2]+0.5*m0JacIonB[1]+0.5*m0JacIonB[0]>0.)) {
    phiS_qp[3] = (log(GammaJacIonB[3]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[2]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[1]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0])+GammaJacIonB[0]/(x3HalfMomJacElcB[3]+x3HalfMomJacElcB[2]+x3HalfMomJacElcB[1]+x3HalfMomJacElcB[0]))*T_e)/q_e;
  } else {
    phiS_qp[3] = 0.0;
  }

  // Sheath potential
  out[8] = 0.7071067811865475*(phiS_qp[3]+phiS_qp[2]+phiS_qp[1]+phiS_qp[0]); 
  out[9] = 0.7071067811865475*(phiS_qp[3]+phiS_qp[2]-1.0*(phiS_qp[1]+phiS_qp[0])); 
  out[10] = 0.7071067811865475*(phiS_qp[3]-1.0*phiS_qp[2]+phiS_qp[1]-1.0*phiS_qp[0]); 
  out[12] = 0.7071067811865475*(phiS_qp[3]-1.0*(phiS_qp[2]+phiS_qp[1])+phiS_qp[0]); 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_3x_ser_p1(double q_e, double T_e, const double *jacInv, const double *m0JacIon, const double *sheathvals, double *phi) 
{ 
  // q_e:        electron change.
  // T_e:        electron temperature.
  // jacInv:     reciprocal of the geometry Jacobian (1/J).
  // m0JacIon:   ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi:        electrostatic potential in domain volume.

  double m0Ion[8];
  m0Ion[0] = 0.3535533905932737*jacInv[5]*m0JacIon[5]+0.3535533905932737*jacInv[3]*m0JacIon[3]+0.3535533905932737*jacInv[1]*m0JacIon[1]+0.3535533905932737*jacInv[0]*m0JacIon[0]; 
  m0Ion[1] = 0.3535533905932737*jacInv[3]*m0JacIon[5]+0.3535533905932737*m0JacIon[3]*jacInv[5]+0.3535533905932737*jacInv[0]*m0JacIon[1]+0.3535533905932737*m0JacIon[0]*jacInv[1]; 
  m0Ion[2] = 0.3535533905932737*jacInv[5]*m0JacIon[7]+0.3535533905932737*jacInv[3]*m0JacIon[6]+0.3535533905932737*jacInv[1]*m0JacIon[4]+0.3535533905932737*jacInv[0]*m0JacIon[2]; 
  m0Ion[3] = 0.3535533905932737*jacInv[1]*m0JacIon[5]+0.3535533905932737*m0JacIon[1]*jacInv[5]+0.3535533905932737*jacInv[0]*m0JacIon[3]+0.3535533905932737*m0JacIon[0]*jacInv[3]; 
  m0Ion[4] = 0.3535533905932737*jacInv[3]*m0JacIon[7]+0.3535533905932737*jacInv[5]*m0JacIon[6]+0.3535533905932737*jacInv[0]*m0JacIon[4]+0.3535533905932737*jacInv[1]*m0JacIon[2]; 
  m0Ion[5] = 0.3535533905932737*jacInv[0]*m0JacIon[5]+0.3535533905932737*m0JacIon[0]*jacInv[5]+0.3535533905932737*jacInv[1]*m0JacIon[3]+0.3535533905932737*m0JacIon[1]*jacInv[3]; 
  m0Ion[6] = 0.3535533905932737*jacInv[1]*m0JacIon[7]+0.3535533905932737*jacInv[0]*m0JacIon[6]+0.3535533905932737*m0JacIon[4]*jacInv[5]+0.3535533905932737*m0JacIon[2]*jacInv[3]; 
  m0Ion[7] = 0.3535533905932737*jacInv[0]*m0JacIon[7]+0.3535533905932737*jacInv[1]*m0JacIon[6]+0.3535533905932737*m0JacIon[2]*jacInv[5]+0.3535533905932737*jacInv[3]*m0JacIon[4]; 

  double phi_qp[8];
  phi_qp[0] = -((1.0*log(-((1.0*m0Ion[7])/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))+m0Ion[6]/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[5]/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[4]/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[3])/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[2])/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[7])+sheathvals[6]+sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.3535533905932737*sheathvals[15]+0.3535533905932737*sheathvals[14]+0.3535533905932737*sheathvals[13]+0.3535533905932737*sheathvals[12]-0.3535533905932737*sheathvals[11]-0.3535533905932737*sheathvals[10]-0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[1] = -((1.0*log(m0Ion[7]/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[6])/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[5])/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[4]/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[3]/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[2])/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[7]-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.3535533905932737*sheathvals[15]-0.3535533905932737*sheathvals[14]-0.3535533905932737*sheathvals[13]+0.3535533905932737*sheathvals[12]+0.3535533905932737*sheathvals[11]-0.3535533905932737*sheathvals[10]-0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[2] = -((1.0*log(m0Ion[7]/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[6])/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[5]/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[4])/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[3])/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[2]/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[7]-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.3535533905932737*sheathvals[15]-0.3535533905932737*sheathvals[14]+0.3535533905932737*sheathvals[13]-0.3535533905932737*sheathvals[12]-0.3535533905932737*sheathvals[11]+0.3535533905932737*sheathvals[10]-0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[3] = -((1.0*log(-((1.0*m0Ion[7])/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))+m0Ion[6]/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[5])/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[4])/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[3]/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[2]/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[7])+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.3535533905932737*sheathvals[15]+0.3535533905932737*sheathvals[14]-0.3535533905932737*sheathvals[13]-0.3535533905932737*sheathvals[12]+0.3535533905932737*sheathvals[11]+0.3535533905932737*sheathvals[10]-0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[4] = -((1.0*log(m0Ion[7]/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[6]/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[5])/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[4])/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[3])/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[2])/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[7]+sheathvals[6]-1.0*sheathvals[5]-1.0*sheathvals[4]-1.0*sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.3535533905932737*sheathvals[15]+0.3535533905932737*sheathvals[14]-0.3535533905932737*sheathvals[13]-0.3535533905932737*sheathvals[12]-0.3535533905932737*sheathvals[11]-0.3535533905932737*sheathvals[10]+0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[5] = -((1.0*log(-((1.0*m0Ion[7])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))-(1.0*m0Ion[6])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[5]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[4])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[3]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[2])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]+sheathvals[5]-1.0*sheathvals[4]+sheathvals[3]-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.3535533905932737*sheathvals[15]-0.3535533905932737*sheathvals[14]+0.3535533905932737*sheathvals[13]-0.3535533905932737*sheathvals[12]+0.3535533905932737*sheathvals[11]-0.3535533905932737*sheathvals[10]+0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[6] = -((1.0*log(-((1.0*m0Ion[7])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0]))-(1.0*m0Ion[6])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[5])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[4]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])-(1.0*m0Ion[3])/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[2]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[7])-1.0*sheathvals[6]-1.0*sheathvals[5]+sheathvals[4]-1.0*sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.3535533905932737*sheathvals[15]-0.3535533905932737*sheathvals[14]-0.3535533905932737*sheathvals[13]+0.3535533905932737*sheathvals[12]-0.3535533905932737*sheathvals[11]+0.3535533905932737*sheathvals[10]+0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 
  phi_qp[7] = -((1.0*log(m0Ion[7]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[6]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[5]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[4]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[3]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[2]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[7]+sheathvals[6]+sheathvals[5]+sheathvals[4]+sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.3535533905932737*sheathvals[15]+0.3535533905932737*sheathvals[14]+0.3535533905932737*sheathvals[13]+0.3535533905932737*sheathvals[12]+0.3535533905932737*sheathvals[11]+0.3535533905932737*sheathvals[10]+0.3535533905932737*sheathvals[9]+0.3535533905932737*sheathvals[8]; 

  phi[0] = 0.3535533905932737*(phi_qp[7]+phi_qp[6]+phi_qp[5]+phi_qp[4]+phi_qp[3]+phi_qp[2]+phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.3535533905932737*(phi_qp[7]+phi_qp[6]+phi_qp[5]+phi_qp[4])-0.3535533905932737*(phi_qp[3]+phi_qp[2]+phi_qp[1]+phi_qp[0]); 
  phi[2] = 0.3535533905932737*(phi_qp[7]+phi_qp[6])-0.3535533905932737*(phi_qp[5]+phi_qp[4])+0.3535533905932737*(phi_qp[3]+phi_qp[2])-0.3535533905932737*(phi_qp[1]+phi_qp[0]); 
  phi[3] = 0.3535533905932737*phi_qp[7]-0.3535533905932737*phi_qp[6]+0.3535533905932737*phi_qp[5]-0.3535533905932737*phi_qp[4]+0.3535533905932737*phi_qp[3]-0.3535533905932737*phi_qp[2]+0.3535533905932737*phi_qp[1]-0.3535533905932737*phi_qp[0]; 
  phi[4] = 0.3535533905932737*(phi_qp[7]+phi_qp[6])-0.3535533905932737*(phi_qp[5]+phi_qp[4]+phi_qp[3]+phi_qp[2])+0.3535533905932737*(phi_qp[1]+phi_qp[0]); 
  phi[5] = 0.3535533905932737*phi_qp[7]-0.3535533905932737*phi_qp[6]+0.3535533905932737*phi_qp[5]-0.3535533905932737*(phi_qp[4]+phi_qp[3])+0.3535533905932737*phi_qp[2]-0.3535533905932737*phi_qp[1]+0.3535533905932737*phi_qp[0]; 
  phi[6] = 0.3535533905932737*phi_qp[7]-0.3535533905932737*(phi_qp[6]+phi_qp[5])+0.3535533905932737*(phi_qp[4]+phi_qp[3])-0.3535533905932737*(phi_qp[2]+phi_qp[1])+0.3535533905932737*phi_qp[0]; 
  phi[7] = 0.3535533905932737*phi_qp[7]-0.3535533905932737*(phi_qp[6]+phi_qp[5])+0.3535533905932737*phi_qp[4]-0.3535533905932737*phi_qp[3]+0.3535533905932737*(phi_qp[2]+phi_qp[1])-0.3535533905932737*phi_qp[0]; 
}

