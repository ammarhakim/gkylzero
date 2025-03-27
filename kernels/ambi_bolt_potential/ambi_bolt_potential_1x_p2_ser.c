#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
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
  out[0] = 0.031943828249996996*((53.75872022286246*jacInv[2]-34.292856398964496*jacInv[1]+49.49747468305833*jacInv[0])*m0JacIon[2]+(49.49747468305833*m0JacIon[0]-34.292856398964496*m0JacIon[1])*jacInv[2]+(66.40783086353598*jacInv[1]-38.34057902536163*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(22.135943621178658*jacInv[0]-38.34057902536163*jacInv[1])); 

  double GammaJacIonB[1];
  GammaJacIonB[0] = GammaJacIonB[0]; 

  double x3HalfMomJacElcB[1];
  x3HalfMomJacElcB[0] = x3HalfMomJacElcB[0]; 

  double m0JacIonB[1];
  m0JacIonB[0] = m0JacIonB[0]; 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/x3HalfMomJacElcB[0])*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[3] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_ser_p2(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *cmag, const double *jacobtotInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
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
  out[0] = 0.031943828249996996*((53.75872022286246*jacInv[2]+34.292856398964496*jacInv[1]+49.49747468305833*jacInv[0])*m0JacIon[2]+(34.292856398964496*m0JacIon[1]+49.49747468305833*m0JacIon[0])*jacInv[2]+(66.40783086353598*jacInv[1]+38.34057902536163*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(38.34057902536163*jacInv[1]+22.135943621178658*jacInv[0])); 

  double GammaJacIonB[1];
  GammaJacIonB[0] = GammaJacIonB[0]; 

  double x3HalfMomJacElcB[1];
  x3HalfMomJacElcB[0] = x3HalfMomJacElcB[0]; 

  double m0JacIonB[1];
  m0JacIonB[0] = m0JacIonB[0]; 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (log(GammaJacIonB[0]/x3HalfMomJacElcB[0])*T_e)/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[3] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_ser_p2(double q_e, double T_e, const double *jacInv, const double *m0JacIon, const double *sheathvals, double *phi) 
{ 
  // q_e:        electron change.
  // T_e:        electron temperature.
  // jacInv:     reciprocal of the geometry Jacobian (1/J).
  // m0JacIon:   ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi:        electrostatic potential in domain volume.

  double m0Ion[3];
  m0Ion[0] = 0.7071067811865475*jacInv[2]*m0JacIon[2]+0.7071067811865475*jacInv[1]*m0JacIon[1]+0.7071067811865475*jacInv[0]*m0JacIon[0]; 
  m0Ion[1] = 0.6324555320336759*jacInv[1]*m0JacIon[2]+0.6324555320336759*m0JacIon[1]*jacInv[2]+0.7071067811865475*jacInv[0]*m0JacIon[1]+0.7071067811865475*m0JacIon[0]*jacInv[1]; 
  m0Ion[2] = 0.45175395145262565*jacInv[2]*m0JacIon[2]+0.7071067811865475*jacInv[0]*m0JacIon[2]+0.7071067811865475*m0JacIon[0]*jacInv[2]+0.6324555320336759*jacInv[1]*m0JacIon[1]; 

  double phi_qp[3];
  phi_qp[0] = -((1.0*log(-((3.0*m0Ion[1])/(2.0*sheathvals[2]-3.0*sheathvals[1]+2.23606797749979*sheathvals[0]))+(1.4142135623730951*m0Ion[2])/(1.4142135623730951*sheathvals[2]-2.1213203435596424*sheathvals[1]+1.5811388300841895*sheathvals[0])+m0Ion[0]/(0.8944271909999159*sheathvals[2]-1.3416407864998738*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.6324555320336759*sheathvals[5]-0.9486832980505137*sheathvals[4]+0.7071067811865475*sheathvals[3]; 
  phi_qp[1] = -((1.0*log(m0Ion[0]/(sheathvals[0]-1.118033988749895*sheathvals[2])-(2.23606797749979*m0Ion[2])/(2.0*sheathvals[0]-2.23606797749979*sheathvals[2]))*T_e)/q_e)-0.7905694150420947*sheathvals[5]+0.7071067811865475*sheathvals[3]; 
  phi_qp[2] = -((1.0*log((3.0*m0Ion[1])/(2.0*sheathvals[2]+3.0*sheathvals[1]+2.23606797749979*sheathvals[0])+(1.4142135623730951*m0Ion[2])/(1.4142135623730951*sheathvals[2]+2.1213203435596424*sheathvals[1]+1.5811388300841895*sheathvals[0])+m0Ion[0]/(0.8944271909999159*sheathvals[2]+1.3416407864998738*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.6324555320336759*sheathvals[5]+0.9486832980505137*sheathvals[4]+0.7071067811865475*sheathvals[3]; 

  phi[0] = 0.39283710065919303*phi_qp[2]+0.6285393610547091*phi_qp[1]+0.39283710065919303*phi_qp[0]; 
  phi[1] = 0.5270462766947298*phi_qp[2]-0.5270462766947298*phi_qp[0]; 
  phi[2] = 0.35136418446315326*phi_qp[2]-0.7027283689263066*phi_qp[1]+0.35136418446315326*phi_qp[0]; 
}

