#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_1x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.6123724356957944*GammaJac_i[1]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 

  double m0JacIonB[1];
  m0JacIonB[0] = 0.7071067811865475*m0JacIon[0]-1.224744871391589*m0JacIon[1]; 

  // Particle number density evaluate at the sheath entrance
  out[0] = 0.5*((1.4142135623730951*jacInv[1]-2.4494897427831783*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(1.4142135623730951*jacInv[0]-2.4494897427831783*jacInv[1])); 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (T_e*log((2.5066282746310007*GammaJacIonB[0])/(m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[2] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_1x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  double GammaJacIonB[1];
  GammaJacIonB[0] = 0.3535533905932737*GammaJac_i[0]*sheathDirDx-0.6123724356957944*GammaJac_i[1]*sheathDirDx; 

  double m0JacIonB[1];
  m0JacIonB[0] = 1.224744871391589*m0JacIon[1]+0.7071067811865475*m0JacIon[0]; 

  // Particle number density evaluate at the sheath entrance
  out[0] = 0.5*((1.4142135623730951*jacInv[1]+2.4494897427831783*jacInv[0])*m0JacIon[1]+m0JacIon[0]*(2.4494897427831783*jacInv[1]+1.4142135623730951*jacInv[0])); 

  double phiS_qp[1];
  if ((isfinite(GammaJacIonB[0]/2)) && (GammaJacIonB[0]/2>0.) && (m0JacIonB[0]/2>0.)) {
    phiS_qp[0] = (T_e*log((2.5066282746310007*GammaJacIonB[0])/(m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }

  // Sheath potential
  out[2] = 1.4142135623730951*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_1x_ser_p1(double q_e, double T_e, const double *jacInv, const double *m0JacIon, const double *sheathvals, double *phi) 
{ 
  // q_e:        electron change.
  // T_e:        electron temperature.
  // jacInv:     reciprocal of the geometry Jacobian (1/J).
  // m0JacIon:   ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi:        electrostatic potential in domain volume.

  double m0Ion[2];
  m0Ion[0] = 0.7071067811865475*jacInv[1]*m0JacIon[1]+0.7071067811865475*jacInv[0]*m0JacIon[0]; 
  m0Ion[1] = 0.7071067811865475*jacInv[0]*m0JacIon[1]+0.7071067811865475*m0JacIon[0]*jacInv[1]; 

  double phi_qp[2];
  phi_qp[0] = -((1.0*log(m0Ion[0]/(sheathvals[0]-1.0*sheathvals[1])-(1.0*m0Ion[1])/(sheathvals[0]-1.0*sheathvals[1]))*T_e)/q_e)-0.7071067811865475*sheathvals[3]+0.7071067811865475*sheathvals[2]; 
  phi_qp[1] = -((1.0*log(m0Ion[1]/(sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.7071067811865475*sheathvals[3]+0.7071067811865475*sheathvals[2]; 

  phi[0] = 0.7071067811865475*(phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.7071067811865475*phi_qp[1]-0.7071067811865475*phi_qp[0]; 
}

