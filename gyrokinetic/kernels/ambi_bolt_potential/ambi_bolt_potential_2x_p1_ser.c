#include <gkyl_ambi_bolt_potential_kernels.h>

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_lower_2x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  double GammaJacIonB[2];
  GammaJacIonB[0] = 0.6123724356957944*GammaJac_i[2]*sheathDirDx+0.3535533905932737*GammaJac_i[0]*sheathDirDx; 
  GammaJacIonB[1] = 0.6123724356957944*GammaJac_i[3]*sheathDirDx+0.3535533905932737*GammaJac_i[1]*sheathDirDx; 

  double m0JacIonB[2];
  m0JacIonB[0] = 0.7071067811865475*m0JacIon[0]-1.224744871391589*m0JacIon[2]; 
  m0JacIonB[1] = 0.7071067811865475*m0JacIon[1]-1.224744871391589*m0JacIon[3]; 

  // Particle number density evaluate at the sheath entrance
  out[0] = 0.5*((jacInv[3]-1.7320508075688772*jacInv[1])*m0JacIon[3]-1.7320508075688772*m0JacIon[1]*jacInv[3]+(jacInv[2]-1.7320508075688772*jacInv[0])*m0JacIon[2]-1.7320508075688772*m0JacIon[0]*jacInv[2]+jacInv[1]*m0JacIon[1]+jacInv[0]*m0JacIon[0]); 
  out[1] = 0.5*((jacInv[2]-1.7320508075688772*jacInv[0])*m0JacIon[3]+(m0JacIon[2]-1.7320508075688772*m0JacIon[0])*jacInv[3]-1.7320508075688772*jacInv[1]*m0JacIon[2]+m0JacIon[1]*(jacInv[0]-1.7320508075688772*jacInv[2])+m0JacIon[0]*jacInv[1]); 

  double phiS_qp[2];
  if ((isfinite(0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1])) && (0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1]>0.) && (0.7071067811865475*m0JacIonB[0]-0.7071067811865475*m0JacIonB[1]>0.)) {
    phiS_qp[0] = (T_e*log((1.7724538509055159*GammaJacIonB[0])/(0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e)-0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e))-(1.7724538509055159*GammaJacIonB[1])/(0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e)-0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }
  if ((isfinite(0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0])) && (0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0]>0.) && (0.7071067811865475*m0JacIonB[1]+0.7071067811865475*m0JacIonB[0]>0.)) {
    phiS_qp[1] = (T_e*log((1.7724538509055159*GammaJacIonB[1])/(0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e)+0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e))+(1.7724538509055159*GammaJacIonB[0])/(0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e)+0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[1] = 0.0;
  }

  // Sheath potential
  out[4] = phiS_qp[1]+phiS_qp[0]; 
  out[5] = phiS_qp[1]-1.0*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_sheath_calc_upper_2x_ser_p1(double sheathDirDx, double q_e, double m_e, double T_e, const double *jacInv, const double *GammaJac_i, const double *m0JacIon, double *out) 
{ 
  // sheathDirDx: cell length in direction of the sheath.
  // q_e:         electron change.
  // m_e:         electron mass.
  // T_e:         electron temperature.
  // jacInv:      reciprocal of the geometry Jacobian (1/J).
  // GammaJac_i:  ion particle flux (times the Jacobian) through sheath entrance.
  // m0JacIon:    ion density (times the geometry Jacobian).
  // out:         ion density and electrostatic potential at the sheath entrance.

  double GammaJacIonB[2];
  GammaJacIonB[0] = 0.3535533905932737*GammaJac_i[0]*sheathDirDx-0.6123724356957944*GammaJac_i[2]*sheathDirDx; 
  GammaJacIonB[1] = 0.3535533905932737*GammaJac_i[1]*sheathDirDx-0.6123724356957944*GammaJac_i[3]*sheathDirDx; 

  double m0JacIonB[2];
  m0JacIonB[0] = 1.224744871391589*m0JacIon[2]+0.7071067811865475*m0JacIon[0]; 
  m0JacIonB[1] = 1.224744871391589*m0JacIon[3]+0.7071067811865475*m0JacIon[1]; 

  // Particle number density evaluate at the sheath entrance
  out[0] = 0.5*((jacInv[3]+1.7320508075688772*jacInv[1])*m0JacIon[3]+1.7320508075688772*m0JacIon[1]*jacInv[3]+(jacInv[2]+1.7320508075688772*jacInv[0])*m0JacIon[2]+1.7320508075688772*m0JacIon[0]*jacInv[2]+jacInv[1]*m0JacIon[1]+jacInv[0]*m0JacIon[0]); 
  out[1] = 0.5*((jacInv[2]+1.7320508075688772*jacInv[0])*m0JacIon[3]+(m0JacIon[2]+1.7320508075688772*m0JacIon[0])*jacInv[3]+1.7320508075688772*jacInv[1]*m0JacIon[2]+m0JacIon[1]*(1.7320508075688772*jacInv[2]+jacInv[0])+m0JacIon[0]*jacInv[1]); 

  double phiS_qp[2];
  if ((isfinite(0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1])) && (0.7071067811865475*GammaJacIonB[0]-0.7071067811865475*GammaJacIonB[1]>0.) && (0.7071067811865475*m0JacIonB[0]-0.7071067811865475*m0JacIonB[1]>0.)) {
    phiS_qp[0] = (T_e*log((1.7724538509055159*GammaJacIonB[0])/(0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e)-0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e))-(1.7724538509055159*GammaJacIonB[1])/(0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e)-0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[0] = 0.0;
  }
  if ((isfinite(0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0])) && (0.7071067811865475*GammaJacIonB[1]+0.7071067811865475*GammaJacIonB[0]>0.) && (0.7071067811865475*m0JacIonB[1]+0.7071067811865475*m0JacIonB[0]>0.)) {
    phiS_qp[1] = (T_e*log((1.7724538509055159*GammaJacIonB[1])/(0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e)+0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e))+(1.7724538509055159*GammaJacIonB[0])/(0.7071067811865475*m0JacIonB[1]*sqrt(T_e/m_e)+0.7071067811865475*m0JacIonB[0]*sqrt(T_e/m_e))))/q_e;
  } else {
    phiS_qp[1] = 0.0;
  }

  // Sheath potential
  out[4] = phiS_qp[1]+phiS_qp[0]; 
  out[5] = phiS_qp[1]-1.0*phiS_qp[0]; 

}

GKYL_CU_DH void ambi_bolt_potential_phi_calc_2x_ser_p1(double q_e, double T_e, const double *jacInv, const double *m0JacIon, const double *sheathvals, double *phi) 
{ 
  // q_e:        electron change.
  // T_e:        electron temperature.
  // jacInv:     reciprocal of the geometry Jacobian (1/J).
  // m0JacIon:   ion density.
  // sheathvals: ion density and electrostatic potential at the sheath entrance.
  // phi:        electrostatic potential in domain volume.

  double m0Ion[4];
  m0Ion[0] = 0.5*jacInv[3]*m0JacIon[3]+0.5*jacInv[2]*m0JacIon[2]+0.5*jacInv[1]*m0JacIon[1]+0.5*jacInv[0]*m0JacIon[0]; 
  m0Ion[1] = 0.5*jacInv[2]*m0JacIon[3]+0.5*m0JacIon[2]*jacInv[3]+0.5*jacInv[0]*m0JacIon[1]+0.5*m0JacIon[0]*jacInv[1]; 
  m0Ion[2] = 0.5*jacInv[1]*m0JacIon[3]+0.5*m0JacIon[1]*jacInv[3]+0.5*jacInv[0]*m0JacIon[2]+0.5*m0JacIon[0]*jacInv[2]; 
  m0Ion[3] = 0.5*jacInv[0]*m0JacIon[3]+0.5*m0JacIon[0]*jacInv[3]+0.5*jacInv[1]*m0JacIon[2]+0.5*m0JacIon[1]*jacInv[2]; 

  double phi_qp[4];
  phi_qp[0] = -((1.0*log(m0Ion[3]/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[2])/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[3]-1.0*sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.5*sheathvals[7]-0.5*sheathvals[6]-0.5*sheathvals[5]+0.5*sheathvals[4]; 
  phi_qp[1] = -((1.0*log(-((1.0*m0Ion[3])/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))+m0Ion[2]/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])-(1.0*m0Ion[1])/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[3])+sheathvals[2]-1.0*sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.5*sheathvals[7]+0.5*sheathvals[6]-0.5*sheathvals[5]+0.5*sheathvals[4]; 
  phi_qp[2] = -((1.0*log(-((1.0*m0Ion[3])/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))-(1.0*m0Ion[2])/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(-(1.0*sheathvals[3])-1.0*sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)-0.5*sheathvals[7]-0.5*sheathvals[6]+0.5*sheathvals[5]+0.5*sheathvals[4]; 
  phi_qp[3] = -((1.0*log(m0Ion[3]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[2]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[1]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0])+m0Ion[0]/(sheathvals[3]+sheathvals[2]+sheathvals[1]+sheathvals[0]))*T_e)/q_e)+0.5*sheathvals[7]+0.5*sheathvals[6]+0.5*sheathvals[5]+0.5*sheathvals[4]; 

  phi[0] = 0.5*(phi_qp[3]+phi_qp[2]+phi_qp[1]+phi_qp[0]); 
  phi[1] = 0.5*(phi_qp[3]+phi_qp[2])-0.5*(phi_qp[1]+phi_qp[0]); 
  phi[2] = 0.5*phi_qp[3]-0.5*phi_qp[2]+0.5*phi_qp[1]-0.5*phi_qp[0]; 
  phi[3] = 0.5*phi_qp[3]-0.5*(phi_qp[2]+phi_qp[1])+0.5*phi_qp[0]; 
}

