#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_lower_1x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_lower_1x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_upper_1x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_upper_1x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_lower_2x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_lower_2x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_upper_2x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_upper_2x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_lower_3x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_lower_3x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_upper_3x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_upper_3x_ser_p1(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_lower_1x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_lower_1x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_upper_1x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_upper_1x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_lower_2x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_lower_2x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_upper_2x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_upper_2x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_lower_3x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_lower_3x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 

GKYL_CU_DH void sheath_rarefaction_phi_mod_kinetic_elc_upper_3x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_boltzmann_elc_upper_3x_ser_p2(double elem_q, double mElc, double tempElcBoltz, const double *momsElc, double mIon, const double *momsIon, const double *phiWall, double *phi); 


EXTERN_C_END

