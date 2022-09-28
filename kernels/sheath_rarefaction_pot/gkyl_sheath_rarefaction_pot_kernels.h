#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_1x_ser_p1(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_1x_ser_p1(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_1x_ser_p2(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_1x_ser_p2(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_1x_tensor_p1(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_1x_tensor_p1(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_lower_1x_tensor_p2(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 
GKYL_CU_DH void sheath_rarefaction_phi_mod_upper_1x_tensor_p2(double elem_q, double mElc, const double *momsElc, const double *m2parElc, double mIon, const double *momsIon, const double *m2parIon, const double *phiWall, double *phi); 

EXTERN_C_END

