#pragma once 
#include <math.h> 
#include <gkyl_util.h> 

GKYL_CU_DH void vlasov_alpha_gen_geo_3x3v_ser_p1(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_geo);
