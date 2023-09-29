#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_1x1v_ser_p1(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_1x1v_ser_p2(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_1x2v_ser_p1(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_1x2v_ser_p2(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_2x2v_ser_p1(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_2x2v_ser_p2(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_3x2v_ser_p1(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

GKYL_CU_DH void gyrokinetic_mom_cross_bgk_3x2v_ser_p2(const double beta, const double m_self, const double *moms_self, const double m_other, const double *moms_other, const double nu_sr, const double nu_rs, double *moms_cross); 

EXTERN_C_END 
