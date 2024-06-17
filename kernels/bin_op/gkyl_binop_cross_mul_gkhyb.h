// Tue Jul 12 14:48:26 2022
#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG

GKYL_CU_DH void binop_cross_mul_1x1v_gkhyb_p1(const double *f, const double *g, double *fg );

GKYL_CU_DH void binop_cross_mul_1x2v_gkhyb_p1(const double *f, const double *g, double *fg );

GKYL_CU_DH void binop_cross_mul_2x2v_gkhyb_p1(const double *f, const double *g, double *fg );

GKYL_CU_DH void binop_cross_mul_3x2v_gkhyb_p1(const double *f, const double *g, double *fg );

GKYL_CU_DH void binop_cross_mul_comp_par_1x1v_gkhyb_p1(const double *f, const double *g, double *fg, int linc2 );

GKYL_CU_DH void binop_cross_mul_comp_par_1x2v_gkhyb_p1(const double *f, const double *g, double *fg, int linc2 );

GKYL_CU_DH void binop_cross_mul_comp_par_2x2v_gkhyb_p1(const double *f, const double *g, double *fg, int linc2 );

GKYL_CU_DH void binop_cross_mul_comp_par_3x2v_gkhyb_p1(const double *f, const double *g, double *fg, int linc2 );

EXTERN_C_END
