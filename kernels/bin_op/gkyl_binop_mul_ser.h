// Tue Jul 27 09:22:21 2021
#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG

GKYL_CU_DH void binop_mul_1d_ser_p0(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p0(void);

GKYL_CU_DH void binop_mul_1d_ser_p1(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p1(void);

GKYL_CU_DH void binop_mul_1d_ser_p2(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p2(void);

GKYL_CU_DH void binop_mul_1d_ser_p3(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p3(void);

GKYL_CU_DH void binop_mul_2d_ser_p0(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_2d_ser_p0(void);

GKYL_CU_DH void binop_mul_2d_ser_p1(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_2d_ser_p1(void);

GKYL_CU_DH void binop_mul_2d_ser_p2(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_2d_ser_p2(void);

GKYL_CU_DH void binop_mul_2d_ser_p3(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_2d_ser_p3(void);

GKYL_CU_DH void binop_mul_3d_ser_p0(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_3d_ser_p0(void);

GKYL_CU_DH void binop_mul_3d_ser_p1(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_3d_ser_p1(void);

GKYL_CU_DH void binop_mul_3d_ser_p2(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_3d_ser_p2(void);

GKYL_CU_DH void binop_mul_3d_ser_p3(const double *f, const double *g, double *fg );
struct gkyl_kern_op_count op_count_binop_mul_3d_ser_p3(void);
EXTERN_C_END
