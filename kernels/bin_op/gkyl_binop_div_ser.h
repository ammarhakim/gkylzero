// Mon Jul 19 01:21:18 2021
#pragma once
#include <gkyl_util.h>
#include <gkyl_mat.h>
EXTERN_C_BEG
void binop_div_1d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_1d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_1d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_1d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_2d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_2d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_2d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_2d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_3d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_3d_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_3d_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
void binop_div_3d_ser_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g, double *fdivg );
EXTERN_C_END
