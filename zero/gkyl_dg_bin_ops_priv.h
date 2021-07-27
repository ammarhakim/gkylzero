// Private header: not for direct use
#pragma once

#include <gkyl_util.h>

// Function pointer type for multiplication
typedef void (*mul_op_t)(const double *f, const double *g, double *fg);

// Serendipity multiplication kernels
GKYL_CU_D
static struct { mul_op_t mul[4] } ser_mul_list[] = {
  { NULL, NULL, NULL, NULL }, // No 0D basis functions
  { binop_mul_1d_ser_p0, binop_mul_1d_ser_p1, binop_mul_1d_ser_p2, binop_mul_1d_ser_p3 },
  { binop_mul_2d_ser_p0, binop_mul_2d_ser_p1, binop_mul_2d_ser_p2, binop_mul_2d_ser_p3 },
  { binop_mul_3d_ser_p0, binop_mul_3d_ser_p1, binop_mul_3d_ser_p2, binop_mul_3d_ser_p3 }
};

static mul_op_t
choose_ser_mul_kern(int dim, int poly_order)
{
  return ser_mul_list[dim].mul[poly_order];
}


