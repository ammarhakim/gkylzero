#include <gkyl_binop_mul_ser.h>
#include <gkyl_dg_bin_ops.h>

#include <assert.h>

typedef void (*mul_op_t)(const double *f, const double *g, double *fg);

// Serendipity multiplication kernels
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

void
gkyl_dg_mul_op(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  mul_op_t mul_op = choose_ser_mul_kern(basis.ndim, basis.polyOrder);
  int num_basis = basis.numBasis;

  assert( (out->size == lop->size) && (out->size == rop->size) );

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);
    double *out_d = gkyl_array_fetch(out, i);

    mul_op(lop_d, rop_d, out_d);
  }
}
