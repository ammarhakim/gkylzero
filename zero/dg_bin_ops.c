#include <assert.h>

#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_mat.h>

// multiplication
void
gkyl_dg_mul_op(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  mul_op_t mul_op = choose_ser_mul_kern(basis.ndim, basis.poly_order);
  int num_basis = basis.num_basis;

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);
    double *out_d = gkyl_array_fetch(out, i);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }
}

void gkyl_dg_mul_op_range(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  mul_op_t mul_op = choose_ser_mul_kern(basis.ndim, basis.poly_order);
  int num_basis = basis.num_basis;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *lop_d = gkyl_array_cfetch(lop, loc);
    const double *rop_d = gkyl_array_cfetch(rop, loc);
    double *out_d = gkyl_array_fetch(out, loc);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);    
  }
}

struct gkyl_kern_op_count
gkyl_dg_mul_op_count(struct gkyl_basis basis)
{
  mul_op_count_t mul_op = choose_ser_mul_op_count_kern(basis.ndim, basis.poly_order);
  return mul_op();
}

// division
void
gkyl_dg_div_op(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  div_op_t div_op = choose_ser_div_kern(basis.ndim, basis.poly_order);
  int num_basis = basis.num_basis;

  // allocate memory for use in kernels
  struct gkyl_mat *A = gkyl_mat_new(num_basis, num_basis, 0.0);
  struct gkyl_mat *x = gkyl_mat_new(num_basis, 1, 0.0);

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);
    double *out_d = gkyl_array_fetch(out, i);

    gkyl_mat_clear(A, 0.0); gkyl_mat_clear(x, 0.0);
    div_op(A, x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }

  gkyl_mat_release(A);
  gkyl_mat_release(x);
}

void gkyl_dg_div_op_range(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  div_op_t div_op = choose_ser_div_kern(basis.ndim, basis.poly_order);
  int num_basis = basis.num_basis;

  // allocate memory for use in kernels
  struct gkyl_mat *A = gkyl_mat_new(num_basis, num_basis, 0.0);
  struct gkyl_mat *x = gkyl_mat_new(num_basis, 1, 0.0);  

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *lop_d = gkyl_array_cfetch(lop, loc);
    const double *rop_d = gkyl_array_cfetch(rop, loc);
    double *out_d = gkyl_array_fetch(out, loc);

    gkyl_mat_clear(A, 0.0); gkyl_mat_clear(x, 0.0);    
    div_op(A, x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }

  gkyl_mat_release(A);
  gkyl_mat_release(x);  
}
