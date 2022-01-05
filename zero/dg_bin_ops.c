#include <assert.h>

#include <gkyl_util.h>
#include <gkyl_array_ops_priv.h>
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
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_dg_mul_op_cu(basis, c_oop, out, c_lop, lop, c_rop, rop); }
#endif

  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

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
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_dg_mul_op_range_cu(basis, c_oop, out, c_lop, lop, c_rop, rop, range); }
#endif

  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

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

// division
void
gkyl_dg_div_op(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_dg_div_op_cu(basis, c_oop, out, c_lop, lop, c_rop, rop); }
#endif

  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);

  // allocate memory for use in kernels
  struct gkyl_nmat *As = gkyl_nmat_new(out->size, num_basis, num_basis);
  struct gkyl_nmat *xs = gkyl_nmat_new(out->size, num_basis, 1);

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);

    struct gkyl_mat A = gkyl_nmat_get(As, i);
    struct gkyl_mat x = gkyl_nmat_get(xs, i);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);
    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }

  bool status = gkyl_nmat_linsolve_lu(As, xs);

  for (size_t i=0; i<out->size; ++i) {
    double *out_d = gkyl_array_fetch(out, i);
    struct gkyl_mat x = gkyl_nmat_get(xs, i);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }

  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);
}

void gkyl_dg_div_op_range(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_dg_div_op_range_cu(basis, c_oop, out, c_lop, lop, c_rop, rop, range); }
#endif

  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);

  // allocate memory for use in kernels
  struct gkyl_nmat *As = gkyl_nmat_new(range.volume, num_basis, num_basis);
  struct gkyl_nmat *xs = gkyl_nmat_new(range.volume, num_basis, 1);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *lop_d = gkyl_array_cfetch(lop, loc);
    const double *rop_d = gkyl_array_cfetch(rop, loc);

    struct gkyl_mat A = gkyl_nmat_get(As, loc);
    struct gkyl_mat x = gkyl_nmat_get(xs, loc);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);  
    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }

  bool status = gkyl_nmat_linsolve_lu(As, xs);

  for (size_t i=0; i<out->size; ++i) {
    double *out_d = gkyl_array_fetch(out, i);
    struct gkyl_mat x = gkyl_nmat_get(xs, i);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }

  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);  
}
