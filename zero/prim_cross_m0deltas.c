#include <gkyl_prim_cross_m0deltas.h>
#include <gkyl_prim_cross_m0deltas_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_alloc.h>

gkyl_prim_cross_m0deltas*
gkyl_prim_cross_m0deltas_new(const struct gkyl_basis *basis, const struct gkyl_range *range,
  double betap1, bool use_gpu)
{
  gkyl_prim_cross_m0deltas *up = gkyl_malloc(sizeof(gkyl_prim_cross_m0deltas));

  // MF 2022/11/19: hardcoded arrays for a max of 3x p2 Ser basis below.
  assert(basis->num_basis <= 20);

  up->betap1 = betap1;
  up->use_gpu = use_gpu;

  // Preallocate memory for the weak division.
  up->mem = use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(range->volume, basis->num_basis)
	             : gkyl_dg_bin_op_mem_new(range->volume, basis->num_basis);

  return up;
}

void
gkyl_prim_cross_m0deltas_advance(gkyl_prim_cross_m0deltas *up, struct gkyl_basis basis,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  const struct gkyl_array* prem0s, const struct gkyl_range *range, struct gkyl_array* out)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_prim_cross_m0deltas_advance_cu(up, basis, massself, m0self, nuself,
      massother, m0other, nuother, prem0s, range, out);

#endif

  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  // Allocate memory for use in kernels.
  struct gkyl_nmat *As = up->mem->As;
  struct gkyl_nmat *xs = up->mem->xs;
  // MF 2022/11/19: Hardcoded to a max number of basis for 3x p2 ser.
  double denom[20], numer[20];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *m0self_d = gkyl_array_cfetch(m0self, loc);
    const double *nuself_d = gkyl_array_cfetch(nuself, loc);
    const double *m0other_d = gkyl_array_cfetch(m0other, loc);
    const double *nuother_d = gkyl_array_cfetch(nuother, loc);
    const double *prem0s_d = gkyl_array_cfetch(prem0s, loc);

    // compute the numerator and denominator in:
    //   m0_s*delta_s = m0_s*2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs)

    mul_op(nuself_d, m0self_d, denom);
    mul_op(nuother_d, m0other_d, numer);

    for (int k=0; k<num_basis; k++) {
      denom[k] *= massself;
      numer[k] *= 2.*up->betap1*massother;
    }

    array_acc1(num_basis, denom, 0.5/up->betap1, numer);

    mul_op(prem0s_d, numer, numer);

    struct gkyl_mat A = gkyl_nmat_get(As, count);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);

    div_set_op(&A, &x, numer, denom);

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem->lu_mem, As, xs);
  assert(status);

  gkyl_range_iter_init(&iter, range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    double *out_d = gkyl_array_fetch(out, loc);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    binop_div_copy_sol(&x, out_d);

    count += 1;
  }
}

void
gkyl_prim_cross_m0deltas_release(gkyl_prim_cross_m0deltas* up)
{
  gkyl_dg_bin_op_mem_release(up->mem);
  gkyl_free(up);
}
