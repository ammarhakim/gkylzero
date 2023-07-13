#include <gkyl_vtsq.h>
#include <gkyl_vtsq_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_alloc.h>

struct gkyl_vtsq*
gkyl_vtsq_new(const struct gkyl_basis *basis, const struct gkyl_range *range,
  int m1comps, int vdim_phys, bool use_gpu)
{
  gkyl_vtsq *up = gkyl_malloc(sizeof(struct gkyl_vtsq));

  // MF 2022/11/19: hardcoded arrays for a max of 3x p2 Ser basis below.
  assert(basis->num_basis <= 20);

  up->m1comps = m1comps;
  up->vdim_phys = vdim_phys;
  up->use_gpu = use_gpu;

  // Preallocate memory for the weak division.
  up->mem = use_gpu? gkyl_dg_bin_op_mem_cu_dev_new(range->volume, basis->num_basis)
                    : gkyl_dg_bin_op_mem_new(range->volume, basis->num_basis);

  return up;
}

void
gkyl_vtsq_advance(struct gkyl_vtsq *up, struct gkyl_basis basis,
  const struct gkyl_array *moms, const struct gkyl_range *range,
  struct gkyl_array *out)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_vtsq_advance_cu(up, basis, moms, range, out);
#endif

  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  // Allocate memory for division.
  struct gkyl_nmat *As = up->mem->As;
  struct gkyl_nmat *xs = up->mem->xs;
  // MF 2022/11/19: Hardcoded to a max number of basis for 3x p2 ser.
  double denom[20] = {0.}, numer[20] = {0.};

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, loc);

    // Compute the numerator and denominator in:
    //   (m0*m2 - m1 . m1)/(d_v*m0*m0)

    const double *m0_d = &moms_d[0];
    const double *m2_d = &moms_d[(up->m1comps+1)*num_basis];
    mul_op(m0_d, m2_d, numer);
    for (size_t d=0; d<up->m1comps; d++) {
       const double *m1i_d = &moms_d[(d+1)*num_basis];
       mul_op(m1i_d, m1i_d, denom);  // reuse denom to save memory.
       for (int k=0; k<num_basis; k++) numer[k] -= denom[k];
    }

    mul_op(m0_d, m0_d, denom);
    for (int k=0; k<num_basis; k++) denom[k] *= up->vdim_phys;

    // Place numerator and denominator in mem allocated for weak division.
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
gkyl_vtsq_release(struct gkyl_vtsq *up)
{
  gkyl_dg_bin_op_mem_release(up->mem);
  gkyl_free(up);
}
