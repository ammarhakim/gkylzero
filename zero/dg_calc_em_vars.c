#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>

// Methods for computing b_hat and ExB from the basis_inv method, which stores the inverse operation
void gkyl_calc_em_vars_bvar_basis_inv(struct gkyl_basis basis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* bvar)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(bvar)) {
    return gkyl_calc_em_vars_bvar_basis_inv_cu(basis, range, em, bvar);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  em_basis_inv_t em_bvar_basis_inv = choose_ser_em_bvar_basis_inv_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);

    double *bvar_d = gkyl_array_fetch(bvar, loc);
    em_bvar_basis_inv(em_d, bvar_d);
  }
}

void gkyl_calc_em_vars_ExB_basis_inv(struct gkyl_basis basis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* ExB)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(ExB)) {
    return gkyl_calc_em_vars_ExB_basis_inv_cu(basis, range, em, ExB);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  em_basis_inv_t em_ExB_basis_inv = choose_ser_em_ExB_basis_inv_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);

    double *ExB_d = gkyl_array_fetch(ExB, loc);
    em_ExB_basis_inv(em_d, ExB_d);
  }
}

// Methods for computing b_hat and ExB from linsolve methods setting a matrix 
// to compute 1/|B|^2 and then performing the necessary weak multiplications
void gkyl_calc_em_vars_bvar(gkyl_dg_bin_op_mem *mem, 
  struct gkyl_basis basis, const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* bvar)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(bvar)) {
    return gkyl_calc_em_vars_bvar_cu(mem, basis, range, em, bvar);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  // fetch matrix memory for use in kernels
  struct gkyl_nmat *As = mem->As;
  struct gkyl_nmat *xs = mem->xs;

  em_set_t em_set_magB2;
  em_copy_t em_copy_bvar;
  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      em_set_magB2 = choose_ser_em_set_magB2_kern(cdim, poly_order);
      em_copy_bvar = choose_ser_em_copy_bvar_kern(cdim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      em_set_magB2 = choose_ten_em_set_magB2_kern(cdim, poly_order);
      em_copy_bvar = choose_ten_em_copy_bvar_kern(cdim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *em_d = gkyl_array_cfetch(em, loc);
    int* cell_avg_magB2_d = gkyl_array_fetch(cell_avg_magB2, loc);

    struct gkyl_mat A = gkyl_nmat_get(As, count);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0); 

    cell_avg_magB2_d[0] = em_set_magB2(&A, &x, em_d);

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(mem->lu_mem, As, xs);
  assert(status);

  gkyl_range_iter_init(&iter, range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = gkyl_array_cfetch(cell_avg_magB2, loc);
    double *bvar_d = gkyl_array_fetch(bvar, loc);

    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    em_copy_bvar(&x, em_d, cell_avg_magB2_d, bvar_d);

    count += 1;
  }  
}

void gkyl_calc_em_vars_ExB(gkyl_dg_bin_op_mem *mem, 
  struct gkyl_basis basis, const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* ExB)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(ExB)) {
    return gkyl_calc_em_vars_ExB_cu(mem, basis, range, em, ExB);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  // fetch matrix memory for use in kernels
  struct gkyl_nmat *As = mem->As;
  struct gkyl_nmat *xs = mem->xs;

  em_set_t em_set_magB2;
  em_copy_t em_copy_ExB;
  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      em_set_magB2 = choose_ser_em_set_magB2_kern(cdim, poly_order);
      em_copy_ExB = choose_ser_em_copy_ExB_kern(cdim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      em_set_magB2 = choose_ten_em_set_magB2_kern(cdim, poly_order);
      em_copy_ExB = choose_ten_em_copy_ExB_kern(cdim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *em_d = gkyl_array_cfetch(em, loc);
    int* cell_avg_magB2_d = gkyl_array_fetch(cell_avg_magB2, loc);

    struct gkyl_mat A = gkyl_nmat_get(As, count);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0); 

    cell_avg_magB2_d[0] = em_set_magB2(&A, &x, em_d);

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(mem->lu_mem, As, xs);
  assert(status);

  gkyl_range_iter_init(&iter, range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = gkyl_array_cfetch(cell_avg_magB2, loc);
    double *ExB_d = gkyl_array_fetch(ExB, loc);

    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    em_copy_ExB(&x, em_d, cell_avg_magB2_d, ExB_d);

    count += 1;
  }  
}
