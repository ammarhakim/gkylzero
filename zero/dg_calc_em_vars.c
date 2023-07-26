#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_em_vars*
gkyl_dg_calc_em_vars_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *conf_range, 
  bool is_ExB, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_em_vars_cu_dev_new(cbasis, conf_range, is_ExB);
  } 
#endif     
  gkyl_dg_calc_em_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_em_vars));
  int nc = cbasis->num_basis;
  enum gkyl_basis_type b_type = cbasis->b_type;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->poly_order = poly_order;
  if (is_ExB) {
    up->Ncomp = 3;
    up->em_calc_temp = choose_em_calc_num_ExB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_ExB_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_ExB_kern(b_type, cdim, poly_order);
  }
  else {
    up->Ncomp = 6;
    up->em_calc_temp = choose_em_calc_BB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_bvar_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_bvar_kern(b_type, cdim, poly_order);    
  }
  // There are Ncomp more linear systems to be solved 
  // 6 components of bb and 3 components of E x B
  up->As = gkyl_nmat_new(up->Ncomp*conf_range->volume, nc, nc);
  up->xs = gkyl_nmat_new(up->Ncomp*conf_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);
  // 6 component temporary variable for either storing B_i B_j (for computing bb) 
  // or (E x B)_i and B_i^2 (for computing E x B/|B|^2)
  up->temp_var = gkyl_array_new(GKYL_DOUBLE, 6*nc, conf_range->volume);
  up->conf_range = *conf_range;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_em_vars_advance(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_calc_em_vars_advance_cu(up, em, cell_avg_magB2, out);
  }
#endif
  gkyl_array_clear(up->temp_var, 0.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->conf_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->conf_range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);
    int* cell_avg_magB2_d = gkyl_array_fetch(cell_avg_magB2, loc);

    up->em_calc_temp(em_d, gkyl_array_fetch(up->temp_var, loc));
    cell_avg_magB2_d[0] = up->em_set(count, up->As, up->xs, gkyl_array_fetch(up->temp_var, loc));

    count += up->Ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->conf_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->conf_range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = gkyl_array_fetch(cell_avg_magB2, loc);
    double *out_d = gkyl_array_fetch(out, loc);

    up->em_copy(count, up->xs, em_d, cell_avg_magB2_d, out_d);

    count += up->Ncomp;
  }  
}

void gkyl_dg_calc_em_vars_release(gkyl_dg_calc_em_vars *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);
  gkyl_array_release(up->temp_var);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
