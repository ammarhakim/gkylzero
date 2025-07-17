#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_gk_neut_fluid_prim_vars.h>
#include <gkyl_gk_neut_fluid_prim_vars_priv.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_euler.h>
#include <gkyl_util.h>

gkyl_gk_neut_fluid_prim_vars*
gkyl_gk_neut_fluid_prim_vars_new(double gas_gamma, const struct gkyl_basis* cbasis,
  const struct gkyl_range *mem_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_gk_neut_fluid_prim_vars_cu_dev_new(gas_gamma, cbasis, mem_range);
  }
#endif
  gkyl_gk_neut_fluid_prim_vars *up = gkyl_malloc(sizeof(gkyl_gk_neut_fluid_prim_vars));

  up->gas_gamma = gas_gamma;

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->num_basis = cbasis->num_basis;
  up->udrift_ncomp = 3;
  up->mem_range = *mem_range;

  up->udrift_set_prob_ker = choose_udrift_set_prob_ker(b_type, cdim, poly_order);
  up->udrift_get_sol_ker = choose_udrift_get_sol_ker(b_type, cdim, poly_order);
  up->pressure_ker = choose_pressure_ker(b_type, cdim, poly_order);
//  up->driftE_ker = choose_driftE_ker(b_type, cdim, poly_order);
//  up->thermalE_ker = choose_thermalE_ker(b_type, cdim, poly_order);

  // There are udrift_ncomp*range->volume linear systems to be solved
  // for 3 components of u: ux, uy, uz.
  up->As = gkyl_nmat_new(up->udrift_ncomp*mem_range->volume, up->num_basis, up->num_basis);
  up->xs = gkyl_nmat_new(up->udrift_ncomp*mem_range->volume, up->num_basis, 1);
  if (up->poly_order > 1) {
    up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void gkyl_gk_neut_fluid_prim_vars_udrift_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *udrift)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(u)) {
    return gkyl_gk_neut_fluid_prim_vars_udrift_advance_cu(up, moms, udrift);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_set_prob_ker(count, up->As, up->xs, moms_d);

    count += up->udrift_ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    double* udrift_d = gkyl_array_fetch(udrift, linidx);

    up->udrift_get_sol_ker(count, up->xs, udrift_d);

    count += up->udrift_ncomp;
  }
}

void gkyl_gk_neut_fluid_prim_vars_pressure_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *pressure)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(u)) {
    return gkyl_gk_neut_fluid_prim_vars_pressure_advance_cu(up, moms, udrift);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_set_prob_ker(count, up->As, up->xs, moms_d);

    count += up->udrift_ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    double udrift_d[up->num_basis];
    for (int i=0; i<up->num_basis; i++)
      udrift_d[i] = 0.0;

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* pressure_d = gkyl_array_fetch(pressure, linidx);

    up->udrift_get_sol_ker(count, up->xs, udrift_d);
    up->pressure_ker(up->gas_gamma, moms_d, udrift_d, pressure_d);

    count += up->udrift_ncomp;
  }
}

void gkyl_gk_neut_fluid_prim_vars_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *prim_vars)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(u)) {
    return gkyl_gk_neut_fluid_prim_vars_advance_cu(up, moms, udrift);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_set_prob_ker(count, up->As, up->xs, moms_d);

    count += up->udrift_ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* prim_vars_d = gkyl_array_fetch(prim_vars, linidx);

    double* udrift_d = prim_vars_d;
    double* pressure_d = &prim_vars_d[up->udrift_ncomp*up->num_basis];

    up->udrift_get_sol_ker(count, up->xs, udrift_d);
    up->pressure_ker(up->gas_gamma, moms_d, udrift_d, pressure_d);

    count += up->udrift_ncomp;
  }
}

void gkyl_gk_neut_fluid_prim_vars_release(gkyl_gk_neut_fluid_prim_vars *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  if (up->poly_order > 1) {
    gkyl_nmat_linsolve_lu_release(up->mem);
  }

  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
