#include <gkyl_alloc.h>
//#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_calc_priv.h>
#include <gkyl_prim_lbo_kernels.h> 
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_with_fluid.h>
#include <gkyl_mat.h>
#include <assert.h>

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim)
{
  gkyl_prim_lbo_calc *up = gkyl_malloc(sizeof(gkyl_prim_lbo_calc));
  up->grid = *grid;
  up->prim = gkyl_prim_lbo_type_acquire(prim);

  up->is_first = true;
  up->As = up->xs = 0;
  up->mem = 0;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void
gkyl_prim_lbo_calc_advance(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array *uout, struct gkyl_array *vtSqout)
{
  struct gkyl_range_iter conf_iter;

  // allocate memory for use in kernels
  int nc = cbasis.num_basis;
  int udim = calc->prim->udim;
  int N = nc*(udim + 1);

  if (calc->is_first) {
    calc->As = gkyl_nmat_new(conf_rng.volume, N, N);
    calc->xs = gkyl_nmat_new(conf_rng.volume, N, 1);
    calc->mem = gkyl_nmat_linsolve_lu_new(calc->As->num, calc->As->nr);
    calc->is_first = false;
  }

  gkyl_array_clear_range(uout, 0.0, conf_rng);
  gkyl_array_clear_range(vtSqout, 0.0, conf_rng);

  // loop over configuration space cells.
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  long count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);

    struct gkyl_mat lhs = gkyl_nmat_get(calc->As, count);
    struct gkyl_mat rhs = gkyl_nmat_get(calc->xs, count);
    gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

    calc->prim->self_prim(calc->prim, &lhs, &rhs, conf_iter.idx,
      gkyl_array_cfetch(moms, midx), gkyl_array_cfetch(boundary_corrections, midx)
    );

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);
    
    struct gkyl_mat out = gkyl_nmat_get(calc->xs, count);
    double *u = gkyl_array_fetch(uout, midx);
    double *vtSq = gkyl_array_fetch(vtSqout, midx);
    prim_lbo_copy_sol(&out, nc, udim, u, vtSq);
    count += 1;
  }
}

const struct gkyl_prim_lbo_type* gkyl_prim_lbo_calc_get_prim(gkyl_prim_lbo_calc* calc)
{
  return calc->prim;
}

void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* up)
{
  gkyl_prim_lbo_type_release(up->prim);

  if (up->As)
    gkyl_nmat_release(up->As);
  if (up->xs)
    gkyl_nmat_release(up->xs);
  if (up->mem)
    gkyl_nmat_linsolve_lu_release(up->mem);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

// "derived" class constructors
gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_vlasov_new(cbasis, pbasis);
  return gkyl_prim_lbo_calc_new(grid, prim);
}

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_with_fluid_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_vlasov_with_fluid_new(cbasis, pbasis, conf_rng);
  return gkyl_prim_lbo_calc_new(grid, prim);
}

gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_gyrokinetic_new(cbasis, pbasis);
  return gkyl_prim_lbo_calc_new(grid, prim);
}

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_vlasov_cu_dev_new(cbasis, pbasis);
  return gkyl_prim_lbo_calc_cu_dev_new(grid, prim);
}

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_with_fluid_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_vlasov_with_fluid_cu_dev_new(cbasis, pbasis, conf_rng);
  return gkyl_prim_lbo_calc_cu_dev_new(grid, prim);
}

gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis)
{
  struct gkyl_prim_lbo_type *prim; // LBO primitive moments type
  prim = gkyl_prim_lbo_gyrokinetic_cu_dev_new(cbasis, pbasis);
  return gkyl_prim_lbo_calc_cu_dev_new(grid, prim);
}

#ifndef GKYL_HAVE_CUDA

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim)
{
  assert(false);
  return 0;
}

void
gkyl_prim_lbo_calc_advance_cu(gkyl_prim_lbo_calc* calc, const struct gkyl_basis cbasis,
  struct gkyl_range conf_rng,
  const struct gkyl_array* moms, const struct gkyl_array* boundary_corrections,
  struct gkyl_array* uout, struct gkyl_array* vtSqout)
{
  assert(false);
}

#endif
