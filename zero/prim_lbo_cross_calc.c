#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_cross_calc_priv.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
#include <gkyl_mat.h>
#include <assert.h>

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim, int nspecies)
{
  gkyl_prim_lbo_cross_calc *up = gkyl_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  up->grid = *grid;
  up->prim = gkyl_prim_lbo_type_acquire(prim);
  up->nspecies = nspecies;

  up->is_first = true;
  up->As = up->xs = 0;
  up->mem = 0;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void
gkyl_prim_lbo_cross_calc_advance(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range conf_rng,
  struct gkyl_array *greene[GKYL_MAX_SPECIES], const double self_m,
  const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  const double cross_m[GKYL_MAX_SPECIES], struct gkyl_array *cross_u[GKYL_MAX_SPECIES],
  struct gkyl_array *cross_vtsq[GKYL_MAX_SPECIES], const struct gkyl_array *moms,
  const struct gkyl_array *boundary_corrections, struct gkyl_array *u_out[GKYL_MAX_SPECIES],
  struct gkyl_array *vtsq_out[GKYL_MAX_SPECIES])
{
  struct gkyl_range_iter conf_iter;

  int nspecies = calc->nspecies;
  double *cross_us[GKYL_MAX_SPECIES];
  double *cross_vtsqs[GKYL_MAX_SPECIES];
  double *u_outs[GKYL_MAX_SPECIES];
  double *vtsq_outs[GKYL_MAX_SPECIES];
  double *greenes[GKYL_MAX_SPECIES];
  
  // allocate memory for use in kernels
  int nc = cbasis.num_basis;
  int vdim = calc->prim->pdim - calc->prim->cdim;
  int N = nc*(vdim + 1);

  if (calc->is_first) {
    calc->As = gkyl_nmat_new(nspecies*conf_rng.volume, N, N);
    calc->xs = gkyl_nmat_new(nspecies*conf_rng.volume, N, 1);
    calc->mem = gkyl_nmat_linsolve_lu_new(calc->As->num, calc->As->nr);
    calc->is_first = false;
  }

  // loop over configuration space cells.
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  long count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);
    for (int n=0; n<nspecies; ++n) {
      cross_us[n] = gkyl_array_fetch(cross_u[n], midx);
      cross_vtsqs[n] = gkyl_array_fetch(cross_vtsq[n], midx);
      greenes[n] = gkyl_array_fetch(greene[n], midx);

      struct gkyl_mat lhs = gkyl_nmat_get(calc->As, count);
      struct gkyl_mat rhs = gkyl_nmat_get(calc->xs, count);

      gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

      calc->prim->cross_prim(calc->prim, &lhs, &rhs, conf_iter.idx, greenes[n], self_m,
	gkyl_array_cfetch(self_u, midx), gkyl_array_cfetch(self_vtsq, midx),
        cross_m[n], cross_us[n], cross_vtsqs[n], gkyl_array_cfetch(moms, midx),
        gkyl_array_cfetch(boundary_corrections, midx)
      );

      count += 1;
    }
  }

  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);
    for (int n=0; n<nspecies; ++n) {
      struct gkyl_mat out = gkyl_nmat_get(calc->xs, count);
      u_outs[n] = gkyl_array_fetch(u_out[n], midx);
      vtsq_outs[n] = gkyl_array_fetch(vtsq_out[n], midx);
      prim_lbo_copy_sol(&out, nc, vdim, u_outs[n], vtsq_outs[n]);
      count += 1;
    }
  }
}

void gkyl_prim_lbo_cross_calc_release(gkyl_prim_lbo_cross_calc* up)
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

#ifndef GKYL_HAVE_CUDA

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim, int nspecies)
{
  assert(false);
  return 0;
}

void
gkyl_prim_lbo_cross_calc_advance_cu(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range conf_rng,
  struct gkyl_array *greene[GKYL_MAX_SPECIES], const double self_m,
  const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  const double cross_m[GKYL_MAX_SPECIES], struct gkyl_array *cross_u[GKYL_MAX_SPECIES],
  struct gkyl_array *cross_vtsq[GKYL_MAX_SPECIES], const struct gkyl_array *moms,
  const struct gkyl_array *boundary_corrections, struct gkyl_array *u_out[GKYL_MAX_SPECIES],
  struct gkyl_array *vtsq_out[GKYL_MAX_SPECIES])
{
  assert(false);
}

#endif
