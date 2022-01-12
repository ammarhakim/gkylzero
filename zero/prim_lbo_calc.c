#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
#include <gkyl_mat.h>
#include <assert.h>

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo *prim)
{
  gkyl_prim_lbo_calc *up = gkyl_malloc(sizeof(gkyl_prim_lbo_calc));
  up->grid = *grid;
  up->prim = prim;
  return up;
}

void
gkyl_prim_lbo_calc_advance(gkyl_prim_lbo_calc* calc, const struct gkyl_basis cbasis,
  const struct gkyl_range conf_rng, const struct gkyl_array *m0, const struct gkyl_array *m1,
  const struct gkyl_array *m2, const struct gkyl_array *cM, const struct gkyl_array *cE,
  struct gkyl_array *uout, struct gkyl_array *vtSqout)
{
  struct gkyl_range_iter conf_iter;

  // allocate memory for use in kernels
  int nc = cbasis.num_basis;
  int vdim = calc->prim->pdim - calc->prim->cdim;
  int N = nc*(vdim + 1);
  struct gkyl_nmat *As = gkyl_nmat_new(conf_rng.volume, N, N);
  struct gkyl_nmat *xs = gkyl_nmat_new(conf_rng.volume, N, 1);

  gkyl_array_clear_range(uout, 0.0, conf_rng);
  gkyl_array_clear_range(vtSqout, 0.0, conf_rng);

  // loop over configuration space cells.
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  long count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);

    struct gkyl_mat lhs = gkyl_nmat_get(As, count);
    struct gkyl_mat rhs = gkyl_nmat_get(xs, count);
    gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

    calc->prim->self_prim(calc->prim, &lhs, &rhs, gkyl_array_cfetch(m0, midx), gkyl_array_cfetch(m1, midx),
      gkyl_array_cfetch(m2, midx), gkyl_array_cfetch(cM, midx), gkyl_array_cfetch(cE, midx)
    );

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu(As, xs);
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);
    
    struct gkyl_mat out = gkyl_nmat_get(xs, count);
    double *u = gkyl_array_fetch(uout, midx);
    double *vtSq = gkyl_array_fetch(vtSqout, midx);
    prim_lbo_copy_sol(&out, nc, vdim, u, vtSq);
    count += 1;
  }

  gkyl_nmat_release(As);
  gkyl_nmat_release(xs);
}

void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* up)
{
  gkyl_prim_lbo_release(up->prim);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo *prim)
{
  assert(false);
  return 0;
}

void
gkyl_prim_lbo_calc_advance_cu(gkyl_prim_lbo_calc* mcalc, const struct gkyl_basis cbasis,
  const struct gkyl_range conf_rng, const struct gkyl_array* m0, const struct gkyl_array* m1,
  const struct gkyl_array* m2, const struct gkyl_array* cM, const struct gkyl_array* cE,
  struct gkyl_array* uout, struct gkyl_array* vtSqout)
{
  assert(false);
}

#endif
