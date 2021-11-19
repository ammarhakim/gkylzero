#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_prim_vlasov_calc.h>
#include <gkyl_mat.h>

#include <assert.h>

gkyl_prim_vlasov_calc*
gkyl_prim_vlasov_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_prim_vlasov *prim)
{
  gkyl_prim_vlasov_calc *up = gkyl_malloc(sizeof(gkyl_prim_vlasov_calc));
  up->grid = *grid;
  up->prim = prim;
  return up;
}

void
gkyl_prim_vlasov_calc_advance(const gkyl_prim_vlasov_calc* calc, const struct gkyl_basis cbasis,
  const struct gkyl_range conf_rng, const struct gkyl_array *m0, const struct gkyl_array *m1,
  const struct gkyl_array *m2, const struct gkyl_array *cM, const struct gkyl_array *cE,
  struct gkyl_array *uout, struct gkyl_array *vtSqout)
{
  struct gkyl_range_iter conf_iter;

  // allocate memory for use in kernels
  int N = cbasis.num_basis*(calc->prim->pdim - calc->prim->cdim);
  struct gkyl_mat *A = gkyl_mat_new(N, N, 0.0);
  struct gkyl_mat *rhs = gkyl_mat_new(N, 1, 0.0);

  gkyl_array_clear_range(uout, 0.0, conf_rng);
  gkyl_array_clear_range(vtSqout, 0.0, conf_rng);

  // loop over configuration space cells.
  gkyl_range_iter_init(&conf_iter, &conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(&conf_rng, conf_iter.idx);
    gkyl_mat_clear(A, 0.0); gkyl_mat_clear(rhs, 0.0);

    calc->prim->self_prim(A, rhs, gkyl_array_cfetch(m0, midx), gkyl_array_cfetch(m1, midx),
      gkyl_array_cfetch(m2, midx), gkyl_array_cfetch(cM, midx), gkyl_array_cfetch(cE, midx),
      gkyl_array_fetch(uout, midx), gkyl_array_fetch(vtSqout, midx)
    );
  }

  gkyl_mat_release(A);
  gkyl_mat_release(rhs);
}

void gkyl_prim_vlasov_calc_release(gkyl_prim_vlasov_calc* up)
{
  gkyl_prim_vlasov_release(up->prim);
  gkyl_free(up);
}
