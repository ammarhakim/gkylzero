#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_exp_priv.h>
#include <gkyl_range.h>

struct gkyl_dg_pkpm_harm_exp *
gkyl_dg_pkpm_harm_exp_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_harm_exp_cu_dev_new(inp);
  }
#endif
  assert(inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  struct gkyl_dg_pkpm_harm_exp *up = gkyl_malloc(sizeof(*up));
  up->cdim = inp->conf_basis->ndim;
  up->num_conf_nodes = 1 << up->cdim;
  up->num_phase_nodes = 3 * up->num_conf_nodes;
  up->qbym = inp->qbym;
  up->use_gpu = false;
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

int
gkyl_dg_pkpm_harm_exp_num_coef(const struct gkyl_dg_pkpm_harm_exp *up)
{
  return up->num_conf_nodes * PKPM_HX_NC;
}

int
gkyl_dg_pkpm_harm_exp_num_stage(const struct gkyl_dg_pkpm_harm_exp *up)
{
  return up->num_phase_nodes * 6;
}

void
gkyl_dg_pkpm_harm_exp_coef(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range, double h,
  const struct gkyl_array *em, struct gkyl_array *coef
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(coef)) {
    gkyl_dg_pkpm_harm_exp_coef_cu(up, conf_range, h, em, coef);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    pkpm_hx_coef_cell(up, h, gkyl_array_cfetch(em, loc), gkyl_array_fetch(coef, loc));
  }
}

void
gkyl_dg_pkpm_harm_exp_accum(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, struct gkyl_array *S, struct gkyl_array *A
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(S)) {
    gkyl_dg_pkpm_harm_exp_accum_cu(up, stage, h, conf_range, phase_range, coef, fin, rhs, S, A);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    pkpm_hx_accum_cell(up, stage, h, gkyl_array_cfetch(coef, cloc), gkyl_array_cfetch(fin, ploc),
      gkyl_array_cfetch(rhs, ploc), gkyl_array_fetch(S, ploc), gkyl_array_fetch(A, ploc));
  }
}

void
gkyl_dg_pkpm_harm_exp_state(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  const struct gkyl_array *S, const struct gkyl_array *A, struct gkyl_array *out
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    gkyl_dg_pkpm_harm_exp_state_cu(up, stage, h, conf_range, phase_range, coef, f0, S, A, out);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(conf_range, iter.idx);
    long ploc = gkyl_range_idx(phase_range, iter.idx);
    pkpm_hx_state_cell(up, stage, h, gkyl_array_cfetch(coef, cloc), gkyl_array_cfetch(f0, ploc),
      gkyl_array_cfetch(S, ploc), gkyl_array_cfetch(A, ploc), gkyl_array_fetch(out, ploc));
  }
}

void
gkyl_dg_pkpm_harm_exp_release(struct gkyl_dg_pkpm_harm_exp *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_harm_exp *
gkyl_dg_pkpm_harm_exp_cu_dev_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
  assert(false);
  return 0;
}

#endif
