/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_exp_priv.h>
}

// One thread per configuration cell.
__global__ static void
gkyl_dg_pkpm_harm_exp_coef_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, struct gkyl_range conf_range, double h,
  const struct gkyl_array *em, struct gkyl_array *coef
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    pkpm_hx_coef_cell(up, h, (const double *)gkyl_array_cfetch(em, loc), (double *)gkyl_array_fetch(coef, loc));
  }
}

void
gkyl_dg_pkpm_harm_exp_coef_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, const struct gkyl_range *conf_range, double h,
  const struct gkyl_array *em, struct gkyl_array *coef
)
{
  gkyl_dg_pkpm_harm_exp_coef_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, h, em->on_dev, coef->on_dev
  );
}

// One thread per phase cell.
__global__ static void
gkyl_dg_pkpm_harm_exp_accum_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, struct gkyl_range conf_range,
  struct gkyl_range phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, struct gkyl_array *S, struct gkyl_array *A
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < phase_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc, idx);
    long cloc = gkyl_range_idx(&conf_range, idx);
    long ploc = gkyl_range_idx(&phase_range, idx);
    pkpm_hx_accum_cell(up, stage, h, (const double *)gkyl_array_cfetch(coef, cloc),
      (const double *)gkyl_array_cfetch(fin, ploc), (const double *)gkyl_array_cfetch(rhs, ploc),
      (double *)gkyl_array_fetch(S, ploc), (double *)gkyl_array_fetch(A, ploc));
  }
}

void
gkyl_dg_pkpm_harm_exp_accum_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  const struct gkyl_array *rhs, struct gkyl_array *S, struct gkyl_array *A
)
{
  gkyl_dg_pkpm_harm_exp_accum_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, stage, h, *conf_range, *phase_range, coef->on_dev, fin->on_dev, rhs->on_dev,
    S->on_dev, A->on_dev
  );
}

__global__ static void
gkyl_dg_pkpm_harm_exp_state_cu_kernel(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, struct gkyl_range conf_range,
  struct gkyl_range phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  const struct gkyl_array *S, const struct gkyl_array *A, struct gkyl_array *out
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < phase_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc, idx);
    long cloc = gkyl_range_idx(&conf_range, idx);
    long ploc = gkyl_range_idx(&phase_range, idx);
    pkpm_hx_state_cell(up, stage, h, (const double *)gkyl_array_cfetch(coef, cloc),
      (const double *)gkyl_array_cfetch(f0, ploc), (const double *)gkyl_array_cfetch(S, ploc),
      (const double *)gkyl_array_cfetch(A, ploc), (double *)gkyl_array_fetch(out, ploc));
  }
}

void
gkyl_dg_pkpm_harm_exp_state_cu(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *f0,
  const struct gkyl_array *S, const struct gkyl_array *A, struct gkyl_array *out
)
{
  gkyl_dg_pkpm_harm_exp_state_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, stage, h, *conf_range, *phase_range, coef->on_dev, f0->on_dev, S->on_dev,
    A->on_dev, out->on_dev
  );
}

struct gkyl_dg_pkpm_harm_exp *
gkyl_dg_pkpm_harm_exp_cu_dev_new(const struct gkyl_dg_pkpm_harm_exp_inp *inp)
{
  struct gkyl_dg_pkpm_harm_exp *up =
    (struct gkyl_dg_pkpm_harm_exp *)gkyl_malloc(sizeof(struct gkyl_dg_pkpm_harm_exp));
  up->cdim = inp->conf_basis->ndim;
  up->num_conf_nodes = 1 << up->cdim;
  up->num_phase_nodes = 3 * up->num_conf_nodes;
  up->qbym = inp->qbym;
  up->use_gpu = true;
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  struct gkyl_dg_pkpm_harm_exp *up_cu =
    (struct gkyl_dg_pkpm_harm_exp *)gkyl_cu_malloc(sizeof(struct gkyl_dg_pkpm_harm_exp));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_dg_pkpm_harm_exp), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;
  return up;
}
