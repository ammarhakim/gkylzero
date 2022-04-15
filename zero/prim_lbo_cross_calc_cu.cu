/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_cross_calc_priv.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_prim_lbo_cross_calc_set_cu_ker(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis cbasis, struct gkyl_range conf_rng,
  int nspecies, const struct gkyl_array* greene[GKYL_MAX_SPECIES], const double self_m,
  const struct gkyl_array* self_u, const struct gkyl_array* self_vtsq,
  const double cross_m[GKYL_MAX_SPECIES], const struct gkyl_array* cross_u[GKYL_MAX_SPECIES],
  const struct gkyl_array* cross_vtsq[GKYL_MAX_SPECIES],
  const struct gkyl_array* moms, const struct gkyl_array* boundary_corrections)
{
  long count = 0;
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_rng.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_rng, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&conf_rng, idx);

    for (int n=1; n<nspecies; ++n) {
      struct gkyl_mat lhs = gkyl_nmat_get(As, count);
      struct gkyl_mat rhs = gkyl_nmat_get(xs, count);
      const double *self_u_d = (const double*) gkyl_array_cfetch(self_u, linc);
      const double *greene_d = (const double*) gkyl_array_cfetch(greene[n], linc);
      const double *self_vtsq_d = (const double*) gkyl_array_cfetch(self_vtsq, linc);
      const double *cross_u_d = (const double*) gkyl_array_cfetch(cross_u[n], linc);
      const double *cross_vtsq_d = (const double*) gkyl_array_cfetch(cross_vtsq[n], linc);
      const double *moms_d = (const double*) gkyl_array_cfetch(moms, linc);
      const double *boundary_corrections_d = (const double*) gkyl_array_cfetch(boundary_corrections, linc);

      gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

      calc->prim->cross_prim(calc->prim, &lhs, &rhs, greene_d, self_m,
        self_u_d, self_vtsq_d,
        cross_m[n], cross_u_d, cross_vtsq_d, moms_d,
        boundary_corrections_d
      );
      count += 1;
    }
  }
}

__global__ static void
gkyl_prim_lbo_copy_sol_cu_ker(struct gkyl_nmat *xs,
  struct gkyl_basis cbasis, struct gkyl_range conf_rng,
  int nc, int vdim, int nspecies,
  struct gkyl_array* u_out[GKYL_MAX_DIM], struct gkyl_array* vtsq_out[GKYL_MAX_DIM])
{
  long count = 0;
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_rng.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_rng, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&conf_rng, idx);

    for (int n=0; n<nspecies; ++n) {
      struct gkyl_mat out_d = gkyl_nmat_get(xs, count);
      double *u_d = (double*) gkyl_array_fetch(u_out[n], linc);
      double *vtsq_d = (double*) gkyl_array_fetch(vtsq_out[n], linc);
    
      prim_lbo_copy_sol(&out_d, nc, vdim, u_d, vtsq_d);
      count += 1;
    }
  }
}

void
gkyl_prim_lbo_cross_calc_advance_cu(gkyl_prim_lbo_cross_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range conf_rng, struct gkyl_array* greene[GKYL_MAX_SPECIES], const double self_m,
  const struct gkyl_array* self_u, const struct gkyl_array* self_vtsq,
  const double cross_m[GKYL_MAX_SPECIES], struct gkyl_array* cross_u[GKYL_MAX_SPECIES],
  struct gkyl_array* cross_vtsq[GKYL_MAX_SPECIES], const struct gkyl_array* moms,
  const struct gkyl_array* boundary_corrections, struct gkyl_array* u_out[GKYL_MAX_SPECIES],
  struct gkyl_array* vtsq_out[GKYL_MAX_SPECIES])
{
  int nspecies = calc->nspecies;

  const struct gkyl_array *cross_us[GKYL_MAX_SPECIES];
  const struct gkyl_array *cross_vtsqs[GKYL_MAX_SPECIES];
  struct gkyl_array *u_outs[GKYL_MAX_SPECIES];
  struct gkyl_array *vtsq_outs[GKYL_MAX_SPECIES];
  const struct gkyl_array *greenes[GKYL_MAX_SPECIES];
  
  for (int n=0; n<nspecies; ++n) {
    cross_us[n] = cross_u[n]->on_dev;
    cross_vtsqs[n] = cross_vtsq[n]->on_dev;
    u_outs[n] = u_out[n]->on_dev;
    vtsq_outs[n] = vtsq_out[n]->on_dev;
    greenes[n] = greene[n]->on_dev;
  }

  int nc = cbasis.num_basis;
  int vdim = calc->prim->pdim - calc->prim->cdim;
  int N = nc*(vdim + 1);
  
  if (calc->is_first) {
    calc->As = gkyl_nmat_cu_dev_new(nspecies*conf_rng.volume, N, N);
    calc->xs = gkyl_nmat_cu_dev_new(nspecies*conf_rng.volume, N, 1);
    calc->mem = gkyl_nmat_linsolve_lu_cu_dev_new(calc->As->num, calc->As->nr);
    calc->is_first = false;
  }

  gkyl_prim_lbo_cross_calc_set_cu_ker<<<conf_rng.nblocks, conf_rng.nthreads>>>(calc->on_dev,
    calc->As->on_dev, calc->xs->on_dev, cbasis, conf_rng,
    nspecies, greenes, self_m, self_u->on_dev, self_vtsq->on_dev,
    cross_m, cross_us, cross_vtsqs,
    moms->on_dev, boundary_corrections->on_dev);
  
  bool status = gkyl_nmat_linsolve_lu_pa(calc->mem, calc->As, calc->xs);

  gkyl_prim_lbo_copy_sol_cu_ker<<<conf_rng.nblocks, conf_rng.nthreads>>>(calc->xs->on_dev,
    cbasis, conf_rng, nc, vdim, nspecies,
    u_outs, vtsq_outs);
}

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim, int nspecies)
{
  gkyl_prim_lbo_cross_calc *up = (gkyl_prim_lbo_cross_calc*) gkyl_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  up->grid = *grid;
  up->prim = prim;
  up->nspecies = nspecies;

  up->is_first = true;
  up->As = up->xs = 0;
  up->mem = 0;

  struct gkyl_prim_lbo_type *pt = gkyl_prim_lbo_type_acquire(prim);
  up->prim = pt->on_dev; // so memcpy below gets dev copy

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_prim_lbo_cross_calc *up_cu = (gkyl_prim_lbo_cross_calc*) gkyl_cu_malloc(sizeof(gkyl_prim_lbo_cross_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_prim_lbo_cross_calc), GKYL_CU_MEMCPY_H2D);

  up->prim = pt; // host portion of struct should have host copy
  up->on_dev = up_cu; // host pointer
  
  return up;
}
