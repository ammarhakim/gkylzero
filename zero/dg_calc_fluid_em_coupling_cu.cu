/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_fluid_em_coupling.h>
#include <gkyl_dg_calc_fluid_em_coupling_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_fluid_em_coupling_set_cu_kernel(gkyl_dg_calc_fluid_em_coupling* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], struct gkyl_array* em)
{
  int num_fluids = up->num_fluids;
  double *fluids[GKYL_MAX_SPECIES];
  const double *app_accels[GKYL_MAX_SPECIES];  
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    for (int n=0; n<num_fluids; ++n) {
      fluids[n] = (double*) gkyl_array_fetch(fluid[n], loc);
      app_accels[n] = (const double*) gkyl_array_cfetch(app_accel[n], loc);
    }
    const double *ext_em_d = (const double*) gkyl_array_cfetch(ext_em, loc);
    const double *app_current_d = (const double*) gkyl_array_cfetch(app_current, loc);
    double *em_d = (double*) gkyl_array_fetch(em, loc);

    up->fluid_em_coupling_set(linc1, up->num_fluids, up->qbym, up->epsilon0, dt, 
      up->As, up->xs, 
      app_accels, ext_em_d, app_current_d,
      fluids, em_d);
  }
}

__global__ static void
gkyl_dg_calc_fluid_em_coupling_copy_cu_kernel(gkyl_dg_calc_fluid_em_coupling* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], struct gkyl_array* em)
{
  int num_fluids = up->num_fluids;
  double *fluids[GKYL_MAX_SPECIES];  
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    for (int n=0; n<num_fluids; ++n) {
      fluids[n] = (double*) gkyl_array_fetch(fluid[n], loc);
    }
    double *em_d = (double*) gkyl_array_fetch(em, loc);

    up->fluid_copy(count, up->num_fluids, up->qbym, up->epsilon0, 
      up->xs, fluids, em_d);
  }
}

// Host-side wrapper for primitive variable calculation
void gkyl_dg_calc_fluid_em_coupling_advance_cu(struct gkyl_dg_calc_fluid_em_coupling *up, double dt, 
  const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], struct gkyl_array* em)
{
  struct gkyl_range conf_range = up->mem_range;

  int num_fluids = up->num_fluids;
  struct gkyl_array* fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array* app_accels[GKYL_MAX_SPECIES];
  for (int n=0; n<num_fluids; ++n) {
    fluids[n] = fluid[n]->on_dev;
    app_accels[n] = app_accel[n]->on_dev;
  }

  gkyl_dg_calc_fluid_em_coupling_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev, dt
    up->As->on_dev, up->xs->on_dev, conf_range, 
    app_accels, ext_em->on_dev, app_current->on_dev,
    fluids, em->on_dev);

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  gkyl_dg_calc_fluid_em_coupling_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, fluids, em->on_dev);
}

// CUDA kernel to set device pointers to fluid vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_fluid_em_coupling_set_cu_dev_ptrs(struct gkyl_dg_calc_fluid_em_coupling *up, const struct gkyl_wv_eqn *wv_eqn, 
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  up->fluid_em_coupling_set = choose_fluid_em_coupling_set_kern(b_type, cdim, poly_order);
  up->fluid_em_coupling_copy = choose_fluid_em_coupling_copy_kern(b_type, cdim, poly_order);
  up->fluid_em_coupling_energy = choose_fluid_em_coupling_energy_kern(b_type, cdim, poly_order);
}

gkyl_dg_calc_fluid_em_coupling*
gkyl_dg_calc_fluid_em_coupling_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *mem_range, 
  int num_fluids, double qbym[GKYL_MAX_SPECIES], double epsilon)
{
  struct gkyl_dg_calc_fluid_em_coupling *up = (struct gkyl_dg_calc_fluid_em_coupling*) gkyl_malloc(sizeof(gkyl_dg_calc_fluid_em_coupling));

  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->mem_range = *mem_range;

  // Linear system size is nc*(3*num_fluids + 3)
  up->num_fluids = num_fluids;
  up->As = gkyl_nmat_cu_dev_new(mem_range->volume, nc*(3*up->num_fluids + 3), nc*(3*up->num_fluids + 3));
  up->xs = gkyl_nmat_cu_dev_new(mem_range->volume, nc*(3*up->num_fluids + 3), 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  // Needed constants for the source solve
  up->epsilon0 = epsilon0;
  for (int n = 0; n < num_fluids; ++n) {
    up->qbym[n] = qbym[n];
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_fluid_em_coupling *up_cu = (struct gkyl_dg_calc_fluid_em_coupling*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_fluid_em_coupling));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_fluid_em_coupling), GKYL_CU_MEMCPY_H2D);

  dg_calc_fluid_em_coupling_set_cu_dev_ptrs<<<1,1>>>(up_cu, wv_eqn->on_dev, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;
}
