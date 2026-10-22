#include <gkyl_alloc.h>
#include <gkyl_flash_interp_vars_priv.h>

struct gkyl_flash_interp_vars*
gkyl_flash_interp_vars_new(struct gkyl_basis *basis, 
  struct gkyl_rect_grid *grid_kinetic, struct gkyl_rect_grid *grid_fluid,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid) 
{
  struct gkyl_flash_interp_vars *up = gkyl_malloc(sizeof(struct gkyl_flash_interp_vars));

  up->cdim = grid_fluid->ndim;
  up->basis = basis;

  up->grid_fluid = grid_fluid;
  up->grid_kinetic = grid_kinetic;

  // Interpolation updater
  up->interp = gkyl_dg_fv_interpolate_new(basis, grid_kinetic, grid_fluid);

  return up;
}

// Interpolate momentum density from finite volume grid to kinetic grid.
void
gkyl_flash_interp_vars_rhou_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_fluid_in, struct gkyl_array *dg_rhou_out)
{
  struct gkyl_array *tmp = gkyl_array_new(GKYL_DOUBLE, 3, range_fluid->volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range_fluid);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(range_fluid, iter.idx);
    const double *fluid_d = (const double *)gkyl_array_cfetch(fv_fluid_in, linc);
    double *tmp_d = (double *)gkyl_array_fetch(tmp, linc);

    for (int d=0; d<3; ++d) {
      tmp_d[d] = fluid_d[d+1]/fluid_d[0];
    }
  }

  // Interpolate momentum density vector onto DG grid.
  for (int i=0; i<3; ++i) {
    gkyl_dg_fv_interpolate_fv_to_modal_recov_advance(up->interp,
    range_kinetic, range_fluid,
    tmp, i, dg_rhou_out, i);
    }

  gkyl_array_release(tmp);
}

// Interpolate fluid velocity vector from finite volume grid to kinetic grid.
void
gkyl_flash_interp_vars_vel_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_fluid_vel_in, struct gkyl_array *dg_vel_out)
{
  // Interpolate momentum density vector onto DG grid.
  for (int i=0; i<3; ++i) {
    gkyl_dg_fv_interpolate_fv_to_modal_recov_advance(up->interp,
    range_kinetic, range_fluid,
    fv_fluid_vel_in, i, dg_vel_out, i);
  }
}

// Interpolate initial condition arrays [n, ux, uy, uz, T/m] onto kinetic grid.
void
gkyl_flash_interp_vars_init_dist_moms_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_moms_in, struct gkyl_array *dg_moms_out)
{
  // Interpolate moments onto DG grid.
  for (int i=0; i<5; ++i) {
    gkyl_dg_fv_interpolate_fv_to_modal_recov_advance(up->interp,
    range_kinetic, range_fluid,
    fv_moms_in, i, dg_moms_out, i);
  }
}

void
gkyl_flash_interp_vars_coupling_fluxes_to_fv(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *dg_coupling_fluxes, struct gkyl_array *fv_coupling_fluxes) 
{
  for (int i=0; i<3; ++i) {
    gkyl_dg_fv_interpolate_modal_to_fv_advance(up->interp,
      range_kinetic, range_fluid,
      dg_coupling_fluxes, i, fv_coupling_fluxes, i);
  }
}

void
gkyl_flash_interp_vars_release(struct gkyl_flash_interp_vars *up)
{
  gkyl_dg_fv_interpolate_release(up->interp);
  gkyl_free(up);
}