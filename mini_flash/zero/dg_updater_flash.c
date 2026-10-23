
#include <gkyl_alloc.h>
#include <gkyl_dg_updater_flash.h>
#include <gkyl_dg_updater_flash_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

gkyl_dg_updater_flash*
gkyl_dg_updater_flash_new(const struct gkyl_rect_grid *conf_grid,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_basis *conf_basis, 
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range, 
  const struct gkyl_range *phase_range, const bool *is_zero_flux_dir,
  struct gkyl_dg_vlasov_pkpm_auxfields *pkpm_aux)
{
  gkyl_dg_updater_flash *up = gkyl_malloc(sizeof(gkyl_dg_updater_flash));
  bool use_gpu = 0; // HARDCODE: No GPUs

  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;

  int up_dirs_conf[GKYL_MAX_DIM], zero_flux_flags_conf[GKYL_MAX_DIM];
  int up_dirs_phase[GKYL_MAX_DIM], zero_flux_flags_phase[GKYL_MAX_DIM];
    for (int d=0; d<cdim; ++d) {
    up_dirs_conf[d] = d;
    up_dirs_phase[d] = d;
    zero_flux_flags_conf[d] = zero_flux_flags_conf[d+cdim] = is_zero_flux_dir[d] ? 1 : 0;
    zero_flux_flags_phase[d] = zero_flux_flags_phase[d+pdim] = is_zero_flux_dir[d] ? 1 : 0;
  }
  int num_up_dirs_conf = cdim;

  for (int d=cdim; d<pdim; ++d) {
    up_dirs_phase[d] = d;
    zero_flux_flags_phase[d] = zero_flux_flags_phase[d+pdim] = 1; // zero-flux BCs in vel-space
  }

  // Flow-frame Vlasov equation from PKPM.
  up->eqn_vlasov = gkyl_dg_vlasov_pkpm_new(conf_basis, phase_basis, conf_range, phase_range, use_gpu);
  gkyl_vlasov_pkpm_set_auxfields(up->eqn_vlasov, *pkpm_aux);

  up->up_vlasov = gkyl_hyper_dg_new(phase_grid, phase_basis, up->eqn_vlasov, pdim,
    up_dirs_phase, zero_flux_flags_phase, 1, use_gpu);

  up->vlasov_tm = 0.0;

  // Set kernels for computing timestep separately from the update
  up->compute_cfl_kernels.boundary_accel = choose_flash_compute_cfl_vlasov_pkpm_boundary_accel_kernel(cdim);
  up->compute_cfl_kernels.surf_accel = choose_flash_compute_cfl_vlasov_pkpm_surf_accel_kernel(cdim);
  up->compute_cfl_kernels.surf_stream = choose_flash_compute_cfl_vlasov_pkpm_surf_kernel(cdim);
  up->compute_cfl_kernels.vol = choose_flash_compute_cfl_vlasov_pkpm_vol_kernel(cdim);
  up->compute_cfl_kernels.lbo_vol = choose_flash_compute_cfl_lbo_pkpm_vol_kernel(cdim);

  return up;
}

void
gkyl_dg_updater_flash_compute_cflrate(gkyl_dg_updater_flash *up,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_range *update_phase_rng,
  const struct gkyl_range *update_conf_rng, const struct gkyl_array *pkpm_accel, 
  const struct gkyl_array *pkpm_lax, const struct gkyl_array *div_b, 
  const struct gkyl_array *max_b, struct gkyl_array *cflrate)
{
  // Compute maximum stable timestep without actually performing the update
  int pdim = update_phase_rng->ndim;
  int cdim = update_conf_rng->ndim;

  double xcc[GKYL_MAX_DIM];
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  int edge;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_phase_rng);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(phase_grid, idxc, xcc);

    long linc = gkyl_range_idx(update_phase_rng, idxc);
    long linc_conf = gkyl_range_idx(update_conf_rng, idxc);

    const double *pkpm_accel_c = (const double *)gkyl_array_cfetch(pkpm_accel, linc_conf);
    const double *div_b_c = (const double *)gkyl_array_cfetch(div_b, linc_conf);
    const double *max_b_c = (const double *)gkyl_array_cfetch(max_b, linc_conf);

    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, linc);

    double cflr, cfls = 0.0;
 
    // Volume term contribution
    cflr = up->compute_cfl_kernels.vol(xcc, phase_grid->dx, pkpm_accel_c);

    for (int d=0; d<pdim; ++d) {
      if (d < cdim) {
        // Stream update in configuration space
        gkyl_copy_int_arr(pdim, idxc, idxl);
        gkyl_copy_int_arr(pdim, idxc, idxr);
        idxl[d] = idxl[d] - 1;
        idxr[d] = idxr[d] + 1;
        long linl_conf = gkyl_range_idx(update_conf_rng, idxl);
        long linr_conf = gkyl_range_idx(update_conf_rng, idxr);
        const double *pkpm_lax_l = (const double *)gkyl_array_cfetch(pkpm_lax, linl_conf);
        const double *pkpm_lax_r = (const double *)gkyl_array_cfetch(pkpm_lax, linr_conf);

        cfls = up->compute_cfl_kernels.surf_stream(xcc, phase_grid->dx, 
          max_b_c,  pkpm_lax_l, pkpm_lax_r, pkpm_accel_c);
      }
      else {
        // Accel update in vpar
        if ((idxc[d] == update_phase_rng->lower[d]) ||
            (idxc[d] == update_phase_rng->upper[d])) {
          edge = (idxc[d] == update_phase_rng->lower[d]) ? -1 : 1;

          cfls = up->compute_cfl_kernels.boundary_accel(xcc, phase_grid->dx, div_b_c, pkpm_accel_c, edge);
        }
        else {
          // Accel update
          cfls = up->compute_cfl_kernels.surf_accel(xcc, phase_grid->dx, pkpm_accel_c);
        }
      }
      cflrate_d[0] += cflr + cfls;
    }
  }
}

void
gkyl_dg_updater_flash_compute_cflrate_lbo(gkyl_dg_updater_flash *up,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_range *update_phase_rng,
  const struct gkyl_range *update_conf_rng, const struct gkyl_array *nuSum, 
  const struct gkyl_array *nuPrimMomsSum, struct gkyl_array *cflrate)
{
  int pdim = update_phase_rng->ndim;
  int cdim = update_conf_rng->ndim;

  double xcc[GKYL_MAX_DIM];
  int idxc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_phase_rng);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(phase_grid, idxc, xcc);

    long linc = gkyl_range_idx(update_phase_rng, idxc);
    long linc_conf = gkyl_range_idx(update_conf_rng, idxc);

    const double *nuSum_c = (const double *)gkyl_array_cfetch(nuSum, linc_conf);
    const double *nuPrimMomsSum_c = (const double *)gkyl_array_cfetch(nuPrimMomsSum, linc_conf);

    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, linc);
 
    // Volume term contribution
    double cflr = up->compute_cfl_kernels.lbo_vol(xcc, phase_grid->dx, nuSum_c, nuPrimMomsSum_c);
    cflrate_d[0] += cflr;
  }
}

void
gkyl_dg_updater_flash_advance(gkyl_dg_updater_flash *up, 
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng, 
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate_f,
  struct gkyl_array* GKYL_RESTRICT rhs_f)
{
  // Compute RHS of Vlasov update
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->up_vlasov, update_phase_rng, fIn, cflrate_f, rhs_f);
  up->vlasov_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_dg_updater_flash_release(gkyl_dg_updater_flash *up)
{
  gkyl_dg_eqn_release(up->eqn_vlasov);
  gkyl_hyper_dg_release(up->up_vlasov);
  gkyl_free(up);
}