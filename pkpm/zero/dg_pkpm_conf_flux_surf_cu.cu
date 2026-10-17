/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_conf_flux_surf_priv.h>
#include <gkyl_util.h>
}

// One thread per interior phase-space cell: the cell's lower faces in every
// configuration direction, plus the upper ghost-owned face at the domain edge
// (written by exactly one interior cell per direction, so no atomics are needed).
__global__ void
gkyl_dg_pkpm_conf_flux_surf_advance_cu_kernel(
  struct gkyl_dg_pkpm_conf_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, struct gkyl_range phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *mass_flux
)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int nflux = 2 * up->num_surf_basis;
  int nmf = up->num_face_basis;
  const struct gkyl_array *th = up->acoustic_floor ? theta : 0;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xc_r[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *b_c = (const double *)gkyl_array_cfetch(bvar, cidx);
    const double *u_c = (const double *)gkyl_array_cfetch(u, cidx);
    const double *th_c = th ? (const double *)gkyl_array_cfetch(th, cidx) : 0;
    const double *f_c = (const double *)gkyl_array_cfetch(fin, pidx);
    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, pidx);
    double *flux = (double *)gkyl_array_fetch(conf_flux_surf, pidx);
    double *mf = (double *)gkyl_array_fetch(mass_flux, pidx);

    for (int dir = 0; dir < cdim; ++dir) {
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long cidx_l = gkyl_range_idx(&conf_range, idx_l);
      long pidx_l = gkyl_range_idx(&phase_range, idx_l);
      cflrate_d[0] += up->conf_flux[dir](
        xc, up->phase_grid.dx, (const double *)gkyl_array_cfetch(bvar, cidx_l), b_c,
        (const double *)gkyl_array_cfetch(u, cidx_l), u_c,
        th ? (const double *)gkyl_array_cfetch(th, cidx_l) : 0, th_c,
        (const double *)gkyl_array_cfetch(fin, pidx_l), f_c, &flux[dir * nflux], &mf[dir * nmf]
      );

      if (idx[dir] == phase_range.upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(&conf_range, idx_r);
        long pidx_r = gkyl_range_idx(&phase_range_ext, idx_r);
        gkyl_rect_grid_cell_center(&up->phase_grid, idx_r, xc_r);
        double *flux_r = (double *)gkyl_array_fetch(conf_flux_surf, pidx_r);
        double *mf_r = (double *)gkyl_array_fetch(mass_flux, pidx_r);
        cflrate_d[0] += up->conf_flux[dir](
          xc_r, up->phase_grid.dx, b_c, (const double *)gkyl_array_cfetch(bvar, cidx_r), u_c,
          (const double *)gkyl_array_cfetch(u, cidx_r), th_c,
          th ? (const double *)gkyl_array_cfetch(th, cidx_r) : 0, f_c,
          (const double *)gkyl_array_cfetch(fin, pidx_r), &flux_r[dir * nflux], &mf_r[dir * nmf]
        );
      }
    }
  }
}

// One thread per interior configuration cell: sums the per-cell kinetic mass fluxes
// over the velocity column (deterministic, no atomics), including the upper
// ghost-owned faces.
__global__ void
gkyl_dg_pkpm_conf_flux_surf_reduce_cu_kernel(
  struct gkyl_dg_pkpm_conf_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, struct gkyl_range phase_range_ext,
  const struct gkyl_array *mass_flux, struct gkyl_array *mass_flux_face
)
{
  int cdim = up->cdim;
  int nmf = up->num_face_basis;
  int idx[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    double *face = (double *)gkyl_array_fetch(mass_flux_face, cidx);
    for (int iv = phase_range.lower[cdim]; iv <= phase_range.upper[cdim]; ++iv) {
      idx[cdim] = iv;
      const double *mf = (const double *)gkyl_array_cfetch(mass_flux, gkyl_range_idx(&phase_range, idx));
      for (int k = 0; k < cdim * nmf; ++k) {
        face[k] += mf[k];
      }
    }
    for (int dir = 0; dir < cdim; ++dir) {
      if (idx[dir] == conf_range.upper[dir]) {
        gkyl_copy_int_arr(cdim, idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        double *face_r = (double *)gkyl_array_fetch(mass_flux_face, gkyl_range_idx(&conf_range, idx_r));
        for (int iv = phase_range.lower[cdim]; iv <= phase_range.upper[cdim]; ++iv) {
          idx_r[cdim] = iv;
          const double *mf_r =
            (const double *)gkyl_array_cfetch(mass_flux, gkyl_range_idx(&phase_range_ext, idx_r));
          for (int k = 0; k < nmf; ++k) {
            face_r[dir * nmf + k] += mf_r[dir * nmf + k];
          }
        }
      }
    }
  }
}

// One thread per interior phase-space cell: corrects the stored fluxes of its lower
// faces and of the upper ghost-owned faces (same ownership as the flux kernel).
__global__ void
gkyl_dg_pkpm_conf_flux_surf_corr_cu_kernel(
  struct gkyl_dg_pkpm_conf_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, struct gkyl_range phase_range_ext, double mass,
  const struct gkyl_array *moms, const struct gkyl_array *roe_mass_flux,
  const struct gkyl_array *mass_flux_face, const struct gkyl_array *fin,
  struct gkyl_array *conf_flux_surf
)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int nflux = 2 * up->num_surf_basis;
  int nmf = up->num_face_basis;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);

    const double *moms_c = (const double *)gkyl_array_cfetch(moms, cidx);
    const double *roe_c = (const double *)gkyl_array_cfetch(roe_mass_flux, cidx);
    const double *kin_c = (const double *)gkyl_array_cfetch(mass_flux_face, cidx);
    const double *f_c = (const double *)gkyl_array_cfetch(fin, pidx);
    double *flux = (double *)gkyl_array_fetch(conf_flux_surf, pidx);

    for (int dir = 0; dir < cdim; ++dir) {
      bool skip = up->skip_domain_faces[dir];
      if (!(skip && idx[dir] == phase_range.lower[dir])) {
        gkyl_copy_int_arr(pdim, idx, idx_l);
        idx_l[dir] = idx_l[dir] - 1;
        long cidx_l = gkyl_range_idx(&conf_range, idx_l);
        long pidx_l = gkyl_range_idx(&phase_range, idx_l);
        up->conf_flux_corr[dir](
          mass, (const double *)gkyl_array_cfetch(moms, cidx_l), moms_c, &roe_c[dir * nmf],
          &kin_c[dir * nmf], (const double *)gkyl_array_cfetch(fin, pidx_l), f_c, &flux[dir * nflux]
        );
      }

      if (idx[dir] == phase_range.upper[dir] && !skip) {
        gkyl_copy_int_arr(pdim, idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(&conf_range, idx_r);
        long pidx_r = gkyl_range_idx(&phase_range_ext, idx_r);
        const double *roe_r = (const double *)gkyl_array_cfetch(roe_mass_flux, cidx_r);
        const double *kin_r = (const double *)gkyl_array_cfetch(mass_flux_face, cidx_r);
        double *flux_r = (double *)gkyl_array_fetch(conf_flux_surf, pidx_r);
        up->conf_flux_corr[dir](
          mass, moms_c, (const double *)gkyl_array_cfetch(moms, cidx_r), &roe_r[dir * nmf],
          &kin_r[dir * nmf], f_c, (const double *)gkyl_array_cfetch(fin, pidx_r), &flux_r[dir * nflux]
        );
      }
    }
  }
}

void
gkyl_dg_pkpm_conf_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *mass_flux
)
{
  gkyl_dg_pkpm_conf_flux_surf_advance_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, *phase_range_ext, bvar->on_dev, u->on_dev,
    theta ? theta->on_dev : 0, fin->on_dev, cflrate->on_dev, conf_flux_surf->on_dev,
    mass_flux->on_dev
  );
}

void
gkyl_dg_pkpm_conf_flux_surf_constrain_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext, double mass,
  const struct gkyl_array *moms, const struct gkyl_array *roe_mass_flux,
  const struct gkyl_array *mass_flux, struct gkyl_array *mass_flux_face,
  const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf
)
{
  gkyl_array_clear(mass_flux_face, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_reduce_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, *phase_range_ext, mass_flux->on_dev,
    mass_flux_face->on_dev
  );
  gkyl_dg_pkpm_conf_flux_surf_corr_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, *phase_range_ext, mass, moms->on_dev,
    roe_mass_flux->on_dev, mass_flux_face->on_dev, fin->on_dev, conf_flux_surf->on_dev
  );
}

__global__ static void
dg_pkpm_conf_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_pkpm_conf_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  for (int d = 0; d < cdim; ++d) {
    up->conf_flux[d] = choose_pkpm_conf_flux_kern(b_type, cdim, poly_order, d);
    up->conf_flux_corr[d] = choose_pkpm_conf_flux_corr_kern(b_type, cdim, poly_order, d);
  }
}

struct gkyl_dg_pkpm_conf_flux_surf *
gkyl_dg_pkpm_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_conf_flux_surf_inp *inp)
{
  struct gkyl_dg_pkpm_conf_flux_surf *up =
    (struct gkyl_dg_pkpm_conf_flux_surf *)gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->acoustic_floor = inp->acoustic_floor;
  for (int d = 0; d < GKYL_MAX_CDIM; ++d) {
    up->skip_domain_faces[d] = inp->skip_domain_faces[d];
  }
  up->use_gpu = inp->use_gpu;
  up->num_surf_basis = inp->phase_basis->num_basis / 2;
  up->num_face_basis = inp->conf_basis->num_basis / 2;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_pkpm_conf_flux_surf *up_cu =
    (struct gkyl_dg_pkpm_conf_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_conf_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  up->on_dev = up_cu;
  return up;
}
