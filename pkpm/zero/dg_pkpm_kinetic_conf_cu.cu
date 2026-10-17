/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_dg_pkpm_kinetic_conf_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set pointers to auxiliary fields.
__global__ static void
gkyl_dg_pkpm_kinetic_conf_set_auxfields_cu_kernel(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin
)
{
  struct dg_pkpm_kinetic_conf *pkpm = container_of(eqn, struct dg_pkpm_kinetic_conf, eqn);
  pkpm->auxfields.bvar = auxin.bvar;
  pkpm->auxfields.grad_b = auxin.grad_b;
  pkpm->auxfields.u = auxin.u;
  pkpm->auxfields.grad_u = auxin.grad_u;
  pkpm->auxfields.theta = auxin.theta;
  pkpm->auxfields.nu_sum = auxin.nu_sum;
  pkpm->auxfields.nu_vthsq_sum = auxin.nu_vthsq_sum;
  pkpm->auxfields.conf_flux_surf = auxin.conf_flux_surf;
}

// Host-side wrapper for set_auxfields_cu_kernel: the arrays must be device arrays.
void
gkyl_dg_pkpm_kinetic_conf_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin
)
{
  struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin_dev = {
    .bvar = auxin.bvar ? auxin.bvar->on_dev : 0,
    .grad_b = auxin.grad_b ? auxin.grad_b->on_dev : 0,
    .u = auxin.u ? auxin.u->on_dev : 0,
    .grad_u = auxin.grad_u ? auxin.grad_u->on_dev : 0,
    .theta = auxin.theta ? auxin.theta->on_dev : 0,
    .nu_sum = auxin.nu_sum ? auxin.nu_sum->on_dev : 0,
    .nu_vthsq_sum = auxin.nu_vthsq_sum ? auxin.nu_vthsq_sum->on_dev : 0,
    .conf_flux_surf = auxin.conf_flux_surf ? auxin.conf_flux_surf->on_dev : 0,
  };
  gkyl_dg_pkpm_kinetic_conf_set_auxfields_cu_kernel<<<1, 1>>>(eqn, auxin_dev);
}

// CUDA kernel to set device function pointers.
__global__ static void
dg_pkpm_kinetic_conf_set_cu_dev_ptrs(
  struct dg_pkpm_kinetic_conf *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;
  pkpm->vol_kern = choose_pkpm_kinetic_conf_vol_kern(b_type, cdim, poly_order);
  for (int d = 0; d < cdim; ++d) {
    pkpm->stream_surf[d] = choose_pkpm_stream_surf_kern(b_type, cdim, poly_order, d);
  }
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_kinetic_conf_cu_dev_inew(const struct gkyl_dg_pkpm_kinetic_conf_inp *inp)
{
  struct dg_pkpm_kinetic_conf *pkpm =
    (struct dg_pkpm_kinetic_conf *)gkyl_malloc(sizeof(struct dg_pkpm_kinetic_conf));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  pkpm->cdim = cdim;
  pkpm->pdim = pdim;
  pkpm->num_surf_basis = inp->phase_basis->num_basis / 2;
  pkpm->eqn.num_equations = 2;
  pkpm->conf_range = *inp->conf_range;
  pkpm->phase_range = *inp->phase_range;
  pkpm->auxfields = inp->auxfields;

  pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(pkpm->eqn.flags);
  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_kinetic_conf_free);

  struct dg_pkpm_kinetic_conf *pkpm_cu =
    (struct dg_pkpm_kinetic_conf *)gkyl_cu_malloc(sizeof(struct dg_pkpm_kinetic_conf));
  gkyl_cu_memcpy(pkpm_cu, pkpm, sizeof(struct dg_pkpm_kinetic_conf), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_kinetic_conf_set_cu_dev_ptrs<<<1, 1>>>(
    pkpm_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  pkpm->eqn.on_dev = &pkpm_cu->eqn;
  return &pkpm->eqn;
}
