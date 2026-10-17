/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_kinetic_vel_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set pointers to auxiliary fields.
__global__ static void
gkyl_dg_pkpm_kinetic_vel_set_auxfields_cu_kernel(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_vel_auxfields auxin
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  pkpm->auxfields.a_P = auxin.a_P;
  pkpm->auxfields.vel_flux_surf = auxin.vel_flux_surf;
}

// Host-side wrapper for set_auxfields_cu_kernel: the arrays must be device arrays.
void
gkyl_dg_pkpm_kinetic_vel_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_vel_auxfields auxin
)
{
  struct gkyl_dg_pkpm_kinetic_vel_auxfields auxin_dev = {
    .a_P = auxin.a_P ? auxin.a_P->on_dev : 0,
    .vel_flux_surf = auxin.vel_flux_surf ? auxin.vel_flux_surf->on_dev : 0,
  };
  gkyl_dg_pkpm_kinetic_vel_set_auxfields_cu_kernel<<<1, 1>>>(eqn, auxin_dev);
}

// CUDA kernel to set device function pointers.
__global__ static void
dg_pkpm_kinetic_vel_set_cu_dev_ptrs(
  struct dg_pkpm_kinetic_vel *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;
  choose_pkpm_kinetic_vel_kerns(pkpm, b_type, cdim, poly_order);
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_kinetic_vel_cu_dev_inew(const struct gkyl_dg_pkpm_kinetic_vel_inp *inp)
{
  struct dg_pkpm_kinetic_vel *pkpm =
    (struct dg_pkpm_kinetic_vel *)gkyl_malloc(sizeof(struct dg_pkpm_kinetic_vel));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  pkpm->cdim = cdim;
  pkpm->pdim = pdim;
  pkpm->eqn.num_equations = 2;
  pkpm->conf_range = *inp->conf_range;
  pkpm->phase_range = *inp->phase_range;
  pkpm->auxfields = inp->auxfields;

  pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(pkpm->eqn.flags);
  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_kinetic_vel_free);

  struct dg_pkpm_kinetic_vel *pkpm_cu =
    (struct dg_pkpm_kinetic_vel *)gkyl_cu_malloc(sizeof(struct dg_pkpm_kinetic_vel));
  gkyl_cu_memcpy(pkpm_cu, pkpm, sizeof(struct dg_pkpm_kinetic_vel), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_kinetic_vel_set_cu_dev_ptrs<<<1, 1>>>(
    pkpm_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  pkpm->eqn.on_dev = &pkpm_cu->eqn;
  return &pkpm->eqn;
}
