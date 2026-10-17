/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set pointers to auxiliary fields.
__global__ static void
gkyl_dg_pkpm_fluid_set_auxfields_cu_kernel(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_fluid_auxfields auxin
)
{
  struct dg_pkpm_fluid *pkpm = container_of(eqn, struct dg_pkpm_fluid, eqn);
  pkpm->auxfields.moms = auxin.moms;
  pkpm->auxfields.bvar = auxin.bvar;
  pkpm->auxfields.u = auxin.u;
  pkpm->auxfields.fluid_flux_surf = auxin.fluid_flux_surf;
}

// Host-side wrapper for set_auxfields_cu_kernel: the arrays must be device arrays.
void
gkyl_dg_pkpm_fluid_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_fluid_auxfields auxin
)
{
  struct gkyl_dg_pkpm_fluid_auxfields auxin_dev = {
    .moms = auxin.moms ? auxin.moms->on_dev : 0,
    .bvar = auxin.bvar ? auxin.bvar->on_dev : 0,
    .u = auxin.u ? auxin.u->on_dev : 0,
    .fluid_flux_surf = auxin.fluid_flux_surf ? auxin.fluid_flux_surf->on_dev : 0,
  };
  gkyl_dg_pkpm_fluid_set_auxfields_cu_kernel<<<1, 1>>>(eqn, auxin_dev);
}

// CUDA kernel to set device function pointers.
__global__ static void
dg_pkpm_fluid_set_cu_dev_ptrs(
  struct dg_pkpm_fluid *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  pkpm->eqn.vol_term = vol;
  pkpm->eqn.surf_term = surf;
  pkpm->eqn.boundary_surf_term = boundary_surf;
  choose_pkpm_fluid_kerns(pkpm, b_type, cdim, poly_order);
}

struct gkyl_dg_eqn *
gkyl_dg_pkpm_fluid_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_inp *inp)
{
  struct dg_pkpm_fluid *pkpm = (struct dg_pkpm_fluid *)gkyl_malloc(sizeof(struct dg_pkpm_fluid));

  int cdim = inp->conf_basis->ndim;
  pkpm->cdim = cdim;
  pkpm->num_surf_basis = inp->conf_basis->num_basis / 2;
  pkpm->eqn.num_equations = 3;
  pkpm->conf_range = *inp->conf_range;
  pkpm->auxfields = inp->auxfields;

  pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(pkpm->eqn.flags);
  pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_dg_pkpm_fluid_free);

  struct dg_pkpm_fluid *pkpm_cu = (struct dg_pkpm_fluid *)gkyl_cu_malloc(sizeof(struct dg_pkpm_fluid));
  gkyl_cu_memcpy(pkpm_cu, pkpm, sizeof(struct dg_pkpm_fluid), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_fluid_set_cu_dev_ptrs<<<1, 1>>>(
    pkpm_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  pkpm->eqn.on_dev = &pkpm_cu->eqn;
  return &pkpm->eqn;
}
