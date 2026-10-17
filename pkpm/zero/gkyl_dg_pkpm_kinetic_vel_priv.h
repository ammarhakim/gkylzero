// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef double (*pkpm_kinetic_vel_vol_t)(
  const double *w, const double *dxv, const double *a_P, const double *f,
  double *GKYL_RESTRICT out
);

typedef double (*pkpm_accel_surf_t)(
  const double *w, const double *dxv, const double *flux_l, const double *flux_r,
  double *GKYL_RESTRICT out
);

typedef double (*pkpm_accel_boundary_surf_t)(
  const double *w, const double *dxv, const int edge, const double *flux,
  double *GKYL_RESTRICT out
);

struct dg_pkpm_kinetic_vel {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Configuration-space dimensions
  int pdim; // Phase-space dimensions
  pkpm_kinetic_vel_vol_t vol_kern; // Volume kernel (pressure-force term)
  pkpm_accel_surf_t accel_surf; // Surface kernel scattering stored v_par fluxes
  pkpm_accel_boundary_surf_t accel_boundary_surf; // Velocity domain edge version
  struct gkyl_range conf_range; // Configuration-space range
  struct gkyl_range phase_range; // Phase-space range
  struct gkyl_dg_pkpm_kinetic_vel_auxfields auxfields; // Auxiliary fields.
};

// Kernel tables (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_kinetic_vel_vol_t ser_pkpm_kinetic_vel_vol_kernels[] = {
  pkpm_kinetic_vel_vol_1x1v_ser_p1,
  pkpm_kinetic_vel_vol_2x1v_ser_p1,
  pkpm_kinetic_vel_vol_3x1v_ser_p1,
};

GKYL_CU_D static const pkpm_accel_surf_t ser_pkpm_accel_surf_kernels[] = {
  pkpm_accel_surfvpar_1x1v_ser_p1,
  pkpm_accel_surfvpar_2x1v_ser_p1,
  pkpm_accel_surfvpar_3x1v_ser_p1,
};

GKYL_CU_D static const pkpm_accel_boundary_surf_t ser_pkpm_accel_boundary_surf_kernels[] = {
  pkpm_accel_boundary_surfvpar_1x1v_ser_p1,
  pkpm_accel_boundary_surfvpar_2x1v_ser_p1,
  pkpm_accel_boundary_surfvpar_3x1v_ser_p1,
};

GKYL_CU_D
static void
choose_pkpm_kinetic_vel_kerns(
  struct dg_pkpm_kinetic_vel *pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      pkpm->vol_kern = ser_pkpm_kinetic_vel_vol_kernels[cdim - 1];
      pkpm->accel_surf = ser_pkpm_accel_surf_kernels[cdim - 1];
      pkpm->accel_boundary_surf = ser_pkpm_accel_boundary_surf_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
}

// Volume term: the pressure-force term in v_par.
GKYL_CU_D
static double
vol(
  const struct gkyl_dg_eqn *eqn, const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  long cidx = gkyl_range_idx(&pkpm->conf_range, idx);
  return pkpm->vol_kern(
    xc, dx, (const double *)gkyl_array_cfetch(pkpm->auxfields.a_P, cidx), qIn, qRhsOut
  );
}

// Surface term in v_par: each cell owns its lower face; the upper face is the lower
// face of the cell above.
GKYL_CU_D
static double
surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR, const int *idxL, const int *idxC,
  const int *idxR, const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  if (dir == pkpm->cdim) {
    long pidxC = gkyl_range_idx(&pkpm->phase_range, idxC);
    long pidxR = gkyl_range_idx(&pkpm->phase_range, idxR);
    const double *flux_l = (const double *)gkyl_array_cfetch(pkpm->auxfields.vel_flux_surf, pidxC);
    const double *flux_r = (const double *)gkyl_array_cfetch(pkpm->auxfields.vel_flux_surf, pidxR);
    return pkpm->accel_surf(xcC, dxC, flux_l, flux_r, qRhsOut);
  }
  return 0.0;
}

// Velocity domain edges (no ghost cells): the lowest cell uses the lower-face flux
// of the cell above it (edge = -1), the uppermost cell its own lower-face flux
// (edge = +1); the domain-edge face itself carries zero flux.
GKYL_CU_D
static double
boundary_surf(
  const struct gkyl_dg_eqn *eqn, int dir, const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin, const int *idxEdge, const int *idxSkin,
  const int edge, const double *qInEdge, const double *qInSkin, double *GKYL_RESTRICT qRhsOut
)
{
  struct dg_pkpm_kinetic_vel *pkpm = container_of(eqn, struct dg_pkpm_kinetic_vel, eqn);
  if (dir == pkpm->cdim) {
    long pidx = gkyl_range_idx(&pkpm->phase_range, edge == -1 ? idxEdge : idxSkin);
    const double *flux = (const double *)gkyl_array_cfetch(pkpm->auxfields.vel_flux_surf, pidx);
    return pkpm->accel_boundary_surf(xcSkin, dxSkin, edge, flux, qRhsOut);
  }
  return 0.0;
}
