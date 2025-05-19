#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_rect_decomp.h>

struct gkyl_mirror_grid_gen_x {
  
};

// context for use in root finder
struct psirz_ctx {
  struct gkyl_basis_ops_evalf *evcub; // cubic eval functions
  double Z; // local Z value
  double psi; // psi to match
};

static
double psirz(double R, void *ctx)
{
  struct psirz_ctx *rctx = ctx;
  double Z = rctx->Z;
  double xn[2] = { R, Z };
  double fout[1];
  rctx->evcub->eval_cubic(0, xn, fout, rctx->evcub->ctx);
  return fout[0] - rctx->psi;
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  struct gkyl_mirror_grid_gen *geo = gkyl_malloc(sizeof *geo);
  geo->gg_x = gkyl_malloc(sizeof *geo->gg_x);

  int nr = inp->nrnodes, nz = inp->nznodes;
  int cells[] = { nr, nz };
  double lower[2] = { inp->R[0], inp->Z[0] };
  double upper[2] = { inp->R[1], inp->Z[1] };

  struct gkyl_rect_grid gridRZ;
  gkyl_rect_grid_init(&gridRZ, 2, lower, upper, cells);

  struct gkyl_basis_ops_evalf *evcub =
    gkyl_dg_basis_ops_evalf_new(&gridRZ, inp->psiRZ);

  do {
    const char *fname = inp->psi_cubic_fname ? inp->psi_cubic_fname : "psi_cubic.gkyl";
    if (inp->write_psi_cubic)
      gkyl_dg_basis_ops_evalf_write_cubic(evcub, fname);
  } while (0);

  // construct grid in RZ plane
  
  enum { NPSI, NZ };
  int nc[2];
  nc[NPSI] = inp->comp_grid->cells[0]+1;
  nc[NZ] = inp->comp_grid->cells[2]+1;
  
  long nctot = nc[NPSI]*nc[NZ];
  geo->nodesrz = gkyl_array_new(GKYL_DOUBLE, 2, nctot);

  struct gkyl_range node_rng;
  gkyl_range_init_from_shape(&node_rng, 2, nc);

  double zlow = inp->comp_grid->lower[2], zup = inp->comp_grid->upper[2];
  double dz = (zup-zlow)/(nc[NZ]-1);

  double psi_lo = inp->comp_grid->lower[0];
  double psi_up = inp->comp_grid->upper[0];

  bool inc_axis = inp->include_axis;

  // adjust if we are using sqrt(psi) as radial coordinate
  double psic_lo = psi_lo, psic_up = psi_up;
  if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // when we include axis the psi_lo is ignored
    psic_lo = inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }
  
  double dpsi = (psic_up-psic_lo)/(nc[NPSI]-1);

  double rlow = lower[0], rup = upper[0];
  double rmin = rlow + 1e-6*(rup-rlow);

  struct psirz_ctx pctx = { .evcub = evcub };

  bool status = true;
  for (int iz=0; iz<nc[NZ]; ++iz) {
    double zcurr = zlow + iz*dz;

    double psi_min[1], psi_max[1];
    evcub->eval_cubic(0.0, (double[2]) { rmin, zcurr }, psi_min, evcub->ctx);
    evcub->eval_cubic(0.0, (double[2]) { rup, zcurr }, psi_max, evcub->ctx);

    for (int ipsi=0; ipsi<nc[NPSI]; ++ipsi) {

      if (inc_axis && (ipsi == 0)) {
        int idx[2] = { ipsi, iz };
        double *rz = gkyl_array_fetch(geo->nodesrz, gkyl_range_idx(&node_rng, idx));
        rz[0] = 0.0; rz[1] = zcurr;
      }
      else {
        double psic_curr = psic_lo + ipsi*dpsi;
        pctx.Z = zcurr;

        // we continue to do root-finding for psi and not sqrt(psi)
        double psi_curr = psic_curr;
        if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          psi_curr = psic_curr*psic_curr;
        
        pctx.psi =  psi_curr;
        
        struct gkyl_qr_res root = gkyl_ridders(psirz, &pctx, rmin, rup,
          psi_min[0]-psi_curr, psi_max[0]-psi_curr,
          100, 1e-10);

        if (root.status) {
          status = false;
          goto cleanup;
        }
        
        int idx[2] = { ipsi, iz };
        double *rz = gkyl_array_fetch(geo->nodesrz, gkyl_range_idx(&node_rng, idx));
        rz[0] = root.res; rz[1] = zcurr;
      }
    }
  }

  cleanup:

  if (true != status) {
    gkyl_mirror_grid_gen_release(geo);
    geo = 0;
    fprintf(stderr, "gkyl_mirror_grid_gen_inew failed to generate a grid\n");
  }
  
  gkyl_dg_basis_ops_evalf_release(evcub);
  
  return geo;
}

void
gkyl_mirror_grid_gen_release(struct gkyl_mirror_grid_gen *geom)
{
  gkyl_array_release(geom->nodesrz);
  gkyl_free(geom->gg_x);
  gkyl_free(geom);
}
