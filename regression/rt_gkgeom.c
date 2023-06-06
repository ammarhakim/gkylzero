#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gkgeom.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

// Cerfon equilibrium
struct cerfon_ctx {
  double R0, psi_prefactor;
};

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

void
psi_cerfon(double t, const double *xn, double *fout, void *ctx)
{
  struct cerfon_ctx *s = ctx;
  double R0 = s->R0, psi_prefactor = s->psi_prefactor;
  double R = xn[0], Z = xn[1];
  double x = R/R0, y = Z/R0;

  fout[0] = psi_prefactor*(0.00373804283369699*hex(x)*log(x) - 0.00574955335438162*hex(x) - 0.0448565140043639*qad(x)*sq(y)*log(x) + 0.0503044260840946*qad(x)*sq(y) + 0.017623348727471*qad(x)*log(x) + 0.0956643504553683*qad(x) + 0.0299043426695759*sq(x)*qad(y)*log(x) - 0.0160920841654771*sq(x)*qad(y) - 0.0704933949098842*sq(x)*sq(y)*log(x) + 0.0644725519961135*sq(x)*sq(y) - 7.00898484784405e-5*sq(x)*log(x) - 0.303766642191745*sq(x) - 0.00199362284463839*hex(y) + 0.0117488991516474*qad(y) + 7.00898484784405e-5*sq(y) + 0.0145368720253975);
}

void
cerforn_rt(void)
{
  struct cerfon_ctx ctx = {  .R0 = 2.5, .psi_prefactor = 1.0 };

  // create RZ grid
  double lower[] = { 0.01, -6.0 }, upper[] = { 6.0, 6.0 };
  int cells[] = { 64, 128 };

  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rzpoly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rzpoly_order);

  // allocate psiRZ array, initialize and write it to file
  struct gkyl_array *psiRZ = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid,
    &rzbasis, 1, &psi_cerfon, &ctx);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  gkyl_eval_on_nodes_release(eon);

  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, psiRZ, "cerforn_psi.gkyl");

  // construct geometry calculator
  struct gkyl_gkgeom_inp inp = {
    // psiRZ and related inputs
    .rzgrid = &rzgrid,
    .rzbasis = &rzbasis,
    .psiRZ = psiRZ,
    .rzlocal = &rzlocal

  };

  gkyl_gkgeom *geo = gkyl_gkgeom_new(&inp);

  gkyl_gkgeom_release(geo);
  gkyl_array_release(psiRZ);
}

int
main(int argc, char **argcv)
{
  cerforn_rt();
  
  return 0;
}
