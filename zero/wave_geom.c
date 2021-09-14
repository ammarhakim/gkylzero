#include <gkyl_array.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom_priv.h>

#include <math.h>

void
nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int i=0; i<GKYL_MAX_CDIM; ++i) xp[i] = xc[i];
}

// Geometry information for a single cell: recall a cell "owns" the
// faces on the lower side of cell
struct wave_cell_geom {
  double kappa; // ratio of cell-volume in phy to comp space
  double lenr[GKYL_MAX_CDIM]; // ratio of face-area in phys to comp space for "lower" faces
  double norm[GKYL_MAX_CDIM][GKYL_MAX_CDIM]; // norm[d] is the normal to face perp to direction 'd'
  // tau1[d] X tau2[d] = norm[d] are tangents to face perp to direction 'd'
  double tau1[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  double tau2[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
};

struct gkyl_wave_geom {
  struct gkyl_array *geom; // geometry in each cell
};

// Computes 1D geometry
void
calc_geom_1d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct wave_cell_geom *geo)
{
  double xlc[GKYL_MAX_CDIM], xrc[GKYL_MAX_CDIM];
  double xlp[GKYL_MAX_CDIM], xrp[GKYL_MAX_CDIM];

  xlc[0] = xc[0]-0.5*dx[0]; // left node
  xrc[0] = xc[0]+0.5*dx[0]; // right node

  // compute coordinates of left/right nodes
  mapc2p(0.0, xlc, xlp, ctx);
  mapc2p(0.0, xrc, xrp, ctx);

  geo->kappa = fabs(xrc-xlc)/dx[0];
  geo->lenr[0] = 1.0;

  geo->norm[0][0] = 1.0; geo->norm[0][1] = 0.0; geo->norm[0][2] = 0.0;
  geo->tau1[0][0] = 0.0; geo->tau1[0][1] = 1.0; geo->tau1[0][2] = 0.0;
  geo->tau2[0][0] = 0.0; geo->tau2[0][1] = 0.0; geo->tau2[0][2] = 1.0;
}

void
calc_geom_2d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct wave_cell_geom *geo)
{
}

void
calc_geom_3d(const double *dx, const double *xc, evalf_t mapc2p, void *ctx, struct wave_cell_geom *geo)
{
}

struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  evalf_t mapc2p, void *ctx)
{
  struct gkyl_wave_geom *wg = gkyl_malloc(sizeof(struct gkyl_wave_geom));
  wg->geom = gkyl_array_new(GKYL_USER, sizeof(struct wave_cell_geom), range->volume);

  double xc[GKYL_MAX_CDIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    struct wave_cell_geom *geo = gkyl_array_fetch(wg->geom, gkyl_range_idx(range, iter.idx));
    // compute geometry based on grid dimensions
    switch (grid->ndim) {
      case 1:
        calc_geom_1d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;

      case 2:
        calc_geom_2d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;

      case 3:
        calc_geom_3d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;
    };   
  }

  return wg;
}

void
gkyl_wave_geom_release(struct gkyl_wave_geom *wg)
{
  gkyl_array_release(wg->geom);
  gkyl_free(wg);
}


void
gkyl_wave_geom_set_idx(const struct gkyl_wave_geom *wg, const int *idx)
{
}
