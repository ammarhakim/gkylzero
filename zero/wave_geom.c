#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>

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

struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid, struct gkyl_range *ext_range,
  evalf_t mapc2p, void *ctx)
{
  struct gkyl_wave_geom *wg = gkyl_malloc(sizeof(struct gkyl_wave_geom));
  wg->geom = gkyl_array_new(GKYL_USER, sizeof(struct wave_cell_geom), ext_range->volume);

  return wg;
}

void
gkyl_wave_geom_release(struct gkyl_wave_geom *wg)
{
  gkyl_array_release(wg->geom);
  gkyl_free(wg);
}


void
gkyl_wave_geom_set_idx(const gkyl_wave_geom *wg, const int *idx)
{
}
