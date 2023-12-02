#include <gkyl_gk_geometry.h>

struct gkyl_gk_geometry{
  // stuff for mapc2p and finite differences array
  struct gkyl_range* nrange;
  const struct gkyl_range* range;
  const struct gkyl_range* range_ext;
  const struct gkyl_basis* basis;
  const struct gkyl_rect_grid* grid;
  double* dzc;
  evalf_t mapc2p_func;
  evalf_t bmag_func;

  struct gkyl_array* mc2p;
  struct gkyl_array* mc2p_nodal;
  struct gkyl_array* mc2p_nodal_fd;

  struct gkyl_array* bmag; // bmag
  struct gkyl_array* g_ij;
  struct gkyl_array* jacobgeo;
  struct gkyl_array* jacobgeo_inv;
  struct gkyl_array* gij;
  struct gkyl_array* b_i;
  struct gkyl_array* cmag;
  struct gkyl_array* jacobtot;
  struct gkyl_array* jacobtot_inv;
  struct gkyl_array* bmag_inv;
  struct gkyl_array* bmag_inv_sq;
  struct gkyl_array* gxxj;
  struct gkyl_array* gxyj;
  struct gkyl_array* gyyj;
};
