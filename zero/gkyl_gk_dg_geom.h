#include <gkyl_gk_geometry.h>
#include <gkyl_dg_geom.h>

void gkyl_gk_dg_geom_populate_vol(struct gkyl_dg_geom *dg_geom, struct gk_geometry* gk_geom);

void gkyl_gk_dg_geom_populate_surf(struct gkyl_dg_geom *dg_geom, struct gk_geometry* gk_geom);

void gkyl_gk_dg_geom_write_vol(struct gkyl_dg_geom *dg_geom, struct gk_geometry* gk_geom, const char *name);

void gkyl_gk_dg_geom_write_surf(struct gkyl_dg_geom *dg_geom, struct gk_geometry* gk_geom, const char *name);
