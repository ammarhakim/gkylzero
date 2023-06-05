#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_gk_geom.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

struct gkyl_gkgeom {
  struct gkyl_rect_grid rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  int num_rzbasis; // number of basis functions in RZ

  
};

gkyl_gkgeom*
gkyl_gkgeom_new(const struct gkyl_gkgeom_inp *inp)
{
  struct gkyl_gkgeom *geo = gkyl_malloc(sizeof(*geo));

  geo->rzgrid = *inp->rzgrid;
  geo->psiRZ = gkyl_array_acquire(inp->psiRZ);
  geo->num_rzbasis = inp->rzbasis->num_basis;

  return geo;
}

void
gkyl_gkgeom_release(gkyl_gkgeom *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_free(geo);
}
