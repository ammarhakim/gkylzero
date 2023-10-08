#include <gkyl_efit.h>
#include <gkyl_rect_grid.h>
#include <assert.h>

struct gkyl_efit{
  const struct gkyl_basis *rzbasis;
  struct gkyl_rect_grid *rzgrid;
  bool use_gpu;
  FILE *fp;
  int nr, nz;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;
  struct gkyl_array *psizr;
  struct gkyl_array *psibyrzr;
  struct gkyl_array *psibyr2zr;
  struct gkyl_range *rzlocal;
  struct gkyl_range *rzlocal_ext;
  double rmin, rmax, zmin, zmax;
};

