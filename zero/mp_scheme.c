#include "gkyl_wave_geom.h"
#include <float.h>

#include <gkyl_mp_scheme.h>
#include <gkyl_alloc.h>

struct gkyl_mp_scheme {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_mp_recon mp_recon; // base reconstruction to use
  
  double cfl; // CFL number
  const struct gkyl_wv_eqn *equation; // equation object
  struct gkyl_wave_geom *geom; // geometry object
};

gkyl_mp_scheme *
gkyl_mp_scheme_new(struct gkyl_mp_scheme_inp mpinp)
{
  struct gkyl_mp_scheme *mp = gkyl_malloc(sizeof(*mp));

  mp->grid = *(mpinp.grid);
  mp->ndim = mp->grid.ndim;
  
  mp->num_up_dirs = mpinp.num_up_dirs;
  for (int i=0; i<mpinp.num_up_dirs; ++i)
    mp->update_dirs[i] = mpinp.update_dirs[i];

  mp->mp_recon = mpinp.mp_recon;
  mp->equation = gkyl_wv_eqn_acquire(mpinp.equation);

  mp->geom = gkyl_wave_geom_acquire(mpinp.geom);

  return mp;
}

struct gkyl_mp_scheme_status
gkyl_mp_scheme_advance(gkyl_mp_scheme *mp,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout)
{

  return (struct gkyl_mp_scheme_status) {
    .success = true,
    .dt_suggested = DBL_MAX
  };
}

void
gkyl_mp_scheme_release(gkyl_mp_scheme* mp)
{
  gkyl_wv_eqn_release(mp->equation);
  gkyl_wave_geom_release(mp->geom);
  
  gkyl_free(mp);
}
