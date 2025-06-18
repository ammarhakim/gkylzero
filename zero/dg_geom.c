#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_dg_geom.h>

static bool
dg_geom_is_cu_dev(const struct gkyl_dg_geom* dgg)
{
  return GKYL_IS_CU_ALLOC(dgg->flags);
}

static void
dg_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_geom *dgg = container_of(ref, struct gkyl_dg_geom, ref_count);

  for (int d=0; d<dgg->range.ndim; ++d)
    gkyl_array_release(dgg->surf_geom[d]);
  
  if (dg_geom_is_cu_dev(dgg)) 
    gkyl_cu_free(dgg->on_dev); 

  gkyl_free(dgg);
}


struct gkyl_dg_geom *
gkyl_dg_geom_new(const struct gkyl_dg_geom_inp *inp)
{
  struct gkyl_dg_geom *dgg = gkyl_malloc(sizeof *dgg);

  dgg->range = *inp->range;

  int ndim = dgg->range.ndim;
  int lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; ++d) {
    lower[d] = 0;
    upper[d] = inp->nquad-1; // inclusive
  }
  // NOTE: surfaces are ndim-1 objects
  gkyl_range_init(&dgg->surf_quad_range, ndim-1, lower, upper);
  gkyl_range_init(&dgg->vol_quad_range, ndim, lower, upper);

  for (int d=0; d<ndim; ++d)
    dgg->surf_geom[d] = gkyl_array_new(GKYL_USER,
      sizeof(struct gkyl_dg_surf_geom[dgg->surf_quad_range.volume]), dgg->range.volume);

  dgg->flags = 0;
  GKYL_CLEAR_CU_ALLOC(dgg->flags);
  dgg->ref_count = gkyl_ref_count_init(dg_geom_free);
  dgg->on_dev = dgg; // CPU eqn obj points to itself
  
  return dgg;
}


struct gkyl_dg_geom*
gkyl_dg_geom_acquire(const struct gkyl_dg_geom* dgg)
{
  gkyl_ref_count_inc(&dgg->ref_count);
  return (struct gkyl_dg_geom*) dgg;
}

void
gkyl_dg_geom_release(const struct gkyl_dg_geom *dgg)
{
  gkyl_ref_count_dec(&dgg->ref_count);
}
