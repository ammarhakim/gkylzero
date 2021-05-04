#include <gkyl_alloc.h>
#include <gkyl_kep_scheme.h>

struct gkyl_kep_scheme {
    struct gkyl_rect_grid grid; // grid object
    int ndim; // number of dimensions
    int num_up_dirs; // number of update directions
    int update_dirs[GKYL_MAX_DIM]; // directions to update
    const struct gkyl_wv_eqn *equation; // equation object
};

gkyl_kep_scheme*
gkyl_kep_scheme_new(struct gkyl_kep_scheme_inp inp)
{
  struct gkyl_kep_scheme *up;
  up = gkyl_malloc(sizeof(*up));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = inp.num_up_dirs;
  for (int i=0; i<inp.num_up_dirs; ++i)
    up->update_dirs[i] = inp.update_dirs[i];

  up->equation = gkyl_wv_eqn_aquire(inp.equation);

  return up;
}

void
gkyl_kep_scheme_advance(const gkyl_kep_scheme *kep, const struct gkyl_range *update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  
}

void
gkyl_kep_scheme_release(gkyl_kep_scheme* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_free(up);
}
