#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_wave_prop.h>
#include <gkyl_util.h>

struct gkyl_wave_prop {
    struct gkyl_rect_grid grid; // grid object
    int ndim; // number of dimensions
    int num_up_dirs; // number of update directions
    int update_dirs[GKYL_MAX_DIM]; // directions to update
    const struct gkyl_wv_eqn *equation; // equation object
};

void
gkyl_wave_prop_advance(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{

}

gkyl_wave_prop*
gkyl_wave_prop_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_wv_eqn *equation, int num_up_dirs, int update_dirs[])
{
  gkyl_wave_prop *up = gkyl_malloc(sizeof(gkyl_wave_prop));

  up->grid = *grid;
  up->ndim = grid->ndim;
  up->num_up_dirs = num_up_dirs;

  for (int i=0; i<num_up_dirs; ++i)
    up->update_dirs[i] = update_dirs[i];
  
  up->equation = gkyl_wv_eqn_aquire(equation);

  return up;
}

void
gkyl_wave_prop_release(gkyl_wave_prop* up)
{
  gkyl_wv_eqn_release(up->equation);
  free(up);
}
