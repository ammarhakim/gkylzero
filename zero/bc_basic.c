#include <gkyl_bc_basic.h>
#include <gkyl_bc_basic_priv.h>
#include <gkyl_alloc.h>

struct gkyl_bc_basic*
gkyl_bc_basic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_range* local_range_ext,
  const int *num_ghosts, enum gkyl_bc_basic_type bctype, const struct gkyl_basis *basis,
  int cdim, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_basic *up = gkyl_malloc(sizeof(struct gkyl_bc_basic));

  up->dir = dir;
  up->edge = edge;

  // Create the skin/ghost ranges.
  gkyl_skin_ghost_ranges(&up->skin_r, &up->ghost_r, dir, edge,
                         local_range_ext, num_ghosts);

  // Create function applied to array contents (DG coefficients) when copying
  // to/from buffer.
  up->array_copy_func = gkyl_bc_basic_create_arr_copy_func(dir, cdim, bctype, basis);
  return up;
}

void
gkyl_bc_basic_advance(const struct gkyl_bc_basic *up, struct gkyl_array *buff_arr, struct gkyl_array *f_arr)
{
  // Apply BC in two steps:
  // 1. Copy skin to buffer while applying array_copy_func.
  gkyl_array_flip_copy_to_buffer_fn(buff_arr->data, f_arr, up->dir,
                                    up->skin_r, up->array_copy_func);
  // 2. Copy from buffer to ghost.
  gkyl_array_copy_from_buffer(f_arr, buff_arr->data, up->ghost_r);
}

void gkyl_bc_basic_release(struct gkyl_bc_basic *up)
{
  // Release memory associated with array_copy_func.
  if (gkyl_array_copy_func_is_cu_dev(up->array_copy_func)) {
    gkyl_cu_free(up->array_copy_func->ctx_on_dev);
    gkyl_cu_free(up->array_copy_func->on_dev);
  }
  gkyl_free(up->array_copy_func->ctx);
  gkyl_free(up->array_copy_func);
  // Release updater memory.
  gkyl_free(up);
}
