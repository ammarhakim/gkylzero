#include <gkyl_bc_excision.h>
#include <gkyl_bc_excision_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_bc_excision*
gkyl_bc_excision_new(int tangential_dir, const struct gkyl_rect_grid grid,
  const struct gkyl_basis basis, const struct gkyl_range ghost_r, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_excision *up = gkyl_malloc(sizeof(*up));

  up->tan_dir = tangential_dir;
  up->tan_dir_num_cellsD2 = grid.cells[tangential_dir]/2;
  up->num_basis = basis.num_basis;
  up->ghost_r = ghost_r;
  up->use_gpu = use_gpu;

  gkyl_range_init(&up->buff_r, ghost_r.ndim, ghost_r.lower, ghost_r.upper);

  return up;
}

void
gkyl_bc_excision_advance(const struct gkyl_bc_excision *up,
  const struct gkyl_array *ghost_buffer, struct gkyl_array *distf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_excision_advance_cu(up, ghost_buffer, distf);
    return;
  }
#endif

  int sidx[GKYL_MAX_DIM]; // Shifted index.
  int tan_cellsD2 = up->tan_dir_num_cellsD2; // Number of cells in direction tangential to boundary / 2.

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->ghost_r);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_copy_int_arr(up->ghost_r.ndim, iter.idx, sidx);
    // Shift the index in the direction tangential to the boundary so that we place
    // f from a ghost cell in the ghost cell on the opposite side of the excision.
    int tan_idx = iter.idx[up->tan_dir];
    sidx[up->tan_dir] = tan_idx > tan_cellsD2 ? tan_idx-tan_cellsD2 : tan_idx+tan_cellsD2;

    long buff_loc = gkyl_range_idx(&up->buff_r, sidx);
    long out_loc = gkyl_range_idx(&up->ghost_r, iter.idx);

    const double *buff = (const double*) gkyl_array_cfetch(ghost_buffer, buff_loc);
    double *out = (double*) gkyl_array_fetch(distf, out_loc);

    for (int i=0; i<up->num_basis; i++)
      out[i] = buff[i];
  }
}

void gkyl_bc_excision_release(struct gkyl_bc_excision *up)
{
  // Release memory associated with this updater.
  gkyl_free(up);
}
