#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_bc_excision gkyl_bc_excision;

/**
 * Create a new updater that applies excision BCs. These are BCs that
 * instantaneously transport the field tangential to the boundary direction,
 * like falling through a hole and appearing on the other side.
 *
 * @param tangential_dir Direction tangential to the boundary along which to transport the field.
 * @param grid Grid the field is defined on.
 * @param basis Basis the field is represented with.
 * @param ghost_r Ghost range to fill with BCs.
 * @param use_gpu Whether to use the GPU or not.
 * @return New updater pointer.
 */
struct gkyl_bc_excision*
gkyl_bc_excision_new(int tangential_dir, const struct gkyl_rect_grid grid,
  const struct gkyl_basis basis, const struct gkyl_range ghost_r, bool use_gpu);

/**
 * Apply the excision BC.
 *
 * @param up Excision BC updater.
 * @param ghost_buffer Buffer with the values of the corresponding ghost cells.
 * @param distf Distribution function to apply BC to.
 */
void
gkyl_bc_excision_advance(const struct gkyl_bc_excision *up, const struct gkyl_array *ghost_buffer,
  struct gkyl_array *distf);

/**
 * Free memory associated with excision BC updater.
 *
 * @param up Excision BC updater.
 */
void gkyl_bc_excision_release(struct gkyl_bc_excision *up);
