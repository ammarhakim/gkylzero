// Private header for common functions between species. Do not include
// in any public facing header!
#pragma once

#include "gkyl_gyrokinetic.h"
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_comm.h>

#include <stdio.h>
#include <stdlib.h>

// ranges for use in BCs
struct app_skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// allocate double array (filled with zeros)
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
// allocate integer array (filled with zeros)
static struct gkyl_array*
mk_int_arr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_INT, nc, size);
  else
    a = gkyl_array_new(GKYL_INT, nc, size);
  return a;
}

// Compute out = c1*arr1 + c2*arr2
static inline struct gkyl_array*
array_combine(struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range *rng)
{
  return gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct app_skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}


static struct gkyl_rect_grid augment_grid(struct gkyl_rect_grid grid, struct gkyl_gyrokinetic_geometry geometry)
{
  struct gkyl_rect_grid augmented_grid;
  int cells[3];
  double lower[3];
  double upper[3];

  if (grid.ndim==1) {
    cells[0] = 1;
    cells[1] = 1;
    cells[2] = grid.cells[0];

    lower[0] = geometry.world[0] - 1e-5;
    lower[1] = geometry.world[1] - 1e-1;
    lower[2] = grid.lower[0];

    upper[0] = geometry.world[0] + 1e-5;
    upper[1] = geometry.world[1] + 1e-1;
    upper[2] = grid.upper[0];
  }
  else if (grid.ndim==2) {
    cells[0] = grid.cells[0];
    cells[1] = 1;
    cells[2] = grid.cells[1];

    lower[0] = grid.lower[0];
    lower[1] = geometry.world[0] - 1e-1;
    lower[2] = grid.lower[1];

    upper[0] = grid.upper[0];
    upper[1] = geometry.world[0] + 1e-1;
    upper[2] = grid.upper[1];
  }

  gkyl_rect_grid_init(&augmented_grid, 3, lower, upper, cells);
  return augmented_grid;
}


static void augment_local(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  if (inrange->ndim == 2) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = inrange->lower[0]-nghost[0];
    upper_ext[0] = inrange->upper[0]+nghost[0];
    lower[0] = inrange->lower[0];
    upper[0] = inrange->upper[0];

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[1]-nghost[1];
    upper_ext[2] = inrange->upper[1]+nghost[1];
    lower[2] = inrange->lower[1];
    upper[2] = inrange->upper[1];


    gkyl_range_init(ext_range, inrange->ndim+1, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
  else if (inrange->ndim == 1) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = 1 - 1;
    upper_ext[0] = 1 + 1;
    lower[0] = 1;
    upper[0] = 1;

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[0]-nghost[0];
    upper_ext[2] = inrange->upper[0]+nghost[0];
    lower[2] = inrange->lower[0];
    upper[2] = inrange->upper[0];



    gkyl_range_init(ext_range, inrange->ndim+2, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
}
