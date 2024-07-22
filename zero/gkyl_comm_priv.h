#pragma once

#include <gkyl_range.h>

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  long max_vol; // maximum vol of send/recv region
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
#define G_MAX(a,b) (a)>(b)?(a):(b)
  
  int ndim = parent->ndim;
  long max_vol = 0;
    
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);
    
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

// Create ghost and skin sub-ranges given a parent range: includes
// corners
static void
skin_ghost_ranges_with_corners_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
#define G_MAX(a,b) (a)>(b)?(a):(b)
  
  int ndim = parent->ndim;
  long max_vol = 0;
    
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_with_corners_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);
    
    gkyl_skin_ghost_with_corners_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

// define long -> skin_ghost_ranges ...
#define i_key long
#define i_val struct skin_ghost_ranges
#define i_tag l2sgr
#include <stc/cmap.h>
// ... done with map definition
