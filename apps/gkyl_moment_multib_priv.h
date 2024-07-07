// Private header for use in MB moment app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_moment_priv.h>

// top-level internal App
struct gkyl_moment_multib_app {
  int num_blocks; // total number of blocks
  struct gkyl_moment_app **app; // individual App objects, one per-block

  struct gkyl_moment_stat stat; // statistics
};
