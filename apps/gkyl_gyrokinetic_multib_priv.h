// Private header for use in MB gyrokinetic app: do not include in
// user-facing header files!
#pragma once

#include <gkyl_gyrokinetic_priv.h>

// top-level internal App
struct gkyl_gyrokinetic_multib_app {
  int num_blocks; // total number of blocks
  struct gkyl_gyrokinetic_app **app; // individual App objects, one per-block

  struct gkyl_gyrokinetic_stat stat; // statistics
};
