#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_vlasov.h>

void
gkyl_bc_external_furman_pivi_advance(void *data, const struct gkyl_array *arr,
  int dir, struct gkyl_range range, struct gkyl_array_copy_func *cf);
