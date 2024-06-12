#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

int
idx_to_inup_ker(const int dim, const int* num_cells, const int* idx)
{
  int out_idx = 0;

  for (int d = 0; d < dim; d++) {
    if (idx[d] == num_cells[d]) {
      out_idx += (int)(pow(2.0, d) + 0.5);
    }
  }

  return out_idx;
}

int
idx_to_inloup_ker(const int dim, const int* num_cells, const int* idx)
{
  int out_idx = 0;

  for (int d = 0; d < dim; d++) {
    if (idx[d] == 1) {
      out_idx = (2 * out_idx) + (int)(pow(3.0, d) + 0.5);
    }
    else if (idx[d] == num_cells[d]) {
      out_idx = (2 * out_idx) + (int)(pow(3.0, d) + 0.5) + 1;
    }
  }

  return out_idx;
}