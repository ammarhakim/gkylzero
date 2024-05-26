// Test updater which applies excision BCs.
// 
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_bc_excision.h>

// Allocate array (filled with zeros).
static struct gkyl_array *mkarr(bool use_gpu, long nc, long size) {
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost) {
  int ndim = parent->ndim;
  for (int d = 0; d < ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d], d,
      GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d], d,
      GKYL_UPPER_EDGE, parent, ghost);
  }
}

void evalFunc_1x1v(double t, const double *xn, double *restrict fout,
  void *ctx) {
  double x = xn[0], vx = xn[1];
  fout[0] = (x) * (vx - 0.5) * (vx - 0.5);
}

void evalFunc_1x2v(double t, const double *xn, double *restrict fout,
  void *ctx) {
  double x = xn[0], vx = xn[1], vy = xn[2];
  fout[0] = (x * x) * (vx - 0.5) * (vy - 0.5);
}

void evalFunc_2x2v(double t, const double *xn, double *restrict fout,
  void *ctx) {
  double x = xn[0], y = xn[1];
  double vx = xn[2], vy = xn[3];
  fout[0] = x * y * (vx - 1) * (vy - 2);
}

void test_bc(int cdim, int vdim, int poly_order, bool use_gpu) {
  const int ndim = cdim + vdim;
  double lower[ndim], upper[ndim];
  int cells[ndim];
  for (int i = 0; i < ndim; i++) {
    lower[i] = -2.0;
    upper[i] =  2.0;
    cells[i] = i < cdim ? 8 : 2;
  }
  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int i = 0; i < cdim; i++) {
    lower_conf[i] = lower[i];
    upper_conf[i] = upper[i];
    cells_conf[i] = cells[i];
  }

  // Grids
  struct gkyl_rect_grid grid, grid_conf;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  // Basis functions
  struct gkyl_basis basis, basis_conf;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost[ndim];
  for (int d = 0; d < ndim; d++) ghost[d] = 0;
  for (int d = 0; d < cdim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Projection updater for dist-function
  gkyl_proj_on_basis *projDistf;
  if (cdim == 1 && vdim == 1) {
    projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1,
      evalFunc_1x1v, NULL);
  } else if (cdim == 1 && vdim == 2) {
    projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1,
      evalFunc_1x2v, NULL);
  } else if (cdim == 2 && vdim == 2) {
    projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1,
      evalFunc_2x2v, NULL);
  }

  // Create distribution function array
  struct gkyl_array *distf, *distf_ho, *distf_init_ho;
  distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  distf_init_ho = mkarr(false, distf->ncomp, distf->size);

  // Project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local_ext, distf_ho);
  gkyl_array_copy(distf, distf_ho);
  gkyl_array_copy(distf_init_ho, distf_ho);

  // Determine the size of the BC buffer
  long buff_sz = 0;
  for (int d = 0; d < cdim; ++d) {
    long vol = skin_ghost.lower_ghost[d].volume;
    buff_sz  = buff_sz > vol ? buff_sz : vol;
  }
  struct gkyl_array *bc_buffer = mkarr(use_gpu, basis.num_basis, buff_sz);

  // Create and apply BC at the x-boundary with the direction tangential to the
  // hole being the second direction.
  int bc_dir = 0;
  int tan_dir = 1;
  struct gkyl_bc_excision *bclo = gkyl_bc_excision_new(tan_dir, grid, basis, skin_ghost.lower_ghost[bc_dir], use_gpu);

  // Copy ghost cells into buffer.
  gkyl_array_copy_to_buffer(bc_buffer->data, distf, &skin_ghost.lower_ghost[bc_dir]);
  gkyl_bc_excision_advance(bclo, bc_buffer, distf);

  gkyl_array_copy(distf_ho, distf);
  // Check the lower ghost cells after applying BC.
  int idx_neigh[ndim];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skin_ghost.lower_ghost[bc_dir]);
  printf("\n");
  while (gkyl_range_iter_next(&iter)) {
    for (int d=0; d<ndim; d++) idx_neigh[d] = iter.idx[d];
    idx_neigh[tan_dir] = iter.idx[tan_dir] > cells[tan_dir]/2 ? iter.idx[tan_dir]-cells[tan_dir]/2
                                                              : iter.idx[tan_dir]+cells[tan_dir]/2;

    int linidx = gkyl_range_idx(&skin_ghost.lower_ghost[bc_dir], iter.idx);
    int linidx_neigh = gkyl_range_idx(&skin_ghost.lower_ghost[bc_dir], idx_neigh);

    const double *f = gkyl_array_cfetch(distf_ho, linidx);
    const double *f_neigh = gkyl_array_cfetch(distf_init_ho, linidx_neigh);

    for (int i = 0; i < basis.num_basis; i++) {
      TEST_CHECK(gkyl_compare(f[i], f_neigh[i], 1e-12));
    }
  }

  gkyl_bc_excision_release(bclo);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf_init_ho);
  gkyl_array_release(bc_buffer);
  gkyl_proj_on_basis_release(projDistf);
}

void test_bc_excision_1x1v_p1_ho() { test_bc(1, 1, 1, false); }

#ifdef GKYL_HAVE_CUDA
void test_bc_excision_1x1v_p1_dev() { test_bc(1, 1, 1, true); }
#endif

TEST_LIST = {
  {"test_bc_excision_1x1v_p1_ho", test_bc_excision_1x1v_p1_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_bc_excision_1x1v_p1_dev", test_bc_excision_1x1v_p1_dev},
#endif
  {NULL, NULL},
};
