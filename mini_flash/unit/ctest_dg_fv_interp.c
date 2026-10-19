#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_fv_interpolate.h>
#include <gkyl_fv_proj.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

enum FV_TO_MODAL_TYPE {
  QUAD, RECOV
};

struct interp_ctx {
  double ndim;
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Initialize test function to interpolate.
void
proj_fn(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  struct interp_ctx *app = (struct interp_ctx *)ctx;

  if (app->ndim == 1) {
    fout[0] = sin(M_PI*x);
  }
  else if (app->ndim == 2) {
    fout[0] = sin(M_PI*x)*cos(M_PI*y);
  }
  else if (app->ndim == 3) {
    fout[0] = sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*x+M_PI/3.0);
  }
  else {
    fout[0] = 0.0;
  }
}

void test_interp_modal_to_fv(int ndim, int num_cells_fv)
{

}

void test_interp_fv_to_modal(int ndim, int num_cells_fv, enum FV_TO_MODAL_TYPE type)
{
  struct interp_ctx ctx = {.ndim = ndim};

  // Initialize grids.
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM] = {0.0};
  int cells_dg[GKYL_MAX_DIM], cells_fv[GKYL_MAX_DIM] = {0};
  int ghost[GKYL_MAX_DIM] = {0};

  double L = 1.0;
  for (int d=0; d<ndim; ++d) {
    lower[d] = -L;
    upper[d] = L;
    cells_fv[d] = num_cells_fv;
    cells_dg[d] = (int)(num_cells_fv/2);
  }

  struct gkyl_rect_grid grid_dg, grid_fv;
  gkyl_rect_grid_init(&grid_dg, ndim, lower, upper, cells_dg);
  gkyl_rect_grid_init(&grid_fv, ndim, lower, upper, cells_fv);

  // Construct grid ranges.
  struct gkyl_range range_dg, range_dg_ext, range_fv, range_fv_ext;
  gkyl_create_grid_ranges(&grid_dg, ghost, &range_dg_ext, &range_dg);
  gkyl_create_grid_ranges(&grid_fv, ghost, &range_fv_ext, &range_fv);

  // Initialize basis, always p=1.
  struct gkyl_basis basis;
  int poly_order = 1;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Initialize arrays for FV and DG data.
  struct gkyl_array *fv_data = mkarr(1, range_fv.volume);
  struct gkyl_array *dg_data = mkarr(basis.num_basis, range_dg.volume);
  struct gkyl_array *dg_data_exact = mkarr(basis.num_basis, range_dg.volume);

  // Project initial condition onto finite volume grid.
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid_fv, 2, 1, proj_fn, &ctx);
  gkyl_fv_proj_advance(fv_proj, 0.0, &range_fv, fv_data);
  gkyl_fv_proj_release(fv_proj);

  // Project initial condition onto DG grid for comparison.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid_dg, &basis, 5, 1, proj_fn, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &range_dg, dg_data_exact);
  gkyl_proj_on_basis_release(proj); 
  
  // Interpolate finite volume data onto modal DG grid.
  struct gkyl_dg_fv_interpolate *interp = gkyl_dg_fv_interpolate_new(&basis, &grid_dg, &grid_fv);

  if (type == QUAD) {
    gkyl_dg_fv_interpolate_fv_to_modal_advance(interp, &range_dg, &range_fv, fv_data, 0, dg_data, 0);
  }
  else if (type == RECOV) {
    gkyl_dg_fv_interpolate_fv_to_modal_recov_advance(interp, &range_dg, &range_fv, fv_data, 0, dg_data, 0);
  }
  else {
    gkyl_exit("Unsupported interpolation type. Use enum QUAD or RECOV.");
  }

  // Compute difference beteween projection and interpolation.
  struct gkyl_array *dg_error = mkarr(basis.num_basis, range_dg.volume);
  gkyl_array_accumulate(dg_error, 1.0, dg_data_exact);
  gkyl_array_accumulate(dg_error, -1.0, dg_data);

  char name_base[128];
  if (type == QUAD) snprintf(name_base, 12,"fv_to_modal");
  if (type == RECOV) snprintf(name_base, 18,"fv_to_modal_recov");
  
  char fmt[128] = "%s_%dx_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, name_base, ndim, num_cells_fv);
  char fname[128];
  snprintf(fname, sz+1, fmt, name_base, ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_dg, &range_dg, 0, dg_data, fname);

  char fmt_exact[128] = "%s_exact_%dx_%d.gkyl";
  sz = gkyl_calc_strlen(fmt_exact, name_base, ndim, num_cells_fv);
  char fname_exact[128];
  snprintf(fname_exact, sz+1, fmt_exact, name_base, ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_dg, &range_dg, 0, dg_data_exact, fname_exact);

  char fmt_error[128] = "%s_error_%dx_%d.gkyl";
  sz = gkyl_calc_strlen(fmt_error, name_base, ndim, num_cells_fv);
  char fname_error[128];
  snprintf(fname_error, sz+1, fmt_error, name_base, ndim, num_cells_fv);
  gkyl_grid_sub_array_write(&grid_dg, &range_dg, 0, dg_error, fname_error);

  gkyl_array_release(fv_data);
  gkyl_array_release(dg_data);
  gkyl_array_release(dg_data_exact);
  gkyl_array_release(dg_error);
}

// NOTE: Cell number is for the refined grid.
void test_modal_to_fv_1x_p1_4() { test_interp_modal_to_fv(1, 4); }
void test_modal_to_fv_1x_p1_8() { test_interp_modal_to_fv(1, 8); }
void test_modal_to_fv_1x_p1_16() { test_interp_modal_to_fv(1, 16); }
void test_modal_to_fv_1x_p1_32() { test_interp_modal_to_fv(1, 32); }

void test_fv_to_modal_1x_p1_4() { test_interp_fv_to_modal(1, 4, QUAD); }
void test_fv_to_modal_1x_p1_8() { test_interp_fv_to_modal(1, 8, QUAD); }
void test_fv_to_modal_1x_p1_16() { test_interp_fv_to_modal(1, 16, QUAD); }
void test_fv_to_modal_1x_p1_32() { test_interp_fv_to_modal(1, 32, QUAD); }
void test_fv_to_modal_2x_p1_32() { test_interp_fv_to_modal(2, 32, QUAD); }
void test_fv_to_modal_3x_p1_32() { test_interp_fv_to_modal(3, 32, QUAD); }

void test_fv_to_modal_recov_1x_p1_4() { test_interp_fv_to_modal(1, 4, RECOV); }
void test_fv_to_modal_recov_1x_p1_8() { test_interp_fv_to_modal(1, 8, RECOV); }
void test_fv_to_modal_recov_1x_p1_16() { test_interp_fv_to_modal(1, 16, RECOV); }
void test_fv_to_modal_recov_1x_p1_32() { test_interp_fv_to_modal(1, 32, RECOV); }
void test_fv_to_modal_recov_2x_p1_32() { test_interp_fv_to_modal(2, 32, RECOV); }
void test_fv_to_modal_recov_3x_p1_32() { test_interp_fv_to_modal(3, 32, RECOV); }

TEST_LIST = {
  // { "test_modal_to_fv_1x_p1_4", test_modal_to_fv_1x_p1_4 },
  // { "test_modal_to_fv_1x_p1_8", test_modal_to_fv_1x_p1_8 },

  { "test_fv_to_modal_1x_p1_4", test_fv_to_modal_1x_p1_4 },
  { "test_fv_to_modal_1x_p1_8", test_fv_to_modal_1x_p1_8 },
  { "test_fv_to_modal_1x_p1_16", test_fv_to_modal_1x_p1_16 },
  { "test_fv_to_modal_1x_p1_32", test_fv_to_modal_1x_p1_32 },
  { "test_fv_to_modal_2x_p1_32", test_fv_to_modal_2x_p1_32 },
  { "test_fv_to_modal_3x_p1_32", test_fv_to_modal_3x_p1_32 },

  { "test_fv_to_modal_recov_1x_p1_4", test_fv_to_modal_recov_1x_p1_4 },
  { "test_fv_to_modal_recov_1x_p1_8", test_fv_to_modal_recov_1x_p1_8 },
  { "test_fv_to_modal_recov_1x_p1_16", test_fv_to_modal_recov_1x_p1_16 },
  { "test_fv_to_modal_recov_1x_p1_32", test_fv_to_modal_recov_1x_p1_32 },
  { "test_fv_to_modal_recov_2x_p1_32", test_fv_to_modal_recov_2x_p1_32 },
  { "test_fv_to_modal_recov_3x_p1_32", test_fv_to_modal_recov_3x_p1_32 },

  // { "test_fv_to_modal_recov_1x_p1_4", test_fv_to_modal_recov_1x_p1_4 },
  // { "test_fv_to_modal_recov_1x_p1_8", test_fv_to_modal_recov_1x_p1_8 },
  
  { NULL, NULL }
};