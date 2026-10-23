#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_flash_proj_pkpm_lte_on_basis.h>
#include <gkyl_util.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct flash_ctx {
  double mass;
  double n;
  double ux, uy, uz;
  double temp;

  double L;
  double vmax;
  int NC;
};

struct flash_ctx
init_ctx()
{
  double mass = 1.0;
  double n = 1.0;
  double ux = 0.5;
  double uy = 0.0;
  double uz = 0.0;
  double temp = 3.0;
  double vth = sqrt(temp/mass);
  double vmax = 5.0*vth;

  double L = 1.0;

  return (struct flash_ctx) {
    .mass = mass,
    .n = n,
    .ux = ux,
    .uy = uy,
    .uz = uz,
    .temp = temp,

    .L = L,
    .vmax = vmax,
    .NC = 16,
  };
}

void init_fluid(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct flash_ctx *app = ctx;
  fout[0] = app->n;
  fout[1] = app->ux;
  fout[2] = app->uy;
  fout[3] = app->uz;
  fout[4] = app->temp/app->mass;
}

void
test_proj_lte(int cdim, int vdim, int NV)
{
  struct flash_ctx ctx = init_ctx();
  int poly_order = 1; // Only ever using hybrid phase basis.

  int pdim = cdim+vdim;

  int cells[] = {ctx.NC, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double mass = 1.0;
  double L = 1.0;
  double vmax = 5.0;

  double lower[] = {0.0, -vmax, -vmax, -vmax};
  double upper[] = {L, vmax, vmax, vmax};

  double lower_vel[] = {-vmax, -vmax, -vmax};
  double upper_vel[] = {vmax, vmax, vmax};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis;
  gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);

  // Initialize arrays
  struct gkyl_array *fluidin, *distf;
  fluidin = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  distf = mkarr(2*phase_basis.num_basis, phase_range_ext.volume);

  // Project initial fluid moments
  int num_quad = poly_order+1;
  gkyl_proj_on_basis *proj_fluid = gkyl_proj_on_basis_new(&conf_grid,
    &conf_basis,num_quad, 5, init_fluid, &ctx);
  
  gkyl_proj_on_basis_advance(proj_fluid, 0.0, &conf_range, fluidin);
  gkyl_proj_on_basis_release(proj_fluid);

  // Initialize updater
  struct gkyl_flash_proj_pkpm_lte_on_basis_inp inp = {
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &conf_basis,
    .phase_basis = &phase_basis,
    .conf_range = &conf_range,
    .conf_range_ext = &conf_range_ext,
    .vel_range = &vel_range,
    .phase_range = &phase_range,
    .mass = mass,
    .use_gpu = false,
    .quad_type = 0
  };
  gkyl_flash_proj_pkpm_lte_on_basis *proj = gkyl_flash_proj_pkpm_lte_on_basis_new(&inp);

  // Project distribution function
  gkyl_flash_proj_pkpm_lte_on_basis_advance(proj, &phase_range, &conf_range, fluidin, distf);

  // Write out input moments and distribution functions
  const char *fmt_moms = "ctest_flash_proj_pkpm_lte_moms_%dx_p%d_%d.gkyl";
  const char *fmt_lte = "ctest_flash_proj_pkpm_lte_%dx%dv_p%d_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt_moms, cdim, poly_order, NV);
  char fileNm[sz+1]; // ensures no buffer overflow

  snprintf(fileNm, sizeof(fileNm), fmt_moms, cdim, poly_order, NV);
  gkyl_grid_sub_array_write(&conf_grid, &conf_range, 0, fluidin, fileNm);

  snprintf(fileNm, sizeof(fileNm), fmt_lte, cdim, vdim, poly_order, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, distf, fileNm);
}

void test_1x1v_p1_8() { test_proj_lte(1, 1, 8); }
void test_1x1v_p1_16() { test_proj_lte(1, 1, 16); }

TEST_LIST = {
  // { "test_1x1v_p1_8", test_1x1v_p1_8 },
  { "test_1x1v_p1_16", test_1x1v_p1_16 },

  { NULL, NULL }
};