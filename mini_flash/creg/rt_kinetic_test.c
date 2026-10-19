#include <gkyl_const.h>
#include <gkyl_flash.h>
#include <gkyl_flash_interface.h>
#include <gkyl_proj_on_basis.h>

// Context structure for initializing Gkeyll mini-FLASH app from FLASH.
struct gkyl_flash_ctx {
  int cdim, vdim;
  double mass;
  double nu; // Collision frequency
  double mu0, epsilon0; // Vacuum permeability and permittivity

  int num_cells_block[GKYL_MAX_CDIM], num_cells_vel[GKYL_MAX_VDIM];
  double lower_conf[GKYL_MAX_CDIM], upper_conf[GKYL_MAX_CDIM]; // [x, y, z] grid bounds
  double lower_vel[GKYL_MAX_VDIM], upper_vel[GKYL_MAX_VDIM]; // [vx, vy, vz] grid bounds

  double cfl; // Scaling factor to reduce CFL when computing kinetic dt.
};

/**
 * Initialization function for EM field, called at quadrature points during initialization.
 * Components of the output array are [Ex, Ey, Ez, Bx, By, Bz, phi, psi]
 * (Phi and psi are potentials used for divergence cleaning)
 * 
 * @param t Current time
 * @param xn Physical coordinates of grid point [x, y, z, vx, vy, vz]
 * @param fout Pointer to data in output gkyl_array
 * @param ctx Pointer to context structure with additional data.
 */
void init_field(double t, const double *xn, double *fout, void *ctx) 
{
  struct gkyl_flash_ctx *app = (struct gkyl_flash_ctx *)ctx;

  // Bx0 is *always* 1.0.
  double Bx0 = 1.0;

  // All other components are 0.
  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

  double Bx = Bx0; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).

  // Electric field.
  fout[0] = Ex; fout[1] = Ey; fout[2] = Ez;
  // Magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

/**
 * Initialization function for collision frequency.
 * 
 * @param t Current time
 * @param xn Physical coordinates of grid point [x, y, z, vx, vy, vz]
 * @param fout Pointer to data in output gkyl_array
 * @param ctx Pointer to context structure with additional data.
 */
void 
init_nu(double t, const double *xn, double *fout, void *ctx) 
{
  struct gkyl_flash_ctx *app = (struct gkyl_flash_ctx *)ctx;
  fout[0] = 100.0;
}

void 
initVel(double t, const double *xn, double *fout, void *ctx) 
{
  double ux, uy, uz = 0.0;
  fout[0] = 0.1;
  fout[1] = uy;
  fout[2] = uz;
}

void init_dist(double t, const double *xn, double *fout, void *ctx) {
  struct flash_sod_ctx *app = (struct flash_sod_ctx *)ctx;
  double x = xn[0], vx = xn[1];

  double n, vth;
  if (fabs(x) <= 0.0) {
    n = 1.0;
    vth = sqrt(1.0);
  }
  else {
    n = 0.125;
    vth = sqrt(0.125/0.1);
  }

  double F0 = n/sqrt(2.0*GKYL_PI*vth*vth) * exp(-(vx*vx)/(2.0*vth*vth));

  double G = vth * vth * F0;
  fout[0] = F0;
  fout[1] = G;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_flash_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
      iot->curr += 1;
    }

  gkyl_flash_app_write(app, t_curr, frame);
  }
}

int main(int argc, char **argv)
{
  double Lx = 1.0;
  double vthmax = 5.0;
  int NX = 8;
  int NV = 6;

  int cdim = 1;
  int vdim = 1;
  int pdim = cdim+vdim;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];
  double lower_conf[GKYL_MAX_DIM], upper_conf[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], cells_vel[GKYL_MAX_DIM], cells_conf[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];

  for (int i=0; i<cdim; ++i) {
    lower[i] = -Lx;
    upper[i] = Lx;
    cells[i] = NX;

    lower_conf[i] = -Lx;
    upper_conf[i] = Lx;
    cells_conf[i] = NX;
    ghost[i] = 1;
  }
  // cells[1] = 6;
  // cells_conf[1] = 6;

  for (int i=0; i<vdim; ++i) {
    lower[i+cdim] = -vthmax;
    upper[i+cdim] = vthmax;
    cells[i+cdim] = NV;
    ghost[i+cdim] = 0;

    lower_vel[i] = -vthmax;
    upper_vel[i] = vthmax;
    cells_vel[i] = NV;
  }

  struct gkyl_rect_grid grid, conf_grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&conf_grid, cdim, lower_conf, upper_conf, cells_conf);

  struct gkyl_range range, range_ext, conf_range, conf_range_ext;
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_basis basis, conf_basis;
  gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  gkyl_cart_modal_serendip(&conf_basis, cdim, 0);

  struct gkyl_range lower_skin, lower_ghost, upper_skin, upper_ghost;
  gkyl_skin_ghost_ranges(&lower_skin, &lower_ghost, 0, GKYL_LOWER_EDGE, &range_ext, ghost);
  gkyl_skin_ghost_ranges(&upper_skin, &upper_ghost, 0, GKYL_UPPER_EDGE, &range_ext, ghost);

  struct gkyl_flash_species spec = {
    .mass = 1.0,
    .name_kinetic = "kinetic",
    .name_fluid = "fluid",
    .cells = { NV },
    .lower = { -vthmax, -vthmax },
    .upper = { vthmax, vthmax },


    .collisions = {
      .num_cross_collisions = 0,
      .ctx = 0,
      .self_nu = init_nu,
    },

    .init_dist = init_dist,
    .init_dist_from_arrays = false,
    .init_fluid = initVel,
    .ctx = 0,

    .gas_gamma = 5.0/3.0,
    .limiter = GKYL_MIN_MOD,

    .bcx = { GKYL_SPECIES_COPY, GKYL_SPECIES_COPY },

    .write_pkpm_diags = true,
  };
  // snprintf(spec.name_kinetic, sizeof spec.name_kinetic, fmt_kinetic, (int)NX);
  // snprintf(spec.name_fluid, sizeof spec.name_fluid, fmt_fluid, (int)NX);

  struct gkyl_flash_field field = {
    .mu0 = 1.0, .epsilon0 = 1.0,
    .is_static = true,
    .ctx = 0,
    .init = init_field
  };

  struct gkyl_flash_inp inp = {
    .name = "test",
    .cdim = cdim, .vdim = vdim,
    .lower_conf = { -Lx, -Lx },
    .upper_conf = { Lx, Lx },

    .cells_fluid = { 2*NX, 2*NX },

    .cfl_frac = 0.9,

    .num_species = 1,
    .species = { spec },

    .field = field,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    // .using_flash = true,
  };

  gkyl_flash_app *app = gkyl_flash_app_new(&inp);
  gkyl_flash_app_apply_ic(app, 0.0);

  double dt = 1.0;
  double tcurr = 0.0;
  struct gkyl_update_status status;

  struct gkyl_tm_trigger io_trig = { .curr = 0, .dt = 0.1, .tcurr = tcurr };
  write_data(&io_trig, app, tcurr, false);

  struct gkyl_array *f_skin = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, range_ext.volume);
  gkyl_flash_app_copy_dist_skin(app, f_skin);

  // struct gkyl_range_iter iter;
  // gkyl_range_iter_init(&iter, &upper_skin);
  // while (gkyl_range_iter_next(&iter)) {
  //   long linp = gkyl_range_idx(&range, iter.idx);
  //   long linc = gkyl_range_idx(&conf_range, iter.idx);

  //   double *f_d = (double *)gkyl_array_fetch(f_skin, linp);

  //   printf("(%d, %d)\n", iter.idx[0], iter.idx[1]);
  //   printf("linc: %ld\n", linc);
  //   printf("%.3e\n", f_d[0]);
  //   printf("\n");
  // }

    gkyl_grid_sub_array_write(&grid, &range, 0, f_skin, "f_skin.gkyl");
}