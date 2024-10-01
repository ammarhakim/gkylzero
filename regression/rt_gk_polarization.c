#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_math.h>
#include <gkyl_eval_on_nodes.h>

#include <rt_arg_parse.h>

// Define the context of the simulation. This is basically all the globals
struct gk_polarization_ctx
{
  int cdim, vdim; // Dimensionality.
  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te;
  double n0;
  double beta;
  double tau;
  double Ti;
  double phi_center;
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_min;
  double psi_max;
  // Grid parameters
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int Nx;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
};

// Potential initial condition
void
eval_potential(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_polarization_ctx *app = ctx;
  double psi = xn[0]; // Magnetic flux function psi of field line.
  double z = xn[1];
  double z_min = app->z_min;
  double z_max = app->z_max;
  double psi_min = app->psi_min;
  double psi_max = app->psi_max;
  double z_norm = (z - z_min) / (z_max - z_min);
  double psi_norm = (psi - psi_min) / (psi_max - psi_min);

  double z_dependence = 4.*(pow(0.5,2) - pow(z_norm - 0.5, 2));
  double psi_dependence = 4.*(pow(0.5,2) - pow(psi_norm - 0.5, 2));
  fout[0] = app->phi_center * z_dependence * psi_dependence;
}

// Electrons initial conditions
void
eval_density_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_polarization_ctx *app = ctx;
  fout[0] = app->n0;
}

void
eval_upar_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0;
}

void
eval_temp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_polarization_ctx *app = ctx;
  fout[0] = app->Te;
}

// Ion initial conditions
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_polarization_ctx *app = ctx;
  fout[0] = app->n0;
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0;
}

void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_polarization_ctx *app = ctx;
  fout[0] = app->Ti;
}

struct gk_polarization_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double Ti = Te;

  // Thermal speeds.
  double vti = sqrt(Ti / mi);
  double vte = sqrt(Te / me);
  double c_s = sqrt(Te / mi);

  double phi_center = 1000;

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Geometry parameters.
  double z_min = -M_PI + 1e-1;
  double z_max = M_PI - 1e-1;
  double psi_min = 1e-3;
  double psi_max = 1e-2;

  // Grid parameters
  double vpar_max_elc = 4 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 4 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nx = 8;
  int Nz = 32;
  int Nvpar = 16; // Number of cells in the paralell velocity direction 96
  int Nmu = 16;  // Number of cells in the mu direction 192
  int poly_order = 1;

  struct gk_polarization_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mi = mi,
    .qi = qi,
    .me = me,
    .qe = qe,
    .Te = Te,
    .n0 = n0,
    .beta = beta,
    .Ti = Ti,
    .phi_center = phi_center,
    .z_min = z_min,
    .z_max = z_max,
    .psi_min = psi_min,
    .psi_max = psi_max,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_source_mom(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_polarization_ctx ctx = create_ctx();

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density_elc,
    .ctx_upar = &ctx,
    .upar= eval_upar_elc,
    .ctx_temp = &ctx,
    .temp = eval_temp_elc,
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .lower = {-ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .projection = elc_ic,
    .bcx = {
      .lower = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = elc_ic,
      },
      .upper = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = elc_ic,
      },
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .num_diag_moments = 1,
    .diag_moments = {"MaxwellianMoments"},
  };

  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density_ion,
    .ctx_upar = &ctx,
    .upar= eval_upar_ion,
    .ctx_temp = &ctx,
    .temp = eval_temp_ion,
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .projection = ion_ic,
    .bcx = {
      .lower = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = ion_ic,
      },
      .upper = {
        .type = GKYL_SPECIES_FIXED_FUNC,
        .projection = ion_ic,
      },
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .num_diag_moments = 1,
    .diag_moments = {"MaxwellianMoments"},
  };

  struct gkyl_gyrokinetic_field field =
  {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .poisson_bcs = {
      .lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_DIRICHLET},
      .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_DIRICHLET},
      .lo_value = {0.0, 0.0},
      .up_value = {0.0, 0.0},
    },
    .polarization_potential = eval_potential,
    .polarization_potential_ctx = &ctx,
  };

  struct gkyl_efit_inp efit_inp = {
    .filepath = "./data/eqdsk/wham.geqdsk", // equilibrium to use
    .rz_poly_order = 2,                     // polynomial order for psi(R,Z) used for field line tracing
    .flux_poly_order = 1,                   // polynomial order for fpol(psi)
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary 
  };

  // GK app
  struct gkyl_gk app_inp = {
    .name = "rt_gk_polarization",
    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = {ctx.psi_min, ctx.z_min},
    .upper = {ctx.psi_max, ctx.z_max},
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {0.0},
      .efit_info = efit_inp,
      .mirror_grid_info = grid_inp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 2,
    .species = {elc, ion},
    .field = field,
    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }  

  // Create triggers for IO.
  struct gkyl_tm_trigger trig_write = { .dt = 1, .tcurr = 0, .curr = 0 };

  // Write out ICs (if restart, it overwrites the restart frame).
  write_data(&trig_write, app, t_curr, false);

  struct gkyl_array *phi = gkyl_array_new(GKYL_DOUBLE, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_eval_on_nodes *evalDistf = gkyl_eval_on_nodes_new(&app->grid, &app->confBasis,1, eval_potential, &ctx);
  gkyl_eval_on_nodes_advance(evalDistf, 0.0, &app->local, phi);
  gkyl_eval_on_nodes_release(evalDistf);


    // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &app->local);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&app->local, conf_iter.idx);
    double *phi_d = gkyl_array_fetch(phi, linidx);\
    double *field_d  = gkyl_array_fetch(app->field->phi_host, linidx);
    // Ignore the corners
    if (conf_iter.idx[0] == 1){
      if(conf_iter.idx[1] == 1){
        continue;
      } else if (conf_iter.idx[1] == ctx.Nz){
        continue;
      }
    } else if (conf_iter.idx[0] == ctx.Nx){
      if(conf_iter.idx[1] == 1){
        continue;
      } else if (conf_iter.idx[1] == ctx.Nz){
        continue;
      }
    }
    for (int i=0; i<app->confBasis.num_basis; ++i) {
      assert( fabs(phi_d[i] - field_d[i])/fabs(phi_d[i]) < 3e-2 );
    }
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Passed \n", NULL);
  gkyl_array_release(phi);
  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
