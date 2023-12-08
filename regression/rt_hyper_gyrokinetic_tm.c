#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct kerntm_inp {
  int cdim, vdim, poly_order;
  int ccells[3], vcells[3];
  int nloop;
  bool use_gpu;
};

struct kerntm_inp
get_inp(int argc, char **argv)
{
  int c, cdim = 2, vdim = 2, poly_order = 2, nloop = 10;
  int nx, ny, nz = 8;
  int nvx, nvy, nvz = 16;
  bool use_gpu = false;
  while ((c = getopt(argc, argv, "+hgc:d:p:n:x:y:z:u:v:w:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: app_vlasov_kerntm -c CDIM -d VDIM -p POLYORDER -x NX -y NY -z NZ -u VX -v VY -w VZ -n NLOOP -g\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;        
      
      case 'c':
        cdim = atoi(optarg);
        break;

      case 'd':
        vdim = atoi(optarg);
        break;

      case 'p':
        poly_order = atoi(optarg);
        break;

      case 'n':
        nloop = atoi(optarg);
        break;          

      case 'x':
        nx = atoi(optarg);
        break;          

      case 'y':
        ny = atoi(optarg);
        break;          

      case 'z':
        nz = atoi(optarg);
        break;          

      case 'u':
        nvx = atoi(optarg);
        break;          

      case 'v':
        nvy = atoi(optarg);
        break;          

      case 'w':
        nvz = atoi(optarg);
        break;          

      case '?':
        break;
    }
  }
  
  return (struct kerntm_inp) {
    .cdim = cdim,
    .vdim = vdim,
    .poly_order = poly_order,
    .ccells = { nx, ny, nz },
    .vcells = { nvx, nvy, nvz},
    .nloop = nloop,
    .use_gpu = use_gpu,
  };
}

int
main(int argc, char **argv)
{
  struct kerntm_inp inp = get_inp(argc, argv);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  if (inp.use_gpu) {
    printf("Running kernel timers on GPU with:\n");
    use_gpu = true;
  } else {
    printf("Running kernel timers on CPU with:\n");
  }
#else
  printf("Running kernel timers on CPU with:\n");
#endif

  int cdim = inp.cdim;
  int vdim = inp.vdim;
  int poly_order = inp.poly_order;

  int cells[6];
  int ghost[6];
  double lower[6];
  double upper[6];
  int up_dirs[GKYL_MAX_DIM];
  int zero_flux_flags[GKYL_MAX_DIM];
  
  printf("cdim = %d; vdim = %d; poly_order = %d\n", inp.cdim, inp.vdim, inp.poly_order);
  printf("cells = [");
  for (int d=0; d<inp.cdim; ++d) {
    printf("%d ", inp.ccells[d]);
    cells[d] = inp.ccells[d];
    lower[d] = 0.;
    upper[d] = 1.;
    ghost[d] = 1;
    up_dirs[d] = d;
    zero_flux_flags[d] = 0;
  }
  for (int d=0; d<inp.vdim; ++d) {
    printf("%d ", inp.vcells[d]);
    cells[d+cdim] = inp.vcells[d];
    lower[d+cdim] = 0.;
    upper[d+cdim] = 1.;
    ghost[d+cdim] = 0;
    up_dirs[d+cdim] = d+cdim;
    zero_flux_flags[d+cdim] = 1;
  }
  printf("]\n");
    
  printf("nloop = %d\n", inp.nloop);
  
  // initialize grid and ranges
  int pdim = cdim+vdim;

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Initialize geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&confGrid, &confRange, &confRange_ext, &confBasis, 
    mapc2p, 0, bmag_func, 0, use_gpu);

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  double charge, mass = 1.;
  eqn = gkyl_dg_gyrokinetic_new(&confBasis, &basis, &confRange, charge, mass, gk_geom, use_gpu);

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate;
  struct gkyl_array *fin_h, *rhs_h;
  struct gkyl_array *phi, *apar, *apardot;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);

  phi = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  apar = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);
  apardot = mkarr1(use_gpu, confBasis.num_basis, confRange_ext.volume);

  gkyl_array_clear(phi, 1.);
  gkyl_array_clear(apar, 0.);
  gkyl_array_clear(apardot, 0.);

  struct gkyl_dg_gyrokinetic_auxfields aux = { .phi = phi, .apar = apar, .apardot = apardot };

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  // run hyper_dg_advance
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  struct timespec tm_start = gkyl_wall_clock();
  for(int n=0; n<inp.nloop; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_gyrokinetic_set_auxfields(eqn, aux);
    if (use_gpu) 
      gkyl_hyper_dg_advance_cu(slvr, &phaseRange, fin, cflrate, rhs);
    else
      gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs);
  }

#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  
  double tm_tot = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));
  printf("Avg time for gyrokinetic hyper dg: %g [s]\n", tm_tot/inp.nloop);
  
  return 0;
}
