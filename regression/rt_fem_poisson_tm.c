#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_util.h>

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

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

struct kerntm_inp {
  int dim, poly_order;
  int cells[3];
  int nloop;
  bool use_gpu;
};

struct kerntm_inp
get_inp(int argc, char **argv)
{
  int c, dim = 2, poly_order = 2, nloop = 10;
  int nx, ny, nz = 32;
  bool use_gpu = false;
  while ((c = getopt(argc, argv, "+hgd:p:n:x:y:z:")) != -1) {
    switch (c)
    {
      case 'h':
        printf("Usage: app_vlasov_kerntm -d DIM -p POLYORDER -x NX -y NY -z NZ -n NLOOP -g\n");
        exit(-1);
        break;

      case 'g':
        use_gpu = true;
        break;        
      
      case 'd':
        dim = atoi(optarg);
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

      case '?':
        break;
    }
  }
  
  return (struct kerntm_inp) {
    .dim = dim,
    .poly_order = poly_order,
    .cells = { nx, ny, nz },
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

  int dim = inp.dim;
  int poly_order = inp.poly_order;

  int cells[6];
  int ghost[6];
  double lower[6];
  double upper[6];
  int up_dirs[GKYL_MAX_DIM];
  int zero_flux_flags[GKYL_MAX_DIM];
  
  printf("dim = %d; poly_order = %d\n", inp.dim, inp.poly_order);
  printf("cells = [");
  for (int d=0; d<inp.dim; ++d) {
    printf("%d ", inp.cells[d]);
    cells[d] = inp.cells[d];
    lower[d] = 0.;
    upper[d] = 1.;
    ghost[d] = 1;
  }
  printf("]\n");
    
  printf("nloop = %d\n", inp.nloop);
  
  // initialize grid and ranges

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  // initialize basis
  struct gkyl_basis basis;

  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  // initialize arrays

  // create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr1(use_gpu, basis.num_basis, localRange_ext.volume);
  // create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr1(use_gpu, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *rho_h;

  // create epsilon array
  struct gkyl_array *epsilon = mkarr1(use_gpu, basis.num_basis, localRange_ext.volume);
  gkyl_array_shiftc(epsilon, 1.0*pow(sqrt(2.),dim), 0);
  
  // set initial condition
  int nf = localRange_ext.volume*basis.num_basis;
  double *rho_d;
  if (use_gpu) {
    rho_h = mkarr1(false, basis.num_basis, localRange_ext.volume);
    rho_d = rho_h->data;
  } else {
    rho_d = rho->data;
  }
  for(int i=0; i< nf; i++) {
    rho_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(rho, rho_h);

  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.up_type[0] = GKYL_POISSON_PERIODIC;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;

  // FEM poisson solver.
  gkyl_fem_poisson *poisson = gkyl_fem_poisson_new(&localRange, &grid, basis, &bc_tv, epsilon, NULL, true, use_gpu);

  int nrep = inp.nloop;

  // TIMING SET_RHS
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  struct timespec tm_start = gkyl_wall_clock();
  for(int n=0; n<nrep; n++) {
    // Set the RHS source.
    gkyl_fem_poisson_set_rhs(poisson, rho, NULL);
  }
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  
  double tm_tot = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));
  printf("Avg time for fem_poisson_set_rhs: %g [s]\n", tm_tot/inp.nloop);

  // TIMING SOLVE
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  tm_start = gkyl_wall_clock();
  for(int n=0; n<nrep; n++) {
    // Solve the problem.
    gkyl_fem_poisson_solve(poisson, phi);
  }
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  
  tm_tot = gkyl_time_sec(gkyl_time_diff(tm_start, gkyl_wall_clock()));
  printf("Avg time for fem_poisson_solve: %g [s]\n", tm_tot/inp.nloop);
  
  gkyl_array_release(rho);
  gkyl_array_release(phi);
  gkyl_array_release(epsilon);
  if (use_gpu)
    gkyl_array_release(rho);
  gkyl_fem_poisson_release(poisson);
  
  return 0;
}
