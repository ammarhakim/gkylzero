#include <stdio.h>
#include <math.h>

#include <gkylzero.h>

void evalFluidFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
  fout[1] = 0.0; 
  fout[2] = 0.0; 
  fout[3] = 0.0;
  fout[4] = 3.0/2.0;
}

void evalFieldFunc(double t, const double *xn, double* restrict emout, void *ctx)
{
  emout[0] = 1.0;
  emout[1] = 0.0;
  emout[2] = 0.0;
  emout[3] = 0.0;
  emout[4] = 0.0;
  emout[5] = 0.0;
  emout[6] = 0.0;
  emout[7] = 0.0;
}

void
gkyl_write_field(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_array *em, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = snprintf(0, 0, fmt, "app_sources_1d", frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, "app_sources_1d", frame);
  
  gkyl_grid_array_write(&grid, &range, em, fileNm);
}

void
gkyl_write_species(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_array *f, double tm, int frame)
{
  const char *fmt = "%s-euler_%d.gkyl";
  int sz = snprintf(0, 0, fmt, "app_sources_1d", frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, "app_sources_1d", frame);
  
  gkyl_grid_array_write(&grid, &range, f, fileNm);
}

int
main(void)
{
  int ndim = 1;
  double lower[] = { -1.0 }, upper[] = { 1.0 };
  int cells[] = { 1 };

  // create grid
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // compute local and local extended ranges
  int nghost[] = { 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  // allocate arrays
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, 5, ext_range.volume);
  struct gkyl_array *em = gkyl_array_new(GKYL_DOUBLE, 8, ext_range.volume);

  // initialize simulation
  gkyl_fv_proj *fv_fluid_proj = gkyl_fv_proj_new(&grid, 2, 5, evalFluidFunc, 0);
  gkyl_fv_proj_advance(fv_fluid_proj, 0.0, &range, f);

  gkyl_fv_proj *fv_em_proj = gkyl_fv_proj_new(&grid, 2, 8, evalFieldFunc, 0);
  gkyl_fv_proj_advance(fv_em_proj, 0.0, &range, em);

  gkyl_write_field(grid, range, em, 0.0, 0);
  gkyl_write_species(grid, range, f, 0.0, 0);

  // moment sources object
  struct gkyl_moment_em_coupling *mes =
    gkyl_moment_em_coupling_new( (struct gkyl_moment_em_coupling_inp) {
        .grid = &grid,
        .nfluids = 1,
        .param = {
          { .charge = -1.0, .mass = 1.0, .type = GKYL_ISO_EULER },
        },
        .epsilon0 = 1.0,
      }
    );

  // start, end and initial time-step
  double tcurr = 0.0, tend = 100.0;
  double dt = 0.1;
  int frame = 1;
  struct gkyl_array* fluids[] = {f};
  while (tcurr < tend) {
    printf("Taking time-step at t = %g ...", tcurr);
    gkyl_moment_em_coupling_advance(mes, dt, &range, fluids, em);
    printf(" dt = %g\n", dt); 
    tcurr += dt;
    gkyl_write_field(grid, range, em, tcurr, frame);
    gkyl_write_species(grid, range, f, tcurr, frame);
    frame += 1;
  }

  // free resources
  gkyl_fv_proj_release(fv_fluid_proj);
  gkyl_fv_proj_release(fv_em_proj);
  gkyl_array_release(f);
  gkyl_array_release(em);
  gkyl_moment_em_coupling_release(mes);
  
  return 0;
}
