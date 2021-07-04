#include <assert.h>
#include <stdint.h>

#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>

void
gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells)
{
  grid->ndim = ndim;  
  grid->cellVolume = 1.0;
  for (int i=0; i<ndim; ++i) {
    grid->lower[i] = lower[i];
    grid->upper[i] = upper[i];
    grid->cells[i] = cells[i];
    grid->dx[i] = (upper[i]-lower[i])/cells[i];
    grid->cellVolume *= grid->dx[i];
  }
}

void
gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp)
{
  // dimension and shape are written as 64 bit integers
  uint64_t ndim = grid->ndim;
  uint64_t cells[GKYL_MAX_DIM];
  for (int d=0; d<grid->ndim; ++d)
    cells[d] = grid->cells[d];
  
  fwrite(&ndim, sizeof(uint64_t), 1, fp);
  fwrite(cells, sizeof(uint64_t), grid->ndim, fp);
  fwrite(grid->lower, sizeof(double), grid->ndim, fp);
  fwrite(grid->upper, sizeof(double), grid->ndim, fp);
}

bool
gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp)
{
  uint64_t ndim = grid->ndim;
  uint64_t cells64[GKYL_MAX_DIM];

  if (1 != fread(&ndim, sizeof(uint64_t), 1, fp))
    return false;
  if (ndim != fread(cells64, sizeof(uint64_t), ndim, fp))
    return false;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  if (ndim != fread(lower, sizeof(double), ndim, fp))
    return false;
  if (ndim != fread(upper, sizeof(double), ndim, fp))
    return false;

  // copy into regular int array
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) cells[d] = cells64[d];

  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);

  return true;
}
