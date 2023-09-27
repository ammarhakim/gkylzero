#include <stdlib.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

// Functions to extract ADAS data and project onto DG data
GKYL_CU_DH
static inline struct gkyl_array *
array_from_numpy(FILE *fp, long sz, int Zmax)
{
  struct gkyl_array *arr
    = gkyl_array_new(GKYL_DOUBLE, Zmax, sz);

  long res_sz = fread(arr->data, 1, sizeof(double[Zmax][sz]), fp); //, sizeof(double[sz]), fp);

  if (res_sz != sizeof(double[Zmax][sz])) {
    gkyl_array_release(arr);
    arr = 0;
  }
  return arr;
}

GKYL_CU_DH
static inline double * minmax_from_numpy(FILE *fp, long sz)
{
  double array[sz];
  long res_sz = fread(array, 1, sizeof(double[sz]), fp);
  double min = array[0];
  double max = array[sz-1];
  double *minmax = malloc(2);
  minmax[0] = min;
  minmax[1] = max;
  
  return minmax;
}

// 2d p=1
GKYL_CU_DH
static inline void
nodal_to_modal(const double *f, double *mv)
{
  mv[0] = 0.5*(f[3]+f[2]+f[1]+f[0]);
  mv[1] = 0.2886751345948129*(f[3]-1.0*f[2]+f[1]-1.0*f[0]);
  mv[2] = 0.2886751345948129*(f[3]+f[2]-1.0*(f[1]+f[0]));
  mv[3] = 0.1666666666666667*(f[3]-1.0*(f[2]+f[1])+f[0]);
}

/**
 * @param grid
 * @param range_nodal
 * @param adas_nodal
 * @param adas_dg
 * @param charge_state
 **/
GKYL_CU_DH
static inline void
create_dg_from_nodal(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range_nodal,  const struct gkyl_array *adas_nodal,
  struct gkyl_array *adas_dg, int charge_state)
{
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, grid->cells);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  double nv[4];
  int zi = charge_state - 1;
  
  while (gkyl_range_iter_next(&iter)) {

    int ix = iter.idx[0], iy = iter.idx[1];
    int count = 0;
    for (int j=0; j<2; ++j) {
      for (int i=0; i<2; ++i) {
        long nidx = gkyl_range_idx(range_nodal, (const int[]) { ix+i, iy+j } );
        const double *adas_n = gkyl_array_cfetch(adas_nodal, nidx);
	const double *adas_z_n = &adas_n[zi]; // get data for charge state
	nv[count++] = adas_n[0];
      }
    }
    double *mv = gkyl_array_fetch(adas_dg, gkyl_range_idx(&range, iter.idx));
    nodal_to_modal(nv, mv);
  }
}
