#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef struct adas_field {
  FILE *logData;
  FILE *logT;
  FILE *logN;
  long NT;
  long NN; 
  int Zmax;
  struct gkyl_array fld;
  double Eiz[GKYL_MAX_CHARGE_STATE];
} adas_field;

// Functions to extract ADAS data and project onto DG data
static inline void
array_from_numpy(FILE *fp, long sz, int Zmax, int charge_state, struct gkyl_array *arr)
{
  int zi = charge_state;
  double array[Zmax][sz];
  long res_sz = fread(array, 1, sizeof(double[Zmax][sz]), fp);
		   
  for (int i=0; i<sz; ++i) {
    double *arr_d = (double*) gkyl_array_fetch(arr, i);
    arr_d[0] = array[zi][i]; 
  }
}

static inline void
minmax_from_numpy(FILE *fp, long sz, double minmax[2])
{
  double array[sz];
  long res_sz = fread(array, 1, sizeof(double[sz]), fp);
  double min = array[0];
  double max = array[sz-1];
  minmax[0] = min;
  minmax[1] = max;
}

// 2d p=1
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
  int zi = charge_state;
  
  while (gkyl_range_iter_next(&iter)) {

    int ix = iter.idx[0], iy = iter.idx[1];
    int count = 0;
    for (int j=0; j<2; ++j) {
      for (int i=0; i<2; ++i) {
  	const int array_idx[2] = {ix+i, iy+j};
        long nidx = (long) gkyl_range_idx(range_nodal, array_idx );
        const double *adas_n = (const double*) gkyl_array_cfetch(adas_nodal, nidx);
  	//const double *adas_z_n = &adas_n[zi]; // get data for charge state
  	nv[count++] = adas_n[0];
      }
    }
    double *mv = (double*) gkyl_array_fetch(adas_dg, gkyl_range_idx(&range, iter.idx));
    nodal_to_modal(nv, mv);
  }
}

void
read_adas_field_iz(enum gkyl_ion_type type_ion, struct adas_field *data, const char *base);
 
void
read_adas_field_recomb(enum gkyl_ion_type type_ion, struct adas_field *data, const char *base);
