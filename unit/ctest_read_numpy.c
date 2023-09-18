//#include <acutest.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>

#include <gkylzero.h>

struct gkyl_array *
array_from_numpy(FILE *fp, long sz)
{
  struct gkyl_array *arr
    = gkyl_array_new(GKYL_DOUBLE, 1, sz);
  double array[sz];
  //long res_sz = fread(array, 1, sizeof(double[sz]), fp);
  long res_sz = fread(arr->data, 1, sizeof(double[sz]), fp); //, sizeof(double[sz]), fp);

  if (res_sz != sizeof(double[sz])) {
    gkyl_array_release(arr);
    arr = 0;
  }
  /* for(int k=0;k<sz;k++) { */
  /*   printf("%d:  %g\n",k, array[k]); */
  /* } */
  return arr;
}

// 2d p=1
static void
nodal_to_modal(const double *f, double *mv)
{
  mv[0] = 0.5*(f[3]+f[2]+f[1]+f[0]);
  mv[1] = 0.2886751345948129*(f[3]-1.0*f[2]+f[1]-1.0*f[0]);
  mv[2] = 0.2886751345948129*(f[3]+f[2]-1.0*(f[1]+f[0]));
  mv[3] = 0.1666666666666667*(f[3]-1.0*(f[2]+f[1])+f[0]);
}

void
create_dg_from_nodal(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range_nodal,  const struct gkyl_array *adas_nodal,
  struct gkyl_array *adas_dg)
{
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, grid->cells);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  double nv[4];
  
  while (gkyl_range_iter_next(&iter)) {

    int ix = iter.idx[0], iy = iter.idx[1];
    int count = 0;
    for (int j=0; j<2; ++j) {
      for (int i=0; i<2; ++i) {
        long nidx = gkyl_range_idx(range_nodal, (const int[]) { ix+i, iy+j } );
        const double *adas_n = gkyl_array_cfetch(adas_nodal, nidx);
	//printf("%g\n",adas_n[0]);
	nv[count++] = adas_n[0];
      }
    }
    double *mv = gkyl_array_fetch(adas_dg, gkyl_range_idx(&range, iter.idx));
    nodal_to_modal(nv, mv);
  }
}

int file_isreg(const char *path) {
  struct stat st;
  if (stat(path, &st) < 0){
    return -1;
  }
  return S_ISREG(st.st_mode);
}

int
main(int argc, char **argv)
{
  int NT = 29, NN = 24;
  long sz = NT*NN;
  double logTmin = -0.69897, logTmax = 4., logNmin = 7.69897+6., logNmax = 15.30103+6.;

  //read nodal data from numpy file
  //if(file_isreg("ioniz_h_Z1.dat")==1) printf("found file\n");
  FILE *fp = fopen("ioniz_h_Z1.npy", "rb");
  
  struct gkyl_array *adas_nodal = array_from_numpy(fp, NT*NN);
  fclose(fp);

  if (!adas_nodal) {
    fprintf(stderr, "Unable to read data from nodal numpy file!\n");
    return 0;
  }

  struct gkyl_range range_node;
  gkyl_range_init_from_shape(&range_node, 2, (int[]) { NT, NN } );
  
  // allocate grid and DG array
  struct gkyl_rect_grid tn_grid;
  gkyl_rect_grid_init(&tn_grid, 2,
    (double[]) { logTmin, logNmin},
    (double []) { logTmax, logNmax},
    (int[]) { NT-1, NN-1 }
  );

  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, tn_grid.cells);
  
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1);
  
  struct gkyl_array *adas_dg =
    gkyl_array_new(GKYL_DOUBLE, basis.num_basis, NT*NN);

  create_dg_from_nodal(&tn_grid, &range_node, adas_nodal, adas_dg);

  gkyl_grid_sub_array_write(&tn_grid, &range, adas_dg, "adas_dg.gkyl");

  gkyl_array_release(adas_nodal);
  gkyl_array_release(adas_dg);
  
  return 0;
}

/* void */
/* test_create_gkyl_array() */
/* { */
/*   printf("-1"); */
/*   int NT = 29, NN = 24; */
/*   long sz = NT*NN; */
/*   double logTmin = -0.69897, logTmax = 4., logNmin = 7.69897+6., logNmax = 15.30103+6.; */

/*   printf("0"); */
/*   // read nodal data from numpy file */
/*   FILE *fp = fopen("ioniz_h_Z1.dat", "r"); */
/*   struct gkyl_array *adas_nodal = array_from_numpy(fp, NT*NN); */
/*   fclose(fp); */

/*   if (!adas_nodal) { */
/*     fprintf(stderr, "Unable to read data from nodal numpy file!\n"); */
/*   } */

/*   printf("1"); */
/*   struct gkyl_range range_node; */
/*   gkyl_range_init_from_shape(&range_node, 2, (int[]) { NT, NN } ); */
  
/*   // allocate grid and DG array */
/*   struct gkyl_rect_grid tn_grid; */
/*   gkyl_rect_grid_init(&tn_grid, 2, */
/*     (double[]) { logTmin, logNmin}, */
/*     (double []) { logTmax, logNmax}, */
/*     (int[]) { (NT-1)/2, (NN-1)/2 } */
/*   ); */

/*   printf("2"); */
/*   struct gkyl_range range; */
/*   gkyl_range_init_from_shape(&range, 2, tn_grid.cells); */
  
/*   struct gkyl_basis basis; */
/*   gkyl_cart_modal_serendip(&basis, 2, 2); */
  
/*   struct gkyl_array *adas_dg = */
/*     gkyl_array_new(GKYL_DOUBLE, basis.num_basis, (NT-1)/2*(NN-1)/2); */

/*   create_dg_from_nodal(&tn_grid, &range_node, adas_nodal, adas_dg); */

/*   gkyl_grid_sub_array_write(&tn_grid, &range, adas_dg, "adas_dg.gkyl"); */

/*   gkyl_array_release(adas_nodal); */
/*   gkyl_array_release(adas_dg); */
  
/* } */

/* void create_gkyl_array() {test_create_gkyl_array();} */

/* TEST_LIST = { */
/*   { "create_gkyl_array", create_gkyl_array }, */
/* }; */
