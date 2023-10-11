#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_deflate_geo_priv.h>



gkyl_deflate_geo* gkyl_deflate_geo_new(const struct gkyl_basis *cbasis,const struct gkyl_basis *deflated_cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_rect_grid *deflated_grid, const int *rem_dirs, bool use_gpu){

  gkyl_deflate_geo *up = gkyl_malloc(sizeof(gkyl_deflate_geo));
  up->grid = grid;
  up->deflated_grid = deflated_grid;
  up->basis = cbasis;
  up->deflated_basis = deflated_cbasis;
  up->rem_dirs = gkyl_malloc(3*sizeof(int));
  for(int i=0; i<3; i++)
    up->rem_dirs[i] = rem_dirs[i];
  up->kernel = deflate_geo_choose_kernel(up->rem_dirs, up->deflated_grid->ndim, cbasis->b_type, cbasis->poly_order);

  printf("\n");
  for(int i = 0; i<3; i++){
    printf("in new remdir = %d\n", up->rem_dirs[i]);
  }

  return up;
}


void gkyl_deflate_geo_advance(const gkyl_deflate_geo *up, const struct gkyl_range *range, const struct gkyl_range* deflated_range, const struct gkyl_array *field, struct gkyl_array *deflated_field, int ncomp){
  // Inflated grid will always be 3 long in other directions
  // So deflate reange to the cell at index 2 in the ignored directions (z or z and y)
  // Then loop through and call the kernel on the inflated field and deflated field
  printf("in the advance\n");
  printf("grid ndim is = %d\n", up->grid->ndim);

  struct gkyl_range update_range;
  int loc_dir[3] = {0.};
  printf("allocated locdir\n");

  for(int i = 0; i<3; i++){
    loc_dir[i] = 1;
    if(up->rem_dirs[i] == 1){
      loc_dir[i] = 2;
    }
  }

  int do_idx[3];
  for(int i = 0; i<3; i++)
    do_idx[i] = loc_dir[i];

  //printf("range params:\n");
  //printf("lower[0] = %d, lower[1] = %d, lower[2] = %d \n", range->lower[0], range->lower[1], range->lower[2]);
  //printf("upper[0] = %d, upper[1] = %d, upper[2] = %d \n", range->upper[0], range->upper[1], range->upper[2]);

  //printf("loc_dir = %d %d %d\n", loc_dir[0], loc_dir[1], loc_dir[2]);
  //printf("up->rem_dirs = %d %d %d\n", up->rem_dirs[0], up->rem_dirs[1], up->rem_dirs[2]);
  //printf("made locdir\n");

  //printf("update range params:\n");
  //printf("lower[0] = %d, lower[1] = %d, lower[2] = %d \n", update_range.lower[0], update_range.lower[1], update_range.lower[2]);
  //printf("upper[0] = %d, upper[1] = %d, upper[2] = %d \n", update_range.upper[0], update_range.upper[1], update_range.upper[2]);

  //gkyl_range_deflate(&update_range, range, up->rem_dirs, loc_dir);

  //printf("update range params:\n");
  //printf("lower[0] = %d, lower[1] = %d, lower[2] = %d \n", update_range.lower[0], update_range.lower[1], update_range.lower[2]);
  //printf("upper[0] = %d, upper[1] = %d, upper[2] = %d \n", update_range.upper[0], update_range.upper[1], update_range.upper[2]);

  printf("deflated range params:\n");
  printf("lower[0] = %d, lower[1] = %d, lower[2] = %d \n", deflated_range->lower[0], deflated_range->lower[1], deflated_range->lower[2]);
  printf("upper[0] = %d, upper[1] = %d, upper[2] = %d \n", deflated_range->upper[0], deflated_range->upper[1], deflated_range->upper[2]);


  struct gkyl_range_iter iter;
  //gkyl_range_iter_init(&iter, &update_range);
  gkyl_range_iter_init(&iter, deflated_range);

  printf("before loop iter.idx = %d %d %d\n", iter.idx[0], iter.idx[1], iter.idx[2]);
  int idx_deflated[3];
  while(gkyl_range_iter_next(&iter)){
      printf("iter.idx = %d %d %d\n", iter.idx[0], iter.idx[1], iter.idx[2]);
      int count = 0;
      for(int i = 0; i<3; i++){
        if(up->rem_dirs[i]==0){
          do_idx[i] = iter.idx[count];
          //idx_deflated[count] = do_idx[i];
          count += 1;
        }
      }
      printf("do_idx = %d %d %d\n", do_idx[0], do_idx[1], do_idx[2]);

      //long loc = gkyl_range_idx(range, iter.idx);
      long loc = gkyl_range_idx(range, do_idx);
      const double *fld = gkyl_array_cfetch(field, loc);

      //printf("idx_deflated = %d %d %d\n", idx_deflated[0], idx_deflated[1], idx_deflated[2]);
      //long loc_deflated = gkyl_range_idx(deflated_range, idx_deflated);
      long loc_deflated = gkyl_range_idx(deflated_range, iter.idx);
      double *fld_deflated = gkyl_array_fetch(deflated_field, loc_deflated);
      for(int c = 0; c<ncomp; c++){
        up->kernel(&fld[c*up->basis->num_basis], &fld_deflated[c*up->deflated_basis->num_basis]);
      }


 }

}

void gkyl_deflate_geo_release(gkyl_deflate_geo* up){
  printf("trying to freee\n");
  gkyl_free(up->rem_dirs);
  printf("freed the int*\n");
  gkyl_free(up);

}


