// Test calculation of Vlasov moments of a distribution function.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  fout[0] = (x*x)*(vx)*(vx);
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void test_makeDistFunc()
{
  int poly_order = 1;
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {4, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  
  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
 
  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  
  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local_ext, distf);

  // understanding check on how to use ranges
    /*int i,j;
    for (i=0;i<ndim;i++){
       for (j=0;j<ndim;j++){
         printf("\nlower_ghost[%i].lower[%i] = %i \nlower_ghost[%i].upper[%i] = %i",
             i,j,skin_ghost.lower_ghost[i].lower[j],
             i,j,skin_ghost.lower_ghost[i].upper[j]);
       }
    }

    for (i=0;i<ndim;i++){
       for (j=0;j<ndim;j++){
         printf("\nlower_skin[%i].lower[%i] = %i \nlower_skin[%i].upper[%i] = %i",
             i,j,skin_ghost.lower_skin[i].lower[j],
             i,j,skin_ghost.lower_skin[i].upper[j]);
       }
     }*/
   /* int i;
  for (i = 0; i<12; i++){
    double *dfl = gkyl_array_fetch(distf,i);
    printf("\n Cell %i\n %6.5f\n%6.5f\n",i, dfl[0],dfl[1]);
  }*/

   // Go from 2D index to 1D ranges
   int lgl_idx = gkyl_range_idx(skin_ghost.lower_ghost,skin_ghost.lower_ghost[0].lower);
   int lgu_idx = gkyl_range_idx(skin_ghost.lower_ghost,skin_ghost.lower_ghost[0].upper);
   int lsl_idx = gkyl_range_idx(skin_ghost.lower_skin,skin_ghost.lower_skin[0].lower);
   int lsu_idx = gkyl_range_idx(skin_ghost.lower_skin,skin_ghost.lower_skin[0].upper);
   int usl_idx = gkyl_range_idx(skin_ghost.upper_skin,skin_ghost.upper_skin[0].lower);
   int usu_idx = gkyl_range_idx(skin_ghost.upper_skin,skin_ghost.upper_skin[0].upper);
   int ugl_idx = gkyl_range_idx(skin_ghost.upper_ghost,skin_ghost.upper_ghost[0].lower);
   int ugu_idx = gkyl_range_idx(skin_ghost.upper_ghost,skin_ghost.upper_ghost[0].upper);
    
   int i,j;
   for (i=lgl_idx;i<=lgu_idx;i++){
      double *val = gkyl_array_fetch(distf,i);
      TEST_CHECK(gkyl_compare(16.88888889,val[0],1e-6));
      for (j=0;j<9;j++){
        printf("%10.8f\n",val[j]);
      }
   }
   for (i=ugl_idx;i<=ugu_idx;i++){
      double *val = gkyl_array_fetch(distf,i);
      TEST_CHECK(gkyl_compare(16.88888889,val[0],1e-12));
      for (j=0;j<9;j++){
        printf("%10.8f\n",val[j]);
      }
   }
  // Make a dg equation object
  struct gkyl_dg_eqn* eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &local, GKYL_FIELD_E_B, false);
  
  //Apply boundary condition
  struct gkyl_array_copy_func* bc =  gkyl_vlasov_wall_bc_create(eqn, 0, &basis);
  // struct gkyl_array_copy_func* gkyl_vlasov_absorb_bc_create(eqn, int dir,pbasis)
  long buff_sz = 0;
  for (int d=0; d<cdim; ++d) {
    long vol = skin_ghost.lower_skin[d].volume;
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  struct gkyl_array *bc_buffer;
  bc_buffer = mkarr(basis.num_basis, buff_sz);
  gkyl_array_flip_copy_to_buffer_fn(&bc_buffer, distf, 0, skin_ghost.lower_skin[0], bc);
  gkyl_array_copy_from_buffer(distf, &bc_buffer, skin_ghost.lower_ghost[0]);




  // Check cells after applying boundary condition
   for (i=lgl_idx;i<=lgu_idx;i++){
      double *val = gkyl_array_fetch(distf,i);
      TEST_CHECK(gkyl_compare(2.,val[0],1e-12));
      for (j=0;j<9;j++){
        printf("%10.8f\n",val[j]);
      }
   }
   for (i=lsl_idx;i<=lsu_idx;i++){
      double *val = gkyl_array_fetch(distf,i);
      TEST_CHECK(gkyl_compare(2.,val[0],1e-12));
   }
   for (i=usl_idx;i<=usu_idx;i++){
      double *val = gkyl_array_fetch(distf,i);
      TEST_CHECK(gkyl_compare(2.,val[0],1e-12));
   }
   for (i=ugl_idx;i<=ugu_idx;i++){
      double *val = gkyl_array_fetch(distf,i);
      TEST_CHECK(gkyl_compare(2.,val[0],1e-12));
   }


  // release memory for moment data object
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_bc_release(bc);
}

TEST_LIST = {
  { "makeDistFunc", test_makeDistFunc },
  { NULL, NULL },
};
