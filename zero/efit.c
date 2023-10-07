#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <gkyl_efit.h>
#include <gkyl_efit_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include<gkyl_rect_decomp.h>

#include <gkyl_array.h>
#include <gkyl_range.h>

void nodal_array_to_modal_psi(const struct gkyl_array *nodal_array, struct gkyl_array *modal_array, const struct gkyl_range *update_range, const struct gkyl_range *nrange, const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, int num_ret_vals){
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_basis = basis->num_basis;
  int poly_order = basis->poly_order;
  //initialize the nodes
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, basis->num_basis);
  basis->node_list(gkyl_array_fetch(nodes, 0));
  double fnodal[num_basis]; // to store nodal function values

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  int nidx[3];
  long lin_nidx[num_basis];
  
  while (gkyl_range_iter_next(&iter)) {
     gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {
      const double* temp  = gkyl_array_cfetch(nodes,i);
      for( int j = 0; j < grid->ndim; j++){
        if(poly_order==1){
          //if(j==1 && ginp->bcs[1]==1) // this conversion is valid if ghost cells are included
          //  nidx[j] = iter.idx[j] + (temp[j]+1)/2 ;
          //else // otherwise this conversion is valid
            nidx[j] = iter.idx[j]-1 + (temp[j]+1)/2 ;
        }
        if (poly_order==2)
          nidx[j] = 2*(iter.idx[j]-1) + (temp[j] + 1) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(modal_array, lidx); // pointer to expansion in cell
    double fao[num_basis*num_ret_vals];
  
    for (int i=0; i<num_basis; ++i) {
      const double* temp = gkyl_array_cfetch(nodal_array, lin_nidx[i]);
      for (int j=0; j<num_ret_vals; ++j) {
        fao[i*num_ret_vals + j] = temp[j];
      }
    }

    for (int i=0; i<num_ret_vals; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fnodal[k] = fao[num_ret_vals*k+i];
      // transform to modal expansion
      basis->nodal_to_modal(fnodal, &arr_p[num_basis*i]);
    }
  }

}



gkyl_efit* gkyl_efit_new(char *filepath, const struct gkyl_basis *rzbasis,
  struct gkyl_rect_grid *rzgrid, struct gkyl_range *rzlocal, struct gkyl_range *rzlocal_ext, bool use_gpu)
{

  gkyl_efit *up = gkyl_malloc(sizeof(gkyl_efit));
  up->rzbasis = rzbasis;
  up->rzgrid = rzgrid;
  up->rzlocal = rzlocal;
  up->rzlocal_ext = rzlocal_ext;
  up->use_gpu = use_gpu;

  FILE *ptr = fopen(filepath,"r");
  size_t status;

  // Get the dimensions

  status = fscanf(ptr,"%d%d", &up->nr, &up->nz);
  printf("nr, nz = %d, %d\n",up->nr, up->nz );

  // Read the non-array parameters, all are doubles:
  // rdim,zdim,rcentr,rleft,zmid;
  // rmaxis,zmaxis,simag,sibry,bcentr;
  // current,simag,xdum,rmaxis,xdum;
  // zmaxis,xdum,sibry,xdum,xdum;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;

  status = fscanf(ptr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &rdim, &zdim, &rcentr, &rleft, &zmid, & rmaxis, &zmaxis, &simag, &sibry, &bcentr, & current, &simag, &xdum, &rmaxis, &xdum, & zmaxis, &xdum, &sibry, &xdum, &xdum);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%g sibry=%g bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", rdim, zdim, rcentr, rleft, zmid,  rmaxis, zmaxis, simag, sibry, bcentr,  current, simag, rmaxis,  zmaxis, sibry);


  // Now we need to make the grid
  up->zmin = zmid - zdim/2;
  up->zmax = zmid + zdim/2;
  up->rmin = rleft;
  up->rmax = rleft+rdim;
  double rzlower[2] = {up->rmin, up->zmin };
  double rzupper[2] = {up->rmax, up->zmax};
  int cells[2] = {up->nr, up->nz};
  gkyl_rect_grid_init(up->rzgrid, 2, rzlower, rzupper, cells);


  // Now we 4 of the 1d arrays, all of length nr :
  // fpol, pres, ffprim, pprime
  // I don't actually care about these so just read 4*nr times
  double dummy;
  for(int i = 0; i<4*up->nr; i++){
    status = fscanf(ptr, "%lf", &dummy);
  }

  // Now we are gonna wanna read psi
  int node_nums[2] = {up->nr, up->nz};
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, up->rzgrid->ndim, node_nums);
  struct gkyl_array *psizr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);

  // Now lets loop through
  int idx[2];
  for(int ir = 0; ir < up->nr; ir++){
    idx[1] = ir;
    for(int iz = 0; iz < up->nz; iz++){
      idx[0] = iz;
      double *psi_n = gkyl_array_fetch(psizr_n, gkyl_range_idx(&nrange, idx));
      status = fscanf(ptr,"%lf", psi_n);

    }
  }
  // We filled psizr_nodal
  fclose(ptr);

  // Now we need to create a range for the modal array up->psizr
  // Then we can loop through it with a nodal->modal conversion function
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(up->rzgrid, nghost, up->rzlocal_ext, up->rzlocal);



  up->psizr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal->volume);
  nodal_array_to_modal_psi(psizr_n, up->psizr, up->rzlocal, &nrange, up->rzbasis, up->rzgrid, 1);
  
  return up;
}

