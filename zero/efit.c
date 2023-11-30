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
#include <gkyl_nodal_ops.h>

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



gkyl_efit* gkyl_efit_new(const char *filepath, const struct gkyl_basis *rzbasis, const struct gkyl_basis *fluxbasis, bool use_gpu)
{

  gkyl_efit *up = gkyl_malloc(sizeof(gkyl_efit));
  up->rzbasis = rzbasis;
  up->fluxbasis = fluxbasis;
  //up->rzgrid = rzgrid;
  //up->rzlocal = rzlocal;
  //up->rzlocal_ext = rzlocal_ext;
  up->use_gpu = use_gpu;
  up->filepath = filepath;

  FILE *ptr = fopen(up->filepath,"r");
  size_t status;

  // Get the dimensions

  status = fscanf(ptr,"%d%d", &up->nr, &up->nz);
  printf("nr, nz = %d, %d\n",up->nr, up->nz );

  // Read the non-array parameters, all are doubles:
  // rdim,zdim,rcentr,rleft,zmid;
  // rmaxis,zmaxis,simag,sibry,bcentr;
  // current,simag,xdum,rmaxis,xdum;
  // zmaxis,xdum,sibry,xdum,xdum;
  //double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;

  status = fscanf(ptr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &up->rdim, &up->zdim, &up->rcentr, &up->rleft, &up->zmid, &up-> rmaxis, &up->zmaxis, &up->simag, &up->sibry, &up->bcentr, &up-> current, &up->simag, &up->xdum, &up->rmaxis, &up->xdum, &up-> zmaxis, &up->xdum, &up->sibry, &up->xdum, &up->xdum);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", up->rdim, up->zdim, up->rcentr, up->rleft, up->zmid, up->rmaxis, up->zmaxis, up->simag, up->sibry, up->bcentr, up-> current, up->simag, up->rmaxis, up-> zmaxis, up->sibry);

  fclose(ptr);


  // Now we need to make the grid
  up->zmin = up->zmid - up->zdim/2;
  up->zmax = up->zmid + up->zdim/2;
  up->rmin = up->rleft;
  up->rmax = up->rleft+up->rdim;

  up->rzlower = gkyl_malloc(2*sizeof(double));
  up->rzupper = gkyl_malloc(2*sizeof(double));
  up->rzcells = gkyl_malloc(2*sizeof(double));
  up->rzghost = gkyl_malloc(2*sizeof(double));

  up->rzlower[0] = up->rmin;
  up->rzlower[1] = up->zmin;
  up->rzupper[0] = up->rmax;
  up->rzupper[1] = up->zmax;
  up->rzghost[0] = 1;
  up->rzghost[1] = 1;

  if(up->rzbasis->poly_order==1){
    up->rzcells[0] = up->nr-1;
    up->rzcells[1]= up->nz-1;
  }
  if(up->rzbasis->poly_order==2){
    up->rzcells[0] = (up->nr-1)/2;
    up->rzcells[1] = (up->nz-1)/2;
  }

  // Now we need to make the flux grid
  up->fluxlower = gkyl_malloc(1*sizeof(double));
  up->fluxupper = gkyl_malloc(1*sizeof(double));
  up->fluxcells = gkyl_malloc(1*sizeof(double));
  up->fluxghost = gkyl_malloc(1*sizeof(double));

  up->fluxlower[0] = up->sibry;
  up->fluxupper[0] = up->simag;
  up->fluxghost[0] = 1;

  if(up->fluxbasis->poly_order==1){
    up->fluxcells[0] = up->nr-1;
  }

  if(up->fluxbasis->poly_order==2){
    up->fluxcells[0] = (up->nr-1)/2;
  }

  return up;
}

void gkyl_efit_advance(gkyl_efit* up, struct gkyl_rect_grid* rzgrid, struct gkyl_rect_grid* fluxgrid, struct gkyl_range* rzlocal, struct gkyl_range* rzlocal_ext, struct gkyl_array* psizr, struct gkyl_array* psibyrzr,struct gkyl_array* psibyr2zr, struct gkyl_range* fluxlocal, struct gkyl_range* fluxlocal_ext, struct gkyl_array* fpolflux)
{
  // Do this in g2 now
  //gkyl_rect_grid_init(up->rzgrid, 2, rzlower, rzupper, cells);

  // Skip all the dummy args we already read
  FILE *ptr = fopen(up->filepath,"r");
  size_t status;
  int idummy;
  double ddummy;

  status = fscanf(ptr,"%d%d", &idummy, &idummy);
  status = fscanf(ptr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy, &ddummy);


  // Read pol because we do want that
  int flux_node_nums[1] = {up->nr};
  struct gkyl_range flux_nrange;
  gkyl_range_init_from_shape(&flux_nrange, 1, flux_node_nums);
  struct gkyl_array *fpolflux_n = gkyl_array_new(GKYL_DOUBLE, 1, flux_nrange.volume);
  int fidx[1];
  for(int i = 0; i<up->nr; i++){
      fidx[0] = i;
      double *fpol_n= gkyl_array_fetch(fpolflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", fpol_n);
  }

  int bcs[3] = {0,0,0};
  gkyl_nodal_ops_n2m( up->fluxbasis, fluxgrid, &flux_nrange, fluxlocal, 1, fpolflux_n, fpolflux, bcs);

  // Now we 3 of the 1d arrays, all of length nr :
  // pres, ffprim, pprime
  // I don't actually care about these so just read 4*nr times
  for(int i = 0; i<3*up->nr; i++){
    status = fscanf(ptr, "%lf", &ddummy);
  }

  // Now we are gonna wanna read psi
  int node_nums[2] = {up->nr, up->nz};
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, rzgrid->ndim, node_nums);
  struct gkyl_array *psizr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  struct gkyl_array *psibyrzr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  struct gkyl_array *psibyr2zr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);

  // Now lets loop through
  // Not only do we want psi at the nodes, we also want psi/R and psi/R^2 so we can use them for the magnetc field
  double R = up->rmin;
  double dR = up->rdim/(up->nr-1);
  int idx[2];
  for(int iz = 0; iz < up->nz; iz++){
    idx[1] = iz;
    for(int ir = 0; ir < up->nr; ir++){
      R = up->rmin+ir*dR;
      idx[0] = ir;
      // set psi
      double *psi_n = gkyl_array_fetch(psizr_n, gkyl_range_idx(&nrange, idx));
      status = fscanf(ptr,"%lf", psi_n);
      // set psibyr and psibyr2
      double *psibyr_n = gkyl_array_fetch(psibyrzr_n, gkyl_range_idx(&nrange, idx));
      double *psibyr2_n = gkyl_array_fetch(psibyr2zr_n, gkyl_range_idx(&nrange, idx));
      psibyr_n[0] = psi_n[0]/R;
      psibyr2_n[0] = psi_n[0]/R/R;
    }
  }
  // We filled psizr_nodal
  fclose(ptr);

  // Do the ranges in g2
  // Now we need to create a range for the modal array up->psizr
  // Then we can loop through it with a nodal->modal conversion function
  //int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  //gkyl_create_grid_ranges(rzgrid, nghost, rzlocal_ext, rzlocal);



  // Allocate these in g2 now
  //up->psizr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal_ext->volume);
  //up->psibyrzr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal_ext->volume);
  //up->psibyr2zr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal_ext->volume);

  nodal_array_to_modal_psi(psizr_n, psizr, rzlocal, &nrange, up->rzbasis, rzgrid, 1);
  nodal_array_to_modal_psi(psibyrzr_n, psibyrzr, rzlocal, &nrange, up->rzbasis, rzgrid, 1);
  nodal_array_to_modal_psi(psibyr2zr_n, psibyr2zr, rzlocal, &nrange, up->rzbasis, rzgrid, 1);
}


void gkyl_efit_release(gkyl_efit* up){
  //gkyl_array_release(up->psizr);
  //gkyl_array_release(up->psibyrzr);
  //gkyl_array_release(up->psibyr2zr);
  gkyl_free(up->rzlower);
  gkyl_free(up->rzupper);
  gkyl_free(up->rzcells);
  gkyl_free(up->rzghost);
  gkyl_free(up);
}
