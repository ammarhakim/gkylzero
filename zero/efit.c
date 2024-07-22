#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include<gkyl_rect_decomp.h>
#include <gkyl_efit.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_nodal_ops.h>
#include <assert.h>

gkyl_efit* gkyl_efit_new(const char *filepath, int rz_poly_order, 
  enum gkyl_basis_type rz_basis_type, int flux_poly_order, bool reflect, bool use_gpu)
{
  gkyl_efit *up = gkyl_malloc(sizeof(struct gkyl_efit));
  up->rzbasis = gkyl_malloc(sizeof(struct gkyl_basis));

  up->rzgrid = gkyl_malloc(sizeof( struct gkyl_rect_grid));
  up->rzlocal = gkyl_malloc(sizeof(struct gkyl_range));
  up->rzlocal_ext = gkyl_malloc(sizeof(struct gkyl_range));

  up->fluxbasis = gkyl_malloc(sizeof(struct gkyl_basis));
  up->fluxgrid = gkyl_malloc(sizeof(struct gkyl_rect_grid));
  up->fluxlocal = gkyl_malloc(sizeof(struct gkyl_range));
  up->fluxlocal_ext = gkyl_malloc(sizeof(struct gkyl_range));

  up->use_gpu = use_gpu;
  up->filepath = filepath;

  gkyl_cart_modal_serendip(up->fluxbasis, 1, flux_poly_order);
  switch (rz_basis_type){
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(up->rzbasis, 2, rz_poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(up->rzbasis, 2, rz_poly_order);
      break;
    default:
      assert(false);
      break;
  }

  FILE *ptr = fopen(up->filepath,"r");
  size_t status;

  // Get the dimensions

  status = fscanf(ptr,"%d%d", &up->nr, &up->nz);

  // Read the non-array parameters, all are doubles:
  // rdim,zdim,rcentr,rleft,zmid;
  // rmaxis,zmaxis,simag,sibry,bcentr;
  // current,simag,xdum,rmaxis,xdum;
  // zmaxis,xdum,sibry,xdum,xdum;
  //double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;

  status = fscanf(ptr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
    &up->rdim, &up->zdim, &up->rcentr, &up->rleft, &up->zmid, &up-> rmaxis, &up->zmaxis, 
    &up->simag, &up->sibry, &up->bcentr, &up-> current, &up->simag, &up->xdum, &up->rmaxis, 
    &up->xdum, &up-> zmaxis, &up->xdum, &up->sibry, &up->xdum, &up->xdum);


  // Set zmid to 0 for double null
  if (reflect) {
    up->zmid = 0.0;
    up->zmaxis = 0.0;
  }


  // Now we need to make the grid
  up->zmin = up->zmid - up->zdim/2;
  up->zmax = up->zmid + up->zdim/2;
  up->rmin = up->rleft;
  up->rmax = up->rleft+up->rdim;

  double rzlower[2] = {up->rmin, up->zmin };
  double rzupper[2] = {up->rmax, up->zmax};
  int rzcells[2] = {0};
  int rzghost[2] = {1,1};
  if(up->rzbasis->poly_order==1){
    rzcells[0] = up->nr-1;
    rzcells[1]= up->nz-1;
  }
  if(up->rzbasis->poly_order==2){
    rzcells[0] = (up->nr-1)/2;
    rzcells[1] = (up->nz-1)/2;
  }
  gkyl_rect_grid_init(up->rzgrid, 2, rzlower, rzupper, rzcells);
  gkyl_create_grid_ranges(up->rzgrid, rzghost, up->rzlocal_ext, up->rzlocal);

  double fluxlower[1] = {up->sibry};
  double fluxupper[1] = {up->simag};
  int fluxcells[1] = {0};
  int fluxghost[2] = {1,1};
  if(up->fluxbasis->poly_order==1){
    fluxcells[0] = up->nr-1;
  }
  if(up->fluxbasis->poly_order==2){
    fluxcells[0] = (up->nr-1)/2;
  }

  gkyl_rect_grid_init(up->fluxgrid, 1, fluxlower, fluxupper, fluxcells);
  gkyl_create_grid_ranges(up->fluxgrid, fluxghost, up->fluxlocal_ext, up->fluxlocal);


  // allocate the necessary arrays
  up->psizr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal_ext->volume);
  up->psibyrzr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal_ext->volume);
  up->psibyr2zr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis->num_basis, up->rzlocal_ext->volume);
  up->fpolflux = gkyl_array_new(GKYL_DOUBLE, up->fluxbasis->num_basis, up->fluxlocal_ext->volume);
  up->qflux = gkyl_array_new(GKYL_DOUBLE, up->fluxbasis->num_basis, up->fluxlocal_ext->volume);

  // Read fpol because we do want that
  int flux_node_nums[1] = {up->nr};
  struct gkyl_range flux_nrange;
  gkyl_range_init_from_shape(&flux_nrange, 1, flux_node_nums);
  struct gkyl_array *fpolflux_n = gkyl_array_new(GKYL_DOUBLE, 1, flux_nrange.volume);
  int fidx[1];
  for(int i = up->nr-1; i>=0; i--){
      fidx[0] = i;
      double *fpol_n= gkyl_array_fetch(fpolflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", fpol_n);
  }

  struct gkyl_nodal_ops *n2m_flux = gkyl_nodal_ops_new(up->fluxbasis, up->fluxgrid, false);
  gkyl_nodal_ops_n2m(n2m_flux, up->fluxbasis, up->fluxgrid, 
    &flux_nrange, up->fluxlocal, 1, fpolflux_n, up->fpolflux);

  // Now we 3 of the 1d arrays, all of length nr :
  // pres, ffprim, pprime
  // I don't actually care about these so just read 4*nr times
  for(int i = 0; i<3*up->nr; i++){
    status = fscanf(ptr, "%lf", &up->xdum);
  }

  // Now we are gonna wanna read psi
  int node_nums[2] = {up->nr, up->nz};
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, up->rzgrid->ndim, node_nums);
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
  struct gkyl_nodal_ops *n2m_rz = gkyl_nodal_ops_new(up->rzbasis, up->rzgrid, false);
  gkyl_nodal_ops_n2m(n2m_rz, up->rzbasis, up->rzgrid, &nrange, up->rzlocal, 1, psizr_n, up->psizr);
  gkyl_nodal_ops_n2m(n2m_rz, up->rzbasis, up->rzgrid, &nrange, up->rzlocal, 1, psibyrzr_n, up->psibyrzr);
  gkyl_nodal_ops_n2m(n2m_rz, up->rzbasis, up->rzgrid, &nrange, up->rzlocal, 1, psibyr2zr_n, up->psibyr2zr);
  gkyl_nodal_ops_release(n2m_rz);

  // Reflect psi psi/R and psi/R^2 for double null
  // Reflect DG coeffs rather than nodal data to avoid symmetry errors in n2m conversion
  if (reflect) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, up->rzlocal);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[1] < gkyl_range_shape(up->rzlocal,1)/2 +1 ) {
        int idx_change[2] = {iter.idx[0], gkyl_range_shape(up->rzlocal, 1) - iter.idx[1]+1};
        const double *coeffs_ref = gkyl_array_cfetch(up->psizr, gkyl_range_idx(up->rzlocal, iter.idx));
        double *coeffs  = gkyl_array_fetch(up->psizr, gkyl_range_idx(up->rzlocal, idx_change));
        up->rzbasis->flip_odd_sign( 1, coeffs_ref, coeffs);
        coeffs_ref = gkyl_array_cfetch(up->psibyrzr, gkyl_range_idx(up->rzlocal, iter.idx));
        coeffs  = gkyl_array_fetch(up->psibyrzr, gkyl_range_idx(up->rzlocal, idx_change));
        up->rzbasis->flip_odd_sign( 1, coeffs_ref, coeffs);
        coeffs_ref = gkyl_array_cfetch(up->psibyr2zr, gkyl_range_idx(up->rzlocal, iter.idx));
        coeffs  = gkyl_array_fetch(up->psibyr2zr, gkyl_range_idx(up->rzlocal, idx_change));
        up->rzbasis->flip_odd_sign( 1, coeffs_ref, coeffs);
      }
    }
  }
 
  // Now lets read the q profile
  struct gkyl_array *qflux_n = gkyl_array_new(GKYL_DOUBLE, 1, flux_nrange.volume);
  for(int i = up->nr-1; i>=0; i--){
      fidx[0] = i;
      double *q_n= gkyl_array_fetch(qflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", q_n);
  }
  gkyl_nodal_ops_n2m(n2m_flux, up->fluxbasis, up->fluxgrid, 
    &flux_nrange, up->fluxlocal, 1, qflux_n, up->qflux);
  gkyl_nodal_ops_release(n2m_flux);
  
  // Free nodal arrays
  gkyl_array_release(fpolflux_n);
  gkyl_array_release(psizr_n);
  gkyl_array_release(psibyrzr_n);
  gkyl_array_release(psibyr2zr_n);
  gkyl_array_release(qflux_n);
  // Done, don't care about the rest
  
  fclose(ptr);
  return up;
}

void gkyl_efit_release(gkyl_efit* up){
  gkyl_free(up->rzbasis);
  gkyl_free(up->rzgrid);
  gkyl_free(up->rzlocal);
  gkyl_free(up->rzlocal_ext);
  gkyl_free(up->fluxbasis);
  gkyl_free(up->fluxgrid);
  gkyl_free(up->fluxlocal);
  gkyl_free(up->fluxlocal_ext);
  gkyl_array_release(up->psizr);
  gkyl_array_release(up->psibyrzr);
  gkyl_array_release(up->psibyr2zr);
  gkyl_array_release(up->fpolflux);
  gkyl_array_release(up->qflux);
  gkyl_free(up);
}
