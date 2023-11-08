#include <gkyl_nodal_ops.h>

void gkyl_nodal_ops_n2m(const struct gkyl_basis* cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, const struct gkyl_array* nodal_fld, struct gkyl_array *modal_fld){
  double xc[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, cbasis->num_basis);
  cbasis->node_list(gkyl_array_fetch(nodes, 0));
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
        if(cpoly_order==1){
            nidx[j] = iter.idx[j]-1 + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2)
          nidx[j] = 2*iter.idx[j] + (temp[j] + 1) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *arr_p = gkyl_array_fetch(modal_fld, lidx);
    double fao[num_basis*num_comp];
  
    for (int i=0; i<num_basis; ++i) {
      const double* temp = gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
      for (int j=0; j<num_comp; ++j) {
        fao[i*num_comp + j] = temp[j];
      }
    }

    for (int i=0; i<num_comp; ++i) {
      // copy so nodal values for each return value are contiguous
      // (recall that function can have more than one return value)
      for (int k=0; k<num_basis; ++k)
        fnodal[k] = fao[num_comp*k+i];
      // transform to modal expansion
      cbasis->nodal_to_modal(fnodal, &arr_p[num_basis*i]);
    }
  }
}


void gkyl_nodal_ops_m2n(const struct gkyl_basis* cbasis, const struct gkyl_rect_grid *grid, const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, struct gkyl_array* nodal_fld, const struct gkyl_array *modal_fld){
  double xc[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid->ndim, cbasis->num_basis);
  cbasis->node_list(gkyl_array_fetch(nodes, 0));
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
        if(cpoly_order==1){
            nidx[j] = iter.idx[j]-1 + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2)
          nidx[j] = 2*iter.idx[j] + (temp[j] + 1) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    const double *arr_p = gkyl_array_cfetch(modal_fld, lidx);
    double fao[num_basis*num_comp];

    // already fetched modal coeffs in arr_p
    // Now we are going to need to fill the nodal values
    // So in the loop of num_nodes/num_basis we will fetch the nodal value at lin_nidx[i]
    // at each place do an eval_basis at logical coords
  
    for (int i=0; i<num_basis; ++i) {
      double* temp = gkyl_array_fetch(nodal_fld, lin_nidx[i]);
      const double*  node_i  = gkyl_array_cfetch(nodes,i);
      for (int j=0; j<num_comp; ++j) {
        temp[j] = cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
      }
    }

  }
}
