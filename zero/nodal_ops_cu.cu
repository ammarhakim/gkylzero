/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_nodal_ops.h>
}

__global__ static void
gkyl_nodal_ops_n2m_cu_kernel(const struct gkyl_basis *cbasis, 
  struct gkyl_rect_grid grid, struct gkyl_range nrange, struct gkyl_range update_range, 
  const struct gkyl_array *nodes, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld)
{
  double xc[GKYL_MAX_DIM];
  int idx[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  int nidx[3];
  long lin_nidx[num_basis];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < update_range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idx);
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idx);

    for (int i=0; i<num_basis; ++i) {
      const double *temp  = (const double *) gkyl_array_cfetch(nodes, i);
      for( int j = 0; j < grid->ndim; j++){
        if(cpoly_order==1){
              nidx[j] = iter.idx[j]-1 + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2){
              nidx[j] = 2*(iter.idx[j]-1) + (temp[j]+1) ;
        }
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    double *arr_p = (double *) gkyl_array_fetch(modal_fld, linc);
    double fao[num_basis*num_comp];
  
    for (int i=0; i<num_basis; ++i) {
      const double *temp = (const double *) gkyl_array_cfetch(nodal_fld, lin_nidx[i]);
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

void 
gkyl_nodal_ops_n2m_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  const struct gkyl_array *nodal_fld, struct gkyl_array *modal_fld) 
{
  int nblocks = update_range->nblocks;
  int nthreads = update_range->nthreads;

  gkyl_nodal_ops_n2m_cu_kernel<<<nblocks, nthreads>>>(cbasis->on_dev, *grid, 
    *nrange, *update_range, nodal_ops->nodes->on_dev, num_comp, 
    nodal_fld, modal_fld);
}

__global__ static void
gkyl_nodal_ops_m2n_cu_kernel(const struct gkyl_basis *cbasis, 
  struct gkyl_rect_grid grid, struct gkyl_range nrange, struct gkyl_range update_range, 
  const struct gkyl_array *nodes, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld)
{
  double xc[GKYL_MAX_DIM];
  int idx[GKYL_MAX_DIM];
  int num_basis = cbasis->num_basis;
  int cpoly_order = cbasis->poly_order;
  double fnodal[num_basis]; // to store nodal function values

  int nidx[3];
  long lin_nidx[num_basis];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < update_range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idx);
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idx);

    for (int i=0; i<num_basis; ++i) {
      const double *temp  = (const double *) gkyl_array_cfetch(nodes, i);
      for( int j = 0; j < grid->ndim; j++){
        if(cpoly_order==1){
          nidx[j] = iter.idx[j]-1 + (temp[j]+1)/2 ;
        }
        if (cpoly_order==2)
          nidx[j] = 2*iter.idx[j] + (temp[j] + 1) ;
      }
      lin_nidx[i] = gkyl_range_idx(nrange, nidx);
    }

    const double *arr_p = (const double *) gkyl_array_cfetch(modal_fld, linc);
    double fao[num_basis*num_comp];
  
    // already fetched modal coeffs in arr_p
    // Now we are going to need to fill the nodal values
    // So in the loop of num_nodes/num_basis we will fetch the nodal value at lin_nidx[i]
    // at each place do an eval_basis at logical coords
  
    for (int i=0; i<num_basis; ++i) {
      double *temp = (double *) gkyl_array_fetch(nodal_fld, lin_nidx[i]);
      const double *node_i  = (const double *) gkyl_array_cfetch(nodes, i);
      for (int j=0; j<num_comp; ++j) {
        temp[j] = cbasis->eval_expand(node_i, &arr_p[j*num_basis]);
      }
    }
  }
}

void 
gkyl_nodal_ops_m2n_cu(const struct gkyl_nodal_ops *nodal_ops, 
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid, 
  const struct gkyl_range *nrange, const struct gkyl_range *update_range, int num_comp, 
  struct gkyl_array *nodal_fld, const struct gkyl_array *modal_fld) 
{
  int nblocks = update_range->nblocks;
  int nthreads = update_range->nthreads;

  gkyl_nodal_ops_n2m_cu_kernel<<<nblocks, nthreads>>>(cbasis->on_dev, *grid, 
    *nrange, *update_range, nodal_ops->nodes->on_dev, num_comp, 
    nodal_fld, modal_fld);
}
