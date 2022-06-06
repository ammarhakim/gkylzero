void 
gkyl_fem_poisson_set_rhs_cu(gkyl_fem_poisson *up, struct gkyl_array *rhsin)
{
  gkyl_fem_poisson_set_rhs_kernel<<<dG, dB>>>(up->rhs, rhsin->on_dev, &up->solve_range, up->globalidx_cu, up->bcvals); 
}	

void
gkyl_fem_poisson_solve_cu(gkyl_fem_poisson *up, struct gkyl_array *phiin)
{
  // do linear solve with cusolver
  gkyl_cusolver_solve(up->prob_cu);

  gkyl_fem_poisson_get_sol_kernel<<<dG, dB>>>(up->rhs, rhsin->on_dev, &up->solve_range, up->globalidx_cu, up->bcvals); 
}

__global__ void
gkyl_fem_poisson_set_rhs_kernel(double *rhs_global, struct gkyl_array *rhs_local, struct gkyl_range range, long *globalidx, double *bcvals)
{
  int idx[GKYL_MAX_DIM];
  int idx0[GKYL_MAX_DIM];
  int num_cells[POISSON_MAX_DIM];
  for (int d=0; d<POISSON_MAX_DIM; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);
    
    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);
    
    const double *local_d = (const double*) gkyl_array_cfetch(rhs_local, start);

    int keri = idx_to_inup_ker(range->ndim, num_cells, idx);

    for (size_t d=0; d<range->ndim; d++) idx0[d] = idx[d]-1;

    l2g[keri](num_cells, idx0, globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(range->ndim, num_cells, idx);
    srcker[keri](local_d, bcvals, globalidx, rhs_global);
  }
}
