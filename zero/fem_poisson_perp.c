#include <gkyl_fem_poisson_perp.h>
#include <gkyl_fem_poisson_perp_priv.h>
#include <gkyl_array_average.h>

struct gkyl_fem_poisson_perp*
gkyl_fem_poisson_perp_new(const struct gkyl_range *solve_range, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis basis, struct gkyl_poisson_bc *bcs, struct gkyl_array *epsilon,
  struct gkyl_array *kSq, bool use_gpu)
{

  struct gkyl_fem_poisson_perp *up = gkyl_malloc(sizeof(struct gkyl_fem_poisson_perp));

  up->solve_range = solve_range;
  up->ndim = grid->ndim;
  up->ndim_perp = up->ndim-1;
  up->grid = *grid;
  up->num_basis =  basis.num_basis;
  up->basis_type = basis.b_type;
  up->poly_order = basis.poly_order;
  up->pardir = grid->ndim-1; // Assume parallel direction is always the last.
  up->basis = basis;
  up->use_gpu = use_gpu;
  up->epsilon = epsilon;

  assert(up->ndim > 1);
  assert(up->epsilon->ncomp == (2*(up->ndim-1)-1)*basis.num_basis);

  // We assume epsilon and kSq live on the device, and we create a host-side
  // copies temporarily to compute the LHS matrix. This also works for CPU solves.
  struct gkyl_array *epsilon_ho = gkyl_array_new(GKYL_DOUBLE, up->epsilon->ncomp, up->epsilon->size);
  gkyl_array_copy(epsilon_ho, up->epsilon);
  struct gkyl_array *kSq_ho;
  if (kSq) {
    up->ishelmholtz = true;
    kSq_ho = gkyl_array_new(GKYL_DOUBLE, kSq->ncomp, kSq->size);
    gkyl_array_copy(kSq_ho, kSq);
  } else {
    up->ishelmholtz = false;
    kSq_ho = gkyl_array_new(GKYL_DOUBLE, up->num_basis, 1);
    gkyl_array_clear(kSq_ho, 0.);
  }

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis])); // global index, one for each basis in a cell.

  for (int d=0; d<up->ndim; d++) up->num_cells[d] = up->solve_range->upper[d]-up->solve_range->lower[d]+1;

  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim_perp, up->solve_range->lower, up->solve_range->upper);
  // 1D range of parallel cells.
  int lower1d[] = {up->solve_range->lower[up->pardir]}, upper1d[] = {up->solve_range->upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);

  // Range of perpendicular cells at each parallel location.
  up->perp_range = (struct gkyl_range *) gkyl_malloc(up->par_range1d.volume * sizeof(struct gkyl_range));
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);
    int removeDim[] = {0,0,0},  loc[] = {0,0,0};
    removeDim[up->pardir] = 1;
    loc[up->pardir] = paridx;
    gkyl_range_deflate(&up->perp_range[paridx], up->solve_range, removeDim, loc);
  }
  // Range of parallel cells.
  int sublower[GKYL_MAX_CDIM], subupper[GKYL_MAX_CDIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->lower[d];
  }
  subupper[up->pardir] = up->solve_range->upper[up->pardir];
  gkyl_sub_range_init(&up->par_range, up->solve_range, sublower, subupper);

  // Prepare for periodic domain case.
  for (int d=0; d<up->ndim_perp; d++) {
    // Sanity check.
    if ((bcs->lo_type[d] == GKYL_POISSON_PERIODIC && bcs->up_type[d] != GKYL_POISSON_PERIODIC) ||
        (bcs->lo_type[d] != GKYL_POISSON_PERIODIC && bcs->up_type[d] == GKYL_POISSON_PERIODIC))
      assert(false);
  }
  for (int d=0; d<up->ndim_perp; d++) up->isdirperiodic[d] = bcs->lo_type[d] == GKYL_POISSON_PERIODIC;
  up->isdomperiodic = true;
  for (int d=0; d<up->ndim_perp; d++) up->isdomperiodic = up->isdomperiodic && up->isdirperiodic[d];

  if (up->isdomperiodic) {
    // Create operator to compute the perpendicular average of the rhs.
    struct gkyl_basis perpavg_basis;
    if (basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY)
      gkyl_cart_modal_serendip(&perpavg_basis, 1, basis.poly_order);
    else if (basis.b_type == GKYL_BASIS_MODAL_TENSOR)
      gkyl_cart_modal_tensor(&perpavg_basis, 1, basis.poly_order);

    struct gkyl_range perpavg_local_ext;
    gkyl_range_extend(&perpavg_local_ext, &up->par_range1d, &((int){1}), &((int){1})); // Assume 1 ghost cell.
    gkyl_sub_range_init(&up->perpavg_local, &perpavg_local_ext, &((int){up->par_range1d.lower[0]}), &((int){up->par_range1d.upper[0]}));

    int avg_dim[GKYL_MAX_CDIM] = {1,1,1};
    avg_dim[up->pardir] = 0;

    struct gkyl_array_average_inp avg_inp = {
      .grid = grid,
      .basis = basis,
      .local = solve_range,
      .basis_avg = perpavg_basis,
      .local_avg = &up->perpavg_local,
      .local_avg_ext = &perpavg_local_ext,
      .weight = NULL,
      .avg_dim = avg_dim,
      .use_gpu = use_gpu,
    };
    up->perpavg_op = gkyl_array_average_new(&avg_inp);
    up->perpavg_rhs = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, perpavg_basis.num_basis, perpavg_local_ext.volume)
                             : gkyl_array_new(GKYL_DOUBLE, perpavg_basis.num_basis, perpavg_local_ext.volume);
  }

  // Pack BC values into a single array for easier use in kernels.
  for (int d=0; d<up->ndim_perp; d++) {
    for (int k=0; k<6; k++) up->bcvals[d*2*3+k] = 0.0; // default. Not used in some cases (e.g. periodic).
    if (bcs->lo_type[d] != GKYL_POISSON_PERIODIC) {
      int vnum, voff;
      vnum = 1;
      voff = 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+k] = bcs->lo_value[d].v[k];

      vnum = 1;
      voff = 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+3+k] = bcs->up_value[d].v[k];
    }
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->bcvals_cu = (double *) gkyl_cu_malloc(sizeof(double[PERP_DIM_MAX*3*2]));
    gkyl_cu_memcpy(up->bcvals_cu, up->bcvals, sizeof(double[PERP_DIM_MAX*3*2]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_poisson_perp_global_num_nodes(up->ndim, up->poly_order, basis.b_type, up->num_cells, up->isdirperiodic);
//  up->numnodes_global = 1+gkyl_fem_poisson_perp_global_num_nodes(up->ndim, up->poly_order, basis.b_type, up->num_cells, up->isdirperiodic);

  for (int d=0; d<up->ndim; d++) up->dx[d] = up->grid.dx[d];  // Cell lengths.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->dx_cu = (double *) gkyl_cu_malloc(sizeof(double[GKYL_MAX_CDIM]));
    gkyl_cu_memcpy(up->dx_cu, up->dx, sizeof(double[GKYL_MAX_CDIM]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*up->par_range.volume); // Global right side vector.

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_poisson_perp_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_poisson_perp_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  // Select local-to-global mapping kernels:
  fem_poisson_perp_choose_local2global_kernels(&basis, up->isdirperiodic, up->kernels->l2g);

  // Select lhs kernels:
  fem_poisson_perp_choose_lhs_kernels(&basis, bcs, up->kernels->lhsker);

  // Select rhs src kernels:
  fem_poisson_perp_choose_src_kernels(&basis, bcs, up->kernels->srcker);

  // Select sol kernel:
  up->kernels->solker = fem_poisson_perp_choose_sol_kernels(&basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_poisson_perp_choose_kernels_cu(&basis, bcs, up->isdirperiodic, up->kernels_cu);
#endif

  // Create a linear Ax=B problem for each perp plane. Here A is the discrete (global)
  // matrix representation of the LHS of the perpendiculat Helmholtz equation.
  // cuSolverRF may support for A_i x_i = B_i, so we may revisit this
  // structure for the GPU solve.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->prob_cu = gkyl_culinsolver_prob_new(up->par_range.volume, up->numnodes_global, up->numnodes_global, 1);
  } else {
    up->prob = gkyl_superlu_prob_new(up->par_range.volume, up->numnodes_global, up->numnodes_global, 1);
  }
#else
  up->prob = gkyl_superlu_prob_new(up->par_range.volume, up->numnodes_global, up->numnodes_global, 1);
#endif

  struct gkyl_mat_triples **tri = gkyl_malloc(up->par_range.volume*sizeof(struct gkyl_mat_triples *));
  for (size_t i=0; i<up->par_range.volume; i++) {
    tri[i] = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri[i]);
#endif
  }

  // Assign non-zero elements in A.
  int idx0[GKYL_MAX_CDIM],  idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *eps_p = gkyl_array_fetch(epsilon_ho, linidx);
      double *kSq_p = up->ishelmholtz? gkyl_array_fetch(kSq_ho, linidx) : gkyl_array_fetch(kSq_ho,0);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      // Apply the -nabla . (epsilon*nabla_perp)-kSq stencil.
      keri = idx_to_inloup_ker(up->ndim_perp, up->num_cells, idx1);
      up->kernels->lhsker[keri](eps_p, kSq_p, up->dx, up->bcvals, up->globalidx, tri[paridx]);
    }
  }

//  // Augment the matrix to make it regular (i.e. not singular).
//  // Need to add 1 to up->numnodes_global where it is first assigned.
//  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
//  while (gkyl_range_iter_next(&up->par_iter1d)) {
//    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);
//    for (int i=0; i<up->numnodes_global-1; i++) {
//      gkyl_mat_triples_insert(tri[paridx], i, up->numnodes_global-1, 1.0);
//      gkyl_mat_triples_insert(tri[paridx], up->numnodes_global-1, i, 1.0);
//    }
//  }

  // Set Dirichlet BC for node on the left.
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    for (int i=0; i<up->numnodes_global; i++) {
      gkyl_mat_triples_insert(tri[paridx], 0, i, 0.0);
    }
    gkyl_mat_triples_insert(tri[paridx], 0, 0, 1.0);
    for (int i=0; i<up->numnodes_global; i++) {
      gkyl_mat_triples_insert(tri[paridx], 1, i, 0.0);
    }
    gkyl_mat_triples_insert(tri[paridx], 1, 1, 1.0);
  }

// // Write the matrix elements in a file
//   gkyl_mat_triples_iter *iter2 = gkyl_mat_triples_iter_new(tri[0]);
//     FILE *file = fopen("A.txt", "w");
//     for (size_t i=0; i<gkyl_mat_triples_size(tri[0]); ++i) {
//       gkyl_mat_triples_iter_next(iter2); // bump iterator.
//       struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter2);
//       size_t idx[2] = { mt.row, mt.col };
//    
//       fprintf(file,"a(%zu,%zu) = %g\n",idx[0],idx[1],mt.val);
//     }
//   gkyl_mat_triples_iter_release(iter2);
//   fclose(file); 

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_culinsolver_amat_from_triples(up->prob_cu, tri);
  else
    gkyl_superlu_amat_from_triples(up->prob, tri);
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

//  gkyl_superlu_print_amat(up->prob);

  for (size_t i=0; i<up->par_range.volume; i++)
    gkyl_mat_triples_release(tri[i]);
  gkyl_free(tri);

  gkyl_array_release(epsilon_ho);
  gkyl_array_release(kSq_ho);

  return up;
}

void
gkyl_fem_poisson_perp_set_rhs(gkyl_fem_poisson_perp *up, struct gkyl_array *rhsin)
{

//  if (up->isdomperiodic && !(up->ishelmholtz)) {
//    // Subtract the volume averaged RHS from the RHS.
//    gkyl_array_average_advance(up->perpavg_op, rhsin, up->perpavg_rhs);
//
//    if (!up->use_gpu) {
//      struct gkyl_range_iter iter;
//      gkyl_range_iter_init(&iter, up->solve_range);
//      while (gkyl_range_iter_next(&iter)) {
//        long linidx = gkyl_range_idx(up->solve_range, iter.idx);
//        double *rhsin_p = gkyl_array_fetch(rhsin, linidx);
//
//        int idx_par[] = {iter.idx[up->pardir]};
//        long linidx_par = gkyl_range_idx(&up->perpavg_local, idx_par);
//        double *avg_p = gkyl_array_fetch(up->perpavg_rhs, linidx_par);
//
//        rhsin_p[0] += -pow(sqrt(2.0),up->ndim_perp)*avg_p[0];
//        rhsin_p[up->ndim] += -pow(sqrt(2.0),up->ndim_perp)*avg_p[1];
//      }
//    }
//#ifdef GKYL_HAVE_CUDA
//    if (up->use_gpu) {
//      gkyl_fem_poisson_perp_subtract_perpavg_cu(up, up->perpavg_rhs, rhsin);
//    }
//#endif
//  }


#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_poisson_perp_set_rhs_cu(up, rhsin);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *eps_p = gkyl_array_fetch(up->epsilon, linidx);
      double *rhsin_p = gkyl_array_fetch(rhsin, linidx);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d] - 1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      // Apply the RHS source stencil. It's mostly the mass matrix times a
      // modal-to-nodal operator times the source, modified by BCs in skin cells.
      keri = idx_to_inloup_ker(up->ndim_perp, up->num_cells, idx1);

      long parProbOff = paridx*up->numnodes_global;

      up->kernels->srcker[keri](eps_p, up->dx, rhsin_p, up->bcvals, parProbOff, up->globalidx, brhs_p);
    }

  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_poisson_perp_solve(gkyl_fem_poisson_perp *up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_poisson_perp_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  gkyl_array_clear(phiout, 0.0);

  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);
  
    gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
    while (gkyl_range_iter_next(&up->perp_iter2d)) {
      long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

      for (size_t d=0; d<up->ndim_perp; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *phiout_p = gkyl_array_fetch(phiout, linidx);

      int keri = idx_to_inup_ker(up->ndim_perp, up->num_cells, up->perp_iter2d.idx);
      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d]-1;
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

      long parProbOff = paridx*up->numnodes_global;

      up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), parProbOff, up->globalidx, phiout_p);
    }
  }

}

void gkyl_fem_poisson_perp_release(struct gkyl_fem_poisson_perp *up)
{
  if (up->isdomperiodic) {
    gkyl_array_average_release(up->perpavg_op);
    gkyl_array_release(up->perpavg_rhs);
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
    gkyl_cu_free(up->dx_cu);
    gkyl_cu_free(up->bcvals_cu);
    gkyl_culinsolver_prob_release(up->prob_cu);
  } else {
    gkyl_superlu_prob_release(up->prob);
  }
#else
  gkyl_superlu_prob_release(up->prob);
#endif

  gkyl_array_release(up->brhs);
  gkyl_free(up->kernels);
  gkyl_free(up->perp_range);
  gkyl_free(up->globalidx);
  gkyl_free(up);
}
