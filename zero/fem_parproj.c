#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>

static long
global_num_nodes(const int dim, const int poly_order, const int basis_type, const int parnum_cells)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_parproj_num_nodes_global_1x_ser_p1(parnum_cells);
    } else if (poly_order == 2) {
      return fem_parproj_num_nodes_global_1x_ser_p2(parnum_cells);
    } else if (poly_order == 3) {
      return fem_parproj_num_nodes_global_1x_ser_p3(parnum_cells);
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return fem_parproj_num_nodes_global_3x_ser_p1(parnum_cells);
      } else if (poly_order == 2) {
        return fem_parproj_num_nodes_global_3x_ser_p2(parnum_cells);
      } else if (poly_order == 3) {
        return fem_parproj_num_nodes_global_3x_ser_p3(parnum_cells);
      }
    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
      if (poly_order == 1) {
        return fem_parproj_num_nodes_global_3x_tensor_p1(parnum_cells);
      } else if (poly_order == 2) {
        return fem_parproj_num_nodes_global_3x_tensor_p2(parnum_cells);
      }
    }
  }
  assert(false);  // Other dimensionalities not supported.
  return -1;
}

static void
local_mass(const int dim, const int poly_order, const int basis_type, struct gkyl_mat *massout)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_parproj_mass_1x_ser_p1(massout);
    } else if (poly_order == 2) {
      return fem_parproj_mass_1x_ser_p2(massout);
    } else if (poly_order == 3) {
      return fem_parproj_mass_1x_ser_p3(massout);
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return fem_parproj_mass_3x_ser_p1(massout);
      } else if (poly_order == 2) {
        return fem_parproj_mass_3x_ser_p2(massout);
      } else if (poly_order == 3) {
        return fem_parproj_mass_3x_ser_p3(massout);
      }
    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
      if (poly_order == 1) {
        return fem_parproj_mass_3x_tensor_p1(massout);
      } else if (poly_order == 2) {
        return fem_parproj_mass_3x_tensor_p2(massout);
      }
    }
  }
  assert(false);  // Other dimensionalities not supported.
}

gkyl_fem_parproj*
gkyl_fem_parproj_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis basis,
  const bool isparperiodic, bool use_gpu)
{
  gkyl_fem_parproj *up = gkyl_malloc(sizeof(gkyl_fem_parproj));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_parproj_kernels));
  up->kernels_cu = up->kernels;
 #ifdef GKYL_HAVE_CUDA
   if (use_gpu) {
     up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_parproj_kernels));
   }
 #endif

  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis = basis.num_basis;
  up->basis_type = basis.b_type;
  up->poly_order = basis.poly_order;
  up->pardir = grid->ndim-1; // Assume parallel direction is always the last.
  up->isperiodic = isparperiodic;
  up->use_gpu = use_gpu;

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis]));

  // Local and local-ext ranges for whole-grid arrays.
  int ghost[PARPROJ_MAX_DIM];
  for (int d=0; d<up->ndim; d++) ghost[d] = 1;
  gkyl_create_grid_ranges(grid, ghost, &up->local_range_ext, &up->local_range);

  // Range of cells we'll solve apply the parallel projection
  // operator in, as a sub-range of up->local_range_ext.
  int sublower[PARPROJ_MAX_DIM], subupper[PARPROJ_MAX_DIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->local_range.lower[d];
    subupper[d] = up->local_range.upper[d];
  }
  if (up->isperiodic) {
    // Include ghost cells in parallel direction.
    sublower[up->pardir] = up->local_range_ext.lower[up->pardir];
    subupper[up->pardir] = up->local_range_ext.upper[up->pardir];
  }
  gkyl_sub_range_init(&up->solve_range, &up->local_range_ext, sublower, subupper);

  // Range of parallel cells, as a sub-range of up->solve_range.
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range.lower[d];
    subupper[d] = up->solve_range.lower[d];
  }
  subupper[up->pardir] = up->solve_range.upper[up->pardir];
  gkyl_sub_range_init(&up->par_range, &up->solve_range, sublower, subupper);
  up->parnum_cells = up->par_range.volume;
  // Range of perpendicular cells.
  gkyl_range_shorten(&up->perp_range, &up->solve_range, up->pardir, 1);

  // 1D range of parallel cells.
  int lower1d[] = {up->par_range.lower[up->pardir]}, upper1d[] = {up->par_range.upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);
  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim==3 ? 2 : 1, up->perp_range.lower, up->perp_range.upper);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = global_num_nodes(up->ndim, up->poly_order, basis.b_type, up->par_range.volume);

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*up->perp_range.volume); // Global right side vector.

  // Create local matrices used later.
  up->local_mass = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  local_mass(up->ndim, up->poly_order, basis.b_type, up->local_mass);

  // Select local-to-global mapping kernel:
  up->kernels->l2g = fem_parproj_choose_local2global_kernel(up->ndim, basis.b_type, up->poly_order);

  // Select RHS source kernel:
  up->kernels->srcker = fem_parproj_choose_srcstencil_kernel(up->ndim, basis.b_type, up->poly_order);

  // Select kernel that fetches the solution:
  up->kernels->solker = fem_parproj_choose_solstencil_kernel(up->ndim, basis.b_type, up->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    fem_parproj_choose_kernels_cu(&basis, up->isperiodic, up->kernels_cu);
  }
#endif

  // Create a linear Ax=B problem. We envision two cases:
  //  a) No weight, or weight is a scalar so we can divide the RHS by it. Then
  //     A is the same for every problem, and we can just populate B with a
  //     column for each problem.
  //  b) Weight depends on space. Then we have to create an A matrix and a
  //     separate Ax=B problem for each perpendicular cell.
  // For now we restrict ourselves to a).
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->prob_cu = gkyl_cusolver_prob_new(up->numnodes_global, up->numnodes_global, up->perp_range.volume);
  else
    up->prob = gkyl_superlu_prob_new(up->numnodes_global, up->numnodes_global, up->perp_range.volume);
#else
    up->prob = gkyl_superlu_prob_new(up->numnodes_global, up->numnodes_global, up->perp_range.volume);
#endif

  // Assign non-zero elements in A.
  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri);
#endif
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    up->kernels->l2g(up->parnum_cells, paridx, up->globalidx);
    for (size_t k=0; k<up->numnodes_local; ++k) {
      for (size_t m=0; m<up->numnodes_local; ++m) {
        double val = gkyl_mat_get(up->local_mass,k,m);
        long globalidx_k = up->globalidx[k];
        long globalidx_m = up->globalidx[m];
        gkyl_mat_triples_accum(tri, globalidx_k, globalidx_m, val);
      }
    }
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cusolver_amat_from_triples(up->prob_cu, tri);
  } else {
    gkyl_superlu_amat_from_triples(up->prob, tri);
  }
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

  gkyl_mat_triples_release(tri);

  return up;
}

void
gkyl_fem_parproj_set_rhs(gkyl_fem_parproj* up, const struct gkyl_array *rhsin)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_parproj_set_rhs_cu(up, rhsin);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {

    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);
    const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

    int idx1d[] = {up->solve_iter.idx[up->ndim-1]};
    long paridx = gkyl_range_idx(&up->par_range1d, idx1d);
    up->kernels->l2g(up->parnum_cells, paridx, up->globalidx);

    int idx2d[] = {up->perp_range2d.lower[0], up->perp_range2d.lower[0]};
    for (int d=0; d<up->ndim-1; d++) idx2d[d] = up->solve_iter.idx[d];
    long perpidx2d = gkyl_range_idx(&up->perp_range2d, idx2d);
    long perpProbOff = perpidx2d*up->numnodes_global;

    up->kernels->srcker(rhsin_p, perpProbOff, up->globalidx, brhs_p);

  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_parproj_solve(gkyl_fem_parproj* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_parproj_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {

    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);
    double *phiout_p = gkyl_array_fetch(phiout, linidx);

    int idx1d[] = {up->solve_iter.idx[up->ndim-1]};
    long paridx = gkyl_range_idx(&up->par_range1d, idx1d);
    up->kernels->l2g(up->parnum_cells, paridx, up->globalidx);

    int idx2d[] = {up->perp_range2d.lower[0], up->perp_range2d.lower[0]};
    for (int d=0; d<up->ndim-1; d++) idx2d[d] = up->solve_iter.idx[d];
    long perpidx2d = gkyl_range_idx(&up->perp_range2d, idx2d);
    long perpProbOff = perpidx2d*up->numnodes_global;

    up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), perpProbOff, up->globalidx, phiout_p);

  }

}

void gkyl_fem_parproj_release(gkyl_fem_parproj *up)
{
  gkyl_mat_release(up->local_mass);
#ifdef GKYL_HAVE_CUDA
  gkyl_cu_free(up->kernels_cu);
  if (up->use_gpu) {
    gkyl_cusolver_prob_release(up->prob_cu);
  } else {
    gkyl_superlu_prob_release(up->prob);
  }
#else
  gkyl_superlu_prob_release(up->prob);
#endif
  gkyl_free(up->globalidx);
  gkyl_free(up->brhs);
  gkyl_free(up->kernels);
  gkyl_free(up);
}
