#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>

struct gkyl_fem_parproj*
gkyl_fem_parproj_new(const struct gkyl_range *solve_range, const struct gkyl_range *solve_range_ext, 
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight, bool use_gpu)
{
  struct gkyl_fem_parproj *up = gkyl_malloc(sizeof(struct gkyl_fem_parproj));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_parproj_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_parproj_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  up->solve_range = solve_range;
  up->solve_range_ext = solve_range_ext;
  up->ndim = solve_range->ndim;
  up->num_basis  = basis->num_basis;
  up->basis_type = basis->b_type;
  up->poly_order = basis->poly_order;
  up->pardir = up->ndim-1; // Assume parallel direction is always the last.
  up->isperiodic = bctype == GKYL_FEM_PARPROJ_PERIODIC;
  up->isdirichlet = bctype == GKYL_FEM_PARPROJ_DIRICHLET;
  up->use_gpu = use_gpu;

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis]));

  bool isweighted = false;
  if (weight) isweighted = true;

  // Range of parallel cells, as a sub-range of up->solve_range.
  int sublower[GKYL_MAX_CDIM], subupper[GKYL_MAX_CDIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->lower[d];
  }
  subupper[up->pardir] = up->solve_range->upper[up->pardir];
  gkyl_sub_range_init(&up->par_range, up->solve_range, sublower, subupper);
  up->parnum_cells = up->par_range.volume;

  // Range of perpendicular cells.
  gkyl_range_shorten_from_above(&up->perp_range, up->solve_range, up->pardir, 1);

  // 1D range of parallel cells.
  int lower1d[] = {up->par_range.lower[up->pardir]}, upper1d[] = {up->par_range.upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);
  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim==3 ? 2 : 1, up->perp_range.lower, up->perp_range.upper);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_parproj_global_num_nodes(basis, up->isperiodic, up->par_range.volume);

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*up->perp_range.volume); // Global right side vector.

  // Select local-to-global mapping kernel:
  fem_parproj_choose_local2global_kernel(basis, up->isperiodic, up->kernels->l2g);

  // Select weighted LHS kernel (not always used):
  fem_parproj_choose_lhs_kernel(basis, up->isdirichlet, isweighted, up->kernels->lhsker);

  // Select RHS source kernel:
  fem_parproj_choose_srcstencil_kernel(basis, up->isdirichlet, up->kernels->srcker);

  // Select kernel that fetches the solution:
  up->kernels->solker = fem_parproj_choose_solstencil_kernel(basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_parproj_choose_kernels_cu(basis, up->isperiodic, up->isdirichlet, up->kernels_cu);
#endif

  // MF 2022/08/23: at the moment we only support weight=/1 for cdim=1. For
  // cdim=3 we need to create a separate lhs A matrix for row of z cells.
  if (isweighted) assert(up->ndim==1);

  // Create a linear Ax=B problem. We envision two cases:
  //  a) No weight, or weight is a scalar so we can divide the RHS by it. Then
  //     A is the same for every problem, and we can just populate B with a
  //     column for each problem.
  //  b) Weight depends on space. Then we have to create an A matrix and a
  //     separate Ax=B problem for each perpendicular cell.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->prob_cu = gkyl_culinsolver_prob_new(1, up->numnodes_global, up->numnodes_global, up->perp_range.volume);
  else
    up->prob = gkyl_superlu_prob_new(1, up->numnodes_global, up->numnodes_global, up->perp_range.volume);
#else
    up->prob = gkyl_superlu_prob_new(1, up->numnodes_global, up->numnodes_global, up->perp_range.volume);
#endif

  // Assign non-zero elements in A.
  struct gkyl_mat_triples **tri = gkyl_malloc(sizeof(struct gkyl_mat_triples *));
  tri[0] = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri[0]);
#endif
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);
//    printf("par_iter1d.idx=%d | paridx=%ld | up->parnum_cells=%d\n",up->par_iter1d.idx[0], paridx, up->parnum_cells);

    int keri = up->par_iter1d.idx[0] == up->solve_range->upper[up->pardir]? 1 : 0;
    up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

    const double *wgt_p = isweighted? gkyl_array_cfetch(weight, paridx+1) : NULL;

    // Apply the wgt*phi*basis stencil.
    keri = idx_to_inloup_ker(up->par_range1d.lower[0],
      up->par_range1d.upper[0], up->par_iter1d.idx[0]);
    up->kernels->lhsker[keri](wgt_p, up->globalidx, tri[0]);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_culinsolver_amat_from_triples(up->prob_cu, tri);
  else
    gkyl_superlu_amat_from_triples(up->prob, tri);
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

  gkyl_mat_triples_release(tri[0]);
  gkyl_free(tri);

  return up;
}

void
gkyl_fem_parproj_set_rhs(struct gkyl_fem_parproj* up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));
    if (phibc)
      assert(gkyl_array_is_cu_dev(phibc));

    gkyl_fem_parproj_set_rhs_cu(up, rhsin, phibc);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);
  int skin_idx[GKYL_MAX_CDIM] = {-1};

  gkyl_range_iter_init(&up->solve_iter, up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {

    int idx1d[] = {up->solve_iter.idx[up->pardir]};

    long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);
    const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

    for (int d=0; d<up->ndim-1; d++) skin_idx[d] = up->solve_iter.idx[d];
    skin_idx[up->pardir] = idx1d[0] == up->parnum_cells? idx1d[0] : idx1d[0];
    linidx = gkyl_range_idx(up->solve_range_ext, skin_idx);
    const double *phibc_p = up->isdirichlet? gkyl_array_cfetch(phibc, linidx) : NULL;

    long paridx = gkyl_range_idx(&up->par_range1d, idx1d);
    int keri = idx1d[0] == up->par_range1d.upper[0]? 1 : 0;
    up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

    int idx2d[] = {up->perp_range2d.lower[0], up->perp_range2d.lower[0]};
    for (int d=0; d<up->ndim-1; d++) idx2d[d] = up->solve_iter.idx[d];
    long perpidx2d = gkyl_range_idx(&up->perp_range2d, idx2d);
    long perpProbOff = perpidx2d*up->numnodes_global;

    keri = idx_to_inloup_ker(up->par_range1d.lower[0],
      up->par_range1d.upper[0], idx1d[0]);
    up->kernels->srcker[keri](rhsin_p, phibc_p, perpProbOff, up->globalidx, brhs_p);

  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_parproj_solve(struct gkyl_fem_parproj* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_parproj_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  gkyl_range_iter_init(&up->solve_iter, up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {

    long linidx = gkyl_range_idx(up->solve_range, up->solve_iter.idx);
    double *phiout_p = gkyl_array_fetch(phiout, linidx);

    int idx1d[] = {up->solve_iter.idx[up->ndim-1]};
    long paridx = gkyl_range_idx(&up->par_range1d, idx1d);
    int keri = idx1d[0] == up->par_range1d.upper[0]? 1 : 0;
    up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

    int idx2d[] = {up->perp_range2d.lower[0], up->perp_range2d.lower[0]};
    for (int d=0; d<up->ndim-1; d++) idx2d[d] = up->solve_iter.idx[d];
    long perpidx2d = gkyl_range_idx(&up->perp_range2d, idx2d);
    long perpProbOff = perpidx2d*up->numnodes_global;

    up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), perpProbOff, up->globalidx, phiout_p);

  }

}

void gkyl_fem_parproj_release(struct gkyl_fem_parproj *up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
    gkyl_culinsolver_prob_release(up->prob_cu);
  } else {
    gkyl_superlu_prob_release(up->prob);
  }
#else
  gkyl_superlu_prob_release(up->prob);
#endif
  gkyl_array_release(up->brhs);
  gkyl_free(up->globalidx);
  gkyl_free(up->kernels);
  gkyl_free(up);
}
