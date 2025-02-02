#include <gkyl_fem_parproj_multib.h>
#include <gkyl_fem_parproj_priv.h>
#include <gkyl_fem_parproj_multib_priv.h>

struct gkyl_fem_parproj_multib*
gkyl_fem_parproj_multib_new(int num_blocks, const struct gkyl_range *mbz_range,
  struct gkyl_range *solve_ranges, const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right, bool use_gpu)
{
  struct gkyl_fem_parproj_multib *up = gkyl_malloc(sizeof(*up));

  up->num_blocks = num_blocks;
  up->mbz_range = mbz_range;
  up->ndim = basis->ndim;
  up->num_basis  = basis->num_basis;
  up->basis_type = basis->b_type;
  up->poly_order = basis->poly_order;
  up->pardir = up->ndim-1; // Assume parallel direction is always the last.
  up->isdirichlet = bctype == GKYL_FEM_PARPROJ_DIRICHLET;
  up->use_gpu = use_gpu;

  // Periodic BCs are not allowed.
  assert(!(bctype == GKYL_FEM_PARPROJ_PERIODIC));

  // Must always pass a weight (e.g. the Jacobian).
  assert(weight_left && weight_right);
  up->weight_rhs = gkyl_array_acquire(weight_right);

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis]));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_parproj_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_parproj_kernels));
  else
    up->kernels_cu = up->kernels;
#endif

  // Select local-to-global mapping kernel:
  fem_parproj_choose_local2global_kernel(basis, false, up->kernels->l2g);

  // Select weighted LHS kernel (not always used):
  fem_parproj_multib_choose_lhs_kernel(basis, up->isdirichlet, true, up->kernels->lhsker);

  // Select RHS source kernel:
  fem_parproj_multib_choose_srcstencil_kernel(basis, up->isdirichlet, true, up->kernels->srcker);

  // Select kernel that fetches the solution:
  up->kernels->solker = fem_parproj_choose_solstencil_kernel(basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_parproj_choose_kernels_cu(basis, false, up->isdirichlet, up->kernels_cu);
#endif

  up->solve_range = gkyl_malloc(num_blocks*sizeof(struct gkyl_range));
  for (int bI=0; bI<num_blocks; bI++)
    up->solve_range[bI] = solve_ranges[bI];

  // Range of perpendicular cells (assumed equal for every block).
  const int *lower2d = up->solve_range[0].lower;
  const int *upper2d = up->ndim==1? up->solve_range[0].lower : up->solve_range[0].upper;
  gkyl_range_init(&up->perp_range2d, up->ndim==3 ? 2 : 1, lower2d, upper2d);

  // Number of nodes on perpendicular plane between two blocks.
  assert(up->poly_order=1);
  up->numnodes_perp = pow(2, up->ndim-1); // Using the p=1 assumption here.

  up->par_range1d = gkyl_malloc(num_blocks*sizeof(struct gkyl_range));
  up->parnum_cells = gkyl_malloc(num_blocks*sizeof(int));
  up->numnodes_global = gkyl_malloc(num_blocks*sizeof(long));
  long numnodes_global_tot = 0;
  for (int bI=0; bI<num_blocks; bI++) {
    // 1D range of parallel cells.
    int lower1d[] = {up->solve_range[bI].lower[up->pardir]}, upper1d[] = {up->solve_range[bI].upper[up->pardir]};
    gkyl_range_init(&up->par_range1d[bI], 1, lower1d, upper1d);
  
    up->parnum_cells[bI] = up->par_range1d[bI].volume;
  
    // Compute the number of local and global nodes.
    up->numnodes_global[bI] = gkyl_fem_parproj_global_num_nodes(basis, false, up->par_range1d[bI].volume);

    numnodes_global_tot += up->numnodes_global[bI];
  }
  // Add additional "nodes" (equations really) to impose equality between
  // colocated nodes from neighboring blocks.
  numnodes_global_tot += (num_blocks-1)*up->numnodes_perp;

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, numnodes_global_tot*up->perp_range2d.volume); // Global right side vector.

  // We assume there's a spatially dependent weight. Then we solve A_i x_i=B_i
  // where there's a different A_i for each perp cell and B_i is a single
  // column matrix.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->prob_cu = gkyl_culinsolver_prob_new(up->perp_range2d.volume, numnodes_global_tot, numnodes_global_tot, 1);
  else
    up->prob = gkyl_superlu_prob_new(up->perp_range2d.volume, numnodes_global_tot, numnodes_global_tot, 1);
#else
  up->prob = gkyl_superlu_prob_new(up->perp_range2d.volume, numnodes_global_tot, numnodes_global_tot, 1);
#endif

  // Assign non-zero elements in A.
  struct gkyl_mat_triples **tri = gkyl_malloc(up->perp_range2d.volume*sizeof(struct gkyl_mat_triples *));
  for (size_t i=0; i<up->perp_range2d.volume; i++) {
    tri[i] = gkyl_mat_triples_new(numnodes_global_tot, numnodes_global_tot);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri[i]);
#endif
  }

  int idx1[GKYL_MAX_CDIM];

  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    long block_offset = 0;
    int par_idx_offset = 0;
    for (int bI=0; bI<num_blocks; bI++) {

      gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d[bI]);
      while (gkyl_range_iter_next(&up->par_iter1d)) {
        long paridx = gkyl_range_idx(&up->par_range1d[bI], up->par_iter1d.idx);

        for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = par_idx_offset + up->par_iter1d.idx[0];
        long linidx = gkyl_range_idx(up->mbz_range, idx1);
        const double *wgt_p = gkyl_array_cfetch(weight_left, linidx);

        int keri = up->par_iter1d.idx[0] == up->parnum_cells[bI]? 1 : 0;
        up->kernels->l2g[keri](up->parnum_cells[bI], paridx, up->globalidx);
        // Shift global indices to account for previous blocks.
        for (int k=0; k<up->num_basis; k++) {
          up->globalidx[k] += block_offset;
        }

        // Apply the wgt*phi*basis stencil.
        keri = idx_to_inloup_ker(up->parnum_cells[bI], up->par_iter1d.idx[0]);
        const double *phibc_p = NULL;
        if (up->isdirichlet) {
          // Interior blocks don't apply  Dirichlet BCs.
          if (bI > 0 && keri == 1) keri = 3;
          if (bI < up->num_blocks-1 && keri == 2) keri = 4;
        }
        up->kernels->lhsker[keri](wgt_p, up->globalidx, tri[perpidx]);
      }
      block_offset += up->numnodes_global[bI];

      if (bI < num_blocks-1) {
        // Insert additional conditions setting upper (along z) nodes of this block
        // equal to the lower nodes of the next block.

        long globalidx_lo[up->num_basis];
        idx1[0] = up->par_range1d[bI].upper[0];
        long paridx_lo = gkyl_range_idx(&up->par_range1d[bI], idx1);
        int keri_lo = 1;
        up->kernels->l2g[keri_lo](up->parnum_cells[bI], paridx_lo, globalidx_lo);

        long globalidx_up[up->num_basis];
        idx1[0] = up->par_range1d[bI+1].lower[0];
        long paridx_up = gkyl_range_idx(&up->par_range1d[bI+1], idx1);
        int keri_up = 0;
        up->kernels->l2g[keri_up](up->parnum_cells[bI+1], paridx_up, globalidx_up);

        int local_off = up->numnodes_perp; // Using the p=1 assumption here.
        for (int k=0; k<up->numnodes_perp; k++) {
          // =1 entry for node on lower block.
          long ilo = block_offset+k;
          long jlo = block_offset-up->numnodes_global[bI]+globalidx_lo[local_off+k];
          gkyl_mat_triples_insert(tri[perpidx], ilo, jlo, 1.0);
          gkyl_mat_triples_insert(tri[perpidx], jlo, ilo, 1.0);

          // =-1 entry for node on lower block (corresponding RHS entry =0).
          long iup = block_offset+k;
          long jup = block_offset+up->numnodes_perp+globalidx_up[k];
          gkyl_mat_triples_insert(tri[perpidx], iup, jup, -1.0);
          gkyl_mat_triples_insert(tri[perpidx], jup, iup, -1.0);
        }
        block_offset += up->numnodes_perp;
      }

      par_idx_offset += up->par_range1d[bI].upper[0]-up->par_range1d[bI].lower[0]+1;
    }
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_culinsolver_amat_from_triples(up->prob_cu, tri);
  else
    gkyl_superlu_amat_from_triples(up->prob, tri);
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

  for (size_t i=0; i<up->perp_range2d.volume; i++)
    gkyl_mat_triples_release(tri[i]);
  gkyl_free(tri);

  return up;
}

void
gkyl_fem_parproj_multib_set_rhs(struct gkyl_fem_parproj_multib* up,
  const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));
    if (phibc)
      assert(gkyl_array_is_cu_dev(phibc));

    gkyl_fem_parproj_multib_set_rhs_cu(up, rhsin, phibc);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx1[GKYL_MAX_CDIM];
  long prob_offset = 0;

  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    int par_idx_offset = 0;
    for (int bI=0; bI<up->num_blocks; bI++) {

      gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d[bI]);
      while (gkyl_range_iter_next(&up->par_iter1d)) {
        long paridx = gkyl_range_idx(&up->par_range1d[bI], up->par_iter1d.idx);

        for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = par_idx_offset + up->par_iter1d.idx[0];
        long linidx = gkyl_range_idx(up->mbz_range, idx1);

        const double *wgt_p = gkyl_array_cfetch(up->weight_rhs, linidx);
        const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

        int keri = up->par_iter1d.idx[0] == up->parnum_cells[bI]? 1 : 0;
        up->kernels->l2g[keri](up->parnum_cells[bI], paridx, up->globalidx);

        // Shift global indices by this offset to account for linear
        // problems at other perp cells and for previous blocks.

        keri = idx_to_inloup_ker(up->parnum_cells[bI], up->par_iter1d.idx[0]);
        const double *phibc_p = NULL;
        if (up->isdirichlet) {
          phibc_p = gkyl_array_cfetch(phibc, linidx);
          // Interior blocks don't apply  Dirichlet BCs.
          if (bI > 0 && keri == 1) keri = 3;
          if (bI < up->num_blocks-1 && keri == 2) keri = 4;
        }
        up->kernels->srcker[keri](wgt_p, rhsin_p, phibc_p, prob_offset, up->globalidx, brhs_p);
      }

      par_idx_offset += up->par_range1d[bI].upper[0]-up->par_range1d[bI].lower[0]+1;

      prob_offset += up->numnodes_global[bI];
      if (bI < up->num_blocks-1)
        prob_offset += up->numnodes_perp;
    }
  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_parproj_multib_solve(struct gkyl_fem_parproj_multib* up, struct gkyl_array *phiout)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_parproj_multib_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  int idx1[GKYL_MAX_CDIM];
  long prob_offset = 0;

  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    int par_idx_offset = 0;
    for (int bI=0; bI<up->num_blocks; bI++) {

      gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d[bI]);
      while (gkyl_range_iter_next(&up->par_iter1d)) {
        long paridx = gkyl_range_idx(&up->par_range1d[bI], up->par_iter1d.idx);

        for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = par_idx_offset + up->par_iter1d.idx[0];
        long linidx = gkyl_range_idx(up->mbz_range, idx1);

        double *phiout_p = gkyl_array_fetch(phiout, linidx);

        int keri = up->par_iter1d.idx[0] == up->parnum_cells[bI]? 1 : 0;
        up->kernels->l2g[keri](up->parnum_cells[bI], paridx, up->globalidx);

        up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), prob_offset, up->globalidx, phiout_p);
      }

      par_idx_offset += up->par_range1d[bI].upper[0]-up->par_range1d[bI].lower[0]+1;

      prob_offset += up->numnodes_global[bI];
      if (bI < up->num_blocks-1)
        prob_offset += up->numnodes_perp;
    }
  }

}

void gkyl_fem_parproj_multib_release(struct gkyl_fem_parproj_multib *up)
{
  gkyl_array_release(up->weight_rhs);
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
  gkyl_free(up->solve_range);
  gkyl_free(up->par_range1d);
  gkyl_free(up->parnum_cells);
  gkyl_free(up->numnodes_global);
  gkyl_free(up);
}
