#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_kernels.h>

struct gkyl_fem_parproj {
  void *ctx; // evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int pardir; // parallel (z) direction.
  int parnum_cells; // number of cells in parallel (z) direction.
  bool isperiodic; // =true if parallel direction is periodic.

  struct gkyl_range local_range, local_range_ext;
  struct gkyl_range perp_range; // range of perpendicular cells.
  struct gkyl_range_iter perp_iter;
  struct gkyl_range par_range; // range of parallel cells.
  struct gkyl_range_iter par_iter;
  struct gkyl_range par_range1d; // 1D range of parallel cells.
  struct gkyl_range_iter par_iter1d;
  struct gkyl_range perp_range2d; // 2D range of perpendicular cells.
  struct gkyl_range_iter perp_iter2d;

  int bc_off;
  // These are used with periodic BCs:
  struct gkyl_array *parbc_buff;
  struct gkyl_range parskin_lo, parskin_up, parghost_lo, parghost_up;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_superlu_prob* prob;
  struct gkyl_mat *local_mass; // local mass matrix.
  struct gkyl_mat *local_mass_modtonod; // local mass matrix times modal-to-nodal matrix.
  struct gkyl_mat *local_nodtomod; // local nodal-to-modal matrix.

  long *globalidx;
};

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

static void
local_mass_modtonod(const int dim, const int poly_order, const int basis_type, struct gkyl_mat *mass_mod2nod)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_parproj_mass_times_modtonod_1x_ser_p1(mass_mod2nod);
    } else if (poly_order == 2) {
      return fem_parproj_mass_times_modtonod_1x_ser_p2(mass_mod2nod);
    } else if (poly_order == 3) {
      return fem_parproj_mass_times_modtonod_1x_ser_p3(mass_mod2nod);
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return fem_parproj_mass_times_modtonod_3x_ser_p1(mass_mod2nod);
      } else if (poly_order == 2) {
        return fem_parproj_mass_times_modtonod_3x_ser_p2(mass_mod2nod);
      } else if (poly_order == 3) {
        return fem_parproj_mass_times_modtonod_3x_ser_p3(mass_mod2nod);
      }
    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
      if (poly_order == 1) {
        return fem_parproj_mass_times_modtonod_3x_tensor_p1(mass_mod2nod);
      } else if (poly_order == 2) {
        return fem_parproj_mass_times_modtonod_3x_tensor_p2(mass_mod2nod);
      }
    }
  }
  assert(false);  // Other dimensionalities not supported.
}

static void
local_nodtomod(const int dim, const int poly_order, const int basis_type, struct gkyl_mat *nod2mod)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_parproj_nodtomod_1x_ser_p1(nod2mod);
    } else if (poly_order == 2) {
      return fem_parproj_nodtomod_1x_ser_p2(nod2mod);
    } else if (poly_order == 3) {
      return fem_parproj_nodtomod_1x_ser_p3(nod2mod);
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return fem_parproj_nodtomod_3x_ser_p1(nod2mod);
      } else if (poly_order == 2) {
        return fem_parproj_nodtomod_3x_ser_p2(nod2mod);
      } else if (poly_order == 3) {
        return fem_parproj_nodtomod_3x_ser_p3(nod2mod);
      }
    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
      if (poly_order == 1) {
        return fem_parproj_nodtomod_3x_tensor_p1(nod2mod);
      } else if (poly_order == 2) {
        return fem_parproj_nodtomod_3x_tensor_p2(nod2mod);
      }
    }
  }
  assert(false);  // Other dimensionalities not supported.
}

static void
local_to_global(const int dim, const int poly_order, const int basis_type, const int parnum_cells, const int paridx, long *globalidx)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_parproj_local_to_global_1x_ser_p1(parnum_cells, paridx, globalidx);
    } else if (poly_order == 2) {
      return fem_parproj_local_to_global_1x_ser_p2(parnum_cells, paridx, globalidx);
    } else if (poly_order == 3) {
      return fem_parproj_local_to_global_1x_ser_p3(parnum_cells, paridx, globalidx);
    }
  } else if (dim==3) {
    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      if (poly_order == 1) {
        return fem_parproj_local_to_global_3x_ser_p1(parnum_cells, paridx, globalidx);
      } else if (poly_order == 2) {
        return fem_parproj_local_to_global_3x_ser_p2(parnum_cells, paridx, globalidx);
      } else if (poly_order == 3) {
        return fem_parproj_local_to_global_3x_ser_p3(parnum_cells, paridx, globalidx);
      }
    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
      if (poly_order == 1) {
        return fem_parproj_local_to_global_3x_tensor_p1(parnum_cells, paridx, globalidx);
      } else if (poly_order == 2) {
        return fem_parproj_local_to_global_3x_tensor_p2(parnum_cells, paridx, globalidx);
      }
    }
  }
  assert(false);  // Other dimensionalities not supported.
}

gkyl_fem_parproj*
gkyl_fem_parproj_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const bool isparperiodic, void *ctx)
{
  gkyl_fem_parproj *up = gkyl_malloc(sizeof(gkyl_fem_parproj));

  up->ctx = ctx;
  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis = basis->num_basis;
  up->basis_type = basis->b_type;
  up->poly_order = basis->poly_order;
  up->pardir = grid->ndim-1; // Assume parallel direction is always the last.
  up->isperiodic = isparperiodic;

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis]));

  // Local and local-ext ranges for whole-grid arrays.
  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<up->ndim; d++) ghost[d] = 1;
  gkyl_create_grid_ranges(grid, ghost, &up->local_range_ext, &up->local_range);
  // Range of parallel cells, as a sub-range of up->local_range.
  int sublower[GKYL_MAX_DIM], subupper[GKYL_MAX_DIM];
  for (int d=0; d<up->ndim-1; d++) {sublower[d] = up->local_range.lower[d];  subupper[d] = up->local_range.lower[d];}
  sublower[up->pardir] = up->local_range.lower[up->pardir];
  subupper[up->pardir] = up->local_range.upper[up->pardir];
  up->bc_off = 1;
  if (up->isperiodic) {
    // Include ghost cells in parallel direction.
    sublower[up->pardir] = up->local_range_ext.lower[up->pardir];
    subupper[up->pardir] = up->local_range_ext.upper[up->pardir];
    // Will also need an offset to convert indices in the grid to indices in
    // the linear problem.
    up->bc_off = 0;
  }
  gkyl_sub_range_init(&up->par_range, &up->local_range_ext, sublower, subupper);
  up->parnum_cells = up->par_range.volume;
  // Range of perpendicular cells.
  gkyl_range_shorten(&up->perp_range, &up->local_range, up->pardir, 1);

  // 1D range of parallel cells.
  int lower1d[] = {up->par_range.lower[up->pardir]}, upper1d[] = {up->par_range.upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);
  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim==3 ? 2 : 1, up->perp_range.lower, up->perp_range.upper);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = global_num_nodes(up->ndim, up->poly_order, basis->b_type, up->par_range.volume);

  // Create local matrices used later.
  up->local_mass = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  local_mass(up->ndim, up->poly_order, basis->b_type, up->local_mass);
  up->local_mass_modtonod = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  local_mass_modtonod(up->ndim, up->poly_order, basis->b_type, up->local_mass_modtonod);
  up->local_nodtomod = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  local_nodtomod(up->ndim, up->poly_order, basis->b_type, up->local_nodtomod);

  // Create a linear Ax=B problem. We envision two cases:
  //  a) No weight, or weight is a scalar so we can divide the RHS by it. Then
  //     A is the same for every problem, and we can just populate B with a
  //     column for each problem.
  //  b) Weight depends on space. Then we have to create an A matrix and a
  //     separate Ax=B problem for each perpendicular cell.
  // For now we restrict ourselves to a).
  up->prob = gkyl_superlu_prob_new(up->numnodes_global, up->numnodes_global, up->perp_range.volume);

  // Assign non-zero elements in A.
  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
  gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
  while (gkyl_range_iter_next(&up->par_iter1d)) {
    long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

    local_to_global(up->ndim, up->poly_order, up->basis_type, up->parnum_cells, paridx, up->globalidx);
    for (size_t k=0; k<up->numnodes_local; ++k) {
      for (size_t m=0; m<up->numnodes_local; ++m) {
        double val = gkyl_mat_get(up->local_mass,k,m);
        long globalidx_k = up->globalidx[k];
        long globalidx_m = up->globalidx[m];
        gkyl_mat_triples_accum(tri, globalidx_k, globalidx_m, val);
      }
    }
  }
  gkyl_superlu_amat_from_triples(up->prob, tri);

  gkyl_mat_triples_release(tri);

  return up;
}

void
gkyl_fem_parproj_set_rhs(gkyl_fem_parproj* up, const struct gkyl_array *rhsin)
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->numnodes_global, up->perp_range.volume);

  int cidx[GKYL_MAX_DIM] = {0};

  gkyl_range_iter_init(&up->perp_iter, &up->perp_range);
  while (gkyl_range_iter_next(&up->perp_iter)) {
    long perpidx = gkyl_range_idx(&up->perp_range, up->perp_iter.idx);
    long perpidx2d = gkyl_range_idx(&up->perp_range2d, up->perp_iter.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (int d=0; d<up->ndim-1; d++) cidx[d] = up->perp_iter.idx[d];
      cidx[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(&up->local_range, cidx);
      const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

      local_to_global(up->ndim, up->poly_order, up->basis_type, up->parnum_cells, paridx, up->globalidx);

      for (size_t k=0; k<up->numnodes_local; k++) {
        long globalidx_k = up->globalidx[k];
        double massnod_rhs = 0.;
        for (size_t m=0; m<up->numnodes_local; m++)
          massnod_rhs += gkyl_mat_get(up->local_mass_modtonod,k,m) * rhsin_p[m];
        gkyl_mat_triples_accum(tri, globalidx_k, perpidx2d, massnod_rhs);
      }
    }

  }

  gkyl_superlu_brhs_from_triples(up->prob, tri);
  gkyl_mat_triples_release(tri);

}

void
gkyl_fem_parproj_solve(gkyl_fem_parproj* up, struct gkyl_array *phiout) {
  gkyl_superlu_solve(up->prob);

  int pidx[GKYL_MAX_DIM] = {0};

  gkyl_range_iter_init(&up->perp_iter, &up->perp_range);
  while (gkyl_range_iter_next(&up->perp_iter)) {
    long perpidx = gkyl_range_idx(&up->perp_range, up->perp_iter.idx);
    long perpidx2d = gkyl_range_idx(&up->perp_range2d, up->perp_iter.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (int d=0; d<up->ndim-1; d++) pidx[d] = up->perp_iter.idx[d];
      pidx[up->pardir] = up->par_iter1d.idx[0];

      long linidx = gkyl_range_idx(&up->local_range, pidx);
      double *phiout_p = gkyl_array_fetch(phiout, linidx);

      local_to_global(up->ndim, up->poly_order, up->basis_type, up->parnum_cells, paridx, up->globalidx);

      for (size_t k=0; k<up->numnodes_local; k++) {
        phiout_p[k] = 0.;
        for (size_t m=0; m<up->numnodes_local; m++) {
          long globalidx_m = up->globalidx[m];
          phiout_p[k] += gkyl_mat_get(up->local_nodtomod,k,m) * gkyl_superlu_get_rhs_ij(up->prob, globalidx_m, perpidx2d);
        }
      }

    }
  }

}

void gkyl_fem_parproj_release(gkyl_fem_parproj *up)
{
  if (up->isperiodic) {
    gkyl_array_release(up->parbc_buff);
  }
  gkyl_mat_release(up->local_mass);
  gkyl_mat_release(up->local_mass_modtonod);
  gkyl_mat_release(up->local_nodtomod);
  gkyl_superlu_prob_release(up->prob);
  gkyl_free(up->globalidx);
  gkyl_free(up);
}
