#include <gkyl_fem_poisson.h>
#include <gkyl_fem_poisson_priv.h>

static long
global_num_nodes(const int dim, const int poly_order, const int basis_type, const int *num_cells)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_poisson_num_nodes_global_1x_ser_p1_nonperiodicx(num_cells);
    } else if (poly_order == 2) {
      return fem_poisson_num_nodes_global_1x_ser_p2_nonperiodicx(num_cells);
//    } else if (poly_order == 3) {
//      return fem_poisson_num_nodes_global_1x_ser_p3_nonperiodicx(num_cells);
    }
  }
  assert(false);  // Other dimensionalities not supported.
  return -1;
}

static void
local_stiff(const int dim, const int poly_order, const int basis_type, const double *dx, struct gkyl_mat *stiffout)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_poisson_stiff_1x_ser_p1(dx,stiffout);
    } else if (poly_order == 2) {
      return fem_poisson_stiff_1x_ser_p2(dx,stiffout);
//    } else if (poly_order == 3) {
//      return fem_poisson_stiff_1x_ser_p3(dx,stiffout);
    }
  } else if (dim==2) {
    if (poly_order == 1) {
      return fem_poisson_stiff_2x_ser_p1(dx,stiffout);
    } else if (poly_order == 2) {
      return fem_poisson_stiff_2x_ser_p2(dx,stiffout);
//    } else if (poly_order == 3) {
//      return fem_poisson_stiff_2x_ser_p3(dx,stiffout);
    }
//  } else if (dim==3) {
//    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
//      if (poly_order == 1) {
//        return fem_poisson_stiff_3x_ser_p1(dx,stiffout);
//      } else if (poly_order == 2) {
//        return fem_poisson_stiff_3x_ser_p2(dx,stiffout);
//      } else if (poly_order == 3) {
//        return fem_poisson_stiff_3x_ser_p3(dx,stiffout);
//      }
//    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
//      if (poly_order == 1) {
//        return fem_poisson_stiff_3x_tensor_p1(dx,stiffout);
//      } else if (poly_order == 2) {
//        return fem_poisson_stiff_3x_tensor_p2(dx,stiffout);
//      }
//    }
  }
  assert(false);  // Other dimensionalities not supported.
}

static void
local_mass_modtonod(const int dim, const int poly_order, const int basis_type, struct gkyl_mat *mass_mod2nod)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_poisson_mass_times_modtonod_1x_ser_p1(mass_mod2nod);
    } else if (poly_order == 2) {
      return fem_poisson_mass_times_modtonod_1x_ser_p2(mass_mod2nod);
//    } else if (poly_order == 3) {
//      return fem_poisson_mass_times_modtonod_1x_ser_p3(mass_mod2nod);
    }
  } else if (dim==2) {
    if (poly_order == 1) {
      return fem_poisson_mass_times_modtonod_2x_ser_p1(mass_mod2nod);
    } else if (poly_order == 2) {
      return fem_poisson_mass_times_modtonod_2x_ser_p2(mass_mod2nod);
//    } else if (poly_order == 3) {
//      return fem_poisson_mass_times_modtonod_2x_ser_p3(mass_mod2nod);
    }
//  } else if (dim==3) {
//    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
//      if (poly_order == 1) {
//        return fem_poisson_mass_times_modtonod_3x_ser_p1(mass_mod2nod);
//      } else if (poly_order == 2) {
//        return fem_poisson_mass_times_modtonod_3x_ser_p2(mass_mod2nod);
//      } else if (poly_order == 3) {
//        return fem_poisson_mass_times_modtonod_3x_ser_p3(mass_mod2nod);
//      }
//    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
//      if (poly_order == 1) {
//        return fem_poisson_mass_times_modtonod_3x_tensor_p1(mass_mod2nod);
//      } else if (poly_order == 2) {
//        return fem_poisson_mass_times_modtonod_3x_tensor_p2(mass_mod2nod);
//      }
//    }
  }
  assert(false);  // Other dimensionalities not supported.
}

static void
local_nodtomod(const int dim, const int poly_order, const int basis_type, struct gkyl_mat *nod2mod)
{
  if (dim==1) {
    if (poly_order == 1) {
      return fem_poisson_nodtomod_1x_ser_p1(nod2mod);
    } else if (poly_order == 2) {
      return fem_poisson_nodtomod_1x_ser_p2(nod2mod);
//    } else if (poly_order == 3) {
//      return fem_poisson_nodtomod_1x_ser_p3(nod2mod);
    }
  } else if (dim==2) {
    if (poly_order == 1) {
      return fem_poisson_nodtomod_2x_ser_p1(nod2mod);
    } else if (poly_order == 2) {
      return fem_poisson_nodtomod_2x_ser_p2(nod2mod);
//    } else if (poly_order == 3) {
//      return fem_poisson_nodtomod_2x_ser_p3(nod2mod);
    }
//  } else if (dim==3) {
//    if (basis_type == GKYL_BASIS_MODAL_SERENDIPITY) {
//      if (poly_order == 1) {
//        return fem_poisson_nodtomod_3x_ser_p1(nod2mod);
//      } else if (poly_order == 2) {
//        return fem_poisson_nodtomod_3x_ser_p2(nod2mod);
//      } else if (poly_order == 3) {
//        return fem_poisson_nodtomod_3x_ser_p3(nod2mod);
//      }
//    } if (basis_type == GKYL_BASIS_MODAL_TENSOR) {
//      if (poly_order == 1) {
//        return fem_poisson_nodtomod_3x_tensor_p1(nod2mod);
//      } else if (poly_order == 2) {
//        return fem_poisson_nodtomod_3x_tensor_p2(nod2mod);
//      }
//    }
  }
  assert(false);  // Other dimensionalities not supported.
}

int idx_to_inup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the in-up kernel given the grid index.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == num_cells[d]) iout += (int)(pow(2,d)+0.5);
  }
  return iout;
}

gkyl_fem_poisson*
gkyl_fem_poisson_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis basis,
  struct gkyl_poisson_bc bcs, const double epsilon, void *ctx)
{

  gkyl_fem_poisson *up = gkyl_malloc(sizeof(gkyl_fem_poisson));

  up->ctx = ctx;
  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis =  basis.num_basis;
  up->basis_type = basis.b_type;
  up->poly_order = basis.poly_order;
  up->basis = basis;

  for (int d=0; d<up->ndim; d++) up->isdirperiodic[d] = bcs.lo_type[d] == GKYL_POISSON_PERIODIC;

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis])); // global index, one for each basis in a cell.

  // Local and local-ext ranges for whole-grid arrays.
  int ghost[POISSON_MAX_DIM];
  for (int d=0; d<up->ndim; d++) ghost[d] = 1;
  gkyl_create_grid_ranges(grid, ghost, &up->local_range_ext, &up->local_range);
  // Range of cells we'll solve Poisson in, as
  // a sub-range of up->local_range_ext.
  int sublower[POISSON_MAX_DIM], subupper[POISSON_MAX_DIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->local_range.lower[d];
    subupper[d] = up->local_range.upper[d];
  }
  gkyl_sub_range_init(&up->solve_range, &up->local_range_ext, sublower, subupper);
  for (int d=0; d<up->ndim; d++) up->num_cells[d] = up->solve_range.upper[d]-up->solve_range.lower[d]+1;

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = global_num_nodes(up->ndim, up->poly_order, basis.b_type, &up->num_cells[0]);

  // Prepare for periodic domain case.
  up->isdomperiodic = true;
  for (int d=0; d<up->ndim; d++) up->isdomperiodic = up->isdomperiodic && up->isdirperiodic[d];
  if (up->isdomperiodic) {
    up->rhs_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, up->solve_range.volume);
    gkyl_array_clear(up->rhs_cellavg, 0.0);
    up->mavgfac = -pow(sqrt(2.),up->ndim);
  }

  // Create local matrices used later.
  double dx[POISSON_MAX_DIM];
  for (int d=0; d<up->ndim; d++) dx[d] = up->grid.dx[d];
  up->local_stiff = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  up->local_mass_modtonod = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  up->local_nodtomod = gkyl_mat_new(up->numnodes_local, up->numnodes_local, 0.);
  local_stiff(up->ndim, up->poly_order, basis.b_type, &dx[0], up->local_stiff);
  local_mass_modtonod(up->ndim, up->poly_order, basis.b_type, up->local_mass_modtonod);
  local_nodtomod(up->ndim, up->poly_order, basis.b_type, up->local_nodtomod);

  // Select local-to-global mapping kernels:
  choose_local2global_kernels(&basis, &up->isdirperiodic[0], &up->l2g[0]);

  // Create a linear Ax=B problem. Here A is the discrete (global) stiffness
  // matrix times -epsilon.
  up->prob = gkyl_superlu_prob_new(up->numnodes_global, up->numnodes_global, 1);

  // Assign non-zero elements in A.
  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long idx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);

    int keri = idx_to_inup_ker(up->ndim, &up->num_cells[0], up->solve_iter.idx);
    up->l2g[keri](&up->num_cells[0], &up->solve_iter.idx[0], &up->globalidx[0]);
    for (size_t k=0; k<up->numnodes_local; ++k) {
      for (size_t m=0; m<up->numnodes_local; ++m) {
        double val = -epsilon*gkyl_mat_get(up->local_stiff,k,m);
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
gkyl_fem_poisson_set_rhs(gkyl_fem_poisson* up, const struct gkyl_array *rhsin)
{

  if (up->isdomperiodic) {
    // Subtract the volume averaged RHS.
    gkyl_dg_calc_average_range(up->basis, 0, up->rhs_cellavg, 0, rhsin, up->solve_range);
    gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, up->solve_range);
//    gkyl_array_shiftc0(rhsin, up->mavgfac*up->rhs_avg);
  }

  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->numnodes_global, 1);

  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);

    const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

    int keri = idx_to_inup_ker(up->ndim, &up->num_cells[0], up->solve_iter.idx);
    up->l2g[keri](&up->num_cells[0], &up->solve_iter.idx[0], &up->globalidx[0]);

    for (size_t k=0; k<up->numnodes_local; k++) {
      long globalidx_k = up->globalidx[k];
      double massnod_rhs = 0.;
      for (size_t m=0; m<up->numnodes_local; m++)
        massnod_rhs += gkyl_mat_get(up->local_mass_modtonod,k,m) * rhsin_p[m];
      gkyl_mat_triples_accum(tri, globalidx_k, 0, massnod_rhs);
    }
  }

  gkyl_superlu_brhs_from_triples(up->prob, tri);
  gkyl_mat_triples_release(tri);

}


void
gkyl_fem_poisson_solve(gkyl_fem_poisson* up, struct gkyl_array *phiout) {
  gkyl_superlu_solve(up->prob);

  int pidx[POISSON_MAX_DIM] = {0};

  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);

    double *phiout_p = gkyl_array_fetch(phiout, linidx);

    int keri = idx_to_inup_ker(up->ndim, &up->num_cells[0], up->solve_iter.idx);
    up->l2g[keri](&up->num_cells[0], &up->solve_iter.idx[0], &up->globalidx[0]);

    for (size_t k=0; k<up->numnodes_local; k++) {
      phiout_p[k] = 0.;
      for (size_t m=0; m<up->numnodes_local; m++) {
        long globalidx_m = up->globalidx[m];
        phiout_p[k] += gkyl_mat_get(up->local_nodtomod,k,m) * gkyl_superlu_get_rhs_ij(up->prob, globalidx_m, 0);
      }
    }
  }

}

void gkyl_fem_poisson_release(gkyl_fem_poisson *up)
{
  if (up->isdomperiodic) gkyl_array_release(up->rhs_cellavg);
  gkyl_mat_release(up->local_stiff);
  gkyl_mat_release(up->local_mass_modtonod);
  gkyl_mat_release(up->local_nodtomod);
  gkyl_superlu_prob_release(up->prob);
  gkyl_free(up->globalidx);
  gkyl_free(up);
}
