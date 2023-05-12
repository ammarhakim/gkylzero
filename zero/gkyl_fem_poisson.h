#pragma once

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array_ops.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_dg_bin_ops.h>

#ifndef POISSON_MAX_DIM
# define POISSON_MAX_DIM 3
#endif

// Object type
typedef struct gkyl_fem_poisson gkyl_fem_poisson;

// Boundary condition types.
enum gkyl_poisson_bc_type {
  GKYL_POISSON_PERIODIC=0,
  GKYL_POISSON_DIRICHLET,  // sets the value. 
  GKYL_POISSON_NEUMANN,  // sets the slope normal to the boundary.
  GKYL_POISSON_ROBIN,  // a combination of dirichlet and neumann.
};

// Boundary condition values. Dirichlet and Neumann use only one value,
// Robin uses 3, and periodic ignores the value.
struct gkyl_poisson_bc_value { double v[3]; };

struct gkyl_poisson_bc {
  enum gkyl_poisson_bc_type lo_type[POISSON_MAX_DIM], up_type[POISSON_MAX_DIM];
  struct gkyl_poisson_bc_value lo_value[POISSON_MAX_DIM], up_value[POISSON_MAX_DIM];
};

/**
 * Create new updater to solve the Helmholtz problem
 *   - nabla . (epsilon * nabla phi) - kSq * phi = rho
 * using a FEM to ensure phi is continuous. This solver is also
 * used as a Poisson solver by passing a zero kSq. The input is the
 * DG field rho, which is translated to FEM. The output is the
 * DG field phi, after we've translated the FEM solution to DG.
 * Free using gkyl_fem_poisson_release method.
 *
 * For scalar, spatially constant epsilon use gkyl_fem_poisson_consteps_new,
 * for tensor or heterogenous epsilon use gkyl_fem_poisson_vareps_new.
 *
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param bcs Boundary conditions.
 * @param epsilon_const Constant scalar value of the permittivity.
 * @param epsilon_var Spatially varying permittivity tensor.
 * @param kSq Squared wave number (factor multiplying phi in Helmholtz eq).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_fem_poisson* gkyl_fem_poisson_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis basis, struct gkyl_poisson_bc *bcs,
  double epsilon_const, struct gkyl_array *epsilon_var, struct gkyl_array *kSq, bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM poisson updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_poisson_set_rhs(gkyl_fem_poisson* up, struct gkyl_array *rhsin);

void gkyl_fem_poisson_set_rhs_cu(gkyl_fem_poisson *up, struct gkyl_array *rhsin);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_poisson_solve(gkyl_fem_poisson* up, struct gkyl_array *phiout);

void gkyl_fem_poisson_solve_cu(gkyl_fem_poisson* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_release(gkyl_fem_poisson *up);

GKYL_CU_DH
static inline int idx_to_inup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between upper cells and
  // elsewhere.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == num_cells[d]) iout += (int)(pow(2,d)+0.5);
  }
  return iout;
}

GKYL_CU_DH
static inline int idx_to_inloup_ker(const int dim, const int *num_cells, const int *idx) {
  // Return the index of the kernel (in the array of kernels) needed given the grid index.
  // This function is for kernels that differentiate between lower, interior
  // and upper cells.
  int iout = 0;
  for (int d=0; d<dim; d++) {
    if (idx[d] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[d] == num_cells[d]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}
