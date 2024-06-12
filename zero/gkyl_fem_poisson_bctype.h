#pragma once

// Enumeration of boundary condition types.
enum gkyl_poisson_bc_type {
  GKYL_POISSON_PERIODIC = 0, // Periodic boundary conditions (copied across domain).
  GKYL_POISSON_DIRICHLET, // Dirichlet boundary conditions (impose value at boundary).
  GKYL_POISSON_NEUMANN, // Neumann boundary conditions (impose derivative normal to the boundary).
  GKYL_POISSON_ROBIN, // Robin boundary condition (impose linear combination of value and normal derivative at boundary).
};

// Container for boundary condition values (either a value, a normal derivative, or a linear combination of the two).
// For GKYL_POISSON_PERIODIC, no boundary values are required.
// For GKYL_POISSON_DIRICHLET or GKYL_POISSON_NEUMANN, one value is required (either a boundary value or a derivative normal to the boundary).
// For GKYL_POISSON_ROBIN, three values are required (linear combination of a boundary value and a normal derivative value).
struct gkyl_poisson_bc_value {
  double v[3];
};

// Container for boundary condition types and values in each coordinate direction.
struct gkyl_poisson_bc {
  enum gkyl_poisson_bc_type lo_type[GKYL_MAX_CDIM]; // Type of lower boundary condition.
  struct gkyl_poisson_bc_value lo_value[GKYL_MAX_CDIM]; // Value(s) for lower boundary condition.

  enum gkyl_poisson_bc_type up_type[GKYL_MAX_CDIM]; // Type of upper boundary condition.
  struct gkyl_poisson_bc_value up_value[GKYL_MAX_CDIM]; // Value(s) for upper boundary condition
};

/**
* Compute the kernel index (in the array of kernels) needed for a given grid index, assuming a kernel that distinguishes between upper cells 
* and interior cells.
* @param dim Number of dimensions.
* @param num_cells Array of cell counts.
* @param idx Grid index.
* @return Kernel index (in the array of kernels).
*/
GKYL_CU_DH
int
idx_to_inup_ker(const int dim, const int* num_cells, const int* idx);

/**
* Compute the kernel index (in the array of kernels) needed for a given grid index, assuming a kernel that distinguishes between lower cells,
* upper cells and interior cells.
* @param dim Number of dimensions.
* @param num_cells Array of cell counts.
* @param idx Grid index.
* @return Kernel index (in the array of kernels).
*/
GKYL_CU_DH
int
idx_to_inloup_ker(const int dim, const int* num_cells, const int* idx);