#pragma once

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

/**
* Raise or lower the indices of a rank-1 (spatial) tensor field using the value of the spatial metric tensor field at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param old_indices Array of current tensor indices (true for covariant, false for contravariant).
* @param new_indices Array of new tensor indices (true for covariant, false for contravariant).
* @param tensor Rank-1 (spatial) tensor field.
*/
GKYL_CU_D
double*
gkyl_gr_spatial_tensor_reindex_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double* tensor);

/**
* Raise or lower the indices of a rank-1 (spacetime) tensor field using the value of the spacetime metric tensor field at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param old_indices Array of current tensor indices (true for covariant, false for contravariant).
* @param new_indices Array of new tensor indices (true for covariant, false for contravariant).
* @param tensor Rank-1 (spacetime) tensor field.
*/
GKYL_CU_D
double*
gkyl_gr_spacetime_tensor_reindex_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* old_indices, bool* new_indices, double* tensor);

/**
* Compute the rank-2 covariant derivative of a rank-1 (spatial) tensor field using the values of the spatial Christoffel symbols at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param new_indices Array of tensor indices (true for covariant, false for contravariant).
* @param tensor Rank-1 (spatial) tensor field.
* @param tensor_der Rank-2 (spatial) partal derivative of the tensor field.
*/
GKYL_CU_D
double**
gkyl_gr_spatial_covariant_der_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, bool* indices, double* tensor, double** tensor_der);

/**
* Compute the rank-2 covariant derivative of a rank-1 (spacetime) tensor field using the values of the spacetime Christoffel symbols at a given point in spacetime.
*
* @param spacetime Base spacetime object.
* @param t Time coordinate.
* @param x Spatial coordinate (x-direction).
* @param y Spatial coordinate (y-direction).
* @param z Spatial coordinate (z-direction).
* @param dt Time coordinate spacing.
* @param dx Spatial coordinate spacing (x-direction).
* @param dy Spatial coordinate spacing (y-direction).
* @param dz Spatial coordinate spacing (z-direction).
* @param new_indices Array of tensor indices (true for covariant, false for contravariant).
* @param tensor Rank-1 (spatial) tensor field.
* @param tensor_der Rank-2 (spatial) partal derivative of the tensor field.
*/
GKYL_CU_D
double**
gkyl_gr_spacetime_covariant_der_rank1(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, bool* indices, double* tensor, double** tensor_der);