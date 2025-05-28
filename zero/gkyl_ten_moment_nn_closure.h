#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <kann.h>

struct gkyl_ten_moment_nn_closure_inp
{
  const struct gkyl_rect_grid *grid; // Grid on which to solve equations.
  int poly_order; // Polynomial order of learned DG coefficients.
  kann_t* ann; // Neural network architecture.
};

// Object type.
typedef struct gkyl_ten_moment_nn_closure gkyl_ten_moment_nn_closure;

/**
 * Create a new updater to update the pressure tensor in the ten moment equation system using a machine-learned magnetized closure trained on PKPM simulations,
 * i.e. where the heat flux tensor q_ijk and its derivatives are calculated from predictions for q_par and q_perp on the basis of rho, p_par and p_perp.
 * Returns RHS for accumulation in a forward-Euler method.
 *
 * @param inp Input parameters for updater.
 * @return Pointer to updater.
 */
gkyl_ten_moment_nn_closure*
gkyl_ten_moment_nn_closure_new(struct gkyl_ten_moment_nn_closure_inp inp);

/**
 * Delete updater.
 *
 * @param nnclosure Updater to delete.
 */
 void
 gkyl_ten_moment_nn_closure_release(gkyl_ten_moment_nn_closure *nnclosure);