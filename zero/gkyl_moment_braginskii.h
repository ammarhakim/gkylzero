#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Identifiers for different Braginskii types
enum gkyl_braginskii_type {
  GKYL_UNMAG_BRAG, // Unmagnetized Braginskii
  GKYL_MAG_BRAG, // Magnetized Braginskii
};

struct gkyl_moment_braginskii_data {
    enum gkyl_eqn_type type_eqn; // equation type
    double charge; // species charge
    double mass; // species mass
};

struct gkyl_moment_braginskii_inp {
    const struct gkyl_rect_grid *grid; // grid on which to solve equations
    int nfluids; // number of fluids
    struct gkyl_moment_braginskii_data param[GKYL_MAX_SPECIES]; // species data
    enum gkyl_braginskii_type type_brag; // which Braginskii equations
    double epsilon0;
};

// Object type
typedef struct gkyl_moment_braginskii gkyl_moment_braginskii;

/**
 * Create new updater to update fluid equations with Braginskii 
 * transport terms (viscosity, heat conduction, etc.).
 * Returns RHS for accumulation in a forward Euler method.
 *
 * @param inp Input parameters to updater
 */
gkyl_moment_braginskii* gkyl_moment_braginskii_new(struct gkyl_moment_braginskii_inp inp);

/**
 * Compute RHS contribution from Braginskii transport terms in the
 * multi-fluid system.  The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param bes Braginskii updater object.
 * @param update_rng Range on which to compute.
 * @param fluid Input array of fluid variables (array size: nfluids)
 * @param em Total EM variables (plasma + external)
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output (NOTE: Returns RHS output of all nfluids)
 */

void gkyl_moment_braginskii_advance(const gkyl_moment_braginskii *bes, struct gkyl_range update_rng,
  const struct gkyl_array *fluid[], const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, struct gkyl_array *rhs[]);

/**
 * Delete updater.
 *
 * @param bes Updater to delete.
 */
void gkyl_moment_braginskii_release(gkyl_moment_braginskii *bes);
