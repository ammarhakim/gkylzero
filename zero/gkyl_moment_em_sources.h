#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_moment_em_sources_data {
    double charge;
    double mass;
};

struct gkyl_moment_em_sources_inp {
    const struct gkyl_rect_grid *grid; // grid on which to solve equations
    int nfluids; // number of fluids
    struct gkyl_moment_em_sources_data param[GKYL_MAX_SPECIES]; // species data
    double epsilon0; //
    int has_pressure; // 1 if pressure equation is present
};

// Object type
typedef struct gkyl_moment_em_sources gkyl_moment_em_sources;

/**
 * Create new updater to update electromagnetic sources in fluid
 * equations.  Uses implicit time-stepping (Time-centered
 * Crank-Nicholson).
 *
 * @param inp Input parameters to updater
 */
gkyl_moment_em_sources* gkyl_moment_em_sources_new(struct gkyl_moment_em_sources_inp inp);

/**
 * Compute implicit update of the electromagnetic source terms in the
 * multi-fluid system.  The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param mes Moment electromagnetic sources updater object
 * @param update_rng Range on which to run projection.
 * @param fluid Array of fluid variables (array size: nfluids)
 * @param em EM variables
 */

void gkyl_moment_em_sources_advance(const gkyl_moment_em_sources *mes, double dt,
  const struct gkyl_range *update_rng, struct gkyl_array *fluid[], struct gkyl_array *em);

/**
 * Delete updater.
 *
 * @param mes Updater to delete.
 */
void gkyl_moment_em_sources_release(gkyl_moment_em_sources *mes);
