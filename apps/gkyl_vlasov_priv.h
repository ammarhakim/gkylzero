// Private header for use in Vlasov app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_mom.h>

// Definitions of private structs and APIs attached to these objects
// for use in Vlasov app.

// ranges for use in BCs
struct vm_skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// data for moments
struct vm_species_moment {
  struct gkyl_mom_type *mtype;
  gkyl_mom_calc *mcalc;
  struct gkyl_array *marr;
  struct gkyl_array *marr_host;  
};

// species data
struct vm_species {
  struct gkyl_vlasov_species info; // data for species
    
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct vm_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *f_host; // host copy for use IO and initialization

  struct vm_species_moment m1i; // for computing currents
  struct vm_species_moment *moms; // diagnostic moments

  struct gkyl_dg_eqn *eqn; // Vlasov equation
  gkyl_hyper_dg *slvr; // solver 

  double* omegaCfl_ptr;
};

// field data
struct vm_field {
  struct gkyl_vlasov_field info; // data for field
    
  struct gkyl_array *em, *em1, *emnew; // arrays for updates
  struct gkyl_array *qmem; // array for q/m*(E,B)
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used for both copy and periodic)

  struct gkyl_array *em_host;  // host copy for use IO and initialization

  struct gkyl_dg_eqn *eqn; // Maxwell equation
  gkyl_hyper_dg *slvr; // solver

  double* omegaCfl_ptr;
};

// Vlasov object: used as opaque pointer in user code
struct gkyl_vlasov_app {
  char name[128]; // name of app
  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  double tcurr; // current time
  double cfl; // CFL number

  bool use_gpu; // should we use GPU (if present)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
    
  struct gkyl_rect_grid grid; // config-space grid
  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  struct vm_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  bool has_field; // has field (false if no field is specified)
  struct vm_field field; // field data

  // species data
  int num_species;
  struct vm_species *species; // species data

  struct gkyl_vlasov_stat stat; // statistics
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Compute out = c1*arr1 + c2*arr2
static inline struct gkyl_array*
array_combine(struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range rng)
{
  return gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng),
    c2, arr2, rng);
}


// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct vm_skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

/** vm_species_moment API */

void vm_species_moment_init(struct gkyl_vlasov_app *app, struct vm_species *s,
  struct vm_species_moment *sm, const char *nm);

void vm_species_moment_release(const struct gkyl_vlasov_app *app, const struct vm_species_moment *sm);


/** vm_species API */

void vm_species_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_species *s);

double vm_species_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  const struct gkyl_array *fin, const struct gkyl_array *qmem, struct gkyl_array *rhs);

void vm_species_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f);

void vm_species_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_species *species,
  int dir, struct gkyl_array *f);

void vm_species_apply_bc(gkyl_vlasov_app *app, const struct vm_species *species, struct gkyl_array *f);

void vm_species_release(const gkyl_vlasov_app* app, const struct vm_species *s);

/** vm_field API */

void vm_field_init(struct gkyl_vm *vm, struct gkyl_vlasov_app *app, struct vm_field *f);

double vm_field_rhs(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em, struct gkyl_array *rhs);

void vm_field_apply_periodic_bc(gkyl_vlasov_app *app, const struct vm_field *field, int dir, struct gkyl_array *f);

void vm_field_apply_copy_bc(gkyl_vlasov_app *app, const struct vm_field *field, int dir, struct gkyl_array *f);

void vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f);

void vm_field_release(const gkyl_vlasov_app* app, const struct vm_field *f);
