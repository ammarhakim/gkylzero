// Private header for use in moment app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fv_proj.h>
#include <gkyl_moment.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wv_ten_moment.h>

// Species data
struct moment_species {
  int ndim;
  char name[128]; // species name
  double charge, mass;
  double k0; // closure parameter (default is 0.0, used by 10 moment)

  int evolve; // evolve species? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
    
  struct gkyl_array *fdup, *f[4]; // arrays for updates
  struct gkyl_array *app_accel; // array for applied acceleration/forces
  // pointer to projection operator for applied acceleration/forces function
  gkyl_fv_proj *proj_app_accel;
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  enum gkyl_eqn_type eqn_type; // type ID of equation
  int num_equations; // number of equations in species
  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary condition type
  enum gkyl_species_bc_type lower_bct[3], upper_bct[3];
  // boundary condition solvers on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_q; // integrated conserved quantities
  bool is_first_q_write_call; // flag for dynvec written first time  
};

// Field data
struct moment_field {
  int ndim;
  double epsilon0, mu0;

  int evolve; // evolve species? 1-yes, 0-no
    
  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);    
    
  struct gkyl_array *fdup, *f[4]; // arrays for updates
  struct gkyl_array *app_current; // arrays for applied currents
  // pointer to projection operator for applied current function
  gkyl_fv_proj *proj_app_current;


  bool is_ext_em_static; // flag to indicate if external field is time-independent
  struct gkyl_array *ext_em; // array external fields  
  gkyl_fv_proj *proj_ext_em;   // pointer to projection operator for external fields
  bool was_ext_em_computed; // flag to indicate if we already computed external EM field
  
  struct gkyl_array *bc_buffer; // buffer for periodic BCs

  gkyl_wave_prop *slvr[3]; // solver in each direction

  // boundary condition type
  enum gkyl_field_bc_type lower_bct[3], upper_bct[3];
  // boundary conditions on lower/upper edges in each direction
  gkyl_wv_apply_bc *lower_bc[3], *upper_bc[3];

  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for dynvec written first time
};

// Source data
struct moment_coupling {
  gkyl_moment_em_coupling *slvr; // source solver function
};

// Moment app object: used as opaque pointer in user code
struct gkyl_moment_app {
  char name[128]; // name of app
  int ndim; // space dimensions
  double tcurr; // current time
  double cfl; // CFL number

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int is_dir_skipped[3]; // flags to tell if update in direction are skipped

  enum gkyl_moment_fluid_scheme fluid_scheme; // scheme to update fluid equations
    
  struct gkyl_rect_grid grid; // grid
  struct gkyl_range local, local_ext; // local, local-ext ranges

  bool has_mapc2p; // flag to indicate if we have mapc2p
  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_wave_geom *geom; // geometry needed for species and field solvers

  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  int has_field; // flag to indicate if we have a field
  struct moment_field field; // field data
    
  // species data
  int num_species;
  struct moment_species *species; // species data

  int update_sources; // flag to indicate if sources are to be updated
  struct moment_coupling sources; // sources
    
  struct gkyl_moment_stat stat; // statistics
};

// Function pointer to compute integrated quantities from input
typedef void (*integ_func)(int nc, const double *qin, double *integ_out);
