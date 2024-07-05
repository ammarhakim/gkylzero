#pragma once

#include <gkyl_block_geom.h>
#include <gkyl_moment.h>

typedef struct gkyl_moment_multib_app gkyl_moment_multib_app;

// data common to species on all blocks
struct gkyl_moment_multib_species_info {
  char name[128]; // species name
  double charge, mass; // charge and mass
 
  struct gkyl_wv_eqn *equation; // equation object
  enum gkyl_wave_limiter limiter; // limiter to use
  enum gkyl_wave_split_type split_type; // edge splitting to use

  int evolve; // evolve species? 1-yes, 0-no
  bool force_low_order_flux; // should  we force low-order flux?
};

// data for individual blocks
struct gkyl_moment_multib_species {
  int block_id; // block ID
  
  void *ctx; // context for initial condition init function (and potentially other functions)
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  bool is_app_accel_static; // flag to indicate if applied acceleration is static
  void *app_accel_ctx; // context for applied acceleration
  // pointer to applied acceleration/forces function
  void (*app_accel_func)(double t, const double *xn, double *fout, void *ctx);

  void *nT_source_ctx; // context for nT source
  // pointer to user-defined number density and temperature sources
  void (*nT_source_func)(double t, const double *xn, double *fout, void *ctx);
  bool nT_source_set_only_once;

  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];

  // for function BCs these should be set
  wv_bc_func_t bcx_func[2], bcy_func[2], bcz_func[2];  
};

// data common to species on all blocks
struct gkyl_moment_multib_field_info {
  double epsilon0, mu0;
  double elc_error_speed_fact, mag_error_speed_fact;

  enum gkyl_wave_limiter limiter; // limiter to use

  int evolve; // evolve field? 1-yes, 0-no
  bool use_explicit_em_coupling; // flag to indicate if using explicit em-coupling
};

// data for individual blocks
struct gkyl_moment_multib_field {
  int block_id; // block ID
  
  void *ctx; // context for initial condition init function (and potentially other functions)
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  void *app_current_ctx; // context for applied current
  // pointer to applied current function
  void (*app_current_func)(double t, const double *xn, double *fout, void *ctx);
  double t_ramp_curr; // linear ramp for turning on applied currents

  bool is_ext_em_static; // flag to indicate if external field is time-independent
  void *ext_em_ctx; // context for applied current
  // pointer to external fields
  void (*ext_em_func)(double t, const double *xn, double *fout, void *ctx);
  double t_ramp_ext_em; // linear ramp for turning on external E field

  // boundary conditions
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];
  // for function BCs these should be set
  wv_bc_func_t bcx_func[2], bcy_func[2], bcz_func[2];
};

// Top-level app parameters: this
struct gkyl_moment_multib {
  char name[128]; // name of app

 // geometry and topology of all blocks in simulation  
  struct gkyl_block_geom *block_geom;

  int num_species; // number of species
  // block-independent species info
  struct gkyl_moment_multib_species_info species_info[GKYL_MAX_SPECIES];
  // per-block species info
  struct gkyl_moment_multib_species *species[GKYL_MAX_SPECIES];
 
  // block-independent field info
  struct gkyl_moment_multib_field_info field_info;
  // per-block species info
  struct gkyl_moment_multib_field *field; 
};


/**
 * Construct a new moments multi-block app.
 *
 * @param mbinp Multi-block App inputs. See struct docs.
 * @return New multi-block moment app object.
 */
struct gkyl_moment_multib_app* gkyl_moment_multib_app_new(struct gkyl_moment_multib *mbinp);
                                
  
