// Private header for use in Gyrokinetic app: do not include in user-facing
// header files!
#pragma once

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <stc/cstr.h>

#include <gkyl_alloc.h>
#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_basic.h>
#include <gkyl_bc_emission.h>
#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_boundary_flux.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_gk_neut_hamil.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gyrokinetic_vars.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_recomb.h>
#include <gkyl_dg_updater_diffusion_gyrokinetic.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_dg_updater_lbo_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_vlasov.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_deflated_fem_poisson.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_gk_maxwellian_correct.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_gk_maxwellian_moments.h>
#include <gkyl_tok_geo.h>
#include <gkyl_velocity_map.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_cross_prim_moms_bgk.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_null_pool.h>
#include <gkyl_position_map.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_powsqrt_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_radiation_read.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_tok_geo.h>
#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>

// Definitions of private structs and APIs attached to these objects
// for use in Gyrokinetic app.

// Meta-data for IO
struct gyrokinetic_output_meta {
  int frame; // frame number
  double stime; // output time
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
  char basis_type_nm[64]; // used during read
};

// list of valid moment names for gyrokinetics
static const char *const valid_moment_names[] = {
  "M0",
  "M1",
  "M2",
  "M2par",
  "M2perp",
  "M3par",
  "M3perp",
  "ThreeMoments",
  "FourMoments",
  "MaxwellianMoments", // internal flag for whether we are computing (n, u_par, T/m)
  "BiMaxwellianMoments", // internal flag for whether we are computing (n, u_par, T_par/m, T_perp/m)
  "HamiltonianMoments", // Compute the H moment of f.
  "Integrated", // this is an internal flag, not for passing to moment type
};

// check if name of moment is valid or not for gyrokinetics
static bool
is_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_moment_names)/sizeof(valid_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// list of valid moment names for neutrals
static const char *const valid_neut_moment_names[] = {
  "M0",
  "M1i",
  "M2ij",
  "M2",
  "M3i",
  "M3ijk",
  "FiveMoments",
  "LTEMoments", // this is an internal flag for computing moments (n, V_drift, T/m)
                // of the LTE (local thermodynamic equilibrium) distribution
  "Integrated", // this is an internal flag, not for passing to moment type
};

// check if name of moment is valid or not for neutrals
static bool
is_neut_moment_name_valid(const char *nm)
{
  int n = sizeof(valid_neut_moment_names)/sizeof(valid_neut_moment_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(valid_neut_moment_names[i], nm) == 0)
      return 1;
  return 0;
}

// struct for holding moment correction inputs
struct correct_all_moms_inp {
  bool correct_all_moms; // boolean if we are correcting all the moments or only density
  double iter_eps; // error tolerance for moment fixes (density is always exact)
  int max_iter; // maximum number of iterations
  bool use_last_converged; // use last iteration value regardless of convergence?
};

// data for gyrokinetic moments
struct gk_species_moment {
  struct gk_geometry *gk_geom; // geometry struct for dividing moments by Jacobian
  struct gkyl_dg_bin_op_mem *mem_geo; // memory needed in dividing moments by Jacobian
  bool is_integrated; // boolean for if computing integrated moments 
                      // integrated moments do not need to divide by Jacobian since
                      // the inverse Jacobian is already included in the computation
  int num_mom; // number of moments 

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)

  // Options for moment calculation: 
  // 1. Compute the moment directly with dg_updater_moment_gyrokinetic
  // 2. Compute the moments of the equivalent Maxwellian (n, u_par, T/m)
  // 3. Compute the moments of the equivalent Bi-Maxwellian (n, u_par, T_par/m, T_perp/m)
  //    Latter two options use specialized gkyl_gyrokinetic_maxwellian_moments updater
  union {
    struct {
      struct gkyl_gk_maxwellian_moments *gyrokinetic_maxwellian_moms; 
    };
    struct {
      struct gkyl_vlasov_lte_moments *vlasov_lte_moms; // updater for computing LTE moments
    };
    struct {
      struct gkyl_dg_updater_moment *mcalc; 
    };
  };
  bool is_maxwellian_moms;
  bool is_bimaxwellian_moms;
};

struct gk_rad_drag {  
  enum gkyl_radiation_id radiation_id; // type of radiation
  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  struct gk_neut_species *collide_with_neut[GKYL_MAX_SPECIES]; // pointers to neutral cross-species we collide with
  int collide_with_idx[2*GKYL_MAX_SPECIES]; // index of species we collide with
  bool is_neut_species[2*GKYL_MAX_SPECIES]; // Flag of whether neutral or gk species
  
  // drag coefficients in vparallel and mu for each species being collided with
  struct gkyl_gk_rad_drag *vnu_surf;
  struct gkyl_gk_rad_drag *vnu;
  struct gkyl_gk_rad_drag *vsqnu_surf;
  struct gkyl_gk_rad_drag *vsqnu;
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars;
  struct gkyl_array *rad_fit_ne[2*GKYL_MAX_SPECIES];

  struct gk_species_moment moms[2*GKYL_MAX_SPECIES]; // moments needed in radiation update (need number density)

  struct gk_species_moment m2; // m2 of radiation update (needed for emissivity)
  struct gkyl_array *emissivity[2*GKYL_MAX_SPECIES];
  struct gkyl_array *emissivity_host[2*GKYL_MAX_SPECIES];
  struct gkyl_array *emissivity_rhs;
  struct gkyl_array *emissivity_denominator;

  struct gkyl_array *nvnu_surf; // total vparallel radiation drag surface expansion including density scaling
  struct gkyl_array *nvnu; // total vparallel radiation drag volume expansion including density scaling
  struct gkyl_array *nvsqnu_surf; // total mu radiation drag surface expansion including density scaling
  struct gkyl_array *nvsqnu; // total mu radiation drag volume expansion including density scaling

  struct gkyl_array **vtsq_min_per_species;  // Smallest vtsq that radiation is calculated (one for each fit), divided by configuration space normalization
  struct gk_species_moment prim_moms;
  struct gkyl_array *vtsq;
  struct gkyl_array *m0;

  // host-side copies for I/O
  struct gkyl_array *nvnu_surf_host; 
  struct gkyl_array *nvnu_host; 
  struct gkyl_array *nvsqnu_surf_host; 
  struct gkyl_array *nvsqnu_host; 

  gkyl_dg_updater_collisions *drag_slvr; // radiation solver

  struct gk_species_moment integ_moms; // integrated moments
  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time
};

// forward declare species struct
struct gk_species;

struct gk_lbo_collisions {  
  enum gkyl_collision_id collision_id; // type of collisions
  bool write_diagnostics; // Whether to write diagnostics out.
  struct gkyl_array *boundary_corrections; // LBO boundary corrections
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator
  struct gkyl_array *nu_sum, *prim_moms, *nu_prim_moms; // LBO primitive moments
  struct gkyl_array *nu_sum_host, *prim_moms_host, *nu_prim_moms_host; // LBO primitive moments host-side for I/O
  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  double self_nu_fac; // Self collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double cross_nu_fac[GKYL_MAX_SPECIES]; // Cross collision freqs without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double vtsq_min; // minimum vtsq
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  double betaGreenep1; // value of Greene's factor beta + 1
  double other_m[GKYL_MAX_SPECIES]; // masses of species being collided with
  struct gkyl_array *other_prim_moms[GKYL_MAX_SPECIES]; // self-primitive moments of species being collided with
  struct gkyl_array *cross_prim_moms[GKYL_MAX_SPECIES]; // LBO cross-primitive moments
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // LBO cross-species collision frequencies
  struct gkyl_array *other_nu[GKYL_MAX_SPECIES];
  struct gkyl_array *cross_nu_prim_moms; // weak multiplication of collision frequency and primitive moments
  
  struct gkyl_array *self_nu, *self_nu_prim_moms; // LBO self-primitive moments

  struct gk_species_moment moms; // Moments needed in LBO (M0, M1, M2).
  gkyl_dg_bin_op_mem *dg_div_mem; // Memory needed for weak division.

  struct gkyl_array *m0;
  struct gkyl_array *vtsq;
  struct gkyl_array *m2self; // m2self used for robustness of LBO
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

struct gk_lte {  
  struct gkyl_array *f_lte;

  struct gk_species_moment moms; // moments needed in the equilibrium

  // LTE distribution function projection object
  // also corrects the density of projected distribution function
  union {
    struct {
      struct gkyl_gk_maxwellian_proj_on_basis *proj_max; 
      struct gkyl_gk_maxwellian_correct *corr_max; 
    };
    struct {
      struct gkyl_vlasov_lte_proj_on_basis *proj_lte; 
      struct gkyl_vlasov_lte_correct *corr_lte; 
    };
  };

  long n_iter; // total number of iterations from correcting moments
  long num_corr; // total number of times the correction updater is called
  bool correct_all_moms; // boolean if we are correcting all the moments
  gkyl_dynvec corr_stat;
  bool is_first_corr_status_write_call;
};

struct gk_bgk_collisions {  
  enum gkyl_collision_id collision_id; // type of collisions
  bool write_diagnostics; // Whether to write diagnostics out.
  struct gkyl_array *nu_sum; // BGK collision frequency 
  struct gkyl_array *nu_sum_host; // BGK collision frequency host-side for I/O
  struct gkyl_array *self_nu; // BGK self-collision frequency

  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  double self_nu_fac; // Self collision frequency without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double cross_nu_fac[GKYL_MAX_SPECIES]; // Cross collision freqs without factor of n_r/(v_ts^2+v_tr^2)^(3/2)
  double vtsq_min; // minimum vtsq
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  struct gk_species_moment moms; // Moments needed in BGK (M0, M1, M2).
  struct gkyl_array *m0;
  struct gkyl_array *vtsq;
  
  struct gkyl_array *nu_fmax;

  // Cross collisions inputs, arrays, and updaters
  double betaGreenep1; // value of Greene's factor beta + 1
  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  double other_m[GKYL_MAX_SPECIES]; // masses of species being collided with
  struct gkyl_array *other_moms[GKYL_MAX_SPECIES]; // moments of species being collided with
  struct gkyl_array *other_nu[GKYL_MAX_SPECIES]; // cross-species collision frequencies
  struct gkyl_array *cross_nu[GKYL_MAX_SPECIES]; // cross-species collision frequencies

  struct gkyl_array *cross_moms[GKYL_MAX_SPECIES];
  struct gkyl_array *cross_moms_host[GKYL_MAX_SPECIES];
  struct gkyl_gyrokinetic_cross_prim_moms_bgk *cross_bgk; // cross-species moment computation

  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep)

  bool implicit_step; // whether or not to take an implcit bgk step
  double dt_implicit; // timestep used by the implicit collisions  
};

struct gk_boundary_fluxes {
  struct gk_species_moment gammai[2*GKYL_MAX_CDIM]; // integrated moments
  gkyl_boundary_flux *flux_slvr[2*GKYL_MAX_CDIM]; // boundary flux solver

  // Objects used for boundary flux diagnostics.
  bool allocated_diagnostic; // Signal diagnostic objects were allocated.
  bool allocated_solver; // Signal solver objects were allocated.
  int num_boundaries; // Number of boundaries to compute bfluxes at.
  int boundaries_dir[2*GKYL_MAX_CDIM]; // Direction of bflux boundaries.
  enum gkyl_edge_loc boundaries_edge[2*GKYL_MAX_CDIM]; // Edge of bflux boundaries.
  struct gkyl_range *boundaries_conf_ghost[2*GKYL_MAX_CDIM]; // Ghost range of bflux boundaries.
  bool is_hamiltonian_intmom; // True if need Hamiltonian integrated moments.
  struct gkyl_bc_basic *gfss_bc_op[2*GKYL_MAX_CDIM]; // Applies BCs to bmag and phi.
  struct gkyl_array *bc_buffer; // Buffer used by gfss_bc_op;
  struct gkyl_array **f, **f1, **fnew; // Boundary flux through each boundary (one for each RK stage).
  struct gk_species_moment moms_op; // Moments calculator.
  struct gkyl_array_integrate *integ_op; // Operator that integrates over volume.
  double *int_moms_local, *int_moms_global; // Integrated moments in this time step.
  gkyl_dynvec intmom[2*GKYL_MAX_CDIM]; // Integrated moments of the boundary fluxes.
  double *intmom_cumm_buff; // Cummulative (in time) integrated moments of the boundary fluxes.
  bool is_first_intmom_write_call[2*GKYL_MAX_CDIM]; // Flag 1st writing of blux_intmom.
  bool is_not_first_restart_write_call; // False at first write after restart.
  // Function pointers to various methods.
  void (*bflux_clear_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fin, double val);
  void (*bflux_scale_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fin, double val);
  void (*bflux_step_f_func)(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux, struct gkyl_array **fout,
    double dt, const struct gkyl_array **fin);
  void (*bflux_combine_func)(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
    struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2);
  void (*bflux_copy_func)(gkyl_gyrokinetic_app *app, struct gk_species *gk_s, struct gk_boundary_fluxes *bflux,
    struct gkyl_array **fout, struct gkyl_array **fin);
  void (*bflux_calc_voltime_int_mom_func)(gkyl_gyrokinetic_app* app, const struct gk_species *gk_s,
    struct gk_boundary_fluxes *bflux, double tm);
  void (*bflux_rhs_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gk_boundary_fluxes *bflux,
    const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms);
  void (*bflux_calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, const struct gk_species *gk_s,
    struct gk_boundary_fluxes *bflux, double tm);
  void (*bflux_write_integrated_mom_func)(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
    struct gk_boundary_fluxes *bflux);

  // Additional objects needed for recycle BCs.
  struct gkyl_rect_grid boundary_grid[2*GKYL_MAX_CDIM];
  struct gkyl_rect_grid conf_boundary_grid[2*GKYL_MAX_CDIM];
  struct gkyl_array *flux_arr[2*GKYL_MAX_CDIM];
  struct gkyl_array *mom_arr[2*GKYL_MAX_CDIM];
  struct gkyl_range flux_r[2*GKYL_MAX_CDIM];
  struct gkyl_range conf_r[2*GKYL_MAX_CDIM];
  struct gkyl_dg_updater_moment *integ_moms[2*GKYL_MAX_CDIM];
  gkyl_ghost_surf_calc *ghost_flux_slvr; // boundary flux solver for gk_neut_species
  
};

struct gk_recycle_wall {
  // recycling wall boundary conditions
  int num_species;
  int dir;
  enum gkyl_edge_loc edge;
  double *scale_ptr;
  double t_bound;
  bool elastic;
  struct gkyl_array *f0;
  double delta[GKYL_MAX_SPECIES]; // recycling fraction
  
  gkyl_boundary_flux *f0_flux_slvr[2*GKYL_MAX_CDIM];
  struct gkyl_dg_bin_op_mem *mem_geo; // memory needed in dividing moments by Jacobian
  struct gkyl_spectrum_model *spectrum_model[GKYL_MAX_SPECIES];
  struct gkyl_yield_model *yield_model[GKYL_MAX_SPECIES];
  struct gkyl_elastic_model *elastic_model;
  struct gkyl_bc_emission_ctx *params;

  struct gkyl_bc_emission_spectrum *update[GKYL_MAX_SPECIES];
  struct gkyl_bc_emission_elastic *elastic_update;
  struct gkyl_rect_grid *init_conf_grid;
  struct gkyl_array *f_emit;
  struct gkyl_array *f_diag;
  struct gkyl_array *init_flux;
  struct gkyl_array *init_bflux_arr;
  struct gkyl_array *emit_flux;
  struct gkyl_array *diag_out;
  struct gkyl_array *diag_out_ho;
  struct gkyl_array *buffer;
  struct gkyl_array *elastic_yield;
  struct gkyl_array *yield[GKYL_MAX_SPECIES]; // projected secondary electron yield
  struct gkyl_array *spectrum[GKYL_MAX_SPECIES]; // projected secondary electron spectrum
  struct gkyl_array *weight[GKYL_MAX_SPECIES];
  struct gkyl_array *flux[GKYL_MAX_SPECIES];
  struct gkyl_array *bflux_arr[GKYL_MAX_SPECIES];
  struct gkyl_array *k[GKYL_MAX_SPECIES];
  struct gk_species *impact_species[GKYL_MAX_SPECIES]; // pointers to impacting species
  struct gkyl_range impact_normal_r[GKYL_MAX_SPECIES];
  struct gkyl_dg_updater_moment *flux_slvr[GKYL_MAX_SPECIES]; // moments
  struct gkyl_dg_updater_moment *init_flux_slvr; 

  struct gkyl_rect_grid *impact_conf_grid[GKYL_MAX_SPECIES];
  struct gkyl_rect_grid *impact_grid[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_ghost_r[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_skin_r[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_buff_r[GKYL_MAX_SPECIES];
  struct gkyl_range *impact_cbuff_r[GKYL_MAX_SPECIES];
  
  struct gkyl_rect_grid *emit_grid;
  struct gkyl_range *emit_buff_r;
  struct gkyl_range *emit_cbuff_r;
  struct gkyl_range *emit_ghost_r;
  struct gkyl_range *emit_skin_r;
  struct gkyl_range emit_normal_r;
};

struct gk_react {
  int num_react; // number of reactions
  bool all_gk; // boolean for if reactions are only between gyrokinetic species
  bool write_diagnostics; // Whether to write diagnostics out.
  struct gkyl_gyrokinetic_react_type react_type[GKYL_MAX_REACT]; // input struct for type of reactions

  struct gkyl_array *f_react; // distribution function array which holds update for each reaction
                              // form depends on type_self, e.g., for ionization and type_self == GKYL_SELF_ELC
                              // f_react = n_donor*(fmax1(n_elc, upar_elc, vtiz1^2) + fmax2(n_elc, upar_donor, vtiz2^2) - f_elc)
                              // RHS update is then obtained by incrementing rhs += coeff_react*f_react

  enum gkyl_react_id react_id[GKYL_MAX_REACT]; // what type of reaction (ionization, charge exchange, recombination)
  enum gkyl_react_self_type type_self[GKYL_MAX_REACT]; // what is the role of species in this reaction
  struct gk_species *species_elc[GKYL_MAX_REACT]; // pointers to electron species being reacted with
  struct gk_species *species_ion[GKYL_MAX_REACT]; // pointers to ion species being reacted with
  int elc_idx[GKYL_MAX_REACT]; // integer index of electron species being reacted with 
  int ion_idx[GKYL_MAX_REACT]; // integer index of ion species being reacted with 
  int donor_idx[GKYL_MAX_REACT]; // integer index of donor species being reacted with 
  int partner_idx[GKYL_MAX_REACT]; // integer index of neut species in cx reaction

  struct gkyl_array *coeff_react[GKYL_MAX_REACT]; // reaction rate
  struct gkyl_array *coeff_react_host[GKYL_MAX_REACT]; // reaction rate
  struct gkyl_array *vt_sq_iz1[GKYL_MAX_REACT]; // ionization temperature
  struct gkyl_array *vt_sq_iz2[GKYL_MAX_REACT]; // ionization temperature
  struct gkyl_array *Jm0_elc[GKYL_MAX_REACT]; // J*electron density
  struct gkyl_array *Jm0_ion[GKYL_MAX_REACT]; // J*ion density
  struct gkyl_array *Jm0_donor[GKYL_MAX_REACT]; // J*donor density
  struct gkyl_array *Jm0_partner[GKYL_MAX_REACT]; // J*partner density

  struct gkyl_array *react_lte_moms[GKYL_MAX_REACT]; // LTE/Maxwellian moments used for Maxwellian projection
  struct gkyl_array *u_i[GKYL_MAX_REACT]; // Vector flow velocity for neutrals (ux, uy, uz)
  struct gkyl_array *u_i_dot_b_i[GKYL_MAX_REACT]; // u_i . b_i (Cartesian compoments of magnetic field unit vector)
  struct gkyl_array *vt_sq_donor[GKYL_MAX_REACT]; // Donor thermal velocity 
  struct gkyl_array *upar_ion[GKYL_MAX_REACT]; // Ion vector parallel flow velocity upar b_i
  struct gkyl_array *vt_sq_ion[GKYL_MAX_REACT]; // Ion thermal velocity 
  struct gkyl_array *vt_sq_partner[GKYL_MAX_REACT]; // Neutral (partner) thermal velocity
  
  union {
    // ionization
    struct {
      struct gkyl_dg_iz *iz[GKYL_MAX_REACT];
    };
    // recombination
    struct {
      struct gkyl_dg_recomb *recomb[GKYL_MAX_REACT];
    };
    // charge exchange
    struct {
      struct gkyl_dg_cx *cx[GKYL_MAX_REACT];
    };
  };
};

// Context for c2p function passed to proj_on_basis.
struct gk_proj_on_basis_c2p_func_ctx {
  int cdim, vdim;
  struct gkyl_position_map *pos_map;
  struct gkyl_velocity_map *vel_map;
};

struct gk_proj {
  enum gkyl_projection_id proj_id; // type of projection
  struct gk_proj_on_basis_c2p_func_ctx proj_on_basis_c2p_ctx; // c2p function context.
  // organization of the different projection objects and the required data and solvers
  union {
    // function projection
    struct {
      struct gkyl_proj_on_basis *proj_func; // projection operator for specified function
      struct gkyl_array *proj_host; // array for projection on host-side if running on GPUs
    };
    // Maxwellian and Bi-Maxwellian projection from primitive moments
    struct {
      struct gkyl_array *dens; // host-side density

      struct gkyl_array *prim_moms_host; // host-side prim_moms for initialization with proj_on_basis
      struct gkyl_array *prim_moms; // prim_moms we pass to Maxwellian projection object (potentially on device)

      bool correct_all_moms; // boolean if we are correcting all the moments

      struct gkyl_proj_on_basis *proj_dens; // projection operator for density
      struct gkyl_array *vtsq; // host-side vth^2 = T/m (temperature/mass)
      struct gkyl_proj_on_basis *proj_temp; // projection operator for temperature
      struct gkyl_array *vtsqpar; // host-side vth_par^2 = Tpar/m (parallel temperature/mass)
      struct gkyl_array *vtsqperp; // host-side vth_perp^2 = Tperp/m (perpendicular temperature/mass)
      struct gkyl_proj_on_basis *proj_temppar; // projection operator for parallel temperature
      struct gkyl_proj_on_basis *proj_tempperp; // projection operator for parallel temperature
      union {
        struct { 
          struct gkyl_array *upar; // host-side upar for GK Maxwellian/Bi-Maxwellian projection
          struct gkyl_proj_on_basis *proj_upar; // projection operator for upar for GK Maxwellian/Bi-Maxwellian projection
          struct gkyl_gk_maxwellian_proj_on_basis *proj_max; // Maxwellian projection object for GK
          struct gkyl_gk_maxwellian_correct *corr_max; // Maxwellian correction object for GK
        };
        struct { 
          struct gkyl_array *udrift; // host-side udrift for Vlasov neutrals LTE projection
          struct gkyl_proj_on_basis *proj_udrift; // projection operator for udrift for Vlasov neutrals LTE projection
          struct gkyl_vlasov_lte_proj_on_basis *proj_lte; // Maxwellian projection object for Vlasov neutrals
          struct gkyl_vlasov_lte_correct *corr_lte; // Maxwellian correction object for Vlasov neutrals
        };     
      };  
    };
  };
};

struct gk_source {
  enum gkyl_source_id source_id; // type of source
  bool evolve; // Whether the source is time dependent.
  struct gkyl_array *source; // applied source
  struct gkyl_array *source_host; // host copy for use in IO and projecting
  struct gk_proj proj_source[GKYL_MAX_SOURCES]; // projector for source
  int num_sources; // Number of sources.

  int num_diag_moments; // number of diagnostics moments
  struct gk_species_moment *moms; // diagnostic moments
  struct gk_species_moment integ_moms; // integrated moments
  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time
};

// species data
struct gk_species {
  struct gkyl_gyrokinetic_species info; // data for species

  enum gkyl_gkmodel_id gkmodel_id;
  enum gkyl_gkfield_id gkfield_id;

  struct gkyl_basis basis; // phase-space basis

  // pointer to basis on device
  // (points to host structs if not on GPU)
  struct gkyl_basis *basis_on_dev; 
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_velocity_map *vel_map; // Velocity mapping objects.

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // fixed buffers for time independent BCs 

  struct gkyl_array *f_host; // host copy for use IO and initialization

  struct gkyl_array *alpha_surf; // array for surface phase space flux
  struct gkyl_array *sgn_alpha_surf; // array for the sign of the surface phase space flux at quadrature points
                                     // utilized for numerical flux function
                                     // F = alpha_surf/2 ( (f^+ + f^-) - sign_alpha_surf*(f^+ - f^-) )
  struct gkyl_array *const_sgn_alpha; // boolean array for if the surface phase space flux is single signed
                                      // if true, numerical flux function inside kernels simplifies to
                                      // F = alpha_surf*f^- (if sign_alpha_surf = 1), 
                                      // F = alpha_surf*f^+ (if sign_alpha_surf = -1)
  
  struct gkyl_array *gyro_phi; // Gyroaveraged electrostatic potential.
  // organization of the different equation objects and the required data and solvers
  union {
    // EM GK model
    struct {
      struct gkyl_array *apar; // array for A_parallel
      struct gkyl_array *apardot; // array for d/dt A_parallel
    };
  };

  struct gkyl_dg_calc_gyrokinetic_vars *calc_gk_vars;

  struct gk_species_moment m0; // for computing charge density
  struct gk_species_moment integ_moms; // integrated moments
  struct gk_species_moment *moms; // diagnostic moments
  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // Integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  struct gkyl_array *fdot_mom_old, *fdot_mom_new; // Moments of f_old and f_new.
  gkyl_dynvec fdot_integ_diag; // Integrated moments of Delta f=f_new - f_old..
  bool is_first_fdot_integ_write_call; // flag for integrated moments dynvec written first time

  struct gkyl_array_integrate* integ_wfsq_op; // Operator to integrate w*f^2.
  double *L2norm_local, *L2norm_global; // L2norm in local MPI process and across the communicator.
  gkyl_dynvec L2norm; // L2 norm.
  bool is_first_L2norm_write_call; // flag for L2norm dynvec written first time

  gkyl_dg_updater_gyrokinetic *slvr; // Gyrokinetic solver.
  struct gkyl_dg_eqn *eqn_gyrokinetic; // Gyrokinetic equation object.
  
  int num_periodic_dir; // Number of periodic directions.
  int periodic_dirs[3]; // List of periodic directions.
  bool bc_is_np[3]; // Whether BC is nonperiodic.

  // Boundary conditions on lower/upper edges in each direction.
  struct gkyl_gyrokinetic_bc lower_bc[3], upper_bc[3];
  // gyrokinetic sheath boundary conditions
  struct gkyl_bc_sheath_gyrokinetic *bc_sheath_lo;
  struct gkyl_bc_sheath_gyrokinetic *bc_sheath_up;
  // Pointers to updaters that apply (non-sheath) BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges, valid (i.e. volume>0) in ranks abutting boundaries.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];
  // GK_IWL sims need SOL ghost and skin ranges.
  struct gkyl_range lower_skin_par_sol, lower_ghost_par_sol;
  struct gkyl_range upper_skin_par_sol, upper_ghost_par_sol;
  // GK IWL sims need a core range extended in z, and a TS BC updater.
  struct gkyl_range local_par_ext_core;
  struct gkyl_bc_twistshift *bc_ts_lo, *bc_ts_up;

  struct gk_proj proj_init; // Projector for initial conditions.

  struct gk_source src; // Plasma source.

  // Boundary fluxes used for other solvers.
  struct gk_boundary_fluxes bflux_solver;

  // Boundary fluxes used for diagnostics.
  struct gk_boundary_fluxes bflux_diag;

  struct gk_lte lte; // Object constructing LTE distributions.

  // Collisions.
  union {
    struct {
      struct gk_lbo_collisions lbo; // LBO collisions object
    };
    struct {
      struct gk_bgk_collisions bgk; // BGK collisions object
    };
  }; 

  struct gk_react react; // Object for reactions with charged species.
  struct gk_react react_neut; // Object for reactions with neutral species.

  struct gk_rad_drag rad; // Radiation object.

  // Gyrokinetic diffusion.
  bool has_diffusion; // Flag to indicate there is applied diffusion.
  struct gkyl_array *diffD; // Array for diffusion tensor.
  struct gkyl_dg_updater_diffusion_gyrokinetic *diff_slvr; // Gyrokinetic diffusion equation solver.

  // Updater that enforces positivity by shifting f.
  struct gkyl_positivity_shift_gyrokinetic *pos_shift_op;
  struct gkyl_array *ps_delta_m0; // Number density of the positivity shift.
  struct gkyl_array *ps_delta_m0s_tot; // Density of total positivity shift (like-species).
  struct gkyl_array *ps_delta_m0r_tot; // Density of total positivity shift (other species).
  struct gk_species_moment ps_moms; // Positivity shift diagnostic moments.
  gkyl_dynvec ps_integ_diag; // Integrated moments of the positivity shift.
  bool is_first_ps_integ_write_call; // Flag first time writing ps_integ_diag.

  // Pointer to various functions selected at runtime.
  double (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs);
  double (*rhs_implicit_func)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs, double dt);
  void (*bc_func)(gkyl_gyrokinetic_app *app, const struct gk_species *species,
    struct gkyl_array *f);
  void (*release_func)(const gkyl_gyrokinetic_app* app, const struct gk_species *s);
  void (*step_f_func)(struct gkyl_array* out, double dt, const struct gkyl_array* inp); 
  void (*combine_func)(struct gkyl_array *out, double c1,
    const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
    const struct gkyl_range *rng);
  void (*copy_func)(struct gkyl_array *out, const struct gkyl_array *inp,
    const struct gkyl_range *range);
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);
  void (*write_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks);
  void (*calc_L2norm_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);
  void (*write_L2norm_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks);
  void (*calc_int_mom_dt_func)(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom);

  // Quantities used for FLR model:
  struct gkyl_array *m0_gyroavg; // Gyroaveraged particle density.
  struct gkyl_array *flr_rhoSqD2; // Laplacian weight in FLR operator.
  struct gkyl_array *flr_kSq; // Field multiplying phi in FLR operator.
  struct gkyl_deflated_fem_poisson *flr_op; // Helmholtz solver to invert FLR operator.
  // Pointer to function that performs the gyroaverage.
  void (*gyroaverage)(gkyl_gyrokinetic_app *app, struct gk_species *species,
    struct gkyl_array *field_in, struct gkyl_array *field_gyroavg);

  double *omega_cfl; // Maximum Omega_CFL in this MPI process.
  double *m0_max; // Maximum number density in this MPI process.
};

// neutral species data
struct gk_neut_species {
  struct gkyl_gyrokinetic_neut_species info; // data for neutral species

  struct gkyl_basis basis; // phase-space basis

  // pointer to basis on device
  // (points to host structs if not on GPU)
  struct gkyl_basis *basis_on_dev; 
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_velocity_map *vel_map; // Velocity mapping objects.

  struct gkyl_array *f, *f1, *fnew; // arrays for updates

  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // fixed buffers for time independent BCs
  struct gkyl_array *bc_buffer_lo_recyc, *bc_buffer_up_recyc; // fixed buffers for recycle BCs 

  struct gkyl_array *f_host; // host copy for use IO and initialization

  enum gkyl_field_id field_id; // type of field equation (always GKYL_FIELD_NULL)
  enum gkyl_model_id model_id; // type of Vlasov equation (always GKYL_MODEL_CANONICAL_PB)

  struct gkyl_array *hamil; // Specified hamiltonian function for canonical poisson bracket
  struct gkyl_array *hamil_host; // Host side hamiltonian array for intial projection
  struct gkyl_array *alpha_surf; // array for surface phase space flux (v^i = v . e^i)
  struct gkyl_array *sgn_alpha_surf; // array for the sign of the surface phase space flux at quadrature points
                                     // utilized for numerical flux function
                                     // F = alpha_surf/2 ( (f^+ + f^-) - sign_alpha_surf*(f^+ - f^-) )
  struct gkyl_array *const_sgn_alpha; // boolean array for if the surface phase space flux is single signed
                                      // if true, numerical flux function inside kernels simplifies to
                                      // F = alpha_surf*f^- (if sign_alpha_surf = 1), 
                                      // F = alpha_surf*f^+ (if sign_alpha_surf = -1)

  struct gk_species_moment m0; // for computing density
  struct gk_species_moment integ_moms; // integrated moments
  struct gk_species_moment *moms; // diagnostic moments
  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  gkyl_dg_updater_vlasov *slvr; // Vlasov solver 
  struct gkyl_dg_eqn *eqn_vlasov; // Vlasov equation object

  // boundary fluxes
  struct gk_boundary_fluxes bflux;

  // recycling wall boundary conditions
  struct gk_recycle_wall bc_recycle_lo;
  struct gk_recycle_wall bc_recycle_up;
  bool recyc_lo;
  bool recyc_up;
  
  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
  bool bc_is_np[3]; // whether BC is nonperiodic.

  // boundary conditions on lower/upper edges in each direction  
  struct gkyl_gyrokinetic_bc lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges, valid (i.e. volume>0) in ranks abutting boundaries.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];

  struct gk_proj proj_init; // projector for initial conditions

  struct gk_source src; // applied source

  struct gk_lte lte; // object needed for the lte equilibrium

  // collisions
  union {
    struct {
      struct gk_bgk_collisions bgk; // BGK collisions object
    };
  }; 

  struct gk_react react_neut; // reaction object

  double *omega_cfl;

  // pointer to rhs functions
  double (*rhs_func)(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs);
  double (*rhs_implicit_func)(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
    const struct gkyl_array *fin, struct gkyl_array *rhs, double dt);
  void (*bc_func)(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
    struct gkyl_array *f);
  void (*release_func)(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s);
  void (*step_f_func)(struct gkyl_array* out, double dt, const struct gkyl_array* inp); 
  void (*combine_func)(struct gkyl_array *out, double c1,
    const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
    const struct gkyl_range *rng);
  void (*copy_func)(struct gkyl_array *out, const struct gkyl_array *inp,
    const struct gkyl_range *range);
  void (*write_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);
  void (*write_mom_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);
  void (*calc_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);
  void (*write_integrated_mom_func)(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);
};

// field data
struct gk_field {
  struct gkyl_gyrokinetic_field info; // data for field

  enum gkyl_gkfield_id gkfield_id;

  bool update_field; // Are we updating the field?.
  bool calc_init_field; // Whether to compute the t=0 field.

  struct gkyl_job_pool *job_pool; // Job pool  
  // arrays for local charge density, global charge density, and global smoothed (in z) charge density
  struct gkyl_array *rho_c;
  struct gkyl_array *rho_c_global_dg;
  struct gkyl_array *rho_c_global_smooth; 
  struct gkyl_array *phi_fem, *phi_smooth; // arrays for updates

  struct gkyl_array *phi_host;  // host copy for use IO and initialization

  bool init_phi_pol; // Whether to use the initial user polarization phi.
  struct gkyl_array *phi_pol; // Initial polarization density potential.

  struct gkyl_range global_sub_range; // sub range of intersection of global range and local range
                                      // for solving subset of Poisson solves with parallelization in z

  // organization of the different equation objects and the required data and solvers
  union {
    struct {
      struct gkyl_ambi_bolt_potential *ambi_pot;
      struct gkyl_array *sheath_vals[2*GKYL_MAX_CDIM];
    };
    // EM GK model
    struct {
      struct gkyl_array *apar_fem; // array for A_parallel
      struct gkyl_array *apardot_fem; // array for d/dt A_parallel
    };
  };

  double es_energy_fac_1d; 
  struct gkyl_array *es_energy_fac; 
  bool is_dirichletvar; // Whether user provided spatially varying phi BCs.
  struct gkyl_array *phi_bc; // Spatially varying BC.
  struct gkyl_array *epsilon; // Polarization weight including geometric factors.
  struct gkyl_array *kSq; 

  struct gkyl_fem_parproj *fem_parproj; // FEM smoother for projecting DG functions onto continuous FEM basis
                                        // weight*phi_{fem} = phi_{dg} 
  struct gkyl_fem_parproj *fem_parproj_sol;
  struct gkyl_fem_parproj *fem_parproj_core;

  struct gkyl_deflated_fem_poisson *deflated_fem_poisson; // poisson solver which solves
  struct gkyl_fem_poisson_perp *fem_poisson; // poisson solver which solves
                                             // - nabla . (epsilon * nabla phi) - kSq * phi = rho

  // Objects needed for FLR effects.
  bool use_flr;
  void (*invert_flr)(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi);
  struct gkyl_array *flr_rhoSq_sum; // Laplacian weight in FLR operator.
  struct gkyl_array *flr_kSq; // Field multiplying phi in FLR operator.
  struct gkyl_deflated_fem_poisson *flr_op; // Helmholtz solver to invert FLR operator.

  struct gkyl_array_integrate *calc_em_energy;
  double *em_energy_red, *em_energy_red_global; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components
  bool is_first_energy_write_call; // flag for energy dynvec written first time

  double *em_energy_red_old, *em_energy_red_new; // memory for use in GPU reduction of old EM energy.
  gkyl_dynvec integ_energy_dot; // d/dt of integrated energy components.
  bool is_first_energy_dot_write_call; // flag for d(energy)/dt dynvec written first time

  bool has_phi_wall_lo; // flag to indicate there is biased wall potential on lower wall
  bool phi_wall_lo_evolve; // flag to indicate biased wall potential on lower wall is time dependent
  struct gkyl_array *phi_wall_lo; // biased wall potential on lower wall
  struct gkyl_array *phi_wall_lo_host; // host copy for use in IO and projecting
  gkyl_eval_on_nodes *phi_wall_lo_proj; // projector for biased wall potential on lower wall 

  bool has_phi_wall_up; // flag to indicate there is biased wall potential on upper wall
  bool phi_wall_up_evolve; // flag to indicate biased wall potential on upper wall is time dependent
  struct gkyl_array *phi_wall_up; // biased wall potential on upper wall
  struct gkyl_array *phi_wall_up_host; // host copy for use in IO and projecting
  gkyl_eval_on_nodes *phi_wall_up_proj; // projector for biased wall potential on upper wall 

  // Core and SOL ranges for IWL sims.
  struct gkyl_range global_core, global_ext_core, global_sol, global_ext_sol;


  // Pointer to function that computes the time rate of change of the energy.
  void (*calc_energy_dt_func)(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced);

  // Objects used in IWL simulations and TS BCs.
  struct gkyl_range local_par_ext_core; // Core range extended in parallel direction
  struct gkyl_bc_twistshift *bc_T_LU_lo; // TS BC updater.
  // Objects used by the skin surface to ghost (SSFG) operator.
  struct gkyl_range lower_skin_core, lower_ghost_core;
  struct gkyl_range upper_skin_core, upper_ghost_core;
  struct gkyl_skin_surf_from_ghost *ssfg_lo;
  
  // Pointer to function for the twist-and-shift BCs.
  void (*enforce_zbc) (const gkyl_gyrokinetic_app *app, const struct gk_field *field, struct gkyl_array *finout);
};

// Gyrokinetic object: used as opaque pointer in user code.
struct gkyl_gyrokinetic_app {
  char name[128]; // name of app
  struct gkyl_job_pool *job_pool; // Job pool
  
  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  double tcurr; // current time
  double cfl; // CFL number
  double cfl_omegaH; // CFL number used for omega_H.
  double bmag_ref; // Reference magnetic field

  bool use_gpu; // should we use GPU (if present)

  bool enforce_positivity; // Whether to enforce f>0.

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
    
  struct gkyl_rect_grid grid; // config-space grid
  struct gkyl_range local, local_ext; // local, local-ext conf-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges  
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  // Global skin/ghost ranges, valid (i.e. volume>0) in ranks abutting boundaries.
  struct gkyl_range global_lower_skin[GKYL_MAX_DIM];
  struct gkyl_range global_lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range global_upper_skin[GKYL_MAX_DIM];
  struct gkyl_range global_upper_ghost[GKYL_MAX_DIM];

  struct gkyl_basis basis; // conf-space basis
  
  struct gkyl_rect_decomp *decomp; // Volume decomposition object.
  struct gkyl_comm *comm; // Volume communicator object for conf-space arrays.

  struct gkyl_rect_decomp *decomp_plane[GKYL_MAX_CDIM]; // Volume decomposition object.
  struct gkyl_comm *comm_plane[GKYL_MAX_CDIM]; // Volume communicator object for conf-space arrays.

  // pointer to basis on device
  // (points to host structs if not on GPU)
  struct gkyl_basis *basis_on_dev; 

  struct gk_geometry *gk_geom;
  struct gkyl_array *jacobtot_inv_weak; // 1/(J.B) computed via weak mul and div.
  double omegaH_gf; // Geometry and field model dependent part of omega_H.
  
  struct gkyl_position_map *position_map; // Position mapping object.

  struct gk_field *field; // pointer to field object
  // Pointer to function that computes the fields.
  void (*calc_field_func)(gkyl_gyrokinetic_app* app, double tcurr, const struct gkyl_array *fin[]);

  // species data
  int num_species;
  struct gk_species *species; // data for each species
  struct gkyl_array *ps_delta_m0_ions; // Number density of the total ion positivity shift.
  struct gkyl_array *ps_delta_m0_elcs; // Number density of the total elc positivity shift.
  
  // neutral species data
  int num_neut_species;
  struct gk_neut_species *neut_species; // data for each species

  bool has_implicit_coll_scheme; // Boolean for using implicit bgk scheme (over explicit rk3)

  // pointer to function that takes a single-step of simulation
  struct gkyl_update_status (*update_func)(gkyl_gyrokinetic_app *app, double dt0);

  struct gkyl_gyrokinetic_stat stat; // statistics

  gkyl_dynvec dts; // Record time step over time.
  bool is_first_dt_write_call; // flag for integrated moments dynvec written first time
};

/** gkyl_gyrokinetic_app private API */

/**
 * Create a new array metadata object. It must be freed using
 * gk_array_meta_release.
 *
 * @param meta Gyrokinetic metadata object.
 * @return Array metadata object.
 */
struct gkyl_msgpack_data*
gk_array_meta_new(struct gyrokinetic_output_meta meta);

/**
 * Free memory for array metadata object.
 *
 * @param mt Array metadata object.
 */
void
gk_array_meta_release(struct gkyl_msgpack_data *mt);

/**
 * Return the metadata for outputing gyrokinetic data.
 *
 * @param mt Array metadata object.
 * @return A gyrokinetic metadata object.
 */
struct gyrokinetic_output_meta
gk_meta_from_mpack(struct gkyl_msgpack_data *mt);

/**
 * Allocate a new gyrokinetic app and initialize its conf-space grid and
 * geometry. This method needs to be complemented by
 * gkyl_gyrokinetic_app_new_solver below.
 *
 * @param gk Gyrokinetic input struct.
 * @return A gyrokinetic app object.
 */
gkyl_gyrokinetic_app*
gkyl_gyrokinetic_app_new_geom(struct gkyl_gk *gk);

/**
 * Initialize the rest of the gyrokinetic app solver, after having called
 * the gkyl_gyrokinetic_app_new_geom method.
 *
 * @param gk Gyrokinetic input struct.
 * @param app Gyrokinetic app.
 */
void
gkyl_gyrokinetic_app_new_solver(struct gkyl_gk *gk, gkyl_gyrokinetic_app *app);

/**
 * Find species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Pointer to species with given name. NULL if not found.
 */
struct gk_species* gk_find_species(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Return index of species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of species
 * @return Index of species, -1 if not found
 */
int gk_find_species_idx(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Find neutral species with given name.
 *
 * @param app Top-level app to look into
 * @param nm Name of neutral species
 * @return Pointer to neutral species with given name. NULL if not found.
 */
struct gk_neut_species* gk_find_neut_species(const gkyl_gyrokinetic_app *app, const char *nm);

/**
 * Return index of neutral species in the order it appears in the input.
 *
 * @param app Top-level app to look into
 * @param nm Name of neutral species
 * @return Index of neutral species, -1 if not found
 */
int gk_find_neut_species_idx(const gkyl_gyrokinetic_app *app, const char *nm);

/** gk_species_moment API */

/**
 * Initialize species moment object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param sm Species moment object
 * @param is_integrated Whether to compute the volume integrated moment.
 * @param nm Name string indicating moment type
 */
void gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, const char *nm, bool is_integrated);

/**
 * Calculate moment, given distribution function @a fin.
 * 
 * @param sm Species moment object
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input distribution function array
 */
void gk_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release species moment object.
 *
 * @param app gyrokinetic app object
 * @param sm Species moment object to release
 */
void gk_species_moment_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species_moment *sm);

/** gk_species_radiation API */

/**
 * Initialize species radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param rad Species radiation drag object
 */
void gk_species_radiation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_rad_drag *rad);

/**
 * Compute necessary moments for radiation drag object
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution functions (size num_species)
 * @param fin_neut Input neutral distribution functions (size num_species)
 */
void gk_species_radiation_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species, struct gk_rad_drag *rad, 
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute emissivities 
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution functions (size num_species)
 * @param fin_neut Input neutral distribution functions (size num_species)
 */
void gk_species_radiation_emissivity(gkyl_gyrokinetic_app *app,
  struct gk_species *species, struct gk_rad_drag *rad, 
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute integrated moments of radiation drag object
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution functions (size num_species)
 * @param fin_neut Input neutral distribution functions (size num_species)
 */
void
gk_species_radiation_integrated_moms(gkyl_gyrokinetic_app *app, struct gk_species *species,
				struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute RHS from radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param rad Species radiation drag object
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void gk_species_radiation_rhs(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_rad_drag *rad,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write species radiation drag.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_radiation_write_drag(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Write species radiation emissivity.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_radiation_write_emissivity(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Calculate species radiation integrated moments.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 */
void gk_species_radiation_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Write species radiation integrated moments.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 */
void gk_species_radiation_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Release species radiation drag object.
 *
 * @param app gyrokinetic app object
 * @param rad Species radiation drag object to release
 */
void gk_species_radiation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_rad_drag *rad);

/** gk_species_lbo API */

/**
 * Initialize species LBO collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void gk_species_lbo_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_lbo_collisions *lbo);

/**
 * Initialize species LBO cross-collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param lbo Species LBO object
 */
void gk_species_lbo_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_lbo_collisions *lbo);

/**
 * Compute necessary moments and boundary
 * corrections for LBO collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void gk_species_lbo_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species LBO collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 */
void gk_species_lbo_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lbo_collisions *lbo,
  const struct gkyl_array *fin);

/**
 * Compute RHS from LBO collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param lbo Pointer to LBO
 * @param fin Input distribution function
 * @param rhs On output, the RHS from LBO
 */
void gk_species_lbo_rhs(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lbo_collisions *lbo,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write moments from LBO object.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_lbo_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Release species LBO object.
 *
 * @param app gyrokinetic app object
 * @param lbo Species LBO object to release
 */
void gk_species_lbo_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lbo_collisions *lbo);

/** gk_species_lte API */

/**
 * Initialize species lte object.
 *
 * @param app Gyrokinetic app object
 * @param s Species object 
 * @param lte Species lte object
 * @param corr_inp Input struct with moment correction inputs
 */
void gk_species_lte_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_lte *lte, struct correct_all_moms_inp corr_inp);

/**
 * Compute LTE distribution from input moments
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to species
 * @param lte Pointer to lte object
 * @param moms_lte Input LTE moments
 */
void gk_species_lte_from_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lte *lte,
  const struct gkyl_array *moms_lte);

/**
 * Compute equivalent LTE distribution from input distribution function. 
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to species
 * @param lte Pointer to lte
 * @param fin Input distribution function
 */
void gk_species_lte(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_lte *lte,
  const struct gkyl_array *fin);

/**
 * Write the LTE correction status. 
 *
 * @param app Gyrokinetic app object
 * @param gks Pointer to species
 */
void gk_species_lte_write_max_corr_status(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Release species lte object.
 *
 * @param app gyrokinetic app object
 * @param lte Species lte object to release
 */
void gk_species_lte_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte);

/** gk_species_bgk API */

/**
 * Initialize species BGK collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param bgk Species BGK object
 */
void gk_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_bgk_collisions *bgk);

/**
 * Initialize species BGK cross-collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param bgk Species BGK object
 */
void gk_species_bgk_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_bgk_collisions *bgk);

/**
 * Compute necessary moments for BGK collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void gk_species_bgk_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_bgk_collisions *bgk,
  const struct gkyl_array *fin);

/**
 * Compute necessary moments for cross-species BGK collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void gk_species_bgk_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species,
  struct gk_bgk_collisions *bgk,
  const struct gkyl_array *fin);

/**
 * Compute RHS from BGK collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 * @param rhs On output, the RHS from bgk
 */
void gk_species_bgk_rhs(gkyl_gyrokinetic_app *app,
  struct gk_species *species, struct gk_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write the BGK cross moments. 
 *
 * @param app Gyrokinetic app object
 * @param gks Pointer to species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_bgk_write_cross_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Release species BGK object.
 *
 * @param app gyrokinetic app object
 * @param bgk Species BGK object to release
 */
void gk_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk);

/** gk_species_react API */

/**
 * Initialize species reactions object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param inp Input reaction struct for determining types of reactions
 * @param react Species reaction object
 * @param all_gk Boolean for if the reactions are between only GK species
 */
void gk_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gkyl_gyrokinetic_react inp, struct gk_react *react, bool all_gk);

/**
 * Initialize species reactions "cross-collisions" object
 * for who is reacting with whom
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param react Species react object
 */
void gk_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_react *react);

/**
 * Compute necessary rates and moments for reactions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param react Pointer to react
 * @param f_self Input self distribution function
 * @param fin Input distribution functions (size: num_species)
 * @param fin_neut Input neutral distribution functions (size: num_neut_species)
 */
void gk_species_react_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_species *species, struct gk_react *react,
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute RHS from reactions 
 * (e.g., ionization, charge exchange, recombination, or radiation)
 *
 * @param app gyrokinetic app object
 * @param s Pointer to species
 * @param react Pointer to react
 * @param fin Input distribution function
 * @param rhs On output, the RHS from react (df/dt)
 */
void gk_species_react_rhs(gkyl_gyrokinetic_app *app,
  struct gk_species *s, struct gk_react *react,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write reaction rate.
 *
 * @param app gyrokinetic app object
 * @param gks Pointer to species
 * @param gkr Pointer to react object
 * @param ridx Index for reaction species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_species_react_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, struct gk_react *gkr,
  int ridx, double tm, int frame);

/**
 * Release species react object.
 *
 * @param app gyrokinetic app object
 * @param react Species react object to release
 */
void gk_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react);

/** gk_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Gyrokinetic app object.
 * @param s Species object. 
 * @param bflux Species boundary flux object.
 * @param is_diagnostic Wether this object is used for diagnotics.
 */
void gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_boundary_fluxes *bflux, bool is_diagnostic);

/**
 * Compute boundary flux, either for another solver or for diagnostics. The
 * latter also computes moments of the boundary flux.
 * Note: stores the boundary flux in the ghost cells of rhs
 * The ghost cells are overwritten by apply_bc so computations using
 * boundary fluxes are internal to rhs method (such as integrated flux)
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param bflux Species boundary flux object.
 * @param fin Input distribution function.
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs.
 * @param bflux_out Array of moments of boundary fluxes through every boundary.
 */
void gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs,
  struct gkyl_array **bflux_moms);

/**
 * Clear the boundary fluxes at each boundary.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param bflux Species boundary flux object.
 * @param bflux_in Array of boundary fluxes to clear.
 * @param val Value to set array to.
 */
void
gk_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val);

/**
 * Scale the boundary fluxes at each boundary.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param bflux Species boundary flux object.
 * @param bflux_in Array of boundary fluxes to clear.
 * @param val Value to scale the fluxes by.
 */
void
gk_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_in, double val);

/**
 * Step the diagnotic boundary fluxes forward once.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param dt Time step.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_species_bflux_step_f(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_out, double dt, const struct gkyl_array **bflux_in);

/**
 * Combine the diagnotic boundary fluxes for multi-stage RK stepper.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param fac1 Factor to multiply bflux_in1 by.
 * @param bflux_in1 Array of input boundary fluxes.
 * @param fac2 Factor to multiply bflux_in2 by.
 * @param bflux_in2 Array of input boundary fluxes.
 */
void
gk_species_bflux_combine(gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2);

/**
 * Copy diagnotic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param species Pointer to species.
 * @param bflux Species boundary flux object.
 * @param bflux_out Array of output boundary fluxes.
 * @param bflux_in Array of input boundary fluxes.
 */
void
gk_species_bflux_copy(gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **bflux_fout, struct gkyl_array **bflux_in);

/**
 * Calculate the integrated moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object
 * @param gk_s Species object.
 * @param tm Current simulation time.
 */
void
gk_species_bflux_calc_integrated_mom(gkyl_gyrokinetic_app* app, const struct gk_species *gk_s,
  struct gk_boundary_fluxes *bflux, double tm);

/**
 * Calculate the time integrated, integrated moments
 * of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object
 * @param gk_s Species object.
 * @param tm Current simulation time.
 */
void
gk_species_bflux_calc_voltime_integrated_mom(gkyl_gyrokinetic_app* app,
  const struct gk_species *gk_s, struct gk_boundary_fluxes *bflux, double tm);

/**
 * Write the integrated moments of the diagnostic boundary fluxes.
 *
 * @param app Gyrokinetic app object.
 * @param gks Species object.
 */
void
gk_species_bflux_write_integrated_mom(gkyl_gyrokinetic_app *app,
  const struct gk_species *gks, struct gk_boundary_fluxes *bflux);

/**
 * Release species boundary flux object.
 *
 * @param app Gyrokinetic app object.
 * @param bflux Species boundary flux object to release.
 */
void gk_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux);

/** gk_species_projection API */

/**
 * Initialize species projection object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param inp Input struct for projection (contains functions pointers for type of projection).
 * @param proj Species projection object.
 */
void gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj);

/**
 * Compute species projection.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param proj Species projection object.
 * @param f Output distribution function from projection.
 * @param tm Time for use in projection.
 */
void gk_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_species *species, 
  struct gk_proj *proj, struct gkyl_array *f, double tm);

/**
 * Release species projection object.
 *
 * @param app gyrokinetic app object.
 * @param proj Species projection object to release.
 */
void gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj);

/** gk_species_source API */

/**
 * Initialize species source object.
 *
 * @param app gyrokinetic app object.
 * @param s Species object.
 * @param src Species source object.
 */
void gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_source *src);

/**
 * Compute species applied source term.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param src Species source object.
 * @param tm Time for use in source.
 */
void gk_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_species *species, 
  struct gk_source *src, double tm);

/**
 * Compute RHS contribution from source.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param src Pointer to source.
 * @param fin Input distribution function.
 * @param rhs On output, the distribution function.
 */
void gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write source diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void gk_species_source_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Write source moment diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for source diagnostic.
 * @param frame Output frame.
 */
void gk_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Calc source integrated moment diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 * @param tm Time for source diagnostic.
 */
void gk_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Write source integrated moment diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Pointer to species.
 */
void gk_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Release species source object.
 *
 * @param app gyrokinetic app object.
 * @param src Species source object to release.
 */
void gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_species API */

/**
 * Initialize species.
 *
 * @param gk Input gk data.
 * @param app gyrokinetic app object.
 * @param s On output, initialized species object.
 */
void gk_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param t0 Time for use in ICs.
 */
void gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute the part of the species initial conditions that depends on other
 * species.
 *
 * @param app gyrokinetic app object.
 * @param species Species object.
 * @param t0 Time for use in ICs.
 */
void gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute RHS from species distribution function
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from the species object.
 * @return Maximum stable time-step.
 */
double gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Compute the *implicit* RHS from species distribution function
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param fin Input distribution function.
 * @param rhs On output, the RHS from the species object.
 * @param dt timestep size (used in the implcit coef.).
 * @return Maximum stable time-step.
 */
double gk_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt);

/**
 * Scale and accumulate for forward euler method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param dt Timestep.
 * @param inp Input array.
 */
void gk_species_step_f(struct gk_species *species, struct gkyl_array* out, double dt,
  const struct gkyl_array* inp);

/**
 * Combine for rk3 method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param c1 Scaling factor.
 * @param arr1 Input array.
 * @param c2 Scaling factor.
 * @param arr2 Input array.
 * @param rng Range.
 */
void gk_species_combine(struct gk_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng);

/**
 * Copy for rk3 method.
 *
 * @param species Pointer to species.
 * @param out Output array.
 * @param inp Input array.
 * @param range Range.
 */
void gk_species_copy_range(struct gk_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Apply BCs to dynamic species distribution function.
 *
 * @param app gyrokinetic app object.
 * @param species Pointer to species.
 * @param f Field to apply BCs.
 */
void gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers.
 */
void gk_species_coll_tm(gkyl_gyrokinetic_app *app);

/**
 * Fill stat object in app with total number of iterations
 * used to correct moments in LTE projection object.
 * Also fills stat object with number of times correction object called. 
 *
 * @param app App object to update stat timers.
 */
void gk_species_n_iter_corr(gkyl_gyrokinetic_app *app);

/**
 * Fill stat object in app with collisionless timers.
 *
 * @param app App object to update stat timers.
 */
void gk_species_tm(gkyl_gyrokinetic_app *app);

/**
 * Species write function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_species_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Species moment write function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 * @param frame simulation frame.
 */
void gk_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame);

/**
 * Species calc integrated moment function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 */
void gk_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Species write integrated moment function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 */
void gk_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Species calc L2norm function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param tm simulation time.
 */
void gk_species_calc_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm);

/**
 * Species write L2norm function.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 */
void gk_species_write_L2norm(gkyl_gyrokinetic_app* app, struct gk_species *gks);

/**
 * Calculate the integrated moments divided by dt, for particle/energy balance
 * diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param gks Species object.
 * @param dt Time step.
 * @param fdot_int_mom Integrated moment divided by dt (not yet reduced over comm).
 */
void
gk_species_calc_int_mom_dt(gkyl_gyrokinetic_app* app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom);

/**
 * Delete resources used in species.
 *
 * @param app gyrokinetic app object.
 * @param species Species object to delete.
 */
void gk_species_release(const gkyl_gyrokinetic_app* app, const struct gk_species *s);

/** gk_neut_species_moment API */

/**
 * Initialize neutral species moment object.
 *
 * @param app gyrokinetic app object
 * @param s Neutral species object 
 * @param sm Neutral species moment object
 * @param nm Name string indicating moment type
 */
void gk_neut_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, const char *nm);

/**
 * Calculate neutral species moment, given input neutral distribution function fin.
 * 
 * @param sm Neutral species moment object
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param fin Input neutral distribution function array
 */
void gk_neut_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin);

/**
 * Release neutral species moment object.
 *
 * @param app gyrokinetic app object
 * @param sm Neutral species moment object to release
 */
void gk_neut_species_moment_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species_moment *sm);

/** gk_neut_species_lte API */

/**
 * Initialize species lte object.
 *
 * @param app Gyrokinetic app object
 * @param s Neutral species object 
 * @param lte Neutral species lte object
 * @param corr_inp Input struct with moment correction inputs
 */
void gk_neut_species_lte_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_lte *lte, struct correct_all_moms_inp corr_inp);

/**
 * Compute LTE distribution from input moments
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to neutral species
 * @param lte Pointer to lte object
 * @param moms_lte Input LTE moments
 */
void gk_neut_species_lte_from_moms(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species,
  struct gk_lte *lte,
  const struct gkyl_array *moms_lte);

/**
 * Compute equivalent LTE distribution from input distribution function. 
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to neutral species
 * @param lte Pointer to lte
 * @param fin Input distribution function
 */
void gk_neut_species_lte(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species,
  struct gk_lte *lte,
  const struct gkyl_array *fin);

/**
 * Write the LTE correction status for the neutral species. 
 *
 * @param app Gyrokinetic app object
 * @param gk_ns Pointer to neutral species
 */
void gk_neut_species_lte_write_max_corr_status(gkyl_gyrokinetic_app* app, struct gk_neut_species *gk_ns);

/**
 * Release species lte object.
 *
 * @param app gyrokinetic app object
 * @param lte Neutral species lte object to release
 */
void gk_neut_species_lte_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte);

/** gk_neut_species_bgk API */

/**
 * Initialize neutral species BGK collisions object.
 *
 * @param app gyrokinetic app object
 * @param s Neutral species object 
 * @param bgk Neutral species BGK object
 */
void gk_neut_species_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_bgk_collisions *bgk);

/**
 * Compute necessary moments for BGK collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 */
void gk_neut_species_bgk_moms(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species,
  struct gk_bgk_collisions *bgk,
  const struct gkyl_array *fin);

/**
 * Compute RHS from BGK collisions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param bgk Pointer to BGK
 * @param fin Input distribution function
 * @param rhs On output, the RHS from bgk
 */
void gk_neut_species_bgk_rhs(gkyl_gyrokinetic_app *app,
  struct gk_neut_species *species, struct gk_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species BGK object.
 *
 * @param app gyrokinetic app object
 * @param bgk Neutral species BGK object to release
 */
void gk_neut_species_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_bgk_collisions *bgk);

/** gk_neut_species_boundary_fluxes API */

/**
 * Initialize species boundary flux object.
 *
 * @param app Gyrokinetic app object
 * @param s Species object 
 * @param bflux Species boundary flux object
 */
void gk_neut_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_boundary_fluxes *bflux);

/**
 * Compute boundary flux 
 * Note: stores the boundary flux in the ghost cells of rhs
 * The ghost cells are overwritten by apply_bc so computations using
 * boundary fluxes are internal to rhs method (such as integrated flux)
 *
 * @param app Gyrokinetic app object
 * @param species Pointer to species
 * @param bflux Species boundary flux object
 * @param fin Input distribution function
 * @param rhs On output, the boundary fluxes stored in the ghost cells of rhs
 */
void gk_neut_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species boundary flux object.
 *
 * @param app Gyrokinetic app object
 * @param bflux Species boundary flux object to release
 */
void gk_neut_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux);

/** gk_neut_species_recycle API **/


/**
 * Initialize recycling object.
 *
 * @param app Gyrokinetic app object
 * @param recyc Recycling bc object
 * @param dir Direction for BC (x, y, or z)
 * @param edge Edge for BC (lower/upper)
 * @param ctx Input params for recycling BCs
 * @param f0 Maxwellian distf to scale in ghost
 * @param s Gk_neut_species to apply BCs for
 * @param use_gpu Boolean for using GPUs
 */
void gk_neut_species_recycle_init(struct gkyl_gyrokinetic_app *app, struct gk_recycle_wall *recyc,
  int dir, enum gkyl_edge_loc edge, void *ctx, struct gkyl_array *f0, struct gk_neut_species *s, bool use_gpu);

/**
 * Initialize recycling cross moments.
 *
 * @param app Gyrokinetic app object
 * @param s Gk_neut_species to apply BCs for
 * @param recyc Recycling bc object
 */
void gk_neut_species_recycle_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc);

/**
 * Apply recycling BCs.
 *
 * @param app Gyrokinetic app object
 * @param recyc Recycling bc object
 * @param s Gk_neut_species to apply BCs for
 * @param fout Gk_neut_species distf
 */
void gk_neut_species_recycle_apply_bc(struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc,
  const struct gk_neut_species *s, struct gkyl_array *fout);


/**
 * Write recycle flux diagnostics.
 *
 * @param app Gyrokinetic app object
 * @param s Gk_neut_species to apply BCs for
 * @param recyc Recycling bc object
 * @param tm Simulation time 
 * @param frame Simulation frame
 */
void
gk_neut_species_recycle_write_flux(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc, double tm, int frame);

/**
 * Release recycle BC object.
 *
 * @param app Gyrokinetic app object
 * @param recyc Recycling bc object
 */
void gk_neut_species_recycle_release(const struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc);

/** gk_neut_species_react API */

/**
 * Initialize neutral species reactions object.
 *
 * @param app gyrokinetic app object
 * @param s Neutral species object 
 * @param inp Input reaction struct for determining types of reactions
 * @param react Neutral species reaction object
 */
void gk_neut_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gkyl_gyrokinetic_react inp, struct gk_react *react);

/**
 * Initialize neutral species reactions "cross-collisions" object
 * for who is reacting with whom
 *
 * @param app gyrokinetic app object
 * @param s Neutral species object 
 * @param react Neutral species react object
 */
void gk_neut_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_react *react);

/**
 * Compute necessary rates and moments for reactions
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param react Pointer to react
 * @param fin Input distribution functions (size: num_species)
 * @param fin_neut Input neutral distribution functions (size: num_neut_species)
 */
void gk_neut_species_react_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species,
  struct gk_react *react,
  const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

/**
 * Compute RHS from reactions for neutrals
 * (e.g., ionization, charge exchange, recombination, or radiation)
 *
 * @param app gyrokinetic app object
 * @param s Pointer to neutral species
 * @param react Pointer to react
 * @param fin Input neutral distribution function
 * @param rhs On output, the neutral RHS from react (df/dt)
 */
void gk_neut_species_react_rhs(gkyl_gyrokinetic_app *app,
  struct gk_neut_species *s, struct gk_react *react,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write neutral reaction rate.
 *
 * @param app gyrokinetic app object
 * @param gkns Pointer to neutral species
 * @param gkr Pointer to react object
 * @param ridx Index for reaction species
 * @param tm Simulation time
 * @param frame Simulation output frame
 */
void gk_neut_species_react_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
  int ridx, double tm, int frame);

/**
 * Release neutral species react object.
 *
 * @param app gyrokinetic app object
 * @param react Neutral species react object to release
 */
void gk_neut_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react);

/** gk_neut_species_projection API */

/**
 * Initialize neutral species projection object.
 *
 * @param app gyrokinetic app object
 * @param s Neutral species object 
 * @param inp Input struct for projection (contains functions pointers for type of projection)
 * @param proj Neutral species projection object
 */
void gk_neut_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj);

/**
 * Compute neutral species projection
 *
 * @param app gyrokinetic app object
 * @param species Neutral species object
 * @param proj Neutral species projection object
 * @param f Output Neutral distribution function from projection
 * @param tm Time for use in projection
 */
void gk_neut_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, 
  struct gk_proj *proj, struct gkyl_array *f, double tm);

/**
 * Release neutral species projection object.
 *
 * @param app gyrokinetic app object
 * @param proj Neutral species projection object to release
 */
void gk_neut_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj);

/** gk_neut_species_source API */

/**
 * Initialize neutral species source object.
 *
 * @param app gyrokinetic app object
 * @param s Neutral species object 
 * @param src Neutral species source object
 */
void gk_neut_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_source *src);

/**
 * Compute Neutral species applied source term
 *
 * @param app gyrokinetic app object
 * @param species Neutral species object
 * @param src Neutral species source object
 * @param tm Time for use in source
 */
void gk_neut_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, 
  struct gk_source *src, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app gyrokinetic app object
 * @param species Pointer to Neutral species
 * @param src Pointer to source
 * @param fin Input neutral distribution function
 * @param rhs On output, the incremented rhs (df/dt)
 */
void gk_neut_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Write neutral source diagnostics
 *
 * @param app gyrokinetic app object
 * @param gkns Pointer to species
 * @param tm Time for source diagnostic
 * @param frame Output frame
 */
void gk_neut_species_source_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Write neutral source moment diagnostics
 *
 * @param app gyrokinetic app object
 * @param gkns Pointer to species
 * @param tm Time for source diagnostic
 * @param frame Output frame
 */
void gk_neut_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Calc neutral source integrated moment diagnostics
 *
 * @param app gyrokinetic app object
 * @param gkns Pointer to species
 * @param tm Time for source diagnostic
 */
void gk_neut_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);

/**
 * Write neutral source integrated moment diagnostics
 *
 * @param app gyrokinetic app object
 * @param gkns Pointer to species
 * @param tm Time for source diagnostic
 */
void gk_neut_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

/**
 * Release Neutral species source object.
 *
 * @param app gyrokinetic app object
 * @param src Neutral species source object to release
 */
void gk_neut_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_neut_species API */
/**
 * Initialize neutral species.
 *
 * @param gk Input gk data
 * @param app gyrokinetic app object
 * @param s On output, initialized neutral species object
 */
void gk_neut_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s);

/**
 * Compute neutral species initial conditions.
 *
 * @param app gyrokinetic app object
 * @param species Neutral species object
 * @param t0 Time for use in ICs
 */
void gk_neut_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double t0);

/**
 * Compute RHS from neutral species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param fin Input distribution function
 * @param rhs On output, the RHS from the neutral species object (df/dt)
 * @return Maximum stable time-step
 */
double gk_neut_species_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Compute the *implicit* RHS from neutral species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param fin Input distribution function
 * @param rhs On output, the RHS from the species object
 * @param dt timestep size (used in the implcit coef.)
 * @return Maximum stable time-step
 */
double gk_neut_species_rhs_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, double dt);

/**
 * Apply BCs to neutral species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param f Field to apply BCs
 */
void gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with total number of iterations
 * used to correct moments in LTE projection object for neutral species.
 * Also fills stat object with number of times correction object called. 
 *
 * @param app App object to update stat timers
 */
void gk_neut_species_n_iter_corr(gkyl_gyrokinetic_app *app);

/**
 * Fill stat object in app with collisionless timers for neutral species.
 *
 * @param app App object to update stat timers
 */
void gk_neut_species_tm(gkyl_gyrokinetic_app *app);

/**
 * Scale and accumulate for forward euler method.
 *
 * @param species Pointer to neutral species
 * @param out Output array
 * @param dt Timestep
 * @param inp Input array
 */
void gk_neut_species_step_f(struct gk_neut_species *species, struct gkyl_array* out, double dt,
  const struct gkyl_array* inp);

/**
 * Combine for rk3 method.
 *
 * @param species Pointer to species
 * @param out Output array
 * @param c1 Scaling factor
 * @param arr1 Input array
 * @param c2 Scaling factor
 * @param arr2 Input array
 * @param rng Range
 */
void gk_neut_species_combine(struct gk_neut_species *species, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng);

/**
 * Copy for rk3 method.
 *
 * @param species Pointer to species
 * @param out Output array
 * @param inp Input array
 * @param range Range
 */
void gk_neut_species_copy_range(struct gk_neut_species *species, struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Species write function.
 *
 * @param app gyrokinetic app object
 * @param gkns Neutral species object
 * @param tm simulation time
 * @param frame simulation frame
 */
void gk_neut_species_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Species moment write function.
 *
 * @param app gyrokinetic app object
 * @param gkns Neutral species object
 * @param tm simulation time
 * @param frame simulation frame
 */
void gk_neut_species_write_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

/**
 * Species calc integrated moment function.
 *
 * @param app gyrokinetic app object
 * @param gkns Neutral species object
 * @param tm simulation time
 */
void gk_neut_species_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);

/**
 * Species write integrated moment function.
 *
 * @param app gyrokinetic app object
 * @param gkns Neutral species object
 */
void gk_neut_species_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

/**
 * Delete resources used in neutral species.
 *
 * @param app gyrokinetic app object
 * @param species Neutral species object to delete
 */
void gk_neut_species_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s);

/** gk_field API */

/**
 * Create new field object
 *
 * @param gk Input gk data
 * @param app gyrokinetic app object
 * @return Newly created field
 */
struct gk_field* gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app);

/**
 * Compute biased wall potentials 
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 * @param tm Time to compute biased wall potentials at
 */
void gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm);

/**
 * Accumulate charge density for Poisson solve
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 * @param fin[] Input distribution function (num_species size)
 */
void gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[]);

/**
 * Compute sheath values for use in ambipolar potential model 
 * Computes the ion flux and then the ion density and sheath potential
 * at the entrance of the sheath assuming adiabatic electrons
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 */
void gk_field_calc_ambi_pot_sheath_vals(gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Compute EM field 
 *
 * @param app gyrokinetic app object
 * @param field Pointer to field
 * @param em Output field
 */
void gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Read the field from a file.
 *
 * @param app Gyrokinetic app.
 * @param inp Input struct with importing parameters (and file name).
 */
void gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp);

/**
 * Project the initial field using a user provided function.
 *
 * @param app Gyrokinetic app.
 */
void gk_field_project_init(struct gkyl_gyrokinetic_app *app);

/**
 * Compute field energy diagnostic.
 *
 * @param app gyrokinetic app object.
 * @param tm Time at which diagnostic is computed.
 * @param field Pointer to field.
 */
void gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field);

/**
 * Compute field energy divided by dt for energy balance diagnostics.
 *
 * @param app gyrokinetic app object.
 * @param field Pointer to field.
 * @param dt Time step.
 * @param energy_reduced Integrated field energy (single element double array).
 */
void gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced);

/**
 * Release resources allocated by field
 *
 * @param app gyrokinetic app object
 * @param f Field object to release
 */
void gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f);

/** Time stepping API */

/**
 * Compute the gyrokinetic fields.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param fin Array of distribution functions (one for each species) .
 */
void gyrokinetic_calc_field(gkyl_gyrokinetic_app* app, double tcurr, const struct gkyl_array *fin[]);

/**
 * Compute the gyrokinetic fields and apply boundary conditions.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param distf Array of distribution functions (for each charged species).
 * @param distf_neut Array of distribution functions (for each neutral species).
 */
void gyrokinetic_calc_field_and_apply_bc(gkyl_gyrokinetic_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[]);

/**
 * Compute the RHS of the gyrokinetic equation (df/dt) and the minimum time
 * step it requires for stability based on the CFL constraint.
 *
 * @param app Gyrokinetic app.
 * @param tcurr Current simulation time.
 * @param dt Suggested time step.
 * @param fin Input array of charged-species distribution functions.
 * @param fout Output array of charged-species distribution functions.
 * @param bflux_out Output array of charged-species boundary fluxes.
 * @param fin_neut Input array of neutral-species distribution functions.
 * @param fout_neut Output array of neutral-species distribution functions.
 * @param st Time stepping status object.
 */
void gyrokinetic_rhs(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], struct gkyl_array *fout[], struct gkyl_array **bflux_out[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], 
  struct gkyl_update_status *st); 

/**
 * Take time-step using the RK3 method. Also sets the status object
 * which has the actual and suggested dts used. These can be different
 * from the actual time-step.
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
struct gkyl_update_status gyrokinetic_update_ssp_rk3(gkyl_gyrokinetic_app* app, double dt0);

/**
 * Take time-step of the (BGK) collision operator using a first order implicit method. 
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
void gyrokinetic_update_implicit_coll(gkyl_gyrokinetic_app *app,  double dt0);

/**
 * Take time-step using a first order operator split combining 
 * the RK3 method for the collisionless advection with a first order implicit
 * method for BGK collisions. The first order implicit step occurs *after* the 
 * RK3 step and utilizes the stable RK3 time step. Also sets the status object
 * which has the actual and suggested dts used. These can be different
 * from the actual time-step.
 *
 * @param app Gyrokinetic app.
 * @param dt0 Suggessted time step.
 */
struct gkyl_update_status gyrokinetic_update_op_split(gkyl_gyrokinetic_app *app,  double dt0);
