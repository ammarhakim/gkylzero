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
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_correct_maxwellian_gyrokinetic.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gyrokinetic_vars.h>
#include <gkyl_dg_calc_vlasov_gen_geo_vars.h>
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
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_cross_bgk.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_null_pool.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_cross_calc.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_proj_bimaxwellian_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_radiation_read.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_spitzer_coll_freq.h>
#include <gkyl_tok_geo.h>
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

// data for gyrokinetic moments
struct gk_species_moment {
  struct gk_geometry *gk_geom; // geometry struct for dividing moments by Jacobian
  struct gkyl_dg_bin_op_mem *mem_geo; // memory needed in dividing moments by Jacobian
  bool is_integrated; // boolean for if computing integrated moments 
                      // integrated moments do not need to divide by Jacobian since
                      // the inverse Jacobian is already included in the computation
  int num_mom; // number of moments 

  struct gkyl_dg_updater_moment *mcalc; // moment update

  struct gkyl_array *marr; // array to moment data
  struct gkyl_array *marr_host; // host copy (same as marr if not on GPUs)
};

struct gk_rad_drag {  
  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with
  struct gk_neut_species *collide_with_neut[GKYL_MAX_SPECIES]; // pointers to neutral cross-species we collide with
  int collide_with_idx[2*GKYL_MAX_SPECIES]; // index of species we collide with
  bool is_neut_species[2*GKYL_MAX_SPECIES]; // Flag of whether neutral or gk species
  
  // drag coefficients in vparallel and mu for each species being collided with
  struct gkyl_array *vnu_surf[2*GKYL_MAX_SPECIES]; 
  struct gkyl_array *vnu[2*GKYL_MAX_SPECIES]; 
  struct gkyl_array *vsqnu_surf[2*GKYL_MAX_SPECIES]; 
  struct gkyl_array *vsqnu[2*GKYL_MAX_SPECIES]; 
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars[2*GKYL_MAX_SPECIES]; 

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

  double vtsq_min; // Smallest vtsq that radiation is calculated
  struct gkyl_array *prim_moms;
  struct gkyl_array *boundary_corrections; // boundary corrections
  struct gkyl_array *vtsq;
  
  gkyl_prim_lbo_calc *coll_pcalc; // primitive moment calculator to find te
  struct gkyl_mom_calc_bcorr *bcorr_calc; // LBO boundary corrections calculator for prim_lbo_calc
  struct gk_species_moment lab_moms; // moments needed for te (single array includes Zeroth, First, and Second moment)

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
  struct gkyl_array *integrated_moms_rhs;
};

// forward declare species struct
struct gk_species;

struct gk_lbo_collisions {  
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

  struct gk_species_moment moms; // moments needed in LBO (single array includes Zeroth, First, and Second moment)

  struct gkyl_array *m0;
  struct gkyl_array *vtsq;
  struct gkyl_array *m2self; // m2self used for robustness of LBO
  struct gkyl_array *self_mnu_m0[GKYL_MAX_SPECIES], *self_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *other_mnu_m0[GKYL_MAX_SPECIES], *other_mnu[GKYL_MAX_SPECIES];
  struct gkyl_array *greene_num[GKYL_MAX_SPECIES], *greene_den[GKYL_MAX_SPECIES];
  gkyl_dg_bin_op_mem *greene_factor_mem; // memory needed in computing Greene factor
  struct gkyl_array *greene_factor[GKYL_MAX_SPECIES];

  int num_cross_collisions; // number of species we cross-collide with
  struct gk_species *collide_with[GKYL_MAX_SPECIES]; // pointers to cross-species we collide with

  gkyl_prim_lbo_calc *coll_pcalc; // LBO primitive moment calculator
  gkyl_prim_lbo_cross_calc *cross_calc; // LBO cross-primitive moment calculator
  gkyl_dg_updater_collisions *coll_slvr; // collision solver
};

struct gk_bgk_collisions {  
  struct gkyl_array *nu_sum; // BGK collision frequency 
  struct gkyl_array *nu_sum_host; // BGK collision frequency host-side for I/O
  struct gkyl_array *self_nu; // BGK self-collision frequency

  bool normNu; // Boolean to determine if using Spitzer value
  struct gkyl_array *norm_nu; // Array for normalization factor computed from Spitzer updater n/sqrt(2 vt^2)^3
  struct gkyl_array *nu_init; // Array for initial collisionality when using Spitzer updater
  struct gkyl_spitzer_coll_freq* spitzer_calc; // Updater for Spitzer collisionality if computing Spitzer value

  struct gk_species_moment moms; // moments needed in BGK (single array includes Zeroth, First, and Second moment)

  struct gkyl_array *fmax;
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
  struct gkyl_mom_cross_bgk_gyrokinetic *cross_bgk; // cross-species moment computation

  struct gkyl_correct_maxwellian_gyrokinetic *corr_max; // Maxwellian correction
  struct gkyl_proj_maxwellian_on_basis *proj_max; // Maxwellian projection object
  struct gkyl_bgk_collisions *up_bgk; // BGK updater (also computes stable timestep)
};

struct gk_boundary_fluxes {
  struct gk_species_moment gammai[2*GKYL_MAX_CDIM]; // integrated moments
  gkyl_ghost_surf_calc *flux_slvr; // boundary flux solver
};

struct gk_react {
  int num_react; // number of reactions
  bool all_gk; // boolean for if reactions are only between gyrokinetic species
  struct gkyl_gyrokinetic_react_type react_type[GKYL_MAX_SPECIES]; // input struct for type of reactions

  struct gkyl_array *f_react; // distribution function array which holds update for each reaction
                              // form depend on react->type_self, e.g., for ionization and react->type_self == GKYL_SELF_ELC
                              // f_react = n_elc*coeff_react*(2*fmax(n_elc, upar_donor, vtiz^2) - f_elc)

  struct gkyl_proj_maxwellian_on_basis *proj_max; // Maxwellian projection object

  enum gkyl_react_id react_id[GKYL_MAX_SPECIES]; // what type of reaction (ionization, charge exchange, recombination)
  enum gkyl_react_self_type type_self[GKYL_MAX_SPECIES]; // what is the role of species in this reaction
  struct gk_species *species_elc[GKYL_MAX_SPECIES]; // pointers to electron species being reacted with
  struct gk_species *species_ion[GKYL_MAX_SPECIES]; // pointers to ion species being reacted with
  int elc_idx[GKYL_MAX_SPECIES]; // integer index of electron species being reacted with 
  int ion_idx[GKYL_MAX_SPECIES]; // integer index of ion species being reacted with 
  int donor_idx[GKYL_MAX_SPECIES]; // integer index of donor species being reacted with 

  struct gk_species_moment moms_elc[GKYL_MAX_SPECIES]; // for computing moments of electron species in reaction
  struct gk_species_moment moms_ion[GKYL_MAX_SPECIES]; // for computing moments of ion species in reaction
  struct gk_species_moment moms_donor[GKYL_MAX_SPECIES]; // for computing moments of donor species in reaction

  struct gkyl_array *coeff_react[GKYL_MAX_SPECIES]; // reaction rate
  struct gkyl_array *coeff_react_host[GKYL_MAX_SPECIES]; // reaction rate
  struct gkyl_array *vt_sq_iz[GKYL_MAX_SPECIES]; // ionization temperature
  struct gkyl_array *m0_elc[GKYL_MAX_SPECIES]; // electron density
  struct gkyl_array *m0_ion[GKYL_MAX_SPECIES]; // ion density
  struct gkyl_array *m0_donor[GKYL_MAX_SPECIES]; // donor density
  struct gkyl_array *m0_mod[GKYL_MAX_SPECIES]; // to rescale fmax to have correct density
  struct gkyl_array *prim_vars[GKYL_MAX_SPECIES]; // primitive variables of donor (gk) or ion (vlasov), used for fmax
  union {
    // ionization
    struct {
      struct gkyl_dg_iz *iz[GKYL_MAX_SPECIES];
    };
    // recombination
    struct {
      struct gkyl_dg_recomb *recomb[GKYL_MAX_SPECIES];
    };
  };  
};

struct gk_proj {
  enum gkyl_projection_id proj_id; // type of projection
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
      struct gkyl_array *upar; // host-side upar
      struct gkyl_array *udrift; // host-side udrift
      struct gkyl_array *prim_moms; // host-side prim_moms 

      struct gkyl_array *dens_mod; // array for correcting density

      struct gkyl_array *prim_moms_dev; // device-side prim_moms for GPU simulations
      struct gkyl_array *dens_dev; // device-side density for GPU simulations
      struct gkyl_dg_bin_op_mem *mem; // memory needed in correcting density

      struct gkyl_proj_on_basis *proj_dens; // projection operator for density
      struct gkyl_proj_on_basis *proj_upar; // projection operator for upar
      struct gkyl_proj_on_basis *proj_udrift; // projection operator for upar
      
      union {
        // Maxwellian-specific arrays and functions
        struct {
          struct gkyl_array *vtsq; // host-side vth^2 = T/m (temperature/mass)
          struct gkyl_proj_on_basis *proj_temp; // projection operator for temperature
          struct gkyl_proj_maxwellian_on_basis *proj_max_prim; // Maxwellian projection object
        };
        // Bi-Maxwellian-specific arrays and functions
        struct {
          struct gkyl_array *vtsqpar; // host-side vth_par^2 = Tpar/m (parallel temperature/mass)
          struct gkyl_array *vtsqperp; // host-side vth_perp^2 = Tperp/m (perpendicular temperature/mass)
          struct gkyl_proj_on_basis *proj_temppar; // projection operator for parallel temperature
          struct gkyl_proj_on_basis *proj_tempperp; // projection operator for parallel temperature
          struct gkyl_proj_bimaxwellian_on_basis *proj_bimax; // Bi-Maxwellian projection object
        };
      };
    };
    // Maxwellian from lab moments, includes correction to Maxwellian to produce desired moments
    struct { 
      struct gkyl_array *lab_moms; // lab moms (M0, M1, M2)
      struct gkyl_array *lab_moms_host; // host-side lab moms (M0, M1, M2) for GPU simulations

      struct gkyl_proj_on_basis *proj_lab_moms; // projection operator for (M0, M1, M2)

      struct gkyl_correct_maxwellian_gyrokinetic *corr_max_lab; // Maxwellian correction
      struct gkyl_proj_maxwellian_on_basis *proj_max_lab; // Maxwellian projection object      
    };
  };
};

struct gk_source {
  enum gkyl_source_id source_id; // type of source
  bool write_source; // optional parameter to write out source distribution
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
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    
  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

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
  
  struct gkyl_array *phi; // array for electrostatic potential
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
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  gkyl_dg_updater_gyrokinetic *slvr; // Gyrokinetic solver 
  struct gkyl_dg_eqn *eqn_gyrokinetic; // Gyrokinetic equation object
  
  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions
  bool bc_is_np[3]; // whether BC is nonperiodic.

  // boundary conditions on lower/upper edges in each direction  
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
  // GK_IWL sims need SOL ghost and skin ranges.
  struct gkyl_range lower_skin_par_sol, lower_ghost_par_sol;
  struct gkyl_range upper_skin_par_sol, upper_ghost_par_sol;

  struct gk_proj proj_init; // projector for initial conditions

  enum gkyl_source_id source_id; // type of source
  struct gk_source src; // applied source

  // boundary fluxes
  struct gk_boundary_fluxes bflux;

  // collisions
  struct {
    enum gkyl_collision_id collision_id; // type of collisions
    union {
      struct gk_lbo_collisions lbo; // LBO collisions object
      struct gk_bgk_collisions bgk; // BGK collisions object
    };      
  };

  bool has_reactions; 
  bool has_neutral_reactions; 
  struct gk_react react; // reaction object for reactions with other plasma species
  struct gk_react react_neut; // reaction object for reactions with neutral species

  enum gkyl_radiation_id radiation_id; // type of radiation
  struct gk_rad_drag rad; // radiation object

  // gyrokinetic diffusion
  bool has_diffusion; // flag to indicate there is applied diffusion
  struct gkyl_array *diffD; // array for diffusion tensor
  struct gkyl_dg_updater_diffusion_gyrokinetic *diff_slvr; // gyrokinetic diffusion equation solver

  double *omega_cfl;
};

// neutral species data
struct gk_neut_species {
  struct gkyl_gyrokinetic_neut_species info; // data for neutral species
  
  struct gkyl_job_pool *job_pool; // Job pool
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  struct gkyl_range global, global_ext; // global, global-ext conf-space ranges    
  struct app_skin_ghost_ranges skin_ghost; // conf-space skin/ghost

  struct gkyl_comm *comm;   // communicator object for phase-space arrays
  int nghost[GKYL_MAX_DIM]; // number of ghost-cells in each direction

  struct gkyl_rect_grid grid_vel; // velocity space grid
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext velocity-space ranges

  struct gkyl_array *f, *f1, *fnew; // arrays for updates
  struct gkyl_array *cflrate; // CFL rate in each cell
  struct gkyl_array *bc_buffer; // buffer for BCs (used by bc_basic)
  struct gkyl_array *bc_buffer_lo_fixed, *bc_buffer_up_fixed; // fixed buffers for time independent BCs 

  struct gkyl_array *f_host; // host copy for use IO and initialization

  enum gkyl_field_id field_id; // type of field equation (always GKYL_FIELD_NULL)
  enum gkyl_model_id model_id; // type of Vlasov equation (always GKYL_MODEL_GEN_GEO)

  struct gkyl_array *alpha_surf; // array for surface phase space flux (v^i = v . e^i)
  struct gkyl_array *sgn_alpha_surf; // array for the sign of the surface phase space flux at quadrature points
                                     // utilized for numerical flux function
                                     // F = alpha_surf/2 ( (f^+ + f^-) - sign_alpha_surf*(f^+ - f^-) )
  struct gkyl_array *const_sgn_alpha; // boolean array for if the surface phase space flux is single signed
                                      // if true, numerical flux function inside kernels simplifies to
                                      // F = alpha_surf*f^- (if sign_alpha_surf = 1), 
                                      // F = alpha_surf*f^+ (if sign_alpha_surf = -1)
  struct gkyl_array *cot_vec; // array for cotangent vectors

  struct gk_species_moment m0; // for computing density
  struct gk_species_moment integ_moms; // integrated moments
  struct gk_species_moment *moms; // diagnostic moments
  double *red_integ_diag, *red_integ_diag_global; // for reduction of integrated moments
  gkyl_dynvec integ_diag; // integrated moments reduced across grid
  bool is_first_integ_write_call; // flag for integrated moments dynvec written first time

  gkyl_dg_updater_vlasov *slvr; // Vlasov solver 
  struct gkyl_dg_eqn *eqn_vlasov; // Vlasov equation object
  
  // boundary conditions on lower/upper edges in each direction  
  enum gkyl_species_bc_type lower_bc[3], upper_bc[3];
  // Pointers to updaters that apply BC.
  struct gkyl_bc_basic *bc_lo[3];
  struct gkyl_bc_basic *bc_up[3];
  // To simplify BC application, store local skin and ghost ranges
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  struct gk_proj proj_init; // projector for initial conditions

  enum gkyl_source_id source_id; // type of source
  struct gk_source src; // applied source

  bool has_neutral_reactions;
  struct gk_react react_neut; // reaction object

  double *omega_cfl_ptr;
};

// field data
struct gk_field {
  struct gkyl_gyrokinetic_field info; // data for field

  enum gkyl_gkfield_id gkfield_id;

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

  struct gkyl_array *weight;
  double es_energy_fac_1d; 
  struct gkyl_array *es_energy_fac; 
  struct gkyl_array *epsilon; 
  struct gkyl_array *kSq; 

  struct gkyl_fem_parproj *fem_parproj; // FEM smoother for projecting DG functions onto continuous FEM basis
                                        // weight*phi_{fem} = phi_{dg} 
  struct gkyl_fem_parproj *fem_parproj_sol;
  struct gkyl_fem_parproj *fem_parproj_core;

  struct gkyl_deflated_fem_poisson *deflated_fem_poisson; // poisson solver which solves on lines in x or planes in xy
                                                          // - nabla . (epsilon * nabla phi) - kSq * phi = rho

  struct gkyl_array_integrate *calc_em_energy;
  double *em_energy_red, *em_energy_red_global; // memory for use in GPU reduction of EM energy
  gkyl_dynvec integ_energy; // integrated energy components

  bool is_first_energy_write_call; // flag for energy dynvec written first time

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
};

// gyrokinetic object: used as opaque pointer in user code
struct gkyl_gyrokinetic_app {
  char name[128]; // name of app
  struct gkyl_job_pool *job_pool; // Job pool
  
  int cdim, vdim; // conf, velocity space dimensions
  int poly_order; // polynomial order
  double tcurr; // current time
  double cfl; // CFL number

  bool use_gpu; // should we use GPU (if present)

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

  struct gkyl_basis basis, neut_basis; // phase-space and phase-space basis for neutrals
  struct gkyl_basis confBasis; // conf-space basis
  
  struct gkyl_comm *comm;   // communicator object for conf-space arrays

  // pointers to basis on device (these point to host structs if not
  // on GPU)
  struct {
    struct gkyl_basis *basis, *neut_basis, *confBasis;
  } basis_on_dev;

  struct gk_geometry *gk_geom;

  bool update_field; // are we updating the field?
  struct gk_field *field; // pointer to field object

  // species data
  int num_species;
  struct gk_species *species; // data for each species
  
  // neutral species data
  int num_neut_species;
  struct gk_neut_species *neut_species; // data for each species

  struct gkyl_gyrokinetic_stat stat; // statistics
};

/** gkyl_gyrokinetic_app private API */

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
 * @param nm Name string indicating moment type
 */
void gk_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
  struct gk_species_moment *sm, const char *nm);

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
 * Release species LBO object.
 *
 * @param app gyrokinetic app object
 * @param lbo Species LBO object to release
 */
void gk_species_lbo_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lbo_collisions *lbo);

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
  const struct gk_species *species,
  struct gk_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

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
  const struct gk_species *species,
  struct gk_react *react,
  const struct gkyl_array *f_self, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

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
  const struct gk_species *s, struct gk_react *react,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

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
 * @param app Gyrokinetic app object
 * @param s Species object 
 * @param bflux Species boundary flux object
 */
void gk_species_bflux_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s,
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
void gk_species_bflux_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_boundary_fluxes *bflux, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species boundary flux object.
 *
 * @param app Gyrokinetic app object
 * @param bflux Species boundary flux object to release
 */
void gk_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const struct gk_boundary_fluxes *bflux);

/** gk_species_projection API */

/**
 * Initialize species projection object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param inp Input struct for projection (contains functions pointers for type of projection)
 * @param proj Species projection object
 */
void gk_species_projection_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_projection inp, struct gk_proj *proj);

/**
 * Compute species projection
 *
 * @param app gyrokinetic app object
 * @param species Species object
 * @param proj Species projection object
 * @param f Output distribution function from projection
 * @param tm Time for use in projection
 */
void gk_species_projection_calc(gkyl_gyrokinetic_app *app, const struct gk_species *species, 
  struct gk_proj *proj, struct gkyl_array *f, double tm);

/**
 * Release species projection object.
 *
 * @param app gyrokinetic app object
 * @param proj Species projection object to release
 */
void gk_species_projection_release(const struct gkyl_gyrokinetic_app *app, const struct gk_proj *proj);

/** gk_species_source API */

/**
 * Initialize species source object.
 *
 * @param app gyrokinetic app object
 * @param s Species object 
 * @param src Species source object
 */
void gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_source *src);

/**
 * Compute species applied source term
 *
 * @param app gyrokinetic app object
 * @param species Species object
 * @param src Species source object
 * @param tm Time for use in source
 */
void gk_species_source_calc(gkyl_gyrokinetic_app *app, const struct gk_species *species, 
  struct gk_source *src, double tm);

/**
 * Compute RHS contribution from source
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param src Pointer to source
 * @param fin Input distribution function
 * @param rhs On output, the distribution function
 */
void gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Release species source object.
 *
 * @param app gyrokinetic app object
 * @param src Species source object to release
 */
void gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_species API */

/**
 * Initialize species.
 *
 * @param gk Input gk data
 * @param app gyrokinetic app object
 * @param s On output, initialized species object
 */
void gk_species_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_species *s);

/**
 * Compute species initial conditions.
 *
 * @param app gyrokinetic app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute the part of the species initial conditions that depends on other
 * species.
 *
 * @param app gyrokinetic app object
 * @param species Species object
 * @param t0 Time for use in ICs
 */
void gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *species, double t0);

/**
 * Compute RHS from species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param fin Input distribution function
 * @param rhs On output, the RHS from the species object
 * @return Maximum stable time-step
 */
double gk_species_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

/**
 * Apply BCs to species distribution function.
 *
 * @param app gyrokinetic app object
 * @param species Pointer to species
 * @param f Field to apply BCs
 */
void gk_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with collision timers.
 *
 * @param app App object to update stat timers
 */
void gk_species_coll_tm(gkyl_gyrokinetic_app *app);

/**
 * Fill stat object in app with collisionless timers.
 *
 * @param app App object to update stat timers
 */
void gk_species_tm(gkyl_gyrokinetic_app *app);

/**
 * Delete resources used in species.
 *
 * @param app gyrokinetic app object
 * @param species Species object to delete
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
 * @param f_self Input self distribution function
 * @param fin Input distribution functions (size: num_species)
 * @param fin_neut Input neutral distribution functions (size: num_neut_species)
 */
void gk_neut_species_react_cross_moms(gkyl_gyrokinetic_app *app,
  const struct gk_neut_species *species,
  struct gk_react *react,
  const struct gkyl_array *f_self, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[]);

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
  const struct gk_neut_species *s, struct gk_react *react,
  const struct gkyl_array *fin, struct gkyl_array *rhs);

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
 * Release Neutral species source object.
 *
 * @param app gyrokinetic app object
 * @param src Neutral species source object to release
 */
void gk_neut_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src);

/** gk_neut_species API */

/**
 * Initialize Neutral species.
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
 * Apply BCs to neutral species distribution function
 *
 * @param app gyrokinetic app object
 * @param species Pointer to neutral species
 * @param f Field to apply BCs
 */
void gk_neut_species_apply_bc(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f);

/**
 * Fill stat object in app with collisionless timers for neutral species.
 *
 * @param app App object to update stat timers
 */
void gk_neut_species_tm(gkyl_gyrokinetic_app *app);

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
 * Compute field energy diagnostic
 *
 * @param app gyrokinetic app object
 * @param tm Time at which diagnostic is computed
 * @param field Pointer to field
 */
void gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field);

/**
 * Release resources allocated by field
 *
 * @param app gyrokinetic app object
 * @param f Field object to release
 */
void gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f);
