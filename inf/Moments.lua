-- Thin wrapper around the gkyl_moment app. This tries to preserve
-- compatibility with G2 input files as much as possible. To use this,
-- install G0 (make install) and set LUA_PATH to:
--
-- export LUA_PATH="$HOME/gkylsoft/gkylzero/lib/?.lua;;"
--
-- (Please note the two semi-colons at the end!)
--
-- Then make an input file and run that through luajit. You can
-- install luajit by hand from sources, or using apt get on Ubuntu.
--

local ffi = require "ffi"

-- load shared library
local install_prefix = os.getenv("HOME") .. "/gkylsoft/gkylzero"
local C = ffi.load(install_prefix .. "/lib/libgkylzero.so")

-- set global package paths
package.path = package.path .. ";" .. install_prefix .. "/lib/?.lua"

-- pick default if val is nil, else pick val
local function pickBool(val, default)
   if val == nil then return default else return val end
end

-- declare some top-level things we want to expose
ffi.cdef [[

// This needs to be enum to allow usage below
enum { GKYL_MAX_SPECIES = 8, GKYL_MAX_DIM = 7 };

/**
 * Set the global flag to turn on memory allocation/deallocation
 * tracing.
 *
 * @param flag Flag to set
 */
void gkyl_mem_debug_set(bool flag);

/**
 * Set the global flag to turn on cuda memory allocation/deallocation
 * tracing.
 *
 * @param flag Flag to set
 */
void gkyl_cu_dev_mem_debug_set(bool flag);

]]

-- gkyl_eqn_type.h
ffi.cdef [[

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EQN_EULER, // Euler equations
  GKYL_EQN_SR_EULER, // SR Euler equations
  GKYL_EQN_ISO_EULER, // Isothermal Euler equations
  GKYL_EQN_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_EQN_MAXWELL, // Maxwell equations
  GKYL_EQN_MHD,  // Ideal MHD equations
  GKYL_EQN_BURGERS, // Burgers equations
  GKYL_EQN_ADVECTION, // Scalar advection equation
  GKYL_EQN_EULER_PKPM, // Euler equations with parallel-kinetic-perpendicular-moment (pkpm) model
};

// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_SR_E_B, // Maxwell (E, B) with special relativity
  GKYL_FIELD_PHI, // Poisson (only phi)  
  GKYL_FIELD_PHI_A, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL, // no field is present
  GKYL_FIELD_SR_NULL // no field is present, special relativistic Vlasov
};

// Identifiers for specific collision object types
enum gkyl_collision_id {
  GKYL_NO_COLLISIONS = 0, // No collisions. This is default
  GKYL_BGK_COLLISIONS, // BGK Collision operator
  GKYL_LBO_COLLISIONS // LBO Collision operator
};

// Identifiers for specific source object types
enum gkyl_source_id {
  GKYL_NO_SOURCE = 0, // No source. This is default
  GKYL_FUNC_SOURCE, // Source given by function
  GKYL_BFLUX_SOURCE // Source which scales to boundary fluxes
};

// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD // Gauss-Lobatto quadrature
};
]]

-- gkyl_app.h
ffi.cdef [[

// Update status
struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

// Boundary conditions on particles
enum gkyl_species_bc_type {
  GKYL_SPECIES_COPY = 0, // copy BCs
  GKYL_SPECIES_REFLECT, // perfect reflector
  GKYL_SPECIES_ABSORB, // Absorbing BCs
  GKYL_SPECIES_NO_SLIP, // no-slip boundary conditions
  GKYL_SPECIES_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_SPECIES_FUNC, // Function boundary conditions
  GKYL_SPECIES_FIXED_FUNC, // Fixed function, time-independent, boundary conditions
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // perfect electrical conductor (PEC) BCs
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_FIELD_FUNC, // Function boundary conditions
};

]]


-- gkyl_util.h
ffi.cdef [[
/**
 * Time-trigger. Typical initialization is:
 * 
 * struct gkyl_tm_trigger tmt = { .dt = tend/nframe };
 */
struct gkyl_tm_trigger {
  int curr; // current counter
  double dt, tcurr; // Time-interval, current time
};

/**
 * Check if the tcurr should trigger and bump internal counters if it
 * does. This only works if sequential calls to this method have the
 * tcurr monotonically increasing.
 *
 * @param tmt Time trigger object
 * @param tcurr Current time.
 * @return 1 if triggered, 0 otherwise
 */
int gkyl_tm_trigger_check_and_bump(struct gkyl_tm_trigger *tmt, double tcurr);

/**
 * Compute time in seconds since epoch.
 *
 * @return Time in seconds
 */
double gkyl_time_now(void);
]]

-- gkyl_wave_prop.h
ffi.cdef [[

// Limiters
enum gkyl_wave_limiter {
  GKYL_NO_LIMITER = 1, // to allow default to be 0
  GKYL_MIN_MOD,
  GKYL_SUPERBEE,
  GKYL_VAN_LEER,
  GKYL_MONOTONIZED_CENTERED,
  GKYL_BEAM_WARMING,
  GKYL_ZERO
};
]]

-- gkyl_mp_scheme.h
ffi.cdef [[
// Base reconstruction scheme to use
enum gkyl_mp_recon {
  GKYL_MP_U5 = 0, // upwind-biased 5th order (default)
  GKYL_MP_C2, // centered second-order
  GKYL_MP_C4, // centered fourth-order  
  GKYL_MP_C6, // centered sixth-order
  GKYL_MP_U1, // upwind-biased 1st order
  GKYL_MP_U3, // upwind-biased 3rd order
};
]]

-- gkyl_wv_eqn.h
ffi.cdef [[
// Flux type for use in wave/qfluct methods
enum gkyl_wv_flux_type { GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_LOW_ORDER_FLUX };

// Forward declare for use in function pointers
typedef struct gkyl_wv_eqn gkyl_wv_eqn;

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn* eqn);

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);
]]

-- Various equation objects
ffi.cdef [[

// Type of Rieman problem solver to use
enum gkyl_wv_euler_rp {
  WV_EULER_RP_ROE = 0, // default
  WV_EULER_RP_HLLC,
  WV_EULER_RP_LAX
};

// input packaged as a struct
struct gkyl_wv_euler_inp {
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_euler_rp rp_type; // type of RP to use
};

/**
 * Create a new Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_new(double gas_gamma);

/**
 * Create a new Euler equation object.
 * 
 * @param inp Input parameters
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_inew(const struct gkyl_wv_euler_inp *inp);

/**
 * Create a new isothermal Euler equation object.
 * 
 * @param vt Thermal velocity
 * @return Pointer to isothermal Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_iso_euler_new(double vt);

/**
 * Create a new maxwell equation object.
 * 
 * @param c_speed Speed of light
 * @param e_fact Factor of light-speed for electric field correction
 * @param b_fact Factor of light-speed for magnetic field correction
 * @return Pointer to Maxwell equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_maxwell_new(double c, double e_fact, double b_fact);

/**
 * Create a new SR Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to SR Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_sr_euler_new(double gas_gamma);

/**
 * Create a new advection equation object.
 * 
 * @param c advection speed
 * @return Pointer to Burgers equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_advect_new(double c);

// Wrappers around various eqn object
struct gkyl_wv_euler { struct gkyl_wv_eqn *eqn; };
struct gkyl_wv_iso_euler { struct gkyl_wv_eqn *eqn; };
struct gkyl_wv_advection { struct gkyl_wv_eqn *eqn; };
struct gkyl_wv_ten_moment { struct gkyl_wv_eqn *eqn; };

/**
 * Create a new Ten moment equation object.
 * 
 * @param k0 Closure parameter
 * @return Pointer to Ten moment equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_ten_moment_new(double k0);
]]

-- gkyl_moment_prim_mhd.h
ffi.cdef [[
/**
 * Compute conserved variables from primitive variables.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param pv Primitive variables
 * @param q Conserved variables

 */
void gkyl_mhd_cons_vars(double gas_gamma, const double pv[8], double q[8]);
]]

-- gkyl_wv_mhd.h
ffi.cdef [[

// Type of Rieman problem solver to use
enum gkyl_wv_mhd_rp {
  WV_MHD_RP_ROE = 0, // default
  WV_MHD_RP_HLLD,
  WV_MHD_RP_LAX
};

// flags to indicate which divergence constraint scheme to use
enum gkyl_wv_mhd_div_constraint {
  GKYL_MHD_DIVB_NONE,
  GKYL_MHD_DIVB_EIGHT_WAVES,
  GKYL_MHD_DIVB_GLM
};

struct gkyl_wv_mhd_inp {
  double gas_gamma; // gas adiabatic constant
  enum gkyl_wv_mhd_rp rp_type; // Riemann problem solver
  enum gkyl_wv_mhd_div_constraint divergence_constraint; // divB correction
  double glm_ch; // factor to use in GLM scheme
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
};

/**
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param divb Divergence constraint method
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(const struct gkyl_wv_mhd_inp *inp);

]]

-- gkyl_moment.h
ffi.cdef [[

// Parameters for moment species
struct gkyl_moment_species {
  char name[128]; // species name
  double charge, mass; // charge and mass
  bool has_grad_closure; // has gradient-based closure (only for 10 moment)
  enum gkyl_wave_limiter limiter; // limiter to use
  struct gkyl_wv_eqn *equation; // equation object

  int evolve; // evolve species? 1-yes, 0-no
  bool force_low_order_flux; // should  we force low-order flux?

  void *ctx; // context for initial condition init function (and potentially other functions)
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
  // pointer to applied acceleration/forces function
  void (*app_accel_func)(double t, const double *xn, double *fout, void *ctx);
  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
  // pointer to boundary condition functions along x
  void (*bcx_lower_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  void (*bcx_upper_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  // pointer to boundary condition functions along y
  void (*bcy_lower_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  void (*bcy_upper_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  // pointer to boundary condition functions along z
  void (*bcz_lower_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  void (*bcz_upper_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
};

// Parameter for EM field
struct gkyl_moment_field {
  double epsilon0, mu0;
  double elc_error_speed_fact, mag_error_speed_fact;

  enum gkyl_wave_limiter limiter; // limiter to use

  int evolve; // evolve field? 1-yes, 0-no

  void *ctx; // context for initial condition init function (and potentially other functions)
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);
  // pointer to applied current function
  void (*app_current_func)(double t, const double *xn, double *fout, void *ctx);

  bool is_ext_em_static; // flag to indicate if external field is time-independent
  // pointer to external fields
  void (*ext_em_func)(double t, const double *xn, double *fout, void *ctx);
  
  // boundary conditions
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];
  // pointer to boundary condition functions along x
  void (*bcx_lower_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  void (*bcx_upper_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  // pointer to boundary condition functions along y
  void (*bcy_lower_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  void (*bcy_upper_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  // pointer to boundary condition functions along z
  void (*bcz_lower_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
  void (*bcz_upper_func)(double t, int nc, const double *skin, double *ghost, const double *skin_xc, const double *ghost_xc, void *ctx);
};

// Choices of schemes to use in the fluid solver 
enum gkyl_moment_scheme {
  GKYL_MOMENT_WAVE_PROP = 0, // default, 2nd-order FV
  GKYL_MOMENT_MP, // monotonicity-preserving Suresh-Huynh scheme
  GKYL_MOMENT_KEP // Kinetic-energy preserving scheme
};

/**
 * Range object, representing an N-dimensional integer index
 * set. Lower and upper limits are inclusive.
 */
struct gkyl_range {
  int ndim; // number of dimension
  int lower[GKYL_MAX_DIM]; // lower bound
  int upper[GKYL_MAX_DIM]; // upper bound (inclusive)
  long volume; // total volume of range
    
  // do not access directly
  uint32_t flags; // Flags for internal use
  int ilo[GKYL_MAX_DIM]; // for use in inverse indexer
  long ac[GKYL_MAX_DIM+1]; // coefficients for indexing
  long iac[GKYL_MAX_DIM+1]; // for use in sub-range inverse indexer
  long linIdxZero; // linear index of {0,0,...}
  int nsplit, tid; // number of splits, split ID

  // FOR CUDA ONLY
  int nthreads, nblocks; // CUDA kernel launch specifiers for range-based ops
};

// Forward declaration
struct gkyl_comm;

// Lower-level inputs: in general this does not need to be set by the
// user. It is needed when the App is being created on a sub-range of
// the global range, and is meant for use in higher-level drivers that
// use MPI or other parallel mechanism.
struct gkyl_moment_low_inp {
  // local range over which App operates
  struct gkyl_range local_range;
  // communicator to used
  struct gkyl_comm *comm;
};

// Top-level app parameters
struct gkyl_moment {
  char name[128]; // name of app: used as output prefix

  int ndim; // space dimensions
  double lower[3], upper[3]; // lower, upper bounds
  int cells[3]; // config-space cells

  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  double cfl_frac; // CFL fraction to use

  enum gkyl_moment_scheme scheme_type; // scheme to update fluid and moment eqns
  enum gkyl_mp_recon mp_recon; // reconstruction scheme to use
  bool skip_mp_limiter; // should MP limiter be skipped?
  bool use_hybrid_flux_kep; // should shock-hybrid scheme be used when using KEP?

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_skip_dirs; // number of directions to skip
  int skip_dirs[3]; // directions to skip

  int num_species; // number of species
  struct gkyl_moment_species species[GKYL_MAX_SPECIES]; // species objects
  struct gkyl_moment_field field; // field object

  bool has_collision; // has collisions
  // scaling factors for collision frequencies so that nu_sr=nu_base_sr/rho_s
  // nu_rs=nu_base_rs/rho_r, and nu_base_sr=nu_base_rs
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  // this should not be set by typical user-facing code but only by
  // higher-level drivers
  bool has_low_inp; // should one use low-level inputs?
  struct gkyl_moment_low_inp low_inp; // low-level inputs
};

// Simulation statistics
struct gkyl_moment_stat {
  long nup; // calls to update
  double total_tm; // time for simulation (not including ICs)
  
  long nfail; // number of failed time-steps

  //// wave_prop stuff
  double species_tm; // time to compute species updates
  double field_tm; // time to compute field updates
  double sources_tm; // time to compute source terms

  //// stuff for MP-XX/SSP-RK schemes
  long nfeuler; // calls to forward-Euler method
    
  long nstage_2_fail; // number of failed RK stage-2s
  long nstage_3_fail; // number of failed RK stage-3s

  double stage_2_dt_diff[2]; // [min,max] rel-diff for stage-2 failure
  double stage_3_dt_diff[2]; // [min,max] rel-diff for stage-3 failure
  
  double init_species_tm; // time to initialize all species
  double init_field_tm; // time to initialize fields

  double species_rhs_tm; // time to compute species collisionless RHS
  double species_bc_tm; // time to apply BCs

  double field_rhs_tm; // time to compute field RHS
  double field_bc_tm; // time to apply BCs
};

// Object representing moments app
typedef struct gkyl_moment_app gkyl_moment_app;

/**
 * Construct a new moments app.
 *
 * @param vm App inputs. See struct docs.
 * @return New moment app object.
 */
gkyl_moment_app* gkyl_moment_app_new(struct gkyl_moment *mom);

/**
 * Compute maximum estimated stable dt wtih current app state. Call
 * after app initialized and after initial conditions set.
 *
 * @param app App object.
 * @retuen maximum estimated stable dt
 */
double gkyl_moment_app_max_dt(gkyl_moment_app* app);

/**
 * Initialize species and field.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0);

/**
 * Initialize field.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0);

/**
 * Initialize species.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write(const gkyl_moment_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write_field(const gkyl_moment_app *app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to write
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame);

/**
 * Write field energy to file.
 *
 * @param app App object.
 */
void gkyl_moment_app_write_field_energy(gkyl_moment_app *app);

/**
 * Write integrated moments to file.
 *
 * @param app App object.
 */
void gkyl_moment_app_write_integrated_mom(gkyl_moment_app *app);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_moment_app_stat_write(const gkyl_moment_app* app);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status gkyl_moment_update(gkyl_moment_app *app, double dt);

/**
 * Calculate integrated field energy
 *
 * @param tm Time at which integrated diagnostic are to be computed
 * @param app App object.
 */
void gkyl_moment_app_calc_field_energy(gkyl_moment_app *app, double tm);

/**
 * Calculate integrated moments
 *
 * @param app App object.
 * @param tm Time at which integrated diagnostic are to be computed
 */
void gkyl_moment_app_calc_integrated_mom(gkyl_moment_app *app, double tm);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics.
 */
struct gkyl_moment_stat gkyl_moment_app_stat(gkyl_moment_app *app);

/**
 * Free moment app.
 *
 * @param app App to release.
 */
void gkyl_moment_app_release(gkyl_moment_app* app);

]]

-- name to moment_scheme type
local moment_scheme_tags = {
   ["wave_prop"] = C.GKYL_MOMENT_WAVE_PROP,
   ["mp"] = C.GKYL_MOMENT_MP,
   ["kep"] = C.GKYL_MOMENT_KEP
}

-- name to recovery used in MP scheme
local mp_recon_tags = {
   ["u1"] = C.GKYL_MP_U1,
   ["u3"] = C.GKYL_MP_U3,
   ["u5"] = C.GKYL_MP_U5,
   ["c2"] = C.GKYL_MP_C2,
   ["c4"] = C.GKYL_MP_C4,
   ["c6"] = C.GKYL_MP_C6,
}

-- App container
ffi.cdef [[
// Container to store pointer to app and other data

struct gkyl_moment_app_cont {
  double t0, tend; // start and end times
  int nframe; // number of frames to write
  int nspecies; // number of species

  gkyl_moment_app *app; // pointer to app
};
]]

-- module table
local _M = { }

-- methods to turn on/off memory tracing
_M.mem_debug_set = function(flag)
   C.gkyl_mem_debug_set(flag)
end
_M.cu_dev_mem_debug_set = function(flag)
   C.gkyl_cu_dev_mem_debug_set(flag)
end

-- Time in seconds from epoch
_M.time_now = function()
   return C.gkyl_time_now()
end

-- time-trigger object
local tm_trigger_type = ffi.typeof("struct gkyl_tm_trigger")
local tm_trigger_mt = {
   __new = function(self, dt)
      local tmt = ffi.new(self, { curr = 0, dt = dt, tcurr = 0.0 })
      return tmt
   end,
   __index = {
      checkAndBump = function(self, tcurr)
	 return C.gkyl_tm_trigger_check_and_bump(self, tcurr) == 1
      end
   }
}
_M.TimeTrigger = ffi.metatype(tm_trigger_type, tm_trigger_mt)

-- name to RP-type mappings
local euler_rp_tags = {
   ["roe"] = C.WV_EULER_RP_ROE,
   ["hllc"] = C.WV_EULER_RP_HLLC,
   ["lax"] = C.WV_EULER_RP_LAX
}

-- Euler equation
_M.Euler = function(tbl)
   local einp = ffi.new("struct gkyl_wv_euler_inp")
   einp.gas_gamma = tbl.gasGamma
   einp.rp_type = C.WV_EULER_RP_ROE
   if tbl.rp_type then
      einp.rp_type = euler_rp_tags[tbl.rp_type]
   end
   return ffi.gc(C.gkyl_wv_euler_inew(einp), C.gkyl_wv_eqn_release)
end

-- Isothermal Euler equation
_M.IsoEuler = function(tbl)
   return ffi.gc(C.gkyl_wv_iso_euler_new(tbl.vthermal), C.gkyl_wv_eqn_release)
end

-- Advection
_M.Advection = function(tbl)
   return ffi.gc(C.gkyl_wv_advect_new(tbl.speed), C.gkyl_wv_eqn_release)
end

-- Ten-moment equation
_M.TenMoment = function(tbl)
   return ffi.gc(C.gkyl_wv_ten_moment_new(tbl.k0), C.gkyl_wv_eqn_release)
end

-- name to RP-type mappings
local mhd_rp_tags = {
   ["roe"] = C.WV_MHD_RP_ROE,
   ["hlld"] = C.WV_MHD_RP_HLLD,
   ["lax"] = C.WV_MHD_RP_LAX
}

-- name to divB fix mapping
local mhd_div_tags = {
   ["none"] = C.GKYL_MHD_DIVB_NONE,
   ["eight_waves"] = C.GKYL_MHD_DIVB_EIGHT_WAVES,
   ["glm"] = C.GKYL_MHD_DIVB_GLM
}

-- Ideal MHD
_M.MHD = function(tbl)
   local minp = ffi.new("struct gkyl_wv_mhd_inp")
   minp.gas_gamma = tbl.gasGamma
   minp.rp_type = C.WV_MHD_RP_ROE
   if tbl.rp_type then
      minp.rp_type = mhd_rp_tags[tbl.rp_type]
   end
   minp.divergence_constraint = C.GKYL_MHD_DIVB_EIGHT_WAVES
   if tbl.divergenceConstraint then
      minp.divergence_constraint = mhd_div_tags[tbl.divergenceConstraint]
   end
   minp.glm_ch = 1.0
   minp.glm_alpha = 0.4
   return ffi.gc(C.gkyl_wv_mhd_new(minp), C.gkyl_wv_eqn_release)
end

-- Raw wrapper around MHD function
_M.gkyl_mhd_cons_vars = function(gas_gamma, pv, q)
   return C.gkyl_mhd_cons_vars(gas_gamma, pv, q)
end

-- Wraps user given function in a function that can be passed to the C
-- callback APIs
local function gkyl_eval_moment(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end
      local ret = { func(t, xnl) } -- package return into table
      for i=1,#ret do
         fout[i-1] = ret[i]
      end
   end
end

local function gkyl_eval_applied(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end
      local ux,uy,uz = func(t, xnl)

      fout[0] = ux; fout[1] = uy; fout[2] = uz
   end
end

local function gkyl_eval_bc(func)
   return function(t, nc, skin, ghost, skin_xc, ghost_xc, ctx)
      local ret = { func(t, nc, skin, ctx) } -- package return into table
      for i=1,#ret do
         ghost[i-1] = ret[i]
      end
   end
end

local function gkyl_eval_mapc2p(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end
      local ret = { func(t, xnl) } -- package return into table
      for i=1,#ret do
         fout[i-1] = ret[i]
      end
   end
end

-- Species object: table structure is as follows:
--
-- local elc = {
--    charge = -1.0,
--    mass = 1.0,
--    equation = Moments.Euler { gasGamma = 1.4 },
--    init = function(t, xn)
--       -- return initial conditions
--    end

-- }

local limiter_tags = {
   ["no-limiter"] = C.GKYL_NO_LIMITER,
   ["min-mod"] = C.GKYL_MIN_MOD,
   ["superbee"] = C.GKYL_SUPERBEE,
   ["van-leer"] = C.GKYL_VAN_LEER,
   ["monotonized-centered"] = C.GKYL_MONOTONIZED_CENTERED,
   ["beam-warming"] = C.GKYL_BEAM_WARMING,
   ["zero"] = C.GKYL_ZERO,
}

-- this tables stores pointer to the species equation objects so they
-- are not deleted by the GC while the simulation is being constructed
local species_eqn_tbl = { }

local species_type = ffi.typeof("struct gkyl_moment_species")
local species_mt = {
   __new = function(self, tbl)
      local s = ffi.new(self)
      s.charge = tbl.charge
      s.mass = tbl.mass

      s.limiter = limiter_tags["monotonized-centered"]
      if tbl.limiter then
         s.limiter = limiter_tags[tbl.limiter]
      end

      -- we need to insert equation into species_eqn_tbl to prevent it
      -- from getting GC-ed while sim is being constructed.
      table.insert(species_eqn_tbl, tbl.equation)
      s.equation = tbl.equation

      -- initial conditions
      s.ctx = tbl.ctx 
      if tbl.cinit then
	 -- use C function directly, if specified ...
	 s.init = tbl.cinit
      else
	 -- ... or use Lua function
	 s.init = gkyl_eval_moment(tbl.init)
      end
      
      if tbl.app_accel then
         s.app_accel_func = gkyl_eval_applied(tbl.app_accel)
      end

      s.force_low_order_flux = false
      if tbl.force_low_order_flux then
	 s.force_low_order_flux = tbl.force_low_order_flux
      end

      -- boundary conditions
      if tbl.bcx then
         s.bcx[0], s.bcx[1] = tbl.bcx[1], tbl.bcx[2]
      end
      if tbl.bcy then
         s.bcy[0], s.bcy[1] = tbl.bcy[1], tbl.bcy[2]
      end
      if tbl.bcz then
         s.bcz[0], s.bcz[1] = tbl.bcz[1], tbl.bcz[2]
      end

      if tbl.bcx_lower_func then
         s.bcx_lower_func = gkyl_eval_bc(tbl.bcx_lower_func)
      end
      if tbl.bcy_lower_func then
         s.bcy_lower_func = gkyl_eval_bc(tbl.bcy_lower_func)
      end
      if tbl.bcz_lower_func then
         s.bcz_lower_func = gkyl_eval_bc(tbl.bcz_lower_func)
      end

      if tbl.bcx_upper_func then
         s.bcx_upper_func = gkyl_eval_bc(tbl.bcx_upper_func)
      end
      if tbl.bcy_upper_func then
         s.bcy_upper_func = gkyl_eval_bc(tbl.bcy_upper_func)
      end
      if tbl.bcz_upper_func then
         s.bcz_upper_func = gkyl_eval_bc(tbl.bcz_upper_func)
      end

      return s
   end,
   __index = {
      -- we need this here also to be consistent with G2 App
      bcReflect = C.GKYL_SPECIES_REFLECT,
      bcWall = C.GKYL_SPECIES_REFLECT,
      bcCopy = C.GKYL_SPECIES_COPY,
      bcNoSlip = C.GKYL_SPECIES_NO_SLIP,
      bcWedge = C.GKYL_SPECIES_WEDGE,
      bcFunc = C.GKYL_SPECIES_FUNC
   }
}
_M.Species = ffi.metatype(species_type, species_mt)

-- Wraps user given field initialization function in a function that
-- can be passed to the C callback APIs
local function gkyl_eval_field(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end

      local ex,ey,ez,bx,by,bz = func(t, xnl)

      fout[0] = ex; fout[1] = ey; fout[2] = ez
      fout[3] = bx; fout[4] = by; fout[5] = bz
      fout[6] = 0.0; fout[7] = 0.0
   end
end

local function gkyl_eval_ext_field(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end

      local ex,ey,ez,bx,by,bz = func(t, xnl)

      fout[0] = ex; fout[1] = ey; fout[2] = ez
      fout[3] = bx; fout[4] = by; fout[5] = bz
   end
end

-- Field object: table structure is as follows:
--
-- local field = {
--    epsilon0 = 1.0,  mu0 = 1.0,
--    elcErrorSpeedFactor = 0.0, mgnErrorSpeedFactor = 0.0,
--    init = function(t, xn)
--       -- return EM field: 6 components Ex,Ey,Ez,Bx,By,Bz
--    end
--    evolve = true,
-- }

local field_type = ffi.typeof("struct gkyl_moment_field")
local field_mt = {
   __new = function(self, tbl)
      local f = ffi.new(self)
      f.epsilon0 = tbl.epsilon0
      f.mu0 = tbl.mu0

      f.elc_error_speed_fact = 0.0
      if (tbl.elcErrorSpeedFactor) then
         f.elc_error_speed_fact = tbl.elcErrorSpeedFactor
      end

      f.mag_error_speed_fact = 1.0
      if (tbl.mgnErrorSpeedFactor) then
         f.mag_error_speed_fact = tbl.mgnErrorSpeedFactor
      end

      f.ctx = tbl.ctx -- no need for context
      if tbl.cinit then
	 -- use C function directly, if specified ...
	 f.init = tbl.cinit
      else
	 -- ... or use Lua function
	 f.init = gkyl_eval_field(tbl.init)
      end
      
      if tbl.app_current then
         f.app_current_func = gkyl_eval_applied(tbl.app_current)
      end

      f.is_ext_em_static = false
      if tbl.is_ext_em_static then
	 f.is_ext_em_static = tbl.is_ext_em_static
      end
      if tbl.ext_em_func then
         f.ext_em_func = gkyl_eval_ext_field(tbl.ext_em_func)
      end

      -- boundary conditions
      if tbl.bcx then
         f.bcx[0], f.bcx[1] = tbl.bcx[1], tbl.bcx[2]
      end
      if tbl.bcy then
         f.bcy[0], f.bcy[1] = tbl.bcy[1], tbl.bcy[2]
      end
      if tbl.bcz then
         f.bcz[0], f.bcz[1] = tbl.bcz[1], tbl.bcz[2]
      end

      if tbl.bcx_lower_func then
         f.bcx_lower_func = gkyl_eval_bc(tbl.bcx_lower_func)
      end
      if tbl.bcy_lower_func then
         f.bcy_lower_func = gkyl_eval_bc(tbl.bcy_lower_func)
      end
      if tbl.bcz_lower_func then
         f.bcz_lower_func = gkyl_eval_bc(tbl.bcz_lower_func)
      end

      if tbl.bcx_upper_func then
         f.bcx_upper_func = gkyl_eval_bc(tbl.bcx_upper_func)
      end
      if tbl.bcy_upper_func then
         f.bcy_upper_func = gkyl_eval_bc(tbl.bcy_upper_func)
      end
      if tbl.bcz_upper_func then
         f.bcz_upper_func = gkyl_eval_bc(tbl.bcz_upper_func)
      end

      return f
   end,
   __index = {
      bcOpen = C.GKYL_FIELD_COPY,
      bcCopy = C.GKYL_FIELD_COPY,
      bcReflect = C.GKYL_FIELD_PEC_WALL,
      bcPEC = C.GKYL_FIELD_PEC_WALL,
      bcWedge = C.GKYL_FIELD_WEDGE,
      bcFunc = C.GKYL_FIELD_FUNC
   }
}
_M.Field = ffi.metatype(field_type, field_mt)

-- App
local app_type = ffi.typeof("struct gkyl_moment_app_cont")
local app_mt = {
   __new = function(self, tbl)
      local vm = ffi.new("struct gkyl_moment")

      local species = { }
      local field = nil

      -- first determine all species in system
      for k,v in pairs(tbl) do
	 if ffi.istype(species_type, v) then
	    v.name = k -- assign species name here
	    table.insert(species, v)
	 end
	 if ffi.istype(field_type, v) then
	    field = v -- only one field can be present
	 end
      end
      local num_species = #species

      local name = "moment"
      if GKYL_OUT_PREFIX then
         -- if G0 is being run from gkyl then GKYL_OUT_PREFIX is
         -- defined
         name = GKYL_OUT_PREFIX
      else
         local s, e = string.find(arg[0], ".lua")
         name = string.sub(arg[0], 1, s-1)
      end
	 
      -- set values in input struct
      vm.name = name
      vm.ndim = #tbl.cells

      -- set configuration space grid data
      for d=1, vm.ndim do
         vm.lower[d-1] = tbl.lower[d]
         vm.upper[d-1] = tbl.upper[d]
         vm.cells[d-1] = tbl.cells[d]
      end

      vm.cfl_frac = 0.95 -- for consistency with C app
      if tbl.cflFrac then
         vm.cfl_frac = tbl.cflFrac
      end

      vm.scheme_type = C.GKYL_MOMENT_WAVE_PROP
      if tbl.scheme_type then
	 vm.scheme_type = moment_scheme_tags[tbl.scheme_type]
      end
      if tbl.schemeType then
	 vm.scheme_type = moment_scheme_tags[tbl.schemeType]
      end

      vm.mp_recon = C.GKYL_MP_U5
      if tbl.mp_recon then
	 vm.mp_recon = mp_recon_tags[tbl.mp_recon]
      end
      if tbl.mpRecon then
	 vm.mp_recon = mp_recon_tags[tbl.mpRecon]
      end

      vm.skip_mp_limiter = pickBool(tbl.skipMpLimiter, false)

      vm.use_hybrid_flux_kep = false
      if tbl.useHybridFluxKep then
	 vm.use_hybrid_flux_kep = tbl.useHybridFluxKep
      end

      -- mapc2p
      vm.c2p_ctx = nil -- no need for context
      vm.mapc2p = nil
      if tbl.mapc2p then
         vm.mapc2p = gkyl_eval_mapc2p(tbl.mapc2p)
      end

      -- determine periodic BCs
      vm.num_periodic_dir = 0
      if tbl.periodicDirs then
         vm.num_periodic_dir = #tbl.periodicDirs
         for i=1, #tbl.periodicDirs do
          vm.periodic_dirs[i-1] = tbl.periodicDirs[i]-1 -- note indexing transforms
         end
      end

      -- determine directions to skip, if any
      vm.num_skip_dirs = 0
      if tbl.skip_dirs then
         vm.num_skip_dirs = #tbl.skip_dirs
         for i=1, #tbl.skip_dirs do
          vm.skip_dir[i-1] = tbl.skip_dir[i]
         end
      end

      -- set species
      vm.num_species = #species
      for i=1, #species do
         vm.species[i-1] = species[i]
      end

      -- set field
      if field then
         vm.field = field
      end      

      -- create new Moments app object
      local a = ffi.new(self)

      -- we need to store some stuff in container struct
      a.nspecies = num_species
      if tbl.tStart then
         a.t0 = tbl.tStart
      end
      a.tend = tbl.tEnd
      a.nframe = tbl.nFrame

      -- initialize app from input struct
      a.app = C.gkyl_moment_app_new(vm)
      return a
   end,
   __gc = function(self)
      C.gkyl_moment_app_release(self.app)
   end,
   __index = {
      init = function(self)
         C.gkyl_moment_app_apply_ic(self.app, self.t0)
      end,
      writeField = function(self, tm, frame)
         C.gkyl_moment_app_write_field(self.app, tm, frame)
      end,
      writeSpecies = function(self, tm, frame)
         for i=1, self.nspecies do
          C.gkyl_moment_app_write_species(self.app, i-1, tm, frame)
         end
      end,
      write = function(self, tm, frame)
         C.gkyl_moment_app_write(self.app, tm, frame)
      end,
      writeStat = function(self)
         C.gkyl_moment_app_stat_write(self.app)
      end,
      update = function(self, dt)
         return C.gkyl_moment_update(self.app, dt)
      end,
      calcIntegratedMom = function(self, tcurr)
	 return C.gkyl_moment_app_calc_integrated_mom(self.app, tcurr)
      end,
      calcFieldEnergy = function(self, tcurr)
	 return C.gkyl_moment_app_calc_field_energy(self.app, tcurr)
      end,      
      run = function(self)
	 
	 local frame_trig = _M.TimeTrigger(self.tend/self.nframe)

	 -- function to write data to file
	 local function writeData(tcurr)
	    if frame_trig:checkAndBump(tcurr) then
	       self:write(tcurr, frame_trig.curr-1)
	    end
	 end

	 local p1_trig = _M.TimeTrigger(self.tend/10)
	 -- log messages
	 local function writeLogMessage(tcurr, step, dt)
	    if p1_trig:checkAndBump(tcurr) then
	       io.write(string.format(" Step %6d %.4e. Time-step  %.6e \n", step, tcurr, dt))
	    end
	 end

	 io.write(string.format("Starting GkeyllZero simulation\n"))
	 io.write(string.format("  tstart: %.6e. tend: %.6e\n", 0.0, self.tend))

	 local tinit0 = _M.time_now()
	 self:init()
	 io.write(string.format("  Initialization completed in %g sec\n", _M.time_now() - tinit0))
	 
	 self:calcIntegratedMom(0.0)
	 self:calcFieldEnergy(0.0)
	 writeData(0.0)

	 local tloop0 = _M.time_now()
	 local tcurr, tend = 0.0, self.tend
	 local dt = tend-tcurr
	 local step = 1
	 while tcurr < tend do
	    local status = self:update(dt)
	    tcurr = tcurr + status.dt_actual

	    if status.success == false then
	       io.write(string.format("***** Simulation failed at step %5d at time %.6e\n", step, tcurr))
	       break
	    end

	    self:calcIntegratedMom(tcurr)
	    self:calcFieldEnergy(tcurr)	    

	    writeLogMessage(tcurr, step, status.dt_actual)
	    writeData(tcurr)

	    dt = math.min(status.dt_suggested, (tend-tcurr)*(1+1e-6))
	    step = step + 1
	 end

	 C.gkyl_moment_app_write_integrated_mom(self.app)
	 C.gkyl_moment_app_write_field_energy(self.app)

	 local tloop1 = _M.time_now()
	 
	 io.write(string.format("Completed in %d steps (tend: %.6e). \n", step-1, tcurr))
	 io.write(string.format("Main loop took %g secs to complete\n", tloop1-tloop0))
	 self:writeStat()
	 
      end,
   }
}
_M.App = ffi.metatype(app_type, app_mt)

return _M
