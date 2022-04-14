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
-- set JIT options
if jit.opt then
   jit.opt.start('callunroll=40', 'loopunroll=80', 'maxmcode=40960', 'maxtrace=100000',
		 'maxrecord=40000', 'maxside=1000', 'minstitch=3')
end

-- set global package paths
package.path = package.path .. ";" .. install_prefix .. "/lib/?.lua"

-- declare some top-level things we want to expose
ffi.cdef [[

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

// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_PHI, // Poisson (only phi)  
  GKYL_FIELD_PHI_A, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL // no field is present
};

/* Basis function identifiers */
enum gkyl_basis_type {
  GKYL_BASIS_MODAL_SERENDIPITY,
  GKYL_BASIS_MODAL_TENSOR,
  GKYL_BASIS_MODAL_GK_HYBRID,
};

/* Collision types */
enum gkyl_collision_id {
  GKYL_NO_COLLISIONS = 0, // No collisions. This is default
  GKYL_BGK_COLLISIONS, // BGK Collision operator
  GKYL_LBO_COLLISIONS // LBO Collision operator
};

// Boundary conditions on particles
enum gkyl_species_bc_type {
  GKYL_SPECIES_COPY = 0, // copy BCs
  GKYL_SPECIES_WALL, // perfect reflector
  GKYL_SPECIES_WEDGE, // specialized "wedge" BCs for RZ-theta
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // perfect electrical conductor (PEC) BCs
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
};

// some global constants
enum GKYL_GLOBAL_INTS {
  GKYL_MAX_CDIM = 3, // maximum configuration space dimensions
  GKYL_MAX_DIM = 7, // maximum phase-space dimensions
  GKYL_MAX_SPECIES = 8, // maximum number of species

  GKYL_DEF_ALIGN = 64, // default alignment boundary
};

struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};


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

]]

-- declare functions from Moments C app
ffi.cdef [[

enum gkyl_wave_limiter {
  GKYL_NO_LIMITER = 1, // to allow default to be 0
  GKYL_MIN_MOD,
  GKYL_SUPERBEE,
  GKYL_VAN_LEER,
  GKYL_MONOTONIZED_CENTERED,
  GKYL_BEAM_WARMING,
  GKYL_ZERO
};

// Wave equation object
typedef struct gkyl_wv_eqn gkyl_wv_eqn;

/**
 * Create a new Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_new(double gas_gamma);

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);

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
 * Create a new ideal MHD equation object.
 *
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to mhd equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_mhd_new(
    double gas_gamma, const char *divergence_constraint);

/**
 * Create a new SR Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to SR Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_sr_euler_new(double gas_gamma);

// Wrappers around various eqn object
struct gkyl_wv_euler { struct gkyl_wv_eqn *eqn; };
struct gkyl_wv_iso_euler { struct gkyl_wv_eqn *eqn; };
struct gkyl_wv_ten_moment { struct gkyl_wv_eqn *eqn; };

/**
 * Create a new Ten moment equation object.
 * 
 * @param k0 Closure parameter
 * @return Pointer to Ten moment equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_ten_moment_new(double k0);

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn* eqn);

// Parameters for moment species
struct gkyl_moment_species {
  char name[128]; // species name
  double charge, mass; // charge and mass
  enum gkyl_wave_limiter limiter; // limiter to use
  const struct gkyl_wv_eqn *equation; // equation object

  int evolve; // evolve species? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// Parameter for EM field
struct gkyl_moment_field {
  double epsilon0, mu0;
  double elc_error_speed_fact, mag_error_speed_fact;

  enum gkyl_wave_limiter limiter; // limiter to use

  int evolve; // evolve field? 1-yes, 0-no

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // boundary conditions
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];
};

// Choices of schemes to use in the fluid solver 
enum gkyl_moment_fluid_scheme {
  GKYL_MOMENT_FLUID_WAVE_PROP = 0, // default
  GKYL_MOMENT_FLUID_KEP
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

  enum gkyl_moment_fluid_scheme fluid_scheme; // scheme to update fluid equations

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_skip_dirs; // number of directions to skip
  int skip_dirs[3]; // directions to skip

  int num_species; // number of species
  struct gkyl_moment_species species[GKYL_MAX_SPECIES]; // species objects
  struct gkyl_moment_field field; // field object
};

// Simulation statistics
struct gkyl_moment_stat {
  long nup; // calls to update
  long nfail; // number of failed time-steps

  double total_tm; // time for simulation (not including ICs)
  double species_tm; // time to compute species updates
  double field_tm; // time to compute field updates
  double sources_tm; // time to compute source terms
};

// Object representing moments app
typedef struct gkyl_moment_app gkyl_moment_app;

// Container to store pointer to app and other data
struct gkyl_moment_app_cont {
  double t0, tend; // start and end times
  int nframe; // number of frames to write
  int nspecies; // number of species

  gkyl_moment_app *app; // pointer to app
};

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
void gkyl_moment_app_write_field(const gkyl_moment_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame);

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
struct gkyl_update_status gkyl_moment_update(gkyl_moment_app* app, double dt);

/**
 * Return simulation statistics.
 * 
 * @return Return statistics.
 */
struct gkyl_moment_stat gkyl_moment_app_stat(gkyl_moment_app* app);

/**
 * Free moment app.
 *
 * @param app App to release.
 */
void gkyl_moment_app_release(gkyl_moment_app* app);

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

-- Euler equation
_M.Euler = function(tbl)
   return ffi.gc(C.gkyl_wv_euler_new(tbl.gasGamma), C.gkyl_wv_eqn_release)
end

-- Wraps user given init function in a function that can be passed to
-- the C callback APIs
local function gkyl_eval_moment(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end -- might not be safe?
      local ret = { func(t, xnl) } -- package return into table
      for i=1,#ret do
	 fout[i-1] = ret[i]
      end
   end
end

local function gkyl_eval_mapc2p(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end -- might not be safe?
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
      s.ctx = nil -- no need for a context
      s.init = gkyl_eval_moment(tbl.init)

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

      return s
   end,
   __index = {
      -- we need this here also to be consistent with G2 App
      bcWall = C.GKYL_SPECIES_WALL,
      bcCopy = C.GKYL_SPECIES_COPY,
      bcWedge = C.GKYL_SPECIES_WEDGE
   }
}
_M.Species = ffi.metatype(species_type, species_mt)

-- Wraps user given field initialization function in a function that
-- can be passed to the C callback APIs
local function gkyl_eval_field(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 6 do xnl[i] = xn[i-1] end -- might not be safe?

      local ex,ey,ez,bx,by,bz = func(t, xnl)

      fout[0] = ex; fout[1] = ey; fout[2] = ez
      fout[3] = bx; fout[4] = by; fout[5] = bz
      fout[6] = 0.0; fout[7] = 0.0
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

      f.ctx = nil -- no need for context
      f.init = gkyl_eval_field(tbl.init)

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

      return f
   end,
   __index = {
      bcOpen = C.GKYL_FIELD_COPY,
      bcCopy = C.GKYL_FIELD_COPY,
      bcReflect = C.GKYL_FIELD_PEC_WALL,
      bcPEC = C.GKYL_FIELD_PEC_WALL,
      bcWedge = C.GKYL_FIELD_WEDGE
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

      vm.cfl_frac = 1.0
      if tbl.cflFrac then
	 vm.cfl_frac = tbl.cflFrac
      end

      vm.fluid_scheme = C.GKYL_MOMENT_FLUID_WAVE_PROP

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
	 self:init()
	 writeData(0.0)

	 local tcurr, tend = 0.0, self.tend
	 local dt = tend-tcurr
	 local step = 1
	 while tcurr < tend do
	    local status = self:update(dt);
	    tcurr = tcurr + status.dt_actual

	    writeLogMessage(tcurr, step, status.dt_actual)
	    writeData(tcurr)

	    dt = math.min(status.dt_suggested, (tend-tcurr)*(1+1e-6))
	    step = step + 1
	 end
	 io.write(string.format("Completed in %d steps (tend: %.6e). \n", step-1, tcurr))
	 self:writeStat()
	 
      end,
   }
}
_M.App = ffi.metatype(app_type, app_mt)

return _M
