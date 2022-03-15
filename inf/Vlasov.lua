-- Thin wrapper around the gkyl_vlasov app. This tries to preserve
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

-- declare functions from Vlasov C app
ffi.cdef [[
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
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_PEC_WALL, // perfect electrical conductor (PEC) BCs
};

// This needs to be enum to allow usage below
enum { GKYL_MAX_SPECIES = 8 };

// Parameters for Vlasov species
struct gkyl_vlasov_species {
  char name[128]; // species name

  double charge, mass; // charge and mass
  double lower[3], upper[3]; // lower, upper bounds of velocity-space
  int cells[3]; // velocity-space cells

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  int num_diag_moments; // number of diagnostic moments
  char diag_moments[16][16]; // list of diagnostic moments

  // collision frequency
  void (*nu)(double t, const double *xn, double *fout, void *ctx);
  enum gkyl_collision_id collision_id; // type of collisions (see gkyl_eqn_type.h)

  void *accel_ctx; // context for applied acceleration function
  // pointer to applied acceleration function
  void (*accel)(double t, const double *xn, double *aout, void *ctx);

  // boundary conditions
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
};

// Parameter for EM field
struct gkyl_vlasov_field {
  enum gkyl_field_id field_id; // type of field (see gkyl_eqn_type.h)
  bool is_static; // set to true if field does not change in time

  double epsilon0, mu0;
  double elcErrorSpeedFactor, mgnErrorSpeedFactor;

  void *ctx; // context for initial condition init function
  // pointer to initialization function
  void (*init)(double t, const double *xn, double *fout, void *ctx);

  // boundary conditions
  enum gkyl_field_bc_type bcx[2], bcy[2], bcz[2];
};

// Top-level app parameters
struct gkyl_vm {
  char name[128]; // name of app: used as output prefix

  int cdim, vdim; // conf, velocity space dimensions
  double lower[3], upper[3]; // lower, upper bounds of config-space
  int cells[3]; // config-space cells
  int poly_order; // polynomial order
  enum gkyl_basis_type basis_type; // type of basis functions to use

  double cfl_frac; // CFL fraction to use (default 1.0)

  bool use_gpu; // Flag to indicate if solver should use GPUs

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  int num_species; // number of species
  struct gkyl_vlasov_species species[GKYL_MAX_SPECIES]; // species objects

  bool skip_field; // Skip field update or no field specified
  struct gkyl_vlasov_field field; // field object
};

// Object representing Vlasov app
typedef struct gkyl_vlasov_app gkyl_vlasov_app;

// Container to store pointer to app and other data
struct gkyl_vlasov_app_cont {
  double t0, tend; // start and end times
  int nframe; // number of frames to write
  int nspecies; // number of species

  gkyl_vlasov_app *app; // pointer to app
};

struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

/**
 * Construct a new Vlasov app.
 *
 * @param vm App inputs. See struct docs. All struct params MUST be
 *     initialized
 * @return New vlasov app object.
 */
gkyl_vlasov_app* gkyl_vlasov_app_new(struct gkyl_vm *vm);

/**
 * Initialize species and field by projecting initial conditions on
 * basis functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions.
 */
void gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0);

/**
 * Initialize field by projecting initial conditions on basis
 * functions.
 *
 * @param app App object.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0);

/**
 * Initialize species by projecting initial conditions on basis
 * functions. Species index (sidx) is the same index used to specify
 * the species in the gkyl_vm object used to construct app.
 *
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param t0 Time for initial conditions
 */
void gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0);

/**
 * Calculate diagnostic moments.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app);

/**
 * Write field and species data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write field data to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write species data to file.
 * 
 * @param app App object.
 * @param sidx Index of species to initialize.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame);

/**
 * Write diagnostic moments for species to file.
 * 
 * @param app App object.
 * @param tm Time-stamp
 * @param frame Frame number
 */
void gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame);

/**
 * Write stats to file. Data is written in json format.
 *
 * @param app App object.
 */
void gkyl_vlasov_app_stat_write(const gkyl_vlasov_app* app);

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
struct gkyl_update_status gkyl_vlasov_update(gkyl_vlasov_app* app, double dt);

/**
 * Free Vlasov app.
 *
 * @param app App to release.
 */
void gkyl_vlasov_app_release(gkyl_vlasov_app* app);

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

-- module table
local _M = { }


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

-- Wraps user given distribution function in a function that can be
-- passed to the C callback APIs
local function gkyl_eval_dist(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 6 do xnl[i] = xn[i-1] end -- might not be safe?
      fout[0] = func(t, xnl)
   end
end

-- Wraps user given acceleration function in a function that can be
-- passed to the C callback APIs
local function gkyl_eval_accel(func)
   return function(t, xn, fout, ctx)
      local xnl = ffi.new("double[10]")
      for i=1, 3 do xnl[i] = xn[i-1] end -- might not be safe?
      fout[0], fout[1], fout[2] = func(t, xnl)
   end
end

-- G0 basis types
local basis_type = {}
basis_type["serendipity"] = C.GKYL_BASIS_MODAL_SERENDIPITY
basis_type["tensor"] = C.GKYL_BASIS_MODAL_TENSOR

-- G0 boundary conditions for species: these functions ensure
-- compatibility with G2 boundary specifications
_M.OpenBC = function(tbl) return "bc_copy" end
_M.ReflectBC = function(tbl) return "bc_wall" end

-- G0 BC types for species
local species_bc_type = {}
species_bc_type["bc_copy"] = C.GKYL_SPECIES_COPY
species_bc_type["bc_wall"] = C.GKYL_SPECIES_WALL

-- G0 BC types for field
local field_bc_type = {}
field_bc_type["bc_copy"] = C.GKYL_FIELD_COPY
field_bc_type["bc_pec_wall"] = C.GKYL_FIELD_PEC_WALL

-- Species object: table structure is as follows:
--
-- local elc = {
--    charge = -1.0,
--    mass = 1.0,
--    lower = { -6.0 },
--    upper = { 6.0 },
--    cells = { 32 },
--    init = function(t, xn)
--       -- return distribution function value 
--    end
--    diagnostics = { "M0", "M1i" }
-- }

local species_vdim = 0 -- hack to store species VDIM

local species_type = ffi.typeof("struct gkyl_vlasov_species")
local species_mt = {
   __new = function(self, tbl)
      local s = ffi.new(self)
      s.charge = tbl.charge
      s.mass = tbl.mass

      species_vdim = #tbl.cells -- velocity space dimension
      local vdim = species_vdim
      -- velocity space grid
      for d=1, vdim do
	 s.lower[d-1] = tbl.lower[d]
	 s.upper[d-1] = tbl.upper[d]
	 s.cells[d-1] = tbl.cells[d]
      end

      -- initial conditions
      s.ctx = nil -- no need for a context
      s.init = gkyl_eval_dist(tbl.init)

      s.num_diag_moments = 0
      -- diagnostic moments
      if (tbl.diagnostics) then
	 s.num_diag_moments = #tbl.diagnostics
	 for i=1,#tbl.diagnostics do
	    s.diag_moments[i-1] = tbl.diagnostics[i]
	 end
      end

      -- collisions (at present, this is not how it is done in G2)
      s.nu = nil -- no need for a context
      if tbl.nu then
	 s.collision_id = C.GKYL_LBO_COLLISIONS -- at present, only LBO collisions
	 s.nu = gkyl_eval_dist(tbl.nu)
      end

      -- external acceleration
      s.accel_ctx = nil -- no need for a context
      s.accel = nil
      if tbl.accel then
	 s.accel = gkyl_eval_accel(tbl.accel)
      end

      -- boundary conditions
      if tbl.bcx then
	 s.bcx[0], s.bcx[1] = species_bc_type[tbl.bcx[1]], species_bc_type[tbl.bcx[2]]
      end
      if tbl.bcy then
	 s.bcy[0], s.bcy[1] = species_bc_type[tbl.bcy[1]], species_bc_type[tbl.bcy[2]]
      end
      if tbl.bcz then
	 s.bcz[0], s.bcz[1] = species_bc_type[tbl.bcz[1]], species_bc_type[tbl.bcz[2]]
      end

      return s
   end,
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

local field_type = ffi.typeof("struct gkyl_vlasov_field")
local field_mt = {
   __new = function(self, tbl)
      local f = ffi.new(self)
      f.epsilon0 = tbl.epsilon0
      f.mu0 = tbl.mu0

      f.elcErrorSpeedFactor = 0.0
      if (tbl.elcErrorSpeedFactor) then
	 f.elcErrorSpeedFactor = tbl.elcErrorSpeedFactor
      end

      f.mgnErrorSpeedFactor = 0.0
      if (tbl.mgnErrorSpeedFactor) then
	 f.mgnErrorSpeedFactor = tbl.mgnErrorSpeedFactor
      end

      f.ctx = nil -- no need for context
      f.init = gkyl_eval_field(tbl.init)

      -- boundary conditions
      if tbl.bcx then
	 f.bcx[0], f.bcx[1] = field_bc_type[tbl.bcx[1]], field_bc_type[tbl.bcx[2]]
      end
      if tbl.bcy then
	 f.bcy[0], f.bcy[1] = field_bc_type[tbl.bcy[1]], field_bc_type[tbl.bcy[2]]
      end
      if tbl.bcz then
	 f.bcz[0], f.bcz[1] = field_bc_type[tbl.bcz[1]], field_bc_type[tbl.bcz[2]]
      end

      f.is_static = false
      if (tbl.evolve ~= nil) then
	 f.is_static = not tbl.evolve
      end

      return f
   end,
   __index = {
      bcOpen = "bc_copy",
      bcCopy = "bc_copy",
      bcReflect = "bc_pec_wall"
   }
}
_M.Field = ffi.metatype(field_type, field_mt)

-- App
local app_type = ffi.typeof("struct gkyl_vlasov_app_cont")
local app_mt = {
   __new = function(self, tbl)
      local vm = ffi.new("struct gkyl_vm")

      local species = { }
      local field = nil

      -- first determine all species in system
      for k,v in pairs(tbl) do
	 if ffi.istype(species_type, v) then
	    v.name = k -- assign field name here
	    table.insert(species, v)
	 end
	 if ffi.istype(field_type, v) then
	    field = v -- only one field can be present
	 end
      end
      local num_species = #species

      local name = "vlasov"
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
      vm.cdim = #tbl.cells
      vm.vdim = species_vdim

      -- set configuration space grid data
      for d=1, vm.cdim do
	 vm.lower[d-1] = tbl.lower[d]
	 vm.upper[d-1] = tbl.upper[d]
	 vm.cells[d-1] = tbl.cells[d]
      end

      -- basis functions to use
      vm.basis_type = basis_type[tbl.basis]
      vm.poly_order = tbl.polyOrder

      -- CFL frac
      vm.cfl_frac = 1.0
      if tbl.cflFrac then
	 vm.cfl_frac = tbl.cflFrac
      end

      -- if we should be using GPUs
      vm.use_gpu = false
      if tbl.useGPU ~= nil then
	 vm.use_gpu = tbl.useGPU
      end

      -- determine periodic BCs
      vm.num_periodic_dir = 0
      if tbl.periodicDirs then
	 vm.num_periodic_dir = #tbl.periodicDirs
	 for i=1, #tbl.periodicDirs do
	    vm.periodic_dirs[i-1] = tbl.periodicDirs[i]-1 -- note indexing transforms
	 end
      end

      -- set species
      vm.num_species = #species
      for i=1, #species do
	 vm.species[i-1] = species[i]
      end

      -- set field
      vm.skip_field = true
      if field then
	 vm.skip_field = false
	 vm.field = field
      end

      -- create new Vlasov app object
      local a = ffi.new(self)

      -- we need to store some stuff in container struct
      a.nspecies = num_species
      if tbl.tStart then
	 a.t0 = tbl.tStart
      end
      a.tend = tbl.tEnd
      a.nframe = tbl.nFrame

      -- initialize app from input struct
      a.app = C.gkyl_vlasov_app_new(vm)
      return a
   end,
   __gc = function(self)
      C.gkyl_vlasov_app_release(self.app)
   end,
   __index = {
      init = function(self)
	 C.gkyl_vlasov_app_apply_ic(self.app, self.t0)
      end,
      writeField = function(self, tm, frame)
	 C.gkyl_vlasov_app_write_field(self.app, tm, frame)
      end,
      writeSpecies = function(self, tm, frame)
	 for i=1, self.nspecies do
	    C.gkyl_vlasov_app_write_species(self.app, i-1, tm, frame)
	 end
      end,
      writeMom = function(self, tm, frame)
	 C.gkyl_vlasov_app_write_mom(self.app, tm, frame)
      end,
      write = function(self, tm, frame)
	 C.gkyl_vlasov_app_write(self.app, tm, frame)
      end,
      writeStat = function(self)
	 C.gkyl_vlasov_app_stat_write(self.app)
      end,
      calcMom = function(self)
	 C.gkyl_vlasov_app_calc_mom(self.app)
      end,
      update = function(self, dt)
	 return C.gkyl_vlasov_update(self.app, dt)
      end,
      run = function(self)
	 
	 local frame_trig = _M.TimeTrigger(self.tend/self.nframe)

	 -- function to write data to file
	 local function writeData(tcurr)
	    if frame_trig:checkAndBump(tcurr) then
	       self:write(tcurr, frame_trig.curr-1)
	       self:calcMom()
	       self:writeMom(tcurr, frame_trig.curr-1)
	    end
	 end

	 local count = 0
	 local p1_trig = _M.TimeTrigger(self.tend/99)
	 -- log messages
	 local function writeLogMessage(tcurr, step, dt)
	    if p1_trig:checkAndBump(tcurr) then
	       if count % 10 == 0 then
		  io.write(string.format(" Step %6d %.4e. Time-step  %.6e \n", step, tcurr, dt))
	       end
	       count = count+1
	    end
	 end

	 io.write(string.format("Starting GkeyllZero simulation\n"))
	 self:init()
	 writeData(0.0)

	 local tcurr, tend = 0.0, self.tend
	 local dt = tend-tcurr
	 local step = 1
	 while tcurr < tend do
	    local status = self:update(dt);
	    tcurr = tcurr + status.dt_actual

	    dt = status.dt_suggested
	    writeLogMessage(tcurr, step, dt)
	    writeData(tcurr)

	    step = step + 1
	 end
	 io.write(string.format("Completed in %d steps. Final time-step %.6e\n", step-1, dt))

	 self:writeStat()
	 
      end,
   }
}
_M.App = ffi.metatype(app_type, app_mt)

return _M
