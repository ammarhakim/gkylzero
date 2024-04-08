#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_moment.h>
#include <gkyl_moment_lw.h>
#include <gkyl_moment_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_wv_euler.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Magic IDs for use in distinguishing various species and field types
enum moment_magic_ids {
  MOMENT_SPECIES_DEFAULT = 100,
  MOMENT_FIELD_DEFAULT, // Maxwell equations
  MOMENT_EQN_DEFAULT
};

#define MOMENT_WAVE_EQN_METATABLE_NM "GkeyllZero.App.Moment.Eq"

/** For manipulating gkyl_wv_eqn objects */

// Lua userdata object for constructing field input
struct wv_eqn_lw {
  int magic; // must be first
  struct gkyl_wv_eqn *eqn;
};

// Clean up memory allocated for equation
static int
wv_eqn_lw_gc(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = GKYL_CHECK_UDATA(L, MOMENT_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;

  gkyl_wv_eqn_release(wv_lw->eqn);
  gkyl_free(*l_wv_lw);
  
  return 0;
}

static struct gkyl_wv_eqn*
wv_eqn_get(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = luaL_checkudata(L, -1, MOMENT_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;
  return wv_lw->eqn;
}

/* ****************/
/* Euler Equation */
/* ****************/

static enum gkyl_wv_euler_rp
euler_rp_type_from_str(const char *str)
{
  if (strcmp(str, "roe") == 0)
    return WV_EULER_RP_ROE;
  if (strcmp(str, "hllc") == 0)
    return WV_EULER_RP_HLLC;
  if (strcmp(str, "lax") == 0)
    return WV_EULER_RP_LAX;
  if (strcmp(str, "hll") == 0)
    return WV_EULER_RP_HLL;

  return WV_EULER_RP_ROE;
}

// Eq.Euler.new { gasgamma = 1.4, rpType = "roe" }
// rpType is one of "roe", "hllc", "lax", "hll"
static int
eqn_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_lw = gkyl_malloc(sizeof(*euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  const char *rp_str = glua_tbl_get_string(L, "rpType", "roe");
  enum gkyl_wv_euler_rp rp_type = euler_rp_type_from_str(rp_str);

  euler_lw->magic = MOMENT_EQN_DEFAULT;
  euler_lw->eqn = gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );

  // create Lua userdata ...
  struct wv_eqn_lw **l_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_lw = euler_lw; // ... point it to proper object
  
  // set metatable
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor
static struct luaL_Reg eqn_euler_ctor[] = {
  { "new",  eqn_euler_lw_new },
  { 0, 0 }
};

// Register and load all wave equation objects
static void
eqn_openlibs(lua_State *L)
{
  do {
    luaL_newmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, wv_eqn_lw_gc);
    lua_settable(L, -3);

    luaL_register(L, "G0.Moment.Eq.Euler", eqn_euler_ctor);
    
  } while (0);
}

/* *****************/
/* Species methods */
/* *****************/

// Metatable name for species input struct
#define MOMENT_SPECIES_METATABLE_NM "GkeyllZero.App.Moment.Species"

// Lua userdata object for constructing species input
struct moment_species_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_moment_species mom_species; // input struct to construct species
  bool evolve; // is this species evolved?
  struct lua_func_ctx init_ref; // Lua registery reference to initilization function
};

static int
moment_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_moment_species mom_species = { };

  mom_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  mom_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  bool has_eqn = false;
  with_lua_tbl_key(L, "equation") {
    mom_species.equation = wv_eqn_get(L);
    has_eqn = true;
  }

  if (!has_eqn)
    return luaL_error(L, "Species \"equation\" not specfied or incorrect type!");

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Species must have an \"init\" function for initial conditions!");

  struct moment_species_lw *moms_lw = lua_newuserdata(L, sizeof(*moms_lw));
  moms_lw->magic = MOMENT_SPECIES_DEFAULT;
  moms_lw->evolve = evolve;
  moms_lw->mom_species = mom_species;

  moms_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };
  
  // set metatable
  luaL_getmetatable(L, MOMENT_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg mom_species_ctor[] = {
  {"new", moment_species_lw_new},
  {0, 0}
};

/* *****************/
/* Field methods */
/* *****************/

// Metatable name for field input struct
#define MOMENT_FIELD_METATABLE_NM "GkeyllZero.App.Moment.Field"

// Lua userdata object for constructing field input
struct moment_field_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_moment_field mom_field; // input struct to construct field
  bool evolve; // is this field evolved?
  struct lua_func_ctx init_ref; // Lua registery reference to initilization function
};

static int
moment_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_moment_field mom_field = { };

  mom_field.epsilon0 = glua_tbl_get_number(L, "epsilon0", 1.0);
  mom_field.mu0 = glua_tbl_get_number(L, "mu0", 1.0);
  mom_field.elc_error_speed_fact = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  mom_field.mag_error_speed_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Field must have an \"init\" function for initial conditions!");

  struct moment_field_lw *vmf_lw = lua_newuserdata(L, sizeof(*vmf_lw));

  vmf_lw->magic = MOMENT_FIELD_DEFAULT;
  vmf_lw->evolve = evolve;
  vmf_lw->mom_field = mom_field;
  
  vmf_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // this will be set later
    .nret = 6,
    .L = L,
  };  
  
  // set metatable
  luaL_getmetatable(L, MOMENT_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg mom_field_ctor[] = {
  { "new",  moment_field_lw_new },
  { 0, 0 }
};

/* *************/
/* App methods */
/* *************/

// Metatable name for top-level Moment App
#define MOMENT_APP_METATABLE_NM "GkeyllZero.App.Moment"

// Lua userdata object for holding Moment app and run parameters
struct moment_app_lw {
  gkyl_moment_app *app; // Moment app object
  struct lua_func_ctx species_func_ctx[GKYL_MAX_SPECIES]; // function context for each species
  struct lua_func_ctx field_func_ctx; // function context for field
  
  double tstart, tend; // start and end times of simulation
  int nframe; // number of data frames to write
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct moment_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1};
  
  int curr = 0;
  lua_pushnil(L); // initial key is nil
  while (lua_next(L, TKEY) != 0) {
    // key at TKEY and value at TVAL
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct moment_species_lw *vms = lua_touserdata(L, TVAL);
      if (vms->magic == MOMENT_SPECIES_DEFAULT) {
        
        vms->init_ref.ndim = cdim;
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(vms->mom_species.name, key);
        }
        species[curr++] = vms;
      }
    }
    lua_pop(L, 1);
  }
  return curr;
}

// Create top-level App object
static int
mom_app_new(lua_State *L)
{
  struct moment_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-vlasov" is used.
  const char *sim_name = "g0-vlasov";
  with_lua_global(L, "GKYL_OUT_PREFIX") {
    if (lua_isstring(L, -1))
      sim_name = lua_tostring(L, -1);
  }
  
  // initialize app using table inputs (table is on top of stack)

  app_lw->tstart = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->tend = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->nframe = glua_tbl_get_integer(L, "nFrame", 1);

  struct gkyl_moment mom = { }; // input table for app

  strcpy(mom.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    mom.ndim = cdim = glua_objlen(L);
    for (int d=0; d<cdim; ++d)
      mom.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  int cuts[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) cuts[d] = 1;
  
  with_lua_tbl_tbl(L, "decompCuts") {
    int ncuts = glua_objlen(L);
    for (int d=0; d<ncuts; ++d)
      cuts[d] = glua_tbl_iget_integer(L, d+1, 0);
  }  

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<cdim; ++d)
      mom.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<cdim; ++d)
      mom.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  mom.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);

  mom.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      mom.num_periodic_dir = glua_objlen(L);
      for (int d=0; d<mom.num_periodic_dir; ++d)
        // indexes are off by 1 between Lua and C
        mom.periodic_dirs[d] = glua_tbl_iget_integer(L, d+1, 0)-1;
    }
  }

  struct moment_species_lw *species[GKYL_MAX_SPECIES];
  // set all species input
  mom.num_species = get_species_inp(L, cdim, species);
  for (int s=0; s<mom.num_species; ++s) {
    mom.species[s] = species[s]->mom_species;
    
    app_lw->species_func_ctx[s] = species[s]->init_ref;
    mom.species[s].init = gkyl_lw_eval_cb;
    mom.species[s].ctx = &app_lw->species_func_ctx[s];
  }

  // set field input
  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct moment_field_lw *vmf = lua_touserdata(L, -1);
      if (vmf->magic == MOMENT_FIELD_DEFAULT) {

        vmf->init_ref.ndim = cdim;

        mom.field = vmf->mom_field;

        app_lw->field_func_ctx = vmf->init_ref;
        mom.field.init = gkyl_lw_eval_cb;
        mom.field.ctx = &app_lw->field_func_ctx;
      }
    }
  }

  // create decomp and communicator
  struct gkyl_rect_decomp *decomp
    = gkyl_rect_decomp_new_from_cuts_and_cells(cdim, cuts, mom.cells);

  struct gkyl_comm *comm = 0;

  int rank = 0;

  bool has_mpi = false;
#ifdef GKYL_HAVE_MPI
  with_lua_global(L, "GKYL_MPI_COMM") {
    if (lua_islightuserdata(L, -1)) {
      has_mpi = true;
      MPI_Comm mpi_comm = lua_touserdata(L, -1);
      comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
          .mpi_comm = mpi_comm,
          .decomp = decomp
        }
      );
      MPI_Comm_rank(mpi_comm, &rank);
    }
  }
#endif

  if (!has_mpi) {
    // if there is no proper MPI_Comm specifed, the assume we are a
    // serial sim
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .sync_corners = true
      }
    );
  }

  mom.has_low_inp = true;  
  mom.low_inp = (struct gkyl_app_comm_low_inp) {
    .comm = comm,
    .local_range = decomp->ranges[rank]
  };
  
  app_lw->app = gkyl_moment_app_new(&mom); // create the Moment app

  gkyl_rect_decomp_release(decomp);
  if (comm) gkyl_comm_release(comm);
  
  // create Lua userdata ...
  struct moment_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct moment_app_lw*));
  *l_app_lw = app_lw; // ... point it to the Lua app pointer

  // set metatable
  luaL_getmetatable(L, MOMENT_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool
static int
mom_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_moment_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool
static int
mom_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_moment_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool
static int
mom_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->tstart);
  gkyl_moment_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool
static int
mom_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_moment_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool
static int
mom_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_moment_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool
static int
mom_app_write(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_moment_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool
static int
mom_app_write_field(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_moment_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool
static int
mom_app_write_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_moment_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool
static int
mom_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool
static int
mom_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool
static int
mom_app_stat_write(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file
static void
write_data(struct gkyl_tm_trigger *iot, gkyl_moment_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
}

// Run simulation. (num_steps) -> bool. num_steps is optional
static int
mom_app_run(lua_State *L)
{
  bool ret_status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;
  struct gkyl_moment_app *app = app_lw->app;

  double tcurr = app_lw->tstart;
  double tend = app_lw->tend;
  double dt = tend-tcurr;
  long num_steps = luaL_optinteger(L, 2, INT_MAX);

  int nframe = app_lw->nframe;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_moment_app_calc_integrated_mom(app, tcurr);
  gkyl_moment_app_calc_field_energy(app, tcurr);

  long step = 1;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_moment_app_calc_integrated_mom(app, tcurr);
    gkyl_moment_app_calc_field_energy(app, tcurr);    
    
    if (!status.success) {
      ret_status = false;
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_moment_app_stat_write(app);
  gkyl_moment_app_write_integrated_mom(app);
  gkyl_moment_app_write_field_energy(app);

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation
static int
mom_app_gc(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor
static struct luaL_Reg mom_app_ctor[] = {
  { "new",  mom_app_new },
  { 0, 0 }
};

// App methods
static struct luaL_Reg mom_app_funcs[] = {
  { "apply_ic", mom_app_apply_ic },
  { "apply_ic_field", mom_app_apply_ic_field },
  { "apply_ic_species", mom_app_apply_ic_species },
  { "calc_integrated_mom", mom_app_calc_integrated_mom },
  { "calc_field_energy", mom_app_calc_field_energy },
  { "write", mom_app_write },
  { "write_field", mom_app_write_field },
  { "write_species", mom_app_write_species },
  { "write_integrated_mom", mom_app_write_integrated_mom },
  { "write_field_energy", mom_app_write_field_energy },
  { "stat_write", mom_app_stat_write },
  { "run", mom_app_run },
  { 0, 0 }
};


static void
app_openlibs(lua_State *L)
{
  // Register top-level App
  do {
    luaL_newmetatable(L, MOMENT_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, mom_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, mom_app_funcs);
    
    luaL_register(L, "G0.Moment.App", mom_app_ctor);
    
  } while (0);

  // Register Species input struct
  do {
    luaL_newmetatable(L, MOMENT_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Moment.Species", mom_species_ctor);
  } while (0);

  // Register Field input struct
  do {
    luaL_newmetatable(L, MOMENT_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Moment.Field", mom_field_ctor);
  } while (0);
}

void
gkyl_moment_lw_openlibs(lua_State *L)
{
  eqn_openlibs(L);
  app_openlibs(L);
}

#endif
