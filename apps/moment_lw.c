#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_moment.h>
#include <gkyl_moment_lw.h>
#include <gkyl_moment_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_sr_euler.h>
#include <gkyl_wv_ten_moment.h>

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

// (string, int) pair
struct str_int_pair {
  const char *str;
  int val;
};

// wave limiter -> enum map
static const struct str_int_pair wave_limiter[] = {
  { "no-limiter", GKYL_NO_LIMITER },
  { "min-mod", GKYL_MIN_MOD },
  { "superbee", GKYL_SUPERBEE },
  { "van-leer", GKYL_VAN_LEER },
  { "beam-warming", GKYL_BEAM_WARMING },
  { "zero", GKYL_ZERO },
  { 0, 0 }
};

// edge-splitting -> enum map
static const struct str_int_pair wave_split_type[] = {
  { "qwave", GKYL_WAVE_QWAVE },
  { "fwave", GKYL_WAVE_FWAVE },
  { 0, 0 }
};

// RP -> enum map
static const struct str_int_pair euler_rp_type[] = {
  { "roe", WV_EULER_RP_ROE },
  { "hllc", WV_EULER_RP_HLLC },
  { "lax", WV_EULER_RP_LAX },
  { "hll", WV_EULER_RP_HLL },
  { 0, 0 }
};

// simple linear search in list of pairs
static int
search_str_int_pair_by_str(const struct str_int_pair pairs[], const char *str, int def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (strcmp(pairs[i].str, str) == 0)
      return pairs[i].val;
  }
  return def;
}

// simple linear search in list of pairs
static const char *
search_str_int_pair_by_int(const struct str_int_pair pairs[], int val, const char *def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (pairs[i].val == val)
      return pairs[i].str;
  }
  return def;
}

#define MOMENT_WAVE_EQN_METATABLE_NM "GkeyllZero.App.Moments.Eq"

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

// Euler.new { gasgamma = 1.4, rpType = "roe" }
// rpType is one of "roe", "hllc", "lax", "hll"
static int
eqn_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_lw = gkyl_malloc(sizeof(*euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  const char *rp_str = glua_tbl_get_string(L, "rpType", "roe");
  enum gkyl_wv_euler_rp rp_type = search_str_int_pair_by_str(euler_rp_type, rp_str, WV_EULER_RP_ROE);

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
  {"new", eqn_euler_lw_new},
  {0, 0}
};

/* ***************************/
/* Isothermal Euler Equation */
/* ***************************/

// IsoEuler.new { vThermal = 1.0 }
static int
eqn_iso_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *iso_euler_lw = gkyl_malloc(sizeof(*iso_euler_lw));

  double vt = glua_tbl_get_number(L, "vThermal", -1.0);
  if (vt < 0)
    return luaL_error(L, "Thermal velocity \"vThermal\" not specified properly!");
  
  iso_euler_lw->magic = MOMENT_EQN_DEFAULT;
  iso_euler_lw->eqn = gkyl_wv_iso_euler_new(vt);

  // create Lua userdata ...
  struct wv_eqn_lw **l_iso_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_iso_euler_lw = iso_euler_lw; // ... point it to proper object
  
  // set metatable
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor
static struct luaL_Reg eqn_iso_euler_ctor[] = {
  {"new", eqn_iso_euler_lw_new},
  {0, 0}
};

/* ********************/
/* Tenmoment Equation */
/* ********************/

// Tenmoment.new { k0 = 1, has_grad_closure = false }
static int
eqn_tenmoment_lw_new(lua_State *L)
{
  struct wv_eqn_lw *tenm_lw = gkyl_malloc(sizeof(*tenm_lw));

  double k0 = glua_tbl_get_number(L, "k0", 1.0);
  bool has_grad_closure = glua_tbl_get_bool(L, "has_grad_closure", false);

  tenm_lw->magic = MOMENT_EQN_DEFAULT;
  tenm_lw->eqn = gkyl_wv_ten_moment_new(k0, has_grad_closure);

  // create Lua userdata ...
  struct wv_eqn_lw **l_tenm_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_tenm_lw = tenm_lw; // ... point it to proper object
  
  // set metatable
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor
static struct luaL_Reg eqn_tenmoment_ctor[] = {
  { "new",  eqn_tenmoment_lw_new },
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

    luaL_register(L, "G0.Moments.Eq.Euler", eqn_euler_ctor);
    luaL_register(L, "G0.Moments.Eq.IsoEuler", eqn_iso_euler_ctor);
    luaL_register(L, "G0.Moments.Eq.TenMoment", eqn_tenmoment_ctor);
  } while (0);
}

/* *****************/
/* Species methods */
/* *****************/

// Metatable name for species input struct
#define MOMENT_SPECIES_METATABLE_NM "GkeyllZero.App.Moments.Species"

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

  bool has_eqn = false;
  with_lua_tbl_key(L, "equation") {
    mom_species.equation = wv_eqn_get(L);
    has_eqn = true;
  }

  if (!has_eqn)
    return luaL_error(L, "Species \"equation\" not specfied or incorrect type!");

  const char *lim_str = glua_tbl_get_string(L, "limiter", "monotonized-centered");
  mom_species.limiter = search_str_int_pair_by_str(wave_limiter, lim_str, GKYL_MONOTONIZED_CENTERED);

  const char *split_str = glua_tbl_get_string(L, "split_type", "qwave");
  mom_species.split_type = search_str_int_pair_by_str(wave_split_type, split_str, GKYL_WAVE_QWAVE);

  bool evolve = mom_species.evolve = glua_tbl_get_bool(L, "evolve", true);
  mom_species.force_low_order_flux = glua_tbl_get_bool(L, "forceLowOrderFlux", false);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Species must have an \"init\" function for initial conditions!");

  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);
    for (int i=0; i < (nbc>2 ? 2 : nbc); ++i)
      mom_species.bcx[i] = glua_tbl_iget_integer(L, i+1, 0);
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);
    for (int i=0; i < (nbc>2 ? 2 : nbc); ++i)
      mom_species.bcy[i] = glua_tbl_iget_integer(L, i+1, 0);
  }  

  struct moment_species_lw *moms_lw = lua_newuserdata(L, sizeof(*moms_lw));
  moms_lw->magic = MOMENT_SPECIES_DEFAULT;
  moms_lw->evolve = evolve;
  moms_lw->mom_species = mom_species;

  moms_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // this will be set later
    .nret = mom_species.equation->num_equations,
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
#define MOMENT_FIELD_METATABLE_NM "GkeyllZero.App.Moments.Field"

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
  mom_field.mag_error_speed_fact = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 1.0);

  const char *lim_str = glua_tbl_get_string(L, "limiter", "monotonized-centered");  
  mom_field.limiter = search_str_int_pair_by_str(wave_limiter, lim_str, GKYL_MONOTONIZED_CENTERED);
  
  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Field must have an \"init\" function for initial conditions!");

  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);
    for (int i=0; i < (nbc>2 ? 2 : nbc); ++i)
      mom_field.bcx[i] = glua_tbl_iget_integer(L, i+1, 0);
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);
    for (int i=0; i < (nbc>2 ? 2 : nbc); ++i)
      mom_field.bcy[i] = glua_tbl_iget_integer(L, i+1, 0);
  }  

  struct moment_field_lw *momf_lw = lua_newuserdata(L, sizeof(*momf_lw));

  momf_lw->magic = MOMENT_FIELD_DEFAULT;
  momf_lw->evolve = evolve;
  momf_lw->mom_field = mom_field;
  
  momf_lw->init_ref = (struct lua_func_ctx) {
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

// find index of species in table given its name
static int
app_find_species(const gkyl_moment_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(app->species[i].name, nm) == 0)
      return i;
  return -1;
}

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
      struct moment_field_lw *momf = lua_touserdata(L, -1);
      if (momf->magic == MOMENT_FIELD_DEFAULT) {

        momf->init_ref.ndim = cdim;

        mom.field = momf->mom_field;

        app_lw->field_func_ctx = momf->init_ref;
        mom.field.init = gkyl_lw_eval_cb;
        mom.field.ctx = &app_lw->field_func_ctx;
      }
    }
  }

  // create decomp and communicator
  struct gkyl_rect_decomp *decomp
    = gkyl_rect_decomp_new_from_cuts_and_cells(cdim, cuts, mom.cells);

  struct gkyl_comm *comm = 0;
  bool has_mpi = false;
  
#ifdef GKYL_HAVE_MPI
  with_lua_global(L, "GKYL_MPI_COMM") {
    if (lua_islightuserdata(L, -1)) {
      has_mpi = true;
      MPI_Comm mpi_comm = lua_touserdata(L, -1);
      comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
          .mpi_comm = mpi_comm,
          .decomp = decomp,
          .sync_corners = true
        }
      );

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

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int tot_cuts = 1; for (int d=0; d<cdim; ++d) tot_cuts *= cuts[d];

  if (tot_cuts != comm_sz) {
    printf("tot_cuts = %d (%d)\n", tot_cuts, comm_sz);
    luaL_error(L, "Number of ranks and cuts do not match!");
  }
  
  mom.has_low_inp = true;  
  mom.low_inp = (struct gkyl_app_comm_low_inp) {
    .comm = comm,
    .local_range = decomp->ranges[rank]
  };
  
  app_lw->app = gkyl_moment_app_new(&mom); // create the Moment app

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  
  // create Lua userdata ...
  struct moment_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct moment_app_lw*));
  *l_app_lw = app_lw; // ... point it to the Lua app pointer

  // set metatable
  luaL_getmetatable(L, MOMENT_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Compute maximum time-step () -> double
static int
mom_app_max_dt(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double maxdt = gkyl_moment_app_max_dt(app_lw->app);

  lua_pushnumber(L, maxdt);  
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

// Apply initial conditions to species. (name, time) -> bool
static int
mom_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *sp_name = luaL_checkstring(L, 2);
  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0)
    return luaL_error(L, "Incorrect species name '%s' in apply_ic_species!", sp_name);
  
  double t0 = luaL_optnumber(L, 3, app_lw->tstart);
  gkyl_moment_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// The status table for reads has the following structure
//
// {
//   io_status = true or false.
//   io_status_str = "success", "bad-version", "fopen-failed", "fread-failed", "data-mismatch"
//   frame = frame_num (integer)
//   stime = time in file read
// }

static const struct str_int_pair rio_status[] = {
  { "success", GKYL_ARRAY_RIO_SUCCESS },
  { "bad-version", GKYL_ARRAY_RIO_BAD_VERSION },
  { "fopen-failed", GKYL_ARRAY_RIO_FOPEN_FAILED },
  { "fread-failed", GKYL_ARRAY_RIO_FREAD_FAILED },
  { "data-mismatch", GKYL_ARRAY_RIO_DATA_MISMATCH },
  { 0, 0 }
};

// pushes table with status on stack. Table is left on stack
static void
push_restart_status_table(lua_State *L, struct gkyl_app_restart_status status)
{
  lua_newtable(L);

  lua_pushstring(L, "io_status");
  lua_pushboolean(L, status.io_status == GKYL_ARRAY_RIO_SUCCESS ? true : false);
  lua_rawset(L, -3);

  lua_pushstring(L, "io_status_str");
  lua_pushstring(L, search_str_int_pair_by_int(rio_status, status.io_status, "success"));
  lua_rawset(L, -3);

  lua_pushstring(L, "frame");
  lua_pushinteger(L, status.frame);
  lua_rawset(L, -3);

  lua_pushstring(L, "stime");
  lua_pushnumber(L, status.stime);
  lua_rawset(L, -3);  
}

// Read field from file. (file-name) -> status table. See above
static int
mom_app_from_file_field(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *fname = luaL_checkstring(L, 2);
  struct gkyl_app_restart_status status =
    gkyl_moment_app_from_file_field(app_lw->app, fname);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (file-name, species-name) -> status table. See above
static int
mom_app_from_file_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *fname = luaL_checkstring(L, 2);
  const char *sp_name = luaL_checkstring(L, 3);

  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0)
    return luaL_error(L, "Incorrect species name '%s' in from_file_species!", sp_name);
  
  struct gkyl_app_restart_status status =
    gkyl_moment_app_from_file_species(app_lw->app, sidx, fname);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (frame) -> status table. See above
static int
mom_app_from_frame_field(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int frame = luaL_checkinteger(L, 2);
  struct gkyl_app_restart_status status =
    gkyl_moment_app_from_frame_field(app_lw->app, frame);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (frame, species-name) -> status table. See above
static int
mom_app_from_frame_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int frame = luaL_checkinteger(L, 2);
  const char *sp_name = luaL_checkstring(L, 3);

  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0)
    return luaL_error(L, "Incorrect species name '%s' in from_frame_species!", sp_name);
  
  struct gkyl_app_restart_status status =
    gkyl_moment_app_from_frame_species(app_lw->app, sidx, frame);

  push_restart_status_table(L, status);
  
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

// Write species solution to file (name, time, frame) -> bool
static int
mom_app_write_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *sp_name = luaL_checkstring(L, 2);
  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0)
    return luaL_error(L, "Incorrect species name '%s' in write_species!", sp_name);

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
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
    gkyl_moment_app_write_integrated_mom(app);
    gkyl_moment_app_write_field_energy(app);
  }
}

// Update status table is as follows
// {
//    success = true or false
//    dt_actual = actual time-step taken
//    dt_suggested = suggested stable time-step
// }

// Update the solution by a suggested time-step. (dt) -> update status
// table. See above. For details see the C API doc for this function
static int
mom_app_update(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double dt = luaL_checknumber(L, 2);
  struct gkyl_update_status status = gkyl_moment_update(app_lw->app, dt);

  // return status table on stack
  lua_newtable(L);
  lua_pushstring(L, "success");
  lua_pushboolean(L, status.success);
  lua_rawset(L, -3);  

  lua_pushstring(L, "dt_actual");
  lua_pushnumber(L, status.dt_actual);
  lua_rawset(L, -3);    

  lua_pushstring(L, "dt_suggested");
  lua_pushnumber(L, status.dt_suggested);
  lua_rawset(L, -3);
  
  return 1;
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
  gkyl_moment_app_calc_integrated_mom(app, tcurr);
  gkyl_moment_app_calc_field_energy(app, tcurr);
  
  write_data(&io_trig, app, tcurr);

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
  {"max_dt", mom_app_max_dt},

  {"apply_ic", mom_app_apply_ic},
  {"apply_ic_field", mom_app_apply_ic_field},
  {"apply_ic_species", mom_app_apply_ic_species},
  
  {"from_file_field", mom_app_from_file_field},
  {"from_file_species", mom_app_from_file_species},

  {"from_frame_field", mom_app_from_frame_field},
  {"from_frame_species", mom_app_from_frame_species},

  {"write", mom_app_write},
  {"write_field", mom_app_write_field},
  {"write_species", mom_app_write_species},
  {"write_field_energy", mom_app_write_field_energy},
  {"write_integrated_mom", mom_app_write_integrated_mom},
  {"stat_write", mom_app_stat_write},

  {"calc_field_energy", mom_app_calc_field_energy},
  {"calc_integrated_mom", mom_app_calc_integrated_mom},

  { "update", mom_app_update },
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
    
    luaL_register(L, "G0.Moments.App", mom_app_ctor);
    
  } while (0);

  // Register Species input struct
  do {
    luaL_newmetatable(L, MOMENT_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Moments.Species", mom_species_ctor);
  } while (0);

  // Register Field input struct
  do {
    luaL_newmetatable(L, MOMENT_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Moments.Field", mom_field_ctor);
  } while (0);
}

void
gkyl_moment_lw_openlibs(lua_State *L)
{
  eqn_openlibs(L);
  app_openlibs(L);
}

#endif