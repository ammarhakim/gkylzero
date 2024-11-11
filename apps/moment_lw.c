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
#include <gkyl_wv_coldfluid.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_sr_euler.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_reactive_euler.h>
#include <gkyl_wv_euler_mixture.h>
#include <gkyl_wv_iso_euler_mixture.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Magic IDs for use in distinguishing various species and field types.
enum moment_magic_ids {
  MOMENT_SPECIES_DEFAULT = 100, // Fluid species.
  MOMENT_FIELD_DEFAULT, // Maxwell equations.
  MOMENT_EQN_DEFAULT // Equation object.
};

// Wave limiter -> enum map.
static const struct gkyl_str_int_pair wave_limiter[] = {
  { "no-limiter", GKYL_NO_LIMITER },
  { "min-mod", GKYL_MIN_MOD },
  { "superbee", GKYL_SUPERBEE },
  { "van-leer", GKYL_VAN_LEER },
  { "beam-warming", GKYL_BEAM_WARMING },
  { "zero", GKYL_ZERO },
  { 0, 0 }
};

// Edge-splitting -> enum map.
static const struct gkyl_str_int_pair wave_split_type[] = {
  { "qwave", GKYL_WAVE_QWAVE },
  { "fwave", GKYL_WAVE_FWAVE },
  { 0, 0 }
};

// Ideal MHD Riemann problem -> enum map.
static const struct gkyl_str_int_pair mhd_rp_type[] = {
  { "roe", WV_MHD_RP_ROE },
  { "hlld", WV_MHD_RP_HLLD },
  { "lax", WV_MHD_RP_LAX },
  { 0, 0 }
};

// Ideal divB correction -> enum map.
static const struct gkyl_str_int_pair mhd_divb_type[] = {
  { "none", GKYL_MHD_DIVB_NONE },
  { "glm",  GKYL_MHD_DIVB_GLM },
  { "eight_waves", GKYL_MHD_DIVB_EIGHT_WAVES },
  { 0, 0 }
};

// Reactive Euler Riemann problem -> enum map.
static const struct gkyl_str_int_pair reactive_euler_rp_type[] = {
  { "roe", WV_REACTIVE_EULER_RP_ROE },
  { "lax", WV_REACTIVE_EULER_RP_LAX },
  { 0, 0 }
};

// Euler mixture Riemann problem -> enum map.
static const struct gkyl_str_int_pair euler_mixture_rp_type[] = {
  { "roe", WV_EULER_MIXTURE_RP_ROE },
  { "lax", WV_EULER_MIXTURE_RP_LAX },
  { 0, 0 }
};

// Isothermal Euler mixture Riemann problem -> enum map.
static const struct gkyl_str_int_pair iso_euler_mixture_rp_type[] = {
  { "roe", WV_ISO_EULER_MIXTURE_RP_ROE },
  { "lax", WV_ISO_EULER_MIXTURE_RP_LAX },
  { 0, 0 }
};

// Metatable name for equation object input struct.
#define MOMENT_WAVE_EQN_METATABLE_NM "GkeyllZero.App.Moments.Eq"

// Methods for manipulating gkyl_wv_eqn objects.

// Lua userdata object for constructing wave equation objects.
struct wv_eqn_lw {
  int magic; // This must be the first element in the struct.
  struct gkyl_wv_eqn *eqn; // Equation object.
};

// Clean up memory allocated for equation object.
static int
wv_eqn_lw_gc(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = GKYL_CHECK_UDATA(L, MOMENT_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;

  gkyl_wv_eqn_release(wv_lw->eqn);
  gkyl_free(*l_wv_lw);
  
  return 0;
}

// Acquire equation object.
static struct gkyl_wv_eqn*
wv_eqn_get(lua_State *L)
{
  struct wv_eqn_lw **l_wv_lw = luaL_checkudata(L, -1, MOMENT_WAVE_EQN_METATABLE_NM);
  struct wv_eqn_lw *wv_lw = *l_wv_lw;
  return wv_lw->eqn;
}

/* *************** */
/* Euler Equations */
/* *************** */

// Euler.new { gasGamma = 1.4, rpType = G0.EulerRP.Roe }
// where rpType is one of G0.EulerRP.Roe, G0.EulerRP.Lax, G0.EulerRP.HLL or G0.EulerRP.HLLC.
static int
eqn_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_lw = gkyl_malloc(sizeof(*euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  enum gkyl_wv_euler_rp rp_type = glua_tbl_get_integer(L, "rpType", WV_EULER_RP_ROE);

  euler_lw->magic = MOMENT_EQN_DEFAULT;
  euler_lw->eqn = gkyl_wv_euler_inew( &(struct gkyl_wv_euler_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct wv_eqn_lw **l_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_lw = euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_euler_ctor[] = {
  { "new", eqn_euler_lw_new },
  { 0, 0 }
};

/* ************************** */
/* Isothermal Euler Equations */
/* ************************** */

// IsoEuler.new { vThermal = 1.0 }
static int
eqn_iso_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *iso_euler_lw = gkyl_malloc(sizeof(*iso_euler_lw));

  double vt = glua_tbl_get_number(L, "vThermal", -1.0);
  if (vt < 0) {
    return luaL_error(L, "Thermal velocity \"vThermal\" not specified properly!");
  }
  
  iso_euler_lw->magic = MOMENT_EQN_DEFAULT;
  iso_euler_lw->eqn = gkyl_wv_iso_euler_new(vt);

  // Create Lua userdata.
  struct wv_eqn_lw **l_iso_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_iso_euler_lw = iso_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_iso_euler_ctor[] = {
  {"new", eqn_iso_euler_lw_new},
  {0, 0}
};

/* ************************************ */
/* Special Relativistic Euler Equations */
/* ************************************ */

// SrEuler.new { gasgamma = 1.4 }
static int
eqn_sr_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *sr_euler_lw = gkyl_malloc(sizeof(*sr_euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 5.0 / 3.0);

  sr_euler_lw->magic = MOMENT_EQN_DEFAULT;
  sr_euler_lw->eqn = gkyl_wv_sr_euler_new(gas_gamma);

  // Create Lua userdata.
  struct wv_eqn_lw **l_sr_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_sr_euler_lw = sr_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_sr_euler_ctor[] = {
  { "new", eqn_sr_euler_lw_new },
  { 0, 0 }
};

/* ******************** */
/* Cold Fluid Equations */
/* ******************** */

// ColdFluid.new {  }
static int
eqn_coldfluid_lw_new(lua_State *L)
{
  struct wv_eqn_lw *coldfluid_lw = gkyl_malloc(sizeof(*coldfluid_lw));

  coldfluid_lw->magic = MOMENT_EQN_DEFAULT;
  coldfluid_lw->eqn = gkyl_wv_coldfluid_new();

  // Create Lua userdata.
  struct wv_eqn_lw **l_coldfluid_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_coldfluid_lw = coldfluid_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_coldfluid_ctor[] = {
  { "new", eqn_coldfluid_lw_new },
  { 0, 0 }
};

/* ******************** */
/* Ten-moment Equations */
/* ******************** */

// TenMoment.new { k0 = 1, hasGradClosure = false }
static int
eqn_tenmoment_lw_new(lua_State *L)
{
  struct wv_eqn_lw *tenm_lw = gkyl_malloc(sizeof(*tenm_lw));

  double k0 = glua_tbl_get_number(L, "k0", 1.0);
  bool has_grad_closure = glua_tbl_get_bool(L, "hasGradClosure", false);

  tenm_lw->magic = MOMENT_EQN_DEFAULT;
  tenm_lw->eqn = gkyl_wv_ten_moment_new(k0, has_grad_closure, false);

  // Create Lua userdata.
  struct wv_eqn_lw **l_tenm_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_tenm_lw = tenm_lw; // Point userdata to the equation object.
  
  // Set metatable
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_tenmoment_ctor[] = {
  { "new", eqn_tenmoment_lw_new },
  { 0, 0 }
};

/* ************* */
/* MHD Equations */
/* ************* */

// Mhd.new { gasgamma = 1.4, rpType = "roe", divB = "glm", glmCh = 0.0, glmAlpha = 0.0 }
// rpType is one of "roe", "hlld", "lax"
// divB is "none", "glm", "eight_waves"
static int
eqn_mhd_lw_new(lua_State *L)
{
  struct wv_eqn_lw *mhd_lw = gkyl_malloc(sizeof(*mhd_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  
  const char *rp_str = glua_tbl_get_string(L, "rpType", "roe");
  enum gkyl_wv_mhd_rp rp_type = gkyl_search_str_int_pair_by_str(mhd_rp_type, rp_str, WV_MHD_RP_ROE);

  const char *divb_str = glua_tbl_get_string(L, "divergenceConstraint", "none");
  enum gkyl_wv_mhd_div_constraint divb = gkyl_search_str_int_pair_by_str(
    mhd_divb_type, divb_str, GKYL_MHD_DIVB_NONE);

  double glm_ch = glua_tbl_get_number(L, "glmCh", 1.0);
  double glm_alpha = glua_tbl_get_number(L, "glmAlpha", 0.4);

  mhd_lw->magic = MOMENT_EQN_DEFAULT;
  mhd_lw->eqn = gkyl_wv_mhd_new( &(struct gkyl_wv_mhd_inp) {
      .gas_gamma = gas_gamma,
      .rp_type = rp_type,
      .divergence_constraint = divb,
      .glm_alpha = glm_alpha,
      .glm_ch = glm_ch
    }
  );

  // Create Lua userdata.
  struct wv_eqn_lw **l_mhd_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_mhd_lw = mhd_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_mhd_ctor[] = {
  { "new", eqn_mhd_lw_new },
  { 0, 0 }
};

/* ************************ */
/* Reactive Euler Equations */
/* ************************ */

// ReactiveEuler.new { gasGamma = 1.4, specificHeatCapacity = 2.5, energyOfFormation = 1.0, ignitionTemperature = 0.25, reactionRate = 250.0, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_reactive_euler_lw_new(lua_State *L)
{
  struct wv_eqn_lw *reactive_euler_lw = gkyl_malloc(sizeof(*reactive_euler_lw));

  double gas_gamma = glua_tbl_get_number(L, "gasGamma", 1.4);
  double specific_heat_capacity = glua_tbl_get_number(L, "specificHeatCapacity", 2.5);
  double energy_of_formation = glua_tbl_get_number(L, "energyOfFormation", 1.0);
  double ignition_temperature = glua_tbl_get_number(L, "ignitionTemperature", 0.25);
  double reaction_rate = glua_tbl_get_number(L, "reactionRate", 250.0);

  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_reactive_euler_rp rp_type = gkyl_search_str_int_pair_by_str(reactive_euler_rp_type, rp_str, WV_REACTIVE_EULER_RP_LAX);

  reactive_euler_lw->magic = MOMENT_EQN_DEFAULT;
  reactive_euler_lw->eqn = gkyl_wv_reactive_euler_inew( &(struct gkyl_wv_reactive_euler_inp) {
      .gas_gamma = gas_gamma,
      .specific_heat_capacity = specific_heat_capacity,
      .energy_of_formation = energy_of_formation,
      .ignition_temperature = ignition_temperature,
      .reaction_rate = reaction_rate,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct wv_eqn_lw **l_reactive_euler_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_reactive_euler_lw = reactive_euler_lw; // Point userdata to the equation object.
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Equation constructor.
static struct luaL_Reg eqn_reactive_euler_ctor[] = {
  { "new", eqn_reactive_euler_lw_new },
  { 0, 0 }
};

/* *********************** */
/* Euler Mixture Equations */
/* *********************** */

// EulerMixture.new { numComponents = 2, gasGamma = {1.4, 1.4}, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_euler_mixture_lw_new(lua_State *L)
{
  struct wv_eqn_lw *euler_mixture_lw = gkyl_malloc(sizeof(*euler_mixture_lw));

  int num_components = glua_tbl_get_integer(L, "numComponents", 2);

  double *gas_gamma_s = gkyl_malloc(sizeof(double[num_components]));
  with_lua_tbl_tbl(L, "gasGamma") {
    for (int i = 0; i < num_components; i++) {
      gas_gamma_s[i] = glua_tbl_iget_number(L, i + 1, 1.4);
    }
  }
  
  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_euler_mixture_rp rp_type = gkyl_search_str_int_pair_by_str(euler_mixture_rp_type, rp_str, WV_EULER_MIXTURE_RP_LAX);

  euler_mixture_lw->magic = MOMENT_EQN_DEFAULT;
  euler_mixture_lw->eqn = gkyl_wv_euler_mixture_inew( & (struct gkyl_wv_euler_mixture_inp) {
      .num_species = num_components,
      .gas_gamma_s = gas_gamma_s,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct wv_eqn_lw **l_euler_mixture_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_euler_mixture_lw = euler_mixture_lw; // Point userdata to the equation object.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static const luaL_Reg eqn_euler_mixture_ctor[] = {
  { "new", eqn_euler_mixture_lw_new },
  { 0, 0 }
};

/* ********************************** */
/* Isothermal Euler Mixture Equations */
/* ********************************** */

// IsoEulerMixture.new { numComponents = 2, vThermal = {1.0, 1.0}, rpType = "roe" }
// where rpType is one of "roe" or "lax".
static int
eqn_iso_euler_mixture_lw_new(lua_State *L)
{
  struct wv_eqn_lw *iso_euler_mixture_lw = gkyl_malloc(sizeof(*iso_euler_mixture_lw));

  int num_components = glua_tbl_get_integer(L, "numComponents", 2);

  double *vt_s = gkyl_malloc(sizeof(double[num_components]));
  with_lua_tbl_tbl(L, "vThermal") {
    for (int i = 0; i < num_components; i++) {
      vt_s[i] = glua_tbl_iget_number(L, i + 1, 1.0);
    }
  }
  
  const char *rp_str = glua_tbl_get_string(L, "rpType", "lax");
  enum gkyl_wv_iso_euler_mixture_rp rp_type = gkyl_search_str_int_pair_by_str(iso_euler_mixture_rp_type, rp_str, WV_ISO_EULER_MIXTURE_RP_LAX);

  iso_euler_mixture_lw->magic = MOMENT_EQN_DEFAULT;
  iso_euler_mixture_lw->eqn = gkyl_wv_iso_euler_mixture_inew( & (struct gkyl_wv_iso_euler_mixture_inp) {
      .num_species = num_components,
      .vt_s = vt_s,
      .rp_type = rp_type,
      .use_gpu = false
    }
  );

  // Create Lua userdata.
  struct wv_eqn_lw **l_iso_euler_mixture_lw = lua_newuserdata(L, sizeof(struct wv_eqn_lw*));
  *l_iso_euler_mixture_lw = iso_euler_mixture_lw; // Point userdata to the equation object.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);
  lua_setmetatable(L, -2);

  return 1;
}

// Equation constructor.
static const luaL_Reg eqn_iso_euler_mixture_ctor[] = {
  { "new", eqn_iso_euler_mixture_lw_new },
  { 0, 0 }
};

// Register and load all wave equation objects.
static void
eqn_openlibs(lua_State *L)
{
  luaL_newmetatable(L, MOMENT_WAVE_EQN_METATABLE_NM);

  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, wv_eqn_lw_gc);
  lua_settable(L, -3);

  luaL_register(L, "G0.Moments.Eq.Euler", eqn_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.IsoEuler", eqn_iso_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.SrEuler", eqn_sr_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.ColdFluid", eqn_coldfluid_ctor);
  luaL_register(L, "G0.Moments.Eq.TenMoment", eqn_tenmoment_ctor); 
  luaL_register(L, "G0.Moments.Eq.Mhd", eqn_mhd_ctor);
  luaL_register(L, "G0.Moments.Eq.ReactiveEuler", eqn_reactive_euler_ctor);
  luaL_register(L, "G0.Moments.Eq.EulerMixture", eqn_euler_mixture_ctor);
  luaL_register(L, "G0.Moments.Eq.IsoEulerMixture", eqn_iso_euler_mixture_ctor);
}

/* *************** */
/* Species methods */
/* *************** */

// Metatable name for species input struct.
#define MOMENT_SPECIES_METATABLE_NM "GkeyllZero.App.Moments.Species"

// Lua userdata object for constructing species input.
struct moment_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_moment_species mom_species; // Input struct to construct species.
  bool evolve; // Is this species evolved?

  struct lua_func_ctx init_ctx; // Lua registry reference to initialization function.

  bool has_app_accel; // Is there an applied acceleration function?
  struct lua_func_ctx app_accel_func_ctx; // Lua registry reference to applied acceleration function.

  bool has_nT_source; // Is there a temperature source function?
  struct lua_func_ctx nT_source_func_ctx; // Lua registry reference to temperature source function.
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

  if (has_eqn) {
    // We need to store a reference to equation object in a global
    // table as otherwise it gets garbage-collected while the
    // simulation is being set up.
    lua_getfield(L, -1, "equation");
    int eq_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    
    lua_getglobal(L, "__moment_ref_table");
    int eqtbl_len = glua_objlen(L);
    lua_pushinteger(L, eqtbl_len+1);    
    lua_rawgeti(L, LUA_REGISTRYINDEX, eq_ref);
    lua_rawset(L, -3);

    lua_pop(L, 1);
  }
  else {
    return luaL_error(L, "Species \"equation\" not specfied or incorrect type!");
  }

  const char *lim_str = glua_tbl_get_string(L, "limiter", "monotonized-centered");
  mom_species.limiter = gkyl_search_str_int_pair_by_str(wave_limiter, lim_str, GKYL_MONOTONIZED_CENTERED);

  const char *split_str = glua_tbl_get_string(L, "splitType", "qwave");
  mom_species.split_type = gkyl_search_str_int_pair_by_str(wave_split_type, split_str, GKYL_WAVE_QWAVE);

  bool evolve = mom_species.evolve = glua_tbl_get_bool(L, "evolve", true);
  mom_species.force_low_order_flux = glua_tbl_get_bool(L, "forceLowOrderFlux", false);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  else {
    return luaL_error(L, "Species must have an \"init\" function for initial conditions!");
  }

  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_species.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_species.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_species.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  bool has_app_accel = false;
  int app_accel_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "appliedAcceleration")) {
    has_app_accel = true;
    app_accel_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  mom_species.is_app_accel_static = glua_tbl_get_bool(L, "isAppliedAccelerationStatic", false);

  bool has_nT_source = false;
  int nT_source_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "nTSource")) {
    has_nT_source = true;
    nT_source_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  mom_species.nT_source_set_only_once = glua_tbl_get_bool(L, "nTSourceSetOnlyOnce", false);

  mom_species.has_reactivity = glua_tbl_get_bool(L, "hasReactivity", false);
  if (mom_species.has_reactivity) {
    mom_species.reactivity_gas_gamma = glua_tbl_get_number(L, "reactivityGasGamma", 1.4);
    mom_species.reactivity_specific_heat_capacity = glua_tbl_get_number(L, "reactivitySpecificHeatCapacity", 2.5);
    mom_species.reactivity_energy_of_formation = glua_tbl_get_number(L, "reactivityEnergyOfFormation", 1.0);
    mom_species.reactivity_ignition_temperature = glua_tbl_get_number(L, "reactivityIgnitionTemperature", 0.25);
    mom_species.reactivity_reaction_rate = glua_tbl_get_number(L, "reactivityReactionRate", 250.0);
  }

  mom_species.has_volume_sources = glua_tbl_get_bool(L, "hasVolumeSources", false);
  if (mom_species.has_volume_sources) {
    mom_species.volume_gas_gamma = glua_tbl_get_number(L, "volumeGasGamma", 5.0 / 3.0);
    mom_species.volume_U0 = glua_tbl_get_number(L, "volumeU0", 1.0);
    mom_species.volume_R0 = glua_tbl_get_number(L, "volumeR0", 1.0);
  }

  struct moment_species_lw *moms_lw = lua_newuserdata(L, sizeof(*moms_lw));
  moms_lw->magic = MOMENT_SPECIES_DEFAULT;
  moms_lw->evolve = evolve;
  moms_lw->mom_species = mom_species;

  moms_lw->init_ctx = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = mom_species.equation->num_equations,
    .L = L,
  };

  moms_lw->has_app_accel = has_app_accel;
  moms_lw->app_accel_func_ctx = (struct lua_func_ctx) {
    .func_ref = app_accel_ref,
    .ndim = 0, // This will be set later.
    .nret = GKYL_MOM_APP_NUM_APPLIED_ACCELERATION,
    .L = L,
  };

  moms_lw->has_nT_source = has_nT_source;
  moms_lw->nT_source_func_ctx = (struct lua_func_ctx) {
    .func_ref = nT_source_ref,
    .ndim = 0, // This will be set later.
    .nret = GKYL_MOM_APP_NUM_NT_SOURCE,
    .L = L,
  };    
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg mom_species_ctor[] = {
  { "new", moment_species_lw_new },
  { 0, 0 }
};

/* ************* */
/* Field methods */
/* ************* */

// Metatable name for field input struct.
#define MOMENT_FIELD_METATABLE_NM "GkeyllZero.App.Moments.Field"

// Lua userdata object for constructing field input.
struct moment_field_lw {
  int magic; // This must be first element in the struct.

  bool evolve; // Is this field evolved?  
  struct gkyl_moment_field mom_field; // Input struct to construct field.
  
  struct lua_func_ctx init_ctx; // Lua registry reference to initialization function.

  bool has_app_current; // Is there an applied current function?
  struct lua_func_ctx app_current_ctx; // Lua registry reference to applied current function.

  bool has_ext_em; // Is there an external field function?
  struct lua_func_ctx ext_em_ctx; // Lua registry reference to external EM function.
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
  mom_field.limiter = gkyl_search_str_int_pair_by_str(wave_limiter, lim_str, GKYL_MONOTONIZED_CENTERED);
  
  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init")) {
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  else {
    return luaL_error(L, "Field must have an \"init\" function for initial conditions!");
  }
  
  with_lua_tbl_tbl(L, "bcx") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_field.bcx[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_field.bcy[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    int nbc = glua_objlen(L);

    for (int i = 0; i < (nbc > 2 ? 2 : nbc); i++) {
      mom_field.bcz[i] = glua_tbl_iget_integer(L, i + 1, 0);
    }
  }

  bool has_app_current = false;
  int app_current_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "appliedCurrent")) {
    has_app_current = true;
    app_current_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  mom_field.t_ramp_curr = glua_tbl_get_number(L, "currentRampTime", 0.0);

  bool has_ext_em = false;
  int ext_em_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "externalEm")) {
    has_ext_em = true;
    ext_em_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  mom_field.t_ramp_E = glua_tbl_get_number(L, "externalEmRampTime", 0.0);
  mom_field.is_ext_em_static = glua_tbl_get_bool(L, "isExternalEmStatic", false);

  mom_field.use_explicit_em_coupling = glua_tbl_get_bool(L, "useExplicitEmCoupling", false);

  struct moment_field_lw *momf_lw = lua_newuserdata(L, sizeof(*momf_lw));

  momf_lw->magic = MOMENT_FIELD_DEFAULT;
  momf_lw->evolve = evolve;
  momf_lw->mom_field = mom_field;

  momf_lw->init_ctx = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // This will be set later.
    .nret = 6,
    .L = L,
  };  

  momf_lw->has_app_current = has_app_current;
  momf_lw->app_current_ctx = (struct lua_func_ctx) {
    .func_ref = app_current_ref,
    .ndim = 0, // This will be set later.
    .nret = GKYL_MOM_APP_NUM_APPLIED_CURRENT,
    .L = L,
  };

  momf_lw->has_ext_em = has_ext_em;
  momf_lw->ext_em_ctx = (struct lua_func_ctx) {
    .func_ref = ext_em_ref,
    .ndim = 0, // This will be set later.
    .nret = GKYL_MOM_APP_NUM_EXT_EM,
    .L = L,
  };  
  
  // Set metatable.
  luaL_getmetatable(L, MOMENT_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg mom_field_ctor[] = {
  { "new", moment_field_lw_new },
  { 0, 0 }
};

/* *********** */
/* App methods */
/* *********** */

// Metatable name for top-level Moment App.
#define MOMENT_APP_METATABLE_NM "GkeyllZero.App.Moment"

// Lua userdata object for holding Moment app and simulation parameters.
struct moment_app_lw {
  gkyl_moment_app *app; // Moment app object.
  
  struct lua_func_ctx mapc2p_ctx; // Function context for mapc2p.

  struct lua_func_ctx species_init_ctx[GKYL_MAX_SPECIES]; // Function context for species initial conditions.
  struct lua_func_ctx species_app_accel_func_ctx[GKYL_MAX_SPECIES]; // Function context for applied acceleration.
  struct lua_func_ctx species_nT_source_func_ctx[GKYL_MAX_SPECIES]; // Function context for temperature sources.

  struct lua_func_ctx field_init_ctx; // Function context for field initial conditions.
  struct lua_func_ctx field_app_current_ctx; // Function context for applied current.
  struct lua_func_ctx field_ext_em_ctx; // Function context for external EM field.
  
  double t_start, t_end; // Start and end times of simulation.
  int num_frames; // Number of data frames to write.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Find index of species in table given its name.
static int
app_find_species(const gkyl_moment_app *app, const char *nm)
{
  for (int i = 0; i < app->num_species; i++) {
    if (strcmp(app->species[i].name, nm) == 0) {
      return i;
    }
  }
  
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
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct moment_species_lw *vms = lua_touserdata(L, TVAL);
      if (vms->magic == MOMENT_SPECIES_DEFAULT) {
        
        vms->init_ctx.ndim = cdim;
        vms->app_accel_func_ctx.ndim = cdim;
        vms->nT_source_func_ctx.ndim = cdim;
        
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

// Create top-level App object.
static int
mom_app_new(lua_State *L)
{
  struct moment_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-vlasov" is used.
  const char *sim_name = "g0-vlasov";
  with_lua_global(L, "GKYL_OUT_PREFIX") {
    if (lua_isstring(L, -1)) {
      sim_name = lua_tostring(L, -1);
    }
  }
  
  // Initialize app using table inputs (table is on top of stack).

  app_lw->t_start = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->t_end = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->num_frames = glua_tbl_get_integer(L, "nFrame", 1);
  app_lw->dt_failure_tol = glua_tbl_get_number(L, "dtFailureTol", 1.0e-4);
  app_lw->num_failures_max = glua_tbl_get_integer(L, "numFailuresMax", 20);

  struct gkyl_moment mom = { }; // Input table for app.

  strcpy(mom.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    mom.ndim = cdim = glua_objlen(L);

    for (int d = 0; d < cdim; d++) {
      mom.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  int cuts[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
  
  with_lua_tbl_tbl(L, "decompCuts") {
    int ncuts = glua_objlen(L);

    for (int d = 0; d < ncuts; d++) {
      cuts[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }  

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < cdim; d++) {
      mom.lower[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < cdim; d++) {
      mom.upper[d] = glua_tbl_iget_number(L, d + 1, 0);
    }
  }

  mom.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);

  mom.scheme_type = glua_tbl_get_integer(L, "schemeType", 0);

  mom.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      mom.num_periodic_dir = glua_objlen(L);

      for (int d = 0; d < mom.num_periodic_dir; d++) {
        // Indices are off by 1 between Lua and C.
        mom.periodic_dirs[d] = glua_tbl_iget_integer(L, d + 1, 0) - 1;
      }
    }
  }

  // mapc2p function.
  mom.c2p_ctx = 0;
  mom.mapc2p = 0;
  bool has_mapc2p = false;
  int mapc2p_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "mapc2p")) {
    has_mapc2p = true;
    mapc2p_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  if (has_mapc2p) {
    app_lw->mapc2p_ctx = (struct lua_func_ctx) {
      .func_ref = mapc2p_ref,
      .ndim = cdim,
      .nret = cdim,
      .L = L,
    };
    mom.mapc2p = gkyl_lw_eval_cb;
    mom.c2p_ctx = &app_lw->mapc2p_ctx;
  }

  struct moment_species_lw *species[GKYL_MAX_SPECIES];

  // Set all species input.
  mom.num_species = get_species_inp(L, cdim, species);
  for (int s = 0; s < mom.num_species; s++) {
    mom.species[s] = species[s]->mom_species;
    
    app_lw->species_init_ctx[s] = species[s]->init_ctx;
    mom.species[s].init = gkyl_lw_eval_cb;
    mom.species[s].ctx = &app_lw->species_init_ctx[s];

    if (species[s]->has_app_accel) {
      app_lw->species_app_accel_func_ctx[s] = species[s]->app_accel_func_ctx;
      mom.species[s].app_accel_func = gkyl_lw_eval_cb;
      mom.species[s].app_accel_ctx = &app_lw->species_app_accel_func_ctx[s];
    }

    if (species[s]->has_nT_source) {
      app_lw->species_nT_source_func_ctx[s] = species[s]->nT_source_func_ctx;
      mom.species[s].nT_source_func = gkyl_lw_eval_cb;
      mom.species[s].nT_source_ctx = &app_lw->species_nT_source_func_ctx[s];
    }
  }

  // Set field input.
  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct moment_field_lw *momf = lua_touserdata(L, -1);
      
      if (momf->magic == MOMENT_FIELD_DEFAULT) {
        momf->init_ctx.ndim = cdim;
        
        mom.field = momf->mom_field;

        app_lw->field_init_ctx = momf->init_ctx;
        mom.field.init = gkyl_lw_eval_cb;
        mom.field.ctx = &app_lw->field_init_ctx;

        if (momf->has_app_current) {
          momf->app_current_ctx.ndim = cdim;

          app_lw->field_app_current_ctx = momf->app_current_ctx;
          mom.field.app_current_func = gkyl_lw_eval_cb;
          mom.field.app_current_ctx = &app_lw->field_app_current_ctx;
        }

        if (momf->has_ext_em) {
          momf->ext_em_ctx.ndim = cdim;

          app_lw->field_ext_em_ctx = momf->ext_em_ctx;
          mom.field.ext_em_func = gkyl_lw_eval_cb;
          mom.field.ext_em_ctx = &app_lw->field_app_current_ctx;
        }
      }
    }
  }

  // Create parallelism.
  struct gkyl_comm *comm = 0;
  bool has_mpi = false;

  for (int d = 0; d < cdim; d++) {
    mom.parallelism.cuts[d] = cuts[d]; 
  }

#ifdef GKYL_HAVE_MPI
  with_lua_global(L, "GKYL_MPI_COMM") {
    if (lua_islightuserdata(L, -1)) {
      has_mpi = true;
      MPI_Comm mpi_comm = lua_touserdata(L, -1);
      comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
          .mpi_comm = mpi_comm,
          .sync_corners = true
        }
      );

    }
  }
#endif

  if (!has_mpi) {
    // If there is no proper MPI_Comm specifed, then assume we are a
    // serial sim.
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .sync_corners = true
      }
    );
  }
  mom.parallelism.comm = comm;

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int tot_cuts = 1; for (int d=0; d<cdim; ++d) tot_cuts *= cuts[d];

  if (tot_cuts != comm_sz) {
    printf("tot_cuts = %d (%d)\n", tot_cuts, comm_sz);
    luaL_error(L, "Number of ranks and cuts do not match!");
  }
  
  app_lw->app = gkyl_moment_app_new(&mom); // Create the Moment app.

  gkyl_comm_release(comm);
  
  // Create Lua userdata.
  struct moment_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct moment_app_lw*));
  *l_app_lw = app_lw; // Point it to the Lua app pointer.

  // Set metatable.
  luaL_getmetatable(L, MOMENT_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Return number of species () -> int.
static int
mom_app_num_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  lua_pushnumber(L, app_lw->app->num_species);

  return 1;
}

// Return number of species (int) -> string.
static int
mom_app_species_name(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int idx = luaL_checkinteger(L, 2) - 1;
  lua_pushstring(L, app_lw->app->species[idx].name);

  return 1;
}

// Compute maximum time-step () -> double.
static int
mom_app_max_dt(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double maxdt = gkyl_moment_app_max_dt(app_lw->app);

  lua_pushnumber(L, maxdt);  
  return 1;
}

// Apply initial conditions. (time) -> bool.
static int
mom_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_moment_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool.
static int
mom_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_moment_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (name, time) -> bool.
static int
mom_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *sp_name = luaL_checkstring(L, 2);
  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in apply_ic_species!", sp_name);
  }
  
  double t0 = luaL_optnumber(L, 3, app_lw->t_start);
  gkyl_moment_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// The status table for reads has the following structure:
// {
//   io_status = true or false.
//   io_status_str = "success", "bad-version", "fopen-failed", "fread-failed", "data-mismatch"
//   frame = frame_num (integer)
//   stime = time in file read
// }

static const struct gkyl_str_int_pair rio_status[] = {
  { "success", GKYL_ARRAY_RIO_SUCCESS },
  { "bad-version", GKYL_ARRAY_RIO_BAD_VERSION },
  { "fopen-failed", GKYL_ARRAY_RIO_FOPEN_FAILED },
  { "fread-failed", GKYL_ARRAY_RIO_FREAD_FAILED },
  { "data-mismatch", GKYL_ARRAY_RIO_DATA_MISMATCH },
  { 0, 0 }
};

// Pushes table with status on stack. Table is left on stack.
static void
push_restart_status_table(lua_State *L, struct gkyl_app_restart_status status)
{
  lua_newtable(L);

  lua_pushstring(L, "io_status");
  lua_pushboolean(L, status.io_status == GKYL_ARRAY_RIO_SUCCESS ? true : false);
  lua_rawset(L, -3);

  lua_pushstring(L, "io_status_str");
  lua_pushstring(L, gkyl_search_str_int_pair_by_int(rio_status, status.io_status, "success"));
  lua_rawset(L, -3);

  lua_pushstring(L, "frame");
  lua_pushinteger(L, status.frame);
  lua_rawset(L, -3);

  lua_pushstring(L, "stime");
  lua_pushnumber(L, status.stime);
  lua_rawset(L, -3);  
}

// Read field from file. (file-name) -> status table. See above.
static int
mom_app_from_file_field(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *fname = luaL_checkstring(L, 2);
  struct gkyl_app_restart_status status = gkyl_moment_app_from_file_field(app_lw->app, fname);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (file-name, species-name) -> status table. See above.
static int
mom_app_from_file_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *fname = luaL_checkstring(L, 2);
  const char *sp_name = luaL_checkstring(L, 3);

  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in from_file_species!", sp_name);
  }
  
  struct gkyl_app_restart_status status = gkyl_moment_app_from_file_species(app_lw->app, sidx, fname);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (frame) -> status table. See above.
static int
mom_app_from_frame_field(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int frame = luaL_checkinteger(L, 2);
  struct gkyl_app_restart_status status = gkyl_moment_app_from_frame_field(app_lw->app, frame);

  push_restart_status_table(L, status);
  
  return 1;
}

// Read field from file. (frame, species-name) -> status table. See above.
static int
mom_app_from_frame_species(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int frame = luaL_checkinteger(L, 2);
  const char *sp_name = luaL_checkstring(L, 3);

  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in from_frame_species!", sp_name);
  }
  
  struct gkyl_app_restart_status status = gkyl_moment_app_from_frame_species(app_lw->app, sidx, frame);

  push_restart_status_table(L, status);
  
  return 1;
}

// Compute integrated moments. (tm) -> bool.
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

// get number of field-energy diagnostics stored () -> int.
static int
mom_app_field_energy_ndiag(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  lua_pushinteger(L, gkyl_moment_app_field_energy_ndiag(app_lw->app));
  return 1;
}

// Return the field energy as a table () -> table.
static int
mom_app_get_field_energy(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  int nvals = gkyl_moment_app_field_energy_ndiag(app_lw->app);
  double vals[nvals];
  gkyl_moment_app_get_field_energy(app_lw->app, vals);

  lua_createtable(L, nvals, 0);

  for (int i = 0; i < nvals; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, vals[i]);
    lua_rawset(L, -3);
  }  

  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool.
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

// Write solution (field and species) to file (time, frame) -> bool.
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

// Write field to file (time, frame) -> bool.
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

// Write species solution to file (name, time, frame) -> bool.
static int
mom_app_write_species(lua_State *L)
{
  bool status = true;

  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  const char *sp_name = luaL_checkstring(L, 2);
  int sidx = app_find_species(app_lw->app, sp_name);
  if (sidx < 0) {
    return luaL_error(L, "Incorrect species name '%s' in write_species!", sp_name);
  }

  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_moment_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool.
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

// Write integrated field energy to file () -> bool.
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

// Write simulation statistics to JSON. () -> bool.
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

// Write data from simulation to file.
static void
write_data(struct gkyl_tm_trigger* iot, gkyl_moment_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_moment_app_write(app, t_curr, frame);
  }
}

// Update status table has the following structure:
// {
//    success = true or false
//    dt_actual = actual time-step taken
//    dt_suggested = suggested stable time-step
// }

// Update the solution by a suggested time-step. (dt) -> update status
// table. See above. For details see the C API doc for this function.
static int
mom_app_update(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  double dt = luaL_checknumber(L, 2);
  struct gkyl_update_status status = gkyl_moment_update(app_lw->app, dt);

  // Return status table on stack.
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

// Run simulation. (num_steps) -> bool. num_steps is optional.
static int
mom_app_run(lua_State *L)
{
  bool ret_status = true;

  // Create app object.
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;
  struct gkyl_moment_app *app = app_lw->app;

  // Initial and final simulation times.
  double t_curr = app_lw->t_start, t_end = app_lw->t_end;
  long num_steps = luaL_optinteger(L, 2, INT_MAX);

  // Create trigger for IO.
  int num_frames = app_lw->num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_moment_app_apply_ic(app, t_curr);
  gkyl_moment_app_calc_integrated_mom(app, t_curr);
  gkyl_moment_app_calc_field_energy(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  // Compute estimate of maximum stable time-step.
  double dt = gkyl_moment_app_max_dt(app);

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = app_lw->dt_failure_tol;
  int num_failures = 0, num_failures_max = app_lw->num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    gkyl_moment_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    gkyl_moment_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_moment_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    gkyl_moment_app_calc_integrated_mom(app, t_curr);
    gkyl_moment_app_calc_field_energy(app, t_curr);

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_moment_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_moment_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_moment_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_moment_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_moment_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  gkyl_moment_app_stat_write(app);

  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  gkyl_moment_app_cout(app, stdout, "\n");
  gkyl_moment_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_moment_app_cout(app, stdout, "Number of failed time-steps %ld\n", stat.nfail);
  gkyl_moment_app_cout(app, stdout, "Species updates took %g secs\n", stat.species_tm);
  gkyl_moment_app_cout(app, stdout, "Field updates took %g secs\n", stat.field_tm);
  gkyl_moment_app_cout(app, stdout, "Source updates took %g secs\n", stat.sources_tm);
  gkyl_moment_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  lua_pushboolean(L, ret_status);
  return 1;
}

// Return ghost cell as table of 2 elements.
static int
mom_app_nghost(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;
  struct gkyl_moment_app *app = app_lw->app;

  int nghost[3] = { 0 };
  gkyl_moment_app_nghost(app, nghost);
  
  lua_createtable(L, 3, 0);

  for (int i = 0; i < 3; i++) {
    lua_pushinteger(L, i + 1);
    lua_pushinteger(L, nghost[i]);
    lua_rawset(L, -3);
  }
  
  return 1;
}

// Clean up memory allocated for simulation.
static int
mom_app_gc(lua_State *L)
{
  struct moment_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, MOMENT_APP_METATABLE_NM);
  struct moment_app_lw *app_lw = *l_app_lw;

  gkyl_moment_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor.
static struct luaL_Reg mom_app_ctor[] = {
  { "new",  mom_app_new },
  { 0, 0 }
};

// App methods
static struct luaL_Reg mom_app_funcs[] = {
  { "num_species", mom_app_num_species },
  { "species_name", mom_app_species_name },

  { "max_dt", mom_app_max_dt },  

  { "apply_ic", mom_app_apply_ic },
  { "apply_ic_field", mom_app_apply_ic_field },
  { "apply_ic_species", mom_app_apply_ic_species },
  
  { "from_file_field", mom_app_from_file_field },
  { "from_file_species", mom_app_from_file_species },

  { "from_frame_field", mom_app_from_frame_field },
  { "from_frame_species", mom_app_from_frame_species },

  { "write", mom_app_write },
  { "write_field", mom_app_write_field },
  { "write_species", mom_app_write_species },
  { "write_field_energy", mom_app_write_field_energy },
  { "write_integrated_mom", mom_app_write_integrated_mom },
  { "stat_write", mom_app_stat_write },

  { "calc_field_energy", mom_app_calc_field_energy },
  { "calc_integrated_mom", mom_app_calc_integrated_mom },

  { "update", mom_app_update },
  { "run", mom_app_run },

  // Some low-level functions, typically not used by ordinary users.
  { "nghost", mom_app_nghost },
  { "field_energy_ndiag", mom_app_field_energy_ndiag },
  { "get_field_energy", mom_app_get_field_energy },
  
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App.
  do {
    luaL_newmetatable(L, MOMENT_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, mom_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, mom_app_funcs);
    
    luaL_register(L, "G0.Moments.App", mom_app_ctor);
    
  }
  while (0);

  // Register Species input struct.
  do {
    luaL_newmetatable(L, MOMENT_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Moments.Species", mom_species_ctor);
  }
  while (0);

  // Register Field input struct.
  do {
    luaL_newmetatable(L, MOMENT_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Moments.Field", mom_field_ctor);
  }
  while (0);

  // Add globals and other parameters.
  do {
    // Table to store references so the objects do not get garbage-collected.
    lua_newtable(L);
    lua_setglobal(L, "__moment_ref_table");
  }
  while (0);
}

void
gkyl_moment_lw_openlibs(lua_State *L)
{
  eqn_openlibs(L);
  app_openlibs(L);
}

#endif
