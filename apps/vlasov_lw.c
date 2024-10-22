#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_lw.h>
#include <gkyl_vlasov_priv.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Magic IDs for use in distinguishing various species and field types
enum vlasov_magic_ids {
  VLASOV_SPECIES_DEFAULT = 100, // non-relativistic kinetic species
  VLASOV_FIELD_DEFAULT, // Maxwell equations
};

// Vlasov projection type -> enum map.
static const struct gkyl_str_int_pair projection_type[] = {
  { "func", GKYL_PROJ_FUNC },
  { "maxwellianPrimitive", GKYL_PROJ_MAXWELLIAN_PRIM },
  { "maxwellianLab", GKYL_PROJ_MAXWELLIAN_LAB },
  { "biMaxwellian", GKYL_PROJ_BIMAXWELLIAN },
  { "LTE", GKYL_PROJ_VLASOV_LTE },
  { 0, 0 }
};

// Vlasov model type -> enum map.
static const struct gkyl_str_int_pair model_type[] = {
  { "default", GKYL_MODEL_DEFAULT },
  { "SR", GKYL_MODEL_SR },
  { "generalGeometry", GKYL_MODEL_GEN_GEO },
  { "canonicalPB", GKYL_MODEL_CANONICAL_PB },
  { 0, 0 }
};

// Vlasov collision type -> enum map.
static const struct gkyl_str_int_pair collision_type[] = {
  { "none", GKYL_NO_COLLISIONS },
  { "BGK", GKYL_BGK_COLLISIONS },
  { "LBO", GKYL_LBO_COLLISIONS },
  { "FPO", GKYL_FPO_COLLISIONS },
  { 0, 0 }
};

/* *****************/
/* Species methods */
/* *****************/

// Metatable name for species input struct
#define VLASOV_SPECIES_METATABLE_NM "GkeyllZero.App.Vlasov.Species"

// Lua userdata object for constructing species input
struct vlasov_species_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_vlasov_species vm_species; // input struct to construct species
  int vdim; // velocity dimensions
  bool evolve; // is this species evolved?

  int num_init; // Number of projection objects.
  enum gkyl_projection_id proj_id[GKYL_MAX_PROJ]; // Projection type.

  bool has_init_func[GKYL_MAX_PROJ]; // Is there an initialization function?
  struct lua_func_ctx init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to initialization function.

  bool has_density_init_func[GKYL_MAX_PROJ]; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to density initialization function.

  bool has_V_drift_init_func[GKYL_MAX_PROJ]; // Is there a drift velocity initialiation function?
  struct lua_func_ctx V_drift_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to drift velocity initialization function.

  bool has_temp_init_func[GKYL_MAX_PROJ]; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to temperature initialization function.

  bool correct_all_moms[GKYL_MAX_PROJ]; // Are we correcting all moments in projections, or only density?
  double iter_eps[GKYL_MAX_PROJ]; // Error tolerance for moment fixes in projections (density is always exact).
  int max_iter[GKYL_MAX_PROJ]; // Maximum number of iterations for moment fixes in projections.
  bool use_last_converged[GKYL_MAX_PROJ]; // Use last iteration value in projections regardless of convergence?

  enum gkyl_collision_id collision_id; // Collision type.
  
  bool has_self_nu_func; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ref; // Lua registry reference to self-collision frequency function.

  bool collision_correct_all_moms; // Are we correcting all moments in collisions, or only density?
};

static int
vlasov_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_vlasov_species vm_species = { };

  const char *model_str = glua_tbl_get_string(L, "modelID", "default");
  vm_species.model_id = gkyl_search_str_int_pair_by_str(model_type, model_str, GKYL_MODEL_DEFAULT);
  
  vm_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  vm_species.mass = glua_tbl_get_number(L, "mass", 1.0);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);
    for (int d=0; d<vdim; ++d)
      vm_species.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<vdim; ++d)
      vm_species.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<vdim; ++d)
      vm_species.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  with_lua_tbl_tbl(L, "diagnostics") {
    int num_diag_moments = glua_objlen(L);

    int n = 0;
    for (int i=0; i<num_diag_moments; ++i) {
      const char *mom = glua_tbl_iget_string(L, i+1, "");
      if (is_moment_name_valid(mom))
        strcpy(vm_species.diag_moments[n++], mom);
    }
    vm_species.num_diag_moments = n;
  }
  
  enum gkyl_projection_id proj_id[GKYL_MAX_PROJ];

  bool has_init_func[GKYL_MAX_PROJ];
  int init_func_ref[GKYL_MAX_PROJ];

  bool has_density_init_func[GKYL_MAX_PROJ];
  int density_init_func_ref[GKYL_MAX_PROJ];

  bool has_V_drift_init_func[GKYL_MAX_PROJ];
  int V_drift_init_func_ref[GKYL_MAX_PROJ];

  bool has_temp_init_func[GKYL_MAX_PROJ];
  int temp_init_func_ref[GKYL_MAX_PROJ];
  
  bool correct_all_moms[GKYL_MAX_PROJ];
  double iter_eps[GKYL_MAX_PROJ];
  int max_iter[GKYL_MAX_PROJ];
  bool use_last_converged[GKYL_MAX_PROJ];

  int num_init = glua_tbl_get_integer(L, "numInit", 0);

  with_lua_tbl_tbl(L, "projections") {
    for (int i = 0; i < num_init; i++) {
      if (glua_tbl_iget_tbl(L, i + 1)) {
        const char *projection_str = glua_tbl_get_string(L, "projectionID", "func");
        proj_id[i] = gkyl_search_str_int_pair_by_str(projection_type, projection_str, GKYL_PROJ_FUNC);

        init_func_ref[i] = LUA_NOREF;
        has_init_func[i] = false;
        if (glua_tbl_get_func(L, "init")) {
          init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_init_func[i] = true;
        }

        density_init_func_ref[i] = LUA_NOREF;
        has_density_init_func[i] = false;
        if (glua_tbl_get_func(L, "densityInit")) {
          density_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_density_init_func[i] = true;
        }

        V_drift_init_func_ref[i] = LUA_NOREF;
        has_V_drift_init_func[i] = false;
        if (glua_tbl_get_func(L, "driftVelocityInit")) {
          V_drift_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_V_drift_init_func[i] = true;
        }

        temp_init_func_ref[i] = LUA_NOREF;
        has_temp_init_func[i] = false;
        if (glua_tbl_get_func(L, "temperatureInit")) {
          temp_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
          has_temp_init_func[i] = true;
        }

        correct_all_moms[i] = glua_tbl_get_bool(L, "correctAllMoments", true);
        iter_eps[i] = glua_tbl_get_number(L, "iterationEpsilon", pow(10.0, -12.0));
        max_iter[i] = glua_tbl_get_integer(L, "maxIterations", 100);
        use_last_converged[i] = glua_tbl_get_bool(L, "useLastConverged", true);

        lua_pop(L, 1);
      }
    }
  }

  enum gkyl_collision_id collision_id = GKYL_NO_COLLISIONS;

  bool has_self_nu_func = false;
  int self_nu_func_ref = LUA_NOREF;

  bool collision_correct_all_moms = false;

  with_lua_tbl_tbl(L, "collisions") {
    const char *collision_str = glua_tbl_get_string(L, "collisionID", "none");
    collision_id = gkyl_search_str_int_pair_by_str(collision_type, collision_str, GKYL_NO_COLLISIONS);

    if (glua_tbl_get_func(L, "selfNu")) {
      self_nu_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_self_nu_func = true;
    }

    collision_correct_all_moms = glua_tbl_get_bool(L, "correctAllMoments", true);
  }
  
  struct vlasov_species_lw *vms_lw = lua_newuserdata(L, sizeof(*vms_lw));
  vms_lw->magic = VLASOV_SPECIES_DEFAULT;
  vms_lw->vdim = vdim;
  vms_lw->evolve = evolve;
  vms_lw->vm_species = vm_species;

  vms_lw->num_init = num_init;
  for (int i = 0; i < num_init; i++) {
    vms_lw->proj_id[i] = proj_id[i];

    vms_lw->has_init_func[i] = has_init_func[i];
    vms_lw->init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->has_density_init_func[i] = has_density_init_func[i];
    vms_lw->density_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = density_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->has_V_drift_init_func[i] = has_V_drift_init_func[i];
    vms_lw->V_drift_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = V_drift_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->has_temp_init_func[i] = has_temp_init_func[i];
    vms_lw->temp_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = temp_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    vms_lw->correct_all_moms[i] = correct_all_moms[i];
    vms_lw->iter_eps[i] = iter_eps[i];
    vms_lw->max_iter[i] = max_iter[i];
    vms_lw->use_last_converged[i] = use_last_converged[i];
  }

  vms_lw->collision_id = collision_id;

  vms_lw->has_self_nu_func = has_self_nu_func;
  vms_lw->self_nu_func_ref = (struct lua_func_ctx) {
    .func_ref = self_nu_func_ref,
    .ndim = 0,
    .nret = 1,
    .L = L,
  };

  vms_lw->collision_correct_all_moms = collision_correct_all_moms;
  
  // set metatable
  luaL_getmetatable(L, VLASOV_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg vm_species_ctor[] = {
  {"new", vlasov_species_lw_new},
  {0, 0}
};

/* *****************/
/* Field methods */
/* *****************/

// Metatable name for field input struct
#define VLASOV_FIELD_METATABLE_NM "GkeyllZero.App.Vlasov.Field"

// Lua userdata object for constructing field input
struct vlasov_field_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_vlasov_field vm_field; // input struct to construct field
  bool evolve; // is this field evolved?
  struct lua_func_ctx init_ref; // Lua registery reference to initilization function
};

static int
vlasov_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_vlasov_field vm_field = { };

  vm_field.field_id = GKYL_FIELD_E_B;  
  
  vm_field.epsilon0 = glua_tbl_get_number(L, "epsilon0", 1.0);
  vm_field.mu0 = glua_tbl_get_number(L, "mu0", 1.0);
  vm_field.elcErrorSpeedFactor = glua_tbl_get_number(L, "elcErrorSpeedFactor", 0.0);
  vm_field.mgnErrorSpeedFactor = glua_tbl_get_number(L, "mgnErrorSpeedFactor", 0.0);

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  int init_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "init"))
    init_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Field must have an \"init\" function for initial conditions!");

  struct vlasov_field_lw *vmf_lw = lua_newuserdata(L, sizeof(*vmf_lw));

  vmf_lw->magic = VLASOV_FIELD_DEFAULT;
  vmf_lw->evolve = evolve;
  vmf_lw->vm_field = vm_field;
  
  vmf_lw->init_ref = (struct lua_func_ctx) {
    .func_ref = init_ref,
    .ndim = 0, // this will be set later
    .nret = 6,
    .L = L,
  };  
  
  // set metatable
  luaL_getmetatable(L, VLASOV_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg vm_field_ctor[] = {
  { "new",  vlasov_field_lw_new },
  { 0, 0 }
};

/* *************/
/* App methods */
/* *************/

// Metatable name for top-level Vlasov App
#define VLASOV_APP_METATABLE_NM "GkeyllZero.App.Vlasov"

// Lua userdata object for holding Vlasov app and run parameters
struct vlasov_app_lw {
  gkyl_vlasov_app *app; // Vlasov app object

  int num_init[GKYL_MAX_SPECIES]; // Number of projection objects.
  enum gkyl_projection_id proj_id[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Projection type.

  bool has_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there an initialization function?
  struct lua_func_ctx init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for initialization function.

  bool has_density_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for density initialization function.

  bool has_V_drift_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a drift velocity initialization function?
  struct lua_func_ctx V_drift_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for drift velocity initialziation function.
  
  bool has_temp_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for temperature initialization function.

  bool correct_all_moms[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Are we correcting all moments in projections, or only density?
  double iter_eps[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Error tolerance for moment fixes in projections (density is always exact).
  int max_iter[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Maximum number of iterations for moment fixes in projections.
  bool use_last_converged[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Use last iteration value in projections regardless of convergence?

  enum gkyl_collision_id collision_id[GKYL_MAX_SPECIES]; // Collision type.

  bool has_self_nu_func[GKYL_MAX_SPECIES]; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ctx[GKYL_MAX_SPECIES]; // Context for self-collision frequency function.

  bool collision_correct_all_moms[GKYL_MAX_SPECIES]; // Are we correct all moments in collisions, or only density?

  struct lua_func_ctx field_func_ctx; // function context for field
  
  double tstart, tend; // start and end times of simulation
  int nframe; // number of data frames to write
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct vlasov_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1};
  
  int curr = 0;
  lua_pushnil(L); // initial key is nil
  while (lua_next(L, TKEY) != 0) {
    // key at TKEY and value at TVAL
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct vlasov_species_lw *vms = lua_touserdata(L, TVAL);
      if (vms->magic == VLASOV_SPECIES_DEFAULT) {
        
        for (int i =0; i < vms->num_init; i++) {
          if (vms->has_init_func[i]) {
            vms->init_func_ref[i].ndim = cdim + vms->vdim;
          }

          if (vms->has_density_init_func[i]) {
            vms->density_init_func_ref[i].ndim = cdim + vms->vdim;
          }
          
          if (vms->has_V_drift_init_func[i]) {
            vms->V_drift_init_func_ref[i].ndim = cdim + vms->vdim;
          }

          if (vms->has_temp_init_func[i]) {
            vms->temp_init_func_ref[i].ndim = cdim + vms->vdim;
          }
        }

        if (vms->has_self_nu_func) {
          vms->self_nu_func_ref.ndim = cdim + vms->vdim;
        }
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(vms->vm_species.name, key);
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
vm_app_new(lua_State *L)
{
  struct vlasov_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

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

  struct gkyl_vm vm = { }; // input table for app

  strcpy(vm.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    vm.cdim = cdim = glua_objlen(L);
    for (int d=0; d<cdim; ++d)
      vm.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
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
      vm.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<cdim; ++d)
      vm.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  vm.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);
  vm.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  vm.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  vm.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      vm.num_periodic_dir = glua_objlen(L);
      for (int d=0; d<vm.num_periodic_dir; ++d)
        // indexes are off by 1 between Lua and C
        vm.periodic_dirs[d] = glua_tbl_iget_integer(L, d+1, 0)-1;
    }
  }

  struct vlasov_species_lw *species[GKYL_MAX_SPECIES];
  // set all species input
  vm.num_species = get_species_inp(L, cdim, species);
  for (int s=0; s<vm.num_species; ++s) {
    vm.species[s] = species[s]->vm_species;
    vm.vdim = species[s]->vdim;

    app_lw->num_init[s] = species[s]->num_init;
    for (int i = 0; i < app_lw->num_init[s]; i++) {
      app_lw->proj_id[s][i] = species[s]->proj_id[i];

      app_lw->has_init_func[s][i] = species[s]->has_init_func[i];
      app_lw->init_func_ctx[s][i] = species[s]->init_func_ref[i];

      app_lw->has_density_init_func[s][i] = species[s]->has_density_init_func[i];
      app_lw->density_init_func_ctx[s][i] = species[s]->density_init_func_ref[i];

      app_lw->has_V_drift_init_func[s][i] = species[s]->has_V_drift_init_func[i];
      app_lw->V_drift_init_func_ctx[s][i] = species[s]->V_drift_init_func_ref[i];
      
      app_lw->has_temp_init_func[s][i] = species[s]->has_temp_init_func[i];
      app_lw->temp_init_func_ctx[s][i] = species[s]->temp_init_func_ref[i];

      app_lw->correct_all_moms[s][i] = species[s]->correct_all_moms[i];
      app_lw->iter_eps[s][i] = species[s]->iter_eps[i];
      app_lw->max_iter[s][i] = species[s]->max_iter[i];
      app_lw->use_last_converged[s][i] = species[s]->use_last_converged[i];
    }

    vm.species[s].num_init = app_lw->num_init[s];
    for (int i = 0; i < app_lw->num_init[s]; i++) {
      vm.species[s].projection[i].proj_id = app_lw->proj_id[s][i];

      if (species[s]->has_init_func[i]) {
        vm.species[s].projection[i].func = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_func = &app_lw->init_func_ctx[s][i];
      }

      if (species[s]->has_density_init_func[i]) {
        vm.species[s].projection[i].density = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_density = &app_lw->density_init_func_ctx[s][i];
      }

      if (species[s]->has_V_drift_init_func[i]) {
        vm.species[s].projection[i].V_drift = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_V_drift = &app_lw->V_drift_init_func_ctx[s][i];
      }

      if (species[s]->has_temp_init_func[i]) {
        vm.species[s].projection[i].temp = gkyl_lw_eval_cb;
        vm.species[s].projection[i].ctx_temp = &app_lw->temp_init_func_ctx[s][i];
      }

      vm.species[s].projection[i].correct_all_moms = app_lw->correct_all_moms[s][i];
      vm.species[s].projection[i].iter_eps = app_lw->iter_eps[s][i];
      vm.species[s].projection[i].max_iter = app_lw->max_iter[s][i];
      vm.species[s].projection[i].use_last_converged = app_lw->use_last_converged[s][i];
    }

    app_lw->collision_id[s] = species[s]->collision_id;

    app_lw->has_self_nu_func[s] = species[s]->has_self_nu_func;
    app_lw->self_nu_func_ctx[s] = species[s]->self_nu_func_ref;

    app_lw->collision_correct_all_moms[s] = species[s]->collision_correct_all_moms;

    vm.species[s].collisions.collision_id = app_lw->collision_id[s];

    if (species[s]->has_self_nu_func) {
      vm.species[s].collisions.self_nu = gkyl_lw_eval_cb;
      vm.species[s].collisions.ctx = &app_lw->self_nu_func_ctx[s];
    }

    vm.species[s].collisions.correct_all_moms = app_lw->collision_correct_all_moms[s];
  }

  // Set field input.
  vm.skip_field = glua_tbl_get_bool(L, "skipField", false);

  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct vlasov_field_lw *vmf = lua_touserdata(L, -1);
      if (vmf->magic == VLASOV_FIELD_DEFAULT) {

        vmf->init_ref.ndim = cdim;

        vm.field = vmf->vm_field;
        vm.skip_field = !vmf->evolve;

        app_lw->field_func_ctx = vmf->init_ref;
        vm.field.init = gkyl_lw_eval_cb;
        vm.field.ctx = &app_lw->field_func_ctx;
      }
    }
  }

  // create parallelism
  struct gkyl_comm *comm = 0;
  bool has_mpi = false;

  for (int d=0; d<cdim; ++d)
    vm.parallelism.cuts[d] = cuts[d]; 

#ifdef GKYL_HAVE_MPI
  with_lua_global(L, "GKYL_MPI_COMM") {
    if (lua_islightuserdata(L, -1)) {
      has_mpi = true;
      MPI_Comm mpi_comm = lua_touserdata(L, -1);
      comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
          .mpi_comm = mpi_comm,
        }
      );

    }
  }
#endif

  if (!has_mpi) {
    // if there is no proper MPI_Comm specifed, the assume we are a
    // serial sim
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {} );
  }
  vm.parallelism.comm = comm;

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int tot_cuts = 1; for (int d=0; d<cdim; ++d) tot_cuts *= cuts[d];

  if (tot_cuts != comm_sz) {
    printf("tot_cuts = %d (%d)\n", tot_cuts, comm_sz);
    luaL_error(L, "Number of ranks and cuts do not match!");
  }
  
  app_lw->app = gkyl_vlasov_app_new(&vm);

  gkyl_comm_release(comm);

  // create Lua userdata ...
  struct vlasov_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct vlasov_app_lw*));
  *l_app_lw = app_lw; // ... point it to the Lua app pointer

  // set metatable
  luaL_getmetatable(L, VLASOV_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool
static int
vm_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_vlasov_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to field. (time) -> bool
static int
vm_app_apply_ic_field(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_vlasov_app_apply_ic_field(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool
static int
vm_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->tstart);
  gkyl_vlasov_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute diagnostic moments. () -> bool
static int
vm_app_calc_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_calc_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool
static int
vm_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_vlasov_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated L2 norm of distribution function. (tm) -> bool
static int
vm_app_calc_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_vlasov_app_calc_integrated_L2_f(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool
static int
vm_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_vlasov_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool
static int
vm_app_write(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_vlasov_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool
static int
vm_app_write_field(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_vlasov_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool
static int
vm_app_write_species(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_vlasov_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write diagnostic moments to file (time, frame) -> bool
static int
vm_app_write_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_vlasov_app_write_mom(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool
static int
vm_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated L2 norm of f to file () -> bool
static int
vm_app_write_integrated_L2_f(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_integrated_L2_f(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool
static int
vm_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool
static int
vm_app_stat_write(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file
static void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}

// Run simulation. (num_steps) -> bool. num_steps is optional
static int
vm_app_run(lua_State *L)
{
  bool ret_status = true;

  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;
  struct gkyl_vlasov_app *app = app_lw->app;

  double tcurr = app_lw->tstart;
  double tend = app_lw->tend;
  double dt = tend-tcurr;
  long num_steps = luaL_optinteger(L, 2, INT_MAX);

  int nframe = app_lw->nframe;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
  gkyl_vlasov_app_calc_field_energy(app, tcurr);

  long step = 1;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_vlasov_app_calc_integrated_mom(app, tcurr);
    gkyl_vlasov_app_calc_field_energy(app, tcurr);    
    
    if (!status.success) {
      ret_status = false;
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_vlasov_app_stat_write(app);
  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_write_field_energy(app);

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation
static int
vm_app_gc(lua_State *L)
{
  struct vlasov_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  gkyl_vlasov_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor
static struct luaL_Reg vm_app_ctor[] = {
  { "new",  vm_app_new },
  { 0, 0 }
};

// App methods
static struct luaL_Reg vm_app_funcs[] = {
  { "apply_ic", vm_app_apply_ic },
  { "apply_ic_field", vm_app_apply_ic_field },
  { "apply_ic_species", vm_app_apply_ic_species },
  { "calc_mom", vm_app_calc_mom },
  { "calc_integrated_mom", vm_app_calc_integrated_mom },
  { "calc_integrated_L2_f", vm_app_calc_integrated_L2_f },
  { "calc_field_energy", vm_app_calc_field_energy },
  { "write", vm_app_write },
  { "write_field", vm_app_write_field },
  { "write_species", vm_app_write_species },
  { "write_mom", vm_app_write_mom },
  { "write_integrated_mom", vm_app_write_integrated_mom },
  { "write_integrated_L2_f", vm_app_write_integrated_L2_f },
  { "write_field_energy", vm_app_write_field_energy },
  { "stat_write", vm_app_stat_write },
  { "run", vm_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App
  do {
    luaL_newmetatable(L, VLASOV_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, vm_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, vm_app_funcs);
    
    luaL_register(L, "G0.Vlasov.App", vm_app_ctor);
    
  } while (0);

  // Register Species input struct
  do {
    luaL_newmetatable(L, VLASOV_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Vlasov.Species", vm_species_ctor);
  } while (0);

  // Register Field input struct
  do {
    luaL_newmetatable(L, VLASOV_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Vlasov.Field", vm_field_ctor);
  } while (0);
}

void
gkyl_vlasov_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif