#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_eqn_type.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_lw.h>
#include <gkyl_gyrokinetic_priv.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Magic IDs for use in distinguishing various species and field types.
enum gyrokinetic_magic_ids {
  GYROKINETIC_SPECIES_DEFAULT = 100, // Standard gyrokinetic species.
  GYROKINETIC_FIELD_DEFAULT, // Gyrokinetic Poisson equation.
};

/* *************** */
/* Species methods */
/* *************** */

// Metatable name for species input struct.
#define GYROKINETIC_SPECIES_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Species"

// Lua userdata object for constructing species input.
struct gyrokinetic_species_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_gyrokinetic_species gk_species; // Input struct to construct species.
  int vdim; // Velocity space dimensions.
  bool evolve; // Is this species evolved?

  enum gkyl_projection_id proj_id; // Projection type.

  bool has_init_func; // Is there an initialization function?
  struct lua_func_ctx init_func_ref; // Lua registry reference to initialization function.

  bool has_density_init_func; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ref; // Lua registry reference to density initialization function.

  bool has_Upar_init_func; // Is there a parallel velociy initialiation function?
  struct lua_func_ctx Upar_init_func_ref; // Lua registry reference to parallel velocity initialization function.

  bool has_temp_init_func; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ref; // Lua registry reference to temperature initialization function.

  bool correct_all_moms; // Are we correcting all moments in projection, or only density?

  enum gkyl_collision_id collision_id; // Collision type.
  
  bool has_self_nu_func; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ref; // Lua registry reference to self-collision frequency function.

  int num_cross_collisions; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  bool collision_correct_all_moms; // Are we correcting all moments in collisions, or only density?
  double collision_iter_eps; // Error tolerance for moment fixes in collisions (density is always exact).
  int collision_max_iter; // Maximum number of iterations for moment fixes in collisions.
  bool collision_use_last_converged; // Use last iteration value in collisions regardless of convergence?

  enum gkyl_source_id source_id; // Source type.

  int num_sources; // Number of projection objects in source.
  enum gkyl_projection_id source_proj_id[GKYL_MAX_PROJ]; // Projection type in source.

  bool source_has_init_func[GKYL_MAX_PROJ]; // Is there an initialization function in source?
  struct lua_func_ctx source_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to initialization function in source.

  bool source_has_density_init_func[GKYL_MAX_PROJ]; // Is there a density initialization function in source?
  struct lua_func_ctx source_density_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to density initialization function in source.

  bool source_has_Upar_init_func[GKYL_MAX_PROJ]; // Is there a parallel velocity initialization function in source?
  struct lua_func_ctx source_Upar_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to parallel velocity initialization function in source.

  bool source_has_temp_init_func[GKYL_MAX_PROJ]; // Is there a temperature initialization function in source?
  struct lua_func_ctx source_temp_init_func_ref[GKYL_MAX_PROJ]; // Lua registry reference to temperature initialization function in source.
};

static int
gyrokinetic_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_gyrokinetic_species gk_species = { };
  
  gk_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  gk_species.mass = glua_tbl_get_number(L, "mass", 1.0);
  gk_species.polarization_density = glua_tbl_get_number(L, "polarizationDensity", 0.0);
  gk_species.no_by = glua_tbl_get_bool(L, "noBy", false);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);

    for (int d = 0; d < vdim; d++) {
      gk_species.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
    }
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d = 0; d < vdim; d++) {
      gk_species.lower[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < vdim; d++) {
      gk_species.upper[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  bool evolve = glua_tbl_get_integer(L, "evolve", true);

  with_lua_tbl_tbl(L, "diagnostics") {
    int num_diag_moments = glua_objlen(L);

    int n = 0;
    for (int i = 0; i < num_diag_moments; i ++) {
      const char *mom = glua_tbl_iget_string(L, i+1, "");

      if (is_moment_name_valid(mom)) {
        strcpy(gk_species.diag_moments[n++], mom);
      }
    }

    gk_species.num_diag_moments = n;
  }

  with_lua_tbl_tbl(L, "bcx") {
    with_lua_tbl_tbl(L, "lower") {
      gk_species.bcx.lower.type = glua_tbl_get_integer(L, "type", 0);
    }
    
    with_lua_tbl_tbl(L, "upper") {
      gk_species.bcx.upper.type = glua_tbl_get_integer(L, "type", 0);
    }
  }

  with_lua_tbl_tbl(L, "bcy") {
    with_lua_tbl_tbl(L, "lower") {
      gk_species.bcy.lower.type = glua_tbl_get_integer(L, "type", 0);
    }

    with_lua_tbl_tbl(L, "upper") {
      gk_species.bcy.upper.type = glua_tbl_get_integer(L, "type", 0);
    }
  }

  with_lua_tbl_tbl(L, "bcz") {
    with_lua_tbl_tbl(L, "lower") {
      gk_species.bcz.lower.type = glua_tbl_get_integer(L, "type", 0);
    }

    with_lua_tbl_tbl(L, "upper") {
      gk_species.bcz.upper.type = glua_tbl_get_integer(L, "type", 0);
    }
  }
  
  enum gkyl_projection_id proj_id = GKYL_PROJ_FUNC;

  bool has_init_func = false;
  int init_func_ref = LUA_NOREF;

  bool has_density_init_func = false;
  int density_init_func_ref = LUA_NOREF;

  bool has_Upar_init_func = false;
  int Upar_init_func_ref = LUA_NOREF;

  bool has_temp_init_func = false;
  int temp_init_func_ref = LUA_NOREF;

  bool correct_all_moms = false;

  with_lua_tbl_tbl(L, "projection") {
    proj_id = glua_tbl_get_integer(L, "projectionID", 0);

    if (glua_tbl_get_func(L, "init")) {
      init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_init_func = true;
    }

    if (glua_tbl_get_func(L, "densityInit")) {
      density_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_density_init_func = true;
    }

    if (glua_tbl_get_func(L, "parallelVelocityInit")) {
      Upar_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_Upar_init_func = true;
    }

    if (glua_tbl_get_func(L, "temperatureInit")) {
      temp_init_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_temp_init_func = true;
    }

    correct_all_moms = glua_tbl_get_bool(L, "correctAllMoments", true);
  }

  enum gkyl_collision_id collision_id = GKYL_NO_COLLISIONS;

  bool has_self_nu_func = false;
  int self_nu_func_ref = LUA_NOREF;

  int num_cross_collisions = 0;
  char collide_with[GKYL_MAX_SPECIES][128];

  bool collision_correct_all_moms = false;
  double collision_iter_eps = pow(10.0, -12.0);
  int collision_max_iter = 100;
  bool collision_use_last_converged = true;

  with_lua_tbl_tbl(L, "collisions") {
    collision_id = glua_tbl_get_integer(L, "collisionID", 0);

    if (glua_tbl_get_func(L, "selfNu")) {
      self_nu_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
      has_self_nu_func = true;
    }

    num_cross_collisions = glua_tbl_get_integer(L, "numCrossCollisions", 0);
    with_lua_tbl_tbl(L, "collideWith") {
      for (int i = 0; i < num_cross_collisions; i++) {
        const char* collide_with_char = glua_tbl_iget_string(L, i + 1, "");
        strcpy(collide_with[i], collide_with_char);
      }
    }

    collision_correct_all_moms = glua_tbl_get_bool(L, "correctAllMoments", true);
    collision_iter_eps = glua_tbl_get_number(L, "iterationEpsilon", pow(10.0, -12.0));
    collision_max_iter = glua_tbl_get_integer(L, "maxIterations", 100);
    collision_use_last_converged = glua_tbl_get_bool(L, "useLastConverged", true);
  }

  enum gkyl_source_id source_id = GKYL_NO_SOURCE;

  int num_sources = 0;
  enum gkyl_projection_id source_proj_id[GKYL_MAX_PROJ];

  bool source_has_init_func[GKYL_MAX_PROJ];
  int source_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_density_init_func[GKYL_MAX_PROJ];
  int source_density_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_Upar_init_func[GKYL_MAX_PROJ];
  int source_Upar_init_func_ref[GKYL_MAX_PROJ];

  bool source_has_temp_init_func[GKYL_MAX_PROJ];
  int source_temp_init_func_ref[GKYL_MAX_PROJ];

  with_lua_tbl_tbl(L, "source") {
    source_id = glua_tbl_get_integer(L, "sourceID", 0);
    num_sources = glua_tbl_get_integer(L, "numSources", 0);

    with_lua_tbl_tbl(L, "projections") {
      for (int i = 0; i < num_sources; i++) {
        if (glua_tbl_iget_tbl(L, i + 1)) {
          source_proj_id[i] = glua_tbl_get_integer(L, "projectionID", 0);

          source_init_func_ref[i] = LUA_NOREF;
          source_has_init_func[i] = false;
          if (glua_tbl_get_func(L, "init")) {
            source_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_init_func[i] = true;
          }

          source_density_init_func_ref[i] = LUA_NOREF;
          source_has_density_init_func[i] = false;
          if (glua_tbl_get_func(L, "densityInit")) {
            source_density_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_density_init_func[i] = true;
          }

          source_Upar_init_func_ref[i] = LUA_NOREF;
          source_has_Upar_init_func[i] = false;
          if (glua_tbl_get_func(L, "parallelVelocityInit")) {
            source_Upar_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_Upar_init_func[i] = true;
          }

          source_temp_init_func_ref[i] = LUA_NOREF;
          source_has_temp_init_func[i] = false;
          if (glua_tbl_get_func(L, "temperatureInit")) {
            source_temp_init_func_ref[i] = luaL_ref(L, LUA_REGISTRYINDEX);
            source_has_temp_init_func[i] = true;
          }

          lua_pop(L, 1);
        }
      }
    }
  }
  
  struct gyrokinetic_species_lw *gks_lw = lua_newuserdata(L, sizeof(*gks_lw));
  gks_lw->magic = GYROKINETIC_SPECIES_DEFAULT;
  gks_lw->vdim = vdim;
  gks_lw->evolve = evolve;
  gks_lw->gk_species = gk_species;

  gks_lw->proj_id = proj_id;

  gks_lw->has_init_func = has_init_func;
  gks_lw->init_func_ref = (struct lua_func_ctx) {
    .func_ref = init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_density_init_func = has_density_init_func;
  gks_lw->density_init_func_ref = (struct lua_func_ctx) {
    .func_ref = density_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_Upar_init_func = has_Upar_init_func;
  gks_lw->Upar_init_func_ref = (struct lua_func_ctx) {
    .func_ref = Upar_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->has_temp_init_func = has_temp_init_func;
  gks_lw->temp_init_func_ref = (struct lua_func_ctx) {
    .func_ref = temp_init_func_ref,
    .ndim = 0, // This will be set later.
    .nret = 1,
    .L = L,
  };

  gks_lw->correct_all_moms = correct_all_moms;

  gks_lw->source_id = source_id;
  gks_lw->num_sources = num_sources;
  
  for (int i = 0; i < num_sources; i++) {
    gks_lw->source_proj_id[i] = source_proj_id[i];

    gks_lw->source_has_init_func[i] = source_has_init_func[i];
    gks_lw->source_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    gks_lw->source_has_density_init_func[i] = source_has_density_init_func[i];
    gks_lw->source_density_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_density_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    gks_lw->source_has_Upar_init_func[i] = source_has_Upar_init_func[i];
    gks_lw->source_Upar_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_Upar_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };

    gks_lw->source_has_temp_init_func[i] = source_has_temp_init_func[i];
    gks_lw->source_temp_init_func_ref[i] = (struct lua_func_ctx) {
      .func_ref = source_temp_init_func_ref[i],
      .ndim = 0, // This will be set later.
      .nret = 1,
      .L = L,
    };
  }

  gks_lw->collision_id = collision_id;

  gks_lw->has_self_nu_func = has_self_nu_func;
  gks_lw->self_nu_func_ref = (struct lua_func_ctx) {
    .func_ref = self_nu_func_ref,
    .ndim = 0,
    .nret = 1,
    .L = L,
  };

  gks_lw->num_cross_collisions = num_cross_collisions;
  for (int i = 0; i < num_cross_collisions; i++) {
    strcpy(gks_lw->collide_with[i], collide_with[i]);
  }

  gks_lw->collision_correct_all_moms = collision_correct_all_moms;
  gks_lw->collision_iter_eps = collision_iter_eps;
  gks_lw->collision_max_iter = collision_max_iter;
  gks_lw->collision_use_last_converged = collision_use_last_converged;
  
  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor.
static struct luaL_Reg gk_species_ctor[] = {
  { "new", gyrokinetic_species_lw_new },
  { 0, 0 }
};

/* ************* */
/* Field methods */
/* ************* */

// Metatable name for field input struct.
#define GYROKINETIC_FIELD_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Field"

// Lua userdata object for constructing field input.
struct gyrokinetic_field_lw {
  int magic; // This must be first element in the struct.
  
  struct gkyl_gyrokinetic_field gk_field; // Input struct to construct field.
};

static int
gyrokinetic_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_gyrokinetic_field gk_field = { };

  gk_field.gkfield_id = glua_tbl_get_integer(L, "fieldID", 0);
  gk_field.electron_mass = glua_tbl_get_number(L, "electronMass", 0.0);
  gk_field.electron_charge = glua_tbl_get_number(L, "electronCharge", 0.0);
  gk_field.electron_temp = glua_tbl_get_number(L, "electronTemperature", 0.0);

  gk_field.fem_parbc = glua_tbl_get_integer(L, "femParBc", 0);
  gk_field.kperpSq = glua_tbl_get_number(L, "kPerpSq", 0.0);

  with_lua_tbl_tbl(L, "poissonBcs") {
    with_lua_tbl_tbl(L, "lowerType") {
      int nbc = glua_objlen(L);
      
      for (int i = 0; i < nbc; i++) {
        gk_field.poisson_bcs.lo_type[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }

    with_lua_tbl_tbl(L, "upperType") {
      int nbc = glua_objlen(L);

      for (int i = 0; i < nbc; i++) {
        gk_field.poisson_bcs.up_type[i] = glua_tbl_iget_integer(L, i + 1, 0);
      }
    }

    with_lua_tbl_tbl(L, "lowerValue") {
      int nbc = glua_objlen(L);

      for (int i = 0; i < nbc; i++) {
        struct gkyl_poisson_bc_value lower_bc = {
          .v = { glua_tbl_iget_number(L, i + 1, 0.0) },
        };

        gk_field.poisson_bcs.lo_value[i] = lower_bc;
      }
    }

    with_lua_tbl_tbl(L, "upperValue") {
      int nbc = glua_objlen(L);

      for (int i = 0; i < nbc; i++) {
        struct gkyl_poisson_bc_value upper_bc = {
          .v = { glua_tbl_iget_number(L, i + 1, 0.0) },
        };

        gk_field.poisson_bcs.up_value[i] = upper_bc;
      }
    }
  }

  struct gyrokinetic_field_lw *gkf_lw = lua_newuserdata(L, sizeof(*gkf_lw));

  gkf_lw->magic = GYROKINETIC_FIELD_DEFAULT;
  gkf_lw->gk_field = gk_field;
  
  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Field constructor.
static struct luaL_Reg gk_field_ctor[] = {
  { "new",  gyrokinetic_field_lw_new },
  { 0, 0 }
};

/* *********** */
/* App methods */
/* *********** */

// Metatable name for top-level Gyrokinetic App.
#define GYROKINETIC_APP_METATABLE_NM "GkeyllZero.App.Gyrokinetic"

// Lua userdata object for holding Gyrokinetic app and run parameters.
struct gyrokinetic_app_lw {
  gkyl_gyrokinetic_app *app; // Gyrokinetic app object.

  struct lua_func_ctx mapc2p_ctx; // Function context for mapc2p.
  struct lua_func_ctx bmag_ctx; // Function context for bmag.

  enum gkyl_projection_id proj_id[GKYL_MAX_SPECIES]; // Projection type.

  bool has_init_func[GKYL_MAX_SPECIES]; // Is there an initialization function?
  struct lua_func_ctx init_func_ctx[GKYL_MAX_SPECIES]; // Context for initialization function.

  bool has_density_init_func[GKYL_MAX_SPECIES]; // Is there a density initialization function?
  struct lua_func_ctx density_init_func_ctx[GKYL_MAX_SPECIES]; // Context for density initialization function.

  bool has_Upar_init_func[GKYL_MAX_SPECIES]; // Is there a parallel velocity initialization function?
  struct lua_func_ctx Upar_init_func_ctx[GKYL_MAX_SPECIES]; // Context for parallel velocity initialziation function.
  
  bool has_temp_init_func[GKYL_MAX_SPECIES]; // Is there a temperature initialization function?
  struct lua_func_ctx temp_init_func_ctx[GKYL_MAX_SPECIES]; // Context for temperature initialization function.

  bool correct_all_moms[GKYL_MAX_SPECIES]; // Are we correcting all moments in projection, or only density?

  enum gkyl_collision_id collision_id[GKYL_MAX_SPECIES]; // Collision type.

  bool has_self_nu_func[GKYL_MAX_SPECIES]; // Is there a self-collision frequency function?
  struct lua_func_ctx self_nu_func_ctx[GKYL_MAX_SPECIES]; // Context for self-collision frequency function.

  int num_cross_collisions[GKYL_MAX_SPECIES]; // Number of species that we cross-collide with.
  char collide_with[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES][128]; // Names of species that we cross-collide with.

  bool collision_correct_all_moms[GKYL_MAX_SPECIES]; // Are we correcting all moments in collisions, or only density?
  double collision_iter_eps[GKYL_MAX_SPECIES]; // Error tolerance for moment fixes in collision (density is always exact).
  int collision_max_iter[GKYL_MAX_SPECIES]; // Maximum number of iterations for moment fixes in collisions.
  bool collision_use_last_converged[GKYL_MAX_SPECIES]; // Use last iteration value in collisions regardless of convergence?

  enum gkyl_source_id source_id[GKYL_MAX_SPECIES]; // Source type.

  int num_sources[GKYL_MAX_SPECIES]; // Number of projection objects in source.
  enum gkyl_projection_id source_proj_id[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Projection type in source.

  bool source_has_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there an initialization function in source?
  struct lua_func_ctx source_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for initialization function in source.

  bool source_has_density_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a density initialization function in source?
  struct lua_func_ctx source_density_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for density initialization function in source.

  bool source_has_Upar_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a parallel velocity initialization function in source?
  struct lua_func_ctx source_Upar_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for parallel velocity initialization function in source.

  bool source_has_temp_init_func[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Is there a temperature initialization function in source?
  struct lua_func_ctx source_temp_init_func_ctx[GKYL_MAX_SPECIES][GKYL_MAX_PROJ]; // Context for temperature initialization function in source.
  
  double t_start, t_end; // Start and end times of simulation.
  int num_frames; // Number of data frames to write.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct gyrokinetic_species_lw *species[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1 };
  
  int curr = 0;
  lua_pushnil(L); // Initial key is nil.
  while (lua_next(L, TKEY) != 0) {
    // Key at TKEY and value at TVAL.
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct gyrokinetic_species_lw *gks = lua_touserdata(L, TVAL);

      if (gks->magic == GYROKINETIC_SPECIES_DEFAULT) {
        if (gks->has_init_func) {
          gks->init_func_ref.ndim = cdim + gks->vdim;
        }

        if (gks->has_density_init_func) {
          gks->density_init_func_ref.ndim = cdim + gks->vdim;
        }
        
        if (gks->has_Upar_init_func) {
          gks->Upar_init_func_ref.ndim = cdim + gks->vdim;
        }

        if (gks->has_temp_init_func) {
          gks->temp_init_func_ref.ndim = cdim + gks->vdim;
        }

        if (gks->has_self_nu_func) {
          gks->self_nu_func_ref.ndim = cdim + gks->vdim;
        }

        for (int i = 0; i < gks->num_sources; i++) {
          if (gks->source_has_init_func[i]) {
            gks->source_init_func_ref[i].ndim = cdim + gks->vdim;
          }

          if (gks->source_has_density_init_func[i]) {
            gks->source_density_init_func_ref[i].ndim = cdim + gks->vdim;
          }

          if (gks->source_has_Upar_init_func[i]) {
            gks->source_Upar_init_func_ref[i].ndim = cdim + gks->vdim;
          }

          if (gks->source_has_temp_init_func[i]) {
            gks->source_temp_init_func_ref[i].ndim = cdim + gks->vdim;
          }
        }
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(gks->gk_species.name, key);
        }
        species[curr++] = gks;
      }
    }
    lua_pop(L, 1);
  }

  return curr;
}

// Create top-level App object.
static int
gk_app_new(lua_State *L)
{
  struct gyrokinetic_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-gyrokinetic" is used.
  const char *sim_name = "g0-gyrokinetic";

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

  struct gkyl_gk gk = { }; // Input table for app.

  strcpy(gk.name, sim_name);
  
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    gk.cdim = cdim = glua_objlen(L);

    for (int d = 0; d < cdim; d++) {
      gk.cells[d] = glua_tbl_iget_integer(L, d + 1, 0);
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
      gk.lower[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d = 0; d < cdim; d++) {
      gk.upper[d] = glua_tbl_iget_number(L, d + 1, 0.0);
    }
  }

  gk.cfl_frac = glua_tbl_get_number(L, "cflFrac", 0.95);
  gk.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  gk.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  gk.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      gk.num_periodic_dir = glua_objlen(L);

      for (int d = 0; d < gk.num_periodic_dir; d++) {
        // Indices are off by 1 between Lua and C.
        gk.periodic_dirs[d] = glua_tbl_iget_integer(L, d + 1, 0) - 1;
      }
    }
  }

  with_lua_tbl_tbl(L, "geometry") {
    gk.geometry.geometry_id = glua_tbl_get_integer(L, "geometryID", 0);

    with_lua_tbl_tbl(L, "world") {
      for (int d = 0; d < 3; d++) {
        gk.geometry.world[d] = glua_tbl_iget_number(L, d + 1, 0.0);
      }
    }

    gk.geometry.c2p_ctx = 0;
    gk.geometry.mapc2p = 0;
    bool has_mapc2p = false;
    int mapc2p_ref = LUA_NOREF;
    if (glua_tbl_get_func(L, "mapc2p")) {
      has_mapc2p = true;
      mapc2p_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    }

    gk.geometry.bmag_ctx = 0;
    gk.geometry.bmag_func = 0;
    bool has_bmag_func = false;
    int bmag_func_ref = LUA_NOREF;
    if (glua_tbl_get_func(L, "bmagFunc")) {
      has_bmag_func = true;
      bmag_func_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    }

    if (has_mapc2p) {
      app_lw->mapc2p_ctx = (struct lua_func_ctx) {
        .func_ref = mapc2p_ref,
        .ndim = 3,
        .nret = 3,
        .L = L,
      };
      gk.geometry.mapc2p = gkyl_lw_eval_cb;
      gk.geometry.c2p_ctx = &app_lw->mapc2p_ctx;
    }

    if (has_bmag_func) {
      app_lw->bmag_ctx = (struct lua_func_ctx) {
        .func_ref = bmag_func_ref,
        .ndim = 3,
        .nret = 1,
        .L = L,
      };
      gk.geometry.bmag_func = gkyl_lw_eval_cb;
      gk.geometry.bmag_ctx = &app_lw->bmag_ctx;
    }
  }

  struct gyrokinetic_species_lw *species[GKYL_MAX_SPECIES];

  // Set all species input.
  gk.num_species = get_species_inp(L, cdim, species);
  for (int s = 0; s < gk.num_species; s++) {
    gk.species[s] = species[s]->gk_species;
    gk.vdim = species[s]->vdim;

    app_lw->proj_id[s] = species[s]->proj_id;

    app_lw->has_init_func[s] = species[s]->has_init_func;
    app_lw->init_func_ctx[s] = species[s]->init_func_ref;

    app_lw->has_density_init_func[s] = species[s]->has_density_init_func;
    app_lw->density_init_func_ctx[s] = species[s]->density_init_func_ref;

    app_lw->has_Upar_init_func[s] = species[s]->has_Upar_init_func;
    app_lw->Upar_init_func_ctx[s] = species[s]->Upar_init_func_ref;
    
    app_lw->has_temp_init_func[s] = species[s]->has_temp_init_func;
    app_lw->temp_init_func_ctx[s] = species[s]->temp_init_func_ref;

    app_lw->correct_all_moms[s] = species[s]->correct_all_moms;

    gk.species[s].projection.proj_id = app_lw->proj_id[s];

    if (species[s]->has_init_func) {
      gk.species[s].projection.func = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_func = &app_lw->init_func_ctx[s];
    }

    if (species[s]->has_density_init_func) {
      gk.species[s].projection.density = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_density = &app_lw->density_init_func_ctx[s];
    }

    if (species[s]->has_Upar_init_func) {
      gk.species[s].projection.upar = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_upar = &app_lw->Upar_init_func_ctx[s];
    }

    if (species[s]->has_temp_init_func) {
      gk.species[s].projection.temp = gkyl_lw_eval_cb;
      gk.species[s].projection.ctx_temp = &app_lw->temp_init_func_ctx[s];
    }

    gk.species[s].projection.correct_all_moms = app_lw->correct_all_moms[s];

    app_lw->collision_id[s] = species[s]->collision_id;

    app_lw->has_self_nu_func[s] = species[s]->has_self_nu_func;
    app_lw->self_nu_func_ctx[s] = species[s]->self_nu_func_ref;

    app_lw->num_cross_collisions[s] = species[s]->num_cross_collisions;
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(app_lw->collide_with[s][i], species[s]->collide_with[i]);
    }

    app_lw->collision_correct_all_moms[s] = species[s]->collision_correct_all_moms;
    app_lw->collision_iter_eps[s] = species[s]->collision_iter_eps;
    app_lw->collision_max_iter[s] = species[s]->collision_max_iter;
    app_lw->collision_use_last_converged[s] = species[s]->collision_use_last_converged;

    gk.species[s].collisions.collision_id = app_lw->collision_id[s];

    if (species[s]->has_self_nu_func) {
      gk.species[s].collisions.self_nu = gkyl_lw_eval_cb;
      gk.species[s].collisions.ctx = &app_lw->self_nu_func_ctx[s];
    }

    gk.species[s].collisions.num_cross_collisions = app_lw->num_cross_collisions[s];
    for (int i = 0; i < app_lw->num_cross_collisions[s]; i++) {
      strcpy(gk.species[s].collisions.collide_with[i], app_lw->collide_with[s][i]);
    }

    gk.species[s].collisions.correct_all_moms = app_lw->collision_correct_all_moms[s];
    gk.species[s].collisions.iter_eps = app_lw->collision_iter_eps[s];
    gk.species[s].collisions.max_iter = app_lw->collision_max_iter[s];
    gk.species[s].collisions.use_last_converged = app_lw->collision_use_last_converged[s];

    app_lw->source_id[s] = species[s]->source_id;

    app_lw->num_sources[s] = species[s]->num_sources;
    for (int i = 0; i < app_lw->num_sources[s]; i++) {
      app_lw->source_proj_id[s][i] = species[s]->source_proj_id[i];

      app_lw->source_has_init_func[s][i] = species[s]->source_has_init_func[i];
      app_lw->source_init_func_ctx[s][i] = species[s]->source_init_func_ref[i];

      app_lw->source_has_density_init_func[s][i] = species[s]->source_has_density_init_func[i];
      app_lw->source_density_init_func_ctx[s][i] = species[s]->source_density_init_func_ref[i];

      app_lw->source_has_Upar_init_func[s][i] = species[s]->source_has_Upar_init_func[i];
      app_lw->source_Upar_init_func_ctx[s][i] = species[s]->source_Upar_init_func_ref[i];

      app_lw->source_has_temp_init_func[s][i] = species[s]->source_has_temp_init_func[i];
      app_lw->source_temp_init_func_ctx[s][i] = species[s]->source_temp_init_func_ref[i];
    }

    gk.species[s].source.source_id = app_lw->source_id[s];

    gk.species[s].source.num_sources = app_lw->num_sources[s];
    for (int i = 0; i < app_lw->num_sources[s]; i++) {
      gk.species[s].source.projection[i].proj_id = app_lw->source_proj_id[s][i];

      if (species[s]->source_has_init_func[i]) {
        gk.species[s].source.projection[i].func = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_func = &app_lw->source_init_func_ctx[s][i];
      }

      if (species[s]->source_has_density_init_func[i]) {
        gk.species[s].source.projection[i].density = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_density = &app_lw->source_density_init_func_ctx[s][i];
      }

      if (species[s]->source_has_Upar_init_func[i]) {
        gk.species[s].source.projection[i].upar = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_upar = &app_lw->source_Upar_init_func_ctx[s][i];
      }

      if (species[s]->source_has_temp_init_func[i]) {
        gk.species[s].source.projection[i].temp = gkyl_lw_eval_cb;
        gk.species[s].source.projection[i].ctx_temp = &app_lw->source_temp_init_func_ctx[s][i];
      }
    }
  }

  // Set field input.
  gk.skip_field = glua_tbl_get_bool(L, "skipField", false);

  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_field_lw *gkf = lua_touserdata(L, -1);

      if (gkf->magic == GYROKINETIC_FIELD_DEFAULT) {
        gk.field = gkf->gk_field;
      }
    }
  }

  // Create parallelism.
  struct gkyl_comm *comm = 0;
  bool has_mpi = false;

  for (int d = 0; d < cdim; d++) {
    gk.parallelism.cuts[d] = cuts[d]; 
  }

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
    // If there is no proper MPI_Comm specifed, the assume we are a
    // serial sim.
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {} );
  }
  gk.parallelism.comm = comm;

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int tot_cuts = 1;
  for (int d = 0; d < cdim; d++) {
    tot_cuts *= cuts[d];
  }

  if (tot_cuts != comm_sz) {
    printf("tot_cuts = %d (%d)\n", tot_cuts, comm_sz);
    luaL_error(L, "Number of ranks and cuts do not match!");
  }
  
  app_lw->app = gkyl_gyrokinetic_app_new(&gk);

  gkyl_comm_release(comm);

  // Create Lua userdata.
  struct gyrokinetic_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct gyrokinetic_app_lw*));
  *l_app_lw = app_lw; // Point it to the Lua app pointer.

  // Set metatable.
  luaL_getmetatable(L, GYROKINETIC_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool.
static int
gk_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->t_start);
  gkyl_gyrokinetic_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool.
static int
gk_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->t_start);
  gkyl_gyrokinetic_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool.
static int
gk_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_gyrokinetic_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy (L2 norm of each field
// component). (tm) -> bool.
static int
gk_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_gyrokinetic_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool.
static int
gk_app_write(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_gyrokinetic_app_write(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write field to file (time, frame) -> bool.
static int
gk_app_write_field(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_gyrokinetic_app_write_field(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write species solution to file (sidx, time, frame) -> bool.
static int
gk_app_write_species(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double tm = luaL_checknumber(L, 3);
  int frame = luaL_checkinteger(L, 4);
  gkyl_gyrokinetic_app_write_species(app_lw->app, sidx, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write diagnostic moments to file (time, frame) -> bool.
static int
gk_app_write_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  int frame = luaL_checkinteger(L, 3);
  gkyl_gyrokinetic_app_write_mom(app_lw->app, tm, frame);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated moments to file () -> bool.
static int
gk_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}


// Write integrated field energy to file () -> bool.
static int
gk_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool.
static int
gk_app_stat_write(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file.
static void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_gyrokinetic_app_write(app, t_curr, iot->curr - 1);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, iot->curr - 1);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }
}

// Step message context.
struct step_message_trigs {
  int log_count; // Number of times logging called.
  int tenth, p1c; 
  struct gkyl_tm_trigger log_trig; // 10% trigger.
  struct gkyl_tm_trigger log_trig_1p; // 1% trigger.
};

// Write log message to console.
static void
write_step_message(const struct gkyl_gyrokinetic_app *app, struct step_message_trigs *trigs, int step, double t_curr, double dt_next)
{
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig, t_curr)) {
    if (trigs->log_count > 0) {
      gkyl_gyrokinetic_app_cout(app, stdout, " Step %6d at time %#11.8g.  Time-step  %.6e.  Completed %g%s\n", step, t_curr, dt_next, trigs->tenth * 10.0, "%");
    }
    else {
      trigs->log_count += 1;
    }
    
    trigs->tenth += 1;
  }
  if (gkyl_tm_trigger_check_and_bump(&trigs->log_trig_1p, t_curr)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "%d", trigs->p1c);
    trigs->p1c = (trigs->p1c+1) % 10;
  }
}

// Run simulation. (num_steps) -> bool. num_steps is optional.
static int
gk_app_run(lua_State *L)
{
  bool ret_status = true;

  // Create app object.
  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;
  struct gkyl_gyrokinetic_app *app = app_lw->app;

  // Initial and final simulation times.
  double t_curr = app_lw->t_start, t_end = app_lw->t_end;
  long num_steps = luaL_optinteger(L, 2, INT_MAX);

  gkyl_gyrokinetic_app_cout(app, stdout, "Initializing Gyrokinetic Simulation ...\n");

  // Initialize simulation.
  bool is_restart = false;
  lua_getglobal(L, "GKYL_IS_RESTART");
  if (lua_toboolean(L, -1)) {
    is_restart = true;
  }
  lua_pop(L, 1);

  int restart_frame = 0;
  lua_getglobal(L, "GKYL_RESTART_FRAME");
  restart_frame = lua_tointeger(L, -1);
  lua_pop(L, 1);

  int frame_curr = 0;
  if (is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }

  int num_frames = app_lw->num_frames;
  // Triggers for IO and logging.
  struct gkyl_tm_trigger io_trig = { .dt = (t_end - t_curr) / num_frames, .tcurr = t_curr, .curr = frame_curr };

  struct step_message_trigs m_trig = {
    .log_count = 0,
    .tenth = t_curr > 0.0 ? 0.0 : (int) floor(t_curr / t_end * 10.0),
    .p1c = t_curr > 0.0 ? 0.0 : (int) floor(t_curr / t_end * 100.0) % 10,
    .log_trig = { .dt = (t_end - t_curr) / 10.0 },
    .log_trig_1p = { .dt = (t_end - t_curr) / 100.0 },
  };

  struct timespec tm_ic0 = gkyl_wall_clock();
  gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  gkyl_gyrokinetic_app_cout(app, stdout, "Initialization completed in %g sec\n\n", gkyl_time_diff_now_sec(tm_ic0));

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = app_lw->dt_failure_tol;
  int num_failures = 0, num_failures_max = app_lw->num_failures_max;

  bool use_verbose = false;
  lua_getglobal(L, "GKYL_USE_VERBOSE");
  if (lua_toboolean(L, -1)) {
    use_verbose = true;
  }
  lua_pop(L, 1);

  long num_steps_new = -1;
  lua_getglobal(L, "GKYL_NUM_STEPS");
  num_steps_new = lua_tointeger(L, -1);
  lua_pop(L, 1);
  if (num_steps_new != -1) {
    num_steps = num_steps_new;
  }

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (use_verbose) {
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    }
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    if (use_verbose) {
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    }

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    if (!use_verbose) {
      write_step_message(app, &m_trig, step, t_curr, status.dt_suggested);
    }

    step += 1;
  }

  write_data(&io_trig, app, t_curr, false);
  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation.
static int
gk_app_gc(lua_State *L)
{
  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor.
static struct luaL_Reg gk_app_ctor[] = {
  { "new",  gk_app_new },
  { 0, 0 }
};

// App methods.
static struct luaL_Reg gk_app_funcs[] = {
  { "apply_ic", gk_app_apply_ic },
  { "apply_ic_species", gk_app_apply_ic_species },
  { "calc_integrated_mom", gk_app_calc_integrated_mom },
  { "calc_field_energy", gk_app_calc_field_energy },
  { "write", gk_app_write },
  { "write_field", gk_app_write_field },
  { "write_species", gk_app_write_species },
  { "write_mom", gk_app_write_mom },
  { "write_integrated_mom", gk_app_write_integrated_mom },
  { "write_field_energy", gk_app_write_field_energy },
  { "stat_write", gk_app_stat_write },
  { "run", gk_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App.
  do {
    luaL_newmetatable(L, GYROKINETIC_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, gk_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, gk_app_funcs);
    
    luaL_register(L, "G0.Gyrokinetic.App", gk_app_ctor);
    
  }
  while (0);

  // Register Species input struct.
  do {
    luaL_newmetatable(L, GYROKINETIC_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Species", gk_species_ctor);
  }
  while (0);

  // Register Field input struct.
  do {
    luaL_newmetatable(L, GYROKINETIC_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Field", gk_field_ctor);
  }
  while (0);
}

void
gkyl_gyrokinetic_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif
