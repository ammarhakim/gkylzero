#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_lw.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

// Magic IDs for use in distinguishing various species and field types
enum gyrokinetic_magic_ids {
  GYROKINETIC_SPECIES_DEFAULT = 100, // Gyrokinetic model
  GYROKINETIC_FIELD_DEFAULT, // Poisson equation
  GYROKINETIC_COLLISIONS_DEFAULT, // LBO Collisions
  GYROKINETIC_DIFFUSION_DEFAULT, // Diffusion operator
  GYROKINETIC_SOURCE_DEFAULT, // Sources
  GYROKINETIC_GEOMETRY_DEFAULT, // Geometry
  GYROKINETIC_PROJECTION_DEFAULT, // Geometry
};

/* ***************** */
/* Projection methods */
/* ***************** */

// Metatable name for field input struct
#define GYROKINETIC_PROJECTION_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Projection"

// Lua userdata object for constructing field input
struct gyrokinetic_projection_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_projection gyrokinetic_projection; // input struct to construct projection

  struct lua_func_ctx density_ref; // Lua registery reference to density
  struct lua_func_ctx upar_ref; // Lua registery reference to upar
  struct lua_func_ctx temp_ref; // Lua registery reference to temperature
};

static int
gyrokinetic_projection_lw_new(lua_State *L)
{
  struct gkyl_gyrokinetic_projection gyrokinetic_projection = { };

  gyrokinetic_projection.proj_id = GKYL_PROJ_MAXWELLIAN_PRIM; 

  // Register functions for density, upar, and temperature 
  int density_ref = LUA_NOREF;
  int upar_ref = LUA_NOREF;
  int temp_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "density")) 
    density_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Projection must have an \"density\" function if projecting a Maxwellian!");
  if (glua_tbl_get_func(L, "upar"))
    upar_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Projection must have an \"upar\" function if projecting a Maxwellian!");
  if (glua_tbl_get_func(L, "temp")) 
    temp_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Projection must have an \"temp\" function if projecting a Maxwellian!");

  struct gyrokinetic_projection_lw *gyrokinetic_p_lw = lua_newuserdata(L, sizeof(*gyrokinetic_p_lw));

  gyrokinetic_p_lw->magic = GYROKINETIC_PROJECTION_DEFAULT;
  gyrokinetic_p_lw->gyrokinetic_projection = gyrokinetic_projection;

  gyrokinetic_p_lw->density_ref = (struct lua_func_ctx) {
    .func_ref = density_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };

  gyrokinetic_p_lw->upar_ref = (struct lua_func_ctx) {
    .func_ref = upar_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };

  gyrokinetic_p_lw->temp_ref = (struct lua_func_ctx) {
    .func_ref = temp_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };
  
  // set metatable
  luaL_getmetatable(L, GYROKINETIC_PROJECTION_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species diffusion constructor
static struct luaL_Reg gyrokinetic_projection_ctor[] = {
  { "new",  gyrokinetic_projection_lw_new },
  { 0, 0 }
};

/* ***************** */
/* Geometry methods */
/* ***************** */

// Metatable name for field input struct
#define GYROKINETIC_GEOMETRY_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Geometry"

// Lua userdata object for constructing field input
struct gyrokinetic_geometry_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_geometry gyrokinetic_geometry; // input struct to construct geometry

  // Lua registry references for mapc2p and Bmag
  struct lua_func_ctx mapc2p_ref;
  struct lua_func_ctx bmag_ref;
};

static int
gyrokinetic_geometry_lw_new(lua_State *L)
{
  struct gkyl_gyrokinetic_geometry gyrokinetic_geometry = { };

  gyrokinetic_geometry.geometry_id = GKYL_MAPC2P; 

  // Register and initialize mapc2p and bmag
  int mapc2p_ref = LUA_NOREF;
  int bmag_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "mapc2p")) {
    mapc2p_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  if (glua_tbl_get_func(L, "bmag")) {
    bmag_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  struct gyrokinetic_geometry_lw *gyrokinetic_g_lw = lua_newuserdata(L, sizeof(*gyrokinetic_g_lw));

  gyrokinetic_g_lw->magic = GYROKINETIC_GEOMETRY_DEFAULT;
  gyrokinetic_g_lw->gyrokinetic_geometry = gyrokinetic_geometry;

  gyrokinetic_g_lw->mapc2p_ref = (struct lua_func_ctx) {
    .func_ref = mapc2p_ref,
    .ndim = 0, // this will be set later
    .nret = 3,
    .L = L,
  };
  gyrokinetic_g_lw->bmag_ref = (struct lua_func_ctx) {
    .func_ref = bmag_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };
  
  // set metatable
  luaL_getmetatable(L, GYROKINETIC_GEOMETRY_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species geometry constructor
static struct luaL_Reg gyrokinetic_geometry_ctor[] = {
  { "new",  gyrokinetic_geometry_lw_new },
  { 0, 0 }
};

/* ************** */
/* Source methods */
/* ************** */

// Metatable name for field input struct
#define GYROKINETIC_SOURCE_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Source"

// Lua userdata object for constructing field input
struct gyrokinetic_source_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_source gyrokinetic_source; // input struct to construct sources

  struct gyrokinetic_projection_lw *projection_lw; // pointer to Lua projection table
  bool has_projection; // Boolean if projection exists for later initialization
};

static int
gyrokinetic_source_lw_new(lua_State *L)
{
  struct gkyl_gyrokinetic_source gyrokinetic_source = { };
  gyrokinetic_source.source_id = GKYL_PROJ_SOURCE; 

  // Fetch user input table for projection input for sources
  struct gyrokinetic_projection_lw *gyrokinetic_p_lw = 0;
  bool has_projection = false;
  with_lua_tbl_key(L, "projection") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_projection_lw *input_gyrokinetic_projection_lw = lua_touserdata(L, -1);
      if (input_gyrokinetic_projection_lw->magic == GYROKINETIC_PROJECTION_DEFAULT) {
        gyrokinetic_p_lw = input_gyrokinetic_projection_lw;
        has_projection = true;
      }
    }
  }  

  struct gyrokinetic_source_lw *gyrokinetic_src_lw = lua_newuserdata(L, sizeof(*gyrokinetic_src_lw));

  gyrokinetic_src_lw->magic = GYROKINETIC_SOURCE_DEFAULT;
  gyrokinetic_src_lw->gyrokinetic_source = gyrokinetic_source;

  gyrokinetic_src_lw->projection_lw = gyrokinetic_p_lw;
  gyrokinetic_src_lw->has_projection = has_projection;
  
  // set metatable
  luaL_getmetatable(L, GYROKINETIC_SOURCE_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species source constructor
static struct luaL_Reg gyrokinetic_source_ctor[] = {
  { "new",  gyrokinetic_source_lw_new },
  { 0, 0 }
};

/* ***************** */
/* Diffusion methods */
/* ***************** */

// Metatable name for field input struct
#define GYROKINETIC_DIFFUSION_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Diffusion"

// Lua userdata object for constructing field input
struct gyrokinetic_diffusion_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_diffusion gyrokinetic_diffusion; // input struct to construct diffusion operator
};

static int
gyrokinetic_diffusion_lw_new(lua_State *L)
{
  struct gkyl_gyrokinetic_diffusion gyrokinetic_diffusion = { };

  gyrokinetic_diffusion.order = glua_tbl_get_integer(L, "order", 2);
  gyrokinetic_diffusion.num_diff_dir = 0;
  if (glua_tbl_has_key(L, "diffDirs")) {
    with_lua_tbl_tbl(L, "diffDirs") {
      gyrokinetic_diffusion.num_diff_dir = glua_objlen(L);
      for (int d=0; d<gyrokinetic_diffusion.num_diff_dir; ++d) {
        // indexes are off by 1 between Lua and C
        gyrokinetic_diffusion.diff_dirs[d] = glua_tbl_iget_integer(L, d+1, 0)-1;
        gyrokinetic_diffusion.D[d] = glua_tbl_iget_number(L, d+1, 0.0);
      }
    }
  }

  struct gyrokinetic_diffusion_lw *gyrokinetic_d_lw = lua_newuserdata(L, sizeof(*gyrokinetic_d_lw));

  gyrokinetic_d_lw->magic = GYROKINETIC_DIFFUSION_DEFAULT;
  gyrokinetic_d_lw->gyrokinetic_diffusion = gyrokinetic_diffusion;
  
  // set metatable
  luaL_getmetatable(L, GYROKINETIC_DIFFUSION_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species diffusion constructor
static struct luaL_Reg gyrokinetic_diffusion_ctor[] = {
  { "new",  gyrokinetic_diffusion_lw_new },
  { 0, 0 }
};

/* ****************** */
/* Collisions methods */
/* ****************** */

// Metatable name for field input struct
#define GYROKINETIC_COLLISIONS_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Collisions"

// Lua userdata object for constructing field input
struct gyrokinetic_collisions_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_collisions gyrokinetic_collisions; // input struct to construct collisions
  struct lua_func_ctx init_nu_ref; // Lua registery reference to initilization function for (self) collision frequency
};

static int
gyrokinetic_collisions_lw_new(lua_State *L)
{
  struct gkyl_gyrokinetic_collisions gyrokinetic_collisions = { };

  gyrokinetic_collisions.collision_id = GKYL_LBO_COLLISIONS;  

  int init_nu_ref = LUA_NOREF;
  if (glua_tbl_get_func(L, "self_nu"))
    init_nu_ref = luaL_ref(L, LUA_REGISTRYINDEX);
  else
    return luaL_error(L, "Collisions must have an \"self_nu\" function for collision frequency!");

  struct gyrokinetic_collisions_lw *gyrokinetic_c_lw = lua_newuserdata(L, sizeof(*gyrokinetic_c_lw));

  gyrokinetic_c_lw->magic = GYROKINETIC_COLLISIONS_DEFAULT;
  gyrokinetic_c_lw->gyrokinetic_collisions = gyrokinetic_collisions;
  
  gyrokinetic_c_lw->init_nu_ref = (struct lua_func_ctx) {
    .func_ref = init_nu_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };  
  
  // set metatable
  luaL_getmetatable(L, GYROKINETIC_COLLISIONS_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species collisions constructor
static struct luaL_Reg gyrokinetic_collisions_ctor[] = {
  { "new",  gyrokinetic_collisions_lw_new },
  { 0, 0 }
};

/* *****************/
/* Species methods */
/* *****************/

// Metatable name for species input struct
#define GYROKINETIC_SPECIES_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Species"

// Lua userdata object for constructing species input
struct gyrokinetic_species_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_species gyrokinetic_species; // input struct to construct species
  int vdim; // velocity dimensions

  struct gyrokinetic_projection_lw *projection_lw; // pointer to Lua projection table
  bool has_projection; // Boolean if projection exists for later initialization
  struct gyrokinetic_collisions_lw *collisions_lw; // pointer to Lua collisions table
  bool has_collisions; // Boolean if collisions exists for later initialization
  struct gyrokinetic_diffusion_lw *diffusion_lw; // pointer to Lua diffusion table
  bool has_diffusion; // Boolean if diffusion exists for later initialization
  struct gyrokinetic_source_lw *source_lw; // pointer to Lua source table
  bool has_source; // Boolean if source exists for later initialization
};

static int
gyrokinetic_species_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_gyrokinetic_species gyrokinetic_species = { };

  gyrokinetic_species.gkmodel_id = GKYL_GK_MODEL_GEN_GEO;
  
  gyrokinetic_species.charge = glua_tbl_get_number(L, "charge", 0.0);
  gyrokinetic_species.mass = glua_tbl_get_number(L, "mass", 1.0);
  gyrokinetic_species.polarization_density = glua_tbl_get_number(L, "polarization_density", 1.0);

  with_lua_tbl_tbl(L, "cells") {
    vdim = glua_objlen(L);
    for (int d=0; d<vdim; ++d)
      gyrokinetic_species.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<vdim; ++d)
      gyrokinetic_species.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<vdim; ++d)
      gyrokinetic_species.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "diagnostics") {
    int num_diag_moments = glua_objlen(L);

    int n = 0;
    for (int i=0; i<num_diag_moments; ++i) {
      const char *mom = glua_tbl_iget_string(L, i+1, "");
      if (is_moment_name_valid(mom))
        strcpy(gyrokinetic_species.diag_moments[n++], mom);
    }
    gyrokinetic_species.num_diag_moments = n;
  }

  // Fetch user input table for projection input
  struct gyrokinetic_projection_lw *gyrokinetic_p_lw = 0;
  bool has_projection = false;
  with_lua_tbl_key(L, "projection") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_projection_lw *input_gyrokinetic_projection_lw = lua_touserdata(L, -1);
      if (input_gyrokinetic_projection_lw->magic == GYROKINETIC_PROJECTION_DEFAULT) {
        gyrokinetic_p_lw = input_gyrokinetic_projection_lw;
        has_projection = true;
      }
    }
  }  

  // Fetch user input table for collision operator 
  struct gyrokinetic_collisions_lw *gyrokinetic_c_lw = 0;
  bool has_collisions = false;
  with_lua_tbl_key(L, "collisions") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_collisions_lw *input_gyrokinetic_collisions_lw = lua_touserdata(L, -1);
      if (input_gyrokinetic_collisions_lw->magic == GYROKINETIC_COLLISIONS_DEFAULT) {
        gyrokinetic_c_lw = input_gyrokinetic_collisions_lw;
        has_collisions = true;
      }
    }
  }  

  // Fetch user input table for diffusion operator 
  struct gyrokinetic_diffusion_lw *gyrokinetic_d_lw = 0;
  bool has_diffusion = false;
  with_lua_tbl_key(L, "diffusion") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_diffusion_lw *input_gyrokinetic_diffusion_lw = lua_touserdata(L, -1);
      if (input_gyrokinetic_diffusion_lw->magic == GYROKINETIC_DIFFUSION_DEFAULT) {
        gyrokinetic_d_lw = input_gyrokinetic_diffusion_lw;
        has_diffusion = true;
      }
    }
  }  

  // Fetch user input table for sources
  struct gyrokinetic_source_lw *gyrokinetic_src_lw = 0;
  bool has_source = false;
  with_lua_tbl_key(L, "source") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_source_lw *input_gyrokinetic_source_lw = lua_touserdata(L, -1);
      if (input_gyrokinetic_source_lw->magic == GYROKINETIC_SOURCE_DEFAULT) {
        gyrokinetic_src_lw = input_gyrokinetic_source_lw;
        has_source = true;
      }
    }
  }  
  
  // Create species lua wrapper struct from inputs
  struct gyrokinetic_species_lw *gyrokinetic_s_lw = lua_newuserdata(L, sizeof(*gyrokinetic_s_lw));
  gyrokinetic_s_lw->magic = GYROKINETIC_SPECIES_DEFAULT;
  gyrokinetic_s_lw->vdim = vdim;
  gyrokinetic_s_lw->gyrokinetic_species = gyrokinetic_species;

  gyrokinetic_s_lw->projection_lw = gyrokinetic_p_lw;
  gyrokinetic_s_lw->has_projection = has_projection;

  gyrokinetic_s_lw->collisions_lw = gyrokinetic_c_lw;
  gyrokinetic_s_lw->has_collisions = has_collisions;

  gyrokinetic_s_lw->diffusion_lw = gyrokinetic_d_lw;
  gyrokinetic_s_lw->has_diffusion = has_diffusion;

  gyrokinetic_s_lw->source_lw = gyrokinetic_src_lw;
  gyrokinetic_s_lw->has_source = has_source;
   
  // set metatable
  luaL_getmetatable(L, GYROKINETIC_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Species constructor
static struct luaL_Reg gyrokinetic_species_ctor[] = {
  {"new", gyrokinetic_species_lw_new},
  {0, 0}
};

/* *****************/
/* Field methods */
/* *****************/

// Metatable name for field input struct
#define GYROKINETIC_FIELD_METATABLE_NM "GkeyllZero.App.Gyrokinetic.Field"

// Lua userdata object for constructing field input
struct gyrokinetic_field_lw {
  int magic; // this must be first element in the struct
  
  struct gkyl_gyrokinetic_field gyrokinetic_field; // input struct to construct field

  struct lua_func_ctx phi_wall_lo_ref; // Lua registery reference to biasing potential at lower wall
  bool has_phi_wall_lo; // Boolean if biasing potential at lower wall exists for later initialization
  struct lua_func_ctx phi_wall_up_ref; // Lua registery reference to biasing potential at upper wall
  bool has_phi_wall_up; // Boolean if biasing potential at upper wall exists for later initialization
};

static int
gyrokinetic_field_lw_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_gyrokinetic_field gyrokinetic_field = { };

  gyrokinetic_field.gkfield_id = GKYL_GK_FIELD_ES;  
  
  gyrokinetic_field.bmag_fac = glua_tbl_get_number(L, "bmag_fac", 1.0);

  int phi_wall_lo_ref = LUA_NOREF;
  bool has_phi_wall_lo = false;
  if (glua_tbl_get_func(L, "phi_wall_lo")) {
    phi_wall_lo_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_phi_wall_lo = true;
  }
  // By default, biasing potential at lower wall is not time dependent
  gyrokinetic_field.phi_wall_lo_evolve = glua_tbl_get_integer(L, "phi_wall_lo_evolve", false);

  int phi_wall_up_ref = LUA_NOREF;
  bool has_phi_wall_up = false;
  if (glua_tbl_get_func(L, "phi_wall_up")) {
    phi_wall_up_ref = luaL_ref(L, LUA_REGISTRYINDEX);
    has_phi_wall_up = true;
  }
  // By default, biasing potential at upper wall is not time dependent
  gyrokinetic_field.phi_wall_up_evolve = glua_tbl_get_integer(L, "phi_wall_up_evolve", false);

  struct gyrokinetic_field_lw *gyrokinetic_f_lw = lua_newuserdata(L, sizeof(*gyrokinetic_f_lw));

  gyrokinetic_f_lw->magic = GYROKINETIC_FIELD_DEFAULT;
  gyrokinetic_f_lw->gyrokinetic_field = gyrokinetic_field;

  gyrokinetic_f_lw->phi_wall_lo_ref = (struct lua_func_ctx) {
    .func_ref = phi_wall_lo_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };  
  gyrokinetic_f_lw->has_phi_wall_lo = has_phi_wall_lo;

  gyrokinetic_f_lw->phi_wall_up_ref = (struct lua_func_ctx) {
    .func_ref = phi_wall_up_ref,
    .ndim = 0, // this will be set later
    .nret = 1,
    .L = L,
  };  
  gyrokinetic_f_lw->has_phi_wall_up = has_phi_wall_up;

  // set metatable
  luaL_getmetatable(L, GYROKINETIC_FIELD_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Field constructor
static struct luaL_Reg gyrokinetic_field_ctor[] = {
  { "new",  gyrokinetic_field_lw_new },
  { 0, 0 }
};

/* *************/
/* App methods */
/* *************/

// Metatable name for top-level Gyrokinetic App
#define GYROKINETIC_APP_METATABLE_NM "GkeyllZero.App.Gyrokinetic"

// Lua userdata object for holding Gyrokinetic app and run parameters
struct gyrokinetic_app_lw {
  gkyl_gyrokinetic_app *app; // Gyokinetic app object

  // Function context for mapc2p and Bmag for geometry
  struct lua_func_ctx mapc2p_func_ctx;
  struct lua_func_ctx bmag_func_ctx;

  // Function contexts for ICs, either initialization of distribution function or Maxwellian ICs
  struct lua_func_ctx species_density_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' density 
  struct lua_func_ctx species_upar_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' upar
  struct lua_func_ctx species_temp_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' temperature

  // Function context for each species' (self) collision frequency
  struct lua_func_ctx collisions_func_ctx[GKYL_MAX_SPECIES]; 

  // Function contexts for sources, either sourcing of distribution function or Maxwellian sources
  struct lua_func_ctx src_density_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' source density 
  struct lua_func_ctx src_upar_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' source upar
  struct lua_func_ctx src_temp_func_ctx[GKYL_MAX_SPECIES]; // function context for each species' source temperature

  struct lua_func_ctx phi_wall_lo_func_ctx; // function context for biasing potential on lower wall (if present)
  struct lua_func_ctx phi_wall_up_func_ctx; // function context for biasing potential on upper wall (if present)
  
  double tstart, tend; // start and end times of simulation
  int nframe; // number of data frames to write
};

// Gets all species objects from the App table, which must on top of
// the stack. The number of species is returned and the appropriate
// pointers set in the species pointer array.
static int
get_species_inp(lua_State *L, int cdim, struct gyrokinetic_species_lw *gyrokinetic_s_lw[GKYL_MAX_SPECIES])
{
  enum { TKEY = -2, TVAL = -1};
  
  int curr = 0;
  lua_pushnil(L); // initial key is nil
  while (lua_next(L, TKEY) != 0) {
    // key at TKEY and value at TVAL
    if (lua_type(L, TVAL) == LUA_TUSERDATA) {
      struct gyrokinetic_species_lw *input_gyrokinetic_species_lw = lua_touserdata(L, TVAL);
      if (input_gyrokinetic_species_lw->magic == GYROKINETIC_SPECIES_DEFAULT) {

        if (input_gyrokinetic_species_lw->has_projection) {
          input_gyrokinetic_species_lw->projection_lw->density_ref.ndim = cdim;
          input_gyrokinetic_species_lw->projection_lw->upar_ref.ndim = cdim;
          input_gyrokinetic_species_lw->projection_lw->temp_ref.ndim = cdim;
        }
        if (input_gyrokinetic_species_lw->has_collisions) {
          input_gyrokinetic_species_lw->collisions_lw->init_nu_ref.ndim = cdim;
        }
        if (input_gyrokinetic_species_lw->has_source) {
          if (input_gyrokinetic_species_lw->source_lw->has_projection) {
            input_gyrokinetic_species_lw->source_lw->projection_lw->density_ref.ndim = cdim;
            input_gyrokinetic_species_lw->source_lw->projection_lw->upar_ref.ndim = cdim;
            input_gyrokinetic_species_lw->source_lw->projection_lw->temp_ref.ndim = cdim;            
          }
        }
        
        if (lua_type(L,TKEY) == LUA_TSTRING) {
          const char *key = lua_tolstring(L, TKEY, 0);
          strcpy(input_gyrokinetic_species_lw->gyrokinetic_species.name, key);
        }
        gyrokinetic_s_lw[curr++] = input_gyrokinetic_species_lw;
      }
    }
    lua_pop(L, 1);
  }
  return curr;
}

// Create top-level App object
static int
gyrokinetic_app_new(lua_State *L)
{
  struct gyrokinetic_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // The output prefix to use is stored in the global
  // GKYL_OUT_PREFIX. If this is not found then "g0-gyrokinetic" is used.
  const char *sim_name = "g0-gyrokinetic";
  with_lua_global(L, "GKYL_OUT_PREFIX") {
    if (lua_isstring(L, -1))
      sim_name = lua_tostring(L, -1);
  }
  
  // initialize app using table inputs (table is on top of stack)

  app_lw->tstart = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->tend = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->nframe = glua_tbl_get_integer(L, "nFrame", 1);

  struct gkyl_gk gyrokinetic = { }; // input table for app

  strcpy(gyrokinetic.name, sim_name);
  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    gyrokinetic.cdim = cdim = glua_objlen(L);
    for (int d=0; d<cdim; ++d)
      gyrokinetic.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<cdim; ++d)
      gyrokinetic.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<cdim; ++d)
      gyrokinetic.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  gyrokinetic.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  gyrokinetic.basis_type = get_basis_type(
    glua_tbl_get_string(L, "basis", "serendipity")
  );

  gyrokinetic.num_periodic_dir = 0;
  if (glua_tbl_has_key(L, "periodicDirs")) {
    with_lua_tbl_tbl(L, "periodicDirs") {
      gyrokinetic.num_periodic_dir = glua_objlen(L);
      for (int d=0; d<gyrokinetic.num_periodic_dir; ++d)
        // indexes are off by 1 between Lua and C
        gyrokinetic.periodic_dirs[d] = glua_tbl_iget_integer(L, d+1, 0)-1;
    }
  }

  // set geometry input
  with_lua_tbl_key(L, "geometry") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_geometry_lw *gyrokinetic_g_lw = lua_touserdata(L, -1);
      if (gyrokinetic_g_lw->magic == GYROKINETIC_GEOMETRY_DEFAULT) {
        // assign the App's geometry struct to user input geometry struct 
        gyrokinetic.geometry = gyrokinetic_g_lw->gyrokinetic_geometry;
        // get context for mapc2p and bmag functions
        gyrokinetic_g_lw->mapc2p_ref.ndim = cdim;
        app_lw->mapc2p_func_ctx = gyrokinetic_g_lw->mapc2p_ref;
        gyrokinetic.geometry.mapc2p = gkyl_lw_eval_cb;
        gyrokinetic.geometry.c2p_ctx = &app_lw->mapc2p_func_ctx;

        gyrokinetic_g_lw->bmag_ref.ndim = cdim;
        app_lw->bmag_func_ctx = gyrokinetic_g_lw->bmag_ref;
        gyrokinetic.geometry.bmag_func = gkyl_lw_eval_cb;
        gyrokinetic.geometry.bmag_ctx = &app_lw->bmag_func_ctx;
      }
    }
  }

  struct gyrokinetic_species_lw *gyrokinetic_s_lw[GKYL_MAX_SPECIES];
  // set all species input
  gyrokinetic.num_species = get_species_inp(L, cdim, gyrokinetic_s_lw);
  for (int s=0; s<gyrokinetic.num_species; ++s) {
    gyrokinetic.species[s] = gyrokinetic_s_lw[s]->gyrokinetic_species;
    gyrokinetic.vdim = gyrokinetic_s_lw[s]->vdim;
    
    // Note: user input projection, collisions, diffusion, and source structs
    // initialized in gyrokinetic_species_lw_new as part of species' initialization
    if (gyrokinetic_s_lw[s]->has_projection) {
      // assign the App's species object projection struct to user input projection struct 
      gyrokinetic.species[s].projection = gyrokinetic_s_lw[s]->projection_lw->gyrokinetic_projection;
      // get context for density, upar, and temperature
      app_lw->species_density_func_ctx[s] = gyrokinetic_s_lw[s]->projection_lw->density_ref;
      gyrokinetic.species[s].projection.density = gkyl_lw_eval_cb;
      gyrokinetic.species[s].projection.ctx_density = &app_lw->species_density_func_ctx[s];

      app_lw->species_upar_func_ctx[s] = gyrokinetic_s_lw[s]->projection_lw->upar_ref;
      gyrokinetic.species[s].projection.upar = gkyl_lw_eval_cb;
      gyrokinetic.species[s].projection.ctx_upar = &app_lw->species_upar_func_ctx[s];

      app_lw->species_temp_func_ctx[s] = gyrokinetic_s_lw[s]->projection_lw->temp_ref;
      gyrokinetic.species[s].projection.temp = gkyl_lw_eval_cb;
      gyrokinetic.species[s].projection.ctx_temp = &app_lw->species_temp_func_ctx[s];
    }
    if (gyrokinetic_s_lw[s]->has_collisions) {
      // assign the App's species object collision struct to user input collision struct 
      gyrokinetic.species[s].collisions = gyrokinetic_s_lw[s]->collisions_lw->gyrokinetic_collisions;
      // get context for (self) collision frequency
      app_lw->collisions_func_ctx[s] = gyrokinetic_s_lw[s]->collisions_lw->init_nu_ref;
      gyrokinetic.species[s].collisions.self_nu = gkyl_lw_eval_cb;
      gyrokinetic.species[s].collisions.ctx = &app_lw->collisions_func_ctx[s];
    }
    if (gyrokinetic_s_lw[s]->has_diffusion) {
      // assign the App's species object diffusion struct to user input diffusion struct 
      gyrokinetic.species[s].diffusion = gyrokinetic_s_lw[s]->diffusion_lw->gyrokinetic_diffusion;
    }
    if (gyrokinetic_s_lw[s]->has_source) {
      // assign the App's species object source struct to user input source struct 
      gyrokinetic.species[s].source = gyrokinetic_s_lw[s]->source_lw->gyrokinetic_source;
      // get context for distribution functions or density, upar, and temperature for Maxwellian initialization
      if (gyrokinetic_s_lw[s]->source_lw->has_projection) {
        app_lw->src_density_func_ctx[s] = gyrokinetic_s_lw[s]->source_lw->projection_lw->density_ref;
        gyrokinetic.species[s].source.projection[0].density = gkyl_lw_eval_cb;
        gyrokinetic.species[s].source.projection[0].ctx_density = &app_lw->species_density_func_ctx[s];

        app_lw->src_upar_func_ctx[s] = gyrokinetic_s_lw[s]->source_lw->projection_lw->upar_ref;
        gyrokinetic.species[s].source.projection[0].upar = gkyl_lw_eval_cb;
        gyrokinetic.species[s].source.projection[0].ctx_upar = &app_lw->src_upar_func_ctx[s];

        app_lw->src_temp_func_ctx[s] = gyrokinetic_s_lw[s]->source_lw->projection_lw->temp_ref;
        gyrokinetic.species[s].source.projection[0].temp = gkyl_lw_eval_cb;
        gyrokinetic.species[s].source.projection[0].ctx_temp = &app_lw->src_temp_func_ctx[s];
      }   
    }      
  }

  // set field input
  with_lua_tbl_key(L, "field") {
    if (lua_type(L, -1) == LUA_TUSERDATA) {
      struct gyrokinetic_field_lw *gyrokinetic_f_lw = lua_touserdata(L, -1);
      if (gyrokinetic_f_lw->magic == GYROKINETIC_FIELD_DEFAULT) {
        gyrokinetic.field = gyrokinetic_f_lw->gyrokinetic_field;

        // get context for biasing potential on lower wall if it exists
        if (gyrokinetic_f_lw->has_phi_wall_lo) {
          gyrokinetic_f_lw->phi_wall_lo_ref.ndim = cdim;
          app_lw->phi_wall_lo_func_ctx = gyrokinetic_f_lw->phi_wall_lo_ref;
          gyrokinetic.field.phi_wall_lo = gkyl_lw_eval_cb;
          gyrokinetic.field.phi_wall_lo_ctx = &app_lw->phi_wall_lo_func_ctx;
        }
        // get context for biasing potential on upper wall if it exists
        if (gyrokinetic_f_lw->has_phi_wall_up) {
          gyrokinetic_f_lw->phi_wall_up_ref.ndim = cdim;
          app_lw->phi_wall_up_func_ctx = gyrokinetic_f_lw->phi_wall_up_ref;
          gyrokinetic.field.phi_wall_up = gkyl_lw_eval_cb;
          gyrokinetic.field.phi_wall_up_ctx = &app_lw->phi_wall_up_func_ctx;
        }
      }
    }
  }
  
  app_lw->app = gkyl_gyrokinetic_app_new(&gyrokinetic);
  
  // create Lua userdata ...
  struct gyrokinetic_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct gyrokinetic_app_lw*));
  *l_app_lw = app_lw; // ... point it to the Lua app pointer

  // set metatable
  luaL_getmetatable(L, GYROKINETIC_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Apply initial conditions. (time) -> bool
static int
gyrokinetic_app_apply_ic(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double t0 = luaL_optnumber(L, 2, app_lw->tstart);
  gkyl_gyrokinetic_app_apply_ic(app_lw->app, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Apply initial conditions to species. (sidx, time) -> bool
static int
gyrokinetic_app_apply_ic_species(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  int sidx = luaL_checkinteger(L, 2);
  double t0 = luaL_optnumber(L, 3, app_lw->tstart);
  gkyl_gyrokinetic_app_apply_ic_species(app_lw->app, sidx, t0);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated moments. (tm) -> bool
static int
gyrokinetic_app_calc_integrated_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_gyrokinetic_app_calc_integrated_mom(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Compute integrated field energy. (tm) -> bool
static int
gyrokinetic_app_calc_field_energy(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  double tm = luaL_checknumber(L, 2);
  gkyl_gyrokinetic_app_calc_field_energy(app_lw->app, tm);

  lua_pushboolean(L, status);  
  return 1;
}

// Write solution (field and species) to file (time, frame) -> bool
static int
gyrokinetic_app_write(lua_State *L)
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

// Write field to file (time, frame) -> bool
static int
gyrokinetic_app_write_field(lua_State *L)
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

// Write species solution to file (sidx, time, frame) -> bool
static int
gyrokinetic_app_write_species(lua_State *L)
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

// Write diagnostic moments to file (time, frame) -> bool
static int
gyrokinetic_app_write_mom(lua_State *L)
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

// Write integrated moments to file () -> bool
static int
gyrokinetic_app_write_integrated_mom(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_write_integrated_mom(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write integrated field energy to file () -> bool
static int
gyrokinetic_app_write_field_energy(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_write_field_energy(app_lw->app);

  lua_pushboolean(L, status);  
  return 1;
}

// Write simulation statistics to JSON. () -> bool
static int
gyrokinetic_app_stat_write(lua_State *L)
{
  bool status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_stat_write(app_lw->app);

  lua_pushboolean(L, status);
  return 1;  
}

// Write data from simulation to file
static void
write_data(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) 
    gkyl_gyrokinetic_app_write(app, tcurr, iot->curr-1);
}

// Run simulation. (num_steps) -> bool. num_steps is optional
static int
gyrokinetic_app_run(lua_State *L)
{
  bool ret_status = true;

  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;
  struct gkyl_gyrokinetic_app *app = app_lw->app;

  double tcurr = app_lw->tstart;
  double tend = app_lw->tend;
  double dt = tend-tcurr;
  long num_steps = luaL_optinteger(L, 2, INT_MAX);

  int nframe = app_lw->nframe;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_gyrokinetic_app_calc_integrated_mom(app, tcurr);
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  long step = 1;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, tcurr);
    gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);    
    
    if (!status.success) {
      ret_status = false;
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);

    step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);
  gkyl_gyrokinetic_app_write_integrated_mom(app);
  gkyl_gyrokinetic_app_write_field_energy(app);

  lua_pushboolean(L, ret_status);
  return 1;
}

// Clean up memory allocated for simulation
static int
gyrokinetic_app_gc(lua_State *L)
{
  struct gyrokinetic_app_lw **l_app_lw = GKYL_CHECK_UDATA(L, GYROKINETIC_APP_METATABLE_NM);
  struct gyrokinetic_app_lw *app_lw = *l_app_lw;

  gkyl_gyrokinetic_app_release(app_lw->app);
  gkyl_free(*l_app_lw);
  
  return 0;
}

// App constructor
static struct luaL_Reg gyrokinetic_app_ctor[] = {
  { "new",  gyrokinetic_app_new },
  { 0, 0 }
};

// App methods
static struct luaL_Reg gyrokinetic_app_funcs[] = {
  { "apply_ic", gyrokinetic_app_apply_ic },
  { "apply_ic_species", gyrokinetic_app_apply_ic_species },
  { "calc_integrated_mom", gyrokinetic_app_calc_integrated_mom },
  { "calc_field_energy", gyrokinetic_app_calc_field_energy },
  { "write", gyrokinetic_app_write },
  { "write_field", gyrokinetic_app_write_field },
  { "write_species", gyrokinetic_app_write_species },
  { "write_mom", gyrokinetic_app_write_mom },
  { "write_integrated_mom", gyrokinetic_app_write_integrated_mom },
  { "write_field_energy", gyrokinetic_app_write_field_energy },
  { "stat_write", gyrokinetic_app_stat_write },
  { "run", gyrokinetic_app_run },
  { 0, 0 }
};

static void
app_openlibs(lua_State *L)
{
  // Register top-level App
  do {
    luaL_newmetatable(L, GYROKINETIC_APP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, gyrokinetic_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, gyrokinetic_app_funcs);
    
    luaL_register(L, "G0.Gyrokinetic.App", gyrokinetic_app_ctor);
    
  } while (0);

  // Register Geometry input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_GEOMETRY_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Geometry", gyrokinetic_geometry_ctor);
  } while (0);

  // Register Species input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_SPECIES_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Species", gyrokinetic_species_ctor);
  } while (0);

  // Register Projection input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_PROJECTION_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Projection", gyrokinetic_projection_ctor);
  } while (0);

  // Register Collisions input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_COLLISIONS_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Collisions", gyrokinetic_collisions_ctor);
  } while (0);

  // Register Diffusion input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_DIFFUSION_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Diffusion", gyrokinetic_diffusion_ctor);
  } while (0);

  // Register Source input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_SOURCE_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Source", gyrokinetic_source_ctor);
  } while (0);

  // Register Field input struct
  do {
    luaL_newmetatable(L, GYROKINETIC_FIELD_METATABLE_NM);
    luaL_register(L, "G0.Gyrokinetic.Field", gyrokinetic_field_ctor);
  } while (0);
}

void
gkyl_gyrokinetic_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif
