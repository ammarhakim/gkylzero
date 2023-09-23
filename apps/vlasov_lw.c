#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_lua_utils.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_lw.h>
#include <gkyl_vlasov_priv.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

// Check and fetch user-data based on metatable name
#define check_meta(L, mnm) luaL_checkudata(L, 1, mnm)

// For debugging
#define trace_stack_top(L, fnm) do { \
      printf("Inside function %s\n", fnm);                              \
      printf("--> Top of stack is %s\n", lua_typename(L, lua_type(L, -1))); \
    } while (0);

// Get basis type from string
static enum gkyl_basis_type
get_basis_type(const char *bnm)
{
  if (strcmp(bnm, "serendipity") == 0)
    return GKYL_BASIS_MODAL_SERENDIPITY;
  if (strcmp(bnm, "tensor") == 0)
    return GKYL_BASIS_MODAL_TENSOR;
  if (strcmp(bnm, "hybrid") == 0)
    return GKYL_BASIS_MODAL_HYBRID;

  return GKYL_BASIS_MODAL_SERENDIPITY;
}

/* *****************/
/* Species methods */
/* *****************/

// Metatable name for species input struct
#define VLASOV_SPECIES_METATABLE_NM "GkeyllZero.Vlasov.Species"

// Lua userdata object for constructing species input
struct vlasov_species_lw {
  struct gkyl_vlasov_species vm_species; // input struct to construct species
  int vdim; // velocity dimensions
  bool evolve; // is this species evolved?
};

static int
vm_species_new(lua_State *L)
{
  int vdim  = 0;
  struct gkyl_vlasov_species vm_species;
  
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

  struct vlasov_species_lw *vms_lw = lua_newuserdata(L, sizeof(*vms_lw));
  vms_lw->vdim = vdim;
  vms_lw->evolve = evolve;
  vms_lw->vm_species = vm_species;
  
  // set metatable
  luaL_getmetatable(L, VLASOV_SPECIES_METATABLE_NM);
  lua_setmetatable(L, -2);  
  
  return 1;
}

// Species constructor
static struct luaL_Reg vm_species_ctor[] = {
  { "new",  vm_species_new },
  { 0, 0 }
};

/* *************/
/* App methods */
/* *************/

// Metatable name for top-level Vlasov App
#define VLASOV_APP_METATABLE_NM "GkeyllZero.Vlasov.App"

// Lua userdata object for holding Vlasov app and run parameters
struct vlasov_app_lw {
  gkyl_vlasov_app *app; // Vlasov app object
  double tstart, tend; // start and end times of simulation
  int nframes; // number of data frames to write
};

// Create top-level App object
static int
vm_app_new(lua_State *L)
{
  struct vlasov_app_lw *app_lw = gkyl_malloc(sizeof(*app_lw));

  // initialize app using table inputs (table is on top of stack)

  app_lw->tstart = glua_tbl_get_number(L, "tStart", 0.0);
  app_lw->tend = glua_tbl_get_number(L, "tEnd", 1.0);
  app_lw->nframes = glua_tbl_get_integer(L, "nFrame", 1);

  struct gkyl_vm vm = { }; // input table for app

  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    vm.cdim = cdim = glua_objlen(L);
    for (int d=0; d<cdim; ++d)
      vm.cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "lower") {
    for (int d=0; d<cdim; ++d)
      vm.lower[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "upper") {
    for (int d=0; d<cdim; ++d)
      vm.upper[d] = glua_tbl_iget_number(L, d+1, 0);
  }

  vm.poly_order = glua_tbl_get_integer(L, "polyOrder", 1);

  vm.vdim = 1;
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

  app_lw->app = gkyl_vlasov_app_new(&vm);
  
  // create Lua userdata ...
  struct vlasov_app_lw **l_app_lw = lua_newuserdata(L, sizeof(struct vlasov_app_lw*));
  *l_app_lw = app_lw; // ... point it to the Lua app pointer

  // set metatable
  luaL_getmetatable(L, VLASOV_APP_METATABLE_NM);
  lua_setmetatable(L, -2);
  
  return 1;
}

// Run simulation
static int
vm_app_run(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app_lw = check_meta(L, VLASOV_APP_METATABLE_NM);
  struct vlasov_app_lw *app_lw = *l_app_lw;

  printf("Running simulation from %g to %g with %d frames!\n",
    app_lw->tstart, app_lw->tend, app_lw->nframes
  );  

  lua_pushboolean(L, status);
  return 1;
}

// Clean up memory allocated for simulation
static int
vm_app_gc(lua_State *L)
{
  struct vlasov_app_lw **l_app_lw = check_meta(L, VLASOV_APP_METATABLE_NM);
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
    
    luaL_register(L, "Plasma.App", vm_app_ctor);
    
  } while (0);

  // Register Species input struct
  do {
    luaL_newmetatable(L, VLASOV_SPECIES_METATABLE_NM);
    luaL_register(L, "Plasma.Species", vm_species_ctor);
  } while (0);  
}

void
gkyl_vlasov_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif
