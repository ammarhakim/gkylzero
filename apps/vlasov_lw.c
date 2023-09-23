#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_lua_utils.h>
#include <gkyl_vlasov.h>
#include <gkyl_vlasov_lw.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

// Check and fetch user-data based on metatable name
#define check_meta(L, mnm) luaL_checkudata(L, 1, mnm)

/* *************/
/* App methods */
/* *************/

// Metatable name for top-level Vlasov App
#define VM_METATABLE_NM "GkeyllZero.Vlasov.App"

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
  struct vlasov_app_lw *app = gkyl_malloc(sizeof(*app));

  // initialize app using table inputs (table is on top of stack)

  app->tstart = glua_tbl_get_number(L, "tStart", 0.0);
  app->tend = glua_tbl_get_number(L, "tEnd", 1.0);
  app->nframes = glua_tbl_get_integer(L, "nFrame", 1);

  struct gkyl_vm vm; // input table for app

  int cdim = 0;
  with_lua_tbl_tbl(L, "cells") {
    cdim = glua_objlen(L);
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

  // create Lua userdata
  struct vlasov_app_lw **l_app = lua_newuserdata(L, sizeof(struct vlasov_app_lw*));
  *l_app = app;

  // set metatable
  luaL_getmetatable(L, VM_METATABLE_NM);
  lua_setmetatable(L, -2);  
  
  return 1;
}

// Run simulation
static int
vm_app_run(lua_State *L)
{
  bool status = true;

  struct vlasov_app_lw **l_app = check_meta(L, VM_METATABLE_NM);
  struct vlasov_app_lw *app = *l_app;

  printf("Running simulation from %g to %g with %d frames!\n",
    app->tstart, app->tend, app->nframes
  );  

  lua_pushboolean(L, status);
  return 1;
}

// Clean up memory allocated for simulation
static int
vm_app_gc(lua_State *L)
{
  struct vlasov_app_lw **l_app = check_meta(L, VM_METATABLE_NM);
  gkyl_free(*l_app);
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
    luaL_newmetatable(L, VM_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, vm_app_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, vm_app_funcs);
    
    luaL_register(L, "Plasma.App", vm_app_ctor);
    
  } while (0);
}

void
gkyl_vlasov_lw_openlibs(lua_State *L)
{
  app_openlibs(L);
}

#endif
