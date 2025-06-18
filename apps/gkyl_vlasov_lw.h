#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the Vlasov-Maxwell Lua wrapper module into the Lua
 * interpreter.
 *
 * @param L Lua state
 */
void gkyl_vlasov_lw_openlibs(lua_State *L);

#endif
