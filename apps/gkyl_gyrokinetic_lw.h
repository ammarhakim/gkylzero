#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the gyrokinetic Lua wrapper module into the Lua
 * interpreter.
 *
 * @param L Lua state
 */
void gkyl_gyrokinetic_lw_openlibs(lua_State *L);

#endif
