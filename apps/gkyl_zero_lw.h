#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the zero-level Lua interface into the Lua interpreter.
 *
 * @param L Lua state
 */
void gkyl_zero_lw_openlibs(lua_State *L);

#endif
