#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>

/**
 * Load the Moments Lua wrapper module into the Lua interpreter.
 *
 * @param L Lua state
 */
void gkyl_moment_lw_openlibs(lua_State *L);

#endif
