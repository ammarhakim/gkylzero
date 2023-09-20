#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <stdio.h>
#include <stdbool.h>

/**
 * Check if table has specifed key in it. Table must be on top of the
 * stack
 *
 * @param L Lua state
 * @param key Key to check
 * @return true if key is present, false otherwise
 */
bool glua_tbl_has_key(lua_State *L, const char *key);

/**
 * Return number from table, keyed by @a key. Table must be on top of
 * the stack. If number does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return number corresponding to key, or def
 */
double glua_tbl_get_number(lua_State *L, const char *key, double def);

/**
 * Return integer from table, keyed by @a key. Table must be on top of
 * the stack. If number does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return number corresponding to key, or def
 */
long glua_tbl_get_integer(lua_State *L, const char *key, long def);

/**
 * Run Lua code stored in @a str buffer. The size of the buffer is
 * sz. This buffer can either be Lua byte-code or plain-text Lua code.
 * On error the error message is written to the @a err file
 * stream. Pass NULL to not print any error.
 *
 * @param L Lua state
 * @param str Buffer with Lua code (text or bytecode)
 * @param sz Size of buffer
 * @param err Error stream to write errors (if any)
 * @return 0 on success.
 */
int glua_run_lua(lua_State *L, const char *str, long sz, FILE *err);

#endif
