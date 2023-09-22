#pragma once

#ifdef GKYL_HAVE_LUA

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <stdio.h>
#include <stdbool.h>

// Gets a table with given name from a table, pushing it on the
// stack. Table is popped when the scope ends
#define with_lua_tbl_tbl(L, tname)                            \
    for (bool _break = glua_tbl_get_tbl(L, tname); _break;    \
         _break = false, lua_pop(L, 1))

// This macro pushes the named global on the stack and restores the
// stack when the scope is complete.
#define with_lua_global(L, name)                                        \
    for (bool _break = (lua_getglobal(L, name), !lua_isnil(L,-1)), _isfun = lua_isfunction(L,-1); \
         _break;                                                        \
         _break = false, _isfun ? 0 : lua_pop(L, 1))

/**
 * Get length of object on top of stack. (Table size, string length
 * and memory allocated)
 *
 * @return Lenght of object on top of stack.
 */
static inline size_t glua_objlen(lua_State *L) { return lua_objlen(L, -1); }

/**
 * Check if table has specifed key in it. Table must be on top of the
 * stack
 *
 * @param L Lua state
 * @param key Key to check
 * @return true if key is present, false otherwise
 */
bool glua_tbl_has_key(lua_State *L, const char *key);

// In the following table access methods, the get_XYZ methods take a
// string as a key while the iget_XYZ takes a integer index as key.

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
double glua_tbl_iget_number(lua_State *L, long key, double def);

/**
 * Return integer from table, keyed by @a key. Table must be on top of
 * the stack. If number does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return integer corresponding to key, or def
 */
long glua_tbl_get_integer(lua_State *L, const char *key, long def);
long glua_tbl_iget_integer(lua_State *L, long key, long def);

/**
 * Return string from table, keyed by @a key. Table must be on top of
 * the stack. If string does not exist, @a def is returned.
 *
 * @param L Lua state
 * @param key Key
 * @param def Default value if key is not present in table
 * @return string corresponding to key, or def
 */
const char *glua_tbl_get_string(lua_State *L, const char *key, const char *def);
const char *glua_tbl_iget_string(lua_State *L, long key, const char *def);

/**
 * Fetches table named @a key from table on top of stack and pushes
 * that on stack. Returns false if no such table was found. You must
 * pop the table youself once you are done with it.
 *
 * @param L Lua state.
 * @param key Name of table to fetch
 * @return true if table exists, false otherwise
 */
bool glua_tbl_get_tbl(lua_State *L, const char *key);

/**
 * Fetches function named @a key from table on top of stack and pushes
 * that on stack. Returns false if no such function was found.
 *
 * @param L Lua state.
 * @param key Name of function to fetch
 * @return true if function exists, false otherwise 
 */
bool glua_tbl_get_func(lua_State *L, const char *key);

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
