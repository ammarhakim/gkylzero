#include <gkyl_lua_utils.h>

bool
glua_tbl_has_key(lua_State *L, const char *key)
{
  lua_getfield(L, -1, key);
  bool has_key = !lua_isnil(L, -1);
  lua_pop(L, 1);
  return has_key;  
}

double
glua_tbl_get_number(lua_State *L, const char *key, double def)
{
  double out = def;
  lua_getfield(L, -1, key);
  if (!lua_isnil(L, -1) && lua_isnumber(L, -1))
    out = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return out;
}

long
glua_tbl_get_integer(lua_State *L, const char *key, long def)
{
  long out = def;
  lua_getfield(L, -1, key);
  if (!lua_isnil(L, -1) && lua_isnumber(L, -1))
    out = lua_tointeger(L, -1);
  lua_pop(L, 1);
  return out;
}

const char *
glua_tbl_get_string(lua_State *L, const char *key, const char *def)
{
  const char *out = def;
  lua_getfield(L, -1, key);
  if (!lua_isnil(L, -1) && lua_isstring(L, -1))
    out = lua_tostring(L, -1);
  lua_pop(L, 1);
  return out;
}

int
glua_run_lua(lua_State *L, const char *str, long sz, FILE *err)
{
  if (luaL_loadbuffer(L, str, sz, "gkyl_run_lua-inp") || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    const char* ret = lua_tostring(L, -1);
    if (err)
      fprintf(err, "*** ERROR: %s\n", ret);
    return 1;
  }
  return 0;
}
