#include <acutest.h>

#ifdef GKYL_HAVE_LUA

#include <gkyl_lua_utils.h>

static lua_State *
new_lua_State(void)
{
  lua_State *L = luaL_newstate();
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  lua_gc(L, LUA_GCRESTART, -1);
  return L;
}

void
test_1(void)
{
  lua_State *L = new_lua_State();
  const char *lcode =
    "val = { x1 = 25.5, n1 = 2345, x2 = 32.5, n2 = 1234, bc = \"periodic\" }";

  glua_run_lua(L, lcode, strlen(lcode), stderr);

  // push table on top of stack
  lua_getglobal(L, "val");

  TEST_CHECK( glua_tbl_has_key(L, "x1") );
  TEST_CHECK( glua_tbl_has_key(L, "x3") == false );
  
  TEST_CHECK( 25.5 == glua_tbl_get_number(L, "x1", 0.0) );
  TEST_CHECK( 32.5 == glua_tbl_get_number(L, "x2", 0.0) );
  TEST_CHECK( 12.5 == glua_tbl_get_number(L, "x3", 12.5) );

  TEST_CHECK( glua_tbl_has_key(L, "n1") );
  TEST_CHECK( glua_tbl_has_key(L, "n3") == false );  

  TEST_CHECK( 2345 == glua_tbl_get_integer(L, "n1", 0) );
  TEST_CHECK( 1234 == glua_tbl_get_integer(L, "n2", 0) );
  TEST_CHECK( 7890 == glua_tbl_get_integer(L, "n3", 7890) );

  TEST_CHECK( glua_tbl_has_key(L, "bc") );
  TEST_CHECK( strcmp("periodic", glua_tbl_get_string(L, "bc", "N")) == 0 );
  
  lua_close(L);
}

TEST_LIST = {
  { "test_1", test_1 },
  {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  {NULL, NULL},
};

#endif
