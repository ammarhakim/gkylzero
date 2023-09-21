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
  
  const char *lcode1 =
    "kvpairs = { x1 = 25.5, n1 = 2345, x2 = 32.5, n2 = 1234, bc = \"periodic\" }";
  glua_run_lua(L, lcode1, strlen(lcode1), stderr);

  const char *lcode2 =
    "nums = { 0.5, 10.5, 20.5, 30.5 }";
  glua_run_lua(L, lcode2, strlen(lcode2), stderr);

  // push table on top of stack
  lua_getglobal(L, "kvpairs");

  TEST_CHECK( 0 == glua_objlen(L) );

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

  lua_pop(L, 1);

  /* // push table on top of stack */
  /* lua_getglobal(L, "nums"); */

  /* TEST_CHECK( 4 == glua_objlen(L) ); */

  /* lua_pushinteger(L, 2); */
  /* lua_gettable(L, -2); */
  /* double v = lua_tonumber(L, -1); */
  /* lua_pop(L, 1); */
  /* printf("Num %lg\n", v); */
  
  lua_close(L);
}

void
test_2(void)
{
  lua_State *L = new_lua_State();
  
  const char *lcode1 =
    "kvpairs = { tEnd = 101.1, nums = { 11, 12, 13, v = 2222 }, names = { \"Vlasov\", \"Maxwell\"}  }";
  glua_run_lua(L, lcode1, strlen(lcode1), stderr);

  // push table on top of stack
  lua_getglobal(L, "kvpairs");
  
  TEST_CHECK( glua_tbl_has_key(L, "nums") );

  with_lua_table(L, "nums") {
    TEST_CHECK( glua_tbl_has_key(L, "v") );
    TEST_CHECK( 2222 == glua_tbl_get_integer(L, "v", 0) );

    TEST_CHECK( 3 == glua_objlen(L) );

    for (int i=1; i<=glua_objlen(L); ++i)
      TEST_CHECK( 10+i == glua_tbl_iget_integer(L, i, 0) );

    for (int i=1; i<=glua_objlen(L); ++i)
      TEST_CHECK( 10+i == glua_tbl_iget_number(L, i, 0) );
  }

  TEST_CHECK( glua_tbl_has_key(L, "names") );

  with_lua_table(L, "names") {
    TEST_CHECK( 2 == glua_objlen(L) );
    TEST_CHECK( strcmp("Vlasov", glua_tbl_iget_string(L, 1, "Einstein")) == 0 );
    TEST_CHECK( strcmp("Maxwell", glua_tbl_iget_string(L, 2, "Einstein")) == 0 );
  }  

  TEST_CHECK( glua_tbl_has_key(L, "tEnd") );
  TEST_CHECK( 101.1 == glua_tbl_get_number(L, "tEnd", 0.0) );

  lua_close(L);  
}

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },
  {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  {NULL, NULL},
};

#endif
