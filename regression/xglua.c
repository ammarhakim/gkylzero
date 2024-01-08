#ifdef GKYL_HAVE_LUA

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_lua_utils.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_lw.h>

static int
calc_output_prefix_len(const char *fn)
{
  const char *suff = strrchr(fn, '.');
  return strlen(fn) - (suff ? strlen(suff) : 0);
}

static const char*
get_fname(const char *fn)
{
  return strrchr(fn, '/');
}

int
main(int argc, char **argv)
{
  lua_State *L = luaL_newstate();
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  gkyl_vlasov_lw_openlibs(L);
  lua_gc(L, LUA_GCRESTART, -1);

  if (argc > 1) {
    if (gkyl_check_file_exists(argv[1])) {
      // set file prefix as a global
      const char *suff = get_fname(argv[1]);
      const char *suff1 = suff ? suff+1 : argv[1];
      lua_pushlstring(L, suff1, calc_output_prefix_len(suff1));
      lua_setglobal(L, "GKYL_OUT_PREFIX");
      
      int64_t sz;
      char *buff = gkyl_load_file(argv[1], &sz);
      glua_run_lua(L, buff, sz, stderr);
      gkyl_free(buff);
    }
  }
  lua_close(L);  
  return 0;
}

#else

#include <stdio.h>

int
main(int argc, char **argv)
{
  fprintf(stderr, "GkeyllZero built without Lua support!\n");
  return 0;
}

#endif
