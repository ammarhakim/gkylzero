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

#include <gkyl_moment_lw.h>
#include <gkyl_vlasov_lw.h>
#include <gkyl_zero_lw.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#endif

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

// This is a basic Lua interpreter that sets some global values and
// then simply calls the input file specified by the -i flag. In
// general, all flags except -M (use MPI), -g (use GPU) -m (trace
// memory and -i <inp_file> are ignored.
int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
  
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  lua_State *L = luaL_newstate();
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  gkyl_zero_lw_openlibs(L);
  gkyl_vlasov_lw_openlibs(L);
  gkyl_moment_lw_openlibs(L);
  lua_gc(L, LUA_GCRESTART, -1);

#ifdef GKYL_HAVE_MPI
  lua_pushboolean(L, 1);
  lua_setglobal(L, "GKYL_HAVE_MPI");
  
  if (app_args.use_mpi) {
    lua_pushlightuserdata(L, MPI_COMM_WORLD);
    lua_setglobal(L, "GKYL_MPI_COMM");
  }
  else {
    lua_pushboolean(L, 0);
    lua_setglobal(L, "GKYL_MPI_COMM");
  }
#else
  lua_pushboolean(L, 0);
  lua_setglobal(L, "GKYL_HAVE_MPI");

  lua_pushboolean(L, 0);
  lua_setglobal(L, "GKYL_MPI_COMM");
#endif 

  if (strcmp(app_args.file_name, APP_ARGS_DEFAULT_FILE_NAME) != 0) {
    const char *inp_name = app_args.file_name;
    
    if (gkyl_check_file_exists(inp_name)) {
      // set file prefix as a global
      const char *suff = get_fname(inp_name);
      const char *suff1 = suff ? suff+1 : inp_name;
      lua_pushlstring(L, suff1, calc_output_prefix_len(suff1));
      lua_setglobal(L, "GKYL_OUT_PREFIX");
      
      int64_t sz;
      char *buff = gkyl_load_file(inp_name, &sz);
      glua_run_lua(L, buff, sz, stderr);
      gkyl_free(buff);
    }
  }
  lua_close(L);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
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
