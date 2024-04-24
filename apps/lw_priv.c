#ifdef GKYL_HAVE_LUA

#include <gkyl_app.h>
#include <gkyl_lw_priv.h>

struct reg_str_int_pair {
  const char *str;
  int val;
};

// species BC strings -> enum map
static const struct reg_str_int_pair species_bcs[] = {
  {"bcCopy", GKYL_SPECIES_COPY},
  {"bcWall", GKYL_SPECIES_REFLECT},
  {"bcReflect", GKYL_SPECIES_REFLECT},
  {"bcAbsorb", GKYL_SPECIES_ABSORB},
  {"bcNoSlip", GKYL_SPECIES_NO_SLIP},
  {"bcWedge", GKYL_SPECIES_WEDGE},
  {"bcFunc", GKYL_SPECIES_FUNC},
  {"bcFixedFunc", GKYL_SPECIES_FIXED_FUNC},
  {"bcZeroFlux", GKYL_SPECIES_ZERO_FLUX},
  {"bcGkSheath", GKYL_SPECIES_GK_SHEATH},
  {"bcRecycle", GKYL_SPECIES_RECYCLE},
  {"bcGkIWL", GKYL_SPECIES_GK_IWL},
  {0, 0},
};

// field BC strings -> enum map
static const struct reg_str_int_pair field_bcs[] = {
  { "bcCopy", GKYL_FIELD_COPY }, 
  { "bcWall", GKYL_FIELD_PEC_WALL },
  { "bcPECWall", GKYL_FIELD_PEC_WALL },
  { "bcSymWall", GKYL_FIELD_SYM_WALL },
  { "bcReservoir", GKYL_FIELD_RESERVOIR },
  { "bcWedge", GKYL_FIELD_WEDGE },
  { "bcFunc", GKYL_FIELD_FUNC },
  { 0, 0 },
};

static void
register_types(lua_State *L, const struct reg_str_int_pair types[], const char *nm)
{
  lua_getglobal(L, "G0"); // push in a table inside global G0 table
  lua_pushstring(L, nm);
  
  lua_newtable(L);
  for (int i=0; types[i].str != 0; ++i) {
    lua_pushstring(L, types[i].str);
    lua_pushinteger(L, types[i].val);
    lua_rawset(L, -3);
  }
  
  lua_rawset(L, -3);
}

void
gkyl_register_species_bc_types(lua_State *L)
{
  register_types(L, species_bcs, "SpeciesBc");
}


void
gkyl_register_field_bc_types(lua_State *L)
{
  register_types(L, field_bcs, "FieldBc");
}

void
gkyl_lw_eval_cb(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct lua_func_ctx *fr = ctx;
  lua_State *L = fr->L;

  int ndim = fr->ndim;
  int nret = fr->nret;
  lua_rawgeti(L, LUA_REGISTRYINDEX, fr->func_ref);
  lua_pushnumber(L, t);
  lua_createtable(L, GKYL_MAX_DIM, 0);

  for (int i=0; i<ndim; ++i) {
    lua_pushnumber(L, xn[i]);
    lua_rawseti(L, -2, i+1); 
  }

  if (lua_pcall(L, 2, nret, 0)) {
    const char* ret = lua_tostring(L, -1);
    luaL_error(L, "*** gkyl_lw_eval_cb ERROR: %s\n", ret);
  }

  for (int i=nret-1; i>=0; --i) { // need to fetch in reverse order
    fout[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }  
}

#endif
