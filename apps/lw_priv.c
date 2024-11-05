#ifdef GKYL_HAVE_LUA

#include <gkyl_app.h>
#include <gkyl_moment.h>
#include <gkyl_wv_euler.h>
#include <gkyl_vlasov.h>
#include <gkyl_lw_priv.h>
#include <gkyl_util.h>

// Species boundary conditions -> enum map.
static const struct gkyl_str_int_pair species_bcs[] = {
  { "bcCopy", GKYL_SPECIES_COPY },
  { "bcWall", GKYL_SPECIES_REFLECT },
  { "bcReflect", GKYL_SPECIES_REFLECT },
  { "bcAbsorb", GKYL_SPECIES_ABSORB },
  { "bcNoSlip", GKYL_SPECIES_NO_SLIP },
  { "bcWedge", GKYL_SPECIES_WEDGE },
  { "bcFunc", GKYL_SPECIES_FUNC },
  { "bcFixedFunc", GKYL_SPECIES_FIXED_FUNC },
  { "bcZeroFlux", GKYL_SPECIES_ZERO_FLUX },
  { "bcGkSheath", GKYL_SPECIES_GK_SHEATH },
  { "bcRecycle", GKYL_SPECIES_RECYCLE },
  { "bcGkIWL", GKYL_SPECIES_GK_IWL },
  { 0, 0 },
};

// Field boundary conditions -> enum map.
static const struct gkyl_str_int_pair field_bcs[] = {
  { "bcCopy", GKYL_FIELD_COPY }, 
  { "bcWall", GKYL_FIELD_PEC_WALL },
  { "bcPECWall", GKYL_FIELD_PEC_WALL },
  { "bcSymWall", GKYL_FIELD_SYM_WALL },
  { "bcReservoir", GKYL_FIELD_RESERVOIR },
  { "bcWedge", GKYL_FIELD_WEDGE },
  { "bcFunc", GKYL_FIELD_FUNC },
  { 0, 0 },
};

// Moment scheme type -> enum map.
static const struct gkyl_str_int_pair moment_scheme_type[] = {
  { "WaveProp", GKYL_MOMENT_WAVE_PROP },
  { "MP", GKYL_MOMENT_MP },
  { "KEP", GKYL_MOMENT_KEP },
  { 0, 0 }
};

// Euler Riemann problem -> enum map.
static const struct gkyl_str_int_pair euler_rp_type[] = {
  { "Roe", WV_EULER_RP_ROE },
  { "HLLC", WV_EULER_RP_HLLC },
  { "Lax", WV_EULER_RP_LAX },
  { "HLL", WV_EULER_RP_HLL },
  { 0, 0 }
};

// Vlasov projection type -> enum map.
static const struct gkyl_str_int_pair projection_type[] = {
  { "Func", GKYL_PROJ_FUNC },
  { "MaxwellianPrimitive", GKYL_PROJ_MAXWELLIAN_PRIM },
  { "MaxwellianLab", GKYL_PROJ_MAXWELLIAN_LAB },
  { "BiMaxwellian", GKYL_PROJ_BIMAXWELLIAN },
  { "LTE", GKYL_PROJ_VLASOV_LTE },
  { 0, 0 }
};

// Vlasov model type -> enum map.
static const struct gkyl_str_int_pair model_type[] = {
  { "Default", GKYL_MODEL_DEFAULT },
  { "SR", GKYL_MODEL_SR },
  { "GeneralGeometry", GKYL_MODEL_GEN_GEO },
  { "CanonicalPB", GKYL_MODEL_CANONICAL_PB },
  { 0, 0 }
};

// Vlasov collision type -> enum map.
static const struct gkyl_str_int_pair collision_type[] = {
  { "None", GKYL_NO_COLLISIONS },
  { "BGK", GKYL_BGK_COLLISIONS },
  { "LBO", GKYL_LBO_COLLISIONS },
  { "FPO", GKYL_FPO_COLLISIONS },
  { 0, 0 }
};

static void
register_types(lua_State *L, const struct gkyl_str_int_pair types[], const char *nm)
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
gkyl_register_moment_scheme_types(lua_State *L)
{
  register_types(L, moment_scheme_type, "SchemeType");
}

void
gkyl_register_euler_rp_types(lua_State *L)
{
  register_types(L, euler_rp_type, "EulerRP");
}

void
gkyl_register_vlasov_projection_types(lua_State *L)
{
  register_types(L, projection_type, "Projection");
}

void
gkyl_register_vlasov_model_types(lua_State *L)
{
  register_types(L, model_type, "Model");
}

void
gkyl_register_vlasov_collision_types(lua_State *L)
{
  register_types(L, collision_type, "Collisions");
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
