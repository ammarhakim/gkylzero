#ifdef GKYL_HAVE_LUA

#include <gkyl_app.h>
#include <gkyl_moment.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_vlasov.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_lw_priv.h>
#include <gkyl_util.h>

// Define options for moments of a distribution function.
static const struct gkyl_str_int_pair distribution_moms[] = {
  { "M0",                 GKYL_F_MOMENT_M0 }, // Number density.
  { "M1",                 GKYL_F_MOMENT_M1 }, // Momentum density.
  { "M2",                 GKYL_F_MOMENT_M2 }, // Kinetic energy density.
  { "M2par",              GKYL_F_MOMENT_M2PAR }, // Parallel kinetic energy density.
  { "M2perp",             GKYL_F_MOMENT_M2PERP }, // Perpendicular kinetic energy density.
  { "M2ij",               GKYL_F_MOMENT_M2IJ }, // Kinetic energy tensor..
  { "M3",                 GKYL_F_MOMENT_M3 }, // Heat flux.
  { "M3par",              GKYL_F_MOMENT_M3PAR }, // Parallel energy flux.
  { "M3perp",             GKYL_F_MOMENT_M3PERP }, // Perpendicular energy flux.
  { "M3ijk",              GKYL_F_MOMENT_M3IJK }, // Heat flux in lab frame.
  { "MaxwellianMoments",  GKYL_F_MOMENT_MAXWELLIAN }, // M0, drift speed, T/m.
  { "BiMaxwellianMoments",GKYL_F_MOMENT_BIMAXWELLIAN }, // M0, drift speed, Tpar/m, Tperp/m.
  { "LTEMoments",         GKYL_F_MOMENT_LTE }, // Maxwellian or Maxwell-Juttner moments.
  { "M0M1M2",             GKYL_F_MOMENT_M0M1M2 },  // M0, M1, M2.
  { "M0M1M2parM2perp",    GKYL_F_MOMENT_M0M1M2PARM2PERP },  // M0, M1, M2par, M2perp.
  { "HamiltonianMoments", GKYL_F_MOMENT_HAMILTONIAN },  // M0, mass*M1, H moments.
  { "M1_from_H",          GKYL_F_MOMENT_M1_FROM_H }, // dH/dv / m moment.
  { "EnergyMoment",       GKYL_F_MOMENT_ENERGY }, // H moment.
  { "M0EnergyM3",         GKYL_F_MOMENT_M0ENERGYM3 }, // M0, Energy (H) and M3 moments.
  { "Ni",                 GKYL_F_MOMENT_NI }, // M0, M1i for-vector.
  { "Tij",                GKYL_F_MOMENT_TIJ }, // Stress-energy tensor.
  { 0, 0 }
};

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
  { 0, 0 }
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
  { 0, 0 }
};

// Poisson boundary conditions -> enum map.
static const struct gkyl_str_int_pair poisson_bcs[] = {
  { "bcPeriodic", GKYL_POISSON_PERIODIC },
  { "bcDirichlet", GKYL_POISSON_DIRICHLET },
  { "bcNeumann", GKYL_POISSON_NEUMANN },
  { "bcRobin", GKYL_POISSON_ROBIN },
  { 0, 0 }
};

// Moment scheme type -> enum map.
static const struct gkyl_str_int_pair moment_scheme_type[] = {
  { "WaveProp", GKYL_MOMENT_WAVE_PROP },
  { "MP", GKYL_MOMENT_MP },
  { "KEP", GKYL_MOMENT_KEP },
  { 0, 0 }
};

// Wave limiter -> enum map.
static const struct gkyl_str_int_pair wave_limiter[] = {
  { "NoLimiter", GKYL_NO_LIMITER },
  { "MonotonizedCentered", GKYL_MONOTONIZED_CENTERED },
  { "MinMod", GKYL_MIN_MOD },
  { "SuperBee", GKYL_SUPERBEE },
  { "VanLeer", GKYL_VAN_LEER },
  { "BeamWarming", GKYL_BEAM_WARMING },
  { "Zero", GKYL_ZERO },
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

// MHD Riemann problem -> enum map.
static const struct gkyl_str_int_pair mhd_rp_type[] = {
  { "Roe", WV_MHD_RP_ROE },
  { "HLLD", WV_MHD_RP_HLLD },
  { "Lax", WV_MHD_RP_LAX },
  { 0, 0 }
};

// MHD divergence correction -> enum map.
static const struct gkyl_str_int_pair mhd_divb_type[] = {
  { "None", GKYL_MHD_DIVB_NONE },
  { "GLM", GKYL_MHD_DIVB_GLM },
  { "EightWaves", GKYL_MHD_DIVB_EIGHT_WAVES },
  { 0, 0 }
};

// Braginskii type -> enum map.
static const struct gkyl_str_int_pair braginskii_type[] = {
  { "Mag", GKYL_BRAG_MAG },
  { "Visc", GKYL_BRAG_VISC },
  { "HeatFlux", GKYL_BRAG_HEATFLUX },
  { "UnmagFull", GKYL_BRAG_UNMAG_FULL },
  { "MagFull", GKYL_BRAG_MAG_FULL },
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
  { "CanonicalPB", GKYL_MODEL_CANONICAL_PB },
  { "CanonicalPBGR", GKYL_MODEL_CANONICAL_PB_GR },
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

// Vlasov source type -> enum map.
static const struct gkyl_str_int_pair source_type[] = {
  { "None", GKYL_NO_SOURCE },
  { "Func", GKYL_FUNC_SOURCE },
  { "Proj", GKYL_PROJ_SOURCE },
  { "BoundaryFlux", GKYL_BFLUX_SOURCE },
  { 0, 0 }
};

// Gyrokinetic FEM boundary conditions -> enum map.
static const struct gkyl_str_int_pair parproj_type[] = {
  { "None", GKYL_FEM_PARPROJ_NONE },
  { "Periodic", GKYL_FEM_PARPROJ_PERIODIC },
  { "Dirichlet", GKYL_FEM_PARPROJ_DIRICHLET },
  { 0, 0 }
};

// Gyrokinetic geometry type -> enum map.
static const struct gkyl_str_int_pair geometry_type[] = {
  { "Tokamak", GKYL_TOKAMAK },
  { "Mirror", GKYL_MIRROR },
  { "MapC2P", GKYL_MAPC2P },
  { "FromFile", GKYL_GEOMETRY_FROMFILE },
  { 0, 0 }
};

// Gyrokinetic position map type -> enum map.
static const struct gkyl_str_int_pair position_map_type[] = {
  { "UserInput", GKYL_PMAP_USER_INPUT },
  { "ConstantPolynomial", GKYL_PMAP_CONSTANT_DB_POLYNOMIAL },
  { "ConstantNumeric", GKYL_PMAP_CONSTANT_DB_NUMERIC },
  { 0, 0 }
};

// Gyrokinetic field type -> enum map.
static const struct gkyl_str_int_pair gk_field_type[] = {
  { "Electrostatic", GKYL_GK_FIELD_ES },
  { "Boltzmann", GKYL_GK_FIELD_BOLTZMANN },
  { "Adiabatic", GKYL_GK_FIELD_ADIABATIC },
  { "ElectrostaticIWL", GKYL_GK_FIELD_ES_IWL },
  { "Electromagnetic", GKYL_GK_FIELD_EM },
  { 0, 0 }
};

// Gyrokinetic radiation type -> enum map.
static const struct gkyl_str_int_pair gk_radiation_type[] = {
  { "None", GKYL_NO_RADIATION },
  { "GKRadiation", GKYL_GK_RADIATION },
  { "VMComptonRadiation", GKYL_VM_COMPTON_RADIATION },
  { 0, 0 }
};

// Gyrokinetic radiation Te model type -> enum map.
static const struct gkyl_str_int_pair gk_radiation_te_type[] = {
  { "Conservative", GKYL_VARY_TE_CONSERVATIVE },
  { "Aggressive", GKYL_VARY_TE_AGGRESSIVE },
  { "Const", GKYL_CONST_TE },
  { 0, 0 }
};

// Gyrokinetic reaction type -> enum map.
static const struct gkyl_str_int_pair gk_react_type[] = {
  { "None", GKYL_NO_REACT },
  { "Ionization", GKYL_REACT_IZ },
  { "ChargeExchange", GKYL_REACT_CX },
  { "Recombination", GKYL_REACT_RECOMB },
  { 0, 0 }
};

// Gyrokinetic ion type -> enum map.
static const struct gkyl_str_int_pair gk_ion_type[] = {
  { "Hydrogen", GKYL_ION_H },
  { "Deuterium", GKYL_ION_D },
  { "Helium", GKYL_ION_HE },
  { "Lithium", GKYL_ION_LI },
  { "Beryllium", GKYL_ION_BE },
  { "Boron", GKYL_ION_B },
  { "Carbon", GKYL_ION_C },
  { "Nitrogen", GKYL_ION_N },
  { "Oxygen", GKYL_ION_O },
  { "Neon", GKYL_ION_NE },
  { "Argon", GKYL_ION_AR },
  { 0, 0 }
};

// Gyrokinetic self-reaction type -> enum map.
static const struct gkyl_str_int_pair gk_react_self_type[] = {
  { "Electron", GKYL_SELF_ELC },
  { "Ion", GKYL_SELF_ION },
  { "Donor", GKYL_SELF_DONOR },
  { "Receiver", GKYL_SELF_RECVR },
  { "Partner", GKYL_SELF_PARTNER },
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
gkyl_register_distribution_moment_types(lua_State *L)
{
  register_types(L, distribution_moms, "DistributionMoment");
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
gkyl_register_poisson_bc_types(lua_State *L)
{
  register_types(L, poisson_bcs, "PoissonBc");
}

void
gkyl_register_moment_scheme_types(lua_State *L)
{
  register_types(L, moment_scheme_type, "SchemeType");
}

void
gkyl_register_wave_limiter_types(lua_State *L)
{
  register_types(L, wave_limiter, "WaveLimiter");
}

void
gkyl_register_euler_rp_types(lua_State *L)
{
  register_types(L, euler_rp_type, "EulerRP");
}

void
gkyl_register_mhd_rp_types(lua_State *L)
{
  register_types(L, mhd_rp_type, "MHDRP");
}

void
gkyl_register_mhd_divb_types(lua_State *L)
{
  register_types(L, mhd_divb_type, "DivB");
}

void
gkyl_register_braginskii_types(lua_State *L)
{
  register_types(L, braginskii_type, "Braginskii");
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
gkyl_register_vlasov_source_types(lua_State *L)
{
  register_types(L, source_type, "Source");
}

void
gkyl_register_gyrokinetic_fem_bc_types(lua_State *L)
{
  register_types(L, parproj_type, "ParProjBc");
}

void
gkyl_register_gyrokinetic_geometry_types(lua_State *L)
{
  register_types(L, geometry_type, "Geometry");
}

void
gkyl_register_gyrokinetic_position_map_types(lua_State *L)
{
  register_types(L, position_map_type, "PositionMap");
}

void
gkyl_register_gyrokinetic_field_types(lua_State *L)
{
  register_types(L, gk_field_type, "GKField");
}

void
gkyl_register_gyrokinetic_radiation_types(lua_State *L)
{
  register_types(L, gk_radiation_type, "Radiation");
}

void
gkyl_register_gyrokinetic_radiation_Te_types(lua_State *L)
{
  register_types(L, gk_radiation_te_type, "TeMinModel");
}

void
gkyl_register_gyrokinetic_reaction_types(lua_State *L)
{
  register_types(L, gk_react_type, "Reaction");
}

void
gkyl_register_gyrokinetic_ion_types(lua_State *L)
{
  register_types(L, gk_ion_type, "Ion");
}

void
gkyl_register_gyrokinetic_self_reaction_types(lua_State *L)
{
  register_types(L, gk_react_self_type, "Self");
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
