#pragma once

#ifdef GKYL_HAVE_LUA

#include <gkyl_basis.h>
#include <gkyl_util.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

// Get basis type from string
static enum gkyl_basis_type
get_basis_type(const char *bnm)
{
  if (strcmp(bnm, "serendipity") == 0)
    return GKYL_BASIS_MODAL_SERENDIPITY;
  if (strcmp(bnm, "tensor") == 0)
    return GKYL_BASIS_MODAL_TENSOR;
  if (strcmp(bnm, "hybrid") == 0)
    return GKYL_BASIS_MODAL_HYBRID;

  return GKYL_BASIS_MODAL_SERENDIPITY;
}

// Used in call back passed to the initial conditions
struct lua_func_ctx {
  int func_ref; // reference to Lua function in registery
  int ndim, nret; // dimensions of function, number of return values
  lua_State *L; // Lua state
};

/**
 * Add boundary condition flags for species into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_species_bc_types(lua_State *L);

/**
 * Add boundary condition flags for field into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_field_bc_types(lua_State *L);

/**
 * Add boundary condition flags for Poisson solver into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_poisson_bc_types(lua_State *L);

/**
* Add momentum reconstruction type flags for fluid solvers into Lua interpreter.
*
* @param L Lua state to use.
 */
void
gkyl_register_moment_scheme_types(lua_State *L);

/**
* Add Riemann problem type flags for Euler equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void
gkyl_register_euler_rp_types(lua_State *L);

/**
* Add Braginskii type flags for moment equations into Lua interpreter.
*
* @param L Lua state to use.
 */
void
gkyl_register_braginskii_types(lua_State *L);

/**
 * Add projection type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_projection_types(lua_State *L);

/**
 * Add model type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_model_types(lua_State *L);

/**
 * Add collision type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_collision_types(lua_State *L);

/**
 * Add source type flags for Vlasov species initialization into Lua interpreter.
 *
 * @param L Lua state to use.
 */
void
gkyl_register_vlasov_source_types(lua_State *L);

/**
 * Wrapper around Lua function for use in eval callbacks.
 */
void
gkyl_lw_eval_cb(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx);

#endif
