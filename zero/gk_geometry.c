#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_evalf_def.h>
#include <gkyl_util.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_derived_geo.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

static void dumpstack (lua_State *L) {
  int top=lua_gettop(L);
  for (int i=1; i <= top; i++) {
    printf("%d\t%s\t", i, luaL_typename(L,i));
    switch (lua_type(L, i)) {
      case LUA_TNUMBER:
        printf("%g\n",lua_tonumber(L,i));
        break;
      case LUA_TSTRING:
        printf("%s\n",lua_tostring(L,i));
        break;
      case LUA_TBOOLEAN:
        printf("%s\n", (lua_toboolean(L, i) ? "true" : "false"));
        break;
      case LUA_TNIL:
        printf("%s\n", "nil");
        break;
      default:
        printf("%p\n",lua_topointer(L,i));
        break;
    }
  }
}


// write out nodal coordinates 
static void
write_nodal_coordinates(const char *nm, struct gkyl_range *nrange,
  struct gkyl_array *nodes)
{
  double lower[3] = { 0.0, 0.0, 0.0 };
  double upper[3] = { 1.0, 1.0, 1.0 };
  int cells[3];
  for (int i=0; i<nrange->ndim; ++i)
    cells[i] = gkyl_range_shape(nrange, i);
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  gkyl_grid_sub_array_write(&grid, nrange, nodes, nm);
}
bool
gkyl_gk_geometry_is_cu_dev(const struct gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

void
gkyl_gk_geometry_free(const struct gkyl_ref_count *ref)
{
  struct gk_geometry *up = container_of(ref, struct gk_geometry, ref_count);
  gkyl_array_release(up->mc2p_nodal_fd);
  gkyl_array_release(up->mc2p_nodal);
  gkyl_array_release(up->mc2p);

  gkyl_array_release(up->bmag);
  gkyl_array_release(up->g_ij);
  gkyl_array_release(up->jacobgeo);
  gkyl_array_release(up->jacobgeo_inv);
  gkyl_array_release(up->gij);
  gkyl_array_release(up->b_i);
  gkyl_array_release(up->cmag);
  gkyl_array_release(up->jacobtot);
  gkyl_array_release(up->jacobtot_inv);
  gkyl_array_release(up->bmag_inv);
  gkyl_array_release(up->bmag_inv_sq);
  gkyl_array_release(up->gxxj);
  gkyl_array_release(up->gxyj);
  gkyl_array_release(up->gyyj);
  if (gkyl_gk_geometry_is_cu_dev(up)) 
    gkyl_cu_free(up->on_dev); 

  gkyl_free(up);
}

struct gk_geometry*
gkyl_gk_geometry_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_cu_dev_new(grid, range, range_ext, basis, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx);
  } 
#endif 
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = basis;
  up->range = range;
  up->range_ext = range_ext;
  up->grid = grid;
  up->nrange = gkyl_malloc(sizeof(struct gkyl_range));
  double *dzc = gkyl_malloc(3*sizeof(double));

  int poly_order = basis->poly_order;
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<grid->ndim; ++d)
      nodes[d] = grid->cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<grid->ndim; ++d)
      nodes[d] = 2*(grid->cells[d]) + 1;
  }

  gkyl_range_init_from_shape(up->nrange, up->grid->ndim, nodes);
  up->mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim*13, up->nrange->volume);
  up->mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, up->nrange->volume);
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, up->range_ext->volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis->num_basis, up->range_ext->volume);
  up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);


  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  return up;
}

void gkyl_gk_geometry_advance(struct gk_geometry* up, 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, struct gkyl_array *mc2p, 
  struct gkyl_array *bmag, struct gkyl_array *g_ij, 
  struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_inv, 
  struct gkyl_array *gij, struct gkyl_array *b_i, struct gkyl_array *cmag, struct gkyl_array *jacobtot, 
  struct gkyl_array *jacobtot_inv, struct gkyl_array *bmag_inv, struct gkyl_array *bmag_inv_sq, 
  struct gkyl_array *gxxj, struct gkyl_array *gxyj, struct gkyl_array *gyyj) 
{
  lua_State* L;
  L = luaL_newstate();
  luaL_dofile(L, "/home/akash/pppl/scripts/new_g0g2_geom/test_g2/lapd/func.lua");
  luaL_openlibs(L);

  //printf("About to call the bmag func in g0\n");
  //gkyl_eval_on_nodes *eval_bmag = gkyl_eval_on_nodes_new(up->grid, up->basis, 1, bmag_func, bmag_ctx);
  //gkyl_eval_on_nodes_advance(eval_bmag, 0.0, up->range, bmag);
  //Now project mapc2p and the FD array
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid->dx[TH_IDX],
    dpsi = up->grid->dx[PH_IDX],
    dalpha = up->grid->dx[AL_IDX];
  
  double theta_lo = up->grid->lower[TH_IDX],
    phi_lo = up->grid->lower[PH_IDX],
    alpha_lo = up->grid->lower[AL_IDX];

  double dx_fact = up->basis->poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-4;
  double delta_psi = dpsi*1e-6;
  double delta_theta = dtheta*1e-4;
  up->dzc[0] = delta_psi;
  up->dzc[1] = delta_alpha;
  up->dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};
  int cidx[3] = { 0 };
  for(int ia=up->nrange->lower[AL_IDX]; ia<=up->nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 5; ia_delta++){ // should be <5
      if(ia == up->nrange->lower[AL_IDX]){
        if(ia_delta == 1 || ia_delta == 3)
          continue; // want to use one sided stencils at edge
      }
      else if(ia == up->nrange->upper[AL_IDX]){
          if(ia_delta == 2 || ia_delta == 4)
            continue; // want to use one sided stencils at edge
      }
      else{ //interior
        if( ia_delta == 3 || ia_delta == 4)
          continue; //dont do two away
      }
      double alpha_curr = alpha_lo + ia*dalpha + modifiers[ia_delta]*delta_alpha;

      for (int ip=up->nrange->lower[PH_IDX]; ip<=up->nrange->upper[PH_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if(ip == up->nrange->lower[PH_IDX]){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if(ip == up->nrange->upper[PH_IDX]){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }
          double psi_curr = phi_lo + ip*dpsi + modifiers[ip_delta]*delta_psi;
          cidx[PH_IDX] = ip;
          for (int it=up->nrange->lower[TH_IDX]; it<=up->nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 5; // should be 5
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              if(it == up->nrange->lower[TH_IDX]){
                if(it_delta == 1 || it_delta == 3)
                  continue; // want to use one sided stencils at edge
              }
              else if(it == up->nrange->upper[TH_IDX]){
                if(it_delta == 2 || it_delta == 4)
                  continue; // want to use one sided stencils at edge
              }
              else{
                if( it_delta == 3 || it_delta == 4)
                  continue; //dont do two away
              }
              double theta_curr = theta_lo + it*dtheta + modifiers[it_delta]*delta_theta;
              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double *mc2p_fd_n = (double *) gkyl_array_fetch(up->mc2p_nodal_fd, gkyl_range_idx(up->nrange, cidx));
              double *mc2p_n = (double *) gkyl_array_fetch(up->mc2p_nodal, gkyl_range_idx(up->nrange, cidx));

              double xyz[3] = {psi_curr, alpha_curr, theta_curr};
              double XYZ[3] = {0.};
              lua_getglobal(L, "mapc2p");
              lua_pushnumber(L, psi_curr);
              lua_pushnumber(L, alpha_curr);
              lua_pushnumber(L, theta_curr);
              lua_call(L, 3, 3);
              XYZ[2] = lua_tonumber(L, -1);
              XYZ[1] = lua_tonumber(L, -2);
              XYZ[0] = lua_tonumber(L, -3);
              lua_pop(L, 3);

              // set node coordinates
              mc2p_fd_n[lidx+X_IDX] = XYZ[X_IDX];
              mc2p_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
              mc2p_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_n[X_IDX] = XYZ[X_IDX];
                mc2p_n[Y_IDX] = XYZ[Y_IDX];
                mc2p_n[Z_IDX] = XYZ[Z_IDX];
              }
            }
          }
        }
      }
    }
  }

  lua_close(L);

  int bcs[3] = {0};
  gkyl_nodal_ops_n2m(up->basis, up->grid, up->nrange, up->range, 3, up->mc2p_nodal, mc2p, bcs);

  char str1[50] = "xyz.gkyl";
  char str2[50] = "allxyz.gkyl";
  write_nodal_coordinates(str1, up->nrange, up->mc2p_nodal);
  write_nodal_coordinates(str2, up->nrange, up->mc2p_nodal_fd);

  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(up->basis, up->grid, bcs, false);
  gkyl_calc_metric_advance(mcalc, up->nrange, up->mc2p_nodal_fd, up->dzc, g_ij, up->range);

  // calculate the derived geometric quantities
  gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(up->basis, up->grid, false);
  gkyl_calc_derived_geo_advance( jcalculator, up->range, g_ij, bmag, 
    jacobgeo, jacobgeo_inv, gij, b_i, cmag, jacobtot, jacobtot_inv, 
    bmag_inv, bmag_inv_sq, gxxj, gxyj, gyyj);
}


struct gk_geometry*
gkyl_gk_geometry_acquire(const struct gk_geometry* up)
{
  gkyl_ref_count_inc(&up->ref_count);
  return (struct gk_geometry*) up;
}

void
gkyl_gk_geometry_release(const struct gk_geometry *up)
{
  gkyl_ref_count_dec(&up->ref_count);
}



