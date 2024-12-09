#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

struct gk_geometry_mapc2p_advance_ctx {
  struct gk_geometry* gk_geom;

  evalf_t bmag_func;
  void *bmag_ctx;

  evalf_t mapc2p_func;
  void *mapc2p_ctx;
  struct gkyl_array *mc2p_nodal, *mc2p_nodal_fd, *mc2p;

  struct gkyl_nonuniform_position_map_info *nonuniform_map_info;
  struct gkyl_array *mc2nu_pos_nodal, *mc2nu_pos_nodal_fd, *mu2nu_pos;

  struct gkyl_range nrange;
  double dtheta, dpsi, dalpha;
  double theta_lo, phi_lo, alpha_lo;
  double delta_alpha, delta_psi, delta_theta;
  int modifiers[5];
  double dzc[3];
  int cidx[3];

  double psi_curr, alpha_curr, theta_curr;
  int ip_delta, ia_delta, it_delta;
  int ip, ia, it;
};


enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates

void
gk_geometry_mapc2p_advance_init(struct gk_geometry_mapc2p_advance_ctx *up){
  // First just do bmag
  struct gkyl_eval_on_nodes *eval_bmag = gkyl_eval_on_nodes_new(&up->gk_geom->grid, &up->gk_geom->basis, 1, up->bmag_func, up->bmag_ctx);
  gkyl_eval_on_nodes_advance(eval_bmag, 0.0, &up->gk_geom->local, up->gk_geom->bmag);
  gkyl_eval_on_nodes_release(eval_bmag);

  //Now project mapc2p and the FD array
                                
  up->dtheta = up->gk_geom->grid.dx[TH_IDX];
  up->dpsi = up->gk_geom->grid.dx[PH_IDX];
  up->dalpha = up->gk_geom->grid.dx[AL_IDX];

  up->theta_lo = up->gk_geom->grid.lower[TH_IDX] + (up->gk_geom->local.lower[TH_IDX] - up->gk_geom->global.lower[TH_IDX])*up->gk_geom->grid.dx[TH_IDX];
  up->phi_lo = up->gk_geom->grid.lower[PH_IDX] + (up->gk_geom->local.lower[PH_IDX] - up->gk_geom->global.lower[PH_IDX])*up->gk_geom->grid.dx[PH_IDX];
  up->alpha_lo = up->gk_geom->grid.lower[AL_IDX] + (up->gk_geom->local.lower[AL_IDX] - up->gk_geom->global.lower[AL_IDX])*up->gk_geom->grid.dx[AL_IDX];

  double dx_fact = up->gk_geom->basis.poly_order == 1 ? 1 : 0.5;
  up->dtheta *= dx_fact; up->dpsi *= dx_fact; up->dalpha *= dx_fact;

  // used for finite differences 
  up->delta_alpha = up->dalpha*1e-4;
  up->delta_psi = up->dpsi*1e-6;
  up->delta_theta = up->dtheta*1e-4;
  up->dzc[0] = up->delta_psi;
  up->dzc[1] = up->delta_alpha;
  up->dzc[2] = up->delta_theta;

  int modifiers_temp[5] = {0, -1, 1, -2, 2};
  for (int i = 0; i < 5; ++i) {
    up->modifiers[i] = modifiers_temp[i];
  }
  up->cidx[0] = 0;
}

void
gk_geometry_mapc2p_advance_calc(struct gk_geometry_mapc2p_advance_ctx *up)
{
  if (up->nonuniform_map_info->mapping != 0)
  {
    double coords[3] = {up->psi_curr, up->alpha_curr, up->theta_curr};
    up->nonuniform_map_info->mapping(0.0, coords, coords, up->nonuniform_map_info->ctx);
    up->psi_curr = coords[0];
    up->alpha_curr = coords[1];
    up->theta_curr = coords[2];
  }

  up->cidx[TH_IDX] = up->it;
  int lidx = 0;
  if (up->ip_delta != 0)
    lidx = 3 + 3*(up->ip_delta-1);
  if (up->ia_delta != 0)
    lidx = 15 + 3*(up->ia_delta-1);
  if (up->it_delta != 0)
    lidx = 27 + 3*(up->it_delta-1);

  double *mc2p_fd_n = (double *) gkyl_array_fetch(up->mc2p_nodal_fd, gkyl_range_idx(&up->nrange, up->cidx));
  double *mc2p_n = (double *) gkyl_array_fetch(up->mc2p_nodal, gkyl_range_idx(&up->nrange, up->cidx));

  double xyz[3] = {up->psi_curr, up->alpha_curr, up->theta_curr};
  double XYZ[3] = {0.};
  up->mapc2p_func(0.0, xyz, XYZ, up->mapc2p_ctx);

  mc2p_fd_n[lidx+X_IDX] = XYZ[X_IDX];
  mc2p_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
  mc2p_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

  if(up->ip_delta==0 && up->ia_delta==0 && up->it_delta==0){
    mc2p_n[X_IDX] = XYZ[X_IDX];
    mc2p_n[Y_IDX] = XYZ[Y_IDX];
    mc2p_n[Z_IDX] = XYZ[Z_IDX];
  }

  double *c2fa_fd_n = gkyl_array_fetch(up->mc2nu_pos_nodal_fd, gkyl_range_idx(&up->nrange, up->cidx));
  double *c2fa_n = gkyl_array_fetch(up->mc2nu_pos_nodal, gkyl_range_idx(&up->nrange, up->cidx));
  c2fa_fd_n[lidx+X_IDX] = up->psi_curr;
  c2fa_fd_n[lidx+Y_IDX] = up->alpha_curr;
  c2fa_fd_n[lidx+Z_IDX] = up->theta_curr;
  if(up->ip_delta==0 && up->ia_delta==0 && up->it_delta==0) {
    c2fa_n[X_IDX] = up->psi_curr;
    c2fa_n[Y_IDX] = up->alpha_curr;
    c2fa_n[Z_IDX] = up->theta_curr;
  }
}

void
gk_geometry_mapc2p_advance_set(struct gk_geometry_mapc2p_advance_ctx *up)
{
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->gk_geom->basis, &up->gk_geom->grid, false);
  gkyl_nodal_ops_n2m(n2m, &up->gk_geom->basis, &up->gk_geom->grid, &up->nrange, &up->gk_geom->local, 3, up->mc2p_nodal, up->mc2p);
  gkyl_nodal_ops_n2m(n2m, &up->gk_geom->basis, &up->gk_geom->grid, &up->nrange, &up->gk_geom->local, 3, up->mc2nu_pos_nodal, up->mu2nu_pos);
  gkyl_nodal_ops_release(n2m);

  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->gk_geom->basis, &up->gk_geom->grid, &up->gk_geom->global, &up->gk_geom->global_ext, &up->gk_geom->local, &up->gk_geom->local_ext, false);
  gkyl_calc_metric_advance(mcalc, &up->nrange, up->mc2p_nodal_fd, up->dzc, up->gk_geom->g_ij, up->gk_geom->dxdz, up->gk_geom->dzdx, up->gk_geom->dualmag, up->gk_geom->normals, &up->gk_geom->local);
  
  // calculate the derived geometric quantities
  struct gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->gk_geom->basis, &up->gk_geom->grid, false);
  gkyl_calc_derived_geo_advance(jcalculator, &up->gk_geom->local, up->gk_geom->g_ij, up->gk_geom->bmag, 
    up->gk_geom->jacobgeo, up->gk_geom->jacobgeo_inv, up->gk_geom->gij,  up->gk_geom->b_i,  up->gk_geom->cmag, up->gk_geom->jacobtot, up->gk_geom->jacobtot_inv, 
    up->gk_geom->bmag_inv, up->gk_geom->bmag_inv_sq,  up->gk_geom->gxxj, up->gk_geom->gxyj, up->gk_geom->gyyj, up->gk_geom->gxzj, up->gk_geom->eps2);
  gkyl_calc_derived_geo_release(jcalculator);
  gkyl_calc_metric_advance_bcart(mcalc, &up->nrange, up->gk_geom->b_i, up->gk_geom->dzdx, up->gk_geom->bcart, &up->gk_geom->local);
  gkyl_calc_metric_release(mcalc);
}

static
void gkyl_gk_geometry_mapc2p_advance(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p,
  struct gkyl_array* mc2nu_pos_nodal_fd, struct gkyl_array* mc2nu_pos_nodal, struct gkyl_array* mu2nu_pos,
  struct gkyl_nonuniform_position_map_info *nonuniform_map_info)
{
  struct gk_geometry_mapc2p_advance_ctx ctx = {
    .gk_geom = up,
    .nrange = *nrange,
    .dzc = *dzc,
    .mapc2p_func = mapc2p_func,
    .mapc2p_ctx = mapc2p_ctx,
    .bmag_func = bmag_func,
    .bmag_ctx = bmag_ctx,
    .mc2p_nodal_fd = mc2p_nodal_fd,
    .mc2p_nodal = mc2p_nodal,
    .mc2p = mc2p,
    .mc2nu_pos_nodal_fd = mc2nu_pos_nodal_fd,
    .mc2nu_pos_nodal = mc2nu_pos_nodal,
    .mu2nu_pos = mu2nu_pos,
    .nonuniform_map_info = nonuniform_map_info,
  };
  gk_geometry_mapc2p_advance_init(&ctx);

  // set alpha coordinates
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    ctx.cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 5; ia_delta++){ // should be <5
      if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) ){
        if(ia_delta == 1 || ia_delta == 3)
          continue; // want to use one sided stencils at edge
      }
      else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])){
          if(ia_delta == 2 || ia_delta == 4)
            continue; // want to use one sided stencils at edge
      }
      else{ //interior
        if( ia_delta == 3 || ia_delta == 4)
          continue; //dont do two away
      }
      ctx.alpha_curr = ctx.alpha_lo + ia*ctx.dalpha + ctx.modifiers[ia_delta]*ctx.delta_alpha;

      // set psi coordinates
      for (int ip=nrange->lower[PH_IDX]; ip<=nrange->upper[PH_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if((ip == nrange->lower[PH_IDX]) && (up->local.lower[PH_IDX]== up->global.lower[PH_IDX]) ){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if((ip == nrange->upper[PH_IDX]) && (up->local.upper[PH_IDX]== up->global.upper[PH_IDX])){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }
          ctx.psi_curr = ctx.phi_lo + ip*ctx.dpsi + ctx.modifiers[ip_delta]*ctx.delta_psi;
          ctx.cidx[PH_IDX] = ip;

          // set theta coordinates
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 5; // should be 5
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])){
                if(it_delta == 1 || it_delta == 3)
                  continue; // want to use one sided stencils at edge
              }
              else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])){
                if(it_delta == 2 || it_delta == 4)
                  continue; // want to use one sided stencils at edge
              }
              else{
                if( it_delta == 3 || it_delta == 4)
                  continue; //dont do two away
              }
              ctx.theta_curr = ctx.theta_lo + it*ctx.dtheta + ctx.modifiers[it_delta]*ctx.delta_theta;

              ctx.ip = ip;
              ctx.ia = ia;
              ctx.it = it;
              ctx.ip_delta = ip_delta;
              ctx.ia_delta = ia_delta;
              ctx.it_delta = it_delta;
              gk_geometry_mapc2p_advance_calc(&ctx);
            }
          }
        }
      }
    }
  }
  gk_geometry_mapc2p_advance_set(&ctx);
}


struct gk_geometry*
gkyl_gk_geometry_mapc2p_new(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;

  struct gkyl_range nrange;
  double dzc[3] = {0.0};

  int poly_order = up->basis.poly_order;
  int nodes[GKYL_MAX_DIM];
  if (poly_order == 1) {
    for (int d=0; d<up->grid.ndim; ++d)
      nodes[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->grid.ndim; ++d)
      nodes[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }

  gkyl_range_init_from_shape(&nrange, up->grid.ndim, nodes);
  int num_fd_nodes = 13;
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  struct gkyl_array* map_c2fa_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange.volume);
  struct gkyl_array* map_c2fa_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  up->mu2nu_pos = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->bmag_mid = gkyl_array_new(GKYL_DOUBLE, 1, 1);

  gkyl_gk_geometry_mapc2p_advance(up, &nrange, dzc, geometry_inp->mapc2p, geometry_inp->c2p_ctx,
    geometry_inp->bmag_func, geometry_inp->bmag_ctx, mc2p_nodal_fd, mc2p_nodal, up->mc2p,
    map_c2fa_nodal_fd, map_c2fa_nodal, up->mu2nu_pos, &geometry_inp->nonuniform_map_info);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  gkyl_array_release(map_c2fa_nodal_fd);
  gkyl_array_release(map_c2fa_nodal);
                   
  gkyl_array_release(mc2p_nodal_fd);
  gkyl_array_release(mc2p_nodal);

  return up;
}

