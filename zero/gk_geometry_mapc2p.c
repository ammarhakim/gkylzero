#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_comm.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

static
void gk_geometry_mapc2p_advance(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p,
  struct gkyl_array* mc2nu_nodal, struct gkyl_array* mc2nu_pos,
  struct gkyl_position_map *position_map)
{
  // First just do bmag
  struct gkyl_eval_on_nodes *eval_bmag = gkyl_eval_on_nodes_new(&up->grid, &up->basis, 1, bmag_func, bmag_ctx);
  gkyl_eval_on_nodes_advance(eval_bmag, 0.0, &up->local, up->bmag);
  gkyl_eval_on_nodes_release(eval_bmag);

  //Now project mapc2p and the FD array
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid.dx[TH_IDX],
    dpsi = up->grid.dx[PSI_IDX],
    dalpha = up->grid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-4;
  double delta_psi = dpsi*1e-6;
  double delta_theta = dtheta*1e-4;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  position_map->constB_ctx->psi_max   = up->grid.upper[PSI_IDX];
  position_map->constB_ctx->psi_min   = up->grid.lower[PSI_IDX];
  position_map->constB_ctx->alpha_max = up->grid.upper[AL_IDX];
  position_map->constB_ctx->alpha_min = up->grid.lower[AL_IDX];
  position_map->constB_ctx->theta_max = up->grid.upper[TH_IDX];
  position_map->constB_ctx->theta_min = up->grid.lower[TH_IDX];
  position_map->constB_ctx->N_theta_boundaries = up->global.upper[TH_IDX] - up->global.lower[TH_IDX];
  gkyl_position_map_optimize(position_map);
                                
  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
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
      double alpha_curr = alpha_lo + ia*dalpha + modifiers[ia_delta]*delta_alpha;

      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) ){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX])){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }
          double psi_curr = psi_lo + ip*dpsi + modifiers[ip_delta]*delta_psi;
          cidx[PSI_IDX] = ip;
          // set node coordinates
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
              double theta_curr = theta_lo + it*dtheta + modifiers[it_delta]*delta_theta;

              position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
              position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);
              position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double *mc2p_fd_n = (double *) gkyl_array_fetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_n = (double *) gkyl_array_fetch(mc2p_nodal, gkyl_range_idx(nrange, cidx));

              double xyz[3] = {psi_curr, alpha_curr, theta_curr};
              double XYZ[3] = {0.};
              mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);

              mc2p_fd_n[lidx+X_IDX] = XYZ[X_IDX];
              mc2p_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
              mc2p_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_n[X_IDX] = XYZ[X_IDX];
                mc2p_n[Y_IDX] = XYZ[Y_IDX];
                mc2p_n[Z_IDX] = XYZ[Z_IDX];
              }

              double *mc2nu_n = gkyl_array_fetch(mc2nu_nodal, gkyl_range_idx(nrange, cidx));
              if(ip_delta==0 && ia_delta==0 && it_delta==0) {
                mc2nu_n[X_IDX] = psi_curr;
                mc2nu_n[Y_IDX] = alpha_curr;
                mc2nu_n[Z_IDX] = theta_curr;
              }
            }
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 3, mc2p_nodal, mc2p, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 3, mc2nu_nodal, mc2nu_pos, false);
  gkyl_nodal_ops_release(n2m);
}

static
void gk_geometry_mapc2p_advance_interior(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_quad_fd, struct gkyl_array *mc2p_nodal_quad, struct gkyl_array *mc2p_quad,
  struct gkyl_array *bmag_nodal, struct gkyl_position_map *position_map)
{

  //Now project mapc2p and the FD array
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid.dx[TH_IDX],
    dpsi = up->grid.dx[PSI_IDX],
    dalpha = up->grid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo = theta_lo + dels[1]*dtheta/2.0;
  psi_lo = psi_lo + dels[1]*dpsi/2.0;
  alpha_lo = alpha_lo + dels[1]*dalpha/2.0;

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  position_map->constB_ctx->psi_max   = up->grid.upper[PSI_IDX];
  position_map->constB_ctx->psi_min   = up->grid.lower[PSI_IDX];
  position_map->constB_ctx->alpha_max = up->grid.upper[AL_IDX];
  position_map->constB_ctx->alpha_min = up->grid.lower[AL_IDX];
  position_map->constB_ctx->theta_max = up->grid.upper[TH_IDX];
  position_map->constB_ctx->theta_min = up->grid.lower[TH_IDX];
  position_map->constB_ctx->N_theta_boundaries = up->global.upper[TH_IDX] - up->global.lower[TH_IDX];
  gkyl_position_map_optimize(position_map);
                                
  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 3; ia_delta++){ // interior stencil
      double alpha_curr = calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha) + modifiers[ia_delta]*delta_alpha;

      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
        int ip_delta_max = 3;// interior
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          double psi_curr = calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) + modifiers[ip_delta]*delta_psi;
          cidx[PSI_IDX] = ip;
          // set node coordinates
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 3; // interior
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              double theta_curr = calc_running_coord(theta_lo, it-nrange->lower[TH_IDX], dtheta) + modifiers[it_delta]*delta_theta;

              position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
              position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);
              position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double *mc2p_quad_fd_n = (double *) gkyl_array_fetch(mc2p_nodal_quad_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_quad_n = (double *) gkyl_array_fetch(mc2p_nodal_quad, gkyl_range_idx(nrange, cidx));
              double *bmag_n = (double *) gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));

              double xyz[3] = {psi_curr, alpha_curr, theta_curr};
              double XYZ[3] = {0.};
              double B[1] = {0.};
              mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);
              bmag_func(0.0, xyz, B, bmag_ctx);

              mc2p_quad_fd_n[lidx+X_IDX] = XYZ[X_IDX];
              mc2p_quad_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
              mc2p_quad_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_quad_n[X_IDX] = XYZ[X_IDX];
                mc2p_quad_n[Y_IDX] = XYZ[Y_IDX];
                mc2p_quad_n[Z_IDX] = XYZ[Z_IDX];
                bmag_n[0] = B[0];
              }
            }
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, false);
  gkyl_nodal_ops_n2m(n2m, &up->basis, &up->grid, nrange, &up->local, 3, mc2p_nodal_quad, mc2p_quad, true);
  gkyl_nodal_ops_release(n2m);

  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext, false);
  gkyl_calc_metric_advance_interior(mcalc, nrange, mc2p_nodal_quad_fd, dzc, up->g_ij, up->dxdz, up->dzdx, up->dualmag, up->normals, &up->local);
  gkyl_array_copy(up->g_ij_neut, up->g_ij);
  
  // calculate the derived geometric quantities
  struct gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
  gkyl_calc_derived_geo_advance(jcalculator, &up->local, up->g_ij, up->bmag, 
    up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, 
    up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj, up->gxzj, up->eps2);
  gkyl_array_copy(up->gij_neut, up->gij);
  gkyl_calc_derived_geo_release(jcalculator);
  gkyl_calc_metric_advance_bcart(mcalc, nrange, up->b_i, up->dzdx, up->bcart, &up->local);
  gkyl_calc_metric_release(mcalc);
}

void gk_geometry_mapc2p_advance_surface(struct gk_geometry* up, int dir, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *bmag_nodal,
  struct gkyl_position_map *position_map)
{

  //Now project mapc2p and the FD array
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
                                
  double dtheta = up->grid.dx[TH_IDX],
    dpsi = up->grid.dx[PSI_IDX],
    dalpha = up->grid.dx[AL_IDX];

  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo += dir == 2 ? 0.0 : dels[1]*dtheta/2.0;
  psi_lo += dir == 0 ? 0.0 : dels[1]*dpsi/2.0;
  alpha_lo += dir == 1 ? 0. : dels[1]*dalpha/2.0;

  double dx_fact = up->basis.poly_order == 1 ? 1 : 0.5;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-4;
  double delta_psi = dpsi*1e-6;
  double delta_theta = dtheta*1e-4;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  position_map->constB_ctx->psi_max   = up->grid.upper[PSI_IDX];
  position_map->constB_ctx->psi_min   = up->grid.lower[PSI_IDX];
  position_map->constB_ctx->alpha_max = up->grid.upper[AL_IDX];
  position_map->constB_ctx->alpha_min = up->grid.lower[AL_IDX];
  position_map->constB_ctx->theta_max = up->grid.upper[TH_IDX];
  position_map->constB_ctx->theta_min = up->grid.lower[TH_IDX];
  position_map->constB_ctx->N_theta_boundaries = up->global.upper[TH_IDX] - up->global.lower[TH_IDX];
  gkyl_position_map_optimize(position_map);
                                
  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    for(int ia_delta = 0; ia_delta < 5; ia_delta++){ // should be <5
      if((ia == nrange->lower[AL_IDX]) && (up->local.lower[AL_IDX]== up->global.lower[AL_IDX]) && dir==1){
        if(ia_delta == 1 || ia_delta == 3)
          continue; // want to use one sided stencils at edge
      }
      else if((ia == nrange->upper[AL_IDX])  && (up->local.upper[AL_IDX]== up->global.upper[AL_IDX])&& dir==1){
          if(ia_delta == 2 || ia_delta == 4)
            continue; // want to use one sided stencils at edge
      }
      else{ //interior
        if( ia_delta == 3 || ia_delta == 4)
          continue; //dont do two away
      }
      double alpha_curr = dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);
      alpha_curr += modifiers[ia_delta]*delta_alpha;

      for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
        int ip_delta_max = 5;// should be 5
        if(ia_delta != 0)
          ip_delta_max = 1;
        for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
          if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0){
            if(ip_delta == 1 || ip_delta == 3)
              continue; // want to use one sided stencils at edge
          }
          else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0){
            if(ip_delta == 2 || ip_delta == 4)
              continue; // want to use one sided stencils at edge
          }
          else{ // interior 
            if( ip_delta == 3 || ip_delta == 4)
              continue; //dont do two away
          }
          double psi_curr = dir == 0 ? psi_lo + ip*dpsi : calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) ;
          psi_curr += modifiers[ip_delta]*delta_psi;
          cidx[PSI_IDX] = ip;
          // set node coordinates
          for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
            int it_delta_max = 5; // should be 5
            if(ia_delta != 0 || ip_delta != 0 )
              it_delta_max = 1;
            for(int it_delta = 0; it_delta < it_delta_max; it_delta++){
              if((it == nrange->lower[TH_IDX]) && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX]) && dir==2){
                if(it_delta == 1 || it_delta == 3)
                  continue; // want to use one sided stencils at edge
              }
              else if((it == nrange->upper[TH_IDX]) && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX]) && dir==2){
                if(it_delta == 2 || it_delta == 4)
                  continue; // want to use one sided stencils at edge
              }
              else{
                if( it_delta == 3 || it_delta == 4)
                  continue; //dont do two away
              }
              double theta_curr = dir==2 ? theta_lo + it*dtheta: calc_running_coord(theta_lo, it-nrange->lower[TH_IDX], dtheta);
              theta_curr += modifiers[it_delta]*delta_theta;

              position_map->maps[0](0.0, &psi_curr,   &psi_curr,   position_map->ctxs[0]);
              position_map->maps[1](0.0, &alpha_curr, &alpha_curr, position_map->ctxs[1]);
              position_map->maps[2](0.0, &theta_curr, &theta_curr, position_map->ctxs[2]);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              double *mc2p_fd_n = (double *) gkyl_array_fetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_n = (double *) gkyl_array_fetch(mc2p_nodal, gkyl_range_idx(nrange, cidx));
              double *bmag_n = (double *) gkyl_array_fetch(bmag_nodal, gkyl_range_idx(nrange, cidx));

              double xyz[3] = {psi_curr, alpha_curr, theta_curr};
              double XYZ[3] = {0.};
              double B[1] = {0.};
              mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);
              bmag_func(0.0, xyz, B, bmag_ctx);

              mc2p_fd_n[lidx+X_IDX] = XYZ[X_IDX];
              mc2p_fd_n[lidx+Y_IDX] = XYZ[Y_IDX];
              mc2p_fd_n[lidx+Z_IDX] = XYZ[Z_IDX];

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_n[X_IDX] = XYZ[X_IDX];
                mc2p_n[Y_IDX] = XYZ[Y_IDX];
                mc2p_n[Z_IDX] = XYZ[Z_IDX];
                bmag_n[0] = B[0];
              }
            }
          }
        }
      }
    }
  }

  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext, false);
  gkyl_calc_metric_advance_surface(mcalc, dir, nrange, up->geo_surf[dir].mc2p_nodal_fd, dzc, 
    up->geo_surf[dir].bmag_nodal, up->geo_surf[dir].jacobgeo_nodal, up->geo_surf[dir].b_i_nodal,
    up->geo_surf[dir].cmag_nodal, up->geo_surf[dir].jacobtot_inv_nodal, &up->local);
  gkyl_calc_metric_release(mcalc);
  gk_geometry_surf_calc_expansions(up, dir, *nrange);
}



struct gk_geometry*
gk_geometry_mapc2p_init(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;

  struct gkyl_range nrange;
  struct gkyl_range nrange_quad;
  struct gkyl_range nrange_quad_surf[3];
  double dzc[3] = {0.0};

  int poly_order = up->basis.poly_order;

  // nodes tensor
  int num_nodes_corners[GKYL_MAX_CDIM];
  if (poly_order == 1) {
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_corners[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_corners[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }

  int num_quad_points = poly_order+1;

  int num_nodes_quad_interior[GKYL_MAX_CDIM];
  for (int d=0; d<up->grid.ndim; ++d)
    num_nodes_quad_interior[d] = gkyl_range_shape(&up->local, d)*num_quad_points;

  int num_nodes_quad_surf_in_dir[up->grid.ndim][GKYL_MAX_CDIM];
  for (int dir=0; dir<up->grid.ndim; ++dir)
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_quad_surf_in_dir[dir][d] = d == dir ? gkyl_range_shape(&up->local, d)+1 : gkyl_range_shape(&up->local, d)*num_quad_points;

  gkyl_range_init_from_shape(&nrange, up->grid.ndim, num_nodes_corners);
  gkyl_range_init_from_shape(&nrange_quad, up->grid.ndim, num_nodes_quad_interior);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gkyl_range_init_from_shape(&nrange_quad_surf[dir], up->grid.ndim, num_nodes_quad_surf_in_dir[dir]);

  // Initialize surface basis abd allocate surface geo
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, poly_order);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_nodal(up, dir, nrange_quad_surf[dir]);
    gk_geometry_surf_alloc_expansions(up, dir);
  }

  int num_fd_nodes = 13;
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange.volume);
  struct gkyl_array* mc2p_nodal_quad_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange_quad.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  struct gkyl_array* mc2p_nodal_quad = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange_quad.volume);
  struct gkyl_array *mc2p_quad = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  struct gkyl_array* mc2nu_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  up->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange_quad.volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_ghost = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
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

  // calculate mapc2p in cartesian coords at corner nodes for
  // getting cell coordinates (used only for plotting)
  gk_geometry_mapc2p_advance(up, &nrange, dzc, geometry_inp->mapc2p, geometry_inp->c2p_ctx,
    geometry_inp->bmag_func, geometry_inp->bmag_ctx, mc2p_nodal_fd, mc2p_nodal, up->mc2p,
    mc2nu_nodal, up->mc2nu_pos, geometry_inp->position_map);
  // calculate mapc2p in cartesian coords at interior nodes for
  // calculating geo quantity volume expansions 
  gk_geometry_mapc2p_advance_interior(up, &nrange_quad, dzc, geometry_inp->mapc2p, geometry_inp->c2p_ctx,
    geometry_inp->bmag_func, geometry_inp->bmag_ctx, mc2p_nodal_quad_fd, mc2p_nodal_quad, mc2p_quad, bmag_nodal, geometry_inp->position_map);
  // calculate mapc2p in cylindrical coords at surfaces
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    gk_geometry_mapc2p_advance_surface(up, dir, &nrange_quad_surf[dir], dzc, geometry_inp->mapc2p, geometry_inp->c2p_ctx,
      geometry_inp->bmag_func, geometry_inp->bmag_ctx, up->geo_surf[dir].mc2p_nodal_fd, up->geo_surf[dir].mc2p_nodal, up->geo_surf[dir].bmag_nodal, geometry_inp->position_map);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  // Release interior nodal data
  gkyl_array_release(mc2nu_nodal);
  gkyl_array_release(mc2p_nodal_quad_fd);
  gkyl_array_release(mc2p_nodal);
  gkyl_array_release(mc2p_nodal_quad);
  gkyl_array_release(mc2p_quad);
  gkyl_array_release(bmag_nodal);
  // Release surface nodal data
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gk_geometry_surf_release_nodal(up, dir);

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_mapc2p_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;


  if (geometry_inp->position_map == 0){
    geometry_inp->position_map = gkyl_position_map_new((struct gkyl_position_map_inp) {}, \
      geometry_inp->grid, geometry_inp->local, geometry_inp->local_ext, geometry_inp->local, \
      geometry_inp->local_ext, geometry_inp->basis);
    gk_geom_3d = gk_geometry_mapc2p_init(geometry_inp);
    gkyl_position_map_release(geometry_inp->position_map);
  }
  else {
    // First construct the uniform 3d geometry
    gk_geom_3d = gk_geometry_mapc2p_init(geometry_inp);
    if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
        geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
      // The array mc2nu is computed using the uniform geometry, so we need to deflate it
      // Must deflate the 3D uniform geometry in order for the allgather to work
      if(geometry_inp->grid.ndim < 3)
        gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
      else
        gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

      geometry_inp->position_map->to_optimize = true;
      gkyl_comm_array_allgather_host(geometry_inp->comm, &geometry_inp->local, \
      &geometry_inp->global, gk_geom->bmag, (struct gkyl_array*) geometry_inp->position_map->bmag_ctx->bmag);

      gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
      gkyl_gk_geometry_release(gk_geom); // release 3d geometry

      // Construct the non-uniform grid
      gk_geom_3d = gk_geometry_mapc2p_init(geometry_inp);
    }
  }
  return gk_geom_3d;
}
