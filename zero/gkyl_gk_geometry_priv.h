#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry.h>

#include <gkyl_nodal_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_derived_geo.h>



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

void  gkyl_gk_geometry_advance(struct gkyl_gk_geometry* up, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx){
  //First just do bmag
  gkyl_eval_on_nodes *eval_bmag = gkyl_eval_on_nodes_new(up->grid, up->basis, 1, bmag_func, bmag_ctx);
  gkyl_eval_on_nodes_advance(eval_bmag, 0.0, up->range, up->bmag);
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
          // set node coordinates
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
              mapc2p_func(0.0, xyz, XYZ, mapc2p_ctx);

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

  gkyl_nodal_ops_n2m(up->basis, up->grid, up->nrange, up->range, 3, up->mc2p_nodal, up->mc2p);

  char str1[50] = "xyz.gkyl";
  char str2[50] = "allxyz.gkyl";
  write_nodal_coordinates(str1, up->nrange, up->mc2p_nodal);
  write_nodal_coordinates(str2, up->nrange, up->mc2p_nodal_fd);

  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(up->basis, up->grid, false);
  gkyl_calc_metric_advance(mcalc, up->nrange, up->mc2p_nodal_fd, up->dzc, up->g_ij, up->range);

  // calculate the derived geometric quantities
  gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(up->basis, up->grid, false);
  gkyl_calc_derived_geo_advance( jcalculator, up->range, up->g_ij, up->bmag, up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj);
}

