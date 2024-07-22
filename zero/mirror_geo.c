#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_efit.h>
#include <gkyl_mirror_geo_priv.h>

#include <math.h>
#include <string.h>




double
mirror_plate_psi_func(double s, void *ctx){
  // uses a pointer to the plate function to get R(s), Z(s)
  // Then calculates psi(R, Z)
  // will be used by ridders later
  
  struct plate_ctx *gc = ctx;
  double RZ[2];
  if(gc->lower==true)
    gc->geo->plate_func_lower(s, RZ);
  else
    gc->geo->plate_func_upper(s, RZ);

  double R = RZ[0];
  double Z = RZ[1];

  // Now find the cell where this R and Z is
  int rzidx[2];
  rzidx[0] = fmin(gc->geo->rzlocal.lower[0] + (int) floor((R - gc->geo->rzgrid.lower[0])/gc->geo->rzgrid.dx[0]), gc->geo->rzlocal.upper[0]);
  rzidx[1] = fmin(gc->geo->rzlocal.lower[1] + (int) floor((Z - gc->geo->rzgrid.lower[1])/gc->geo->rzgrid.dx[1]), gc->geo->rzlocal.upper[1]);
  long loc = gkyl_range_idx(&gc->geo->rzlocal, rzidx);
  const double *coeffs = gkyl_array_cfetch(gc->geo->psiRZ,loc);

  double xc[2];
  gkyl_rect_grid_cell_center(&gc->geo->rzgrid, rzidx, xc);
  double xy[2];
  xy[0] = (R-xc[0])/(gc->geo->rzgrid.dx[0]*0.5);
  xy[1] = (Z-xc[1])/(gc->geo->rzgrid.dx[1]*0.5);
  double psi = gc->geo->rzbasis.eval_expand(xy, coeffs);
  return psi - gc->psi_curr;
}



// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;
  ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
  return ival;
}

gkyl_mirror_geo*
gkyl_mirror_geo_new(const struct gkyl_mirror_geo_efit_inp *inp)
{
  struct gkyl_mirror_geo *geo = gkyl_malloc(sizeof(*geo));

  geo->efit = gkyl_efit_new(inp->filepath, inp->rzpoly_order, inp->rz_basis_type, inp->fluxpoly_order, false);

  geo->plate_spec = inp->plate_spec;
  geo->plate_func_lower = inp->plate_func_lower;
  geo->plate_func_upper = inp->plate_func_upper;

  geo->rzbasis= *geo->efit->rzbasis;
  geo->rzgrid = *geo->efit->rzgrid;
  geo->psiRZ = gkyl_array_acquire(geo->efit->psizr);
  geo->psibyrRZ = gkyl_array_acquire(geo->efit->psibyrzr);
  geo->psibyr2RZ = gkyl_array_acquire(geo->efit->psibyr2zr);

  geo->num_rzbasis = geo->efit->rzbasis->num_basis;
  geo->rzlocal = *geo->efit->rzlocal;
  geo->rzlocal_ext = *geo->efit->rzlocal_ext;
  geo->fgrid = *geo->efit->fluxgrid;
  geo->fbasis = *geo->efit->fluxbasis;
  geo->frange = *geo->efit->fluxlocal;
  geo->frange_ext = *geo->efit->fluxlocal_ext;
  geo->fpoldg= gkyl_array_acquire(geo->efit->fpolflux);
  geo->qdg= gkyl_array_acquire(geo->efit->qflux);
  geo->psisep = geo->efit->sibry;
  geo->zmaxis = geo->efit->zmaxis;

  geo->root_param.eps =
    inp->root_param.eps > 0 ? inp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    inp->root_param.max_iter > 0 ? inp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    inp->quad_param.max_levels > 0 ? inp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    inp->quad_param.eps > 0 ? inp->quad_param.eps : 1e-10;

  if (geo->efit->rzbasis->poly_order == 1)
    geo->calc_roots = calc_RdR_p1;
  else if (geo->efit->rzbasis->poly_order == 2)
    geo->calc_roots = calc_RdR_p2;

  geo->stat = (struct gkyl_mirror_geo_stat) { };

  
  return geo;
}

double
gkyl_mirror_geo_integrate_psi_contour(const gkyl_mirror_geo *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_mirror_geo_R_psiZ(const gkyl_mirror_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  return R_psiZ(geo, psi, Z, nmaxroots, R, dR);
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

  gkyl_grid_sub_array_write(&grid, nrange, 0, nodes, nm);
}


void gkyl_mirror_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p)
{

  struct gkyl_mirror_geo *geo = mapc2p_ctx;
  struct gkyl_mirror_geo_grid_inp *inp = bmag_ctx;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
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

  double rclose = inp->rclose;

  int nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };

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
      // This is the convention described in Noah Mandell's Thesis Eq 5.104. comp coord y = -alpha.
      alpha_curr*=-1.0;

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
          double zmin = inp->zmin, zmax = inp->zmax;
          double darcL, arcL_curr, arcL_lo;

          mirror_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr, &zmin, &zmax, arc_memo);
          double arcL = arc_ctx.arcL_tot;
          darcL = arcL/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;

          // at the beginning of each theta loop we need to reset things
          cidx[PSI_IDX] = ip;
          arcL_curr = 0.0;
          arcL_lo = (theta_lo + M_PI)/2/M_PI*arcL;
          double ridders_min, ridders_max;
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
              arcL_curr = arcL_lo + it*darcL + modifiers[it_delta]*delta_theta*(arcL/2/M_PI);
              double theta_curr = arcL_curr*(2*M_PI/arcL) - M_PI ; 

              mirror_set_ridders(inp, &arc_ctx, psi_curr, arcL, arcL_curr, zmin, zmax, &rclose, &ridders_min, &ridders_max);

              struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
                arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
                geo->root_param.max_iter, 1e-10);
              double z_curr = res.res;
              ((gkyl_mirror_geo *)geo)->stat.nroot_cont_calls += res.nevals;
              double R[4] = { 0 }, dR[4] = { 0 };
              int nr = R_psiZ(geo, psi_curr, z_curr, 4, R, dR);
              double r_curr = choose_closest(rclose, R, R, nr);

              cidx[TH_IDX] = it;
              int lidx = 0;
              if (ip_delta != 0)
                lidx = 3 + 3*(ip_delta-1);
              if (ia_delta != 0)
                lidx = 15 + 3*(ia_delta-1);
              if (it_delta != 0)
                lidx = 27 + 3*(it_delta-1);

              if(ia_delta==0 && ip_delta==0 && it_delta==0)
                lidx = 0;

              double phi_curr = alpha_curr;
              double *mc2p_fd_n = gkyl_array_fetch(mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
              double *mc2p_n = gkyl_array_fetch(mc2p_nodal, gkyl_range_idx(nrange, cidx));
              mc2p_fd_n[lidx+X_IDX] = r_curr*cos(phi_curr);
              mc2p_fd_n[lidx+Y_IDX] = r_curr*sin(phi_curr);
              mc2p_fd_n[lidx+Z_IDX] = z_curr;

              if(ip_delta==0 && ia_delta==0 && it_delta==0){
                mc2p_n[X_IDX] = r_curr*cos(phi_curr);
                mc2p_n[Y_IDX] = r_curr*sin(phi_curr);
                mc2p_n[Z_IDX] = z_curr;
              }
            }
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, mc2p_nodal, mc2p);
  gkyl_nodal_ops_release(n2m);

  char str1[50] = "xyz";
  char str2[50] = "allxyz";
  if (inp->write_node_coord_array){
    write_nodal_coordinates(strcat(str1, inp->node_file_nm), nrange, mc2p_nodal);
    write_nodal_coordinates(strcat(str2, inp->node_file_nm), nrange, mc2p_nodal_fd);
  }
  gkyl_free(arc_memo);
}


struct gkyl_mirror_geo_stat
gkyl_mirror_geo_get_stat(const gkyl_mirror_geo *geo)
{
  return geo->stat;
}

void
gkyl_mirror_geo_release(gkyl_mirror_geo *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->psibyrRZ);
  gkyl_array_release(geo->psibyr2RZ);
  gkyl_array_release(geo->fpoldg);
  gkyl_array_release(geo->qdg);
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
