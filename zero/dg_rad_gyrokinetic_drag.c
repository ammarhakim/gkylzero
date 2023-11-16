#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_rad_gyrokinetic_drag_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_const.h>

void
gkyl_rad_gyrokinetic_drag_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag  = container_of(base, struct dg_rad_gyrokinetic_drag, eqn);

  if (GKYL_IS_CU_ALLOC(rad_gyrokinetic_drag->eqn.flags))
    gkyl_cu_free(rad_gyrokinetic_drag->eqn.on_dev);
  
  gkyl_free(rad_gyrokinetic_drag);
}

void
gkyl_rad_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_drag_auxfields auxin)
{

  // TO DO:
#ifdef GKYL_HAVE_CUDA
  //if (gkyl_array_is_cu_dev(auxin.bmag_inv) && gkyl_array_is_cu_dev(auxin.nuSum) &&
  //  gkyl_array_is_cu_dev(auxin.nuPrimMomsSum) && gkyl_array_is_cu_dev(auxin.m2self)) {
  // gkyl_rad_gyrokinetic_drag_set_auxfields_cu(eqn->on_dev, auxin);
  // return;
 }
#endif

  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  rad_gyrokinetic_drag->auxfields.bmag = auxin.bmag;
  rad_gyrokinetic_drag->auxfields.vnu = auxin.vnu;
  rad_gyrokinetic_drag->auxfields.vsqnu = auxin.vsqnu;
  rad_gyrokinetic_drag->auxfields.nI = auxin.nI;
}

/* Done first pass
 * To do:
 *  Add GPU functionality
 *  Possibly remove auxfields.m2self
 */
struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
				 const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, const struct gkyl_array *bmag, const struct gkyl_array *fit_params, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    //Return to do
    //return gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(cbasis, pbasis, conf_range, pgrid, mass);
#endif
    printf("In drag new\n");
  struct dg_rad_gyrokinetic_drag* rad_gyrokinetic_drag = gkyl_malloc(sizeof(struct dg_rad_gyrokinetic_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  rad_gyrokinetic_drag->cdim = cdim;
  rad_gyrokinetic_drag->pdim = pdim;

  rad_gyrokinetic_drag->eqn.num_equations = 1;
  rad_gyrokinetic_drag->eqn.surf_term = surf;
  rad_gyrokinetic_drag->eqn.boundary_surf_term = boundary_surf;

  rad_gyrokinetic_drag->vparMax = pgrid->upper[cdim];
  rad_gyrokinetic_drag->vparMaxSq = pow(pgrid->upper[cdim],2);
  rad_gyrokinetic_drag->num_cbasis = cbasis->num_basis;

  const gkyl_dg_rad_gyrokinetic_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_rad_gyrokinetic_drag_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_rad_gyrokinetic_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  printf("set pointers\n");
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      surf_mu_kernels = ser_surf_mu_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      boundary_surf_mu_kernels = ser_boundary_surf_mu_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  rad_gyrokinetic_drag->eqn.vol_term = CK(vol_kernels, cdim, vdim, poly_order);

  rad_gyrokinetic_drag->surf[0] = CK(surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    rad_gyrokinetic_drag->surf[1] = CK(surf_mu_kernels, cdim, vdim, poly_order);

  rad_gyrokinetic_drag->boundary_surf[0] = CK(boundary_surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    rad_gyrokinetic_drag->boundary_surf[1] = CK(boundary_surf_mu_kernels, cdim, vdim, poly_order);
  printf("Chose kernel\n");
  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) assert(rad_gyrokinetic_drag->surf[i]);
  for (int i=0; i<vdim; ++i) assert(rad_gyrokinetic_drag->boundary_surf[i]);

  /*  int num_v_cells;
  if (vdim>1) {
    num_v_cells = grid->cells[cdim] * grid->cells[cdim+1];
  } else {
    num_v_cells = grid->cells[cdim];
    }*/

  rad_gyrokinetic_drag->auxfields.bmag = bmag;
  rad_gyrokinetic_drag->auxfields.nI = 0;
  rad_gyrokinetic_drag->auxfields.fit_params = fit_params;
  
  struct gkyl_range prange, prange_ext;
  const int nghost[GKYL_MAX_DIM]={0};
  printf("Set AuxFields\n");
  //gkyl_create_global_range(pdim, pgrid->cells, prange);
  gkyl_create_grid_ranges(pgrid, nghost, &prange_ext, &prange);
  printf("assign vnu_params\n");
  struct gkyl_dg_rad_vnu_params vnu_params;// = gkyl_malloc(sizeof(gkyl_dg_rad_vnu_params));
  vnu_params.cdim = cdim;
  vnu_params.vdim = vdim;
  vnu_params.bmag = bmag;
  vnu_params.fit_params = fit_params;
  vnu_params.pgrid = pgrid;
  vnu_params.conf_range = conf_range;
  printf("assigned vnu_params, vnu_params.cdim=%i\n",vnu_params.cdim);
  struct gkyl_array *vnu_temp = gkyl_array_new(GKYL_DOUBLE, pbasis->num_basis, prange.volume);
  printf("Now do eval on nodes\n");
  gkyl_eval_on_nodes *evnu = gkyl_eval_on_nodes_new(pgrid, pbasis, 1, &vnu_calc, &vnu_params);
  gkyl_eval_on_nodes_advance(evnu, 0.0, &prange, vnu_temp);//Not on ghosts?
  for (int i=0; i<32; i++){
    printf("i=%i, vnu[%i]=%e\n",i,i,vnu_temp[i]);
  }

  printf("After evnu advance\n");
  gkyl_eval_on_nodes_release(evnu);
  printf("vnu calculated\n");
  rad_gyrokinetic_drag->auxfields.vnu = vnu_temp;

  gkyl_eval_on_nodes *evsqnu = gkyl_eval_on_nodes_new(pgrid, pbasis, 1, &vsqnu_calc, &vnu_params);
  gkyl_eval_on_nodes_advance(evsqnu, 0.0, &prange, vnu_temp);//Not on ghosts?
  gkyl_eval_on_nodes_release(evsqnu);
  rad_gyrokinetic_drag->auxfields.vsqnu = vnu_temp;
  for (int i=0; i<16; i++){
    printf("i=%i, vsqnu[%i]=%e\n",i,i,vnu_temp[i]);
  }

  rad_gyrokinetic_drag->conf_range = *conf_range;

  rad_gyrokinetic_drag->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(rad_gyrokinetic_drag->eqn.flags);
  rad_gyrokinetic_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_rad_gyrokinetic_drag_free);
  rad_gyrokinetic_drag->eqn.on_dev = &rad_gyrokinetic_drag->eqn;
  
  return &rad_gyrokinetic_drag->eqn;
}

void vnu_calc(double t, const double *xn, double *fout, void *ctx){
  double vpar, mu, vmag;
  const struct gkyl_dg_rad_vnu_params *my_ctx = ctx;
  const double *fit_params = (const double*) gkyl_array_cfetch(my_ctx->fit_params, 0);
  const double a = fit_params[0];
  const double alpha = fit_params[1];
  const double beta = fit_params[2];
  const double gamma = fit_params[3];
  const double v0 = fit_params[4];
  const double scaled_v0 = v0/sqrt(GKYL_ELECTRON_MASS/(2*GKYL_ELEMENTARY_CHARGE));
  int idx[3];
  for (int i = 0; i<my_ctx->cdim; i++) 
    idx[i] = (xn[i]-my_ctx->pgrid->lower[i])/my_ctx->pgrid->dx[i];

  //printf("a=%f,alpha=%f,beta=%f,gamma=%f,v0=%f\n",a,alpha,beta,gamma,scaled_v0); all non-zero
  const double bmag = *((const double*) gkyl_array_cfetch(my_ctx->bmag,gkyl_range_idx(my_ctx->conf_range,idx)));
  const double const_mult = a*(alpha+beta)*8*bmag*sqrt(GKYL_PI)*pow(GKYL_ELEMENTARY_CHARGE,5.0/2.0)/GKYL_ELECTRON_MASS;

  
  vpar = xn[my_ctx->cdim];
  if (my_ctx->vdim>1) {
    mu = xn[my_ctx->cdim+1];
  } else {
    mu = 0;
  }
  vmag = sqrt(vpar*vpar+2*bmag*mu/GKYL_ELECTRON_MASS);
  if (vmag == 0) {
    fout[0] = 0;
  } else {
    fout[0] = 2/GKYL_PI*const_mult*pow(vmag,gamma/2+1)/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta));
  }
  printf("Bmag=%f, C=%e, xn[0]=%e, xn[1]=%e, vmag=%e, vnu=%e\n",bmag, const_mult, xn[0],xn[1],vmag,fout[0]);
}


void vsqnu_calc(double t, const double *xn, double *fout, void *ctx){
  double vpar, mu, vmag, vnu;
  const struct gkyl_dg_rad_vnu_params *my_ctx = ctx;
  const double* fit_params = (const double*) gkyl_array_cfetch(my_ctx->fit_params, 0);
  const double a = fit_params[0];
  const double alpha = fit_params[1];
  const double beta = fit_params[2];
  const double gamma = fit_params[3];
  const double v0 = fit_params[4];
  const double scaled_v0 = v0/sqrt(GKYL_ELECTRON_MASS/(2*GKYL_ELEMENTARY_CHARGE));
  int idx[3];
  for (int i = 0; i<my_ctx->cdim; i++) 
    idx[i] = (xn[i]-my_ctx->pgrid->lower[i])/my_ctx->pgrid->dx[i];
  
  const double bmag = *((const double*) gkyl_array_cfetch(my_ctx->bmag,gkyl_range_idx( my_ctx->conf_range,idx)));
  const double const_mult = a*(alpha+beta)*8*bmag*sqrt(GKYL_PI)*pow(GKYL_ELEMENTARY_CHARGE,5.0/2.0)/GKYL_ELECTRON_MASS;

  vpar = xn[my_ctx->cdim];
  if (my_ctx->vdim > 1) {
    mu = xn[my_ctx->cdim+1];
  } else {
    mu = 0;
  }
  vmag = sqrt(vpar*vpar+2*bmag*mu/GKYL_ELECTRON_MASS);
  if (vmag == 0) {
    fout[0] = 0;
  } else {
    vnu = 2/GKYL_PI*const_mult*pow(vmag,gamma/2+1)/(beta*pow(vmag/scaled_v0,-alpha)+alpha*pow(vmag/scaled_v0,beta));
    fout[0] = vmag*sqrt(mu)*GKYL_PI*pow(GKYL_ELECTRON_MASS/(2*bmag),3./2)*vnu;
  }
}


#ifndef GKYL_HAVE_CUDA
//TO DO:
struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, double mass)
{
  assert(false);
  return 0;
}

#endif
