/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion.h>    
#include <gkyl_dg_diffusion_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_diffusion_set_auxfields_cu_kernel(const struct gkyl_dg_eqn* eqn, const struct gkyl_array* D)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  diffusion->auxfields.D = D;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_auxfields auxin)
{
  gkyl_diffusion_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.D->on_dev);
}

__global__ void static
dg_diffusion_set_cu_dev_ptrs(struct dg_diffusion *diffusion, enum gkyl_basis_type b_type, int cdim, int vdim, int poly_order, enum gkyl_diffusion_id diffusion_id, int diff_order, int diffdirs_linidx)
{
  diffusion->auxfields.D = 0; 

  const gkyl_dg_diffusion_vol_kern_list *vol_kernels;
  const gkyl_dg_diffusion_surf_kern_list *surfx_kernels;
  const gkyl_dg_diffusion_surf_kern_list *surfy_kernels;
  const gkyl_dg_diffusion_surf_kern_list *surfz_kernels;
  const gkyl_dg_diffusion_boundary_surf_kern_list *boundary_surfx_kernels;
  const gkyl_dg_diffusion_boundary_surf_kern_list *boundary_surfy_kernels;
  const gkyl_dg_diffusion_boundary_surf_kern_list *boundary_surfz_kernels;

  if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_VLASOV)) {
    switch (cbasis->b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        vol_kernels            = diffusion->const_coeff? ser_vol_kernels_constcoeff                   : ser_vol_kernels_varcoeff                  ;
        surfx_kernels          = diffusion->const_coeff? ser_vlasov_surfx_kernels_constcoeff          : ser_vlasov_surfx_kernels_varcoeff         ;
        surfy_kernels          = diffusion->const_coeff? ser_vlasov_surfy_kernels_constcoeff          : ser_vlasov_surfy_kernels_varcoeff         ;
        surfz_kernels          = diffusion->const_coeff? ser_vlasov_surfz_kernels_constcoeff          : ser_vlasov_surfz_kernels_varcoeff         ;
        boundary_surfx_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfx_kernels_constcoeff : ser_vlasov_boundary_surfx_kernels_varcoeff;
        boundary_surfy_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfy_kernels_constcoeff : ser_vlasov_boundary_surfy_kernels_varcoeff;
        boundary_surfz_kernels = diffusion->const_coeff? ser_vlasov_boundary_surfz_kernels_constcoeff : ser_vlasov_boundary_surfz_kernels_varcoeff;
        break;

      default:
        assert(false);
        break;
    }
  } else if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_GYROKINETIC) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_GYROKINETIC)) {
    switch (cbasis->b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        vol_kernels            = diffusion->const_coeff? ser_vol_kernels_constcoeff                        : ser_vol_kernels_varcoeff                       ;
        surfx_kernels          = diffusion->const_coeff? ser_gyrokinetic_surfx_kernels_constcoeff          : ser_gyrokinetic_surfx_kernels_varcoeff         ;
        surfy_kernels          = diffusion->const_coeff? ser_gyrokinetic_surfy_kernels_constcoeff          : ser_gyrokinetic_surfy_kernels_varcoeff         ;
        surfz_kernels          = diffusion->const_coeff? ser_gyrokinetic_surfz_kernels_constcoeff          : ser_gyrokinetic_surfz_kernels_varcoeff         ;
        boundary_surfx_kernels = diffusion->const_coeff? ser_gyrokinetic_boundary_surfx_kernels_constcoeff : ser_gyrokinetic_boundary_surfx_kernels_varcoeff;
        boundary_surfy_kernels = diffusion->const_coeff? ser_gyrokinetic_boundary_surfy_kernels_constcoeff : ser_gyrokinetic_boundary_surfy_kernels_varcoeff;
        boundary_surfz_kernels = diffusion->const_coeff? ser_gyrokinetic_boundary_surfz_kernels_constcoeff : ser_gyrokinetic_boundary_surfz_kernels_varcoeff;
        break;

      default:
        assert(false);
        break;
    }
  } else {
    switch (cbasis->b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        vol_kernels            = diffusion->const_coeff? ser_vol_kernels_constcoeff            : ser_vol_kernels_varcoeff           ;
        surfx_kernels          = diffusion->const_coeff? ser_surfx_kernels_constcoeff          : ser_surfx_kernels_varcoeff         ;
        surfy_kernels          = diffusion->const_coeff? ser_surfy_kernels_constcoeff          : ser_surfy_kernels_varcoeff         ;
        surfz_kernels          = diffusion->const_coeff? ser_surfz_kernels_constcoeff          : ser_surfz_kernels_varcoeff         ;
        boundary_surfx_kernels = diffusion->const_coeff? ser_boundary_surfx_kernels_constcoeff : ser_boundary_surfx_kernels_varcoeff;
        boundary_surfy_kernels = diffusion->const_coeff? ser_boundary_surfy_kernels_constcoeff : ser_boundary_surfy_kernels_varcoeff;
        boundary_surfz_kernels = diffusion->const_coeff? ser_boundary_surfz_kernels_constcoeff : ser_boundary_surfz_kernels_varcoeff;
        break;

      default:
        assert(false);
        break;
    }
  }

  diffusion->eqn.num_equations = diffusion->num_equations;
  diffusion->eqn.surf_term = surf;
  diffusion->eqn.boundary_surf_term = boundary_surf;

  diffusion->eqn.vol_term = CKVOL(vol_kernels, cdim, diff_order, poly_order, diffdirs_linidx);

  diffusion->surf[0] = CKSURF(surfx_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>1)
    diffusion->surf[1] = CKSURF(surfy_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>2)
    diffusion->surf[2] = CKSURF(surfz_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);

  diffusion->boundary_surf[0] = CKSURF(boundary_surfx_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>1)
    diffusion->boundary_surf[1] = CKSURF(boundary_surfy_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
  if (cdim>2)
    diffusion->boundary_surf[2] = CKSURF(boundary_surfz_kernels, diff_order, cdim, vdim, poly_order, diffusion_id);
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  enum gkyl_diffusion_id diffusion_id, const bool *diff_in_dir, int diff_order,
  const struct gkyl_range *conf_range)
{
  struct dg_diffusion* diffusion = (struct dg_diffusion*) gkyl_malloc(sizeof(struct dg_diffusion));

  int cdim = cbasis->ndim;
  int vdim = basis->ndim - cdim;
  int poly_order = cbasis->poly_order;

  diffusion->num_equations = numeq_from_diffid(diffusion_id);
  diffusion->const_coeff = constcoeff_from_diffid(diffusion_id);
  diffusion->num_basis = basis->num_basis;
  for (int d=0; d<cdim; d++) diffusion->diff_in_dir[d] = diff_in_dir[d];

  int dirs_linidx = diffdirs_linidx(diff_in_dir, cdim);

  diffusion->conf_range = *conf_range;

  diffusion->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion->eqn.flags);
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_free);

  diffusion->D = D;

  // copy the host struct to device struct
  struct dg_diffusion* diffusion_cu = (struct dg_diffusion*) gkyl_cu_malloc(sizeof(struct dg_diffusion));
  gkyl_cu_memcpy(diffusion_cu, diffusion, sizeof(struct dg_diffusion), GKYL_CU_MEMCPY_H2D);
  dg_diffusion_set_cu_dev_ptrs<<<1,1>>>(diffusion_cu, cbasis->b_type, cdim, vdim, poly_order, diffusion_id, diff_order, dirs_linidx);

  // set parent on_dev pointer
  diffusion->eqn.on_dev = &diffusion_cu->eqn;

  return &diffusion->eqn;
}
