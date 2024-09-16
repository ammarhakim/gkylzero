// extern "C" {
//   #include <gkyl_range.h>
//   #include <gkyl_util.h>
//   #include <gkyl_alloc.h>
//   #include <gkyl_array_ops.h>
//   #include <gkyl_array_ops_priv.h>
//   #include <gkyl_dg_fpo_vlasov_diff_coeff.h>
//   #include <gkyl_dg_fpo_vlasov_diff_coeff_priv.h>
// } 
//
// __global__ void
// gkyl_calc_fpo_diff_coeff_recovery_cu_kernel(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery,
//   const struct gkyl_rect_grid grid, struct gkyl_basis pbasis, 
//   const struct gkyl_range phase_range, const struct gkyl_range conf_range, 
//   const struct gkyl_array *gamma, const struct gkyl_array *fpo_g,
//   const struct gkyl_array *fpo_g_surf, const struct gkyl_array *fpo_dgdv_surf, 
//   const struct gkyl_array *fpo_d2gdv2_surf, 
//   struct gkyl_array *fpo_diff_coeff, struct gkyl_array *fpo_diff_coeff_surf)
// {
//   int pdim = pbasis.ndim;
//   int vdim = 3;
//   int cdim = pdim-vdim; 
//
//   int poly_order = pbasis.poly_order;
//
//   // Indices in each direction
//   int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], conf_idxc[GKYL_MAX_DIM];
//   int idx_edge[GKYL_MAX_DIM], idx_skin[GKYL_MAX_DIM];
//   int edge;
//
//   for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
//       tid < phase_range.volume;
//       tid += gridDim.x*blockDim.x)
//   {
//     gkyl_sub_range_inv_idx(&phase_range, tid, idxc);  
//     long linp = gkyl_range_idx(&phase_range, idxc);
//     long linc = gkyl_range_idx(&conf_range, idxc);
//
//     const double *fpo_dgdv_surf_c = (const double *)gkyl_array_cfetch(fpo_dgdv_surf, linp);
//     const double *fpo_d2gdv2_surf_c = (const double *)gkyl_array_cfetch(fpo_d2gdv2_surf, linp);
//     double *fpo_diff_coeff_c = (doubl e*)gkyl_array_fetch(fpo_diff_coeff, linp);
//
//     const double *gamma_c = (const double *)gkyl_array_cfetch(gamma, linc);
//   }
// }
//
// void 
// gkyl_calc_fpo_diff_coeff_recovery_cu(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery,
//   const struct gkyl_rect_grid *grid, struct gkyl_basis pbasis,
//   const struct gkyl_range *phase_range, const struct gkyl_range *conf_range, 
//   const struct gkyl_array *gamma, const struct gkyl_array *fpo_g,
//   const struct gkyl_array *fpo_g_surf, const struct gkyl_array *fpo_dgdv_surf, 
//   const struct gkyl_array *fpo_d2gdv2_surf, 
//   struct gkyl_array *fpo_diff_coeff, struct gkyl_array *fpo_diff_coeff_surf)
// {
//   int nblocks = phase_range->nblocks;
//   int nthreads = phase_range->nthreads;
//
//   gkyl_calc_fpo_diff_coeff_recovery_cu_kernel<<<nblocks, nthreads>>>(coeff_recovery->on_dev, 
//     *grid, pbasis, *phase_range, *conf_range, gamma->on_dev, 
//     fpo_g->on_dev, fpo_g_surf->on_dev, fpo_dgdv_surf->on_dev,
//     fpo_d2gdv2_surf->on_dev, fpo_diff_coeff->on_dev, fpo_diff_coeff_surf->on_dev);
// }
