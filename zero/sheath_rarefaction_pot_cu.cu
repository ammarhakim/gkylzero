/* -*- c++ -*- */

extern "C" {
#include <gkyl_sheath_rarefaction_pot.h>
#include <gkyl_sheath_rarefaction_pot_priv.h>
}

// CUDA kernel to set device pointers to kernel that computes the reflected f.
__global__ static void
gkyl_sheath_rare_pot_set_cu_ker_ptrs(int dim, enum gkyl_basis_type b_type, int poly_order,
  enum gkyl_edge_loc edge, struct gkyl_sheath_rarefaction_pot_kernels *kers)
{

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->phimod = ser_sheath_rarepot_list[edge].list[dim-1].kernels[poly_order-1];
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      kers->phimod = tensor_sheath_rarepot_list[edge].list[dim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gkyl_sheath_rarepot_choose_phimod_kernel_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_sheath_rarefaction_pot_kernels *kers)
{
  gkyl_sheath_rare_pot_set_cu_ker_ptrs<<<1,1>>>(basis->ndim, basis->b_type, basis->poly_order, edge, kers);
}

__global__ static void
gkyl_sheath_rarefaction_pot_advance_cu_ker(struct gkyl_range skin_r,
  struct gkyl_sheath_rarefaction_pot_kernels *kers,
  double elem_charge,
  double mass_e, const struct gkyl_array *moms_e, const struct gkyl_array *m2par_e,
  double mass_i, const struct gkyl_array *moms_i, const struct gkyl_array *m2par_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi)
{
  int cidx[3];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < skin_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange

    gkyl_sub_range_inv_idx(&skin_r, linc, cidx);

    long loc = gkyl_range_idx(&skin_r, cidx);

    const double *momse_p = (const double*) gkyl_array_cfetch(moms_e, loc);
    const double *m2pare_p = (const double*) gkyl_array_cfetch(m2par_e, loc);
    const double *momsi_p = (const double*) gkyl_array_cfetch(moms_i, loc);
    const double *m2pari_p = (const double*) gkyl_array_cfetch(m2par_i, loc);
    const double *phiwall_p = (const double*) gkyl_array_cfetch(phi_wall, loc);
    double *phi_p = (double*) gkyl_array_fetch(phi, loc);

    // Modify the potential at the boundary to account for rarefaction wave.
    kers->phimod(elem_charge, mass_e, momse_p, m2pare_p, mass_i, momsi_p, m2pari_p, phiwall_p, phi_p);
  }
}

void
gkyl_sheath_rarefaction_pot_advance_cu(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_array *moms_e, const struct gkyl_array *m2par_e,
  const struct gkyl_array *moms_i, const struct gkyl_array *m2par_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi)
{
  int nblocks = up->skin_r.nblocks, nthreads = up->skin_r.nthreads;

  gkyl_sheath_rarefaction_pot_advance_cu_ker<<<nblocks, nthreads>>>(up->skin_r, 
    up->kernels_cu, up->elem_charge, up->mass_e, moms_e->on_dev, m2par_e->on_dev,
    up->mass_i, moms_i->on_dev, m2par_i->on_dev, phi_wall->on_dev, phi->on_dev);
}
