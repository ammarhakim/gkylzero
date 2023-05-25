/* -*- c++ -*- */

extern "C" {
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
}

// CUDA kernel to set device pointers to kernel that computes react rate.
__global__ static void
gkyl_iz_set_cu_ker_ptrs(const struct gkyl_basis *basis,
  struct gkyl_bc_sheath_gyrokinetic_kernels *kers)
{
  int cdim = basis->ndim;
  int poly_order = basis->poly_order;
  kers->react_ratef = ser_iz_react_rate_list[cdim-1].kernels[poly_order]; 
};

void
gkyl_iz_choose_react_ratef_kernel_cu(const struct gkyl_basis *basis,
  struct gkyl_iz_kernels *kers)
{
  gkyl_iz_set_cu_ker_ptrs<<<1,1>>>(basis, kers);
}

__global__ static void
gkyl_iz_react_rate_cu_ker(struct gkyl_range conf_rng, double elem_charge, double mass_elc, double E,
  double A, double K, double P, double X, const struct gkyl_array *moms_neut,
  struct gkyl_iz_kernels *kers, struct gkyl_array *coef_iz) {

  int cidx[3];
  
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < conf_rng.volume; linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&conf_rng, linc, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);
   
    const double *moms_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *m0_neut_d = &moms_neut_d[0]; 

    const double *vth_sq_neut_d = (const double*) gkyl_array_cfetch(vth_sq_neut, loc);
    const double *vth_sq_elc_d = (const double*) gkyl_array_cfetch(vth_sq_elc, loc);
    double *coef_iz_d = (double*) gkyl_array_fetch(coef_iz, loc); 

    // put kernels in separate object
    kers->react_ratef(elem_charge, mass_elc, E, A, K, P, X,
      m0_neut_d, vth_sq_neut_d, vth_sq_elc_d, coef_iz_d);
  }
}

void
gkyl_iz_react_rate_cu(const struct gkyl_dg_iz *up, const struct gkyl_array *moms_neut,
  const struct gkyl_array *coef_iz){
  int nblocks = up->conf_rng.nblocks, nthreads = up->conf_rng.nthreads;
  
  gkyl_iz_react_rate_cu_ker<<<nblocks, nthreads>>>gkyl_iz_react_rate_cu_ker(up->conf_rng, up->elem_charge,
    up->mass_elc, up->E, up->A, up->K, up->P, up->X, moms_neut->on_dev, up->vth_sq_neut->on_dev, up->vth_sq_ion->on_dev,
    up->kernels_cu, coef_iz->on_dev);
}					
