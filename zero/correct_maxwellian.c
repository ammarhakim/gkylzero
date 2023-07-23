#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>
#include <math.h>

// Allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct gkyl_correct_maxwellian 
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;
  
  struct gkyl_array *m0, *m1, *m2;
  struct gkyl_array *dm1, *dm2;
  struct gkyl_array *ddm1, *ddm2;
  struct gkyl_array *m0scl;
  //struct gkyl_array *fout;
  
  gkyl_proj_maxwellian_on_basis *maxwellian; // maxwellian projector
  gkyl_mom_calc *m0calc, *m1calc, *m2calc; // moment calculators
  gkyl_dg_bin_op_mem *weak_divide, *weak_multiply, *weak_multiply_confPhase; // memory for weak operators
};

gkyl_correct_maxwellian *
gkyl_correct_maxwellian_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext, const struct gkyl_array *bmag, double mass, bool use_gpu)
{
  gkyl_correct_maxwellian *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;

  // Allocate memory
  up->m0scl = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->m0 = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->m1 = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->m2 = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->dm1 = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->dm2 = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->ddm1 = mkarr(conf_basis->num_basis, conf_local_ext->volume);
  up->ddm2 = mkarr(conf_basis->num_basis, conf_local_ext->volume);

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(conf_basis, phase_basis, conf_local, mass, "M0", use_gpu);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(conf_basis, phase_basis, conf_local, mass, "M1", use_gpu);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(conf_basis, phase_basis, conf_local, mass, "M2", use_gpu);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  up->m0calc = gkyl_mom_calc_new(grid, M0_t, use_gpu);
  up->m1calc = gkyl_mom_calc_new(grid, M1_t, use_gpu);
  up->m2calc = gkyl_mom_calc_new(grid, M2_t, use_gpu); 
 
  up->weak_divide = gkyl_dg_bin_op_mem_new(conf_local->volume, conf_basis->num_basis);
  up->weak_multiply = gkyl_dg_bin_op_mem_new(conf_local->volume, conf_basis->num_basis);
  up->weak_multiply_confPhase = gkyl_dg_bin_op_mem_new(conf_local->volume, conf_basis->num_basis);

  return up;
}

void gkyl_correct_maxwellian_gyrokinetic(gkyl_correct_maxwellian *cmax, struct gkyl_array *fM, const struct gkyl_array *m0_corr, const struct gkyl_array *m1_corr, const struct gkyl_array *m2_corr, const struct gkyl_array *jacob_tot, const struct gkyl_array *bmag, double mass, double err_max, int iter_max, const struct gkyl_range *conf_local, const struct gkyl_range *phase_local, const struct gkyl_range *conf_local_ext, int poly_order, bool use_gpu)
  
{
  int vdim = cmax->phase_basis.ndim - cmax->conf_basis.ndim;
  int vdim_phys = vdim==1?1:3;
  
  double err1[1];
  double err2[1];
  struct gkyl_array *mvals1 = mkarr(1, conf_local->volume);
  struct gkyl_array *mvals2 = mkarr(1, conf_local->volume);
  gkyl_array_clear(mvals1, 0.0);
  gkyl_array_clear(mvals2, 0.0);

  // Copy the input field into the output field
  //gkyl_array_clear(cmax->fout, 0.0);
  //gkyl_array_accumulate(cmax->fout, 1.0, fin);

  // Make the maxwellian projection routine
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&cmax->grid, &cmax->conf_basis, &cmax->phase_basis, poly_order+1, use_gpu); 

  // Rescale the maxwellian
  gkyl_mom_calc_advance(cmax->m0calc, phase_local, conf_local, fM, cmax->m0);
  gkyl_dg_div_op_range(cmax->weak_divide, cmax->conf_basis, 0, cmax->m0scl, 0, m0_corr, 0, cmax->m0, conf_local);
  gkyl_dg_mul_conf_phase_op_range(&cmax->conf_basis, &cmax->phase_basis, fM, cmax->m0scl, fM, conf_local, phase_local);

  // Calculate the moments
  gkyl_mom_calc_advance(cmax->m0calc, phase_local, conf_local, fM, cmax->m0);
  gkyl_mom_calc_advance(cmax->m1calc, phase_local, conf_local, fM, cmax->m1);
  gkyl_mom_calc_advance(cmax->m2calc, phase_local, conf_local, fM, cmax->m2);

  // Initialize the error
  gkyl_array_clear(cmax->ddm1, 0.0);
  gkyl_array_accumulate(cmax->ddm1, -1.0, cmax->m1);
  gkyl_array_accumulate(cmax->ddm1, 1.0, m1_corr);
  gkyl_dg_calc_l2_range(cmax->conf_basis, 0, mvals1, 0, cmax->ddm1, *conf_local);
  gkyl_array_reduce_range(err1, mvals1, GKYL_SUM, *conf_local);

  gkyl_array_clear(cmax->ddm2, 0.0);
  gkyl_array_accumulate(cmax->ddm2, -1.0, cmax->m2);
  gkyl_array_accumulate(cmax->ddm2, 1.0, m2_corr);
  gkyl_dg_calc_l2_range(cmax->conf_basis, 0, mvals2, 0, cmax->ddm2, *conf_local);
  gkyl_array_reduce_range(err2, mvals2, GKYL_SUM, *conf_local);

  // Main iteration loop
  int i = 0;
  while ((i<iter_max) && (err1[0]>err_max) || (err2[0]>err_max))
  {
    printf("Iteration %d\n", i);
    // Correct the moments
    gkyl_array_clear(cmax->ddm1, 0.0);
    gkyl_array_accumulate(cmax->ddm1, -1.0, cmax->m1);
    gkyl_array_accumulate(cmax->ddm1, 1.0, m1_corr);
    gkyl_array_accumulate(cmax->dm1, 1.0, cmax->ddm1);
    gkyl_array_clear(cmax->m1, 0.0);
    gkyl_array_accumulate(cmax->m1, 1.0, m1_corr);
    gkyl_array_accumulate(cmax->m1, 1.0, cmax->dm1);

    gkyl_dg_calc_l2_range(cmax->conf_basis, 0, mvals1, 0, cmax->ddm1, *conf_local);
    gkyl_array_scale_range(mvals1, cmax->grid.cellVolume, *conf_local);
    gkyl_array_reduce_range(err1, mvals1, GKYL_SUM, *conf_local);
    err1[0] = sqrt(err1[0]/cmax->grid.cellVolume/conf_local->volume); 

    gkyl_array_clear(cmax->ddm2, 0.0);
    gkyl_array_accumulate(cmax->ddm2, -1.0, cmax->m2);
    gkyl_array_accumulate(cmax->ddm2, 1.0, m2_corr);
    gkyl_array_accumulate(cmax->dm2, 1.0, cmax->ddm2);
    gkyl_array_clear(cmax->m2, 0.0);
    gkyl_array_accumulate(cmax->m2, 1.0, m2_corr);
    gkyl_array_accumulate(cmax->m2, 1.0, cmax->dm2);

    gkyl_dg_calc_l2_range(cmax->conf_basis, 0, mvals2, 0, cmax->ddm2, *conf_local);
    gkyl_array_scale_range(mvals2, cmax->grid.cellVolume, *conf_local);
    gkyl_array_reduce_range(err2, mvals2, GKYL_SUM, *conf_local);
    err2[0] = sqrt(err2[0]/cmax->grid.cellVolume/conf_local->volume); 

    printf("err1=%10.8e, err2=%10.8e\n", err1[0], err2[0]);

    // Project the maxwellian
    // (1) proj_maxwellian expects the moments as a single array
    struct gkyl_array *moms_ho = mkarr(3*cmax->conf_basis.num_basis, conf_local_ext->volume);
    gkyl_array_set_offset(moms_ho, 1., cmax->m0, 0*cmax->conf_basis.num_basis);
    gkyl_array_set_offset(moms_ho, 1., cmax->m1, 1*cmax->conf_basis.num_basis);
    gkyl_array_set_offset(moms_ho, 1., cmax->m2, 2*cmax->conf_basis.num_basis);
    struct gkyl_array *moms;
    if (use_gpu)
    { // (2) copy host array to device
      moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*cmax->conf_basis.num_basis, conf_local_ext->volume);
      gkyl_array_copy(moms, moms_ho);
    } 
    else 
    {
     moms = moms_ho;
    }
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, phase_local, conf_local, moms, bmag, jacob_tot, mass, fM);
    // Rescale the maxwellian
    gkyl_mom_calc_advance(cmax->m0calc, phase_local, conf_local, fM, cmax->m0);
    gkyl_dg_div_op_range(cmax->weak_divide, cmax->conf_basis, 0, cmax->m0scl, 0, m0_corr, 0, cmax->m0, conf_local);
    gkyl_dg_mul_conf_phase_op_range(&cmax->conf_basis, &cmax->phase_basis, fM, cmax->m0scl, fM, conf_local, phase_local);
    // Calculate the moments
    gkyl_mom_calc_advance(cmax->m0calc, phase_local, conf_local, fM, cmax->m0);
    gkyl_mom_calc_advance(cmax->m1calc, phase_local, conf_local, fM, cmax->m1);
    gkyl_mom_calc_advance(cmax->m2calc, phase_local, conf_local, fM, cmax->m2);
    // Increment i
    i += 1;
  } // Main iteration loop ends
  gkyl_array_release(mvals1);
  gkyl_array_release(mvals2);
}

void
gkyl_correct_maxwellian_release(gkyl_correct_maxwellian* cmax)
{
  gkyl_free(cmax);
}
