#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_correct_maxwellian_gyrokinetic.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>
#include <math.h>

// Allocate cu_dev array
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) : gkyl_array_new(GKYL_DOUBLE, nc, size);  
  return a;
}

struct gkyl_correct_maxwellian_gyrokinetic 
{
  bool use_gpu;
  double mass;
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  struct gkyl_array *bmag, *jacob_tot;  
  struct gkyl_array *m0_in, *m0, *m0scl;
  struct gkyl_array *m12_in, *m12, *dm12, *ddm12;
  struct gkyl_array *moms;
  
  gkyl_proj_maxwellian_on_basis *proj_maxwellian; // maxwellian projector
  gkyl_mom_calc *m0calc, *momsCalc; // moment calculators
  gkyl_dg_bin_op_mem *weak_divide, *weak_multiply, *weak_multiply_confPhase; // memory for weak operators

  double *err1_cu, *err2_cu;
  struct gkyl_array *mvals1, *mvals2;
};

gkyl_correct_maxwellian_gyrokinetic *
gkyl_correct_maxwellian_gyrokinetic_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_local, const struct gkyl_range *conf_local_ext, const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass, bool use_gpu)
{
  gkyl_correct_maxwellian_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->mass = mass;
  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;

  up->bmag = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->jacob_tot = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  if (use_gpu) {
    gkyl_array_copy(up->bmag, bmag);
    gkyl_array_copy(up->jacob_tot, jacob_tot);
  } else {
    up->bmag = bmag;
    up->jacob_tot = jacob_tot;
  }

  // Allocate memory
  up->m0_in = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->m0 = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->m0scl = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->m12_in = mkarr(2*conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->m12 = mkarr(2*conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->dm12 = mkarr(2*conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->ddm12 = mkarr(2*conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->moms = mkarr(3*conf_basis->num_basis, conf_local_ext->volume, use_gpu);

  up->proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(grid, conf_basis, phase_basis, conf_basis->num_basis+1, use_gpu); 

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(conf_basis, phase_basis, conf_local, mass, "M0", use_gpu);
  struct gkyl_mom_type *MOMS_t = gkyl_mom_gyrokinetic_new(conf_basis, phase_basis, conf_local, mass, "ThreeMoments", use_gpu);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(MOMS_t, bmag);
  up->m0calc = gkyl_mom_calc_new(grid, M0_t, use_gpu);
  up->momsCalc = gkyl_mom_calc_new(grid, MOMS_t, use_gpu);

  if (use_gpu) {
    up->weak_divide = gkyl_dg_bin_op_mem_cu_dev_new(conf_local->volume, conf_basis->num_basis);
    up->weak_multiply = gkyl_dg_bin_op_mem_cu_dev_new(conf_local->volume, conf_basis->num_basis);
    up->weak_multiply_confPhase = gkyl_dg_bin_op_mem_cu_dev_new(conf_local->volume, conf_basis->num_basis);
    up->err1_cu = gkyl_cu_malloc(sizeof(double[1]));
    up->err2_cu = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    up->weak_divide = gkyl_dg_bin_op_mem_new(conf_local->volume, conf_basis->num_basis);
    up->weak_multiply = gkyl_dg_bin_op_mem_new(conf_local->volume, conf_basis->num_basis);
    up->weak_multiply_confPhase = gkyl_dg_bin_op_mem_new(conf_local->volume, conf_basis->num_basis);
  }

  up->mvals1 = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);
  up->mvals2 = mkarr(conf_basis->num_basis, conf_local_ext->volume, use_gpu);

  return up;
}

void gkyl_correct_maxwellian_gyrokinetic_fix(gkyl_correct_maxwellian_gyrokinetic *up, struct gkyl_array *fM, const struct gkyl_array *moms_in, double err_max, int iter_max, const struct gkyl_range *conf_local, const struct gkyl_range *phase_local, const struct gkyl_range *conf_local_ext)
{
  double err1[1], err2[1];
  
  // Decompose the input moments
  gkyl_array_set_offset(up->m0_in, 1., moms_in, 0*up->conf_basis.num_basis);
  gkyl_array_set_offset(up->m12_in, 1., moms_in, 1*up->conf_basis.num_basis);

  // Rescale the maxwellian
  (up->use_gpu) ? gkyl_mom_calc_advance_cu(up->m0calc, phase_local, conf_local, fM, up->m0) : gkyl_mom_calc_advance(up->m0calc, phase_local, conf_local, fM, up->m0);
  gkyl_dg_div_op_range(up->weak_divide, up->conf_basis, 0, up->m0scl, 0, up->m0_in, 0, up->m0, conf_local);
  gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis, fM, up->m0scl, fM, conf_local, phase_local);

  // Calculate the moments
  (up->use_gpu) ? gkyl_mom_calc_advance_cu(up->momsCalc, phase_local, conf_local, fM, up->moms) : gkyl_mom_calc_advance(up->momsCalc, phase_local, conf_local, fM, up->moms);
  gkyl_array_set_offset(up->m0, 1., up->moms, 0*up->conf_basis.num_basis);
  gkyl_array_set_offset(up->m12, 1., up->moms, 1*up->conf_basis.num_basis);

  // Initialize the error
  gkyl_array_clear(up->ddm12, 0.0);
  gkyl_array_accumulate(up->ddm12, -1.0, up->m12);
  gkyl_array_accumulate(up->ddm12, 1.0, up->m12_in);
  
  gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals1, 0, up->ddm12, *conf_local);
  gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals2, 1, up->ddm12, *conf_local);
  gkyl_array_scale_range(up->mvals1, up->grid.cellVolume, *conf_local);
  gkyl_array_scale_range(up->mvals2, up->grid.cellVolume, *conf_local);
  if (up->use_gpu) {
    gkyl_array_reduce_range(up->err1_cu, up->mvals1, GKYL_SUM, *conf_local);
    gkyl_array_reduce_range(up->err2_cu, up->mvals2, GKYL_SUM, *conf_local);
    gkyl_cu_memcpy(err1, up->err1_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(err2, up->err2_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(err1, up->mvals1, GKYL_SUM, *conf_local);
    gkyl_array_reduce_range(err2, up->mvals2, GKYL_SUM, *conf_local);
  }

  // Main iteration loop
  int i = 0;
  while ((i<iter_max) && (err1[0]>err_max) || (err2[0]>err_max))
  {
    // Correct the moments
    gkyl_array_clear(up->ddm12, 0.0);
    gkyl_array_accumulate(up->ddm12, -1.0, up->m12);
    gkyl_array_accumulate(up->ddm12, 1.0, up->m12_in);
    gkyl_array_accumulate(up->dm12, 1.0, up->ddm12);
    gkyl_array_clear(up->m12, 0.0);
    gkyl_array_accumulate(up->m12, 1.0, up->m12_in);
    gkyl_array_accumulate(up->m12, 1.0, up->dm12);

    gkyl_array_clear(up->mvals1, 0.0);
    gkyl_array_clear(up->mvals2, 0.0);
    gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals1, 0, up->ddm12, *conf_local);
    gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals2, 1, up->ddm12, *conf_local);
    gkyl_array_scale_range(up->mvals1, up->grid.cellVolume, *conf_local);
    gkyl_array_scale_range(up->mvals2, up->grid.cellVolume, *conf_local);
    if (up->use_gpu) {
      gkyl_array_reduce_range(up->err1_cu, up->mvals1, GKYL_SUM, *conf_local);
      gkyl_array_reduce_range(up->err2_cu, up->mvals2, GKYL_SUM, *conf_local);
      gkyl_cu_memcpy(err1, up->err1_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
      gkyl_cu_memcpy(err2, up->err2_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      gkyl_array_reduce_range(err1, up->mvals1, GKYL_SUM, *conf_local);
      gkyl_array_reduce_range(err2, up->mvals2, GKYL_SUM, *conf_local);
    }
    err1[0] = sqrt(err1[0]/up->grid.cellVolume/conf_local->volume); 
    err2[0] = sqrt(err2[0]/up->grid.cellVolume/conf_local->volume); 

    // Project the maxwellian
    gkyl_array_set_offset(up->moms, 1., up->m0, 0*up->conf_basis.num_basis);
    gkyl_array_set_offset(up->moms, 1., up->m12, 1*up->conf_basis.num_basis);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(up->proj_maxwellian, phase_local, conf_local, up->moms, up->bmag, up->jacob_tot, up->mass, fM);
    // Rescale the maxwellian
    (up->use_gpu) ? gkyl_mom_calc_advance_cu(up->m0calc, phase_local, conf_local, fM, up->m0) : gkyl_mom_calc_advance(up->m0calc, phase_local, conf_local, fM, up->m0);
    gkyl_dg_div_op_range(up->weak_divide, up->conf_basis, 0, up->m0scl, 0, up->m0_in, 0, up->m0, conf_local);
    gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis, fM, up->m0scl, fM, conf_local, phase_local);
    // Calculate the moments
    (up->use_gpu) ? gkyl_mom_calc_advance_cu(up->momsCalc, phase_local, conf_local, fM, up->moms) : gkyl_mom_calc_advance(up->momsCalc, phase_local, conf_local, fM, up->moms);
    gkyl_array_set_offset(up->m0, 1., up->moms, 0*up->conf_basis.num_basis);
    gkyl_array_set_offset(up->m12, 1., up->moms, 1*up->conf_basis.num_basis);
    
    i += 1;
  } // Main iteration loop ends
}

void
gkyl_correct_maxwellian_gyrokinetic_release(gkyl_correct_maxwellian_gyrokinetic* up)
{
  if (up->use_gpu) {
    gkyl_cu_free(up->err1_cu);
    gkyl_cu_free(up->err2_cu);
  }
  gkyl_array_release(up->bmag);
  gkyl_array_release(up->jacob_tot);
  gkyl_array_release(up->m0_in);
  gkyl_array_release(up->m0);
  gkyl_array_release(up->m0scl);
  gkyl_array_release(up->m12_in);
  gkyl_array_release(up->m12);
  gkyl_array_release(up->dm12);
  gkyl_array_release(up->ddm12);
  gkyl_array_release(up->mvals1);
  gkyl_array_release(up->mvals2);
  gkyl_array_release(up->moms);
  gkyl_mom_calc_release(up->m0calc);
  gkyl_mom_calc_release(up->momsCalc);
  gkyl_proj_maxwellian_on_basis_release(up->proj_maxwellian);
  gkyl_dg_bin_op_mem_release(up->weak_divide);
  gkyl_dg_bin_op_mem_release(up->weak_multiply);
  gkyl_dg_bin_op_mem_release(up->weak_multiply_confPhase);
  gkyl_free(up);
}
