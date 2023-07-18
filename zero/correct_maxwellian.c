#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>

// Allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Determine the Jacobian for the maxwellian projection
void eval_jacob_tot(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

struct gkyl_correct_maxwellian 
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;
  
  struct gkyl_array *m0, *m1, *m2;
  struct gkyl_array *dm0, *dm1, *dm2;
  struct gkyl_array *ddm0, *ddm1, *ddm2;
  struct gkyl_array *m0scl;
  struct gkyl_array *fout;
  
  gkyl_proj_maxwellian_on_basis *maxwellian; // maxwellian projector
  gkyl_mom_calc *calc_m0, *calc_m1, *calc_m2; // moment calculators
  gkyl_dg_bin_op_mem *weak_divide, *weak_multiply, *weak_multiply_confPhase; // memory for weak operators
};

gkyl_correct_maxwellian *
gkyl_correct_maxwellian_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_local, const struct gkyl_range *phase_local double mass bool use_gpu)
{
  gkyl_correct_maxwellian *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;

  // Allocate memory
  up->m0scl = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->m0 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->m1 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->m2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->dm0 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->dm1 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->dm2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->ddm0 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->ddm1 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->ddm2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&conf_basis, &phase_basis, &conf_local, mass, "M0", use_gpu);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&conf_basis, &phase_basis, &conf_local, mass, "M1", use_gpu);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&conf_basis, &phase_basis, &conf_local, mass, "M2", use_gpu);

  gkyl_mom_calc up->calc_m0 = gkyl_mom_calc_new(&grid, M0_t, use_gpu);
  gkyl_mom_calc up->calc_m1 = gkyl_mom_calc_new(&grid, M1_t, use_gpu);
  gkyl_mom_calc up->calc_m1 = gkyl_mom_calc_new(&grid, M2_t, use_gpu); 
 
  up->weak_divide = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);
  up->weak_multiply = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);
  up->weak_multiply_confPhase = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);

  return up;
}

void gkyl_correct_maxwellian_gyrokinetic(gkyl_correct_maxwellian *cmax, struct gkyl_array *fin, const struct gkyl_array *m0_corr, const struct gkyl_array *m1_corr, const struct gkyl_array *m2_corr, const struct gkyl_array *bmag, double mass, double err_max, int iter_max, const struct gkyl_range *conf_local, const struct gkyl_range *phase_local, const struct gkyl_range *conf_local_ext, int poly_order, bool use_gpu)
  
{
  int vdim = cmax->phase_basis.ndim - cmax->conf_basis.ndim;
  int vdim_phys = vdim==1?1:3;
  
  struct gkyl_array *jacob_tot;
  jacob_tot = mkarr(cmax->conf_basis.num_basis, conf_local_ext.volume);

  struct gkyl_array *err = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);

  // Copy the input field into the output field
  gkyl_array_clear(cmj->fout, 0.0);
  gkyl_array_accumulate(cmj->fout, 1.0, fin);

  // Make the maxwellian projection routine
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&cmax->grid, &cmax->conf_basis, &cmax->phase_basis, poly_order+1, use_gpu); 

  // Rescale the maxwellian
  gkyl_mom_calc_advance(cmax->calc_m0, fout, cmax->m0);
  gkyl_dg_div_op(cmax->weak_divide, cmax->conf_basis, 0, cmax->m0scl, 0, m0_corr, 0, cmax->m0);
  gkyl_dg_mul_conf_phase_op_range(cmax->conf_basis, cmax->phase_basis, cmax->fout, cmax->m0scl, fout, conf_range, phase_range); // ranges to be decided

  // Calculate the moments
  gkyl_mom_calc_advance(cmax->calc_m0, &phase_local, &conf_local, fout, cmax->m0);
  gkyl_mom_calc_advance(cmax->calc_m1, &phase_local, &conf_local, fout, cmax->m1);
  gkyl_mom_calc_advance(cmax->calc_m2, &phase_local, &conf_local, fout, cmax->m2);
  
  // Main iteration loop
  while ((i<iter_max) && (err[0]>err_max) && (err[1]>err_max))
  {
    gkyl_array_clear(l2, 0.0);
    // Correct the moments
    gkyl_array_clear(cmax->ddm1, 0.0);
    gkyl_array_accumulate(cmax->ddm1, -1.0, cmax->m1);
    gkyl_array_accumulate(cmax->ddm1, 1.0, m1_corr);
    gkyl_array_accumulate(cmax->dm1, 1.0, cmax->ddm1);
    gkyl_array_clear(cmax->m1, 0.0);
    gkyl_array_accumulate(cmax->m1, 1.0, m1_corr);
    gkyl_array_accumulate(cmax->m1, 1.0, cmax->dm1);

    gkyl_dg_calc_l2_range(conf_basis, 0, err, 0, cmax->ddm1, conf_range);

    gkyl_array_clear(cmax->ddm2, 0.0);
    gkyl_array_accumulate(cmax->ddm2, -1.0, cmax->m2);
    gkyl_array_accumulate(cmax->ddm2, 1.0, m2_corr);
    gkyl_array_accumulate(cmax->dm2, 1.0, cmax->ddm2);
    gkyl_array_clear(cmax->m2, 0.0);
    gkyl_array_accumulate(cmax->m2, 1.0, m2_corr);
    gkyl_array_accumulate(cmax->m2, 1.0, cmax->dm2);

    gkyl_dg_calc_l2_range(conf_basis, 1, err, 0, cmax_>ddm2, conf_range);

    // Project the maxwellian
    // (1) proj_maxwellian expects the moments as a single array
    struct gkyl_array *moms_ho = mkarr((vdim+2)*cmax->conf_basis.num_basis, conf_local_ext.volume);
    gkyl_array_set_offset(moms_ho, 1., m0, 0*cmax->conf_basis.num_basis);
    gkyl_array_set_offset(moms_ho, 1., m1, 1*cmax->conf_basis.num_basis);
    gkyl_array_set_offset(moms_ho, 1., m2, (vdim+1)*cmax->conf_basis.num_basis);
    struct gkyl_array *moms;
    if (use_gpu)
    { // (2) copy host array to device
      moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim+2)*confBasis.num_basis, confLocal_ext.volume);
      gkyl_array_copy(moms, moms_ho);
    } 
    else 
    {
     moms = moms_ho;
    }
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &phase_local, &conf_local, moms, bmag, jacob_tot, mass, fout);
    // Rescale the maxwellian
    gkyl_mom_calc_advance(cmax->calc_m0, fout, cmax->m0);
    gkyl_dg_div_op(cmax->weak_divide, cmax->conf_basis, 0, cmax->m0scl, 0, m0_corr, 0, cmax->m0);
    gkyl_dg_mul_conf_phase_op_range(cmax->conf_basis, cmax->phase_basis, cmax->fout, cmax->m0scl, fout, conf_range, phase_range); // ranges to be decided
    // Increment i
    i += 1;
  } // Main iteration loop ends
}

void
gkyl_correct_maxwellian_release(gkyl_correct_maxwellian* cmax)
{
  gkyl_mom_calc_release(cmax->m0calc);
  gkyl_array_release(cmax->num_ratio);
  gkyl_dg_bin_op_mem_release(cmax->mem);
  
  gkyl_free(cmax);
}
