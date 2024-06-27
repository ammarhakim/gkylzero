#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_fpo_vlasov_coeffs_correct.h>
#include <gkyl_fpo_vlasov_coeffs_correct_priv.h>

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct gkyl_fpo_coeffs_correct*
gkyl_fpo_coeffs_correct_new(bool use_gpu, const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_rng)
{
  struct gkyl_fpo_coeffs_correct *up = gkyl_malloc(sizeof(struct gkyl_fpo_coeffs_correct)); 

  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  int vdim = 3;
  
  up->grid = grid;
  up->conf_basis = conf_basis;

  up->is_first = true;

  // Matrices and memory for linear solve
  up->As = 0;
  up->xs = 0;
  up->mem = 0;

  // Kernels for setting linear system matrices
  const gkyl_fpo_coeffs_correct_mat_set_kern_list *mat_set_kern_list;
  const gkyl_fpo_coeffs_correct_accum_kern_list *accum_kern_list;

  switch (conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mat_set_kern_list = ser_fpo_coeffs_correct_mat_set_kernels;
      accum_kern_list = ser_fpo_coeffs_correct_accum_kernels;
      break;

    default:
      assert(false);
      break;
  }
  up->mat_set_kernel = mat_set_kern_list[cdim].kernels[poly_order];
  up->accum_kernel = accum_kern_list[cdim].kernels[poly_order];

  return up;
}

void gkyl_fpo_coeffs_correct_advance(struct gkyl_fpo_coeffs_correct *up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *fpo_moms, const struct gkyl_array *boundary_corrections,
  const struct gkyl_array *moms, struct gkyl_array *drag_diff_coeff_corrs,
  struct gkyl_array *drag_coeff, struct gkyl_array *drag_coeff_surf,
  struct gkyl_array *diff_coeff, struct gkyl_array *diff_coeff_surf)
{
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  // allocate memory for use in kernels
  int nc = up->conf_basis->num_basis;
  int vdim = 3;
  int N = nc*(vdim + 1);

  // Initialize matrices if this is the first call to advance
  if (up->is_first) {
    up->As = gkyl_nmat_new(conf_rng->volume, N, N);
    up->xs = gkyl_nmat_new(conf_rng->volume, N, 1);
    up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);
    up->is_first = false;
  }

  // Loop over configuration space cells
  gkyl_range_iter_init(&conf_iter, conf_rng);
  long count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long linc = gkyl_range_idx(conf_rng, conf_iter.idx);

    struct gkyl_mat lhs = gkyl_nmat_get(up->As, count);
    struct gkyl_mat rhs = gkyl_nmat_get(up->xs, count);
    gkyl_mat_clear(&lhs, 0.0);
    gkyl_mat_clear(&rhs, 0.0);

    const double *fpo_moms_d = gkyl_array_cfetch(fpo_moms, linc);
    const double *boundary_corrections_d = gkyl_array_cfetch(boundary_corrections, linc);
    const double *moms_d = gkyl_array_cfetch(moms, linc); 

    // Set matrix elements
    up->mat_set_kernel(&lhs, &rhs, fpo_moms_d, boundary_corrections_d, moms_d);
    count += 1;
  }

  // Solve linear system to determine corrections
  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  // Loop over configuration space cells, retrieve solutions,
  // and accumulate corrections onto drag/diffusion coefficients
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<conf_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, conf_rng);
  count = 0;
  while (gkyl_range_iter_next(&conf_iter)) {
    long linc = gkyl_range_idx(conf_rng, conf_iter.idx);

    const struct gkyl_mat out = gkyl_nmat_get(up->xs, count);
    double *drag_diff_coeff_corrs_d = gkyl_array_fetch(drag_diff_coeff_corrs, linc);
    for (size_t i=0; i<N; ++i) {
      drag_diff_coeff_corrs_d[i] = gkyl_mat_get(&out, i, 0);
    }
    count += 1;

    // Inner loop over velocity space to accumulate corrections onto drag and diff coeffs
    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      long linp = gkyl_range_idx(&vel_rng, vel_iter.idx);
      
      const double *drag_diff_coeff_corrs_d = gkyl_array_cfetch(drag_diff_coeff_corrs, linc);
      double *drag_coeff_d = gkyl_array_fetch(drag_coeff, linp);
      double *drag_coeff_surf_d = gkyl_array_fetch(drag_coeff_surf, linp);
      double *diff_coeff_d = gkyl_array_fetch(diff_coeff, linp);
      double *diff_coeff_surf_d = gkyl_array_fetch(diff_coeff_surf, linp);

      // Call to kernels to accumulate corrections
      up->accum_kernel(drag_diff_coeff_corrs_d, drag_coeff_d, drag_coeff_surf_d, 
        diff_coeff_d, diff_coeff_surf_d); 
    }
  }
}

