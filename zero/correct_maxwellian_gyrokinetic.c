#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_correct_maxwellian_gyrokinetic.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>
#include <math.h>
#include <stdio.h>

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
  long max_iter;
  double eps_err;
  struct gkyl_rect_grid phase_grid, conf_grid;
  struct gkyl_basis conf_basis, phase_basis;

  const struct gk_geometry *gk_geom; // Pointer to geometry struct
  struct gkyl_array *m0_tar, *m0_num, *m0_scl;
  struct gkyl_array *m12, *m12_tar, *m12_num, *delta_m12;
  struct gkyl_array *moms, *moms_num;
  
  struct gkyl_proj_maxwellian_on_basis *proj_maxwellian; // maxwellian projector
  struct gkyl_dg_updater_moment *m0_calc, *moms_calc; // moment calculators
  struct gkyl_dg_bin_op_mem *weak_divide; // memory for weak operators

  double *err1_cu, *err2_cu;
  struct gkyl_array *mvals, *mvals1, *mvals2;
  struct gkyl_array *mvals_host;

  struct gkyl_velocity_map *vel_map;
};

gkyl_correct_maxwellian_gyrokinetic*
gkyl_correct_maxwellian_gyrokinetic_new(const struct gkyl_correct_maxwellian_gyrokinetic_inp *inp)
{
  gkyl_correct_maxwellian_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = inp->use_gpu;
  up->max_iter = inp->max_iter;
  up->eps_err = inp->eps_err;
  up->mass = inp->mass;
  up->phase_grid = *inp->phase_grid;
  up->conf_grid = *inp->conf_grid;
  up->phase_basis = *inp->phase_basis;
  up->conf_basis = *inp->conf_basis;

  // Ccquire pointer to geometry and v-space jacobian.
  up->gk_geom = gkyl_gk_geometry_acquire(inp->gk_geom);
  up->vel_map = gkyl_velocity_map_acquire(inp->vel_map);

  // Allocate memory
  up->m0_tar = mkarr(up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->m0_num = mkarr(up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->m0_scl = mkarr(up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->m12 = mkarr(2*up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->m12_tar = mkarr(2*up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->m12_num = mkarr(2*up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->delta_m12 = mkarr(2*up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->moms = mkarr(3*up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);
  up->moms_num = mkarr(3*up->conf_basis.num_basis, inp->conf_local_ext->volume, up->use_gpu);

  struct gkyl_proj_maxwellian_on_basis_inp proj_max_inp = {
    .grid = &up->phase_grid,
    .conf_basis = &up->conf_basis,
    .phase_basis = &up->phase_basis,
    .num_quad = up->phase_basis.poly_order+1,
    .vel_map = inp->vel_map,
    .use_gpu = up->use_gpu,
  };
  up->proj_maxwellian = gkyl_proj_maxwellian_on_basis_inew(&proj_max_inp);

  up->m0_calc = gkyl_dg_updater_moment_gyrokinetic_new(&up->phase_grid, &up->conf_basis, 
    &up->phase_basis, inp->conf_local, up->mass, inp->vel_map, inp->gk_geom,
    "M0", false, up->use_gpu);

  up->moms_calc = gkyl_dg_updater_moment_gyrokinetic_new(&up->phase_grid, &up->conf_basis, 
    &up->phase_basis, inp->conf_local, up->mass, inp->vel_map, inp->gk_geom,
    "ThreeMoments", false, up->use_gpu); 

  if (up->use_gpu) {
    up->weak_divide = gkyl_dg_bin_op_mem_cu_dev_new(inp->conf_local->volume, up->conf_basis.num_basis);
    up->err1_cu = gkyl_cu_malloc(sizeof(double[1]));
    up->err2_cu = gkyl_cu_malloc(sizeof(double[1]));
  }
  else {
    up->weak_divide = gkyl_dg_bin_op_mem_new(inp->conf_local->volume, up->conf_basis.num_basis);
  }

  up->mvals = mkarr(3, inp->conf_local_ext->volume, up->use_gpu);
  up->mvals1 = mkarr(1, inp->conf_local_ext->volume, up->use_gpu);
  up->mvals2 = mkarr(1, inp->conf_local_ext->volume, up->use_gpu);

  // Compute the mean on the host-side
  up->mvals_host = up->mvals;
  if (up->use_gpu)
    up->mvals_host = mkarr(3, inp->conf_local_ext->volume, false);

  return up;
}

void gkyl_correct_maxwellian_gyrokinetic_advance(gkyl_correct_maxwellian_gyrokinetic *up, 
  struct gkyl_array *fM, const struct gkyl_array *moms_tar, 
  const struct gkyl_range *conf_local, const struct gkyl_range *phase_local)
{
  double epsilon = 1.00;
  double err1[1], err2[1], mean[3];
  // Decompose the input moments
  gkyl_array_set_offset(up->m0_tar, 1., moms_tar, 0*up->conf_basis.num_basis);
  gkyl_array_set_offset(up->m12_tar, 1., moms_tar, 1*up->conf_basis.num_basis);

  // Rescale the maxwellian
  gkyl_dg_updater_moment_gyrokinetic_advance(up->m0_calc, phase_local, conf_local, fM, up->m0_num);
  gkyl_dg_div_op_range(up->weak_divide, up->conf_basis, 0, up->m0_scl, 0, up->m0_tar, 0, up->m0_num, conf_local);
  gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis, fM, up->m0_scl, fM, conf_local, phase_local);

  // Calculate the moments
  gkyl_dg_updater_moment_gyrokinetic_advance(up->moms_calc, phase_local, conf_local, fM, up->moms_num);
  gkyl_array_set_offset(up->m0_num, 1., up->moms_num, 0*up->conf_basis.num_basis);
  gkyl_array_set_offset(up->m12_num, 1., up->moms_num, 1*up->conf_basis.num_basis);

  // Initialize the m12 and delta_m12
  gkyl_array_set(up->m12, 1.0, up->m12_tar);
  gkyl_array_clear(up->delta_m12, 0.0);
  gkyl_array_accumulate(up->delta_m12, 1.0, up->m12_tar);
  gkyl_array_accumulate(up->delta_m12, -1.0, up->m12_num);
  // Calculate the average
  gkyl_array_clear(up->mvals, 0.0);
  gkyl_dg_calc_average_range(up->conf_basis, 0, up->mvals, 0, moms_tar, *conf_local);
  gkyl_dg_calc_average_range(up->conf_basis, 1, up->mvals, 1, moms_tar, *conf_local);
  gkyl_dg_calc_average_range(up->conf_basis, 2, up->mvals, 2, moms_tar, *conf_local);
  // Compute the reduction on the host-side for ease of subsequent divisions 
  // This reduction is only done once
  if (up->use_gpu) 
    gkyl_array_copy(up->mvals_host, up->mvals);
  gkyl_array_reduce_range(mean, up->mvals_host, GKYL_SUM, conf_local);
  for (int j=0; j<3; j++) {
    mean[j] = mean[j] / conf_local->volume;
  }
  mean[1] = epsilon * sqrt((mean[2]-mean[1]*mean[1]/mean[0])/mean[0])*mean[0];

  // Calculate the absolute error
  gkyl_array_clear(up->mvals1, 0.0);
  gkyl_array_clear(up->mvals2, 0.0);
  gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals1, 0, up->delta_m12, *conf_local);
  gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals2, 1, up->delta_m12, *conf_local);
  gkyl_array_scale_range(up->mvals1, up->conf_grid.cellVolume, conf_local);
  gkyl_array_scale_range(up->mvals2, up->conf_grid.cellVolume, conf_local);
  if (up->use_gpu) {
    gkyl_array_reduce_range(up->err1_cu, up->mvals1, GKYL_SUM, conf_local);
    gkyl_array_reduce_range(up->err2_cu, up->mvals2, GKYL_SUM, conf_local);
    gkyl_cu_memcpy(err1, up->err1_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(err2, up->err2_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(err1, up->mvals1, GKYL_SUM, conf_local);
    gkyl_array_reduce_range(err2, up->mvals2, GKYL_SUM, conf_local);
  }
  err1[0] = sqrt(err1[0]/up->conf_grid.cellVolume/conf_local->volume) / mean[1]; 
  err2[0] = sqrt(err2[0]/up->conf_grid.cellVolume/conf_local->volume) / mean[2]; 

  // Main iteration loop
  int i = 0;
  while ( (i<up->max_iter) && ((err1[0]>up->eps_err) || (err2[0]>up->eps_err)) )
  {
    // Correct the moments
    gkyl_array_clear(up->delta_m12, 0.0);
    gkyl_array_accumulate(up->delta_m12, 1.0, up->m12_tar);
    gkyl_array_accumulate(up->delta_m12, -1.0, up->m12_num);
    gkyl_array_accumulate(up->m12, 1.0, up->delta_m12);   

    gkyl_array_clear(up->mvals1, 0.0);
    gkyl_array_clear(up->mvals2, 0.0);
    gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals1, 0, up->delta_m12, *conf_local);
    gkyl_dg_calc_l2_range(up->conf_basis, 0, up->mvals2, 1, up->delta_m12, *conf_local);
    gkyl_array_scale_range(up->mvals1, up->conf_grid.cellVolume, conf_local);
    gkyl_array_scale_range(up->mvals2, up->conf_grid.cellVolume, conf_local);
    if (up->use_gpu) {
      gkyl_array_reduce_range(up->err1_cu, up->mvals1, GKYL_SUM, conf_local);
      gkyl_array_reduce_range(up->err2_cu, up->mvals2, GKYL_SUM, conf_local);
      gkyl_cu_memcpy(err1, up->err1_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
      gkyl_cu_memcpy(err2, up->err2_cu, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    } else {
      gkyl_array_reduce_range(err1, up->mvals1, GKYL_SUM, conf_local);
      gkyl_array_reduce_range(err2, up->mvals2, GKYL_SUM, conf_local);
    }
    err1[0] = sqrt(err1[0]/up->conf_grid.cellVolume/conf_local->volume) / mean[1]; 
    err2[0] = sqrt(err2[0]/up->conf_grid.cellVolume/conf_local->volume) / mean[2]; 

    // Project the maxwellian
    gkyl_array_set_offset(up->moms, 1., up->m0_num, 0*up->conf_basis.num_basis);
    gkyl_array_set_offset(up->moms, 1., up->m12, 1*up->conf_basis.num_basis);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(up->proj_maxwellian, phase_local, conf_local, up->moms, 
      up->gk_geom->bmag, up->gk_geom->bmag, up->mass, fM);
    gkyl_array_scale_by_cell(fM, up->vel_map->jacobvel);
    // Rescale the maxwellian
    gkyl_dg_updater_moment_gyrokinetic_advance(up->m0_calc, phase_local, conf_local, fM, up->m0_num);
    gkyl_dg_div_op_range(up->weak_divide, up->conf_basis, 0, up->m0_scl, 0, up->m0_tar, 0, up->m0_num, conf_local);
    gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis, fM, up->m0_scl, fM, conf_local, phase_local);
    // Calculate the moments
    gkyl_dg_updater_moment_gyrokinetic_advance(up->moms_calc, phase_local, conf_local, fM, up->moms_num);
    gkyl_array_set_offset(up->m0_num, 1., up->moms_num, 0*up->conf_basis.num_basis);
    gkyl_array_set_offset(up->m12_num, 1., up->moms_num, 1*up->conf_basis.num_basis);
    
    i += 1;
  } // Main iteration loop ends
  
  // Project maxwellian with the target moments if it fails to converge
  if (i==up->max_iter) {
    gkyl_array_set_offset(up->moms, 1., moms_tar, 0*up->conf_basis.num_basis);
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(up->proj_maxwellian, phase_local, conf_local, up->moms, 
      up->gk_geom->bmag, up->gk_geom->bmag, up->mass, fM);
    gkyl_array_scale_by_cell(fM, up->vel_map->jacobvel);
  }
}

void
gkyl_correct_maxwellian_gyrokinetic_release(gkyl_correct_maxwellian_gyrokinetic* up)
{
  if (up->use_gpu) {
    gkyl_cu_free(up->err1_cu);
    gkyl_cu_free(up->err2_cu);
    gkyl_array_release(up->mvals_host);
  }

  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_velocity_map_release(up->vel_map);

  gkyl_array_release(up->m0_tar);
  gkyl_array_release(up->m0_num);
  gkyl_array_release(up->m0_scl);
  gkyl_array_release(up->m12);
  gkyl_array_release(up->m12_tar);
  gkyl_array_release(up->m12_num);
  gkyl_array_release(up->delta_m12);
  gkyl_array_release(up->moms);
  gkyl_array_release(up->moms_num);
  gkyl_array_release(up->mvals);
  gkyl_array_release(up->mvals1);
  gkyl_array_release(up->mvals2);
  
  gkyl_dg_updater_moment_gyrokinetic_release(up->m0_calc);
  gkyl_dg_updater_moment_gyrokinetic_release(up->moms_calc);

  gkyl_proj_maxwellian_on_basis_release(up->proj_maxwellian);
  gkyl_dg_bin_op_mem_release(up->weak_divide);
  gkyl_free(up);
}
