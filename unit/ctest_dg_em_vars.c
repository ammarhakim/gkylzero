#include "math.h"

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <time.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_field_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double Lx = 2.0*M_PI;
  double kx = 2.0*M_PI/Lx;

  pcg64_random_t rng = gkyl_pcg64_init(0); // RNG for use in IC

  double Ex = 1.0;
  double Ey = 0.0;
  double Ez = 0.0;
  double Bx = 1.0;
  double By = 0.0;
  double Bz = 0.0;
  for (int i=0; i<4; ++i) {
    Ey += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    Ez += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    By += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    Bz += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
  }

  fout[0] = Ex;
  fout[1] = Ey;
  fout[2] = Ez;
  fout[3] = Bx;
  fout[4] = By;
  fout[5] = Bz;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void eval_field_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double Lx = 2.0*M_PI;
  double Ly = 2.0*M_PI;
  double kx = 2.0*M_PI/Lx;
  double ky = 2.0*M_PI/Ly;

  pcg64_random_t rng = gkyl_pcg64_init(0); // RNG for use in IC

  double Ex = 0.0;
  double Ey = 0.0;
  double Ez = 0.0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      Ex += gkyl_pcg64_rand_double(&rng)*sin(j*kx*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      Ey += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      Ez += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(j*ky*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      Bx += gkyl_pcg64_rand_double(&rng)*sin(j*kx*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      By += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      Bz += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(i*kx*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
    }
  }

  fout[0] = Ex;
  fout[1] = Ey;
  fout[2] = Ez;
  fout[3] = Bx;
  fout[4] = By;
  fout[5] = Bz;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void eval_field_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  double Lx = 2.0*M_PI;
  double Ly = 2.0*M_PI;
  double Lz = 2.0*M_PI;
  double kx = 2.0*M_PI/Lx;
  double ky = 2.0*M_PI/Ly;
  double kz = 2.0*M_PI/Ly;

  pcg64_random_t rng = gkyl_pcg64_init(0); // RNG for use in IC

  double Ex = 0.0;
  double Ey = 0.0;
  double Ez = 0.0;
  double Bx = 0.0;
  double By = 0.0;
  double Bz = 0.0;
  for (int i=0; i<2; ++i) {
    for (int j=0; j<2; ++j) {
      for (int k=0; k<2; ++k) {
        Ex += gkyl_pcg64_rand_double(&rng)*sin(j*kx*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(k*kz*z + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
        Ey += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(k*kz*z + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
        Ez += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(j*ky*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
        Bx += gkyl_pcg64_rand_double(&rng)*sin(j*kx*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(k*kz*z + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
        By += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(k*kz*z + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
        Bz += gkyl_pcg64_rand_double(&rng)*sin(i*kx*x + 2.0*M_PI*gkyl_pcg64_rand_double(&rng))*sin(i*kx*y + 2.0*M_PI*gkyl_pcg64_rand_double(&rng));
      }
    }
  }

  fout[0] = Ex;
  fout[1] = Ey;
  fout[2] = Ez;
  fout[3] = Bx;
  fout[4] = By;
  fout[5] = Bz;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

void
test(int ndim, int Nx, int poly_order, double eps, bool use_gpu)
{
  double L = 2.*M_PI;

  double lower[ndim], upper[ndim];
  int cells[ndim], ghost[ndim];
  for (int n=0; n<ndim; ++n) {
    lower[n] = -L/2.0;
    upper[n] = L/2.0;
    cells[n] = Nx;
    ghost[n] = 1;
  }
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Local, local-ext phase-space ranges.
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_proj_on_basis *proj_field;
  if (ndim == 1) {
    proj_field = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 8, eval_field_1x, NULL);
  }
  else if (ndim == 2) {
    proj_field = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 8, eval_field_2x, NULL);    
  }
  else {
    proj_field = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 8, eval_field_3x, NULL);       
  }

  // Create EM, bvar, and ExB arrays.
  struct gkyl_array *field, *bvar, *ExB;
  field = mkarr(8*basis.num_basis, local_ext.volume);
  bvar = mkarr(9*basis.num_basis, local_ext.volume);
  ExB = mkarr(3*basis.num_basis, local_ext.volume);

  gkyl_proj_on_basis_advance(proj_field, 0.0, &local, field);

  struct gkyl_array *field_cu, *bvar_cu, *ExB_cu;
  if (use_gpu) { // Create device copies
    field_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, 8*basis.num_basis, local_ext.volume);
    bvar_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
    ExB_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, local_ext.volume);
  }

  if (use_gpu) {
    // Copy host array to device.
    gkyl_array_copy(field_cu , field );
  }

  struct timespec tm = gkyl_wall_clock();

  // Updaters to compute bvar and ExB with a single loop over the grid
  // Updaters compute |B|^2, 1/|B|^2, and other quantities simultaneously
  // and then construct the final bvar
  // magnetic field unit vector (b_i = B_i/|B|, three components) 
  // and magnetic unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components)
  // and ExB: ExB velocity (E x B/|B|^2, 3 components)
  if (use_gpu) {
    gkyl_calc_em_vars_bvar(basis, &local, field_cu, bvar_cu);
    gkyl_calc_em_vars_ExB(basis, &local, field_cu, ExB_cu);
    // Copy host array to device.
    gkyl_array_copy(bvar , bvar_cu );
    gkyl_array_copy(ExB , ExB_cu );
  }
  else {
    gkyl_calc_em_vars_bvar(basis, &local, field, bvar);
    gkyl_calc_em_vars_ExB(basis, &local, field, ExB);
  }

  double em_tm = gkyl_time_diff_now_sec(tm);

  /* printf("\nEM variable computation on (%d)^%d took %g sec\n", cells[0], ndim,  */
  /*   em_tm); */

  // Check if b . b = 1 from EM vars computation
  struct gkyl_array *bibj_check, *b_dot_b;
  bibj_check = mkarr(3*basis.num_basis, local_ext.volume);
  b_dot_b = mkarr(basis.num_basis, local_ext.volume);
  for (int i=0; i<3; ++i) {
    gkyl_dg_mul_op_range(basis, i, bibj_check, i, bvar, i, bvar, &local);
  }
  gkyl_array_accumulate_offset_range(b_dot_b, 1.0, bibj_check, 0*basis.num_basis, local);
  gkyl_array_accumulate_offset_range(b_dot_b, 1.0, bibj_check, 1*basis.num_basis, local);
  gkyl_array_accumulate_offset_range(b_dot_b, 1.0, bibj_check, 2*basis.num_basis, local);

  // Create intermediate arrays and dg_bin_op_memory to construct bvar
  // and ExB by the relevant sequence of operations
  struct gkyl_array *magB2, *int_BiBj, *int_ExB1, *int_ExB2, *alt_bibj, *alt_ExB;
  struct gkyl_array *alt_bibj_cu, *alt_ExB_cu;
  struct gkyl_dg_bin_op_mem *magB2_mem;

  alt_bibj = mkarr(9*basis.num_basis, local_ext.volume);
  alt_ExB = mkarr(3*basis.num_basis, local_ext.volume);

  struct timespec tm2 = gkyl_wall_clock();
  if (use_gpu) {
    alt_bibj_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 9*basis.num_basis, local_ext.volume);
    alt_ExB_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, local_ext.volume);

    magB2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    int_BiBj = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*basis.num_basis, local_ext.volume);
    int_ExB1 = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, local_ext.volume);
    int_ExB2 = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*basis.num_basis, local_ext.volume);

    magB2_mem = gkyl_dg_bin_op_mem_cu_dev_new(local.volume, basis.num_basis);
    int ctr = 0;
    for (int i=0; i<3; ++i) {
      for (int j=i; j<3; ++j) {
        gkyl_dg_mul_op_range(basis, ctr, int_BiBj, i+3, field_cu, j+3, field_cu, &local);
        ctr += 1;
      }
    }
    gkyl_array_accumulate_offset_range(magB2, 1.0, int_BiBj, 0*basis.num_basis, local);
    gkyl_array_accumulate_offset_range(magB2, 1.0, int_BiBj, 3*basis.num_basis, local);
    gkyl_array_accumulate_offset_range(magB2, 1.0, int_BiBj, 5*basis.num_basis, local);

    for (int i=0; i<6; ++i) {
      gkyl_dg_div_op_range(magB2_mem, basis, 3+i, alt_bibj_cu, i, int_BiBj, 0, magB2, &local);
    }

    gkyl_dg_mul_op_range(basis, 0, int_ExB1, 1, field_cu, 5, field_cu, &local);
    gkyl_dg_mul_op_range(basis, 0, int_ExB2, 2, field_cu, 4, field_cu, &local);

    gkyl_dg_mul_op_range(basis, 0, int_ExB1, 2, field_cu, 3, field_cu, &local);
    gkyl_dg_mul_op_range(basis, 0, int_ExB2, 0, field_cu, 5, field_cu, &local);

    gkyl_dg_mul_op_range(basis, 0, int_ExB1, 0, field_cu, 4, field_cu, &local);
    gkyl_dg_mul_op_range(basis, 0, int_ExB2, 1, field_cu, 3, field_cu, &local);

    gkyl_array_accumulate_range(int_ExB1, -1.0, int_ExB2, local);
    for (int i=0; i<3; ++i) {
      gkyl_dg_div_op_range(magB2_mem, basis, i, alt_ExB_cu, i, int_ExB1, 0, magB2, &local);
    }    

    // copy from device and check if things are ok
    gkyl_array_copy(alt_bibj, alt_bibj_cu);
    gkyl_array_copy(alt_ExB, alt_ExB_cu);
  }
  else {
    magB2 = mkarr(basis.num_basis, local_ext.volume);
    int_BiBj = mkarr(6*basis.num_basis, local_ext.volume);
    int_ExB1 = mkarr(3*basis.num_basis, local_ext.volume);
    int_ExB2 = mkarr(3*basis.num_basis, local_ext.volume);

    magB2_mem = gkyl_dg_bin_op_mem_new(local.volume, basis.num_basis);

    int ctr = 0;
    for (int i=0; i<3; ++i) {
      for (int j=i; j<3; ++j) {
        gkyl_dg_mul_op_range(basis, ctr, int_BiBj, i+3, field, j+3, field, &local);
        ctr += 1;
      }
    }
    gkyl_array_accumulate_offset_range(magB2, 1.0, int_BiBj, 0*basis.num_basis, local);
    gkyl_array_accumulate_offset_range(magB2, 1.0, int_BiBj, 3*basis.num_basis, local);
    gkyl_array_accumulate_offset_range(magB2, 1.0, int_BiBj, 5*basis.num_basis, local);

    for (int i=0; i<6; ++i) {
      gkyl_dg_div_op_range(magB2_mem, basis, 3+i, alt_bibj, i, int_BiBj, 0, magB2, &local);
    }

    gkyl_dg_mul_op_range(basis, 0, int_ExB1, 1, field, 5, field, &local);
    gkyl_dg_mul_op_range(basis, 0, int_ExB2, 2, field, 4, field, &local);

    gkyl_dg_mul_op_range(basis, 1, int_ExB1, 2, field, 3, field, &local);
    gkyl_dg_mul_op_range(basis, 1, int_ExB2, 0, field, 5, field, &local);

    gkyl_dg_mul_op_range(basis, 2, int_ExB1, 0, field, 4, field, &local);
    gkyl_dg_mul_op_range(basis, 2, int_ExB2, 1, field, 3, field, &local);

    gkyl_array_accumulate_range(int_ExB1, -1.0, int_ExB2, local);
    for (int i=0; i<3; ++i) {
      gkyl_dg_div_op_range(magB2_mem, basis, i, alt_ExB, i, int_ExB1, 0, magB2, &local);
    }    
  }

  double em_2_tm = gkyl_time_diff_now_sec(tm2);

  /* printf("dg_bin_op EM variable computation on (%d)^%d took %g sec\n", cells[0], ndim,  */
  /*   em_2_tm); */

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);

    // Check magnetic field unit tensor
    const double *bvar_p = gkyl_array_cfetch(bvar, linidx);
    const double *alt_bvar_p = gkyl_array_cfetch(alt_bibj, linidx);
    for (int m=3*basis.num_basis; m<9*basis.num_basis; ++m) {
      TEST_CHECK( gkyl_compare(alt_bvar_p[m], bvar_p[m], eps) );
      if (ndim == 1)
        TEST_MSG("Expected: %.13e, coefficient (%d) in cell (%d)", alt_bvar_p[m], m, iter.idx[0]);
      else if (ndim == 2)
        TEST_MSG("Expected: %.13e, coefficient (%d) in cell (%d, %d)", alt_bvar_p[m], m, iter.idx[0], iter.idx[1]);
      else
        TEST_MSG("Expected: %.13e, coefficient (%d) in cell (%d, %d, %d)", alt_bvar_p[m], m, iter.idx[0], iter.idx[1], iter.idx[2]);
      TEST_MSG("Produced: %.13e, coefficient (%d)", bvar_p[m], m);
    }
    
    // Check if B . B/|B|^2 = 1
    TEST_CHECK( gkyl_compare(alt_bvar_p[3*basis.num_basis] + alt_bvar_p[6*basis.num_basis] + alt_bvar_p[8*basis.num_basis], 
      bvar_p[3*basis.num_basis] + bvar_p[6*basis.num_basis] + bvar_p[8*basis.num_basis], 1.0e-14) );
      if (ndim == 1) 
        TEST_MSG("Expected: %.13e in cell (%d)", sqrt(2.0), iter.idx[0]);
      else if (ndim == 2)
        TEST_MSG("Expected: %.13e in cell (%d, %d)", 2.0, iter.idx[0], iter.idx[1]);
      else
        TEST_MSG("Expected: %.13e in cell (%d, %d, %d)", 2.0*sqrt(2.0), iter.idx[0], iter.idx[1], iter.idx[2]);

      TEST_MSG("Cell average B . B/|B|^2 produced by EM vars computation: %.13e", bvar_p[3*basis.num_basis] + bvar_p[6*basis.num_basis] + bvar_p[8*basis.num_basis]);
      TEST_MSG("Cell average B . B/|B|^2 Produced by dg_bin_op: %.13e", alt_bvar_p[3*basis.num_basis] + alt_bvar_p[6*basis.num_basis] + alt_bvar_p[8*basis.num_basis]);

    // Check b . b = 1
    const double *b_dot_b_p = gkyl_array_cfetch(b_dot_b, linidx);
    TEST_CHECK( gkyl_compare(b_dot_b_p[0], pow(2.0, ndim/2.0), 1.0e-14) );
    TEST_MSG("b . b cell average from EM vars computation: %.13e", b_dot_b_p[0]);

    // Check E x B velocity
    const double *ExB_p = gkyl_array_cfetch(ExB, linidx);
    const double *alt_ExB_p = gkyl_array_cfetch(alt_ExB, linidx);
    for (int m=0; m<3*basis.num_basis; ++m) {
      TEST_CHECK( gkyl_compare(alt_ExB_p[m], ExB_p[m], eps) );
      if (ndim == 1)
        TEST_MSG("Expected: %.13e, coefficient (%d) in cell (%d)", alt_ExB_p[m], m, iter.idx[0]);
      else if (ndim == 2)
        TEST_MSG("Expected: %.13e, coefficient (%d) in cell (%d, %d)", alt_ExB_p[m], m, iter.idx[0], iter.idx[1]);
      else
        TEST_MSG("Expected: %.13e, coefficient (%d) in cell (%d, %d, %d)", alt_ExB_p[m], m, iter.idx[0], iter.idx[1], iter.idx[2]);
      TEST_MSG("Produced: %.13e, coefficient (%d)", ExB_p[m], m);
    }
  }

  gkyl_array_release(field);
  gkyl_array_release(bvar);
  gkyl_array_release(ExB);

  gkyl_array_release(bibj_check);
  gkyl_array_release(b_dot_b);

  gkyl_array_release(alt_bibj);
  gkyl_array_release(alt_ExB);
  gkyl_array_release(magB2);
  gkyl_array_release(int_BiBj);
  gkyl_array_release(int_ExB1);
  gkyl_array_release(int_ExB2);

  gkyl_dg_bin_op_mem_release(magB2_mem);
  if (use_gpu) {
    gkyl_array_release(field_cu);
    gkyl_array_release(bvar_cu);
    gkyl_array_release(ExB_cu);

    gkyl_array_release(alt_bibj_cu);
    gkyl_array_release(alt_ExB_cu);
  }  

  gkyl_proj_on_basis_release(proj_field);
}

void test_1x_p1() { test(1, 8, 1, 1.0e-12, false); }
void test_2x_p1() { test(2, 8, 1, 1.0e-11, false); }
void test_3x_p1() { test(3, 8, 1, 1.0e-10, false); }

void test_2x_p1_big() { test(2, 512, 1, 1.0e-9, false); }
void test_3x_p1_big() { test(3, 128, 1, 1.0e-8, false); }

void test_1x_p2() { test(1, 512, 2, 1.0e-4, false); }
void test_2x_p2() { test(2, 512, 2, 1.0e-2, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x_p1_gpu() { test(1, 8, 1, 1.0e-9, true); }
void test_2x_p1_gpu() { test(2, 8, 1, 1.0e-9, true); }
void test_3x_p1_gpu() { test(3, 8, 1, 1.0e-9, true); }

void test_2x_p1_big_gpu() { test(2, 512, 1, 1.0e-9, true); }
void test_3x_p1_big_gpu() { test(3, 128, 1, 1.0e-9, true); }

void test_1x_p2_gpu() { test(1, 512, 2, 1.0e-4, true); }
void test_2x_p2_gpu() { test(2, 512, 2, 1.0e-4, true); }


#endif

TEST_LIST = {
  { "test_1x_p1", test_1x_p1 },
  { "test_2x_p1", test_2x_p1 },
  { "test_3x_p1", test_3x_p1 },

  // { "test_2x_p1_big", test_2x_p1_big },
  // { "test_3x_p1_big", test_3x_p1_big },

  { "test_1x_p2", test_1x_p2 },
  // { "test_2x_p2", test_2x_p2 },

#ifdef GKYL_HAVE_CUDA
  { "test_1x_p1_gpu", test_1x_p1_gpu },
  { "test_2x_p1_gpu", test_2x_p1_gpu },
  { "test_3x_p1_gpu", test_3x_p1_gpu },

  // { "test_2x_p1_big_gpu", test_2x_p1_big_gpu },
  // { "test_3x_p1_big_gpu", test_3x_p1_big_gpu },

  { "test_1x_p2_gpu", test_1x_p2_gpu },
  // { "test_2x_p2_gpu", test_2x_p2_gpu },

#endif
  { NULL, NULL },
};
