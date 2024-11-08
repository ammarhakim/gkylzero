#include <acutest.h>
#include <math.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_zsurf.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_fem_poisson.h>

#include <gkyl_fem_parproj.h>

#include <gkyl_deflated_dg_bin_ops.h>


void
proj_jac(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = 1.0 + z*z*cos(x);
}

void
proj_rho(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double z = xn[1];
  fout[0] = z*sin((2.*M_PI/(2.*M_PI))*x);
}

// Check continuity along last dim in 2x
void check_continuity_2x(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field)
{
  struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, grid.ndim, basis.num_basis);
  basis.node_list(gkyl_array_fetch(nodes, 0));
  // Check continuity
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  int nidx;
  long lin_nidx;
  int idx[3];
  const double *node_i  ;

  while (gkyl_range_iter_next(&iter)) {
    if (iter.idx[1] != range.upper[1]) {
      long lidx = gkyl_range_idx(&range, iter.idx);
      idx[0] = iter.idx[0];
      idx[1] = iter.idx[1] + 1;
      long lidx_up = gkyl_range_idx(&range, idx);
      double *arr = gkyl_array_fetch(field, lidx);
      double *arr_up = gkyl_array_fetch(field, lidx_up);
      node_i  = gkyl_array_cfetch(nodes, 2);
      double temp1 = basis.eval_expand(node_i, arr);
      node_i  = gkyl_array_cfetch(nodes, 3);
      double temp2 = basis.eval_expand(node_i, arr);
      node_i  = gkyl_array_cfetch(nodes, 0);
      double temp_up1 = basis.eval_expand(node_i, arr_up);
      node_i  = gkyl_array_cfetch(nodes, 1);
      double temp_up2 = basis.eval_expand(node_i, arr_up);
      TEST_CHECK( gkyl_compare(temp1, temp_up1, 1e-12) );
      TEST_CHECK( gkyl_compare(temp2, temp_up2, 1e-12) );
    }
  }
  gkyl_array_release(nodes);
}

void check_same(struct gkyl_range range, struct gkyl_basis basis, struct gkyl_array *field1, struct gkyl_array* field2)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&range, iter.idx);
    const double *f1 = gkyl_array_cfetch(field1, lidx);
    const double *f2 = gkyl_array_cfetch(field2, lidx);
    for(int i = 0; i< basis.num_basis; i++)
      TEST_CHECK( gkyl_compare(f1[i], f2[i], 1e-10) );
  }
}



void
test_aop()
{
  int c_lop = 0;
  int c_rop = 0;
  int c_oop = 0;
  // create the 2d field
  // create xz grid
  double lower[] = { -M_PI, 0.0 }, upper[] = { M_PI, 1.0 };
  int cells[] = { 12, 8 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  //ranges
  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);
  struct gkyl_range global_sub_range = local;

  // basis function
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, poly_order);

  bool use_gpu = false;
#ifdef GKYL_HAVE_CUDA
  use_gpu = true;
  struct gkyl_basis *basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis_on_dev, 2, poly_order);
#else
  struct gkyl_basis *basis_on_dev = &basis;
#endif

  // project initial functions on 2d field
  struct gkyl_array *rho = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&grid, &basis, 2, 1, &proj_rho, 0);
  gkyl_proj_on_basis_advance(proj, 0.0, &local, rho);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_array *jac= gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&grid, &basis, 1, &proj_jac, 0);
  gkyl_eval_on_nodes_advance(eon, 0.0, &local, jac);
  gkyl_eval_on_nodes_release(eon);

  gkyl_grid_sub_array_write(&grid, &local, 0, rho, "rho.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, 0, jac, "jac.gkyl");

  // Get C = rho.J
  struct gkyl_array *Cxz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_mul_op(basis, 0, Cxz, 0,rho, 0, jac);
  gkyl_grid_sub_array_write(&grid, &local, 0, Cxz, "Cxz.gkyl");

#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *Cxz_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(Cxz_dev, Cxz);
  struct gkyl_array *jac_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_copy(jac_dev, jac);
#else
  struct gkyl_array *Cxz_dev = Cxz;
  struct gkyl_array *jac_dev = jac;
#endif

  // Allocate Exz which will store C/J, the smoothed version of it, and
  // Fxz = Exz*J which will hold the final result
  struct gkyl_array *Exz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *Exz_smooth = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *Fxz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#ifdef GKYL_HAVE_CUDA
  struct gkyl_array *Exz_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *Exz_smooth_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *Fxz_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
#else
  struct gkyl_array *Exz_dev = Exz;
  struct gkyl_array *Exz_smooth_dev = Exz_smooth;
  struct gkyl_array *Fxz_dev = Fxz;
#endif

  // Make deflated array operator
  struct gkyl_deflated_dg_bin_ops* operator = gkyl_deflated_dg_bin_ops_new(grid, basis_on_dev, basis, local, use_gpu);

  // Divide
  gkyl_deflated_dg_bin_ops_div(operator, c_oop, Exz_dev, c_lop, Cxz_dev, c_rop, jac_dev);

#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(Exz, Exz_dev);
#endif
  gkyl_grid_sub_array_write(&grid, &local, 0, Exz, "Exz.gkyl");

  // Smooth E
  struct gkyl_fem_parproj *fem_parproj = gkyl_fem_parproj_new(&local, &local_ext, &basis, GKYL_FEM_PARPROJ_NONE, NULL, use_gpu);
  gkyl_fem_parproj_set_rhs(fem_parproj, Exz_dev, Exz_dev);
  gkyl_fem_parproj_solve(fem_parproj, Exz_smooth_dev);
  gkyl_fem_parproj_release(fem_parproj);

#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(Exz_smooth, Exz_smooth_dev);
#endif
  gkyl_grid_sub_array_write(&grid, &local, 0, Exz_smooth, "Exz_smooth.gkyl");

  // Multiply
  gkyl_deflated_dg_bin_ops_mul(operator, c_oop, Fxz_dev, c_lop, Exz_smooth_dev, c_rop, jac_dev);
  
#ifdef GKYL_HAVE_CUDA
  gkyl_array_copy(Fxz, Fxz_dev);
#endif
  gkyl_grid_sub_array_write(&grid, &local, 0, Fxz, "Fxz.gkyl");

  check_continuity_2x(grid, local, basis, Fxz);
  
#ifdef GKYL_HAVE_CUDA 
  gkyl_cu_free(basis_on_dev);
  gkyl_array_release(Cxz_dev);
  gkyl_array_release(jac_dev);
  gkyl_array_release(Exz_dev);
  gkyl_array_release(Exz_smooth_dev);
  gkyl_array_release(Fxz_dev);
#endif
  gkyl_array_release(rho);
  gkyl_array_release(jac);
  gkyl_array_release(Cxz);
  gkyl_array_release(Exz);
  gkyl_array_release(Exz_smooth);
  gkyl_array_release(Fxz);
  gkyl_deflated_dg_bin_ops_release(operator);

}

TEST_LIST = {
  { "test_aop", test_aop},
  { NULL, NULL },
};
