#include <acutest.h>

#include <gkyl_app.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>

#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars_priv.h>
#include <gkyl_canonical_pb_kernels.h> 


void 
info_h_ij_inv(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  // [h^{xx},h^{xy},h^{yy}]
  double R = 1.0;
  double q_theta = xn[0], q_phi = xn[1];
  const double q[2] = {q_theta, q_phi};

  // [h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0 / pow(R, 2);
  fout[1] = 0.0;
  fout[2] = 1.0 / pow(R * sin(q[0]), 2);
}

void 
info_det_h(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  double R = 1.0;
  double q_theta = xn[0], q_phi = xn[1];
  const double q[2] = {q_theta, q_phi};
  fout[0] = pow(R, 2)*sin(q[0]);
}

void 
info_hamil(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double w0 = xn[2], w1 = xn[3];
  const double w[2] = {w0, w1};
  struct kh_2d_ctx *app = (struct kh_2d_ctx *)ctx;
  double *h_inv = malloc(3 * sizeof(double));
  info_h_ij_inv(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * h_inv[2] * w[1] * w[1];
  free(h_inv);
}


void 
test_2x2v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = {pi/4,pi/4, -5.0, -5.0}, upper[] = {(1.01)*pi/4,(1.01)*pi/4, 5.0, 5.0};
  int cells[] = {4, 4, 4, 4};
  int vdim = 2, cdim = 2;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  double velLower[] = {lower[2], lower[3]}, velUpper[] = {upper[2], upper[3]};
  int confCells[] = {cells[0], cells[1]};
  int velCells[] = {cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
    //printf("Choosing hybrid basis for phase space\n");
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
    //printf("Choosing serindipty basis for phase space\n");
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);
    //printf("Choosing tensor basis for velocity space\n");
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);


  int confGhost[] = {1, 1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], confGhost[1], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*cdim*(cdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, cdim*(cdim+1)/2, info_h_ij_inv, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction.
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always 2*cdim
  int alpha_surf_sz = (2*cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (2*cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (2*cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);



  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[4] = { 0.0, 0.0, 0.0, 0.0 };
  double dxv[4] = {confGrid.dx[0], confGrid.dx[1], vel_grid.dx[0], vel_grid.dx[1]};
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

  // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);


    int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    const_sgn_alpha_surf = canonical_pb_alpha_surfvx_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    const_sgn_alpha_surf = canonical_pb_alpha_surfvy_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local);

    //const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfvx_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local);  
    //const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfvy_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local);
  
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);

    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        // if(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis] != 0.0 && alpha_surf_comp_R_local[n + dir*surf_basis.num_basis] != 0.0){
        //   printf("(%1.16e) == (%1.16e)\n",alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis]);
        // }
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-12));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}


// Check the 2x2v_p continuity
void test_2x2v_p1_continuity_tensor() { test_2x2v(1, GKYL_BASIS_MODAL_TENSOR); }
void test_2x2v_p1_continuity_ser() { test_2x2v(1, GKYL_BASIS_MODAL_SERENDIPITY); }
void test_2x2v_p2_continuity_ser() { test_2x2v(2, GKYL_BASIS_MODAL_SERENDIPITY); }

TEST_LIST = {
  {"test_2x2v_p1_continuity_tensor", test_2x2v_p1_continuity_tensor},
  {"test_2x2v_p1_continuity_ser", test_2x2v_p1_continuity_ser},
  {"test_2x2v_p2_continuity_ser", test_2x2v_p2_continuity_ser},
  {NULL, NULL},
};
