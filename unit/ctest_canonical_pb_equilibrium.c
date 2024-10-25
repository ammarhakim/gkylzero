#include <acutest.h>

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

// Check the moments, correction routine, work on the surface of a sphere
// This ctest is intended to test:
// 1. The initial projection of a LTE moment-matches the densired values
// 2. The correction routine matches moments properly

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void 
eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void 
eval_M2(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
  fout[1] = 1.0;
}

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
test_2x2v(int poly_order)
{
  double pi = 3.14159265359;
  double lower[] = {pi/4,pi/4, -5.0, -5.0}, upper[] = {(1.01)*pi/4,(1.01)*pi/4, 5.0, 5.0};
  int cells[] = {2, 2, 16, 16};
  int vdim = 2, cdim = 2;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  double velLower[] = {lower[2], lower[3]}, velUpper[] = {upper[2], upper[3]};
  int confCells[] = {cells[0], cells[1]};
  int velCells[] = {cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1, 1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], confGhost[1], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(confBasis.num_basis*cdim*(cdim+1)/2, local_ext.volume);
  det_h = mkarr(confBasis.num_basis, local_ext.volume);
  hamil = mkarr(confBasis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, cdim*(cdim+1)/2, info_h_ij_inv, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_hamil, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &confLocal, hamil);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .h_ij_inv = h_ij_inv,
    .det_h = det_h,
    .hamil = hamil,
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .h_ij_inv = h_ij_inv,
    .det_h = det_h,
    .hamil = hamil,
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .use_gpu = false,
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  struct gkyl_vlasov_lte_correct_status status_corr;
  status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Correct the distribution function
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .h_ij_inv = h_ij_inv,
    .det_h = det_h,
    .hamil = hamil,
    .model_id = GKYL_MODEL_CANONICAL_PB,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  gkyl_array_set_offset_range(m0, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_can_pb_eq_2x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // // Write the output (moments)
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_n_corr.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,m0_corr,fname);
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_vb_corr.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,m1i_corr,fname);
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_T_corr.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,m2_corr,fname);

  // // Write the output (moments)
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_n.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,m0,fname);
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_vb.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,m1i,fname);
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_T.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,m2,fname);

  // // Write the h^{ij}, det(h_ij)
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_h_ij_inv.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,h_ij_inv,fname);
  // sprintf(fname, "ctest_can_pb_eq_2x2v_p%d_det_h.gkyl", poly_order);
  // gkyl_grid_sub_array_write(&confGrid,&confLocal,det_h,fname);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {7.4375205660623764e-02, 1.7502197646465149e-04, 
    5.1094379072269073e-19, 1.7196663881036155e-02, 3.2823939815277384e-02, 
    1.9095742813612671e-18, 4.0467707945451574e-05, -1.2193251617523399e-18, 
    1.0688988222497508e-05, 1.0430430692729768e-18, 7.5893875551799181e-03, 
    -2.2385893900551390e-07, -3.2456911333379994e-17, 7.5101204021213246e-04, 
    4.7573706725506440e-03, 2.1246392245945747e-19, 3.3851391014577862e-19, 
    2.4714545130704096e-06, -3.2130900226042711e-19, 1.9001628988953240e-18, 
    4.4692284022917053e-18, -5.1759546703890748e-08, -8.7646888322434962e-18, 
    1.7673041769548388e-06, 1.2766139686089529e-19, -1.2313454577087634e-07, 
    -1.0559775696066436e-17, 3.3144343991404806e-04, -1.3064935169987437e-05, 
    3.7101538691845707e-19, 1.0999755050985389e-03, -1.9925087386218194e-19, 
    2.1498701757639759e-18, 1.1074385582705253e-18, 1.1143514926488228e-19, 
    -8.1911466330328986e-19, 1.4840074476861625e-18, -2.8470555165622034e-08, 
    -1.1253012634920306e-18, 1.0793326594515037e-07, 1.0706453651527448e-19, 
    -2.8632875454162582e-19, -3.0208091090346236e-06, 1.1166823777084870e-19, 
    9.9107783923514388e-20, 1.1310712903906065e-18, 3.9889112700602783e-19, 
    -2.8395875891910364e-20 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[4]){1, 1, 8, 8}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-10));
      //printf("p2_vals = %1.16e fv = %1.16e\n", p2_vals[i], fv[i]);
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
}


// special note, the p1 basis does not function
void test_2x2v_p2() { test_2x2v(2); }

TEST_LIST = {
  {"test_2x2v_p2", test_2x2v_p2},
  {NULL, NULL},
};
