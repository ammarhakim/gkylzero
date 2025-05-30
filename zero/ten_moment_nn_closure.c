#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_ten_moment_nn_closure.h>
#include <gkyl_moment_non_ideal_priv.h>

// Makes indexing cleaner.
static const unsigned Q111 = 0;
static const unsigned Q112 = 1;
static const unsigned Q113 = 2;
static const unsigned Q122 = 3;
static const unsigned Q123 = 4;
static const unsigned Q133 = 5;
static const unsigned Q222 = 6;
static const unsigned Q223 = 7;
static const unsigned Q233 = 8;
static const unsigned Q333 = 9;

// 1D stencil locations (L: lower, U: upper).
enum loc_1d {
  L_1D, U_1D
};

// 2D stencil locations (L: lower, U: upper).
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

struct gkyl_ten_moment_nn_closure
{
  struct gkyl_rect_grid grid; // Grid on which to solve equations.
  int ndim; // Number of dimensions.
  int poly_order; // Polynomial order of learned DG coefficients.
  double k0; // Damping coefficient.
  kann_t* ann; // Neural network architecture.
};

static void
create_offsets_vertices(const struct gkyl_range* range, long offsets[])
{
  // Box-spanning stencil.
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // Construct a list of offsets.
  int count = 0;
  while (gkyl_range_iter_next(&iter3)) {
    offsets[count] = gkyl_range_offset(range, iter3.idx);
    count += 1;
  }
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // Box-spanning stencil.
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // Construct a list of offsets.
  int count = 0;
  while (gkyl_range_iter_next(&iter3)) {
    offsets[count] = gkyl_range_offset(range, iter3.idx);
    count += 1;
  }
}

static void
var_setup(const gkyl_ten_moment_nn_closure *nnclosure, int start, int end, const double *fluid_d[], double rho[], double p[][6])
{
  for (int j = start; j <= end; j++) {
    rho[j] = fluid_d[j][RHO];
    p[j][0] = fluid_d[j][P11] - ((fluid_d[j][MX] * fluid_d[j][MX]) / fluid_d[j][RHO]);
    p[j][1] = fluid_d[j][P12] - ((fluid_d[j][MX] * fluid_d[j][MY]) / fluid_d[j][RHO]);
    p[j][2] = fluid_d[j][P13] - ((fluid_d[j][MX] * fluid_d[j][MZ]) / fluid_d[j][RHO]);
    p[j][3] = fluid_d[j][P22] - ((fluid_d[j][MY] * fluid_d[j][MY]) / fluid_d[j][RHO]);
    p[j][4] = fluid_d[j][P23] - ((fluid_d[j][MY] * fluid_d[j][MZ]) / fluid_d[j][RHO]);
    p[j][5] = fluid_d[j][P33] - ((fluid_d[j][MZ] * fluid_d[j][MZ]) / fluid_d[j][RHO]);
  }
}

static void
calc_mag_heat_flux(const gkyl_ten_moment_nn_closure *nnclosure, const double *fluid_d[], const double *em_tot_d[], double *heat_flux_d)
{
  const int ndim = nnclosure->ndim;
  const int poly_order = nnclosure->poly_order;
  kann_t* ann = nnclosure->ann;

  double drho_dx = 0.0, drho_dy = 0.0, drho_dz = 0.0;
  double dp_dx[6] = { 0.0 };
  double dp_dy[6] = { 0.0 };
  double dp_dz[6] = { 0.0 };
  double rho_avg = 0.0;
  double p_avg[6] = { 0.0 };
  double B_avg[3] = { 0.0 };
  double q_unmag[6] = { 0.0 };

  float **input_data = gkyl_malloc(sizeof(float*[3]));
  float **output_data = gkyl_malloc(sizeof(float*[3]));
  for (int i = 0; i < 3; i++) {
    if (ndim == 1) {
      if (poly_order == 1) {
        input_data[i] = gkyl_malloc(sizeof(float[6]));
        output_data[i] = gkyl_malloc(sizeof(float[4]));
      }
    }
    else if (ndim == 2) {
      if (poly_order == 1) {
        input_data[i] = gkyl_malloc(sizeof(float[12]));
        output_data[i] = gkyl_malloc(sizeof(float[8]));
      }
    }
  }

  if (ndim == 1) {
    if (poly_order == 1) {
      const double dx = nnclosure->grid.dx[0];
      double rho[2] = { 0.0 };
      double p[2][6] = { 0.0 };
      var_setup(nnclosure, L_1D, U_1D, fluid_d, rho, p);

      drho_dx = calc_sym_grad_1D(dx, rho[L_1D], rho[U_1D]);
      dp_dx[0] = calc_sym_grad_1D(dx, p[L_1D][0], p[U_1D][0]);
      dp_dx[1] = calc_sym_grad_1D(dx, p[L_1D][1], p[U_1D][1]);
      dp_dx[2] = calc_sym_grad_1D(dx, p[L_1D][2], p[U_1D][2]);
      dp_dx[3] = calc_sym_grad_1D(dx, p[L_1D][3], p[U_1D][3]);
      dp_dx[4] = calc_sym_grad_1D(dx, p[L_1D][4], p[U_1D][4]);
      dp_dx[5] = calc_sym_grad_1D(dx, p[L_1D][5], p[U_1D][5]);

      rho_avg = calc_arithm_avg_1D(rho[L_1D], rho[U_1D]);
      p_avg[0] = calc_arithm_avg_1D(p[L_1D][0], p[U_1D][0]);
      p_avg[1] = calc_arithm_avg_1D(p[L_1D][1], p[U_1D][1]);
      p_avg[2] = calc_arithm_avg_1D(p[L_1D][2], p[U_1D][2]);
      p_avg[3] = calc_arithm_avg_1D(p[L_1D][3], p[U_1D][3]);
      p_avg[4] = calc_arithm_avg_1D(p[L_1D][4], p[U_1D][4]);
      p_avg[5] = calc_arithm_avg_1D(p[L_1D][5], p[U_1D][5]);

      B_avg[0] = calc_arithm_avg_1D(em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
      B_avg[1] = calc_arithm_avg_1D(em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
      B_avg[2] = calc_arithm_avg_1D(em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

      for (int i = 0; i < 3; i++) {
        input_data[i][0] = rho_avg;
        input_data[i][1] = drho_dx;
      }

      // p_par = P_11, p_perp = P_12.
      input_data[0][2] = p_avg[0];
      input_data[0][3] = dp_dx[0];
      input_data[0][4] = p_avg[1];
      input_data[0][5] = dp_dx[1];

      // p_par = P_33, p_perp = P_13.
      input_data[1][2] = p_avg[5];
      input_data[1][3] = dp_dx[5];
      input_data[1][4] = p_avg[2];
      input_data[1][5] = dp_dx[2];

      // p_par = P_22, p_perp = P_23.
      input_data[2][2] = p_avg[3];
      input_data[2][3] = dp_dx[3];
      input_data[2][4] = p_avg[4];
      input_data[2][5] = dp_dx[4];

      for (int i = 0; i < 3; i++) {
        const float *output_data_predicted = kann_apply1(ann, input_data[i]);

        for (int j = 0; j < 4; j++) {
          output_data[i][j] = output_data_predicted[j];
        }
      }

      // Q_11 = q_par, Q_12 = q_perp.
      q_unmag[0] = output_data[0][0];
      q_unmag[1] = output_data[0][2];

      // Q_33 = q_par, Q_13 = q_perp.
      q_unmag[5] = output_data[1][0];
      q_unmag[2] = output_data[1][2];

      // Q_22 = q_par, Q_23 = q_perp.
      q_unmag[3] = output_data[2][0];
      q_unmag[4] = output_data[2][2];
    }
  }
  else if (ndim == 2) {
    if (poly_order == 1) {
      const double dx = nnclosure->grid.dx[0];
      const double dy = nnclosure->grid.dx[1];
      double rho[4] = { 0.0 };
      double p[4][6] = { 0.0 };
      var_setup(nnclosure, LL_2D, UU_2D, fluid_d, rho, p);

      drho_dx = calc_sym_gradx_2D(dx, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dx[0] = calc_sym_gradx_2D(dx, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dx[1] = calc_sym_gradx_2D(dx, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dx[2] = calc_sym_gradx_2D(dx, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dx[3] = calc_sym_gradx_2D(dx, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dx[4] = calc_sym_gradx_2D(dx, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dx[5] = calc_sym_gradx_2D(dx, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);

      drho_dy = calc_sym_grady_2D(dy, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dy[0] = calc_sym_grady_2D(dy, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dy[1] = calc_sym_grady_2D(dy, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dy[2] = calc_sym_grady_2D(dy, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dy[3] = calc_sym_grady_2D(dy, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dy[4] = calc_sym_grady_2D(dy, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dy[5] = calc_sym_grady_2D(dy, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);

      rho_avg = calc_arithm_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      p_avg[0] = calc_arithm_avg_2D(p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      p_avg[1] = calc_arithm_avg_2D(p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      p_avg[2] = calc_arithm_avg_2D(p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      p_avg[3] = calc_arithm_avg_2D(p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      p_avg[4] = calc_arithm_avg_2D(p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      p_avg[5] = calc_arithm_avg_2D(p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);

      B_avg[0] = calc_arithm_avg_2D(em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      B_avg[1] = calc_arithm_avg_2D(em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      B_avg[2] = calc_arithm_avg_2D(em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      for (int i = 0; i < 3; i++) {
        input_data[i][0] = rho_avg;
        input_data[i][1] = drho_dx;
        input_data[i][2] = drho_dy;
        input_data[i][3] = 0.0; // Set cross derivative to zero, for now.
      }

      // p_par = P_11, p_perp = P_12.
      input_data[0][4] = p_avg[0];
      input_data[0][5] = dp_dx[0];
      input_data[0][6] = dp_dy[0];
      input_data[0][7] = 0.0; // Set cross derivative to zero, for now.
      input_data[0][8] = p_avg[1];
      input_data[0][9] = dp_dx[1];
      input_data[0][10] = dp_dy[1];
      input_data[0][11] = 0.0; // Set cross derivative to zero, for now.

      // p_par = P_33, p_perp = P_13.
      input_data[1][4] = p_avg[5];
      input_data[1][5] = dp_dx[5];
      input_data[1][6] = dp_dy[5];
      input_data[1][7] = 0.0; // Set cross derivative to zero, for now.
      input_data[1][8] = p_avg[2];
      input_data[1][9] = dp_dx[2];
      input_data[1][10] = dp_dy[2];
      input_data[1][11] = 0.0; // Set cross derivative to zero, for now.

      // p_par = P_22, p_perp = P_23.
      input_data[2][4] = p_avg[3];
      input_data[2][5] = dp_dx[3];
      input_data[2][6] = dp_dy[3];
      input_data[2][7] = 0.0; // Set cross derivative to zero, for now.
      input_data[2][8] = p_avg[4];
      input_data[2][9] = dp_dx[4];
      input_data[2][10] = dp_dy[4];
      input_data[2][11] = 0.0; // Set cross derivative to zero, for now.

      for (int i = 0; i < 3; i++) {
        const float *output_data_predicted = kann_apply1(ann, input_data[i]);

        for (int j = 0; j < 8; j++) {
          output_data[i][j] = output_data_predicted[j];
        }
      }

      // Q_11 = q_par, Q_12 = q_perp.
      q_unmag[0] = output_data[0][0];
      q_unmag[1] = output_data[0][4];

      // Q_33 = q_par, Q_13 = q_perp.
      q_unmag[5] = output_data[1][0];
      q_unmag[2] = output_data[1][4];

      // Q_22 = q_par, Q_23 = q_perp.
      q_unmag[3] = output_data[2][0];
      q_unmag[4] = output_data[2][4];
    }
  }

  double k0 = nnclosure->k0;
  double alpha = 1.0 / k0;

  if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
    B_avg[0] = 1.0;
  }

  heat_flux_d[Q111] = alpha * q_unmag[0] * B_avg[0];
  heat_flux_d[Q112] = alpha * q_unmag[0] * B_avg[1];
  heat_flux_d[Q113] = alpha * q_unmag[0] * B_avg[2];
  heat_flux_d[Q122] = alpha * q_unmag[1] * B_avg[1];
  heat_flux_d[Q123] = alpha * q_unmag[1] * B_avg[2];
  heat_flux_d[Q133] = alpha * q_unmag[2] * B_avg[2];
  heat_flux_d[Q222] = alpha * q_unmag[3] * B_avg[1];
  heat_flux_d[Q223] = alpha * q_unmag[3] * B_avg[2];
  heat_flux_d[Q233] = alpha * q_unmag[4] * B_avg[2];
  heat_flux_d[Q333] = alpha * q_unmag[5] * B_avg[2];

  for (int i = 0; i < 3; i++) {
    gkyl_free(input_data[i]);
    gkyl_free(output_data[i]);
  }
  gkyl_free(input_data);
  gkyl_free(output_data);
}

static void
calc_nn_closure_update(const gkyl_ten_moment_nn_closure *nnclosure, const double *fluid_d[], const double *em_tot_d[], double *rhs)
{
  const int ndim = nnclosure->ndim;
  const int poly_order = nnclosure->poly_order;
  kann_t* ann = nnclosure->ann;

  double drho_dx = 0.0, drho_dy = 0.0, drho_dz = 0.0;
  double dp_dx[6] = { 0.0 };
  double dp_dy[6] = { 0.0 };
  double dp_dz[6] = { 0.0 };
  double rho_avg = 0.0;
  double p_avg[6] = { 0.0 };
  double B_avg[3] = { 0.0 };
  double divQx[6] = { 0.0 };
  double divQy[6] = { 0.0 };
  double divQz[6] = { 0.0 };

  float **input_data = gkyl_malloc(sizeof(float*[3]));
  float **output_data = gkyl_malloc(sizeof(float*[3]));
  for (int i = 0; i < 3; i++) {
    input_data[i] = gkyl_malloc(sizeof(float[6]));
    output_data[i] = gkyl_malloc(sizeof(float[4]));
  }

  if (ndim == 1) {
    if (poly_order == 1 ) {
      const double dx = nnclosure->grid.dx[0];
      double rho[2] = { 0.0 };
      double p[2][6] = { 0.0 };
      var_setup(nnclosure, L_1D, U_1D, fluid_d, rho, p);

      drho_dx = calc_sym_grad_1D(dx, rho[L_1D], rho[U_1D]);
      dp_dx[0] = calc_sym_grad_1D(dx, p[L_1D][0], p[U_1D][0]);
      dp_dx[1] = calc_sym_grad_1D(dx, p[L_1D][1], p[U_1D][1]);
      dp_dx[2] = calc_sym_grad_1D(dx, p[L_1D][2], p[U_1D][2]);
      dp_dx[3] = calc_sym_grad_1D(dx, p[L_1D][3], p[U_1D][3]);
      dp_dx[4] = calc_sym_grad_1D(dx, p[L_1D][4], p[U_1D][4]);
      dp_dx[5] = calc_sym_grad_1D(dx, p[L_1D][5], p[U_1D][5]);

      rho_avg = calc_arithm_avg_1D(rho[L_1D], rho[U_1D]);
      p_avg[0] = calc_arithm_avg_1D(p[L_1D][0], p[U_1D][0]);
      p_avg[1] = calc_arithm_avg_1D(p[L_1D][1], p[U_1D][1]);
      p_avg[2] = calc_arithm_avg_1D(p[L_1D][2], p[U_1D][2]);
      p_avg[3] = calc_arithm_avg_1D(p[L_1D][3], p[U_1D][3]);
      p_avg[4] = calc_arithm_avg_1D(p[L_1D][4], p[U_1D][4]);
      p_avg[5] = calc_arithm_avg_1D(p[L_1D][5], p[U_1D][5]);
      B_avg[0] = calc_arithm_avg_1D(em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
      B_avg[1] = calc_arithm_avg_1D(em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
      B_avg[2] = calc_arithm_avg_1D(em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

      for (int i = 0; i < 3; i++) {
        input_data[i][0] = rho_avg;
        input_data[i][1] = drho_dx;
      }

      // p_par = P_11, p_perp = P_12.
      input_data[0][2] = p_avg[0];
      input_data[0][3] = dp_dx[0];
      input_data[0][4] = p_avg[1];
      input_data[0][5] = dp_dx[1];

      // p_par = P_33, p_perp = P_13.
      input_data[1][2] = p_avg[5];
      input_data[1][3] = dp_dx[5];
      input_data[1][4] = p_avg[2];
      input_data[1][5] = dp_dx[2];

      // p_par = P_22, p_perp = P_23.
      input_data[2][2] = p_avg[3];
      input_data[2][3] = dp_dx[3];
      input_data[2][4] = p_avg[4];
      input_data[2][5] = dp_dx[4];

      for (int i = 0; i < 3; i++) {
        const float *output_data_predicted = kann_apply1(ann, input_data[i]);

        for (int j = 0; j < 4; j++) {
          output_data[i][j] = output_data_predicted[j];
        }
      }

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
      }

      // div(Q_11) = div(q_par), div(Q_12) = q_perp.
      divQx[0] = output_data[0][1] * B_avg[0];
      divQx[1] = output_data[0][3] * B_avg[0];

      // div(Q_33) = div(q_par), div(Q_13) = q_perp.
      divQx[5] = output_data[1][1] * B_avg[0];
      divQx[2] = output_data[1][3] * B_avg[0];

      // div(Q_22) = div(q_par), div(Q_23) = q_perp.
      divQx[3] = output_data[2][1] * B_avg[0];
      divQx[4] = output_data[2][3] * B_avg[0];
    }
  }
  else if (ndim == 2) {
    if (poly_order == 1) {
      const double dx = nnclosure->grid.dx[0];
      const double dy = nnclosure->grid.dx[1];
      double rho[4] = { 0.0 };
      double p[4][6] = { 0.0 };
      var_setup(nnclosure, LL_2D, UU_2D, fluid_d, rho, p);

      drho_dx = calc_sym_gradx_2D(dx, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dx[0] = calc_sym_gradx_2D(dx, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dx[1] = calc_sym_gradx_2D(dx, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dx[2] = calc_sym_gradx_2D(dx, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dx[3] = calc_sym_gradx_2D(dx, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dx[4] = calc_sym_gradx_2D(dx, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dx[5] = calc_sym_gradx_2D(dx, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);

      drho_dy = calc_sym_grady_2D(dy, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dy[0] = calc_sym_grady_2D(dy, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dy[1] = calc_sym_grady_2D(dy, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dy[2] = calc_sym_grady_2D(dy, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dy[3] = calc_sym_grady_2D(dy, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dy[4] = calc_sym_grady_2D(dy, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dy[5] = calc_sym_grady_2D(dy, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);

      rho_avg = calc_arithm_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      p_avg[0] = calc_arithm_avg_2D(p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      p_avg[1] = calc_arithm_avg_2D(p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      p_avg[2] = calc_arithm_avg_2D(p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      p_avg[3] = calc_arithm_avg_2D(p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      p_avg[4] = calc_arithm_avg_2D(p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      p_avg[5] = calc_arithm_avg_2D(p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);

      B_avg[0] = calc_arithm_avg_2D(em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      B_avg[1] = calc_arithm_avg_2D(em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      B_avg[2] = calc_arithm_avg_2D(em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      for (int i = 0; i < 3; i++) {
        input_data[i][0] = rho_avg;
        input_data[i][1] = drho_dx;
        input_data[i][2] = drho_dy;
        input_data[i][3] = 0.0; // Set cross derivative to zero, for now.
      }

      // p_par = P_11, p_perp = P_12.
      input_data[0][4] = p_avg[0];
      input_data[0][5] = dp_dx[0];
      input_data[0][6] = dp_dy[0];
      input_data[0][7] = 0.0; // Set cross derivative to zero, for now.
      input_data[0][8] = p_avg[1];
      input_data[0][9] = dp_dx[1];
      input_data[0][10] = dp_dy[1];
      input_data[0][11] = 0.0; // Set cross derivative to zero, for now.

      // p_par = P_33, p_perp = P_13.
      input_data[1][4] = p_avg[5];
      input_data[1][5] = dp_dx[5];
      input_data[1][6] = dp_dy[5];
      input_data[1][7] = 0.0; // Set cross derivative to zero, for now.
      input_data[1][8] = p_avg[2];
      input_data[1][9] = dp_dx[2];
      input_data[1][10] = dp_dy[2];
      input_data[1][11] = 0.0; // Set cross derivative to zero, for now.

      // p_par = P_22, p_perp = P_23.
      input_data[2][4] = p_avg[3];
      input_data[2][5] = dp_dx[3];
      input_data[2][6] = dp_dy[3];
      input_data[2][7] = 0.0; // Set cross derivative to zero, for now.
      input_data[2][8] = p_avg[4];
      input_data[2][9] = dp_dx[4];
      input_data[2][10] = dp_dy[4];
      input_data[2][11] = 0.0; // Set cross derivative to zero, for now.

      for (int i = 0; i < 3; i++) {
        const float *output_data_predicted = kann_apply1(ann, input_data[i]);

        for (int j = 0; j < 8; j++) {
          output_data[i][j] = output_data_predicted[j];
        }
      }

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
      }

      // div(Q_11) = div(q_par), div(Q_12) = q_perp.
      divQx[0] = output_data[0][1] * B_avg[0];
      divQy[0] = output_data[0][2] * B_avg[1];
      divQx[1] = output_data[0][5] * B_avg[0];
      divQy[1] = output_data[0][6] * B_avg[1];

      // div(Q_33) = div(q_par), div(Q_13) = q_perp.
      divQx[5] = output_data[1][1] * B_avg[0];
      divQy[5] = output_data[1][2] * B_avg[1];
      divQx[2] = output_data[1][5] * B_avg[0];
      divQy[2] = output_data[1][6] * B_avg[1];

      // div(Q_22) = div(q_par), div(Q_23) = q_perp.
      divQx[3] = output_data[2][1] * B_avg[0];
      divQy[3] = output_data[2][2] * B_avg[1];
      divQx[4] = output_data[2][5] * B_avg[0];
      divQy[4] = output_data[2][6] * B_avg[1];
    }
  }

  double k0 = nnclosure->k0;
  double alpha = 1.0 / k0;

  alpha = 0.0;
  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = alpha * (divQx[0] + divQy[0] + divQz[0]);
  rhs[P12] = alpha * (divQx[1] + divQy[1] + divQz[1]);
  rhs[P13] = alpha * (divQx[2] + divQy[2] + divQz[2]);
  rhs[P22] = alpha * (divQx[3] + divQy[3] + divQz[3]);
  rhs[P23] = alpha * (divQx[4] + divQy[4] + divQz[4]);
  rhs[P33] = alpha * (divQx[5] + divQy[5] + divQz[5]);

  for (int i = 0; i < 3; i++) {
    gkyl_free(input_data[i]);
    gkyl_free(output_data[i]);
  }
  gkyl_free(input_data);
  gkyl_free(output_data);
}

void
gkyl_ten_moment_nn_closure_advance(const gkyl_ten_moment_nn_closure *nnclosure, const struct gkyl_range *heat_flux_rng, const struct gkyl_range *update_rng,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot, struct gkyl_array *heat_flux, struct gkyl_array *rhs)
{
  int ndim = update_rng->ndim;
  long sz[] = { 2, 4, 8 };

  long offsets_vertices[sz[ndim - 1]];
  create_offsets_vertices(update_rng, offsets_vertices);

  long offsets_centers[sz[ndim - 1]];
  create_offsets_centers(heat_flux_rng, offsets_centers);

  const double *fluid_d[sz[ndim - 1]];
  const double *em_tot_d[sz[ndim - 1]];
  double *heat_flux_d;
  const double *heat_flux_up[sz[ndim - 1]];
  double *rhs_d;

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_rng);
  while (gkyl_range_iter_next(&iter_vertex)) {
    long linc_vertex = gkyl_range_idx(heat_flux_rng, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_rng, iter_vertex.idx);

    for (int i = 0; i < sz[ndim - 1]; i++) {
      em_tot_d[i] = gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices[i]);
    }

    heat_flux_d = gkyl_array_fetch(heat_flux, linc_vertex);

    // Unnecessary to calculate heat flux separately for now, so commenting this out for efficiency.
    // calc_mag_heat_flux(nnclosure, fluid_d, em_tot_d, heat_flux_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, update_rng);
  while (gkyl_range_iter_next(&iter_center)) {
    long linc_vertex = gkyl_range_idx(heat_flux_rng, iter_center.idx);
    long linc_center = gkyl_range_idx(update_rng, iter_center.idx);

    for (int i = 0; i < sz[ndim - 1]; i++) {
      em_tot_d[i] = gkyl_array_cfetch(em_tot, linc_vertex + offsets_centers[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_vertex + offsets_centers[i]);
    }

    rhs_d = gkyl_array_fetch(rhs, linc_center);

    calc_nn_closure_update(nnclosure, fluid_d, em_tot_d, rhs_d);
  }
}

gkyl_ten_moment_nn_closure*
gkyl_ten_moment_nn_closure_new(struct gkyl_ten_moment_nn_closure_inp inp)
{
  gkyl_ten_moment_nn_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_nn_closure));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp.k0;
  up->poly_order = inp.poly_order;
  up->ann = inp.ann;

  return up;
}

void
gkyl_ten_moment_nn_closure_release(gkyl_ten_moment_nn_closure *nnclosure)
{
  free(nnclosure);
}