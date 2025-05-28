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
enum loc1d {
  L_1D, U_1D
};

struct gkyl_ten_moment_nn_closure
{
  struct gkyl_rect_grid grid; // Grid on which to solve equations.
  int ndim; // Number of dimensions.
  int poly_order; // Polynomial order of learned DG coefficients.
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
    input_data[i] = gkyl_malloc(sizeof(float[6]));
    output_data[i] = gkyl_malloc(sizeof(float[4]));
  }

  for (int i = 0; i < 3; i++) {
    input_data[i][0] = rho_avg;
    input_data[i][1] = drho_dx;
  }

  if (ndim == 1) {
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

    rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
    p_avg[0] = calc_harmonic_avg_1D(p[L_1D][0], p[U_1D][0]);
    p_avg[1] = calc_harmonic_avg_1D(p[L_1D][1], p[U_1D][1]);
    p_avg[2] = calc_harmonic_avg_1D(p[L_1D][2], p[U_1D][2]);
    p_avg[3] = calc_harmonic_avg_1D(p[L_1D][3], p[U_1D][3]);
    p_avg[4] = calc_harmonic_avg_1D(p[L_1D][4], p[U_1D][4]);
    p_avg[5] = calc_harmonic_avg_1D(p[L_1D][5], p[U_1D][5]);
    B_avg[0] = calc_harmonic_avg_1D(em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
    B_avg[1] = calc_harmonic_avg_1D(em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
    B_avg[2] = calc_harmonic_avg_1D(em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

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

  if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
    B_avg[0] = 1.0;
  }

  heat_flux_d[Q111] = q_unmag[0] * B_avg[0];
  heat_flux_d[Q112] = q_unmag[0] * B_avg[1];
  heat_flux_d[Q113] = q_unmag[0] * B_avg[2];
  heat_flux_d[Q122] = q_unmag[1] * B_avg[1];
  heat_flux_d[Q123] = q_unmag[1] * B_avg[2];
  heat_flux_d[Q133] = q_unmag[2] * B_avg[2];
  heat_flux_d[Q222] = q_unmag[3] * B_avg[1];
  heat_flux_d[Q223] = q_unmag[3] * B_avg[2];
  heat_flux_d[Q233] = q_unmag[4] * B_avg[2];
  heat_flux_d[Q333] = q_unmag[5] * B_avg[2];

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

  for (int i = 0; i < 3; i++) {
    input_data[i][0] = rho_avg;
    input_data[i][1] = drho_dx;
  }

  if (ndim == 1) {
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

    rho_avg = calc_harmonic_avg_1D(rho[L_1D], rho[U_1D]);
    p_avg[0] = calc_harmonic_avg_1D(p[L_1D][0], p[U_1D][0]);
    p_avg[1] = calc_harmonic_avg_1D(p[L_1D][1], p[U_1D][1]);
    p_avg[2] = calc_harmonic_avg_1D(p[L_1D][2], p[U_1D][2]);
    p_avg[3] = calc_harmonic_avg_1D(p[L_1D][3], p[U_1D][3]);
    p_avg[4] = calc_harmonic_avg_1D(p[L_1D][4], p[U_1D][4]);
    p_avg[5] = calc_harmonic_avg_1D(p[L_1D][5], p[U_1D][5]);
    B_avg[0] = calc_harmonic_avg_1D(em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
    B_avg[1] = calc_harmonic_avg_1D(em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
    B_avg[2] = calc_harmonic_avg_1D(em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

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

    // div(Q_11) = div(q_par), div(Q_12) = q_perp.
    divQx[0] = output_data[0][1];
    divQx[1] = output_data[0][3];

    // div(Q_33) = div(q_par), div(Q_13) = q_perp.
    divQx[5] = output_data[1][1];
    divQx[2] = output_data[1][3];

    // div(Q_22) = div(q_par), div(Q_23) = q_perp.
    divQx[3] = output_data[2][1];
    divQx[4] = output_data[2][3];
  }

  if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
    B_avg[0] = 1.0;
  }

  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = (divQx[0] * B_avg[0]) + (divQy[0] * B_avg[1]) + (divQz[0] * B_avg[2]);
  rhs[P12] = (divQx[1] * B_avg[0]) + (divQy[1] * B_avg[1]) + (divQz[1] * B_avg[2]);
  rhs[P13] = (divQx[2] * B_avg[0]) + (divQy[2] * B_avg[1]) + (divQz[2] * B_avg[2]);
  rhs[P22] = (divQx[3] * B_avg[0]) + (divQy[3] * B_avg[1]) + (divQz[3] * B_avg[2]);
  rhs[P23] = (divQx[4] * B_avg[0]) + (divQy[4] * B_avg[1]) + (divQz[4] * B_avg[2]);
  rhs[P33] = (divQx[5] * B_avg[0]) + (divQy[5] * B_avg[1]) + (divQz[5] * B_avg[2]);
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

    calc_mag_heat_flux(nnclosure, fluid_d, em_tot_d, heat_flux_d);
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
  up->ann = inp.ann;

  return up;
}

void
gkyl_ten_moment_nn_closure_release(gkyl_ten_moment_nn_closure *nnclosure)
{
  free(nnclosure);
}