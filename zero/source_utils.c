#include <gkyl_source_utils.h>

void
pressure_tensor_rotate(double q_over_m, double dt, const double* em, const double* ext_em, double p_tensor_old[6], double p_tensor_rhs[6],
  double p_tensor_new[6])
{
  double Bx = em[3] + ext_em[3];
  double By = em[4] + ext_em[4];
  double Bz = em[5] + ext_em[5];

  double dt1 = 0.5 * dt;
  double dt1_sq = dt1 * dt1;
  double dt1_cu = dt1_sq * dt1;
  double dt1_qu = dt1_cu * dt1;

  double Bx_sq = Bx * Bx;
  double Bx_cu = Bx_sq * Bx;
  double Bx_qu = Bx_cu * Bx;

  double By_sq = By * By;
  double By_cu = By_sq * By;
  double By_qu = By_cu * By;

  double Bz_sq = Bz * Bz;
  double Bz_cu = Bz_sq * Bz;
  double Bz_qu = Bz_cu * Bz;

  double q_over_m_sq = q_over_m * q_over_m;
  double q_over_m_cu = q_over_m_sq * q_over_m;
  double q_over_m_qu = q_over_m_cu * q_over_m;

  double denom = 1.0 + (5.0 * (Bx_sq + By_sq + Bz_sq) * dt1_sq * q_over_m_sq) + (4.0 * (Bx_sq + By_sq + Bz_sq) * (Bx_sq + By_sq + Bz_sq) * dt1_qu * q_over_m_qu);

  p_tensor_new[0] = 2.0 * (p_tensor_rhs[0] + 2.0 * dt1 * (Bz * p_tensor_rhs[1] - By * p_tensor_rhs[2]) * q_over_m + dt1_sq * (5.0 * Bx_sq * p_tensor_rhs[0] +
    2.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + Bz_sq * (3.0 * p_tensor_rhs[0] + 2.0 * p_tensor_rhs[3]) - 4.0 * By * Bz * p_tensor_rhs[4] +
    By_sq * (3.0 * p_tensor_rhs[0] + 2.0 * p_tensor_rhs[5])) * q_over_m_sq + 2.0 * dt1_cu * (4.0 * Bx_sq * (Bz * p_tensor_rhs[1] - By * p_tensor_rhs[2]) -
    (By_sq + Bz_sq) * (-(Bz * p_tensor_rhs[1]) + By * p_tensor_rhs[2]) - 3.0 * Bx * (By_sq * p_tensor_rhs[4] - Bz_sq * p_tensor_rhs[4] + By * Bz *
    (-p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * dt1_qu * (2.0 * Bx_qu * p_tensor_rhs[0] + 4.0 * Bx_cu * (By * p_tensor_rhs[1] +
    Bz * p_tensor_rhs[2]) - 2.0 * Bx * (By_sq + Bz_sq) * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + (By_sq + Bz_sq) * (Bz_sq * (p_tensor_rhs[0] +
    p_tensor_rhs[3]) - 2.0 * By * Bz * p_tensor_rhs[4] + By_sq * (p_tensor_rhs[0] + p_tensor_rhs[5])) + Bx_sq * (4.0 * By * Bz * p_tensor_rhs[4] + By_sq * (3.0 *
    p_tensor_rhs[3] + p_tensor_rhs[5]) + Bz_sq * (p_tensor_rhs[3] + 3.0 * p_tensor_rhs[5]))) * q_over_m_qu) / denom - p_tensor_old[0];

  p_tensor_new[1] = 2.0 * (p_tensor_rhs[1] + dt1 * (Bx * p_tensor_rhs[2] + Bz * (-p_tensor_rhs[0] + p_tensor_rhs[3]) - By * p_tensor_rhs[4]) * q_over_m + dt1_sq *
    (4.0 * Bx_sq * p_tensor_rhs[1] + 4.0 * By_sq * p_tensor_rhs[1] + Bz_sq * p_tensor_rhs[1] + 3.0 * By * Bz * p_tensor_rhs[2] + Bx * (3.0 * Bz * p_tensor_rhs[4] +
    By * (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]))) * q_over_m_sq + dt1_cu * (4.0 * Bx_cu * p_tensor_rhs[2] - 2.0 * Bx * (By_sq + Bz_sq) * 
    p_tensor_rhs[2] + Bz_cu * (-p_tensor_rhs[0] + p_tensor_rhs[3]) - 4.0 * By_cu * p_tensor_rhs[4] + 2.0 * By * Bz_sq * p_tensor_rhs[4] - By_sq * Bz *
    (p_tensor_rhs[0] - 4.0 * p_tensor_rhs[3] + 3.0 * p_tensor_rhs[5]) + Bx_sq * (2.0 * By * p_tensor_rhs[4] + Bz * (-4.0 * p_tensor_rhs[0] + p_tensor_rhs[3] +
    3.0 * p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * Bx * By * dt1_qu * (6.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 6.0 * By * Bz * p_tensor_rhs[4] -
    Bz_sq * (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]) + Bx_sq * (2.0 * p_tensor_rhs[0] - p_tensor_rhs[3] - p_tensor_rhs[5]) - By_sq *
    (p_tensor_rhs[0] - 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m_qu) / denom - p_tensor_old[1];

  p_tensor_new[2] = 2.0 * (p_tensor_rhs[2] + dt1 * (-(Bx * p_tensor_rhs[1]) + Bz * p_tensor_rhs[4] + By * (p_tensor_rhs[0] - p_tensor_rhs[5])) * q_over_m + dt1_sq *
    (3.0 * By * Bz * p_tensor_rhs[1] + 4.0 * Bx_sq * p_tensor_rhs[2] + By_sq * p_tensor_rhs[2] + 4.0 * Bz_sq * p_tensor_rhs[2] + Bx * (3.0 * By * p_tensor_rhs[4] +
    Bz * (p_tensor_rhs[0] - 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_sq + dt1_cu * (-4.0 * Bx_cu * p_tensor_rhs[1] + 2.0 * Bx * (By_sq + Bz_sq) *
    p_tensor_rhs[1] - 2.0 * By_sq * Bz * p_tensor_rhs[4] + 4.0 * Bz_cu * p_tensor_rhs[4] + By * Bz_sq * (p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3] - 4.0 *
    p_tensor_rhs[5]) + By_cu * (p_tensor_rhs[0] - p_tensor_rhs[5]) - Bx_sq * (2.0 * Bz * p_tensor_rhs[4] + By * (-4.0 * p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3] +
    p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * Bx * Bz * dt1_qu * (6.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 6.0 * By * Bz * p_tensor_rhs[4] -
    Bz_sq * (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]) + Bx_sq * (2.0 * p_tensor_rhs[0] - p_tensor_rhs[3] - p_tensor_rhs[5]) - By_sq * (p_tensor_rhs[0] -
    2.0 * p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m_qu) / denom - p_tensor_old[2];
  
  p_tensor_new[3] = 2.0 * (p_tensor_rhs[3] + (-2.0 * Bz * dt1 * p_tensor_rhs[1] + 2.0 * Bx * dt1 * p_tensor_rhs[4]) * q_over_m + dt1_sq * (2.0 * Bx * By *
    p_tensor_rhs[1] + 5.0 * By_sq * p_tensor_rhs[3] + Bz_sq * (2.0 * p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3]) + Bz * (-4.0 * Bx * p_tensor_rhs[2] + 2.0 * By *
    p_tensor_rhs[4]) + Bx_sq * (3.0 * p_tensor_rhs[3] + 2.0 * p_tensor_rhs[5])) * q_over_m_sq + 2.0 * dt1_cu * (Bx_sq * (-(Bz * p_tensor_rhs[1]) + 3.0 * By *
    p_tensor_rhs[2]) - Bz * (4.0 * By_sq * p_tensor_rhs[1] + Bz_sq * p_tensor_rhs[1] + 3.0 * By * Bz * p_tensor_rhs[2]) + Bx_cu * p_tensor_rhs[4] + Bx *
    (4.0 * By_sq * p_tensor_rhs[4] + Bz_sq * p_tensor_rhs[4] + 3.0 * By * Bz * (-p_tensor_rhs[0] + p_tensor_rhs[5]))) * q_over_m_cu + 2.0 * dt1_qu * (-2.0 *
    Bx_cu * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 2.0 * Bx * (2.0 * By_sq - Bz_sq) * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + 2.0 * By_qu *
    p_tensor_rhs[3] + Bz_qu * (p_tensor_rhs[0] + p_tensor_rhs[3]) + 4.0 * By_cu * Bz * p_tensor_rhs[4] - 2.0 * By * Bz_cu * p_tensor_rhs[4] + Bx_qu * (p_tensor_rhs[3] +
    p_tensor_rhs[5]) + By_sq * Bz_sq * (p_tensor_rhs[0] + 3.0 * p_tensor_rhs[5]) + Bx_sq * (-2.0 * By * Bz * p_tensor_rhs[4] + By_sq * (3.0 * p_tensor_rhs[0] +
    p_tensor_rhs[5]) + Bz_sq * (p_tensor_rhs[0] + 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_qu) / denom - p_tensor_old[3];

  p_tensor_new[4] = 2.0 * (p_tensor_rhs[4] + dt1 * (By * p_tensor_rhs[1] - Bz * p_tensor_rhs[2] + Bx * (-p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m + dt1_sq *
    (3.0 * Bx * Bz * p_tensor_rhs[1] + Bx_sq * p_tensor_rhs[4] + 4.0 * By_sq * p_tensor_rhs[4] + 4.0 * Bz_sq * p_tensor_rhs[4] + By * (3.0 * Bx * p_tensor_rhs[2] +
    Bz * (-2.0 * p_tensor_rhs[0] + p_tensor_rhs[3] + p_tensor_rhs[5]))) * q_over_m_sq + dt1_cu * (4.0 * By_cu * p_tensor_rhs[1] - 2.0 * By * Bz_sq * p_tensor_rhs[1] +
    2.0 * By_sq * Bz * p_tensor_rhs[2] - 4.0 * Bz_cu * p_tensor_rhs[2] + Bx_sq * (-2.0 * By * p_tensor_rhs[1] + 2.0 * Bz * p_tensor_rhs[2]) + Bx_cu * (-p_tensor_rhs[3] +
    p_tensor_rhs[5]) + Bx * (-(Bz_sq * (3.0 * p_tensor_rhs[0] + p_tensor_rhs[3] - 4.0 * p_tensor_rhs[5])) + By_sq * (3.0 * p_tensor_rhs[0] - 4.0 * p_tensor_rhs[3] +
    p_tensor_rhs[5]))) * q_over_m_cu - 2.0 * By * Bz * dt1_qu * (-6.0 * Bx * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) - 6.0 * By * Bz * p_tensor_rhs[4] + Bz_sq *
    (p_tensor_rhs[0] + p_tensor_rhs[3] - 2.0 * p_tensor_rhs[5]) + By_sq * (p_tensor_rhs[0] - 2.0 * p_tensor_rhs[3] + p_tensor_rhs[5]) + Bx_sq * (-2.0 * p_tensor_rhs[0] +
    p_tensor_rhs[3] + p_tensor_rhs[5])) * q_over_m_qu) / denom - p_tensor_old[4];

  p_tensor_new[5] = 2.0 * (p_tensor_rhs[5] + 2.0 * dt1 * (By * p_tensor_rhs[2] - Bx * p_tensor_rhs[4]) * q_over_m + dt1_sq * (2.0 * Bx * Bz * p_tensor_rhs[2] + By *
    (-4.0 * Bx * p_tensor_rhs[1] + 2 * Bz * p_tensor_rhs[4]) + 5.0 * Bz_sq * p_tensor_rhs[5] + By_sq * (2.0 * p_tensor_rhs[0] + 3.0 * p_tensor_rhs[5]) + Bx_sq * 
    (2.0 * p_tensor_rhs[3] + 3.0 * p_tensor_rhs[5])) * q_over_m_sq - 2.0 * dt1_cu * (Bx_sq * (3.0 * Bz * p_tensor_rhs[1] - By * p_tensor_rhs[2]) - By * (3.0 * By * Bz *
    p_tensor_rhs[1] + By_sq * p_tensor_rhs[2] + 4.0 * Bz_sq * p_tensor_rhs[2]) + Bx_cu * p_tensor_rhs[4] + Bx * (3.0 * By * Bz * (-p_tensor_rhs[0] + p_tensor_rhs[3]) +
    By_sq * p_tensor_rhs[4] + 4.0 * Bz_sq * p_tensor_rhs[4])) * q_over_m_cu + 2.0 * dt1_qu * (-2.0 * Bx_cu * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) - 2.0 * Bx *
    (By_sq - 2.0 * Bz_sq) * (By * p_tensor_rhs[1] + Bz * p_tensor_rhs[2]) + By_sq * Bz_sq * (p_tensor_rhs[0] + 3.0 * p_tensor_rhs[3]) - 2.0 * By_cu * Bz *
    p_tensor_rhs[4] + 4.0 * By * Bz_cu * p_tensor_rhs[4] + 2.0 * Bz_qu * p_tensor_rhs[5] + By_qu * (p_tensor_rhs[0] + p_tensor_rhs[5]) + Bx_qu * (p_tensor_rhs[3] +
    p_tensor_rhs[5]) + Bx_sq * (Bz_sq * (3.0 * p_tensor_rhs[0] + p_tensor_rhs[3]) - 2.0 * By * Bz * p_tensor_rhs[4] + By_sq * (p_tensor_rhs[0] + p_tensor_rhs[3] + 2.0 *
    p_tensor_rhs[5]))) * q_over_m_qu) / denom - p_tensor_old[5];
}