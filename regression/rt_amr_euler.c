#include <gkyl_amr_core.h>

void initEuler(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];
  double y = xn[1];

  double r = sqrt(((x - 1.5) * (x - 1.5)) + ((y - 1.5) * (y - 1.5)));

  double rho = 0.0;
  double u = 0.0;
  double p = 0.0;

  if (r < 0.75)
  {
    rho = 3.0;
    u = 0.0;
    p = 3.0;
  }
  else
  {
    rho = 1.0;
    u = 0.0;
    p = 1.0;
  }

  fout[0] = rho;
  fout[1] = rho * u; fout[2] = 0.0; fout[3] = 0.0;
  fout[4] = p / (1.4 - 1.0) + 0.5 * rho * u * u;
}

int main(int argc, char **argv)
{
  euler2d_run_level1(argc, argv, initEuler, 128, 128, 2, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0, 3.0, 3.0, 0.95, 0.1);
}