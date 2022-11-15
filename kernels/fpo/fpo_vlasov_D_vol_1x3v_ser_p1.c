#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_D_vol_1x3v_ser_p1(const double* w, const double* dx,
  const double* g, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // g: Rosenbluth potential
  // out: Incremented output

  const double Jvxvx = 4/dx[1]/dx[1];

  out[0 + 0] +=  Jvxvx*(0.0);
  out[0 + 1] +=  Jvxvx*(0.0);
  out[0 + 2] +=  Jvxvx*(0.0);
  out[0 + 3] +=  Jvxvx*(0.0);
  out[0 + 4] +=  Jvxvx*(0.0);
  out[0 + 5] +=  Jvxvx*(0.0);
  out[0 + 6] +=  Jvxvx*(0.0);
  out[0 + 7] +=  Jvxvx*(0.0);
  out[0 + 8] +=  Jvxvx*(0.0);
  out[0 + 9] +=  Jvxvx*(0.0);
  out[0 + 10] +=  Jvxvx*(0.0);
  out[0 + 11] +=  Jvxvx*(0.0);
  out[0 + 12] +=  Jvxvx*(0.0);
  out[0 + 13] +=  Jvxvx*(0.0);
  out[0 + 14] +=  Jvxvx*(0.0);
  out[0 + 15] +=  Jvxvx*(0.0);

  const double Jvxvy = 4/dx[1]/dx[2];

  out[16 + 0] +=  Jvxvy*(0.0);
  out[16 + 1] +=  Jvxvy*(0.0);
  out[16 + 2] +=  Jvxvy*(0.0);
  out[16 + 3] +=  Jvxvy*(0.0);
  out[16 + 4] +=  Jvxvy*(0.0);
  out[16 + 5] +=  Jvxvy*(0.0);
  out[16 + 6] +=  Jvxvy*(0.0);
  out[16 + 7] +=  Jvxvy*(3.0*g[0]);
  out[16 + 8] +=  Jvxvy*(0.0);
  out[16 + 9] +=  Jvxvy*(0.0);
  out[16 + 10] +=  Jvxvy*(0.0);
  out[16 + 11] +=  Jvxvy*(3.0*g[1]);
  out[16 + 12] +=  Jvxvy*(0.0);
  out[16 + 13] +=  Jvxvy*(0.0);
  out[16 + 14] +=  Jvxvy*(3.0*g[4]);
  out[16 + 15] +=  Jvxvy*(3.0*g[8]);

  const double Jvxvz = 4/dx[1]/dx[3];

  out[32 + 0] +=  Jvxvz*(0.0);
  out[32 + 1] +=  Jvxvz*(0.0);
  out[32 + 2] +=  Jvxvz*(0.0);
  out[32 + 3] +=  Jvxvz*(0.0);
  out[32 + 4] +=  Jvxvz*(0.0);
  out[32 + 5] +=  Jvxvz*(0.0);
  out[32 + 6] +=  Jvxvz*(0.0);
  out[32 + 7] +=  Jvxvz*(0.0);
  out[32 + 8] +=  Jvxvz*(0.0);
  out[32 + 9] +=  Jvxvz*(3.0*g[0]);
  out[32 + 10] +=  Jvxvz*(0.0);
  out[32 + 11] +=  Jvxvz*(0.0);
  out[32 + 12] +=  Jvxvz*(3.0*g[1]);
  out[32 + 13] +=  Jvxvz*(0.0);
  out[32 + 14] +=  Jvxvz*(3.0*g[3]);
  out[32 + 15] +=  Jvxvz*(3.0*g[6]);

  const double Jvyvy = 4/dx[2]/dx[2];

  out[48 + 0] +=  Jvyvy*(0.0);
  out[48 + 1] +=  Jvyvy*(0.0);
  out[48 + 2] +=  Jvyvy*(0.0);
  out[48 + 3] +=  Jvyvy*(0.0);
  out[48 + 4] +=  Jvyvy*(0.0);
  out[48 + 5] +=  Jvyvy*(0.0);
  out[48 + 6] +=  Jvyvy*(0.0);
  out[48 + 7] +=  Jvyvy*(0.0);
  out[48 + 8] +=  Jvyvy*(0.0);
  out[48 + 9] +=  Jvyvy*(0.0);
  out[48 + 10] +=  Jvyvy*(0.0);
  out[48 + 11] +=  Jvyvy*(0.0);
  out[48 + 12] +=  Jvyvy*(0.0);
  out[48 + 13] +=  Jvyvy*(0.0);
  out[48 + 14] +=  Jvyvy*(0.0);
  out[48 + 15] +=  Jvyvy*(0.0);

  const double Jvyvz = 4/dx[2]/dx[3];

  out[64 + 0] +=  Jvyvz*(0.0);
  out[64 + 1] +=  Jvyvz*(0.0);
  out[64 + 2] +=  Jvyvz*(0.0);
  out[64 + 3] +=  Jvyvz*(0.0);
  out[64 + 4] +=  Jvyvz*(0.0);
  out[64 + 5] +=  Jvyvz*(0.0);
  out[64 + 6] +=  Jvyvz*(0.0);
  out[64 + 7] +=  Jvyvz*(0.0);
  out[64 + 8] +=  Jvyvz*(0.0);
  out[64 + 9] +=  Jvyvz*(0.0);
  out[64 + 10] +=  Jvyvz*(3.0*g[0]);
  out[64 + 11] +=  Jvyvz*(0.0);
  out[64 + 12] +=  Jvyvz*(0.0);
  out[64 + 13] +=  Jvyvz*(3.0*g[1]);
  out[64 + 14] +=  Jvyvz*(3.0*g[2]);
  out[64 + 15] +=  Jvyvz*(3.0*g[5]);

  const double Jvzvz = 4/dx[3]/dx[3];

  out[80 + 0] +=  Jvzvz*(0.0);
  out[80 + 1] +=  Jvzvz*(0.0);
  out[80 + 2] +=  Jvzvz*(0.0);
  out[80 + 3] +=  Jvzvz*(0.0);
  out[80 + 4] +=  Jvzvz*(0.0);
  out[80 + 5] +=  Jvzvz*(0.0);
  out[80 + 6] +=  Jvzvz*(0.0);
  out[80 + 7] +=  Jvzvz*(0.0);
  out[80 + 8] +=  Jvzvz*(0.0);
  out[80 + 9] +=  Jvzvz*(0.0);
  out[80 + 10] +=  Jvzvz*(0.0);
  out[80 + 11] +=  Jvzvz*(0.0);
  out[80 + 12] +=  Jvzvz*(0.0);
  out[80 + 13] +=  Jvzvz*(0.0);
  out[80 + 14] +=  Jvzvz*(0.0);
  out[80 + 15] +=  Jvzvz*(0.0);

}
