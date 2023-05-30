#include <rt_euler_riem_2d.h>

int
main(int argc, char **argv)
{
  return rt_euler_riem_2d_run(argc, argv, WV_EULER_RP_HLL, GKYL_MOMENT_WAVE_PROP);
}
