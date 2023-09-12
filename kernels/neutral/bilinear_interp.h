#include <gkyl_array.h>

void bilinear_interp(double* x, double* y, double* z,
  const struct gkyl_array* xq, const struct gkyl_array* yq,
  int nx, int ny, int nq, struct gkyl_array* zq);
