#include <gkyl_efit.h>

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol);

void
find_xpts(gkyl_efit* up);

