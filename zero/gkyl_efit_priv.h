#include <gkyl_efit.h>

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol, bool cubics);

void
find_xpts(gkyl_efit* up);

void
find_xpts_cubic(gkyl_efit* up);

