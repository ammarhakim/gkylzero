#include <gkyl_efit_priv.h>
#include <float.h>

static void print_result(int n, double x[], double dx[], double errx, double errf, int niter)
{
  double x0 = x[0];
  double y0 = x[1];
  if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 ) {
    printf("x = ");
    for(int i=0; i<n; i++) printf(" %g ",  x[i]);
    printf("\n");
    printf("dx = ");
    for(int i=0; i<n; i++) printf("= %g ",  dx[i]);
    printf("\n");
    printf("errx = %g\n", errx);
    printf("errf = %g\n", errf);
  }

}

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol, bool cubics)
{
  int n = 2;
  double x[2] = {0.0,0.0};
  double dx[2] = {0.0,0.0};
  double fjac[2][2];
  double fjac_inv[2][2];
  double fvec[2];
  double p[2];
  int ntrial=100;
  double errx = 0.0;
  double errf = 0.0;
  for(int niter = 0; niter < ntrial; niter++) {
    if (cubics) {
      for(int i=0; i<n; i++) fvec[i] = up->rzbasis_cubic.eval_grad_expand(i,x,coeffs);
      fjac[0][0] = up->rzbasis_cubic.eval_laplacian_expand(0,x,coeffs);
      fjac[0][1] = up->rzbasis_cubic.eval_mixedpartial_expand(x,coeffs);
      fjac[1][0] = up->rzbasis_cubic.eval_mixedpartial_expand(x,coeffs);
      fjac[1][1] = up->rzbasis_cubic.eval_laplacian_expand(1,x,coeffs);
    }
    else {
      for(int i=0; i<n; i++) fvec[i] = up->rzbasis.eval_grad_expand(i,x,coeffs);
      fjac[0][0] = up->rzbasis.eval_laplacian_expand(0,x,coeffs);
      fjac[0][1] = up->rzbasis.eval_mixedpartial_expand(x,coeffs);
      fjac[1][0] = up->rzbasis.eval_mixedpartial_expand(x,coeffs);
      fjac[1][1] = up->rzbasis.eval_laplacian_expand(1,x,coeffs);
    }
    errf = 0.0;
    for (int i=0;i<n;i++) errf += fvec[i]*fvec[i];
    errf = sqrt(errf);
    if (errf <= 1e-19) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }
    for (int i=0;i<n;i++) p[i] = -fvec[i];

    double det = fjac[0][0]*fjac[1][1] - fjac[0][1]*fjac[1][0];
    fjac_inv[0][0] = fjac[1][1]/det;
    fjac_inv[1][1] = fjac[0][0]/det;
    fjac_inv[0][1] = -fjac[0][1]/det;
    fjac_inv[1][0] = -fjac[1][0]/det;

    dx[0] = fjac_inv[0][0]*p[0] + fjac_inv[0][1]*p[1];
    dx[1] = fjac_inv[1][0]*p[0] + fjac_inv[1][1]*p[1];

    errx = 0.0;
    for (int i=0;i<n;i++) {
      errx += fabs(dx[i]);
      x[i] += dx[i];
    }

    if(errx<=1e-18) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }

  }
  return false;

}

void
find_xpts(gkyl_efit* up)
{
    bool found_xpt = false;
    double Rsep, Zsep;
    double psisep = DBL_MAX;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &up->rzlocal);
    while (gkyl_range_iter_next(&iter)) {
      if ((iter.idx[1] < gkyl_range_shape(&up->rzlocal,1)/2 + 1) || (!up->reflect)) {
        const double* psi = gkyl_array_cfetch(up->psizr, gkyl_range_idx(&up->rzlocal, iter.idx));
        double xsol[2];
        bool status = newton_raphson(up, psi, xsol, false);
        double x0 = xsol[0];
        double y0 = xsol[1];
        double psi0 = up->rzbasis.eval_expand(xsol, psi);
        if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 && status) {
          found_xpt = true;
          double xc[2];
          gkyl_rect_grid_cell_center(&up->rzgrid, iter.idx, xc);
          double R0 = up->rzgrid.dx[0]*x0/2.0 + xc[0];
          double Z0 = up->rzgrid.dx[1]*y0/2.0 + xc[1];
          if(fabs(psi0 - up->sibry) <= fabs(psisep - up->sibry)) {
              Rsep = R0;
              Zsep = Z0;
              psisep = psi0;
          }
        }
      }
    }

    if (found_xpt) {
      if (up->reflect) {
        up->num_xpts = 2;
        up->Rxpt = gkyl_malloc(sizeof(double)*up->num_xpts);
        up->Zxpt = gkyl_malloc(sizeof(double)*up->num_xpts);
        up->Rxpt[0] = Rsep;
        up->Rxpt[1] = Rsep;
        up->Zxpt[0] = Zsep;
        up->Zxpt[1] = -Zsep;
        up->psisep = psisep;
      }
      else {
        up->num_xpts = 1;
        up->Rxpt = gkyl_malloc(sizeof(double)*up->num_xpts);
        up->Zxpt = gkyl_malloc(sizeof(double)*up->num_xpts);
        up->Rxpt[0] = Rsep;
        up->Zxpt[0] = Zsep;
        up->psisep = psisep;
      }
    }
}

void
find_xpts_cubic(gkyl_efit* up)
{
    bool found_xpt = false;
    double Rsep, Zsep;
    double psisep = DBL_MAX;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &up->rzlocal_cubic);
    while (gkyl_range_iter_next(&iter)) {
      if ((iter.idx[1] < gkyl_range_shape(&up->rzlocal_cubic,1)/2 + 1) || (!up->reflect)) {
        const double* psi = gkyl_array_cfetch(up->psizr_cubic, gkyl_range_idx(&up->rzlocal_cubic, iter.idx));
        double xsol[2];
        bool status = newton_raphson(up, psi, xsol, true);
        double x0 = xsol[0];
        double y0 = xsol[1];
        double psi0 = up->rzbasis_cubic.eval_expand(xsol, psi);
        if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 && status) {
          found_xpt = true;
          double xc[2];
          gkyl_rect_grid_cell_center(&up->rzgrid_cubic, iter.idx, xc);
          double R0 = up->rzgrid_cubic.dx[0]*x0/2.0 + xc[0];
          double Z0 = up->rzgrid_cubic.dx[1]*y0/2.0 + xc[1];
          if(fabs(psi0 - up->sibry) <= fabs(psisep - up->sibry)) {
              Rsep = R0;
              Zsep = Z0;
              psisep = psi0;
          }
        }
      }
    }

    if (found_xpt) {
      if (up->reflect) {
        up->num_xpts_cubic = 2;
        up->Rxpt_cubic = gkyl_malloc(sizeof(double)*up->num_xpts_cubic);
        up->Zxpt_cubic = gkyl_malloc(sizeof(double)*up->num_xpts_cubic);
        up->Rxpt_cubic[0] = Rsep;
        up->Rxpt_cubic[1] = Rsep;
        up->Zxpt_cubic[0] = Zsep;
        up->Zxpt_cubic[1] = -Zsep;
        up->psisep_cubic = psisep;
      }
      else {
        up->num_xpts_cubic = 1;
        up->Rxpt_cubic = gkyl_malloc(sizeof(double)*up->num_xpts_cubic);
        up->Zxpt_cubic = gkyl_malloc(sizeof(double)*up->num_xpts_cubic);
        up->Rxpt_cubic[0] = Rsep;
        up->Zxpt_cubic[0] = Zsep;
        up->psisep_cubic = psisep;
      }
    }
}

