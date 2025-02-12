#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gkgeom.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

static void
calc_exact_quad_1d(double dx, double xc, double coeff[4])
{
/* Exact expansion computed using the following Maxima code:

   load("basis-precalc/basisTensor1x")$
   load("modal-basis")$
   bc : basisC[3]$

   f : -(x1-2.5)^2 + 0.15*x1 + 25$
   fS : subst([x1=xc+dx/2*x],f)$
   proj : calcInnerProdList(varsC,1,bc,fS)$
 */
  coeff[0] = (-1.414213562373095*sq(xc))+7.283199846221438*xc-0.1178511301977579*sq(dx)+26.51650429449552; 
  coeff[1] = 2.102478695888894*dx-0.8164965809277258*dx*xc; 
  coeff[2] = -0.105409255338946*sq(dx); 
  coeff[3] = 0.0; 
}

static void
calc_exact_xquad_2d(double dx, double dy, double xc, double yc, double coeff[16])
{
/* Exact expansion computed using the following Maxima code:

   load("basis-precalc/basisTensor2x")$
   load("modal-basis")$
   bc : basisC[3]$

   f : -(x1-2.5)^2 + 0.15*x1 + 25$   
   fS : subst([x1=xc+dx/2*x, y1=yc+dy/2*y],f)$
   proj : calcInnerProdList(varsC,1,bc,fS)$
 */

coeff[0] = (-2.0*sq(xc))+10.3*xc-0.1666666666666667*sq(dx)+37.5; 
coeff[1] = 2.973353886326573*dx-1.154700538379252*dx*xc; 
coeff[2] = 0.0; 
coeff[3] = 0.0; 
coeff[4] = -0.149071198499986*sq(dx); 
coeff[5] = 0.0; 
coeff[6] = 0.0; 
coeff[7] = 0.0; 
coeff[8] = 0.0; 
coeff[9] = 0.0; 
coeff[10] = 0.0; 
coeff[11] = 0.0; 
coeff[12] = 0.0; 
coeff[13] = 0.0; 
coeff[14] = 0.0; 
coeff[15] = 0.0;
}

static void
calc_exact_yquad_2d(double dx, double dy, double xc, double yc, double coeff[16])
{
/* Exact expansion computed using the following Maxima code:

   load("basis-precalc/basisTensor2x")$
   load("modal-basis")$
   bc : basisC[3]$

   f : -(y1-2.5)^2 + 0.15*y1 + 25$      
   fS : subst([x1=xc+dx/2*x, y1=yc+dy/2*y],f)$
   proj : calcInnerProdList(varsC,1,bc,fS)$
 */

coeff[0] = (-2.0*sq(yc))+10.3*yc-0.1666666666666667*sq(dy)+37.5; 
coeff[1] = 0.0; 
coeff[2] = 2.973353886326573*dy-1.154700538379252*dy*yc; 
coeff[3] = 0.0; 
coeff[4] = 0.0; 
coeff[5] = -0.149071198499986*sq(dy); 
coeff[6] = 0.0; 
coeff[7] = 0.0; 
coeff[8] = 0.0; 
coeff[9] = 0.0; 
coeff[10] = 0.0; 
coeff[11] = 0.0; 
coeff[12] = 0.0; 
coeff[13] = 0.0; 
coeff[14] = 0.0; 
coeff[15] = 0.0; 
}

static void
cubic_1d(void)
{
  double lower[] = { 0.0 }, upper[] = { 10.0 };
  int cells[] = { 16 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0] };
  int nc_cells[] = { cells[0] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, cells[0]+1);
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_1d(cells[0]);
  double xn[1];
   
  do {
    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);

      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = -sq(xn[0]) + 0.15*cub(xn[0]);
    }
    
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_1d_a.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_1d_a.gkyl");
  } while (0);

  do {
    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);

      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = sq(xn[0]) + 10*sq(sin(xn[0]));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_1d_b.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_1d_b.gkyl");
  } while (0);

  do {
    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = -sq(xn[0]-2.5) + 0.15*xn[0] + 25.0;
    }
    
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);

    double err = 0.0;

    // check against exact expansion
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&local, iter.idx);
      
      double xc[1], coeff[4];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      calc_exact_quad_1d(grid.dx[0], xc[0], coeff);
        
      const double *pc = gkyl_array_cfetch(psi_cubic, nidx);
      for (int i=0; i<4; ++i)
        err = fmax(err, fabs(pc[i]-coeff[i]));
    }
    printf("Max error for a 1D quadratic is %lg\n", err);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_1d_c.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_1d_c.gkyl");
  } while (0);

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_dg_basis_op_mem_release(mem);
}

static void
cubic_2d(void)
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 1.0, 1.0 };
  int cells[] = { 8, 8 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);

  double xn[2];
  
  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = sin(2*M_PI*xn[0]);
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_2d_a1.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_2d_a1.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = sin(2*M_PI*xn[1]);
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_2d_a2.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_2d_a2.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = (-sq(xn[1]) + 0.15*cub(xn[1]));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_2d_b.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_2d_b.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);

      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = -sq(xn[0]-2.5) + 0.15*xn[0] + 25.0;
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);

    int ilo = local.lower[0], iup = local.upper[0];
    int jlo = local.lower[1], jup = local.upper[1];

    double err = 0.0;

    // check against exact expansion
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&local, iter.idx);
      
      double xc[2], coeff[16];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      calc_exact_xquad_2d(grid.dx[0], grid.dx[1], xc[0], xc[1], coeff);
        
      const double *pc = gkyl_array_cfetch(psi_cubic, nidx);
      for (int i=0; i<16; ++i)
        err = fmax(err, fabs(pc[i]-coeff[i]));
    }
    printf("Max error for a quadratic in x in 2D is %lg\n", err);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_2d_c1.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_2d_c1.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);

      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = -sq(xn[1]-2.5) + 0.15*xn[1] + 25.0;
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);

    int ilo = local.lower[0], iup = local.upper[0];
    int jlo = local.lower[1], jup = local.upper[1];

    double err = 0.0;

    // check against exact expansion
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&local, iter.idx);
      
      double xc[2], coeff[16];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      calc_exact_yquad_2d(grid.dx[0], grid.dx[1], xc[0], xc[1], coeff);
        
      const double *pc = gkyl_array_cfetch(psi_cubic, nidx);
      for (int i=0; i<16; ++i)
        err = fmax(err, fabs(pc[i]-coeff[i]));
    }
    printf("Max error for a quadratic in y in 2D is %lg\n", err);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_2d_c2.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_2d_c2.gkyl");
  } while (0);  

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = sin(2*M_PI*xn[1])*cub(xn[1]);
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);

    // check continuity of value and gradients

    do {
      int idx[2] = { 2, 4 };
      long nidx = gkyl_range_idx(&local, idx);
      const double *pc = gkyl_array_cfetch(psi_cubic, nidx);

      int idxl[2] = { 1, 4 };
      nidx = gkyl_range_idx(&local, idxl);
      const double *pcl = gkyl_array_cfetch(psi_cubic, nidx);

      double psiL = basis.eval_expand((double[]) { 1.0, 0.0 }, pcl);
      double psiR = basis.eval_expand((double[]) { -1.0, 0.0 }, pc);

      printf("Error in values across cell: %lg\n", fabs(psiL-psiR));

      double dpsixL = basis.eval_grad_expand(0, (double[]) { 1.0, 0.0 }, pcl);
      double dpsixR = basis.eval_grad_expand(0, (double[]) { -1.0, 0.0 }, pc);

      printf("Error in x-gradient across cell: %lg\n", fabs(dpsixL-dpsixR));

      double dpsiyL = basis.eval_grad_expand(1, (double[]) { 1.0, 0.0 }, pcl);
      double dpsiyR = basis.eval_grad_expand(1, (double[]) { -1.0, 0.0 }, pc);

      printf("Error in y-gradient across cell: %lg\n", fabs(dpsiyL-dpsiyR));
    } while (0);
    
    do {
      int idx[2] = { 3, 4 };
      long nidx = gkyl_range_idx(&local, idx);
      const double *pc = gkyl_array_cfetch(psi_cubic, nidx);

      int idxl[2] = { 2, 4 };
      nidx = gkyl_range_idx(&local, idxl);
      const double *pcl = gkyl_array_cfetch(psi_cubic, nidx);

      double psiL = basis.eval_expand((double[]) { 1.0, 0.0 }, pcl);
      double psiR = basis.eval_expand((double[]) { -1.0, 0.0 }, pc);

      printf("Error in values across cell: %lg\n", fabs(psiL-psiR));

      double dpsixL = basis.eval_grad_expand(0, (double[]) { 1.0, 0.0 }, pcl);
      double dpsixR = basis.eval_grad_expand(0, (double[]) { -1.0, 0.0 }, pc);

      printf("Error in x-gradient across cell: %lg\n", fabs(dpsixL-dpsixR));

      double dpsiyL = basis.eval_grad_expand(1, (double[]) { 1.0, 0.0 }, pcl);
      double dpsiyR = basis.eval_grad_expand(1, (double[]) { -1.0, 0.0 }, pc);

      printf("Error in y-gradient across cell: %lg\n", fabs(dpsiyL-dpsiyR));
    } while (0);

    do {
      int idx[2] = { 3, 4 };
      long nidx = gkyl_range_idx(&local, idx);
      const double *pc = gkyl_array_cfetch(psi_cubic, nidx);

      int idxl[2] = { 3, 3 };
      nidx = gkyl_range_idx(&local, idxl);
      const double *pcl = gkyl_array_cfetch(psi_cubic, nidx);

      double psiL = basis.eval_expand((double[]) { 0.5, 1.0 }, pcl);
      double psiR = basis.eval_expand((double[]) { 0.5, -1.0 }, pc);

      printf("Error in values across cell: %lg\n", fabs(psiL-psiR));

      double dpsixL = basis.eval_grad_expand(0, (double[]) { 0.25, 1.0 }, pcl);
      double dpsixR = basis.eval_grad_expand(0, (double[]) { 0.25, -1.0 }, pc);

      printf("Error in x-gradient across cell: %lg\n", fabs(dpsixL-dpsixR));

      double dpsiyL = basis.eval_grad_expand(1, (double[]) { -0.25, 1.0 }, pcl);
      double dpsiyR = basis.eval_grad_expand(1, (double[]) { -0.25, -1.0 }, pc);

      printf("Error in y-gradient across cell: %lg\n", fabs(dpsiyL-dpsiyR));
    } while (0);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_2d_d.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_2d_d.gkyl");
  } while (0);
  
  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_dg_basis_op_mem_release(mem);
}

void
cubic_evalf_2d(void)
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 1.0, 1.0 };
  int cells[] = { 8, 8 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);

  double xn[2];
  
  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);      
      pn[0] = sin(2*M_PI*xn[0])*cub(xn[1]);
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, 0, psi_nodal, "nodal_evf.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic, "cubic_evf.gkyl");
  } while (0);

  // compute evalf function from nodal values
  struct gkyl_basis_ops_evalf *evf = gkyl_dg_basis_ops_evalf_new(&grid, psi_nodal);

  // project the cubic on cubic basis: this should result in the same
  // DG expansions
  gkyl_proj_on_basis *projCub = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .ctx = evf->ctx,
      .eval = evf->eval_cubic
    }
  );

  struct gkyl_array *psi_cubic_DG = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_proj_on_basis_advance(projCub, 0.0, &local, psi_cubic_DG);

  gkyl_grid_sub_array_write(&grid, &local, 0, psi_cubic_DG, "cubic_DG_evf.gkyl");
  
  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_array_release(psi_cubic_DG);
  
  gkyl_dg_basis_op_mem_release(mem);
  gkyl_dg_basis_ops_evalf_release(evf);
  gkyl_proj_on_basis_release(projCub);
}

int
main(int argc, char **argv)
{
  cubic_1d();
  cubic_2d();
  cubic_evalf_2d();
  
  return 0;
}
    
