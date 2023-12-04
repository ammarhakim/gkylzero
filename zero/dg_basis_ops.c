#include <gkyl_dg_basis_ops.h>
#include <gkyl_array.h>
#include <gkyl_alloc.h>

enum dg_basis_op_code {
  GKYL_DG_BASIS_OP_CUBIC_1D,
  GKYL_DG_BASIS_OP_CUBIC_2D
};

struct gkyl_dg_basis_op_mem {
  enum dg_basis_op_code opcode;
  union {
    // data for GKYL_DG_BASIS_OP_CUBIC_1D
    struct {
      struct gkyl_array *grad1dx;
    };
    // data for GKYL_DG_BASIS_OP_CUBIC_2D
    struct {
      struct gkyl_array *grad2dx, *grad2dy, *grad2dxy;
    };
  };
};

void
gkyl_dg_calc_cubic_1d(const double val[2], const double grad[2], double *coeff)
{
  coeff[0] = 0.7071067811865475*val[1]-0.2357022603955158*grad[1]+0.7071067811865475*val[0]+0.2357022603955158*grad[0]; 
  coeff[1] = 0.4898979485566357*val[1]-0.08164965809277261*grad[1]-0.4898979485566357*val[0]-0.08164965809277261*grad[0]; 
  coeff[2] = 0.105409255338946*grad[1]-0.105409255338946*grad[0]; 
  coeff[3] = (-0.05345224838248487*val[1])+0.05345224838248487*grad[1]+0.05345224838248487*val[0]+0.05345224838248487*grad[0];
}

void
gkyl_dg_calc_cubic_2d(const double f[4],
  const double fx[4], const double fy[4], const double fxy[4],
  double *coeff)
{
  coeff[0] = (-0.1666666666666667*fy[3])+0.05555555555555555*fxy[3]-0.1666666666666667*fx[3]+0.5*f[3]+0.1666666666666667*fy[2]-0.05555555555555555*fxy[2]-0.1666666666666667*fx[2]+0.5*f[2]-0.1666666666666667*fy[1]-0.05555555555555555*fxy[1]+0.1666666666666667*fx[1]+0.5*f[1]+0.1666666666666667*fy[0]+0.05555555555555555*fxy[0]+0.1666666666666667*fx[0]+0.5*f[0]; 
  coeff[1] = (-0.1154700538379252*fy[3])+0.01924500897298753*fxy[3]-0.05773502691896259*fx[3]+0.3464101615137755*f[3]+0.1154700538379252*fy[2]-0.01924500897298753*fxy[2]-0.05773502691896259*fx[2]+0.3464101615137755*f[2]+0.1154700538379252*fy[1]+0.01924500897298753*fxy[1]-0.05773502691896259*fx[1]-0.3464101615137755*f[1]-0.1154700538379252*fy[0]-0.01924500897298753*fxy[0]-0.05773502691896259*fx[0]-0.3464101615137755*f[0]; 
  coeff[2] = (-0.05773502691896259*fy[3])+0.01924500897298753*fxy[3]-0.1154700538379252*fx[3]+0.3464101615137755*f[3]-0.05773502691896259*fy[2]+0.01924500897298753*fxy[2]+0.1154700538379252*fx[2]-0.3464101615137755*f[2]-0.05773502691896259*fy[1]-0.01924500897298753*fxy[1]+0.1154700538379252*fx[1]+0.3464101615137755*f[1]-0.05773502691896259*fy[0]-0.01924500897298753*fxy[0]-0.1154700538379252*fx[0]-0.3464101615137755*f[0]; 
  coeff[3] = (-0.04*fy[3])+0.006666666666666667*fxy[3]-0.04*fx[3]+0.24*f[3]-0.04*fy[2]+0.006666666666666667*fxy[2]+0.04*fx[2]-0.24*f[2]+0.04*fy[1]+0.006666666666666667*fxy[1]-0.04*fx[1]-0.24*f[1]+0.04*fy[0]+0.006666666666666667*fxy[0]+0.04*fx[0]+0.24*f[0]; 
  coeff[4] = (-0.02484519974999766*fxy[3])+0.07453559924999298*fx[3]+0.02484519974999766*fxy[2]+0.07453559924999298*fx[2]+0.02484519974999766*fxy[1]-0.07453559924999298*fx[1]-0.02484519974999766*fxy[0]-0.07453559924999298*fx[0]; 
  coeff[5] = 0.07453559924999298*fy[3]-0.02484519974999766*fxy[3]-0.07453559924999298*fy[2]+0.02484519974999766*fxy[2]+0.07453559924999298*fy[1]+0.02484519974999766*fxy[1]-0.07453559924999298*fy[0]-0.02484519974999766*fxy[0]; 
  coeff[6] = (-0.008606629658238702*fxy[3])+0.05163977794943223*fx[3]-0.008606629658238702*fxy[2]-0.05163977794943223*fx[2]+0.008606629658238702*fxy[1]-0.05163977794943223*fx[1]+0.008606629658238702*fxy[0]+0.05163977794943223*fx[0]; 
  coeff[7] = 0.05163977794943223*fy[3]-0.008606629658238702*fxy[3]-0.05163977794943223*fy[2]+0.008606629658238702*fxy[2]-0.05163977794943223*fy[1]-0.008606629658238702*fxy[1]+0.05163977794943223*fy[0]+0.008606629658238702*fxy[0]; 
  coeff[8] = 0.01259881576697424*fy[3]-0.01259881576697424*fxy[3]+0.03779644730092272*fx[3]-0.03779644730092272*f[3]-0.01259881576697424*fy[2]+0.01259881576697424*fxy[2]+0.03779644730092272*fx[2]-0.03779644730092272*f[2]-0.01259881576697424*fy[1]-0.01259881576697424*fxy[1]+0.03779644730092272*fx[1]+0.03779644730092272*f[1]+0.01259881576697424*fy[0]+0.01259881576697424*fxy[0]+0.03779644730092272*fx[0]+0.03779644730092272*f[0]; 
  coeff[9] = 0.03779644730092272*fy[3]-0.01259881576697424*fxy[3]+0.01259881576697424*fx[3]-0.03779644730092272*f[3]+0.03779644730092272*fy[2]-0.01259881576697424*fxy[2]-0.01259881576697424*fx[2]+0.03779644730092272*f[2]+0.03779644730092272*fy[1]+0.01259881576697424*fxy[1]-0.01259881576697424*fx[1]-0.03779644730092272*f[1]+0.03779644730092272*fy[0]+0.01259881576697424*fxy[0]+0.01259881576697424*fx[0]+0.03779644730092272*f[0]; 
  coeff[10] = 0.01111111111111111*fxy[3]-0.01111111111111111*fxy[2]-0.01111111111111111*fxy[1]+0.01111111111111111*fxy[0]; 
  coeff[11] = 0.004364357804719848*fy[3]-0.004364357804719848*fxy[3]+0.02618614682831908*fx[3]-0.02618614682831908*f[3]+0.004364357804719848*fy[2]-0.004364357804719848*fxy[2]-0.02618614682831908*fx[2]+0.02618614682831908*f[2]-0.004364357804719848*fy[1]-0.004364357804719848*fxy[1]+0.02618614682831908*fx[1]+0.02618614682831908*f[1]-0.004364357804719848*fy[0]-0.004364357804719848*fxy[0]-0.02618614682831908*fx[0]-0.02618614682831908*f[0]; 
  coeff[12] = 0.02618614682831908*fy[3]-0.004364357804719848*fxy[3]+0.004364357804719848*fx[3]-0.02618614682831908*f[3]+0.02618614682831908*fy[2]-0.004364357804719848*fxy[2]-0.004364357804719848*fx[2]+0.02618614682831908*f[2]-0.02618614682831908*fy[1]-0.004364357804719848*fxy[1]+0.004364357804719848*fx[1]+0.02618614682831908*f[1]-0.02618614682831908*fy[0]-0.004364357804719848*fxy[0]-0.004364357804719848*fx[0]-0.02618614682831908*f[0]; 
  coeff[13] = (-0.00563436169819011*fy[3])+0.00563436169819011*fxy[3]+0.00563436169819011*fy[2]-0.00563436169819011*fxy[2]+0.00563436169819011*fy[1]+0.00563436169819011*fxy[1]-0.00563436169819011*fy[0]-0.00563436169819011*fxy[0]; 
  coeff[14] = 0.00563436169819011*fxy[3]-0.00563436169819011*fx[3]+0.00563436169819011*fxy[2]+0.00563436169819011*fx[2]-0.00563436169819011*fxy[1]+0.00563436169819011*fx[1]-0.00563436169819011*fxy[0]-0.00563436169819011*fx[0]; 
  coeff[15] = (-0.002857142857142857*fy[3])+0.002857142857142857*fxy[3]-0.002857142857142857*fx[3]+0.002857142857142857*f[3]-0.002857142857142857*fy[2]+0.002857142857142857*fxy[2]+0.002857142857142857*fx[2]-0.002857142857142857*f[2]+0.002857142857142857*fy[1]+0.002857142857142857*fxy[1]-0.002857142857142857*fx[1]-0.002857142857142857*f[1]+0.002857142857142857*fy[0]+0.002857142857142857*fxy[0]+0.002857142857142857*fx[0]+0.002857142857142857*f[0];
}

gkyl_dg_basis_op_mem *
gkyl_dg_alloc_cubic_1d(int cells)
{
  struct gkyl_dg_basis_op_mem *mem = gkyl_malloc(sizeof(*mem));
  mem->opcode = GKYL_DG_BASIS_OP_CUBIC_1D;
  mem->grad1dx = gkyl_array_new(GKYL_DOUBLE, 1, cells+1);
  return mem;
}

gkyl_dg_basis_op_mem *
gkyl_dg_alloc_cubic_2d(int cells[2])
{
  struct gkyl_dg_basis_op_mem *mem = gkyl_malloc(sizeof(*mem));
  mem->opcode = GKYL_DG_BASIS_OP_CUBIC_2D;

  size_t ncells = (cells[0]+1)*(cells[1]+1);
  mem->grad2dx = gkyl_array_new(GKYL_DOUBLE, 1, ncells);
  mem->grad2dy = gkyl_array_new(GKYL_DOUBLE, 1, ncells);
  mem->grad2dxy = gkyl_array_new(GKYL_DOUBLE, 1, ncells);
  
  return mem;
}

void
gkyl_dg_basis_op_mem_release(gkyl_dg_basis_op_mem *mem)
{
  if (mem->opcode == GKYL_DG_BASIS_OP_CUBIC_1D) {
    gkyl_array_release(mem->grad1dx);
  }
  if (mem->opcode == GKYL_DG_BASIS_OP_CUBIC_2D) {
    gkyl_array_release(mem->grad2dx);
    gkyl_array_release(mem->grad2dy);
    gkyl_array_release(mem->grad2dxy);
  }  
  gkyl_free(mem);
}

void
gkyl_dg_calc_cubic_1d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells, double dx,
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic)
{
  enum { I, LL, L, R, RR, XE }; // i, i-2, i-1, i+1, i+2 nodes
  
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, (int[1]) { cells });

  struct gkyl_range nc_range;
  gkyl_range_init_from_shape(&nc_range, 1, (int[1]) { cells+1 });

  long offset[XE];
  offset[I] = 0; // i
  offset[LL] = gkyl_range_offset(&nc_range, (int[]) { -2 } ); // i-1  
  offset[L] = gkyl_range_offset(&nc_range, (int[]) { -1 } ); // i-1
  offset[R] = gkyl_range_offset(&nc_range, (int[]) { 1 } ); // i+1
  offset[RR] = gkyl_range_offset(&nc_range, (int[]) { 2 } ); // i+2

  struct gkyl_array *gradx = mem->grad1dx;

  int ilo = nc_range.lower[0], iup = nc_range.upper[0];

  // Step 1: compute gradients at nodes using differencing
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_range); // loop over node range
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_range, iter.idx);

    if ((iter.idx[0] != ilo) && (iter.idx[0] != iup)) {
      // interior nodes
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      
      double *grad_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      grad_I[0] = (val_R[0]-val_L[0])/(2*dx);
    }

    if (iter.idx[0] == ilo) {
      // left boundary: use second-order one sided differencing
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      const double *val_RR = gkyl_array_cfetch(nodal_vals, nidx+offset[RR]);
      
      double *grad_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      grad_I[0] = -(val_RR[0]-4*val_R[0]+3*val_I[0])/(2*dx);
    }

    if (iter.idx[0] == iup) {
      // right boundary: use second-order one sided differencing
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_LL = gkyl_array_cfetch(nodal_vals, nidx+offset[LL]);
      
      double *grad_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      grad_I[0] = (3*val_I[0]-4*val_L[0]+val_LL[0])/(2*dx);
    }
  }

  // Step 2: compute cubic expansions in each cell
  gkyl_range_iter_init(&iter, &range); // loop is over cells
  while (gkyl_range_iter_next(&iter)) {
    
    long nidx = gkyl_range_idx(&nc_range, iter.idx);

    const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
    const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
    const double *grad_I = gkyl_array_cfetch(gradx, nidx+offset[I]);
    const double *grad_R = gkyl_array_cfetch(gradx, nidx+offset[R]);

    double val[2] = { val_I[0], val_R[0] };
    double grad[2] = { grad_I[0]*dx/2, grad_R[0]*dx/2 };

    long cidx = gkyl_range_idx(&range, iter.idx);    
    gkyl_dg_calc_cubic_1d(val, grad, gkyl_array_fetch(cubic, cidx));
  }
}
