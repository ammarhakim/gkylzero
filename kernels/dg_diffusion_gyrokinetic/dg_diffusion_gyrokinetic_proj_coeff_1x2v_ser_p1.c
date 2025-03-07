#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels.h>
#include <gkyl_basis_ser_1x_p1_inv.h> 

GKYL_CU_DH void dg_diffusion_gyrokinetic_proj_coeff_1x2v_ser_p1_diffdirsx(const double *xc, const double *dx, const double *nu, const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag, const double *vtsq, const double *vmap, const double *vmapSq, double *out) 
{ 
  // xc: cell center coordinates.
  // dx: cell length in each direction.
  // nu: overall multiplicative factor in the coefficient.
  // xi: factor multiplying v^2 in the coefficient.
  // mass: species mass.
  // vtsq_min: Minimum thermal speed supported by the grid (squared).
  // gijJ: contravariant metric multiplied by the Jacobian.
  // bmag: magnetic field amplitude .
  // vtsq: thermal speed squared.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // out: DG coefficients of the diffusion coefficent.

  double vtsq_nod[2];
  vtsq_nod[0] = 0.7071067811865475*vtsq[0]-0.7071067811865475*vtsq[1]; 
  vtsq_nod[1] = 0.7071067811865475*vtsq[1]+0.7071067811865475*vtsq[0]; 

  // If vtsq < 0. at nodes, set it to vtsq_min.
  vtsq_nod[0] = fmax(vtsq_min, vtsq_nod[0]);
  vtsq_nod[1] = fmax(vtsq_min, vtsq_nod[1]);

  double vtsq_floored[2];
  vtsq_floored[0] = 0.7071067811865475*vtsq_nod[1]+0.7071067811865475*vtsq_nod[0]; 
  vtsq_floored[1] = 0.7071067811865475*vtsq_nod[1]-0.7071067811865475*vtsq_nod[0]; 

  double vtsq_inv[2] = {0.0};
  ser_1x_p1_inv(vtsq_floored, vtsq_inv); 

  double vsq[12]; 
  vsq[0] = (2.8284271247461907*bmag[0]*vmap[2])/mass+2.0*vmapSq[0]; 
  vsq[1] = (2.8284271247461907*bmag[1]*vmap[2])/mass; 
  vsq[2] = 2.0*vmapSq[1]; 
  vsq[3] = (2.8284271247461907*bmag[0]*vmap[3])/mass; 
  vsq[5] = (2.8284271247461907*bmag[1]*vmap[3])/mass; 
  vsq[8] = 2.0*vmapSq[2]; 

  out[0] = 0.5*gijJ[0]*nu[0]*xi[0]*vsq[1]*vtsq_inv[1]+0.5*nu[0]*vsq[0]*xi[0]*gijJ[1]*vtsq_inv[1]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[1]+0.5*gijJ[0]*nu[0]*vsq[0]*vtsq_inv[0]*xi[0]+2.0*gijJ[0]*nu[0]; 
  out[1] = 0.9*nu[0]*xi[0]*gijJ[1]*vsq[1]*vtsq_inv[1]+0.5*gijJ[0]*nu[0]*vsq[0]*xi[0]*vtsq_inv[1]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[1]+0.5*nu[0]*vsq[0]*vtsq_inv[0]*xi[0]*gijJ[1]+2.0*nu[0]*gijJ[1]; 
  out[2] = 0.5*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[4]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[4]+0.5*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[2]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[2]; 
  out[3] = 0.5*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[5]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[5]+0.5*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[3]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[3]; 
  out[4] = 0.9*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[4]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[4]+0.5*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[2]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[2]; 
  out[5] = 0.9*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[5]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[5]+0.5*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[3]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[3]; 
  out[6] = 0.5*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[7]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[7]+0.5*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[6]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[6]; 
  out[7] = 0.9*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[7]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[7]+0.5*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[6]+0.5*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[6]; 
  out[8] = 0.5000000000000001*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[9]+0.5000000000000001*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[9]+0.5*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[8]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[8]; 
  out[9] = 0.9*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[9]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[9]+0.5000000000000001*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[8]+0.5000000000000001*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[8]; 
  out[10] = 0.5000000000000001*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[11]+0.5000000000000001*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[11]+0.5*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[10]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[10]; 
  out[11] = 0.9*nu[0]*xi[0]*gijJ[1]*vtsq_inv[1]*vsq[11]+0.5*gijJ[0]*nu[0]*vtsq_inv[0]*xi[0]*vsq[11]+0.5000000000000001*gijJ[0]*nu[0]*xi[0]*vtsq_inv[1]*vsq[10]+0.5000000000000001*nu[0]*vtsq_inv[0]*xi[0]*gijJ[1]*vsq[10]; 

}

