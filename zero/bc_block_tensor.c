#include <gkyl_bc_block_tensor.h>
#include <gkyl_bc_block_tensor_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_util.h>

static inline double dot_product(const double *v1, const double *v2)
{
  double out = 0.0;
  for(int i = 0; i < 3; i++)
    out += v1[i]*v2[i];
  return out;
}


/**
 * Given the modal expansion of the tangent vectors in the block which fluxes are leaving
 * and the duals in the block which the fluxes are entering
 * Compute 
 * simulations.
 *
 * @param inp Input parameters
 * @param New GK geometry updater
 */

struct bc_block_tensor*
gkyl_bc_block_tensor_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, bool use_gpu)
{
  struct bc_block_tensor *up = gkyl_malloc(sizeof(* up));
  up->basis = *basis;
  up->range = *range;
  up->range_ext = *range_ext;
  up->grid = *grid;
  up->cdim = grid->ndim;
  up->poly_order = basis->poly_order;

  struct gkyl_basis surf_basis;
  gkyl_cart_modal_tensor(&surf_basis, up->cdim-1, up->poly_order);
  up->num_surf_nodes = surf_basis.num_basis;

  return up;
}

/**
 * Take in modal expansions of duals of one block and tangents of the other (cartesian components)
 * and calculate T^j'_i = e^j' \dot e_i at the quadrature nodes of the interface
 * @param up bc_block_tensor object
 * @param edge1 edge of block which fluxes leave (0 is lower, 1 is upper)
 * @param edge2 edge of block which fluxes enter(0 is lower, 1 is upper)
 * @param ej duals of block which fluxes enter
 * @param e_i tangent vectors of block which fluxes leave
 */
void calc_tensor(struct bc_block_tensor *up, int dir, int edge1, int edge2, const double *ej, const double *e_i, double *tj_i)
{
  // First evaluate at all the quadrature nodes
  double ej_surf[up->num_surf_nodes][9];
  double e_i_surf[up->num_surf_nodes][9];
  for(int n = 0; n < up->num_surf_nodes; n++) {
    for(int i = 0; i < 9; i++) {
      e_i_surf[n][i] = bc_block_tensor_choose_kernel(up->cdim, edge1, dir, n) (&e_i[i*up->basis.num_basis]);
      ej_surf[n][i] = bc_block_tensor_choose_kernel(up->cdim, edge2, dir, n) (&ej[i*up->basis.num_basis]);
    }
  }

  // Now take the dot prduct at the quadrature nodes and fill the tensor
  // Only Need T11,13,31,33 in 2d
  int jctr=0;
  for (int j = 0; j < 3; j++){
    if(up->cdim==2 &&  j==1)
      continue;
    int ictr=0;
    for (int i = 0; i < 3; i++){
      if(up->cdim==2 &&  i==1)
        continue;
      for(int n = 0; n < up->num_surf_nodes; n++) {
        tj_i[up->cdim*up->num_surf_nodes*jctr + ictr*up->num_surf_nodes + n] = dot_product(&ej_surf[n][3*j], &e_i_surf[n][3*i]);
        //printf("j,i = %d, %d\n", j,i);
        //printf("t = %g\n", tj_i[up->cdim*up->num_surf_nodes*jctr + ictr*up->num_surf_nodes + n]);
        //printf("index = %d\n", up->cdim*up->num_surf_nodes*jctr + ictr*up->num_surf_nodes + n);
        //printf("jctr,ictr,n = %d, %d, %d\n\n", jctr,ictr,n);
      }
      ictr+=1;
    }
    jctr+=1;
  }

}
