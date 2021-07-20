#pragma once

/**
 * Basis function object: treat this as a read-only object.
 */
struct gkyl_basis {
  unsigned ndim, polyOrder, numBasis;
  char id[64]; // "serendipity", "tensor", "maximal-order"
    
/**
 * Evaluate basis in unit cell (i.e. a hypercube with each side
 * [-1,1])
 *
 * @param z Location to evaluate basis. z \in [-1,1]^n
 * @param b On output, value of basis at 'z'
 */
  void (*eval)(const double *z, double *b);

/**
 * Flip-sign function: changes signs of input expansion cofficients by
 * changing monomial terms in specified direction.
 *
 * @param dir Direction to flip sign
 * @param f Input expansion
 * @param fout On output, flipped version of @a f
 */
  void (*flip_sign)(int dir, const double *f, double *fout);
};

/**
 * Create new modal serendipity basis function object.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param polyOrder Polynomial order.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim, int polyOrder);
