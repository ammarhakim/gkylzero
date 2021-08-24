#pragma once

/* Basis function identifiers */
enum gkyl_basis_type {
  GKYL_BASIS_MODAL_SERENDIPITY,
  GKYL_BASIS_MODAL_TENSOR
};

/**
 * Basis function object
 */
struct gkyl_basis {
  unsigned ndim, poly_order, num_basis;
  char id[64]; // "serendipity", "tensor", "maximal-order"
  enum gkyl_basis_type type; // identifier for basis function
    
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
 * @param poly_order Polynomial order.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_serendip(struct gkyl_basis *basis, int ndim, int poly_order);

/**
 * Create new modal tensor-product basis function object.
 *
 * @param basis Basis object to initialize
 * @param ndim Dimension of reference element.
 * @param poly_order Polynomial order.
 * @return Pointer to new basis function.
 */
void gkyl_cart_modal_tensor(struct gkyl_basis *basis, int ndim, int poly_order);
