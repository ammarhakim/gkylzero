local ffi = require "ffi"

ffi.cdef [[
enum gkyl_elem_type { GKYL_INT, GKYL_INT_64, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
  void (*free)(const struct gkyl_ref_count* );
  int count;
};

/**
 * Type of function to project.
 *
 * @param t Time to evaluate function
 * @param xn Coordinates for evaluation
 * @param fout Output vector of 'num_ret_vals'
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

]]
