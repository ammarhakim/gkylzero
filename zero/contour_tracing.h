#ifndef CONTOUR_TRACING
#define CONTOUR_TRACING

typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

typedef struct push_input {
    double *xf, *Bn;
    double t, *xi, top_val;
    int top_idx;
    evalf_t func, gunc, hunc;
    void *ctx;
} push_input;

/**
 * Compute the nodes along a contour that obeys the functions func, gunc, and hunc.
 * The contour is traced from the bottom boundary to the top boundary and nodes are put between
 * the bottom and top boundaries.
 *
 * @param xn 3D array of nodes, including start and end
 * @param N Number of nodes
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
 * @param t time
 * @param xi 3D array of initial position along the bottom boundary
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param top_val value of the contour trace at the top boundary
 * @param top_idx index of the top_val of the contour trace at the top boundary
 */
void find_contour_nodes(double *xn, int N, evalf_t func, evalf_t gunc, evalf_t hunc, double t,
                        double *xi, void *ctx, double top_val, int top_idx);

/**
 * Traces the contour from the bottom boundary to the top boundary and returns the final position
 * 
 * @param xf 3D array of the final position of the contour trace
 * @param sf single number which is the length that the contour has traveled
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
 * @param t time
 * @param xi 3D array of initial position along the bottom boundary
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param max_steps maximum number of steps to take
 * @param top_val value of the contour trace at the top boundary
 * @param top_idx index of the top_val of the contour trace at the top boundary
*/
void trace_vertical_boundary_bottom_to_top(double *xf, double *sf, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, int max_steps, double top_val, int top_idx);

/**
 * Function passed for root finding for the boundary of the contour.
 * Optimizes to find the optimal step size to take in order to land on the top boundary
 * 
 * @param ds step size
 * @param p_ctx pointer to push_input struct that is passed to the function push
*/
double root_top_boundary(double ds, push_input *p_ctx);

/**
 * Pushes the contour trace in the final step onto the boundary
 * 
 * @param xf 3D array of the final position of the contour trace
 * @param Bn 3D array to do the pushing 
 * @param t time
 * @param xi 3D array of initial position
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
 * @param top_val value of the contour trace at the top boundary
 * @param top_idx index of the top_val of the contour trace at the top boundary
 * @param ds_max maximum step size to take
*/
double find_ds_to_top_boundary(double *xf, double *Bn, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, double top_val, int top_idx, double ds_max);

/**
 * Traces the contour to go distance "end" and returns the final position.
 * Reffer to the IAS15 paper for more details.
 * https://arxiv.org/abs/1409.4779
 * This pusher makes some minor changes to the original IAS15 pusher.
 * 
 * @param xf 3D array of the final position of the contour trace
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
 * @param t time
 * @param xi 3D array of initial position along the top boundary
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param end distance to travel
 * @param max_steps maximum number of steps to take
*/
void trace_adaptive(double *xf, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, double end, int maxSteps);

/**
 * Traces the contour over a fixed number of even steps returns the final position
 * 
 * @param xf 3D array of the final position of the contour trace
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
 * @param t time
 * @param xi 3D array of initial position along the top boundary
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param L distance to travel
 * @param N Number of steps to take
*/
void trace(double *xf, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, double L, int N);

/**
 * Pushes the contour trace to the next position
 * 
 * @param xf 3D array of the final position of the contour trace
 * @param Bn 3D array to do the pushing 
 * @param t time
 * @param ds step size
 * @param xi 3D array of initial position
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
*/
void push(double *xf, double *Bn, double t, double *xi, void *ctx, double ds, evalf_t func, evalf_t gunc, evalf_t hunc);

/**
 * Pushes the contour trace to the next positions
 * 
 * @param xh 3D * N array of the final position of the contour trace
 * @param Bn 3D array to do the pushing 
 * @param t time
 * @param xni 3D array of initial position
 * @param ds step size
 * @param hp 3D array of the amounts to push
 * @param len length of the arrays hp
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
*/
void calculate_node_positions(double *xh, double *Bn, double t, double *xni, evalf_t func, evalf_t gunc, evalf_t hunc, double ds, const double *hp, void *ctx, int len);

/**
 * Calculates the derivatives of the contour functions
 * 
 * @param Fn 3D array of the derivatives of the contour functions
 * @param func contour function in the 0th direction
 * @param gunc contour function in the 1st direction
 * @param hunc contour function in the 2nd direction
 * @param t time
 * @param xn 3D array of the position of the contour trace
 * @param ctx Context to pass to functions func, gunc, and hunc
 * @param len length of the arrays hp
*/
void calculate_derivatives(double *Fn, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xn, void *ctx, int len);

/**
 * Calculates the derivatives of the contour functions
 * 
 * @param Bn array of all the B's which are involved in the pusher
 * @param Gn array of all the G's which are involved in the pusher
 * @param ndim number of dimensions of the pusher (1, 2, or 3)
*/
void calculate_Bn_from_Gn(double *Bn, double *Gn, int ndim);

/**
 * Calculates the derivatives of the contour functions
 * 
 * @param Gn array of all the G's which are involved in the pusher
 * @param Fn array of all the Forces which are involved in the pusher
 * @param ndim number of dimensions of the pusher (1, 2, or 3)
*/
void calculate_Gn_from_Fn(double *G, double *F, int ndim);

#endif