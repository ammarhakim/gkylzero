#ifndef HAMILTON
#define HAMILTON

typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

typedef struct push_input {
    double *xf, *Bn;
    const double t, *xi, top_val;
    const int top_idx;
    const evalf_t func, gunc, hunc;
    const void *ctx;
} push_input;

void trace_vertical_boundary_bottom_to_top(double *xf, double *sf, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, int max_steps, double top_val, int top_idx);

double root_top_boundary(double ds, push_input *p_ctx);

double find_ds_to_top_boundary(double *xf, double *Bn, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, double top_val, int top_idx, double ds_max);

void trace_adaptive(double *xf, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, double end, int maxSteps);

void trace(double *xf, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xi, void *ctx, double L, int N);

void push(double *xf, double *Bn, double t, double *xi, void *ctx, double ds, evalf_t func, evalf_t gunc, evalf_t hunc);

void calculate_node_positions(double *xh, double *Bn, double t, double *xni, evalf_t func, evalf_t gunc, evalf_t hunc, double ds, const double *hp, void *ctx, int len);

void calculate_derivatives(double *Fn, evalf_t func, evalf_t gunc, evalf_t hunc, double t, double *xn, void *ctx, int len);

void calculate_Bn_from_Gn(double *Bn, double *Gn, int ndim);

void calculate_Gn_from_Fn(double *G, double *F, int ndim);

#endif