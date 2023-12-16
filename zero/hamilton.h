#ifndef HAMILTON
#define HAMILTON

typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

void trace_adaptive(double *xf, evalf_t func, evalf_t gunc, double t, double *xi, void *ctx, double end, int maxSteps);

void trace(double *xf, evalf_t func, evalf_t gunc, double t, double *xi, void *ctx, double L, int N);

void push(double *xf, double *Bn, double t, double *xi, void *ctx, double ds, evalf_t func, evalf_t gunc);

void calculate_node_positions(double *xh, double *Bn, double t, double *xni, evalf_t func, evalf_t gunc, double ds, const double *hp, void *ctx, int len);

void calculate_derivatives(double *Fn, evalf_t func, evalf_t gunc, double t, double *xn, void *ctx, int len);

void calculate_Bn_from_Gn(double *Bn, double *Gn, int ndim);

void calculate_Gn_from_Fn(double *G, double *F, int ndim);

#endif