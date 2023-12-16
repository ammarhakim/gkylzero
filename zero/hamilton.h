#ifndef HAMILTON
#define HAMILTON

typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

// Define a function pointer type for the derivative function
void push(double *xf, double *Bn, double t, double *xi, void *ctx, double ds, evalf_t func, evalf_t gunc);

void calculateNodePositions(double *xh, double *Bn, double t, double *xni, evalf_t func, evalf_t gunc, double ds, const double *hp, void *ctx, int len);

void calculateDerivatives(double *Fn, evalf_t func, evalf_t gunc, double t, double *xn, void *ctx, int len);

void calculateBn(double *Bn, double *Gn, int ndim);

void calculateGnFromFn(double *G, double *F, int ndim);

#endif