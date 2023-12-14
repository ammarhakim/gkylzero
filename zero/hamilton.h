#ifndef HAMILTON
#define HAMILTON

typedef double (*derivative_func)(double, double);

// Define a function pointer type for the derivative function
static double func_norm(double x, double y, derivative_func func, derivative_func gunc);

static double gunc_norm(double x, double y, derivative_func func, derivative_func gunc);

void adaptiveTrace(double *xf, double *yf, derivative_func func, derivative_func gunc, double xi, double yi, double end, int maxSteps);

void trace(double *xf, double *yf, derivative_func func, derivative_func gunc,
           double xi, double yi, double L, int N);

void push(double *xf, double *yf, double *Bx, double *By, double xi, double yi, double ds, derivative_func func, derivative_func gunc);

void calculateNodePositions(double *xh, double *yh, double *Bx, double *By, double xi, double yi, derivative_func func, derivative_func gunc, double ds, const double *hp, int len);

void calculateDerivatives(double *Fx, double *Fy, derivative_func func, derivative_func gunc, double *x, double *y, int len);

void calculateB(double *B, double *G);

void calculateGFromF(double *G, double *F);

#endif