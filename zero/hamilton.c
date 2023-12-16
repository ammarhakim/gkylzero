// Author: Maxwell Rosen 12/7/2023
// Description: This file contains the C implementation of the first order
//              Gauss Radau integrator. So far I'm just translating it
//              directly from the python implementation. It will be a standalone
//              C file.
#include "hamilton.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Gauss Radau spacings
static const double h[8] = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};
// Other constants
static const double rr[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};

// Define a function pointer type for the derivative function
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

static void
func_norm(double t, const double *xn, double *fout, void *ctx, evalf_t func, evalf_t gunc)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    func(t, xn, func_out, ctx);
    gunc(t, xn, gunc_out, ctx);
    double norm = sqrt(*func_out * *func_out + *gunc_out * *gunc_out);
    *fout = *func_out / norm;
    free(func_out);
    free(gunc_out);
}


static void
gunc_norm(double t, const double *xn, double *fout, void *ctx, evalf_t func, evalf_t gunc)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    func(t, xn, func_out, ctx);
    gunc(t, xn, gunc_out, ctx);
    double norm = sqrt(*func_out * *func_out + *gunc_out * *gunc_out);
    *fout = *gunc_out / norm;
    free(func_out);
    free(gunc_out);
}

// void adaptiveTrace(double *xf, double *yf, evalf_t func, evalf_t gunc, double xi, double yi, double end, int maxSteps)
// {
//     double ds = end;
//     double Bx[7] = {0};
//     double By[7] = {0};
//     double epsilon = 0.1;
//     double minStep = end / maxSteps;
//     int i = 1;
//     double s = 0;
//     while (1)
//     {
//         push(xf, yf, Bx, By, xi, yi, ds, func, gunc);
//         double lastBsum = pow(fabs(Bx[6]), 1.0 / 7) + pow(fabs(By[6]), 1.0 / 7);
//         int tooMuchError = lastBsum > epsilon;
//         int largeEnoughSteps = ds > minStep;
//         if (tooMuchError && largeEnoughSteps)
//         {
//             if (ds / 2 > minStep)
//             {
//                 ds = ds / 2;
//                 continue;
//             }
//         }
//         else
//         {
//             if (s + ds > end)
//             {
//                 ds = end - s;
//                 continue;
//             }
//             s = s + ds;
//             xi = *xf;
//             yi = *yf;
//             ds = ds * 2;
//             if (s >= end)
//             {
//                 break;
//             }
//             i++;
//         }
//     }
// }

// void trace(double *xf, double *yf, evalf_t func, evalf_t gunc,
//            double xi, double yi, double L, int N)
// {
//     double ds = L / N;
//     double Bx[7] = {0};
//     double By[7] = {0};
//     for (int i = 1; i <= N; i++)
//     {
//         push(xf, yf, Bx, By, xi, yi, ds, func, gunc);
//         xi = *xf;
//         yi = *yf;
//     }
// }

void trace(double *xf, evalf_t func, evalf_t gunc,
           double t, double *xi, void *ctx, double L, int N)
{
    double ds = L / N;
    double Bn[21] = {0};
    double xtmp[3];
    xtmp[0] = xi[0];
    xtmp[1] = xi[1];
    xtmp[2] = xi[2];
    for (int i = 1; i <= N; i++)
    {
        push(xf, Bn, t, xtmp, ctx, ds, func, gunc);
        xtmp[0] = xf[0];
        xtmp[1] = xf[1];
        xtmp[2] = xf[2];
    }
}

void push(double *xf, double *Bn,
          double t, double *xi, void *ctx, double ds, evalf_t func, evalf_t gunc)
{
    int len_hF = 8;
    int len_BG = 7;
    double *xh = (double *)malloc(3 * len_hF * sizeof(double));
    double *Fn = (double *)malloc(3 * len_hF * sizeof(double));
    double *Gn = (double *)malloc(3 * len_BG * sizeof(double));
    for (int i = 0; i < 12; i++)
    {
        calculateNodePositions(xh, Bn, t, xi, func, gunc, ds, h, ctx, len_hF);
        calculateDerivatives(Fn, func, gunc, t, xh, ctx, len_hF);
        calculateGnFromFn(Gn, Fn, 3);
        calculateBn(Bn, Gn, 3);
    }
    double hf = 1.0;
    calculateNodePositions(xf, Bn, t, xi, func, gunc, ds, &hf, ctx, 1);
    free(xh);
    free(Fn);
    free(Gn);
}

void calculateNodePositions(double *xh, double *Bn,
                            double t, double *xni, evalf_t func, evalf_t gunc,
                            double ds, const double *hp, void *ctx, int len)
{
    double *F1x = malloc(sizeof(double));
    double *F1y = malloc(sizeof(double));
    func_norm(t, xni, F1x, ctx, func, gunc);
    gunc_norm(t, xni, F1y, ctx, func, gunc);
    for (int n = 0; n < len; n++)
    {
        int i = 0;
        xh[n * 3 + i] = xni[i] + hp[n] * ds * (*F1x + hp[n] * 1 / 2 * (Bn[0 + 7 * i] + hp[n] * 2 / 3 * (Bn[1 + 7 * i] + hp[n] * 3 / 4 * (Bn[2 + 7 * i] + hp[n] * 4 / 5 * (Bn[3 + 7 * i] + hp[n] * 5 / 6 * (Bn[4 + 7 * i] + hp[n] * 6 / 7 * (Bn[5 + 7 * i] + hp[n] * 7 / 8 * Bn[6 + 7 * i])))))));
        i = 1;
        xh[n * 3 + i] = xni[i] + hp[n] * ds * (*F1y + hp[n] * 1 / 2 * (Bn[0 + 7 * i] + hp[n] * 2 / 3 * (Bn[1 + 7 * i] + hp[n] * 3 / 4 * (Bn[2 + 7 * i] + hp[n] * 4 / 5 * (Bn[3 + 7 * i] + hp[n] * 5 / 6 * (Bn[4 + 7 * i] + hp[n] * 6 / 7 * (Bn[5 + 7 * i] + hp[n] * 7 / 8 * Bn[6 + 7 * i])))))));
        i = 2;
        xh[n * 3 + i] = xni[i] + hp[n] * ds * (0 + hp[n] * 1 / 2 * (Bn[0 + 7 * i] + hp[n] * 2 / 3 * (Bn[1 + 7 * i] + hp[n] * 3 / 4 * (Bn[2 + 7 * i] + hp[n] * 4 / 5 * (Bn[3 + 7 * i] + hp[n] * 5 / 6 * (Bn[4 + 7 * i] + hp[n] * 6 / 7 * (Bn[5 + 7 * i] + hp[n] * 7 / 8 * Bn[6 + 7 * i])))))));
    }
    free(F1x);
    free(F1y);
}

void calculateDerivatives(double *Fn, evalf_t func, evalf_t gunc,
                          double t, double *xn, void *ctx, int len)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    for (int i = 0; i < len; i++)
    {
        func_norm(t, &xn[i * 3], func_out, ctx, func, gunc);
        gunc_norm(t, &xn[i * 3], gunc_out, ctx, func, gunc);
        Fn[0 + i] = *func_out;
        Fn[8 + i] = *gunc_out;
    }
    free(func_out);
    free(gunc_out);
}

void calculateBn(double *Bn, double *Gn, int ndim)
{
    for (int i = 0; i < ndim; i++)
    {
        Bn[7 * i + 0] = Gn[7 * i + 0] + c[0]  * Gn[7 * i + 1] + c[1]  * Gn[7 * i + 2] + c[3]  * Gn[7 * i + 3] + c[6]  * Gn[7 * i + 4] + c[10] * Gn[7 * i + 5] + c[15] * Gn[7 * i + 6];
        Bn[7 * i + 1] = Gn[7 * i + 1] + c[2]  * Gn[7 * i + 2] + c[4]  * Gn[7 * i + 3] + c[7]  * Gn[7 * i + 4] + c[11] * Gn[7 * i + 5] + c[16] * Gn[7 * i + 6];
        Bn[7 * i + 2] = Gn[7 * i + 2] + c[5]  * Gn[7 * i + 3] + c[8]  * Gn[7 * i + 4] + c[12] * Gn[7 * i + 5] + c[17] * Gn[7 * i + 6];
        Bn[7 * i + 3] = Gn[7 * i + 3] + c[9]  * Gn[7 * i + 4] + c[13] * Gn[7 * i + 5] + c[18] * Gn[7 * i + 6];
        Bn[7 * i + 4] = Gn[7 * i + 4] + c[14] * Gn[7 * i + 5] + c[19] * Gn[7 * i + 6];
        Bn[7 * i + 5] = Gn[7 * i + 5] + c[20] * Gn[7 * i + 6];
        Bn[7 * i + 6] = Gn[7 * i + 6];
    }
}

void calculateGnFromFn(double *G, double *F, int ndim)
{
    for (int i = 0; i < ndim; i++)
    {
        G[7 * i + 0] =       (F[8 * i + 1] - F[8 * i + 0]) / rr[0];
        G[7 * i + 1] =      ((F[8 * i + 2] - F[8 * i + 0]) / rr[1]  - G[7 * i + 0]) / rr[2];
        G[7 * i + 2] =     (((F[8 * i + 3] - F[8 * i + 0]) / rr[3]  - G[7 * i + 0]) / rr[4]  - G[7 * i + 1]) / rr[5];
        G[7 * i + 3] =    ((((F[8 * i + 4] - F[8 * i + 0]) / rr[6]  - G[7 * i + 0]) / rr[7]  - G[7 * i + 1]) / rr[8]  - G[7 * i + 2]) / rr[9];
        G[7 * i + 4] =   (((((F[8 * i + 5] - F[8 * i + 0]) / rr[10] - G[7 * i + 0]) / rr[11] - G[7 * i + 1]) / rr[12] - G[7 * i + 2]) / rr[13] - G[7 * i + 3]) / rr[14];
        G[7 * i + 5] =  ((((((F[8 * i + 6] - F[8 * i + 0]) / rr[15] - G[7 * i + 0]) / rr[16] - G[7 * i + 1]) / rr[17] - G[7 * i + 2]) / rr[18] - G[7 * i + 3]) / rr[19] - G[7 * i + 4]) / rr[20];
        G[7 * i + 6] = (((((((F[8 * i + 7] - F[8 * i + 0]) / rr[21] - G[7 * i + 0]) / rr[22] - G[7 * i + 1]) / rr[23] - G[7 * i + 2]) / rr[24] - G[7 * i + 3]) / rr[25] - G[7 * i + 4]) / rr[26] - G[7 * i + 5]) / rr[27];
    }
}