/** Author: Maxwell Rosen December 2023
 * Description: This file contains the C implementation of the first order
 * Gauss Radau integrator detailed by Everheart (1985) and the IAS15 paper.
 * https://arxiv.org/abs/1409.4779. These functions are related to contour
 * tracing to follow contours of constant flux function psi
 */ 
#include <hamilton.h>
#include <gkyl_math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Gauss Radau spacings
static const double h[8] = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780,
                            0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558,
                            0.885320946839095768090359771030, 0.977520613561287501891174488626};
// Other constants
static const double rr[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278,
                              0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486,
                              0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106,
                              0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621,
                              0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520,
                              0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223,
                              0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639,
                              0.0921996667221917338008147};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321,
                             -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006,
                             -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460,
                             0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991,
                             0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553,
                             2.9061362593084293014237913, -2.7558127197720458314421588};

// Define a function pointer type for the derivative function
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

static void
func_norm(double t, const double *xn, double *fout, void *ctx, evalf_t func, evalf_t gunc, evalf_t hunc)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    double *hunc_out = malloc(sizeof(double));
    func(t, xn, func_out, ctx);
    gunc(t, xn, gunc_out, ctx);
    hunc(t, xn, hunc_out, ctx);
    double norm = sqrt(*func_out * *func_out + *gunc_out * *gunc_out + *hunc_out * *hunc_out);
    *fout = *func_out / norm;
    free(func_out);
    free(gunc_out);
    free(hunc_out);
}

static void
gunc_norm(double t, const double *xn, double *fout, void *ctx, evalf_t func, evalf_t gunc, evalf_t hunc)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    double *hunc_out = malloc(sizeof(double));
    func(t, xn, func_out, ctx);
    gunc(t, xn, gunc_out, ctx);
    hunc(t, xn, hunc_out, ctx);
    double norm = sqrt(*func_out * *func_out + *gunc_out * *gunc_out + *hunc_out * *hunc_out);
    *fout = *gunc_out / norm;
    free(func_out);
    free(gunc_out);
    free(hunc_out);
}

static void
hunc_norm(double t, const double *xn, double *fout, void *ctx, evalf_t func, evalf_t gunc, evalf_t hunc)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    double *hunc_out = malloc(sizeof(double));
    func(t, xn, func_out, ctx);
    gunc(t, xn, gunc_out, ctx);
    hunc(t, xn, hunc_out, ctx);
    double norm = sqrt(*func_out * *func_out + *gunc_out * *gunc_out + *hunc_out * *hunc_out);
    *fout = *hunc_out / norm;
    free(func_out);
    free(gunc_out);
    free(hunc_out);
}

void find_contour_nodes(double *xn, int N, evalf_t func, evalf_t gunc, evalf_t hunc, double t,
                        double *xi, void *ctx, double top_val, int top_idx)
{
    double *xf = malloc(sizeof(double[3 * (N + 1)]));
    double *sf = malloc(sizeof(double));
    xn[0] = xi[0];
    xn[1] = xi[1];
    xn[2] = xi[2];
    trace_vertical_boundary_bottom_to_top(xf, sf, func, gunc, hunc, t, xi, ctx, 10000, top_val, top_idx);
    xn[3 * (N - 1)] = xf[0];
    xn[3 * (N - 1) + 1] = xf[1];
    xn[3 * (N - 1) + 2] = xf[2];

    for (int i = 1; i < N; i++)
    {
        trace_adaptive(xf, func, gunc, hunc, t, xi, ctx, *sf * i / (N - 1), 10000);
        xn[i * 3 + 0] = xf[0];
        xn[i * 3 + 1] = xf[1];
        xn[i * 3 + 2] = xf[2];
    }
    free(xf);
    free(sf);
}

void trace_vertical_boundary_bottom_to_top(double *xf, double *sf, evalf_t func, evalf_t gunc, evalf_t hunc, double t,
                                           double *xi, void *ctx, int max_steps, double top_val, int top_idx)
{
    double ds = fabs(top_val - xf[top_idx]);
    double Bn[21] = {0};
    double xtmp[3];
    xtmp[0] = xi[0];
    xtmp[1] = xi[1];
    xtmp[2] = xi[2];
    double epsilon = 0.1;
    double min_step = fabs(top_val - xf[top_idx]) / max_steps;
    while (1)
    {
        push(xf, Bn, t, xtmp, ctx, ds, func, gunc, hunc);
        double last_B_sum = pow(fabs(Bn[6]), 1.0 / 7.0) + pow(fabs(Bn[13]), 1.0 / 7.0) + pow(fabs(Bn[20]), 1.0 / 7.0);
        if (isnan(last_B_sum))
        {
            ds = ds / 2;
            continue;
        }
        int too_much_error = last_B_sum > epsilon;
        int large_enough_steps = ds > min_step;
        if (too_much_error && large_enough_steps)
        {
            if (ds / 2 > min_step)
            {
                ds = ds / 2;
            }
        }
        else
        {
            if (xf[top_idx] > top_val)
            {
                double ds_tmp = ds;
                ds = find_ds_to_top_boundary(xf, Bn, func, gunc, hunc, t, xtmp, ctx, top_val, top_idx, ds_tmp);
                push(xf, Bn, t, xtmp, ctx, ds, func, gunc, hunc);
                *sf = *sf + ds;
                break;
            }
            else
            {
                *sf = *sf + ds;
                xtmp[0] = xf[0];
                xtmp[1] = xf[1];
                xtmp[2] = xf[2];
                ds = ds * 2;
            }
        }
    }
}

double root_top_boundary(double ds, push_input *p_ctx)
{
    // I have suspicion that this push is changing the value of p_ctx->xf and it's feeding back
    push(p_ctx->xf, p_ctx->Bn, p_ctx->t, p_ctx->xi, p_ctx->ctx,
         ds, p_ctx->func, p_ctx->gunc, p_ctx->hunc);
    return p_ctx->xf[p_ctx->top_idx] - p_ctx->top_val;
}

double find_ds_to_top_boundary(double *xf, double *Bn, evalf_t func, evalf_t gunc, evalf_t hunc, double t,
                               double *xi, void *ctx, double top_val, int top_idx, double ds_max)
{
    struct push_input p_ctx = {
        .xf = xf,
        .Bn = Bn,
        .t = t,
        .xi = xi,
        .func = func,
        .gunc = gunc,
        .hunc = hunc,
        .ctx = ctx,
        .top_val = top_val,
        .top_idx = top_idx};
    int max_iter = 10000;
    double f_upper = root_top_boundary(ds_max, &p_ctx);
    double f_lower = root_top_boundary(0, &p_ctx);
    struct gkyl_qr_res ds = gkyl_ridders(root_top_boundary, &p_ctx, 0, ds_max, f_lower, f_upper, max_iter, 1e-16);
    if (ds.status != 0)
    {
        printf("Error in root riders in finding the top contour boundary. hamilton.c");
    }
    return ds.res;
}

void trace_adaptive(double *xf, evalf_t func, evalf_t gunc, evalf_t hunc,
                    double t, double *xi, void *ctx, double end, int max_steps)
{
    double ds = end;
    double Bn[21] = {0};
    double xtmp[3];
    xtmp[0] = xi[0];
    xtmp[1] = xi[1];
    xtmp[2] = xi[2];
    double epsilon = 0.1;
    double min_step = end / max_steps;
    int i = 1;
    double s = 0;
    while (1)
    {
        push(xf, Bn, t, xtmp, ctx, ds, func, gunc, hunc);
        double last_B_sum = pow(fabs(Bn[6]), 1.0 / 7.0) + pow(fabs(Bn[13]), 1.0 / 7.0) + pow(fabs(Bn[20]), 1.0 / 7.0);
        if (isnan(last_B_sum))
        {
            ds = ds / 2;
            continue;
        }
        int too_much_error = last_B_sum > epsilon;
        int large_enough_steps = ds > min_step;
        if (too_much_error && large_enough_steps)
        {
            if (ds / 2 > min_step)
            {
                ds = ds / 2;
                continue;
            }
        }
        else
        {
            if (s + ds > end)
            {
                ds = end - s;
                continue;
            }
            s = s + ds;
            xtmp[0] = xf[0];
            xtmp[1] = xf[1];
            xtmp[2] = xf[2];
            ds = ds * 2;
            if (s >= end)
            {
                break;
            }
            i++;
        }
    }
}

void trace(double *xf, evalf_t func, evalf_t gunc, evalf_t hunc,
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
        push(xf, Bn, t, xtmp, ctx, ds, func, gunc, hunc);
        xtmp[0] = xf[0];
        xtmp[1] = xf[1];
        xtmp[2] = xf[2];
    }
}

void push(double *xf, double *Bn, double t, double *xi, void *ctx, double ds,
          evalf_t func, evalf_t gunc, evalf_t hunc)
{
    int len_hF = 8;
    int len_BG = 7;
    double *xh = malloc(sizeof(double[3 * len_hF]));
    double *Fn = malloc(sizeof(double[3 * len_hF]));
    double *Gn = malloc(sizeof(double[3 * len_BG]));
    double Bf[21] = {0};
    double total_diff;
    for (int i = 0; i < 20; i++)
    {
        calculate_node_positions(xh, Bn, t, xi, func, gunc, hunc, ds, h, ctx, len_hF);
        calculate_derivatives(Fn, func, gunc, hunc, t, xh, ctx, len_hF);
        calculate_Gn_from_Fn(Gn, Fn, 3);
        calculate_Bn_from_Gn(Bf, Gn, 3);
        total_diff = 0;
        for (int j = 0; j < 21; j++)
        {
            total_diff += pow(fabs(Bf[j] - Bn[j]), 2);
            Bn[j] = Bf[j];
        }
        if (total_diff < 1e-16)
        {
            break;
        }
    }
    double hf = 1.0;
    calculate_node_positions(xf, Bn, t, xi, func, gunc, hunc, ds, &hf, ctx, 1);
    free(xh);
    free(Fn);
    free(Gn);
}

void calculate_node_positions(double *xh, double *Bn, double t, double *xni, evalf_t func,
                              evalf_t gunc, evalf_t hunc, double ds, const double *hp, void *ctx, int len)
{
    double *F1x = malloc(sizeof(double));
    double *F1y = malloc(sizeof(double));
    double *F1z = malloc(sizeof(double));
    func_norm(t, xni, F1x, ctx, func, gunc, hunc);
    gunc_norm(t, xni, F1y, ctx, func, gunc, hunc);
    hunc_norm(t, xni, F1z, ctx, func, gunc, hunc);
    for (int n = 0; n < len; n++)
    {
        int i = 0;
        xh[n * 3 + i] = xni[i] + hp[n] * ds * (*F1x + hp[n] * 1 / 2 * (Bn[0 + 7 * i] + hp[n] * 2 / 3 * (Bn[1 + 7 * i] + hp[n] * 3 / 4 * (Bn[2 + 7 * i] + hp[n] * 4 / 5 * (Bn[3 + 7 * i] + hp[n] * 5 / 6 * (Bn[4 + 7 * i] + hp[n] * 6 / 7 * (Bn[5 + 7 * i] + hp[n] * 7 / 8 * Bn[6 + 7 * i])))))));
        i = 1;
        xh[n * 3 + i] = xni[i] + hp[n] * ds * (*F1y + hp[n] * 1 / 2 * (Bn[0 + 7 * i] + hp[n] * 2 / 3 * (Bn[1 + 7 * i] + hp[n] * 3 / 4 * (Bn[2 + 7 * i] + hp[n] * 4 / 5 * (Bn[3 + 7 * i] + hp[n] * 5 / 6 * (Bn[4 + 7 * i] + hp[n] * 6 / 7 * (Bn[5 + 7 * i] + hp[n] * 7 / 8 * Bn[6 + 7 * i])))))));
        i = 2;
        xh[n * 3 + i] = xni[i] + hp[n] * ds * (*F1z + hp[n] * 1 / 2 * (Bn[0 + 7 * i] + hp[n] * 2 / 3 * (Bn[1 + 7 * i] + hp[n] * 3 / 4 * (Bn[2 + 7 * i] + hp[n] * 4 / 5 * (Bn[3 + 7 * i] + hp[n] * 5 / 6 * (Bn[4 + 7 * i] + hp[n] * 6 / 7 * (Bn[5 + 7 * i] + hp[n] * 7 / 8 * Bn[6 + 7 * i])))))));
    }
    free(F1x);
    free(F1y);
    free(F1z);
}

void calculate_derivatives(double *Fn, evalf_t func, evalf_t gunc, evalf_t hunc,
                           double t, double *xn, void *ctx, int len)
{
    double *func_out = malloc(sizeof(double));
    double *gunc_out = malloc(sizeof(double));
    double *hunc_out = malloc(sizeof(double));
    for (int i = 0; i < len; i++)
    {
        func_norm(t, &xn[i * 3], func_out, ctx, func, gunc, hunc);
        gunc_norm(t, &xn[i * 3], gunc_out, ctx, func, gunc, hunc);
        hunc_norm(t, &xn[i * 3], hunc_out, ctx, func, gunc, hunc);
        Fn[0 + i] = *func_out;
        Fn[8 + i] = *gunc_out;
        Fn[16 + i] = *hunc_out;
    }
    free(func_out);
    free(gunc_out);
    free(hunc_out);
}

void calculate_Bn_from_Gn(double *Bn, double *Gn, int ndim)
{
    for (int i = 0; i < ndim; i++)
    {
        Bn[7 * i + 0] = Gn[7 * i + 0] + c[0] * Gn[7 * i + 1] + c[1] * Gn[7 * i + 2] + c[3] * Gn[7 * i + 3] + c[6] * Gn[7 * i + 4] + c[10] * Gn[7 * i + 5] + c[15] * Gn[7 * i + 6];
        Bn[7 * i + 1] = Gn[7 * i + 1] + c[2] * Gn[7 * i + 2] + c[4] * Gn[7 * i + 3] + c[7] * Gn[7 * i + 4] + c[11] * Gn[7 * i + 5] + c[16] * Gn[7 * i + 6];
        Bn[7 * i + 2] = Gn[7 * i + 2] + c[5] * Gn[7 * i + 3] + c[8] * Gn[7 * i + 4] + c[12] * Gn[7 * i + 5] + c[17] * Gn[7 * i + 6];
        Bn[7 * i + 3] = Gn[7 * i + 3] + c[9] * Gn[7 * i + 4] + c[13] * Gn[7 * i + 5] + c[18] * Gn[7 * i + 6];
        Bn[7 * i + 4] = Gn[7 * i + 4] + c[14] * Gn[7 * i + 5] + c[19] * Gn[7 * i + 6];
        Bn[7 * i + 5] = Gn[7 * i + 5] + c[20] * Gn[7 * i + 6];
        Bn[7 * i + 6] = Gn[7 * i + 6];
    }
}

void calculate_Gn_from_Fn(double *G, double *F, int ndim)
{
    for (int i = 0; i < ndim; i++)
    {
        G[7 * i + 0] = (F[8 * i + 1] - F[8 * i + 0]) / rr[0];
        G[7 * i + 1] = ((F[8 * i + 2] - F[8 * i + 0]) / rr[1] - G[7 * i + 0]) / rr[2];
        G[7 * i + 2] = (((F[8 * i + 3] - F[8 * i + 0]) / rr[3] - G[7 * i + 0]) / rr[4] - G[7 * i + 1]) / rr[5];
        G[7 * i + 3] = ((((F[8 * i + 4] - F[8 * i + 0]) / rr[6] - G[7 * i + 0]) / rr[7] - G[7 * i + 1]) / rr[8] - G[7 * i + 2]) / rr[9];
        G[7 * i + 4] = (((((F[8 * i + 5] - F[8 * i + 0]) / rr[10] - G[7 * i + 0]) / rr[11] - G[7 * i + 1]) / rr[12] - G[7 * i + 2]) / rr[13] - G[7 * i + 3]) / rr[14];
        G[7 * i + 5] = ((((((F[8 * i + 6] - F[8 * i + 0]) / rr[15] - G[7 * i + 0]) / rr[16] - G[7 * i + 1]) / rr[17] - G[7 * i + 2]) / rr[18] - G[7 * i + 3]) / rr[19] - G[7 * i + 4]) / rr[20];
        G[7 * i + 6] = (((((((F[8 * i + 7] - F[8 * i + 0]) / rr[21] - G[7 * i + 0]) / rr[22] - G[7 * i + 1]) / rr[23] - G[7 * i + 2]) / rr[24] - G[7 * i + 3]) / rr[25] - G[7 * i + 4]) / rr[26] - G[7 * i + 5]) / rr[27];
    }
}