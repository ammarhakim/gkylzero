#include "hamilton.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

// Build using gcc test_hamilton.c hamilton.c -o test_hamilton -lm && ./test_hamilton

// Define a function pointer type for the derivative function
static const double h[8] = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

void dyds_circ(double t, const double *xn, double *fout, void *ctx)
{
    *fout = xn[0];
}

void dxds_circ(double t, const double *xn, double *fout, void *ctx)
{
    *fout = -xn[1];
}

void dyds_double_circ(double t, const double *xn, double *fout, void *ctx)
{
    *fout = 2 * xn[0];
}

void dxds_double_circ(double t, const double *xn, double *fout, void *ctx)
{
    *fout = -2 * xn[1];
}

void dyds_hyperbola(double t, const double *xn, double *fout, void *ctx)
{
    *fout = xn[0];
}

void dxds_hyperbola(double t, const double *xn, double *fout, void *ctx)
{
    *fout = xn[1];
}


void test_1()
{
    // Test calculateGFromF function
    double F_test[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    double G_test[7];
    calculateGnFromFn(G_test, F_test, 1);

    double G_known[7] = {
        17.77380891,
        -53.86059148,
        131.06888234,
        -217.79397728,
        281.20441928,
        -298.65875394,
        340.307225};

    for (int i = 0; i < 7; i++)
    {
        assert(fabs(G_test[i] - G_known[i]) < 1e-8);
    }
}

void test_2()
{
    // Test calculateB function
    double G_test[7] = {1, 2, 3, 4, 5, 6, 7};
    double B_test[7];
    calculateBn(B_test, G_test, 1);

    double B_known[7] = {
        0.91365987,
        1.37249275,
        3.08903225,
        -4.44869317,
        14.12000318,
        -13.29068904,
        7.0};

    // Use assert to check that the test passed
    for (int i = 0; i < 7; i++)
    {
        assert(fabs(B_test[i] - B_known[i]) < 1e-8);
    }
}

void test_3()
{
    // Test calculateGFromF function
    double F_test[16] = {1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8};
    double G_test[14];
    calculateGnFromFn(G_test, F_test, 2);

    double G_known[14] = {
        17.77380891,
        -53.86059148,
        131.06888234,
        -217.79397728,
        281.20441928,
        -298.65875394,
        340.307225,
        17.77380891,
        -53.86059148,
        131.06888234,
        -217.79397728,
        281.20441928,
        -298.65875394,
        340.307225};

    for (int i = 0; i < 14; i++)
    {
        assert(fabs(G_test[i] - G_known[i]) < 1e-8);
    }
}

void test_4()
{
    // Test calculateB function
    double G_test[14] = {1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7};
    double B_test[14];
    calculateBn(B_test, G_test, 2);

    double B_known[14] = {
        0.91365987,
        1.37249275,
        3.08903225,
        -4.44869317,
        14.12000318,
        -13.29068904,
        7.0,
        0.91365987,
        1.37249275,
        3.08903225,
        -4.44869317,
        14.12000318,
        -13.29068904,
        7.0};

    // Use assert to check that the test passed
    for (int i = 0; i < 14; i++)
    {
        assert(fabs(B_test[i] - B_known[i]) < 1e-8);
    }
}


void test_5()
{
    // Test the pusher goes around a circle radius 1 to pi
    double xi[3] = {1, 0, 0};
    double xf[3];
    double s_final = M_PI;
    int nSteps = 30;
    double ds = s_final / nSteps;
    double x_trace[(nSteps + 1) * 3];
    double Bn[21] = {0};
    double t = 0;
    void *ctx = NULL;
    x_trace[0] = xi[0];
    x_trace[1] = xi[1];
    x_trace[2] = xi[2];
    for (int i = 1; i < nSteps + 1; i++)
    {
        xi[0] = x_trace[(i - 1) * 3 + 0];
        xi[1] = x_trace[(i - 1) * 3 + 1];
        xi[2] = x_trace[(i - 1) * 3 + 2];
        push(xf, Bn, t, xi, ctx, ds, dxds_circ, dyds_circ);
        x_trace[i * 3 + 0] = xf[0];
        x_trace[i * 3 + 1] = xf[1];
        x_trace[i * 3 + 2] = xf[2];
    }
    assert(fabs(xf[0] + 1) < 1e-8);
    assert(fabs(xf[1]) < 1e-8);
}


void test_6()
{
    // Test the pusher goes around a circle radius 1 to pi
    // with double derivatives
    double xi[3] = {1, 0, 0};
    double xf[3];
    double s_final = M_PI;
    int nSteps = 100;
    double ds = s_final / nSteps;
    double x_trace[(nSteps + 1) * 3];
    double Bn[21] = {0};
    double t = 0;
    void *ctx = NULL;
    x_trace[0] = xi[0];
    x_trace[1] = xi[1];
    x_trace[2] = xi[2];
    for (int i = 1; i < nSteps + 1; i++)
    {
        xi[0] = x_trace[(i - 1) * 3 + 0];
        xi[1] = x_trace[(i - 1) * 3 + 1];
        xi[2] = x_trace[(i - 1) * 3 + 2];
        push(xf, Bn, t, xi, ctx, ds, dxds_double_circ, dyds_double_circ);
        x_trace[i * 3 + 0] = xf[0];
        x_trace[i * 3 + 1] = xf[1];
        x_trace[i * 3 + 2] = xf[2];
    }
    assert(fabs(xf[0] + 1) < 1e-8);
    assert(fabs(xf[1]) < 1e-8);
}

void test_7()
{
    // Test the pusher goes around a circle radius 2 to pi/2
    double xi[3] = {2, 0, 0};
    double xf[3];
    double s_final = M_PI;
    int nSteps = 30;
    double ds = s_final / nSteps;
    double x_trace[(nSteps + 1) * 3];
    double Bn[21] = {0};
    double t = 0;
    void *ctx = NULL;
    x_trace[0] = xi[0];
    x_trace[1] = xi[1];
    x_trace[2] = xi[2];
    for (int i = 1; i < nSteps + 1; i++)
    {
        xi[0] = x_trace[(i - 1) * 3 + 0];
        xi[1] = x_trace[(i - 1) * 3 + 1];
        xi[2] = x_trace[(i - 1) * 3 + 2];
        push(xf, Bn, t, xi, ctx, ds, dxds_circ, dyds_circ);
        x_trace[i * 3 + 0] = xf[0];
        x_trace[i * 3 + 1] = xf[1];
        x_trace[i * 3 + 2] = xf[2];
    }
    assert(fabs(xf[0]) < 1e-8);
    assert(fabs(xf[1]-2) < 1e-8);
}

void test_8()
{
    // Test pusher around a hyperbola
    double xi[3] = {1, 0, 0};
    double xf[3];
    double s_final = M_PI;
    int nSteps = 100;
    double ds = s_final / nSteps;
    double x_trace[(nSteps + 1) * 3];
    double Bn[21] = {0};
    double t = 0;
    void *ctx = NULL;
    x_trace[0] = xi[0];
    x_trace[1] = xi[1];
    x_trace[2] = xi[2];
    for (int i = 1; i < nSteps + 1; i++)
    {
        xi[0] = x_trace[(i - 1) * 3 + 0];
        xi[1] = x_trace[(i - 1) * 3 + 1];
        xi[2] = x_trace[(i - 1) * 3 + 2];
        push(xf, Bn, t, xi, ctx, ds, dxds_hyperbola, dyds_hyperbola);
        x_trace[i * 3 + 0] = xf[0];
        x_trace[i * 3 + 1] = xf[1];
        x_trace[i * 3 + 2] = xf[2];
    }
    assert(fabs(xf[0] - 2.7401066200785826) < 1e-8);
    assert(fabs(xf[1] - 2.551114323075012) < 1e-8);
}

// // Now these tests are for the same conditions, but with trace
// void test_9()
// {
//     // Test the trace around a circle radius 1 to pi
//     double x0 = 1;
//     double y0 = 0;
//     double s_final = M_PI;
//     int nSteps = 30;
//     double xf;
//     double yf;
//     trace(&xf, &yf, dxds_circ, dyds_circ, x0, y0, s_final, nSteps);
//     assert(fabs(xf+1) < 1e-8);
//     assert(fabs(yf) < 1e-8);


// }

// void test_10()
// {
//     // Test the trace around a circle radius 1 to pi
//     // with double derivatives
//     double x0 = 1;
//     double y0 = 0;
//     double s_final = M_PI;
//     int nSteps = 30;
//     double xf;
//     double yf;
//     trace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, nSteps);
//     assert(fabs(xf+1) < 1e-8);
//     assert(fabs(yf) < 1e-8);

// }

// void test_11()
// {
//     // Test the pusher goes around a circle radius 2 to pi/2
//     double x0 = 2;
//     double y0 = 0;
//     double s_final = M_PI;
//     int nSteps = 30;
//     double xf;
//     double yf;
//     trace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, nSteps);
//     assert(fabs(xf) < 1e-8);
//     assert(fabs(yf-2) < 1e-8);
// }

// void test_12()
// {
//     // Test pusher around a hyperbola
//     double x0 = 1;
//     double y0 = 0;
//     double s_final = M_PI;
//     int nSteps = 30;
//     double xf;
//     double yf;
//     trace(&xf, &yf, dxds_hyperbola, dyds_hyperbola, x0, y0, s_final, nSteps);
//     assert(fabs(xf - 2.7401066200785826) < 1e-8);
//     assert(fabs(yf - 2.551114323075012) < 1e-8);
// }

// // Now these tests are for the same conditions, but with trace
// void test_13()
// {
//     // Test the trace around a circle radius 1 to pi
//     double x0 = 1;
//     double y0 = 0;
//     double s_final = M_PI;
//     int maxSteps = 3000;
//     double xf;
//     double yf;
//     adaptiveTrace(&xf, &yf, dxds_circ, dyds_circ, x0, y0, s_final, maxSteps);
//     assert(fabs(xf+1) < 1e-8);
//     assert(fabs(yf) < 1e-8);
// }

// void test_14()
// {
//     // Test the trace around a circle radius 1 to pi
//     // with double derivatives
//     double x0 = 1;
//     double y0 = 0;
//     double s_final = M_PI;
//     int maxSteps = 3000;
//     double xf;
//     double yf;
//     adaptiveTrace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, maxSteps);
//     assert(fabs(xf+1) < 1e-8);
//     assert(fabs(yf) < 1e-8);

// }

// void test_15()
// {
//     // Test the pusher goes around a circle radius 2 to pi/2
//     double x0 = 2;
//     double y0 = 0;
//     double s_final = M_PI;
//     int maxSteps = 3000;
//     double xf;
//     double yf;
//     adaptiveTrace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, maxSteps);
//     assert(fabs(xf) < 1e-8);
//     assert(fabs(yf-2) < 1e-8);
// }

// void test_16()
// {
//     // Test pusher around a hyperbola
//     double x0 = 1;
//     double y0 = 0;
//     double s_final = M_PI;
//     int maxSteps = 3000;
//     double xf;
//     double yf;
//     adaptiveTrace(&xf, &yf, dxds_hyperbola, dyds_hyperbola, x0, y0, s_final, maxSteps);
//     assert(fabs(xf - 2.7401066200785826) < 1e-8);
//     assert(fabs(yf - 2.551114323075012) < 1e-8);
// }

int main()
{
    test_1();
    printf("Test 1 passed\n");
    test_2();
    printf("Test 2 passed\n");
    test_3();
    printf("Test 3 passed\n");
    test_4();
    printf("Test 4 passed\n");
    test_5();
    printf("Test 5 passed\n");
    test_6();
    printf("Test 6 passed\n");
    test_7();
    printf("Test 7 passed\n");
    test_8();
    printf("Test 8 passed\n");
    // test_9();
    // printf("Test 9 passed\n");
    // test_10();
    // printf("Test 10 passed\n");
    // test_11();
    // printf("Test 11 passed\n");
    // test_12();
    // printf("Test 12 passed\n");
    // test_13();
    // printf("Test 13 passed\n");
    return 0;
}