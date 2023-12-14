#include "hamilton.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

// Build using gcc test_hamilton.c hamilton.c -o test_hamilton -lm && ./test_hamilton

// Define a function pointer type for the derivative function
static const double h[8] = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};

typedef double (*derivative_func)(double, double);

double dyds_circ(double x, double y)
{
    return x;
}

double dxds_circ(double x, double y)
{
    return -y;
}

void test_1()
{
    // Test calculateGFromF function
    double F_test[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    double G_test[7];
    calculateGFromF(G_test, F_test);

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
    calculateB(B_test, G_test);

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
    // Test the pusher goes around a circle radius 1 to pi
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double ds = s_final / nSteps;
    double x[nSteps + 1];
    double y[nSteps + 1];
    x[0] = x0;
    y[0] = y0;
    double Bx[7] = {0};
    double By[7] = {0};
    for (int i = 1; i < nSteps + 1; i++)
    {
        push(&x[i], &y[i], Bx, By, x[i - 1], y[i - 1], ds, dxds_circ, dyds_circ);
    }

    assert(fabs(x[nSteps] + 1) < 1e-8);
    assert(fabs(y[nSteps]) < 1e-8);
}

double dyds_double_circ(double x, double y)
{
    return 2 * x;
}

double dxds_double_circ(double x, double y)
{
    return -2 * y;
}

void test_4()
{
    // Test the pusher goes around a circle radius 1 to pi
    // with double derivatives
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double ds = s_final / nSteps;
    double x[nSteps + 1];
    double y[nSteps + 1];
    x[0] = x0;
    y[0] = y0;
    double Bx[7] = {0};
    double By[7] = {0};
    for (int i = 1; i <= nSteps; i++)
    {
        push(&x[i], &y[i], Bx, By, x[i - 1], y[i - 1], ds, dxds_double_circ, dyds_double_circ);
    }
    assert(fabs(x[nSteps] + 1) < 1e-8);
    assert(fabs(y[nSteps]) < 1e-8);
}

void test_5()
{
    // Test the pusher goes around a circle radius 2 to pi/2
    double x0 = 2;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double ds = s_final / nSteps;
    double x[nSteps + 1];
    double y[nSteps + 1];
    x[0] = x0;
    y[0] = y0;
    double Bx[7] = {0};
    double By[7] = {0};
    for (int i = 1; i <= nSteps; i++)
    {
        push(&x[i], &y[i], Bx, By, x[i - 1], y[i - 1], ds, dxds_circ, dyds_circ);
    }
    assert(fabs(x[nSteps]) < 1e-8);
    assert(fabs(y[nSteps] - 2) < 1e-8);
}

double dyds_hyperbola(double x, double y)
{
    return x;
}

double dxds_hyperbola(double x, double y)
{
    return y;
}

void test_6()
{
    // Test pusher around a hyperbola
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double ds = s_final / nSteps;
    double x[nSteps + 1];
    double y[nSteps + 1];
    x[0] = x0;
    y[0] = y0;
    double Bx[7] = {0};
    double By[7] = {0};
    for (int i = 1; i <= nSteps; i++)
    {
        push(&x[i], &y[i], Bx, By, x[i - 1], y[i - 1], ds, dxds_hyperbola, dyds_hyperbola);
    }
    assert(fabs(x[nSteps] - 2.7401066200785826) < 1e-8);
    assert(fabs(y[nSteps] - 2.551114323075012) < 1e-8);
}

// Now these tests are for the same conditions, but with trace
void test_7()
{
    // Test the trace around a circle radius 1 to pi
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double xf;
    double yf;
    trace(&xf, &yf, dxds_circ, dyds_circ, x0, y0, s_final, nSteps);
    assert(fabs(xf+1) < 1e-8);
    assert(fabs(yf) < 1e-8);
}

void test_8()
{
    // Test the trace around a circle radius 1 to pi
    // with double derivatives
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double xf;
    double yf;
    trace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, nSteps);
    assert(fabs(xf+1) < 1e-8);
    assert(fabs(yf) < 1e-8);
   
}

void test_9()
{
    // Test the pusher goes around a circle radius 2 to pi/2
    double x0 = 2;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double xf;
    double yf;
    trace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, nSteps);
    assert(fabs(xf) < 1e-8);
    assert(fabs(yf-2) < 1e-8);
}

void test_10()
{
    // Test pusher around a hyperbola
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int nSteps = 30;
    double xf;
    double yf;
    trace(&xf, &yf, dxds_hyperbola, dyds_hyperbola, x0, y0, s_final, nSteps);
    assert(fabs(xf - 2.7401066200785826) < 1e-8);
    assert(fabs(yf - 2.551114323075012) < 1e-8);
}

// Now these tests are for the same conditions, but with trace
void test_11()
{
    // Test the trace around a circle radius 1 to pi
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int maxSteps = 3000;
    double xf;
    double yf;
    adaptiveTrace(&xf, &yf, dxds_circ, dyds_circ, x0, y0, s_final, maxSteps);
    assert(fabs(xf+1) < 1e-8);
    assert(fabs(yf) < 1e-8);
}

void test_12()
{
    // Test the trace around a circle radius 1 to pi
    // with double derivatives
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int maxSteps = 3000;
    double xf;
    double yf;
    adaptiveTrace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, maxSteps);
    assert(fabs(xf+1) < 1e-8);
    assert(fabs(yf) < 1e-8);
   
}

void test_13()
{
    // Test the pusher goes around a circle radius 2 to pi/2
    double x0 = 2;
    double y0 = 0;
    double s_final = M_PI;
    int maxSteps = 3000;
    double xf;
    double yf;
    adaptiveTrace(&xf, &yf, dxds_double_circ, dyds_double_circ, x0, y0, s_final, maxSteps);
    assert(fabs(xf) < 1e-8);
    assert(fabs(yf-2) < 1e-8);
}

void test_14()
{
    // Test pusher around a hyperbola
    double x0 = 1;
    double y0 = 0;
    double s_final = M_PI;
    int maxSteps = 3000;
    double xf;
    double yf;
    adaptiveTrace(&xf, &yf, dxds_hyperbola, dyds_hyperbola, x0, y0, s_final, maxSteps);
    assert(fabs(xf - 2.7401066200785826) < 1e-8);
    assert(fabs(yf - 2.551114323075012) < 1e-8);
}

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
    test_9();
    printf("Test 9 passed\n");
    test_10();
    printf("Test 10 passed\n");
    test_11();
    printf("Test 11 passed\n");
    test_12();
    printf("Test 12 passed\n");
    test_13();
    printf("Test 13 passed\n");
    return 0;
}