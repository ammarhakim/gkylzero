#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gkyl_array.h>
#include "bilinear_interp.h"

/*
 * bilinear_interp.c
 *
 *  Created on: June 14, 2023
 *      Author: Jonathan Roeltgen
 */

/* Function bilinear_interp - perform bilinear interpolation
 * on a function with given x, y, z
 * Inputs:
 * 		double* x - array of x abscissa
 * 		double* y - array of y abscissa
 * 		double* z - 2D array of ordinates given as (x,y)
 * 		double* xq - array of query points in x direction
 * 		double* yq - array of query points in y direction
 * 		int nx - length of x array
 * 		int ny - length of y array
 * 		int nq - number of query points
 *
 * Return:
 * 		double* - 2D array of function values interpolated at
 * 			the given points (xq,yq)
 *
 * Note: Need to consider the cases that (xq,yq) is outside of
 * 		domain given by data
 */
void bilinear_interp(double* x, double* y, double* z, const struct gkyl_array* xq,
  const struct gkyl_array* yq, int nx, int ny, int nq, struct gkyl_array* zq){
	int j,k,i11,i12,i21,i22;
	double x2x1, y2y1, x2x, y2y, yy1, xx1;

	/* zq = (double*)calloc(nq,sizeof(double)); */
	/* if(zq == NULL){ */
	/* 	printf("Memory allocation in bilinear_interp failed"); */
	/* 	exit(EXIT_FAILURE); */
	/* } */
	
	//loop through all query points
	for(int i=0; i<nq; i++){
		j = -1;
		k = -1;
		const double *xq_d = gkyl_array_cfetch(xq, i);
		const double *yq_d = gkyl_array_cfetch(yq, i);
		double *zq_d = gkyl_array_fetch(zq, i);
		while(xq_d[0]>=x[j+1]&&j<nx-2){
			j++;
		}
		while(yq_d[0]>=y[k+1]&&k<ny-2){
			k++;
		}
		i11 = j*ny + k;
		i21 = (j + 1)*ny + k;
		i12 = j*ny + k + 1;
		i22 = (j + 1)*ny + k + 1;

	    x2x1 = x[j+1] - x[j];
	    y2y1 = y[k+1] - y[k];
	    x2x = x[j+1] - xq_d[0];
	    y2y = y[k+1] - yq_d[0];
	    yy1 = yq_d[0] - y[k];
	    xx1 = xq_d[0] - x[j];
		//4 points of surrounding square are: x[j],x[j+1],y[k],y[k+1]
		//zq[i] = 1/((x[j+1]-x[j])*(y[k+1]-y[k])) * (z[i11]*(x[j+1]-xq[i])*(y[k+1]-yq[i]) +
		//		z[i21]*(xq[i]-x[j])*(y[k+1]-yq[i]) + z[i12]*(x[j+1]-xq[i])*(yq[i]-y[k]) +
		//		z[i22]*(xq[i]-x[j])*(yq[i]-y[k]));
	    zq_d[0] = pow(10., 1.0 / (x2x1 * y2y1) * (
	        z[i11] * x2x * y2y +
	        z[i21] * xx1 * y2y +
	        z[i12] * x2x * yy1 +
	        z[i22] * xx1 * yy1));
	}
}
