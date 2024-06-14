#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_array.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#ifdef GKYL_HAVE_CUDA
# include <cuda_runtime.h>
# include <cublas_v2.h>


/**
 * @param cuh cublasHandle_t object
 * @param alpha 
*/
void cu_mat_mm_array(cublasHandle_t cuh, double alpha, double beta, enum gkyl_mat_trans transa, struct gkyl_mat *A, 
  enum gkyl_mat_trans transb, struct gkyl_array *B, struct gkyl_array *C);
  
#endif