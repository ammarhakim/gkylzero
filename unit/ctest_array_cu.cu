/* -*- c -*- */

#include <gkylzero.h>
#include <stdio.h>

extern "C" {
    int cu_array_test_and_flip_sign( struct gkyl_array *arr);
}

__global__
void ker_cu_array_test_and_flip_sign( struct gkyl_array *arr, int *nfail)
{
  *nfail = 0;
  
  GKYL_CU_CHECK( arr->type == GKYL_DOUBLE, nfail );
  GKYL_CU_CHECK( arr->elemsz ==sizeof(double), nfail );
  GKYL_CU_CHECK( arr->ncomp == 1, nfail );
  GKYL_CU_CHECK( arr->size == 20, nfail );
  
   double *data = ( double *) arr->data;
  for (unsigned i=0; i<arr->size; ++i) {
    GKYL_CU_CHECK( data[i] == (i+0.5)*0.1, nfail );
    data[i] *= -1;
  }

}

int cu_array_test_and_flip_sign( struct gkyl_array *arr)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_array_test_and_flip_sign<<<1,1>>>(arr, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
