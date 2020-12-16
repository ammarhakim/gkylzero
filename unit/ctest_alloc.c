#include <acutest.h>
#include <gkyl_alloc.h>

void
test_aligned_alloc()
{
  int *d1 = gkyl_aligned_alloc(8, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d1 % 8 == 0 );
  gkyl_aligned_free(d1);

  int *d2 = gkyl_aligned_alloc(16, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d2 % 16 == 0 );
  gkyl_aligned_free(d2);

  int *d3 = gkyl_aligned_alloc(32, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d3 % 32 == 0 );
  gkyl_aligned_free(d3);

  int *d4 = gkyl_aligned_alloc(64, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d4 % 64 == 0 );
  gkyl_aligned_free(d4);
}

void
test_aligned_realloc()
{
  int n = 10;
  int *d = gkyl_aligned_alloc(16, n*sizeof(int));

  for (int i=0; i<n; ++i)
    d[i] = 2*i;

  int *rd = gkyl_aligned_realloc(d, 16, n*sizeof(int), 2*n*sizeof(int));

  TEST_CHECK( (ptrdiff_t) rd % 16 == 0 );

  for (int i=0; i<n; ++i)
    TEST_CHECK( rd[i] == 2*i );

  gkyl_aligned_free(rd);
}

TEST_LIST = {
  { "aligned_alloc", test_aligned_alloc },
  { "aligned_realloc", test_aligned_realloc },
  { NULL, NULL },
};
