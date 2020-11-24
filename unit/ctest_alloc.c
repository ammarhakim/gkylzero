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

TEST_LIST = {
  { "aligned_alloc", test_aligned_alloc },
  { NULL, NULL },
};
