#include <acutest.h>

#include <gkyl_rrobin_decomp.h>

static void
test_1(void)
{
  const struct gkyl_rrobin_decomp *rr = gkyl_rrobin_decomp_new(4, 3,
    (int[]) { 4, 1, 1 });

  TEST_CHECK( 4 == rr->total_ranks );
  TEST_CHECK( 3 == rr->nblocks );

  gkyl_rrobin_decomp_release(rr);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
