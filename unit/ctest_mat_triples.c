#include <acutest.h>
#include <gkyl_mat_triples.h>
#include <gkyl_util.h>

void test_tri_1()
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(5, 5);

  gkyl_mat_triples_release(tri);
}

TEST_LIST = {
  { "tri_1", test_tri_1 },
  { NULL, NULL }
};
