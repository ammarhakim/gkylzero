#include <acutest.h>
#include <gkyl_mat_triples.h>
#include <gkyl_util.h>

void test_tri_1()
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(5, 5);

  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */

  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, 19.0);
  gkyl_mat_triples_insert(tri, 0, 2, 21.0);
  gkyl_mat_triples_insert(tri, 0, 3, 21.0);

  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, 21.0);
  gkyl_mat_triples_insert(tri, 1, 1, 21.0);

  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, 12.0);
  gkyl_mat_triples_insert(tri, 2, 2, 16.0);

  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, 5.0);
  gkyl_mat_triples_insert(tri, 3, 4, 21.0);

  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, 12.0);
  gkyl_mat_triples_insert(tri, 4, 1, 12.0);
  gkyl_mat_triples_insert(tri, 4, 4, 18.0);

  TEST_CHECK( gkyl_mat_triples_get(tri, 4, 4) == 18.0 );
  TEST_CHECK( gkyl_mat_triples_get(tri, 4, 0) == 12.0 );

  TEST_CHECK( gkyl_mat_triples_get(tri, 0, 1) == 0.0 ); // zero element

  TEST_CHECK( gkyl_mat_triples_size(tri) == 12 );

  gkyl_mat_triples_release(tri);
}

TEST_LIST = {
  { "tri_1", test_tri_1 },
  { NULL, NULL }
};
