#include <acutest.h>
#include <gkyl_mat_triples.h>
#include <gkyl_util.h>
#include <gkyl_range.h>

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

void test_tri_2()
{
  gkyl_mat_triples *tri = gkyl_mat_triples_new(3, 3);
  /*  A : matrix([1,0,0],[0,2,0],[0,0,3]); */

  // row 1
  gkyl_mat_triples_insert(tri, 1, 1, 2.0);
  gkyl_mat_triples_insert(tri, 1, 0, 2.1);

  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, 1.0);
  gkyl_mat_triples_insert(tri, 0, 2, 1.1);

  // row 2
  gkyl_mat_triples_insert(tri, 2, 2, 3.0);

  // Loop through keys, obtain coordinates for each, the value,
  // and check.
  long *keys = gkyl_mat_triples_keys(tri);
  for (size_t i=0; i<gkyl_mat_triples_size(tri); i++) {
    int idx[2];
    gkyl_mat_triples_key_to_idx(tri, keys[i], idx);
    TEST_CHECK( gkyl_mat_triples_get(tri, idx[0], idx[1]) == gkyl_mat_triples_val_at_key(tri, keys[i]) );
  }

  // Do it again but in column-major order (colmo).
  long *skeys = gkyl_mat_triples_keys_colmo(tri);
  for (size_t i=0; i<gkyl_mat_triples_size(tri); i++) {
    int idx[2];
    gkyl_mat_triples_key_to_idx(tri, skeys[i], idx);
    TEST_CHECK( gkyl_mat_triples_get(tri, idx[0], idx[1]) == gkyl_mat_triples_val_at_key(tri, skeys[i]) );
  }

  gkyl_mat_triples_release(tri);
  gkyl_free(keys);
  gkyl_free(skeys);
}

TEST_LIST = {
  { "tri_1", test_tri_1 },
  { "tri_2", test_tri_2 },
  { NULL, NULL }
};
