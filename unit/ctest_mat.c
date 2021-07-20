#include <acutest.h>
#include <gkyl_mat.h>

void
test_mat_base()
{
  struct gkyl_mat *m = gkyl_mat_new(10, 20, 0.25);

  TEST_CHECK( 10 == m->nr );
  TEST_CHECK( 20 == m->nc );

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK ( 0.25 == gkyl_mat_get(m, i, j) );

  gkyl_mat_clear(m, 0.1);

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK ( 0.1 == gkyl_mat_get(m, i, j) );

  size_t count = 0;
  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      gkyl_mat_set(m, i, j, count++);

  count = 0;
  for (size_t j=0; j<m->nc; ++j) {
    const double* col = gkyl_mat_get_ccol(m, j);
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( col[i] == count++ );
  }

  count = 0;
  for (size_t i=0; i<m->nr*m->nc; ++i)
    TEST_CHECK( m->data[i] == i );

  gkyl_mat_diag(m, 1.0);
  
  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( gkyl_mat_get(m, i, j) == ( i==j ? 1.0 : 0.0 ) );

  struct gkyl_mat *m2 = gkyl_mat_clone(m);

  TEST_CHECK( m2->nr == m->nr );
  TEST_CHECK( m2->nc == m->nc );

  for (size_t j=0; j<m->nc; ++j)
    for (size_t i=0; i<m->nr; ++i)
      TEST_CHECK( gkyl_mat_get(m, i, j) == gkyl_mat_get(m2, i, j) );

  gkyl_mat_release(m);
  gkyl_mat_release(m2);
}

TEST_LIST = {
  { "mat_base", test_mat_base },
  { NULL, NULL },
};
