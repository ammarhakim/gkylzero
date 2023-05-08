#include <acutest.h>
#include <gkyl_null_comm.h>

void
test_1()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 10 });

  struct gkyl_rect_decomp *rdecomp =
    gkyl_rect_decomp_new_from_cuts(1, (int[]) { 1 }, &range);
  
  struct gkyl_comm *comm = gkyl_null_comm_new( &(struct gkyl_null_comm_inp) {
      .decomp = rdecomp
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  TEST_CHECK( rank == 0 );

  int sz;
  gkyl_comm_get_size(comm, &sz);
  TEST_CHECK( sz == 1 );  

  double out[3], inp[3] = { 2.0, 4.0, 8.0 };
  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_MIN, 3, inp, out);

  for (int i=0; i<3; ++i)
    TEST_CHECK( out[i] == inp[i] );

  gkyl_rect_decomp_release(rdecomp);
  gkyl_comm_release(comm);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
