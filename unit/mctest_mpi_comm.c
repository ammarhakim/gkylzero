#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <math.h>
#include <mpi.h>
#include <stc/cstr.h>

#include <gkyl_mpi_comm.h>

void
test_1()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  TEST_CHECK( rank == m_rank );

  int sz;
  gkyl_comm_get_size(comm, &sz);
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  TEST_CHECK( sz == m_sz );

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
}

void
test_n2()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

  double vals[2];
  if (m_rank == 0) {
    vals[0] = 1.0;
    vals[1] = 2.0;
  }
  if (m_rank == 1) {
    vals[0] = 3.0;
    vals[1] = -1.0;
  }

  double v_max[2], v_min[2], v_sum[2];

  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_MAX, 2, vals, v_max);
  TEST_CHECK( v_max[0] == 3.0 );
  TEST_CHECK( v_max[1] == 2.0 );

  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_MIN, 2, vals, v_min);
  TEST_CHECK( v_min[0] == 1.0 );
  TEST_CHECK( v_min[1] == -1.0 );

  gkyl_comm_all_reduce(comm, GKYL_DOUBLE, GKYL_SUM, 2, vals, v_sum);
  TEST_CHECK( v_sum[0] == 4.0 );
  TEST_CHECK( v_sum[1] == 1.0 );

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
}

void
test_n2_sync_1d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  cstr fileNm = cstr_from_fmt("mctes_mpi_%d.txt", rank);
  FILE *fp = fopen(fileNm.str, "w");
  cstr_drop(&fileNm);
  
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 10 });

  int cuts[] = { 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = false,
    }
  );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  gkyl_print_range(&local, "local", fp);
  gkyl_print_range(&local_ext, "local_ext", fp);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);
    f[0] = iter.idx[0];
  }

  gkyl_comm_gkyl_array_sync(comm, &local, &local_ext, nghost, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_range_intersect(&in_range, &local_ext, &range);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local_ext, iter.idx);
    const double  *f = gkyl_array_cfetch(arr, idx);
    fprintf(fp, "%d: (%d). %lg\n", rank, iter.idx[0], f[0]);
    /* if (gkyl_range_contains_idx(&local, iter.idx) == 0 ) { */
    /*   if ( (iter.idx[0] != f[0]) ) */
    /*     fprintf(fp, "%d: (%d). %lg\n", rank, iter.idx[0], f[0]); */
    /* } */
  }

  fclose(fp);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void
test_n4_sync_2d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  cstr fileNm = cstr_from_fmt("mctes_mpi_%d.txt", rank);
  FILE *fp = fopen(fileNm.str, "w");
  cstr_drop(&fileNm);
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 10, 10 });

  int cuts[] = { 2, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = false,
    }
  );

  int nghost[] = { 2, 2 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  gkyl_print_range(&local, "local", fp);
  gkyl_print_range(&local_ext, "local_ext", fp);
  
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }

  gkyl_comm_gkyl_array_sync(comm, &local, &local_ext, nghost, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_range_intersect(&in_range, &local_ext, &range);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local_ext, iter.idx);
    const double  *f = gkyl_array_cfetch(arr, idx);

    if (gkyl_range_contains_idx(&local, iter.idx) == 0 ) {
      if ( (iter.idx[0] != f[0]) || (iter.idx[1] != f[1]) )
        fprintf(fp, "%d: (%d, %d). %lg %lg\n", rank, iter.idx[0], iter.idx[1], f[0], f[1]);
    }
  }

  fclose(fp);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

TEST_LIST = {
    {"test_1", test_1},
    {"test_n2", test_n2},
    {"test_n2_sync_1d", test_n2_sync_1d},
    {"test_n4_sync_2d", test_n4_sync_2d},
    {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
    {NULL, NULL},
};

#endif
