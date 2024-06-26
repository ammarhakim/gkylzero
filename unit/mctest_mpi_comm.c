#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <math.h>
#include <mpi.h>
#include <stc/cstr.h>
#include <gkyl_util.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_mpi_comm.h>

void
mpi_1()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });

  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  
  int cuts[] = { m_sz, 1 };
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
  TEST_CHECK( sz == m_sz );

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
}

void
mpi_allreduce()
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
mpi_n2_sync_1d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);
    f[0] = iter.idx[0];
  }

  gkyl_comm_array_sync(comm, &local, &local_ext, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &range);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double  *f = gkyl_array_cfetch(arr, idx);
    
    TEST_CHECK( iter.idx[0] == f[0] );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void
mpi_n4_sync_2d(bool use_corners)
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 10, 10 });

  int cuts[] = { 2, 2 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);  
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = use_corners,
    }
  );

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x, local_ext_x, local_y, local_ext_y;
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) {1, 0},
    &local_ext_x, &local_x);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, 1 },
    &local_ext_y, &local_y);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  gkyl_comm_barrier(comm);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }

  gkyl_comm_array_sync(comm, &local, &local_ext, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &range);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double  *f = gkyl_array_cfetch(arr, idx);

    if (use_corners) {
      TEST_CHECK( iter.idx[0] == f[0] );
      TEST_CHECK( iter.idx[1] == f[1] );
    }
    else {
      // excludes corners
      if (gkyl_range_contains_idx(&local_ext_x, iter.idx) || gkyl_range_contains_idx(&local_ext_y, iter.idx)) {
        TEST_CHECK( iter.idx[0] == f[0] );
        TEST_CHECK( iter.idx[1] == f[1] );
      }
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void mpi_n4_sync_2d_no_corner() { mpi_n4_sync_2d(false); }
void mpi_n4_sync_2d_use_corner() { mpi_n4_sync_2d(true); }

void
mpi_n4_sync_1x1v()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 512 });

  int cuts[] = { m_sz };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(1, cuts, &range);

  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
    }
  );

  struct gkyl_range vrange;
  gkyl_range_init(&vrange, 1, (int[]) { 1 }, (int[]) { 64 } );

  struct gkyl_rect_decomp *ext_decomp =
    gkyl_rect_decomp_extended_new(&vrange, decomp);

  int nghost[] = { 1, 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&ext_decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_comm *ext_comm = gkyl_comm_extend_comm(comm, &vrange);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  gkyl_array_clear(arr, 200005);
  gkyl_comm_barrier(comm);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }

  gkyl_comm_array_sync(ext_comm, &local, &local_ext, arr);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &ext_decomp->parent_range);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double  *f = gkyl_array_cfetch(arr, idx);

    TEST_CHECK( iter.idx[0] == f[0] );
    TEST_CHECK( iter.idx[1] == f[1] );
  }

  gkyl_array_release(arr);
  gkyl_comm_release(ext_comm);
  gkyl_rect_decomp_release(ext_decomp);
  gkyl_comm_release(comm);
  gkyl_rect_decomp_release(decomp);
}

void
mpi_n1_per_sync_2d_tests(int num_per_dirs, int *per_dirs)
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 1) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .decomp = decomp,
      .mpi_comm = MPI_COMM_WORLD,
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);

    for (int d=0; d<local.ndim; ++d)
      f[d] = iter.idx[d];
  }

  gkyl_comm_array_per_sync(comm, &local, &local_ext, num_per_dirs, per_dirs, arr );

  int idx[GKYL_MAX_DIM] = { 0 };
  int count = 0;
  
  for (int id=0; id<num_per_dirs; ++id) {
    int d = per_dirs[id];
    int ncell = gkyl_range_shape(&local, d);

    gkyl_range_iter_init(&iter, &local_ext_x[d]);
    while (gkyl_range_iter_next(&iter)) {

      if (!gkyl_range_contains_idx(&local, iter.idx)) {
        long lidx = gkyl_range_idx(&local_ext, iter.idx);
        
        for (int n=0; n<local.ndim; ++n)
          idx[n] = iter.idx[n];
        if (idx[d] > local.upper[d])
          idx[d] = idx[d] - ncell;
        else
          idx[d] = idx[d] + ncell;

        const double  *f = gkyl_array_cfetch(arr, lidx);
        for (int n=0; n<local.ndim; ++n) {
          TEST_CHECK( idx[n] == f[n] );
          TEST_MSG( "rank:%d | At idx=(%d,%d) | Expected: %d | Produced: %.13e", rank, iter.idx[0], iter.idx[1], idx[n], f[n] );
	}
      }
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void
mpi_n1_per_sync_2d()
{
  int per_dirs_0[] = {0};
  int per_dirs_1[] = {1};
  int per_dirs_01[] = {0,1};

  mpi_n1_per_sync_2d_tests(1, per_dirs_0);
  mpi_n1_per_sync_2d_tests(1, per_dirs_1);
  mpi_n1_per_sync_2d_tests(2, per_dirs_01);

  mpi_n1_per_sync_2d_tests(1, per_dirs_0);
  mpi_n1_per_sync_2d_tests(1, per_dirs_1);
  mpi_n1_per_sync_2d_tests(2, per_dirs_01);
}

void
mpi_n2_per_sync_2d_tests(int *cuts, int num_per_dirs, int *per_dirs)
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;
  
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .decomp = decomp,
      .mpi_comm = MPI_COMM_WORLD,
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  // Redefine local_ext_x so it's local shifted in the right direction
  // so it covers the ghost cells of interest.
  int decomp_dir = cuts[0]>1? 0 : 1;
  int delta[] = {0, 0};
  delta[decomp_dir] = 2*rank-1;
  struct gkyl_range local_ext_x_shifted;
  gkyl_range_shift(&local_ext_x_shifted, &local_ext_x[decomp_dir], delta);
  gkyl_sub_range_init(&local_ext_x[decomp_dir], &local_ext, local_ext_x_shifted.lower, local_ext_x_shifted.upper);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);

    for (int d=0; d<local.ndim; ++d)
      f[d] = iter.idx[d];
  }

  gkyl_comm_array_per_sync(comm, &local, &local_ext, num_per_dirs, per_dirs, arr );

  int idx[GKYL_MAX_DIM] = { 0 };
  int count = 0;
  
  for (int id=0; id<num_per_dirs; ++id) {
    int d = per_dirs[id];
    int ncell = gkyl_range_shape(&range, d);

    gkyl_range_iter_init(&iter, &local_ext_x[d]);
    while (gkyl_range_iter_next(&iter)) {

      if (!gkyl_range_contains_idx(&local, iter.idx)) {
        long lidx = gkyl_range_idx(&local_ext, iter.idx);
        
        for (int n=0; n<local.ndim; ++n) idx[n] = iter.idx[n];

        if (idx[d] > local.upper[d])
          idx[d] = idx[d] - ncell;
        else if (idx[d] < local.lower[d])
          idx[d] = idx[d] + ncell;

        const double  *f = gkyl_array_cfetch(arr, lidx);
        for (int n=0; n<local.ndim; ++n) {
          TEST_CHECK( idx[n] == f[n] );
          TEST_MSG( "rank:%d | At idx=(%d,%d) | Expected: %d | Produced: %.13e", rank, iter.idx[0], iter.idx[1], idx[n], f[n] );
	}
      }
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void
mpi_n2_per_sync_2d()
{
  int cuts_21[] = {2,1};
  int cuts_12[] = {1,2};
  int per_dirs_0[] = {0};
  int per_dirs_1[] = {1};
  int per_dirs_01[] = {0,1};

  mpi_n2_per_sync_2d_tests(cuts_21, 1, per_dirs_0);
  mpi_n2_per_sync_2d_tests(cuts_21, 1, per_dirs_1);
  mpi_n2_per_sync_2d_tests(cuts_21, 2, per_dirs_01);

  mpi_n2_per_sync_2d_tests(cuts_12, 1, per_dirs_0);
  mpi_n2_per_sync_2d_tests(cuts_12, 1, per_dirs_1);
  mpi_n2_per_sync_2d_tests(cuts_12, 2, per_dirs_01);
}

void
mpi_n2_array_send_irecv_1d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  double sendval = rank==0? 20005. : 30005.;
  double recvval = rank==0? 30005. : 20005.;

  // Assume the range is not decomposed. 
  struct gkyl_array *arrA = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  struct gkyl_array *arrB = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_clear(arrA, sendval*(1-rank));
  gkyl_array_clear(arrB, sendval*rank);

  struct gkyl_array *recvbuff = rank==0? arrB : arrA;
  struct gkyl_array *sendbuff = rank==0? arrA : arrB;

  struct gkyl_comm_state *cstate = gkyl_comm_state_new(comm);
  int tag = 13;
  // Post irecv before send.
  gkyl_comm_array_irecv(comm, recvbuff, (rank+1) % 2, tag, cstate);
  gkyl_comm_array_send(comm, sendbuff, (rank+1) % 2, tag);
  gkyl_comm_state_wait(comm, cstate);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&range, iter.idx);
    const double *f = gkyl_array_cfetch(recvbuff, idx);
    TEST_CHECK( f[0] == recvval );
  }

  gkyl_comm_barrier(comm);

  gkyl_comm_state_release(comm, cstate);
  gkyl_array_release(arrA);
  gkyl_array_release(arrB);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
}

void
mpi_n2_array_isend_irecv_2d()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 10, 20 });

  int cuts[] = { 1, m_sz };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp,
      .sync_corners = false,
    }
  );

  double sendval = rank==0? 20005. : 30005.;
  double recvval = rank==0? 30005. : 20005.;

  // Assume the range is not decomposed. 
  struct gkyl_array *arrA = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  struct gkyl_array *arrB = gkyl_array_new(GKYL_DOUBLE, 1, range.volume);
  gkyl_array_clear(arrA, sendval*(1-rank));
  gkyl_array_clear(arrB, sendval*rank);

  struct gkyl_array *recvbuff = rank==0? arrB : arrA;
  struct gkyl_array *sendbuff = rank==0? arrA : arrB;

  struct gkyl_comm_state *cstate_s = gkyl_comm_state_new(comm);
  struct gkyl_comm_state *cstate_r = gkyl_comm_state_new(comm);
  int tag = 13;
  // Post irecv before send.
  gkyl_comm_array_irecv(comm, recvbuff, (rank+1) % 2, tag, cstate_r);
  gkyl_comm_array_isend(comm, sendbuff, (rank+1) % 2, tag, cstate_s);

  // Do some other unnecessary work.
  struct gkyl_array *tmp_arr = gkyl_array_new(GKYL_DOUBLE, 3, range.volume);
  gkyl_array_clear(tmp_arr, 13.);

  gkyl_comm_state_wait(comm, cstate_r);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&range, iter.idx);
    const double *f = gkyl_array_cfetch(recvbuff, idx);
    TEST_CHECK( f[0] == recvval );
  }

  gkyl_comm_state_wait(comm, cstate_s);

  gkyl_comm_state_release(comm, cstate_s);
  gkyl_comm_state_release(comm, cstate_r);
  gkyl_array_release(arrA);
  gkyl_array_release(arrB);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
}

void
mpi_n4_multicomm_2d()
{
  // Test the use of two gkyl_comm objects simultaneously, mimicing the case
  // where one is used to decompose space and the other species.
  // We sync across the conf communicator, and send/recv across the species
  // comm.
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 4) return;

  struct gkyl_comm *worldcomm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = 0,
    }
  );

  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 10, 20 });

  int confcuts[] = { 2, 1 };
  struct gkyl_rect_decomp *confdecomp = gkyl_rect_decomp_new_from_cuts(2, confcuts, &range);  
  
  int worldrank;
  gkyl_comm_get_rank(worldcomm, &worldrank);

  int confcolor = floor(worldrank/confdecomp->ndecomp);
  struct gkyl_comm *confcomm = gkyl_comm_split_comm(worldcomm, confcolor, confdecomp);
  int confrank;
  gkyl_comm_get_rank(confcomm, &confrank);

  int speciescolor = worldrank % confdecomp->ndecomp;
  struct gkyl_comm *speciescomm = gkyl_comm_split_comm(worldcomm, speciescolor, confdecomp);
  int speciesrank;
  gkyl_comm_get_rank(speciescomm, &speciesrank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&confdecomp->ranges[confrank], nghost, &local_ext, &local);

  struct gkyl_array *arrA = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  struct gkyl_array *arrB = gkyl_array_new(GKYL_DOUBLE, 2, local_ext.volume);
  gkyl_array_clear(arrA, 0.);
  gkyl_array_clear(arrB, 0.);
  struct gkyl_array *recvbuff = speciesrank==0? arrB : arrA;
  struct gkyl_array *sendbuff = speciesrank==0? arrA : arrB;

  gkyl_comm_barrier(worldcomm);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double *f = gkyl_array_fetch(sendbuff, idx);
    f[0] = iter.idx[0]; f[1] = iter.idx[1];
  }

  // Sync sendbuff array and check results.
  gkyl_comm_array_sync(confcomm, &local, &local_ext, sendbuff);

  struct gkyl_range in_range; // interior, including ghost cells
  gkyl_sub_range_intersect(&in_range, &local_ext, &range);
  struct gkyl_range local_x, local_ext_x, local_y, local_ext_y;
  gkyl_create_ranges(&confdecomp->ranges[confrank], (int[]) {1, 0}, &local_ext_x, &local_x);
  gkyl_create_ranges(&confdecomp->ranges[confrank], (int[]) { 0, 1 }, &local_ext_y, &local_y);
  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double *f = gkyl_array_cfetch(sendbuff, idx);
    // exclude corners
    if (gkyl_range_contains_idx(&local_ext_x, iter.idx) || gkyl_range_contains_idx(&local_ext_y, iter.idx)) {
      TEST_CHECK( iter.idx[0] == f[0] );
      TEST_CHECK( iter.idx[1] == f[1] );
    }
  }

  // Now send/recv across species communicator and check results.
  struct gkyl_comm_state *cstate = gkyl_comm_state_new(speciescomm);
  int tag = 13;
  // Post irecv before send.
  gkyl_comm_array_irecv(speciescomm, recvbuff, (speciesrank+1) % 2, tag, cstate);
  gkyl_comm_array_send(speciescomm, sendbuff, (speciesrank+1) % 2, tag);
  gkyl_comm_state_wait(speciescomm, cstate);

  gkyl_range_iter_init(&iter, &in_range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&in_range, iter.idx);
    const double *f = gkyl_array_cfetch(recvbuff, idx);
    // exclude corners
    if (gkyl_range_contains_idx(&local_ext_x, iter.idx) || gkyl_range_contains_idx(&local_ext_y, iter.idx)) {
      TEST_CHECK( iter.idx[0] == f[0] );
      TEST_CHECK( iter.idx[1] == f[1] );
    }
  }

  gkyl_comm_state_release(speciescomm, cstate);
  gkyl_array_release(arrA);
  gkyl_array_release(arrB);
  gkyl_comm_release(speciescomm);
  gkyl_comm_release(confcomm);
  gkyl_rect_decomp_release(confdecomp);
  gkyl_comm_release(worldcomm);
}
  
TEST_LIST = {
  {"mpi_1", mpi_1},
  {"mpi_allreduce", mpi_allreduce},
  {"mpi_n2_sync_1d", mpi_n2_sync_1d},
  {"mpi_n4_sync_2d_no_corner", mpi_n4_sync_2d_no_corner },
  {"mpi_n4_sync_2d_use_corner", mpi_n4_sync_2d_use_corner},
  {"mpi_n2_sync_1x1v", mpi_n4_sync_1x1v },
  {"mpi_n1_per_sync_2d", mpi_n1_per_sync_2d },
  {"mpi_n2_per_sync_2d", mpi_n2_per_sync_2d },
  {"mpi_n2_array_send_irecv_1d", mpi_n2_array_send_irecv_1d},
  {"mpi_n2_array_isend_irecv_2d", mpi_n2_array_isend_irecv_2d},
  {"mpi_n4_multicomm_2d", mpi_n4_multicomm_2d},
  {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  {NULL, NULL},
};

#endif
