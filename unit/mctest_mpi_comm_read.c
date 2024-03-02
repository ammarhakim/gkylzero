#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <math.h>
#include <mpi.h>
#include <stc/cstr.h>

#include <gkyl_array_rio.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

void
mpi_n2_read()
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != 2) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct gkyl_rect_grid grid;
  struct gkyl_array_header_info hdr;
  gkyl_grid_sub_array_header_read(&grid, &hdr,
    "data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];  
  
  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  int cuts[] = { 2, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);
  
  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);  
}

TEST_LIST = {
  {"mpi_n2_read", mpi_n2_read},
  {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  {NULL, NULL},
};

#endif
