#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_null_comm.h>

#include <string.h>

// Private struct
struct null_comm {
  struct gkyl_comm base; // base communicator
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition
};

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);  
  struct null_comm *null_comm = container_of(comm, struct null_comm, base);
  gkyl_rect_decomp_release(null_comm->decomp);
  gkyl_free(null_comm);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  *rank = 0;
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  *sz = 1;
  return 0;
}

static int
all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  memcpy(out, inp, gkyl_elem_type_size[type]*nelem);
  return 0;
}

static int
array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  const int *nghost, struct gkyl_array *array)
{
  return 0;
}

static int barrier(struct gkyl_comm *comm) { return 0; }

static int
array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  return gkyl_grid_sub_array_write(grid, range, arr, fname);
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, base);
  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, null_comm->decomp);
  struct gkyl_comm *ext_comm = gkyl_null_comm_new( &(struct gkyl_null_comm_inp) {
      .decomp = ext_decomp
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  
  return ext_comm;
}

struct gkyl_comm*
gkyl_null_comm_new(const struct gkyl_null_comm_inp *inp)
{
  struct null_comm *comm = gkyl_malloc(sizeof *comm);
  comm->decomp = gkyl_rect_decomp_acquire(inp->decomp);

  comm->base.get_rank = get_rank;
  comm->base.get_size = get_size;
  comm->base.all_reduce = all_reduce;
  comm->base.gkyl_array_sync = array_sync;
  comm->base.barrier = barrier;
  comm->base.gkyl_array_write = array_write;
  comm->base.extend_comm = extend_comm;

  comm->base.ref_count = gkyl_ref_count_init(comm_free);

  return &comm->base;
}
