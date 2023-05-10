#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_null_comm.h>

#include <string.h>
#include <math.h>

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  long max_vol; // maximum vol of send/recv region
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
#define G_MAX(a,b) (a)>(b)?(a):(b)
  
  int ndim = parent->ndim;
  long max_vol = 0;
    
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);
    
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

// define long -> skin_ghost_ranges ...
#define i_key long
#define i_val struct skin_ghost_ranges
#define i_tag l2sgr
#include <stc/cmap.h>
// ... done with map definition

// Private struct
struct null_comm {
  struct gkyl_comm base; // base communicator
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition

  struct gkyl_range grange; // range to "hash" ghost layout
  cmap_l2sgr l2sgr; // map from int -> skin_ghost_ranges

  gkyl_mem_buff pbuff; // buffer for periodic BCs
};

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);  
  struct null_comm *null_comm = container_of(comm, struct null_comm, base);

  cmap_l2sgr_drop(&null_comm->l2sgr);
  gkyl_rect_decomp_release(null_comm->decomp);
  gkyl_mem_buff_release(null_comm->pbuff);
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
  struct gkyl_array *array)
{
  return 0;
}

// apply periodic BCs
static void
apply_periodic_bc(struct skin_ghost_ranges *sgr, char *data,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(data, f, sgr->lower_skin[dir]);
  gkyl_array_copy_from_buffer(f, data, sgr->upper_ghost[dir]);

  gkyl_array_copy_to_buffer(data, f, sgr->upper_skin[dir]);
  gkyl_array_copy_from_buffer(f, data, sgr->lower_ghost[dir]);
}

static int
array_per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, base);

  int nghost[GKYL_MAX_DIM];
  for (int d=0; d<null_comm->decomp->ndim; ++d)
    nghost[d] = local_ext->upper[d]-local->upper[d];
  
  long lkey = gkyl_range_idx(&null_comm->grange, nghost);

  if (!cmap_l2sgr_contains(&null_comm->l2sgr, lkey)) {
    struct skin_ghost_ranges sgr;
    skin_ghost_ranges_init(&sgr, local_ext, nghost);
    cmap_l2sgr_insert(&null_comm->l2sgr, lkey, sgr);
  }

  cmap_l2sgr_value *val = cmap_l2sgr_get_mut(&null_comm->l2sgr, lkey);
  long max_vol_esnz = val->second.max_vol*array->esznc;
  
  if  (max_vol_esnz > gkyl_mem_buff_size(null_comm->pbuff))
    null_comm->pbuff = gkyl_mem_buff_resize(null_comm->pbuff, max_vol_esnz);

  char *data = gkyl_mem_buff_data(null_comm->pbuff);

  for (int d=0; d<nper_dirs; ++d)
    apply_periodic_bc(&val->second, data, per_dirs[d], array);

  return 0;
}

static int
barrier(struct gkyl_comm *comm)
{
  return 0;
}

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

  // construct range to hash ghost layout
  int lower[GKYL_MAX_DIM] = { 0 };
  int upper[GKYL_MAX_DIM];
  for (int d=0; d<inp->decomp->ndim; ++d) upper[d] = GKYL_MAX_NGHOST;
  gkyl_range_init(&comm->grange, inp->decomp->ndim, lower, upper);

  comm->l2sgr = cmap_l2sgr_init();
  comm->pbuff = gkyl_mem_buff_new(1024); // will be reallocated

  comm->base.get_rank = get_rank;
  comm->base.get_size = get_size;
  comm->base.all_reduce = all_reduce;
  comm->base.gkyl_array_sync = array_sync;
  comm->base.gkyl_array_per_sync = array_per_sync;
  comm->base.barrier = barrier;
  comm->base.gkyl_array_write = array_write;
  comm->base.extend_comm = extend_comm;

  comm->base.ref_count = gkyl_ref_count_init(comm_free);

  return &comm->base;
}
