#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_range.h>

// define map of index -> gkyl_mtriple
#define i_key long
#define i_val struct gkyl_mtriple
#define i_tag triple
#include <stc/cmap.h>
// complete definition of map

// define sorted (column-major) map of index[2] -> key.
struct mat_idx {
  size_t row, col;
};

static inline int
cmp_long(long a, long b)
{
  return a==b ? 0 : a<b ? -1 : 1;
}

static inline int
mat_idx_cmp(const struct mat_idx *mia, const struct mat_idx *mib)
{
  return cmp_long(mia->col, mib->col) == 0 ? cmp_long(mia->row, mib->row) :
    cmp_long(mia->col, mib->col);
}

#define i_key struct mat_idx
#define i_val long 
#define i_cmp mat_idx_cmp
#define i_tag idx2key
#include <stc/csmap.h>
// complete definition of sorted map.

struct gkyl_mat_triples {
  struct gkyl_range range; // range representing matrix shape
  cmap_triple triples; // non-zero values in matrix
};

gkyl_mat_triples*
gkyl_mat_triples_new(size_t nr, size_t nc)
{
  struct gkyl_mat_triples *tri = gkyl_malloc(sizeof(struct gkyl_mat_triples));

  gkyl_range_init_from_shape(&tri->range, 2, (const int[]) { nr, nc} );
  tri->triples = cmap_triple_init();

  return tri;
}

double
gkyl_mat_triples_insert(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  assert(i<gkyl_range_shape(&tri->range, 0) && j<gkyl_range_shape(&tri->range, 1));
  
  long loc = gkyl_ridx(tri->range, i, j);
  cmap_triple_put(&tri->triples, loc,
    (struct gkyl_mtriple) { .row = i, .col = j, .val = val }
  );
  return val;
}

double
gkyl_gkyl_mat_triples_accum(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  assert(i<gkyl_range_shape(&tri->range, 0) && j<gkyl_range_shape(&tri->range, 1));
  
  long loc = gkyl_ridx(tri->range, i, j);
  struct cmap_triple_value *mt = cmap_triple_get_mut(&tri->triples, loc);
  double tot_val = val;
  if (mt) {
    // element exists, add to its current value
    tot_val = (mt->second.val += val);
  }
  else {
    cmap_triple_put(&tri->triples, loc,
      (struct gkyl_mtriple) { .row = i, .col = j, .val = val }
    );
  }
  
  return tot_val;
}

double
gkyl_mat_triples_get(const gkyl_mat_triples *tri, size_t i, size_t j)
{
  long loc = gkyl_ridx(tri->range, i, j);
  const struct cmap_triple_value *mt = cmap_triple_get(&tri->triples, loc);
  return mt ? mt->second.val : 0.0;
}

size_t
gkyl_mat_triples_size(const gkyl_mat_triples *tri)
{
  return cmap_triple_size(tri->triples);
}

long*
gkyl_mat_triples_keys(const gkyl_mat_triples *tri)
{
  long trisize = cmap_triple_size(tri->triples);
  long *keys = gkyl_malloc(sizeof(long[trisize]));
  long i = 0;
  c_foreach(kv, cmap_triple, tri->triples) {
    keys[i] = kv.ref->first;
    i += 1;
  }
  return keys;
}

long*
gkyl_mat_triples_keys_colmo(const gkyl_mat_triples *tri)
{
  // create a sorted map, sorted by row.
  csmap_idx2key ikmap;
  ikmap = csmap_idx2key_init();
  c_foreach(kv, cmap_triple, tri->triples) {
    csmap_idx2key_put(&ikmap, (struct mat_idx) {
        .row = kv.ref->second.row, .col = kv.ref->second.col
      },
      kv.ref->first);
  }
  
  // return array of keys from sorted map.
  long ikmapsize = csmap_idx2key_size(ikmap);
  long *keys = gkyl_malloc(sizeof(long[ikmapsize]));
  long i = 0;
  c_foreach(kv, csmap_idx2key, ikmap) {
    keys[i] = kv.ref->second;
    i += 1;
  }
  csmap_idx2key_drop(&ikmap);
  return keys;
}

double
gkyl_mat_triples_val_at_key(const gkyl_mat_triples *tri, long loc)
{
  const struct cmap_triple_value *mt = cmap_triple_get(&tri->triples, loc);
  return mt ? mt->second.val : 0.0;
}

void
gkyl_mat_triples_key_to_idx(const gkyl_mat_triples *tri, long loc, int idx[2])
{
  const struct cmap_triple_value *mt = cmap_triple_get(&tri->triples, loc);
  if (mt) {
    idx[0] = mt->second.row;
    idx[1] = mt->second.col;
  }
  else {
    // if location is not in triple list, just inverse map it
    gkyl_range_inv_idx(&tri->range, loc, idx);
  }
}

void
gkyl_mat_triples_release(gkyl_mat_triples *tri)
{
  cmap_triple_drop(&tri->triples);
  gkyl_free(tri);
}
