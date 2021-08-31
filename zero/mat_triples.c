#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_range.h>

#include <stc/cmap.h>

struct mat_triple {
  size_t row, col;
  double val;
};

// map of index -> mat_triple
using_cmap(triple, int, struct mat_triple);

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
    (struct mat_triple) { .row = i, .col = j, .val = val }
  );
  return val;
}

double
gkyl_mat_triples_accum(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  assert(i<gkyl_range_shape(&tri->range, 0) && j<gkyl_range_shape(&tri->range, 1));
  
  long loc = gkyl_ridx(tri->range, i, j);
  struct cmap_triple_value_t *mt = cmap_triple_get(&tri->triples, loc);
  double tot_val = val;
  if (mt) {
    // element exists, add to its current value
    tot_val = (mt->second.val += val);
  }
  else {
    cmap_triple_put(&tri->triples, loc,
      (struct mat_triple) { .row = i, .col = j, .val = val }
    );
  }
  
  return tot_val;
}

double
gkyl_mat_triples_get(const gkyl_mat_triples *tri, size_t i, size_t j)
{
  long loc = gkyl_ridx(tri->range, i, j);
  struct cmap_triple_value_t *mt = cmap_triple_get(&tri->triples, loc);
  return mt ? mt->second.val : 0.0;
}

size_t
gkyl_mat_triples_size(const gkyl_mat_triples *tri)
{
  return cmap_triple_size(tri->triples);
}

void
gkyl_mat_triples_release(gkyl_mat_triples *tri)
{
  cmap_triple_del(&tri->triples);
  gkyl_free(tri);
}
