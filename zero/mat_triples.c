#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>

struct gkyl_mat_triples {
  
};

gkyl_mat_triples*
gkyl_mat_triples_new()
{
  struct gkyl_mat_triples *tri = gkyl_malloc(sizeof(struct gkyl_mat_triples));

  return tri;
}

void
gkyl_mat_triples_insert(int i, int j, double val)
{
  // TODO
}

double
gkyl_mat_triples_accum(int i, int j, double val)
{
  // TODO

  return val;
}

void
gkyl_mat_triples_release(gkyl_mat_triples *tri)
{
  gkyl_free(tri);
}
