#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>

#include <stc/cvec.h>

// define vector of doubles
using_cvec(d, double);

struct gkyl_mat_triples {
  cvec_d nzval; // List of non-zero values in matrix
  
};

gkyl_mat_triples*
gkyl_mat_triples_new(size_t nr, size_t nc)
{
  struct gkyl_mat_triples *tri = gkyl_malloc(sizeof(struct gkyl_mat_triples));

  tri->nzval = cvec_d_init();

  return tri;
}

void
gkyl_mat_triples_insert(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  // TODO
}

double
gkyl_mat_triples_accum(gkyl_mat_triples *tri, size_t i, size_t j, double val)
{
  // TODO

  return val;
}

void
gkyl_mat_triples_release(gkyl_mat_triples *tri)
{
  cvec_d_del(&tri->nzval);
  gkyl_free(tri);
}
