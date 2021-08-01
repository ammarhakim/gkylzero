#include <gkyl_alloc.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

#include <string.h>

struct gkyl_mat*
gkyl_mat_new(size_t nr, size_t nc, double val)
{
  struct gkyl_mat *m = gkyl_malloc(sizeof(struct gkyl_mat) + sizeof(double[nr*nc]));
  m->nr = nr; m->nc = nc;
  for (size_t i=0; i<nr*nc; ++i) m->data[i] = val;
  return m;
}

struct gkyl_mat*
gkyl_mat_clone(const struct gkyl_mat *in)
{
  size_t tot = sizeof(struct gkyl_mat) + sizeof(double[in->nr*in->nc]);
  struct gkyl_mat *m = gkyl_malloc(tot);
  memcpy(m, in, tot);
  return m;
}

struct gkyl_mat*
gkyl_mat_clear(struct gkyl_mat *mat, double val)
{
  for (size_t i=0; i<mat->nr*mat->nc; ++i) mat->data[i] = val;
  return mat;
}

struct gkyl_mat*
gkyl_mat_diag(struct gkyl_mat *mat, double val)
{
  gkyl_mat_clear(mat, 0.0);
  for (size_t i=0; i<GKYL_MIN(mat->nr, mat->nc); ++i)
    gkyl_mat_set(mat, i, i, val);
  return mat;
}

void
gkyl_mat_release(struct gkyl_mat *mat)
{
  if (mat) gkyl_free(mat);
  mat = 0;
}
