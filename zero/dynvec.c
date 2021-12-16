#include <gkyl_alloc.h>
#include <gkyl_dynvec.h>
#include <gkyl_ref_count.h>

#include <stdbool.h>
#include <stddef.h>
#include <string.h>

/* Size by which vector grows each time it is reallocated */
static const size_t DYNVEC_ALLOC_SZ = 1024;

// size in bytes for various data-types
static const size_t array_elem_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};

struct gkyl_dynvec_tag {
  enum gkyl_elem_type type; // type of data stored in vector
  size_t elemsz, ncomp; // size of elements, number of 'components'

  size_t cloc; // current location to insert data into
  size_t csize; // current number of elements allocated

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  double *tm_mesh; // time stamps
  
  struct gkyl_ref_count ref_count;  
};

static void
dynvec_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dynvec_tag *dv = container_of(ref, struct gkyl_dynvec_tag, ref_count);
  gkyl_free(dv->data);
  gkyl_free(dv->tm_mesh);
  gkyl_free(dv);
}

gkyl_dynvec
gkyl_dynvec_new(enum gkyl_elem_type type, size_t ncomp)
{
  struct gkyl_dynvec_tag *dv = gkyl_malloc(sizeof(struct gkyl_dynvec_tag));
  dv->type = type;
  dv->elemsz = array_elem_size[type];
  dv->ncomp = ncomp;
  dv->esznc = dv->elemsz*dv->ncomp;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->cloc = 0;
  
  dv->data = gkyl_calloc(dv->csize, dv->esznc);
  dv->tm_mesh = gkyl_calloc(dv->csize, sizeof(double));
  
  dv->ref_count = (struct gkyl_ref_count) { dynvec_free, 1 }; 
  
  return dv;
}

void
gkyl_dynvec_append(gkyl_dynvec dv, double tm, const void *data)
{
  size_t loc = dv->cloc;
  if (loc >= dv->csize) {
    dv->csize += DYNVEC_ALLOC_SZ;
    dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
    dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
  }

  // set data
  dv->tm_mesh[loc] = tm;
  memcpy((char*)dv->data + dv->esznc*loc, data, dv->esznc);
  
  dv->cloc += 1;
}

bool
gkyl_dynvec_get(const gkyl_dynvec dv, size_t idx, void *data)
{
  if (idx >= dv->cloc) return false;
  memcpy(data, (char*)dv->data + dv->esznc*idx, dv->esznc);
  return true;  
}

double
gkyl_dynvec_get_tm(const gkyl_dynvec dv, size_t idx)
{
  if (idx >= dv->cloc) return 0.0;
  return dv->tm_mesh[idx];
}

bool
gkyl_dynvec_getlast(const gkyl_dynvec dv, void *data)
{
  size_t loc = dv->cloc;
  if (loc == 0) return false;
  return gkyl_dynvec_get(dv, loc-1, data);
}

double
gkyl_dynvec_getlast_tm(const gkyl_dynvec dv)
{
  size_t loc = dv->cloc;
  return loc == 0 ? 0.0 : dv->tm_mesh[loc-1];
}

size_t
gkyl_dynvec_size(const gkyl_dynvec vec)
{
  return vec->cloc;
}

void
gkyl_dynvec_clear(gkyl_dynvec dv)
{
  dv->cloc = 0;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
}

gkyl_dynvec
gkyl_dynvec_aquire(const gkyl_dynvec vec)
{
  return vec;
}

void
gkyl_dynvec_release(gkyl_dynvec vec)
{
  if (vec)
    gkyl_ref_count_dec(&vec->ref_count);
}
