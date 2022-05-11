#include <gkyl_alloc.h>
#include <gkyl_dynvec.h>
#include <gkyl_ref_count.h>

#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

/* Size by which vector grows each time it is reallocated */
static const size_t DYNVEC_ALLOC_SZ = 1024;

// code for array datatype for use in IO
static const uint64_t array_data_type[] = {
  [GKYL_INT] = 0,
  [GKYL_FLOAT] = 1,
  [GKYL_DOUBLE] = 2,
  [GKYL_USER] = 32,
};

// size in bytes for various data-types
static const size_t array_elem_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};

// type-id for field data
static const uint64_t dynvec_file_type = 2;

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
  
  dv->ref_count = gkyl_ref_count_init(dynvec_free);
  
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

size_t
gkyl_dynvec_capacity(const gkyl_dynvec vec)
{
  return vec->csize;
}

void
gkyl_dynvec_clear(gkyl_dynvec dv)
{
  dv->cloc = 0;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
}

void
gkyl_dynvec_clear_all_but(gkyl_dynvec dv, size_t num)
{
  if (num>dv->cloc) return;
  
  size_t cloc = dv->cloc;
  dv->csize = DYNVEC_ALLOC_SZ;

  void *data = gkyl_malloc(num*dv->esznc);
  double *tm_mesh = gkyl_malloc(sizeof(double[num]));

  size_t low = num>cloc ? 0 : cloc-num; // lower index to copy from
  size_t ncpy = num>cloc ? cloc : num; // number of elemetns to copy
  dv->cloc = ncpy;

  memcpy(tm_mesh, dv->tm_mesh+low, ncpy*sizeof(double));
  memcpy(data, (char*)dv->data+low*dv->esznc, ncpy*dv->esznc);
  
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));

  memcpy(dv->data, data, ncpy*dv->esznc);
  memcpy(dv->tm_mesh, tm_mesh, ncpy*sizeof(double));

  gkyl_free(data);
  gkyl_free(tm_mesh);
}

gkyl_dynvec
gkyl_dynvec_acquire(const gkyl_dynvec vec)
{
  gkyl_ref_count_inc(&vec->ref_count);
  return (struct gkyl_dynvec_tag*) vec;
}

/*
  IF THIS FORMAT IF MODIFIED, PLEASE COPY AND THEN CHANGE THE
  DESCRIPTION SO WE HAVE THE OLDER VERSIONS DOCUMENTED HERE. UPDATE
  VERSION BY 1 EACH TIME YOU CHANGE THE FORMAT.

  The format of the gkyl binary output is as follows.

  ## Version 1: May 9th 2022. Created by A.H

  Data      Type and meaning
  --------------------------
  gkyl0     5 bytes
  version   uint64_t 
  file_type uint64_t (1: field data, 2: diagnostic data)
  meta_size uint64_t Number of bytes of meta-data
  DATA      meta_size bytes of data. This is in msgpack format

  For file_type = 2 (dynvec) the above header is followed by

  real_type uint64_t. Indicates real type of data 
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  TIME_DATA float64[size] bytes of data
  DATA      size*esznc bytes of data
  
 */

int
gkyl_dynvec_write(const gkyl_dynvec vec, const char *fname)
{
  const char g0[5] = "gkyl0";

  FILE *fp = 0;
  with_file (fp, fname, "w") {  
    // Version 1 header
    fwrite(g0, sizeof(char[5]), 1, fp);
    uint64_t version = 1;
    fwrite(&version, sizeof(uint64_t), 1, fp);
    fwrite(&dynvec_file_type, sizeof(uint64_t), 1, fp);
    uint64_t meta_size = 0; // THIS WILL CHANGE ONCE METADATA IS EMBEDDED
    fwrite(&meta_size, sizeof(uint64_t), 1, fp);
    
    uint64_t real_type = array_data_type[vec->type];
    fwrite(&real_type, sizeof(uint64_t), 1, fp);

    uint64_t esznc = vec->esznc, size = gkyl_dynvec_size(vec);
    fwrite(&esznc, sizeof(uint64_t), 1, fp);
    fwrite(&size, sizeof(uint64_t), 1, fp); 

    fwrite(vec->tm_mesh, sizeof(double)*size, 1, fp);
    fwrite(vec->data, esznc*size, 1, fp);
  }

  return errno;
}

bool
gkyl_dynvec_read(gkyl_dynvec vec, const char *fname)
{

  return false;
}

void
gkyl_dynvec_release(gkyl_dynvec vec)
{
  if (vec)
    gkyl_ref_count_dec(&vec->ref_count);
}
