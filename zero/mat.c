#include <gkyl_alloc.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

// BLAS and LAPACKE includes
#include <cblas.h>
#include <lapacke.h>

#include <assert.h>
#include <string.h>

/** Map Gkyl flags to CBLAS flags */
static int cblas_trans_flags[] = {
  [GKYL_NO_TRANS] = CblasNoTrans,
  [GKYL_TRANS] = CblasTrans,
  [GKYL_CONJ_TRANS] = CblasConjTrans
};

/** Helper functions to determine sizes needed in BLAS/LAPACKE routines */
struct mat_sizes { size_t nr, nc; };

static inline struct mat_sizes
get_mat_sizes(enum gkyl_mat_trans trans, const struct gkyl_mat *A)
{
  if (trans == GKYL_NO_TRANS)
    return (struct mat_sizes) { .nr = A->nr, .nc = A->nc };
  return (struct mat_sizes) { .nr = A->nc, .nc = A->nr };
}

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
gkyl_mat_show(const char *name, FILE *fp, const struct gkyl_mat *mat)
{
  fprintf(fp, "%s : matrix( ", name);

  for (int i=0; i<mat->nr-1; ++i) {
    fprintf(fp, "[");
    for (int j=0; j<mat->nc-1; ++j) {
      fprintf(fp, "%lg, ", gkyl_mat_get(mat,i,j));
    }
    fprintf(fp, "%lg ", gkyl_mat_get(mat,i,mat->nc-1));
    fprintf(fp, "], ");
  }

  fprintf(fp, "[");
  for (int j=0; j<mat->nc-1; ++j) {
    fprintf(fp, "%lg, ", gkyl_mat_get(mat,mat->nr-1,j));
  }
  fprintf(fp, "%lg ", gkyl_mat_get(mat,mat->nr-1,mat->nc-1));
  fprintf(fp, "] ");
  
  fprintf(fp, " )\n");
}

struct gkyl_mat*
gkyl_mat_mm(double alpha, double beta,
  enum gkyl_mat_trans transa, const struct gkyl_mat *A,
  enum gkyl_mat_trans transb, const struct gkyl_mat *B, struct gkyl_mat *C)
{
  // determine matrix sizes
  struct mat_sizes sza = get_mat_sizes(transa, A);
  struct mat_sizes szb = get_mat_sizes(transb, B);
  struct mat_sizes szc = get_mat_sizes(GKYL_NO_TRANS, C);

  // intermediate size
  size_t k = sza.nc; // same as szb.nr
  size_t lda = transa == GKYL_NO_TRANS ? C->nr : k;
  size_t ldb = transb == GKYL_NO_TRANS ? k : C->nc;
  size_t ldc = C->nr;
  
  assert( (sza.nr == szc.nr) && (sza.nc == k) && (szb.nr == k) && (szb.nc == szc.nc) );

  // call BLAS routine to perform matrix-matrix multiply
  cblas_dgemm(CblasColMajor,
    cblas_trans_flags[transa],
    cblas_trans_flags[transb],
    C->nr, C->nc, k,
    alpha,
    A->data, lda,
    B->data, ldb,
    beta, C->data, ldc);

  return C;
}

bool
gkyl_mat_linsolve_lu(struct gkyl_mat *A, struct gkyl_mat *x, struct gkyl_mem_chunk *ipiv)
{
  assert( A->nr == A->nc );
  assert( ipiv->count >= sizeof(lapack_int[A->nr]) );
  
  int info = LAPACKE_dgesv(LAPACK_COL_MAJOR,
    A->nr, x->nc, A->data, A->nr, ipiv->data, x->data, A->nr);
  return info == 0 ? true : false;
}

void
gkyl_mat_release(struct gkyl_mat *mat)
{
  if (mat) gkyl_free(mat);
  mat = 0;
}
