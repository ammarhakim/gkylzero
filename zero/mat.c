#include "gkyl_ref_count.h"
#include <gkyl_alloc.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
#include <stdbool.h>

// BLAS and LAPACKE includes
#ifdef GKYL_USING_FRAMEWORK_ACCELERATE
#include <Accelerate/Accelerate.h>
#else
// On non-Darwin platforms use OpenBLAS
# include <cblas.h>
# include <lapacke.h>
#endif

#include <assert.h>
#include <string.h>

/** Map Gkyl flags to CBLAS flags */
static int cblas_trans_flags[] = {
  [GKYL_NO_TRANS] = CblasNoTrans,
  [GKYL_TRANS] = CblasTrans,
  [GKYL_CONJ_TRANS] = CblasConjTrans
};

// flags and corresponding bit-masks
enum nmat_flags { M_IS_CU_NMAT };
static const uint32_t masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// NV-GPU flags
#define SET_CU_NMAT(flags) (flags) |= masks[M_IS_CU_NMAT]
#define CLEAR_CU_NMAT(flags) (flags) &= ~masks[M_IS_CU_NMAT]
#define IS_CU_NMAT(flags) (((flags) & masks[M_IS_CU_NMAT]) != 0)

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
  struct gkyl_mat *m = gkyl_malloc(sizeof(struct gkyl_mat));
  m->data = gkyl_malloc(sizeof(double[nr*nc]));
  m->nr = nr; m->nc = nc;
  for (size_t i=0; i<nr*nc; ++i) m->data[i] = val;
  return m;
}

struct gkyl_mat*
gkyl_mat_clone(const struct gkyl_mat *in)
{
  struct gkyl_mat *m = gkyl_malloc(sizeof(struct gkyl_mat));
  m->data = gkyl_malloc(sizeof(double[in->nr*in->nc]));
  m->nc = in->nc; m->nr = in->nr;
  size_t tot = sizeof(double[in->nr*in->nc]);
  memcpy(m->data, in->data, tot);
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
gkyl_mat_linsolve_lu(struct gkyl_mat *A, struct gkyl_mat *x, void* ipiv)
{
  assert( A->nr == A->nc );

#ifdef GKYL_USING_FRAMEWORK_ACCELERATE
  // On Darwin need to use old clapack interface. Of course Apple has
  // to do everything "different"
  __CLPK_integer info;
  __CLPK_integer n = A->nr;
  __CLPK_integer nrhs = x->nc;
  __CLPK_integer lda = A->nr;
  __CLPK_integer ldb = A->nr;
  dgesv_(&n, &nrhs, A->data, &lda, ipiv, x->data, &ldb, &info);
#else
  // on non-Darwin platforms modern LAPACKE interface is available
  int info = LAPACKE_dgesv(LAPACK_COL_MAJOR,
    A->nr, x->nc, A->data, A->nr, ipiv, x->data, A->nr);
#endif
  
  return info == 0 ? true : false;
}

void
gkyl_mat_release(struct gkyl_mat *mat)
{
  if (mat) {
    gkyl_free(mat->data);
    gkyl_free(mat);
  }
}

static void
nmat_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_nmat *mat = container_of(ref, struct gkyl_nmat, ref_count);
  gkyl_free(mat->data);
  gkyl_free(mat->mptr);
  gkyl_free(mat);
}

struct gkyl_nmat*
gkyl_nmat_new(size_t num, size_t nr, size_t nc)
{
  struct gkyl_nmat *mat = gkyl_malloc(sizeof(struct gkyl_nmat));
  mat->num = num; mat->nr = nr; mat->nc = nc;

  mat->flags = 0;
  mat->data = gkyl_malloc(sizeof(double[num*nr*nc]));  
  mat->mptr = gkyl_malloc(num*sizeof(double*));
  
  for (size_t i=0; i<num; ++i)
    mat->mptr[i] = mat->data+nr*nc*i;

  mat->on_dev = mat; // on CPU this is a self-reference

  mat->ref_count = gkyl_ref_count_init(nmat_free);

  return mat;
}

bool
gkyl_nmat_is_cu_dev(const struct gkyl_nmat *mat)
{
  return IS_CU_NMAT(mat->flags);
}

struct gkyl_nmat*
gkyl_nmat_acquire(const struct gkyl_nmat *mat)
{
  gkyl_ref_count_inc(&mat->ref_count);
  return (struct gkyl_nmat*) mat;
}

bool
gkyl_nmat_linsolve_lu(struct gkyl_nmat *A, struct gkyl_nmat *x)
{
  size_t num = A->num;
  assert( num <= x->num );

  bool status = true;

  long *ipiv = gkyl_malloc(sizeof(long[A->nr]));
  for (size_t i=0; i<num; ++i) {
    struct gkyl_mat Ai = gkyl_nmat_get(A,i);
    struct gkyl_mat xi = gkyl_nmat_get(x,i);
    status = gkyl_mat_linsolve_lu( &Ai, &xi, ipiv );
    if (!status) break;
  }
  gkyl_free(ipiv);

  return status;
}

void
gkyl_nmat_release(struct gkyl_nmat *mat)
{
  if (mat)
    gkyl_ref_count_dec(&mat->ref_count);
}

// CUDA specific code

#ifdef GKYL_HAVE_CUDA

struct gkyl_nmat*
gkyl_nmat_cu_dev_new(size_t num, size_t nr, size_t nc)
{
  assert(false);
  return 0;
}

#else

struct gkyl_nmat*
gkyl_nmat_cu_dev_new(size_t num, size_t nr, size_t nc)
{
  assert(false);
  return 0;
}

#endif // CUDA specific code
