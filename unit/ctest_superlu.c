#include <acutest.h>
#include <gkyl_util.h>
#include <slu_ddefs.h>

void test_slu_example()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * Users' Guide to illustrate how to call a SuperLU routine, and the
 * matrix data structures used by SuperLU.
 *
 */
  SuperMatrix A, L, U, B;
  double   *a, *rhs;
  double   s, u, p, e, r, l;
  int      *asub, *xa;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      nrhs, info, i, m, n, nnz, permc_spec;
  superlu_options_t options;
  SuperLUStat_t stat;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nnz = 12;
  if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
  if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
  if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
  a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
  asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
  asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
  asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
  xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    
  /* Create right-hand side matrix B. */
  nrhs = 1;
  if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
  /* B : transpose([1,1,1,1,1]);*/
  for (i = 0; i < m; ++i) rhs[i] = 1.0;
  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

  if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = NATURAL;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve linear system */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  /* Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336] */
  TEST_CHECK( gkyl_compare(-1.0/32.0, rhs[0], 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/168.0, rhs[1], 1e-14) );
  TEST_CHECK( gkyl_compare( 3.0/224.0, rhs[2], 1e-14) );
  TEST_CHECK( gkyl_compare( 1.0/16.0, rhs[3], 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/336.0, rhs[4], 1e-14) );

  /* De-allocate storage */
  SUPERLU_FREE (rhs);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
}

TEST_LIST = {
  { "slu_example", test_slu_example },
  { NULL, NULL }
};
