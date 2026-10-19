#include <gkyl_flash_decomp.h>

// Internal struct
struct flash_decomp {
  struct gkyl_flash_decomp decomp;
  int *branks; // Ranks per block
};

const struct gkyl_flash_decomp*
gkyl_flash_decomp_new(int total_ranks, int nblocks, int ndim, const int *branks) 
{
  struct flash_decomp *fdecomp = gkyl_malloc(sizeof(struct flash_decomp));
  fdecomp->decomp.total_ranks = total_ranks;
  fdecomp->decomp.nblocks = nblocks;
  fdecomp->decomp.ndim = ndim;

  fdecomp->branks = gkyl_malloc(nblocks*sizeof(int));
  for (int i=0; i<nblocks; ++i) {
    fdecomp->branks[i] = branks[i];
  }

  int tot_branks = 0;
  for (int i=0; i<nblocks; ++i) {
    tot_branks += branks[i];
  }

  return &fdecomp->decomp;
}

int 
gkyl_flash_decomp_nranks(const struct gkyl_flash_decomp *decomp, int bn)
{
  struct flash_decomp *fdecomp = container_of(decomp, struct flash_decomp, decomp);
  return fdecomp->branks[bn];
}

int 
gkyl_flash_decomp_ndim(const struct gkyl_flash_decomp *decomp)
{
  return decomp->ndim;
}

void
gkyl_flash_decomp_getranks(const struct gkyl_flash_decomp *decomp, int bn, int ranks[])
{ 
  struct flash_decomp *fdecomp = container_of(decomp, struct flash_decomp, decomp);

  int ndim = gkyl_flash_decomp_ndim(decomp);

  // HARDCODE: FLASH uniform grid uses a single rank per block.
  ranks[0] = bn;
}

void
gkyl_flash_decomp_release(const struct gkyl_flash_decomp *decomp)
{
  struct flash_decomp *fdecomp = container_of(decomp, struct flash_decomp, decomp);
  gkyl_free(fdecomp->branks);
  gkyl_free(fdecomp);
}
