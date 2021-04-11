#include <acutest.h>

#include <gkyl_rect_apply_bc.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

void
bcfunc(double t, int dir, const double *skin, double *restrict ghost, void *ctx)
{
  for (int d=0; d<2; ++d) ghost[d] = skin[d];
}

void
test_1()
{
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  gkyl_rect_apply_bc *lbc = gkyl_rect_apply_bc_new(&grid, 0, 0, bcfunc, NULL);
  gkyl_rect_apply_bc *rbc = gkyl_rect_apply_bc_new(&grid, 0, 1, bcfunc, NULL);

  gkyl_rect_apply_bc_release(lbc);
  gkyl_rect_apply_bc_release(rbc);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
