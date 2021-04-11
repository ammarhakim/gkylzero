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
  int ndim = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2 };
  gkyl_rect_apply_bc *lbc = gkyl_rect_apply_bc_new(&grid, 0, GKYL_LOWER_EDGE, nghost, bcfunc, NULL);
  gkyl_rect_apply_bc *rbc = gkyl_rect_apply_bc_new(&grid, 0, GKYL_UPPER_EDGE, nghost, bcfunc, NULL);

  gkyl_rect_apply_bc_release(lbc);
  gkyl_rect_apply_bc_release(rbc);
}

void
test_2()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {16,8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2, 1 };
  gkyl_rect_apply_bc *lbc = gkyl_rect_apply_bc_new(&grid, 0, GKYL_LOWER_EDGE, nghost, bcfunc, NULL);
  gkyl_rect_apply_bc *rbc = gkyl_rect_apply_bc_new(&grid, 0, GKYL_UPPER_EDGE, nghost, bcfunc, NULL);

  gkyl_rect_apply_bc *bbc = gkyl_rect_apply_bc_new(&grid, 1, GKYL_LOWER_EDGE, nghost, bcfunc, NULL);
  gkyl_rect_apply_bc *tbc = gkyl_rect_apply_bc_new(&grid, 1, GKYL_UPPER_EDGE, nghost, bcfunc, NULL);

  gkyl_rect_apply_bc_release(lbc);
  gkyl_rect_apply_bc_release(rbc);
  gkyl_rect_apply_bc_release(bbc);
  gkyl_rect_apply_bc_release(tbc);
}

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },
  { NULL, NULL },
};
