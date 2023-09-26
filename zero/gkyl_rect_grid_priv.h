#pragma once

#include <gkyl_rect_grid.h>

bool isInCell(const struct gkyl_rect_grid *grid, const double *pIn, int *iIn, int *dimTrans[1], const int **knownIdx);

void InDir(const struct gkyl_rect_grid *grid, int *iIn, int *dimTrans[1], const int **knownIdx, double lowerInDir[], double upperInDir[]);
