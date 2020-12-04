#pragma once

#include <stdio.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

/**
 * Write array data to file. Data is written as a binary file.
 *
 * @param arr Array object to write
 * @param fp File handle to write to.
 */
void gkyl_array_write(const struct gkyl_array *arr, FILE *fp);

/**
 * Write part of the array data to file. Data is written as a binary
 * file. The region of the array to write is specified in the range
 * object. This method will fail if range volume is greater than array
 * size.
 *
 * @param range Range describing portion of the array to output.
 * @param arr Array object to write
 * @param fp File handle to write to.
 */
void gkyl_sub_array_write(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp);

