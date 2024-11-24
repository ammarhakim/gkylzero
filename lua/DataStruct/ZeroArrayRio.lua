-- Gkyl ------------------------------------------------------------------------
--
-- Lua interface to gkylzero's gkyl_array structure
---
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

_M = {}

ffi.cdef [[
// Read status flags
enum gkyl_array_rio_status {
  GKYL_ARRAY_RIO_SUCCESS = 0,
  GKYL_ARRAY_RIO_BAD_VERSION,
  GKYL_ARRAY_RIO_FOPEN_FAILED,
  GKYL_ARRAY_RIO_FREAD_FAILED,
  GKYL_ARRAY_RIO_DATA_MISMATCH,
  GKYL_ARRAY_RIO_META_FAILED
};

// Structure to pass meta-data to write methods
struct gkyl_array_meta {
  size_t meta_sz; // size in bytes of meta-data
  char *meta; // meta-data encoded in mpack format
};

/**
 * Write out grid and array data to file in .gkyl format so postgkyl
 * can understand it.
 *
 * @param grid Grid object to write
 * @param range Range describing portion of the array to output.
 * @param meta Meta-data to write. Set to NULL or 0 if no metadata
 * @param arr Array object to write
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, const struct gkyl_array_meta *meta,
  const struct gkyl_array *arr, const char *fname);
]]
