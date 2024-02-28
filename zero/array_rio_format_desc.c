#include <gkyl_array_rio_format_desc.h>

#include <stdint.h>
#include <stdio.h>

size_t
gkyl_base_hdr_size(size_t meta_sz)
{
  size_t sz = 0;
  // magic string
  sz += 5; // "gkyl0"
  // version
  sz += sizeof(uint64_t);
  // file_type
  sz += sizeof(uint64_t);
  // metadata
  sz += sizeof(uint64_t) + meta_sz;

  return sz;
}

size_t
gkyl_file_type_1_partial_hrd_size(int ndim)
{
  size_t sz = 0;
  // real_type
  sz += sizeof(uint64_t);
  // ndim
  sz += sizeof(uint64_t);
  // cells
  sz += sizeof(uint64_t[ndim]);
  // lower, upper
  sz += sizeof(double[2*ndim]);
  return sz;
}

size_t
gkyl_file_type_1_hrd_size(int ndim)
{
  size_t sz = gkyl_file_type_1_partial_hrd_size(ndim);
  // esznc
  sz += sizeof(uint64_t);
  // size (total number of cells)
  sz += sizeof(uint64_t);
  return sz;
}

size_t
gkyl_file_type_2_hrd_size(void)
{
  size_t sz = 0;
  // real_type
  sz += sizeof(uint64_t);
  // esznc
  sz += sizeof(uint64_t);
  // size (total number of cells)
  sz += sizeof(uint64_t);
  return sz;
}

size_t
gkyl_file_type_3_hrd_size(int ndim)
{
  size_t sz = gkyl_file_type_1_hrd_size(ndim);
  sz += sizeof(uint64_t); // nrange
  return sz;
}

size_t
gkyl_file_type_3_range_hrd_size(int ndim)
{
  size_t sz = 0;
  // loidx and upidx
  sz += sizeof(uint64_t[2*ndim]);
  sz += sizeof(uint64_t);
  return sz;
}
