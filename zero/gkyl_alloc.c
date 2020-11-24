#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>

// Compute first 'align' boundary after 'num'
#define align_up(num, align) \
    (((num) + ((align) - 1)) & ~((align) - 1))

typedef uint16_t offset_t; // for storing offset 
#define PTR_OFFSET_SZ sizeof(offset_t)

void*
gkyl_aligned_alloc(size_t align, size_t size)
{
  void *ptr = NULL;
  assert((align & (align - 1)) == 0); // check power of 2

  if (align && size) {
    uint32_t hdr_size = PTR_OFFSET_SZ + (align - 1);
    void *p = gkyl_malloc(size + hdr_size);

    if (p) {
      ptr = (void *) align_up(((uintptr_t)p + PTR_OFFSET_SZ), align);
      // calculate and store offset (just behind our aligned pointer)
      *((offset_t *)ptr - 1) = (offset_t)((uintptr_t)ptr - (uintptr_t)p);
    }
  }
  return ptr;
}

void
gkyl_aligned_free(void* ptr)
{
  assert(ptr);
  offset_t offset = *((offset_t *)ptr - 1);
  void *p = (void *)((uint8_t *)ptr - offset);
  gkyl_free(p);
}
