#ifndef GKYL_REF_COUNT_H
#define GKYL_REF_COUNT_H

// Implementation essentially as presented by Chris Wellons:
// https://nullprogram.com/blog/2015/02/17/

#include <gkyl_util.h>

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
    void (*free)(const struct gkyl_ref_count* );
    int count;
};

/**
 * Increment use count.
 *
 * @param ref Object to increment.
 */
static inline void
gkyl_ref_count_inc(const struct gkyl_ref_count *ref) {
  ((struct gkyl_ref_count *)ref)->count++;
}

/**
 * Decrement use count. If use count is zero the free (destructor)
 * function is called.
 *
 * @param ref Object to decrement.
 */
static inline void
gkyl_ref_count_dec(const struct gkyl_ref_count *ref) {
  if (--((struct gkyl_ref_count *)ref)->count == 0)
    ref->free(ref);
}

#endif // GKYL_REF_COUNT_H
