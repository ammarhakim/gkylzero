#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <time.h>

// random number generator
#include <pcg_basic.h>

// Maximum configuration-space dimensions supported
#ifndef GKYL_MAX_CDIM
# define GKYL_MAX_CDIM 3
#endif

// Maximum dimensions supported
#ifndef GKYL_MAX_DIM
# define GKYL_MAX_DIM 7
#endif

// Maximum number of supported species
#ifndef GKYL_MAX_SPECIES
# define GKYL_MAX_SPECIES 8
#endif


// Default alignment boundary
#ifndef GKYL_DEF_ALIGN
# define GKYL_DEF_ALIGN 64
#endif

// CUDA specific defines etc
#ifdef __NVCC__

#define GKYL_HAVE_CUDA
#define GKYL_CU_DH __device__ __host__

#else

#undef GKYL_HAVE_CUDA
#define GKYL_CU_DH

#endif // CUDA specific defines etc

// This funny looking macro allows getting a pointer to the 'type'
// struct that contains an object 'member' given the 'ptr' to the
// 'member' inside 'type'. (Did I just write this gobbledygook?!)
//
// See https://en.wikipedia.org/wiki/Offsetof
#define container_of(ptr, type, member)                                 \
    ((type *)((char *)(1 ? (ptr) : &((type *)0)->member) - offsetof(type, member)))

// Select type-specific compare function
#define gkyl_compare(a, b, eps)                 \
    _Generic((a),                               \
      float: gkyl_compare_float,                \
      double: gkyl_compare_double)              \
    (a, b, eps)

/**
 * Time-trigger. Typical initialization is:
 * 
 * struct gkyl_tm_trigger tmt = { .dt = tend/nframe };
 */
struct gkyl_tm_trigger {
  int curr; // current counter
  double dt, tcurr; // Time-interval, current time
};

/**
 * Check if the tcurr should trigger and bump internal counters if it
 * does. This only works if sequential calls to this method have the
 * tcurr monotonically increasing.
 *
 * @param tmt Time trigger object
 * @param tcurr Current time.
 * @return 1 if triggered, 0 otherwise
 */
int gkyl_tm_trigger_check_and_bump(struct gkyl_tm_trigger *tmt, double tcurr);

/**
 * Print error message to stderr and exit.
 *
 * @param msg Error message.
 */
void gkyl_exit(const char* msg);

/**
 * Compares two float numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_float(float a, float b, float eps);

/**
 * Compares two double numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_double(double a, double b, double eps);

/**
 * Copy (small) int arrays.
 *
 * @param n Number of elements to copy
 * @param inp Input array
 * @param out Output array
 */
static inline void
gkyl_copy_int_arr(int n, const int *restrict inp, int *restrict out)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

/**
 * Gets wall-clock time in secs/nanoseconds.
 * 
 * @return Time object.
 */
struct timespec gkyl_wall_clock(void);

/**
 * Difference between two timespec objects.
 *
 * @param tstart Start time
 * @param tend End time 
 * @return Time object representing difference
 */
struct timespec gkyl_time_diff(struct timespec tstart, struct timespec tend);

/**
 * Difference between timespec object and "now", returned in seconds.
 * 
 * @param tm Timespec
 * @return Time in seconds
 */
double gkyl_time_diff_now_sec(struct timespec tm);

/**
 * Compute in secs time stored in timespec object.
 *
 * @param tm Timespec object
 * @return Time in seconds
 */
double gkyl_time_sec(struct timespec tm);

/**
 * Initialize 32-bit RNG 
 *
 * @param nd_seed true if to use non-deterministic seed, or false for
 *   a determistic seed.
 */
pcg32_random_t gkyl_pcg32_init(bool nd_seed);

/**
 * Returns a unsigned 32-bit random number
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed 32-bit integer
 */
uint32_t gkyl_pcg32_rand_uint32(pcg32_random_t* rng);

/**
 * Returns a random number in [0,1), rounded to the nearest
 * 1/2^32. This may not be enough resolution for some applications.
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed double in [0,1)
 */
double gkyl_pcg32_rand_double(pcg32_random_t* rng);

/** 64-bit RNG: concatination of two 32-bit RNGs */
typedef struct { pcg32_random_t gen[2]; } pcg64_random_t;

/**
 * Initialize 64-bit RNG 
 *
 * @param nd_seed true if to use non-deterministic seed, or false for
 *   a determistic seed.
 */
pcg64_random_t gkyl_pcg64_init(bool nd_seed);

/**
 * Returns a unsigned 64-bit random number
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed 64-bit integer
 */
uint64_t gkyl_pcg64_rand_uint32(pcg64_random_t* rng);

/**
 * Returns a random number in [0,1), rounded to the nearest
 * 1/2^64.
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed double in [0,1)
 */
double gkyl_pcg64_rand_double(pcg64_random_t* rng);
