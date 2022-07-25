#pragma once

#if defined(__GNUC__) || defined(__GNUG__)
# include <xmmintrin.h>
#endif

/** Disable denormalized floats from occuring */
static void
disable_denorm_float(void)
{
#if defined(__GNUC__) || defined(__GNUG__)
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

#if defined(__clang__)
# if defined(__APPLE__)
  fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);  
# endif
#endif  
}
