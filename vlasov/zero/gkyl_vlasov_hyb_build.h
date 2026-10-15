#pragma once

#include <stddef.h>

// Build gates for the optional heavy kernel sets (tensor p=1 hybrid basis in
// 3x3v, tensor p=2 basis in 2x3v; see the top-level Makefile and ./configure
// --help). Kernel tables wrap the affected entries in these macros so the rows
// are NULL when the set is not built, and the wrapper functions that call those
// kernels are compiled out. The equation/moment constructors assert on NULL
// rows, and the apps check the gates up front with an actionable message. The
// 2x3v hybrid (tensor p=1) kernels are always built.
#ifdef GKYL_BUILD_VLASOV_HYB_3X3V
#define GKYL_HYB_3X3V(kern) kern
#else
#define GKYL_HYB_3X3V(kern) NULL
#endif

// Defined by the build only when both 3x3v knobs are set.
#ifdef GKYL_BUILD_VLASOV_HYB_3X3V_PHASE
#define GKYL_HYB_3X3V_PHASE(kern) kern
#else
#define GKYL_HYB_3X3V_PHASE(kern) NULL
#endif

// Tensor p=2 kernels in 2x3v (velocity-space Hamiltonian Vlasov, moments,
// velocity-map division, LBO/BGK and SR/GR conf-space helpers). The
// phase-space Hamiltonian and triad pieces are p=1 only at 2x3v.
#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
#define GKYL_TEN_2X3V_P2(kern) kern
#else
#define GKYL_TEN_2X3V_P2(kern) NULL
#endif
