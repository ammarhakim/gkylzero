"""LTX shape with a positive q profile preserving the old LCFS q and shear."""
import math

geometry = dict(R_axis=0.406051632, B_axis=0.240606108, R_lcfs=0.61183,
                kappa=1.3, delta=0.4)
R_omp_bounds = (0.58, 0.7)


def q(R_omp):
  # Safety factor at outboard midplane major radius R_omp [m].
  # Positive exponential preserving LCFS q=4.3153848 and shear r*q'/q=2.6899871.
  radius = R_omp-geometry["R_axis"]  # Zero shift in this example.
  a = geometry["R_lcfs"]-geometry["R_axis"]
  return 4.3153848*math.exp(0.5*2.6899871*((radius/a)**2-1))
