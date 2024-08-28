#pragma once

/**
* Rotate the pressure tensor for the 10-moment equations with respect to the magnetic field.
*
* @param q_over_m Ratio of species charge to species mass.
* @param dt Current stable time-step.
* @param em Array of electromagnetic variables.
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param p_tensor_old Old components of the pressure tensor (before rotation).
* @param p_tensor_rhs Source terms appearing on the right-hand-side of the pressure tensor equations.
* @param p_tensor_new New components of the pressure tensor (after rotation).
*/
void pressure_tensor_rotate(double q_over_m, double dt, const double* em, const double* ext_em, double p_tensor_old[6], double p_tensor_rhs[6],
  double p_tensor_new[6]);