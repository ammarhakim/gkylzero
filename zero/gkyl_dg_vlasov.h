#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>

/**
 * Create a new Vlasov equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);

/**
 * Set the EM field needed in updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param qmem Pointer to EM field scaled by q/,
 */
void gkyl_vlasov_em(const struct gkyl_dg_eqn *eqn, struct gkyl_array *qmem);
