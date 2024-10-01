#pragma once
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_fem_poisson_bctype.h>

struct deflated_fem_data {
  struct gkyl_array *deflated_field;
  struct gkyl_array *deflated_phi;
  struct gkyl_array *deflated_epsilon;
  struct gkyl_array *deflated_kSq;
  struct gkyl_array *deflated_nodal_fld;
  struct gkyl_fem_poisson *fem_poisson;
};

// Updater type
struct gkyl_deflated_fem_poisson {
  bool ishelmholtz; // if solving Helmholtz equation (kSq is not zero/NULL).
  struct gkyl_rect_grid grid;
  struct gkyl_rect_grid deflated_grid;
  struct gkyl_basis basis;
  struct gkyl_basis *basis_on_dev;
  struct gkyl_basis deflated_basis;
  struct gkyl_basis *deflated_basis_on_dev;
  struct gkyl_range local;
  struct gkyl_range deflated_local;
  struct gkyl_range global_sub_range;
  struct gkyl_range deflated_local_ext;
  struct gkyl_range nrange;
  struct gkyl_range deflated_nrange;
  struct gkyl_poisson_bc poisson_bc;
  struct deflated_fem_data *d_fem_data;
  struct gkyl_array *nodal_fld;
  int num_solves_z;
  int cdim;
  struct gkyl_nodal_ops *n2m;
  struct gkyl_nodal_ops *n2m_deflated;
  struct gkyl_deflate_zsurf *deflator_lo;
  struct gkyl_deflate_zsurf *deflator_up;
  bool use_gpu;
};
