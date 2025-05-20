#pragma once

#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Forward declare internal private object
struct gkyl_mirror_grid_gen_x;

struct gkyl_mirror_grid_gen {
  struct gkyl_array *nodesrz; // r,z coordinates of corner nodes of cells

  struct gkyl_mirror_grid_gen_x *gg_x;
};  

// flag to indicate what field-line coordinate to use
enum gkyl_mirror_grid_gen_field_line_coord {
  GKYL_MIRROR_GRID_GEN_PSI_CART_Z, // use psi and Cartesian Z coordinate
  GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // use sqrt(psi) and Cartesian Z coordinate
  GKYL_MIRROR_GRID_GEN_FL_LENGTH, // use field-line length NYI
};  

// input struct to construct the mirror geometry
struct gkyl_mirror_grid_gen_inp {
  const struct gkyl_rect_grid *comp_grid; // Computational space grid (psi, phi, z)
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // field-line coordinate to use
  bool include_axis; // add nodes on r=0 axis (the axis is assumed be psi=0)
  
  double R[2], Z[2]; // extents of R,Z grid on which psi(R,Z) is given
  int nrnodes, nznodes; // numer of nodes in R and Z
  const struct gkyl_array *psiRZ; // nodal values of psi(R,Z)

  bool write_psi_cubic; // set to true to write the cubic fit to file
  const char *psi_cubic_fname; // name for cubic fit file
};

/**
 * Create new grid generator for mirror geometries. Return a NULL
 * pointer if the grid generator failed.
 *
 * @param inp Input parameters to the grid generator
 * @return newly create mirror geometry object
 */
struct gkyl_mirror_grid_gen *gkyl_mirror_grid_gen_inew(const struct gkyl_mirror_grid_gen_inp *inp);

/**
 * Release the mirror geometry object.
 *
 * @param geom Geometry object to release
 */
void gkyl_mirror_grid_gen_release(struct gkyl_mirror_grid_gen *geom);
