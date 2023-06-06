#pragma once

// Object type
typedef struct gkyl_gkgeom gkyl_gkgeom;

// Inputs to create a new GK geometry creation object
struct gkyl_gkgeom_inp {
  // psiRZ and related inputs  
  const struct gkyl_rect_grid *rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_basis *rzbasis; // basis functions for R,Z grid
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  const struct gkyl_range *rzlocal; // local range over which psiRZ is defined

  // Parameters for root finder: leave unset to use defaults
  struct {
    int max_iter; // typically 20
    double eps; // typically 1e-10
  } root_param;

  // Parameters for nmumerical quadrature: leave unset to use default
  struct {
   int max_levels; // typically 6-7    
    double eps; // typically 1e-10
  } quad_param;
};

// Some cumulative statistics
struct gkyl_gkgeom_stat {
  long num_roots; // number of root computations
};  

/**
 * Create new updater to compute the geometry (mapc2p) needed in GK
 * simulations.
 *
 * @param inp Input parameters
 * @param New GK geometry updater
 */
gkyl_gkgeom *gkyl_gkgeom_new(const struct gkyl_gkgeom_inp *inp);

/**
 * Integrate along a specified psi countour and return its length. The
 * contour must lie completely inside the RZ domain of the psiRZ DG
 * field. The @a rclose parameter is used to select amongst the
 * multiple possible countours with the same psi. Foe example, to
 * select a flux surface on the outboard side of a double-null
 * configuration choose rclose to be Rmax.
 *
 * @param geo Geometry object
 * @param psi Psi value of contour
 * @param zmin Starting z location
 * @param zmax Ending z location
 * @param rclose Value of radial coordinate to discrimate between multiple
 *    contours
 * @return Length of contour
 */
double gkyl_gkgeom_integrate_psi_contour(const gkyl_gkgeom *geo, double psi,
  double zmin, double zmax, double rclose);

/**
 * Return cumulative statistics from geometry computations
 *
 * @param geo Geometry object
 * @return Cumulative statistics
 */
struct gkyl_gkgeom_stat gkyl_gkgeom_get_stat(const gkyl_gkgeom *geo);

/**
 * Delete updater.
 *
 * @param geo Geometry object to delete
 */
void gkyl_gkgeom_release(gkyl_gkgeom *geo);
