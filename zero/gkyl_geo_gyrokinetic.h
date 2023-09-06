#pragma once

#include <stdbool.h>

// Object type
typedef struct gkyl_geo_gyrokinetic gkyl_geo_gyrokinetic;

// Type of flux surface
enum gkyl_geo_gyrokinetic_type {
  GKYL_SOL_DN, // SOL of double-null configuration
  GKYL_SOL_SN, // SOL of single-null configuration
  GKYL_PF, // Private flux region
  GKYL_CORE // Core (closed flux-surface)
};  

// Inputs to create a new GK geometry creation object
struct gkyl_geo_gyrokinetic_inp {
  // psiRZ and related inputs  
  const struct gkyl_rect_grid *rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_basis *rzbasis; // basis functions for R,Z grid
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  const struct gkyl_range *rzlocal; // local range over which psiRZ is defined
  double B0; // Toroidal Field on axis
  double R0; // Axis

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

// Inputs to create geometry for a specific computational grid
struct gkyl_geo_gyrokinetic_geo_inp {
  const struct gkyl_rect_grid *cgrid;
  const struct gkyl_basis *cbasis;

  enum gkyl_geo_gyrokinetic_type ftype; // type of geometry
  
  double rclose; // closest R to discrimate
  double zmin, zmax; // extents of Z for integration

  bool write_node_coord_array; // set to true if nodal coordinates should be written
  const char *node_file_nm; // name of nodal coordinate file
};

// Some cumulative statistics
struct gkyl_geo_gyrokinetic_stat {
  long nquad_cont_calls; // num calls from quadrature
  long nroot_cont_calls; // num calls from root-finder
};  

/**
 * Create new updater to compute the geometry (mapc2p) needed in GK
 * simulations.
 *
 * @param inp Input parameters
 * @param New GK geometry updater
 */
gkyl_geo_gyrokinetic *gkyl_geo_gyrokinetic_new(const struct gkyl_geo_gyrokinetic_inp *inp);

/**
 * Get R(psi,Z) for a specified psi and Z value. Multiple values may
 * be returned (or none). The R(psi,Z) and dR/dZ are stored in the R
 * and dR arrays which be allocated by the caller.
 *
 * @param geo Geometry object
 * @param psi Psi value
 * @param Z Z value
 * @param nmaxroots Maximum number of roots
 * @param R on output, R(psi,Z)
 * @param dR on output, dR/dZ
 */
int gkyl_geo_gyrokinetic_R_psiZ(const gkyl_geo_gyrokinetic *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR);

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
double gkyl_geo_gyrokinetic_integrate_psi_contour(const gkyl_geo_gyrokinetic *geo, double psi,
  double zmin, double zmax, double rclose);

/**
 * Compute physical coordinates (mapc2p)  given computational coordinates
 *
 * @param geo Geometry object
 * @param xn computational coordinates
 * @param ret physical coordinates
 */
void gkyl_geo_gyrokinetic_mapc2p(const gkyl_geo_gyrokinetic *geo, const struct gkyl_geo_gyrokinetic_geo_inp *inp,
    const double *xn, double *ret);

/**
 * Compute geometry (mapc2p) on a specified computational grid. The
 * output array must be pre-allocated by the caller.
 *
 * @param geo Geometry object
 * @param ginp Input structure for creating mapc2p
 * @param mapc2p On output, the DG representation of mapc2p
 */
void gkyl_geo_gyrokinetic_calcgeom(const gkyl_geo_gyrokinetic *geo,
  const struct gkyl_geo_gyrokinetic_geo_inp *ginp, struct gkyl_array *mapc2p, struct gkyl_range *conversion_range);

/**
 * Return cumulative statistics from geometry computations
 *
 * @param geo Geometry object
 * @return Cumulative statistics
 */
struct gkyl_geo_gyrokinetic_stat gkyl_geo_gyrokinetic_get_stat(const gkyl_geo_gyrokinetic *geo);

/**
 * Delete updater.
 *
 * @param geo Geometry object to delete
 */
void gkyl_geo_gyrokinetic_release(gkyl_geo_gyrokinetic *geo);
