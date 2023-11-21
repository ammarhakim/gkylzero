#pragma once

#include <stdbool.h>

#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

#include <math.h>
#include <string.h>

// Object type
typedef struct gkyl_geo_gyrokinetic gkyl_geo_gyrokinetic;


// Some cumulative statistics
struct gkyl_geo_gyrokinetic_stat {
  long nquad_cont_calls; // num calls from quadrature
  long nroot_cont_calls; // num calls from root-finder
};  

struct gkyl_geo_gyrokinetic {
  const struct gkyl_rect_grid* rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  const struct gkyl_range* rzlocal; // local range over which psiRZ is defined
  int num_rzbasis; // number of basis functions in RZ

  struct { int max_iter; double eps; } root_param;
  struct { int max_level; double eps; } quad_param;

  // pointer to root finder (depends on polyorder)
  struct RdRdZ_sol (*calc_roots)(const double *psi, double psi0, double Z,
    double xc[2], double dx[2]);

  struct gkyl_geo_gyrokinetic_stat stat; 
  double B0;
  double R0;
  struct gkyl_array* mc2p_nodal_fd;
  struct gkyl_range* nrange;
  double* dzc;
};



// Type of flux surface
enum gkyl_geo_gyrokinetic_type {
  GKYL_SOL_DN_OUT, // Outboard SOL of double-null configuration
  GKYL_SOL_DN_IN, // Inboard SOL of double-null configuration
  GKYL_SOL_SN_LO, // SOL of a lower single-null configuration
  GKYL_SOL_SN_UP, // SOL of an upper single-null configuration
  GKYL_PF_UP, // Private flux region at top
  GKYL_PF_LO, // Private flux region at bottom
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
  int* bcs;
  const struct gkyl_basis *cbasis;

  enum gkyl_geo_gyrokinetic_type ftype; // type of geometry
  
  double rclose; // closest R to discrimate
  double rleft; // closest R to discrimate
  double rright; // closest R to discrimate
  double zmin, zmax; // extents of Z for integration
  double zmin_left, zmin_right; // for lower single null and PF cases diff b/t in and outboard side
  double zmax_left, zmax_right; // for upper single null and PF cases diff b/t in and outboard side
  double zmaxis; // z of magnetic axis

  bool write_node_coord_array; // set to true if nodal coordinates should be written
  const char *node_file_nm; // name of nodal coordinate file
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
void gkyl_geo_gyrokinetic_calcgeom(gkyl_geo_gyrokinetic *geo,
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

struct gkyl_range* gkyl_geo_gyrokinetic_get_nrange(gkyl_geo_gyrokinetic* geo);
struct gkyl_array* gkyl_geo_gyrokinetic_get_mc2p_nodal_fd(gkyl_geo_gyrokinetic* geo);
double* gkyl_geo_gyrokinetic_get_dzc(gkyl_geo_gyrokinetic* geo);
