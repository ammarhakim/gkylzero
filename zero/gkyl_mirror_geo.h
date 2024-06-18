#pragma once

#include <stdbool.h>

#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <math.h>
#include <string.h>
#include <gkyl_evalf_def.h>
#include <gkyl_gk_geometry.h>

// Object type
typedef struct gkyl_mirror_geo gkyl_mirror_geo;


// Some cumulative statistics
struct gkyl_mirror_geo_stat {
  long nquad_cont_calls; // num calls from quadrature
  long nroot_cont_calls; // num calls from root-finder
};  

typedef   void (*plate_func)(double s, double* RZ);

struct gkyl_mirror_geo {
  struct gkyl_efit* efit;

  struct gkyl_rect_grid rzgrid; // RZ grid on which psi(R,Z) is defined
  struct gkyl_range rzlocal; // local range over which psiRZ is defined
  struct gkyl_range rzlocal_ext; // extended range
  struct gkyl_basis rzbasis; // basis functions for R,Z grid
  int num_rzbasis; // number of basis functions in RZ
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  const struct gkyl_array *psibyrRZ; // psi(R,Z)/R DG representation
  const struct gkyl_array *psibyr2RZ; // psi(R,Z)/R^2 DG representation
                   
  struct gkyl_rect_grid fgrid; // flux grid for fpol
  struct gkyl_range frange; // flux range
  struct gkyl_range frange_ext; // extended range
  struct gkyl_basis fbasis; // psi basis for fpol
  const struct gkyl_array *fpoldg; // fpol(psi) dg rep
  const struct gkyl_array *qdg; // q(psi) dg rep
                                   
  double psisep; // psi of separatrix
  double zmaxis; // z of magnetic axis

  bool plate_spec;
  plate_func plate_func_lower;
  plate_func plate_func_upper;

  struct { int max_iter; double eps; } root_param;
  struct { int max_level; double eps; } quad_param;

  // pointer to root finder (depends on polyorder)
  struct RdRdZ_sol (*calc_roots)(const double *psi, double psi0, double Z,
    double xc[2], double dx[2]);

  struct gkyl_mirror_geo_stat stat; 
  struct gkyl_array* mc2p_nodal_fd;
  struct gkyl_range* nrange;
  double* dzc;
};




// Inputs to create a new GK geometry creation object
struct gkyl_mirror_geo_efit_inp {
  // Inputs to get psiRZ and related inputs from efit
  char* filepath;
  int rzpoly_order;
  enum gkyl_basis_type rz_basis_type;
  int fluxpoly_order;
  // Specifications for divertor plate
  bool plate_spec;
  plate_func plate_func_lower;
  plate_func plate_func_upper;

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
struct gkyl_mirror_geo_grid_inp {
  struct gkyl_rect_grid cgrid;
  struct gkyl_basis cbasis;
  
  double rclose; // closest R to discrimate
  double rright; // closest R to discrimate
  double zmin, zmax; // extents of Z for integration

  bool write_node_coord_array; // set to true if nodal coordinates should be written
  const char *node_file_nm; // name of nodal coordinate file

  double nonuniform_mapping_fraction; // Zero is uniform mapping, one is fully nonuniform mapping. In between values 
};

// An evaluator object for the non-uniform mapping to be read by initial conditions
// compuational to field aligned context object
struct gkyl_mirror_geo_c2fa_ctx {
  struct gkyl_rect_grid grid;
  struct gkyl_rect_grid cgrid;
  struct gkyl_range range;
  struct gkyl_range crange;
  struct gkyl_range crange_global;
  struct gkyl_basis basis;
  struct gkyl_basis cbasis;
  struct gkyl_array* c2fa; // nonuniform map
  struct gkyl_array* c2fa_deflate; // nonuniform map

  struct gkyl_rect_grid grid_deflate;
  struct gkyl_range range_deflate;
  struct gkyl_range range_global_deflate;
  struct gkyl_basis basis_deflate;
};

/**
 * Updater for advancing the map from computational coordinates to non-uniform coordinates
 * 
 * @param t Time
 * @param xn Computational coordinates as given in deflated geometry
 * @param fout Non-uniform coordinates in full 3D field alligned coordinates
 * @param ctx Context for the map
 */
void gkyl_mirror_geo_comp2fieldalligned_advance(double t, const double *xn, double *fout, void *ctx);
// gkyl_mirror_geo_comp2fieldalligned_advance

/**
 * Create new updater to compute the geometry (mapc2p) needed in GK
 * simulations.
 *
 * @param inp Input parameters
 * @param New GK geometry updater
 */
gkyl_mirror_geo *gkyl_mirror_geo_new(const struct gkyl_mirror_geo_efit_inp *inp);

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
int gkyl_mirror_geo_R_psiZ(const gkyl_mirror_geo *geo, double psi, double Z, int nmaxroots,
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
double gkyl_mirror_geo_integrate_psi_contour(const gkyl_mirror_geo *geo, double psi,
  double zmin, double zmax, double rclose);

/**
 * Compute physical coordinates (mapc2p)  given computational coordinates
 *
 * @param geo Geometry object
 * @param xn computational coordinates
 * @param ret physical coordinates
 */
void gkyl_mirror_geo_mapc2p(const gkyl_mirror_geo *geo, const struct gkyl_mirror_geo_grid_inp *inp,
    const double *xn, double *ret);

/**
 * Compute geometry (mapc2p) on a specified computational grid. The
 * output array must be pre-allocated by the caller.
 *
 * @param geo Geometry object
 * @param ginp Input structure for creating mapc2p
 * @param mapc2p On output, the DG representation of mapc2p
 */
void gkyl_mirror_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
  evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void *bmag_ctx, 
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2p, bool nonuniform,
  struct gkyl_array* map_arcL_nodal_fd, struct gkyl_array* map_arcL_nodal, struct gkyl_array* c2fa);

/**
 * Return cumulative statistics from geometry computations
 *
 * @param geo Geometry object
 * @return Cumulative statistics
 */
struct gkyl_mirror_geo_stat gkyl_mirror_geo_get_stat(const gkyl_mirror_geo *geo);

/**
 * Delete updater.
 *
 * @param geo Geometry object to delete
 */
void gkyl_mirror_geo_release(gkyl_mirror_geo *geo);

struct gkyl_range* gkyl_mirror_geo_get_nrange(gkyl_mirror_geo* geo);
struct gkyl_array* gkyl_mirror_geo_get_mc2p_nodal_fd(gkyl_mirror_geo* geo);
double* gkyl_mirror_geo_get_dzc(gkyl_mirror_geo* geo);
