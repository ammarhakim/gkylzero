#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>
#include <mpi.h>

void shaped_pfunc_lower_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = -8.29;
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = 8.29;
}

void shaped_pfunc_upper_inner(double s, double* RZ){
    RZ[0] = 1.651 + (1.8 - 1.651)*s;
    RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
    RZ[0] = 1.65 + (1.8 - 1.65)*s;
    RZ[1] = -(6.33 + (6.777 - 6.33)*s);
}

struct gkyl_block_geom*
create_block_geom(void)
{
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 12);

  /* Block layout and coordinates

   x  
   ^  
   |
   4  +------------------+------------------+------------------+
   |  |b1                |b2                |b3                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
   |  |                  |                  |                  |
   3  +------------------+------------------+------------------+
   |  |b0               x|o b10            %|$ b4              |
   |  |lower outer PF   x|o outer core     %|$ upper outer PF  |
   |  |                 x|o                %|$                 |
   |  +------------------+------------------+------------------+
   2  +------------------+------------------+------------------+
   |  |b9               x|o b11            %|$ b5              |
   |  |lower inner PF   x|o inner core     %|$ upper inner PF  |
   |  |                 x|o                %|$                 |
   1  +------------------+------------------+------------------+
   |  |b8                |b7                |b6                |
   |  |lower inner SOL   |middle inner SOL  |upper inner SOL   |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+

      0 -----------1------------2------------3 -> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  

  double psisep = 1.5093065418975686;
  double Zxpt_lo = -6.1672666854902927;
  double Zxpt_up = 6.1672666854902927;


  double psi_up_core = 1.8;
  double psi_up_pf = 1.8;
  double psi_lo_outer_sol = 0.934;
  double psi_lo_inner_sol = 1.45;

  int npsi_outer_sol = 4;
  int npsi_core = 4;
  int npsi_inner_sol = 4;
  int npsi_lower_pf = 4;
  int npsi_upper_pf = 4;

  double ntheta_lower  = 8;
  double ntheta_middle = 8;
  double ntheta_upper  = 8;

  double theta_lo = -M_PI + 1e-14, theta_up = M_PI - 1e-14;

  struct gkyl_tok_geo_efit_inp efit_inp_inner= {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
      .rzpoly_order = 2,
      .rz_basis_type = GKYL_BASIS_MODAL_TENSOR,
      .fluxpoly_order = 1,
      .plate_spec = true,
      .quad_param = {  .eps = 1e-10 },
      .reflect = true,
      .plate_func_lower = shaped_pfunc_lower_inner,
      .plate_func_upper = shaped_pfunc_upper_inner,
    };
  
  struct gkyl_tok_geo_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
      .rzpoly_order = 2,
      .rz_basis_type = GKYL_BASIS_MODAL_TENSOR,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 },
      .reflect = true,
    };
  
  struct gkyl_tok_geo_efit_inp efit_inp_outer = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
      .rzpoly_order = 2,
      .rz_basis_type = GKYL_BASIS_MODAL_TENSOR,
      .fluxpoly_order = 1,
      .plate_spec = true,
      .quad_param = {  .eps = 1e-10 },
      .reflect = true,
      .plate_func_lower = shaped_pfunc_lower_outer,
      .plate_func_upper = shaped_pfunc_upper_outer,
    };

  // block 0. Lower outer PF region.
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_lower_pf, ntheta_lower },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_outer,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b0_nodes.gkyl",
          .ftype = GKYL_PF_LO_R,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 2.1,
          .rmax = 6.2,
          .zxpt_lo = Zxpt_lo,
          .zmin = -8.29,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 1. Lower outer SOL.
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_lower},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_outer,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b1_nodes.gkyl",
          .ftype = GKYL_SOL_DN_OUT_LO,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zxpt_lo = Zxpt_lo,
          .zmin = -8.29,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 2. Middle outer SOL.
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_outer,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b2_nodes.gkyl",
          .ftype = GKYL_SOL_DN_OUT_MID,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zxpt_lo = Zxpt_lo,
          .zxpt_up = Zxpt_up,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 10, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_UPPER_POSITIVE}
      }
    }
  );

  // block 3. Upper outer SOL.
  gkyl_block_geom_set_block(bgeom, 3, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_outer,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b3_nodes.gkyl",
          .ftype = GKYL_SOL_DN_OUT_UP,
          .rright = 6.2,
          .rleft = 1.1,
          .rmin = 2.1,
          .rmax = 6.2,
          .zxpt_up = Zxpt_up,
          .zmax = 8.29,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  // block 4. Upper outer PF region.
  gkyl_block_geom_set_block(bgeom, 4, &(struct gkyl_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_upper_pf, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_outer,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b4_nodes.gkyl",
          .ftype = GKYL_PF_UP_R,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 2.1,
          .rmax = 6.2,
          .zxpt_up = Zxpt_up,
          .zmax = 8.29,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );
  
  // block 5. Upper inner PF region.
  gkyl_block_geom_set_block(bgeom, 5, &(struct gkyl_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_upper_pf, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_inner,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b5_nodes.gkyl",
          .ftype = GKYL_PF_UP_L,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 1.6,
          .rmax = 6.2,
          .zxpt_up = Zxpt_up,
          .zmax = 6.34,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 6. Upper inner SOL.
  gkyl_block_geom_set_block(bgeom, 6, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_inner_sol, ntheta_upper},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_inner,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b6_nodes.gkyl",
          .ftype = GKYL_SOL_DN_IN_UP,
          .rleft = 2.0,
          .rright= 6.2,
          .rmin = 1.3,
          .rmax = 6.2,
          .zxpt_up = Zxpt_up,
          .zmax = 6.34,  
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 7, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 7. Middle inner SOL.
  gkyl_block_geom_set_block(bgeom, 7, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_inner,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b7_nodes.gkyl",
          .ftype = GKYL_SOL_DN_IN_MID,
          .rleft = 2.0,
          .rright= 6.2,
          .rmin = 1.3,
          .rmax = 6.2,
          .zxpt_up = Zxpt_up,
          .zxpt_lo = Zxpt_lo,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 11, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 8, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 8. Lower inner SOL.
  gkyl_block_geom_set_block(bgeom, 8, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_outer_sol, ntheta_lower},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_inner,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b8_nodes.gkyl",
          .ftype = GKYL_SOL_DN_IN_LO,
          .rleft = 2.0,
          .rright= 6.2,
          .rmin = 1.3,
          .rmax = 6.2,
          .zmin = -6.34,
          .zxpt_lo = Zxpt_lo,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 9, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  // block 9. Lower inner PF region.
  gkyl_block_geom_set_block(bgeom, 9, &(struct gkyl_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_lower_pf, ntheta_lower },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp_inner,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b9_nodes.gkyl",
          .ftype = GKYL_PF_LO_L,
          .rright = 6.2,
          .rleft = 2.0,
          .rmin = 1.6,
          .rmax = 6.2,
          .zxpt_lo = Zxpt_lo,
          .zmin = -6.34,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 8, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  // block 10. Outer core.
  gkyl_block_geom_set_block(bgeom, 10, &(struct gkyl_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_core, theta_up},
      .cells = { npsi_core, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b10_nodes.gkyl",
          .ftype = GKYL_CORE_R,
          .rclose = 6.2,
          .rleft= 1.1,
          .rright= 6.2,
          .rmin=1.5,
          .rmax=6.2,
          .zxpt_lo = Zxpt_lo,
          .zxpt_up = Zxpt_up,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 11. Outer core.
  gkyl_block_geom_set_block(bgeom, 11, &(struct gkyl_block_geom_info) {
      .lower = { psisep, theta_lo},
      .upper = { psi_up_core, theta_up},
      .cells = { npsi_core, ntheta_middle},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .tok_efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .write_node_coord_array = true,
          .node_file_nm = "b11_nodes.gkyl",
          .ftype = GKYL_CORE_L,
          .rclose = 1.1,
          .rleft= 1.1,
          .rright= 6.2,
          .rmin=1.5,
          .rmax=6.2,
          .zxpt_lo = Zxpt_lo,
          .zxpt_up = Zxpt_up,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  return bgeom;
}

struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massAr; // Argon mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double vtIon;
  double vtElc;
  double vtAr;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double n0Ar; // Argon reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  // Simulation parameters
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_Ar; // Velocity space extents in vparallel for Ar
  double mu_max_Ar; // Velocity space extents in mu for Ar
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};



struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 100*2.8*eV;
  double Ti = 150*2.8*eV;
  double TAr = 40.0*eV;
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/2.8; // Particle density in 1/m^3
  double n0Ar = n0*0.0001/3.0; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(TAr/mAr);

  // Source parameters.
  double nsource = 3.9e23/2.8; // peak source rate in particles/m^3/s 
  double T_source = 285*eV*2.8;
  double cx = 0.0065612*9;
  double cz = 0.4916200;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 18*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 18*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_Ar = 4.0*vtAr;
  double mu_max_Ar = 18.*mAr*vtAr*vtAr/(2.0*B0);

  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 16;
  int Nmu = 8;

  double t_end = 5.0e-7; 
  double num_frames = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti, 
    .TAr = TAr, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .vtAr = vtAr,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .n0Ar = n0Ar,
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar, 
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

static void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app,
  double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }
    gkyl_gyrokinetic_multib_app_write(app, t_curr, frame);
  }
}

void
initDensity(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double n = input->n0;
  fout[0] = n;
}

void
initDensityImpurity(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = 1.0e5;
}


void
initDensityNeutral(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double n = input->n0Ar;
  fout[0] = n;
}

void
initTempElc(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->Te;
  fout[0] = T;
}

void
initTempIon(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->Ti;
  fout[0] = T;
}

void
initTempAr(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->TAr;
  fout[0] = T;
}

void
initUpar(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
initUDriftNeutral(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
sourceDensity(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double n = input->nsource;
  fout[0] = n;
}

void
sourceUpar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
sourceTemp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  double T = input->T_source;
  fout[0] = T;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuIon;
}


int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

struct gkyl_comm *comm = 0;
  if (app_args.use_mpi) {
#ifdef GKYL_HAVE_MPI
    MPI_Init(&argc, &argv);
        comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#endif
  }
  if (comm == 0)
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  // construct block geometry
  struct gkyl_block_geom *bgeom = create_block_geom();
  int nblocks = gkyl_block_geom_num_blocks(bgeom);


  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.
  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[1];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempElc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };


  struct gkyl_block_physical_bcs elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = {  ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .enforce_positivity = false,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 2,
      .collide_with = { "ion", "Ar1" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .radiation = {
      .radiation_id = GKYL_GK_RADIATION, 
      .num_cross_collisions = 1, 
      .collide_with = { "Ar1" },
      .z = 18,
      .charge_state = 1,
      .num_of_densities = 1, // Must be 1 for now
    },

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1", // ion is always the higher charge state
          .donor_nm = "Ar0", // interacts with elc to give up charge
          .charge_state = 0, // corresponds to lower charge state (donor)
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .recvr_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    }, 

    .duplicate_across_blocks = true,
    .blocks = elc_blocks,
    .num_physical_bcs = 20,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[1];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = sourceDensity,
        .ctx_upar = &ctx,
        .upar = sourceUpar,
        .ctx_temp = &ctx,
        .temp = sourceTemp,      
      }, 
    },

  };

  struct gkyl_block_physical_bcs ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = {  ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .enforce_positivity = false,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "Ar1" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 
  
    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 20,
    .bcs = ion_phys_bcs,
  };

  // Ar+1 Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb Ar1_blocks[1];
  Ar1_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .polarization_density = ctx.n0Ar,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensityImpurity,
      .ctx_upar = &ctx,
      .upar = initUpar,
      .ctx_temp = &ctx,
      .temp = initTempAr,
    },

  };

  struct gkyl_block_physical_bcs Ar1_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_GK_SHEATH},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = {  ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { cells_v[0], cells_v[1] },
    .no_by = true,
    .enforce_positivity = false,
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0Ar, // Density used to calculate coulomb logarithm
      .T_ref = ctx.TAr, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 2,
      .collide_with = { "elc", "ion" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .donor_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .recvr_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },

  
    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 20,
    .bcs = ion_phys_bcs,
  };

  // Neutral Ar0 Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_neut_species_pb Ar0_blocks[1];
  Ar0_blocks[0] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = initDensityNeutral,
      .ctx_udrift = &ctx,
      .upar = initUDriftNeutral,
      .ctx_temp = &ctx,
      .temp = initTempAr,
    },

  };

  struct gkyl_block_physical_bcs Ar0_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 10 BCs
    //{ .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
    // block 11 BCs
    //{ .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_FIXED_FUNC},
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_neut_species Ar0 = {
    .name = "Ar0",
    .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = {  ctx.vpar_max_Ar,  ctx.vpar_max_Ar,  ctx.vpar_max_Ar },
    .cells = { cells_v[0], cells_v[0], cells_v[0] },
    .is_static = true,
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},

    .duplicate_across_blocks = true,
    .blocks = Ar0_blocks,
    .num_physical_bcs = 20,
    .bcs = Ar0_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
  };

  struct gkyl_block_physical_bcs field_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_NONE},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_NONE},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_NONE },
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_NONE},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_NONE },
    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET },
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_NONE },
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET},
    // block 8 BCs
    { .bidx = 8, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET },
    { .bidx = 8, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_NONE},
    // block 9 BCs
    { .bidx = 9, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_DIRICHLET },
    { .bidx = 9, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_NONE},
    // block 10 BCs
    { .bidx = 10, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_NEUMANN},
    // block 11 BCs
    { .bidx = 11, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_FIELD_NEUMANN},
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 20,
    .bcs = field_phys_bcs,
  };



  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "multib_step_2x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,

    .block_geom = bgeom,
    .cfl_frac = 0.9,

    .num_species = 3,
    .species = { elc, ion, Ar1},

    .num_neut_species = 1,
    .neut_species = { Ar0 },

    .field = field,

    .comm = comm
  };

  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = 0.6;

  // Create trigger for IO.
  int num_frames = 4;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames };

  // Initialize simulation.
  gkyl_gyrokinetic_multib_app_apply_ic(app, t_curr);
  write_data(&io_trig, app, t_curr, false);

  gkyl_comm_release(comm);
  gkyl_block_geom_release(bgeom);
  gkyl_gyrokinetic_multib_app_release(app);
  

  finish:

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;



}
