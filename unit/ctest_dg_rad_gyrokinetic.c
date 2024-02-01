// A test for the vlasov LBO collision updater
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <math.h>


struct gk_rad_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double c_s; // sound speed
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double B0; // reference magnetic field
  double n0; // reference density
  // Simulation parameters
  double Lx; // Box size in x
  double kperp; // perpendicular wave number used in Poisson solve
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double finalTime; // end time
  int numFrames; // number of output frames
};


void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_rad_ctx *app = ctx;
  fout[0] = app->B0;
}

struct gk_rad_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 3*GKYL_PROTON_MASS; // (Lithium) ion mass
  double me = GKYL_ELECTRON_MASS; // electron mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double Te = 30.3027*eV;//7.003256*eV;//
  double Ti = Te;//30.01373*eV;
  double B0 = 1.0; // Magnetic field magnitude in Tesla
  double n0 = 1.0e19; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nuFrac = 0.1;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq
  
  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lx = 100.0*rho_s;

  // Perpendicular wavenumber in SI units:
  double kperpRhos = 0.1;
  double kperp = kperpRhos / rho_s;

  double vpar_max_elc = 4.0*vtElc;
  double mu_max_elc = 0.75*me*(4.0*vtElc)*(4.0*vtElc)/(2.0*B0);

  double vpar_max_ion = 4.0*vtIon;
  double mu_max_ion = 0.75*mi*(4.0*vtIon)*(4.0*vtIon)/(2.0*B0);

  double finalTime = 6.0e-8; // Because it crashes at about 5e-8 
  double numFrames = 10;

  struct gk_rad_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .B0 = B0, 
    .n0 = n0, 
    .Lx = Lx, 
    .kperp = kperp, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .finalTime = finalTime, 
    .numFrames = numFrames, 
  };
  return ctx;
}



// allocate double array (filled with zeros)
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
mapc2p_1x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; 
}

void
bmag_func_1x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0]; 
  fout[0] = 1.0;
}

void
mapc2p_2x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; 
}

void
bmag_func_2x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1];
  fout[0] = 1.0;
}

void
mapc2p_3x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func_3x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 29.999913327161483*GKYL_ELEMENTARY_CHARGE;
}

void
test_1x(int poly_order, int vdim, bool use_gpu)
{
  double mass = 1.0;
  double charge = -1.0*GKYL_ELEMENTARY_CHARGE;
  int cdim = 1;
  int ndim = cdim + vdim;
  double lower[ndim], upper[ndim];
  int cells[ndim];
  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];

  struct gk_rad_ctx ctx = create_ctx(); // context for init functions

  // Phase space and Configuration space extents and resolution
  lower[0] = 0.0;
  upper[0] = 1.0;
  cells[0] = 2;
  lower[1] = -ctx.vpar_max_elc;
  upper[1] = ctx.vpar_max_elc;
  cells[1] = 32;
  if (vdim == 2) {
    lower[2] = 0.0;
    upper[2] = ctx.mu_max_elc;
    cells[2] = 16;
  }
  confLower[0] = lower[0]; 
  confUpper[0] = upper[0];
  confCells[0] = cells[0];

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, surf_vpar_basis, surf_mu_basis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
    if (vdim>1)
      gkyl_cart_modal_serendip(&surf_mu_basis, ndim-1, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
    if (vdim>1) // constant mu surface
      gkyl_cart_modal_gkhybrid(&surf_mu_basis, cdim, poly_order);
  }
  // constant vpar surface
  gkyl_cart_modal_serendip(&surf_vpar_basis, ndim-1, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  printf("Next: initialize geo\n");
  // Initialize geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&confGrid, &confLocal, &confLocal_ext, &confBasis, 
    mapc2p_1x, 0, bmag_func_1x, 0, use_gpu);
  /*  struct gk_geometry *gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&confGrid, &confLocal, &confLocal_ext, &confBasis, 
    mapc2p, &ctx, bmag_func, &ctx, use_gpu);
  printf("Initialized geo\n");
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &confGrid, &confLocal, &confLocal_ext, &confBasis, use_gpu);*/
  printf("Geometry deflated\n");
  // allocate drag coefficients in vparallel and mu for each collision
  // vnu = v_par*nu(v)
  // vsqnu = 2*mu*nu(v)
  // Note that through the spatial variation of B, both these drag coefficients depend on the full phase space
  struct gkyl_array *vnu, *vsqnu, *vnu_surf, *vsqnu_surf, *nvnu, *nvsqnu, *nvnu_surf, *nvsqnu_surf,
    *nvnu_host, *nvsqnu_host, *nvnu_surf_host, *nvsqnu_surf_host;
  vnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vsqnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vnu_surf = mkarr(use_gpu, surf_vpar_basis.num_basis, local_ext.volume);
  vsqnu_surf = mkarr(use_gpu, surf_mu_basis.num_basis, local_ext.volume);

  double a, alpha, beta, gamma, v0;
  a = 0.153650876536253;
  alpha = 8000.006932403581;
  beta = 0.892102642790662;
  gamma = -3.923194017288736;
  v0 = 3.066473173090881;
  printf("Next: call calc_gk_rad_vars\n");
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&grid, &confBasis, &basis, 
		  charge, mass, gk_geom, a, alpha, beta, gamma, v0, use_gpu);

  gkyl_dg_calc_gk_rad_vars_nu_advance(calc_gk_rad_vars, &confLocal, &local, vnu_surf, vnu, vsqnu_surf, vsqnu);

  nvnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvsqnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvnu_surf = mkarr(use_gpu, surf_vpar_basis.num_basis, local_ext.volume);
  nvsqnu_surf = mkarr(use_gpu, surf_mu_basis.num_basis, local_ext.volume);
  
  nvnu_host = nvnu;
  nvsqnu_host = nvsqnu;
  nvnu_surf_host = nvnu_surf;
  nvsqnu_surf_host = nvsqnu_surf;
  if (use_gpu) {
    nvnu_host = mkarr(false, basis.num_basis, local_ext.volume);
    nvsqnu_host = mkarr(false, basis.num_basis, local_ext.volume);
    nvnu_surf_host = mkarr(false, surf_vpar_basis.num_basis, local_ext.volume);
    nvsqnu_surf_host = mkarr(false, surf_mu_basis.num_basis, local_ext.volume);
  }
  struct gkyl_dg_rad_gyrokinetic_auxfields drag_inp = { .nvnu_surf = nvnu_surf, .nvnu = nvnu,
							     .nvsqnu_surf = nvsqnu_surf, .nvsqnu = nvsqnu};

  printf("created auxfields\n");
  if (use_gpu) {
    gkyl_array_copy(nvnu, nvnu_host);
    gkyl_array_copy(nvsqnu, nvsqnu_host);
    gkyl_array_copy(nvnu_surf, nvnu_surf_host);
    gkyl_array_copy(nvsqnu_surf, nvsqnu_surf_host);
  }

  // Initialize distribution function with proj_gkmaxwellian_on_basis
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  // Project n, udrift, and vt^2 based on input functions
  struct gkyl_array *m0 = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_density, 0);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_upar, 0);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_temp_elc, 0);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal_ext, m0); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal_ext, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal_ext, vtsq);
  gkyl_array_scale(vtsq, 1.0/mass);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms = mkarr(false, 2*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms, 1.0, udrift, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, vtsq  , 1*confBasis.num_basis);

  // Initialize Maxwellian projection object
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
      &confBasis, &basis, poly_order+1, use_gpu);

  // If on GPUs, need to copy n, udrift, and vt^2 onto device
  struct gkyl_array *prim_moms_dev, *m0_dev;
  if (use_gpu) {
    prim_moms_dev = mkarr(use_gpu, 2*confBasis.num_basis, confLocal_ext.volume);
    m0_dev = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);

    gkyl_array_copy(prim_moms_dev, prim_moms);
    gkyl_array_copy(m0_dev, m0);
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, m0_dev, prim_moms_dev,
      gk_geom->bmag, gk_geom->bmag, mass, f);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, m0, prim_moms,
      gk_geom->bmag, gk_geom->bmag, mass, f);
  }

  // initialize solver 
  struct gkyl_dg_updater_collisions *slvr;
  slvr = gkyl_dg_updater_rad_gyrokinetic_new(&grid, &confBasis, &basis, &local, &drag_inp, use_gpu);
  
  struct gkyl_array *cflrate, *rhs, *fin, *nI, *fmax;
  cflrate = mkarr(use_gpu, 1, local_ext.volume);
  rhs = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  struct gkyl_array *m0_ion = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  struct gkyl_dg_updater_moment *m0_ion_up = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis, &confLocal, 0, GKYL_PROTON_MASS, gk_geom, "M0", false, use_gpu);
  gkyl_dg_updater_moment_gyrokinetic_advance(m0_ion_up, &local, &confLocal, f, m0_ion);
  printf("Before hyper dg\n");
  // run hyper_dg_advance
  int nrep = 1;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_array_clear(nvnu_surf, 0.0);
    gkyl_array_clear(nvnu, 0.0);
    gkyl_array_clear(nvsqnu_surf, 0.0);
    gkyl_array_clear(nvsqnu, 0.0);

    gkyl_dg_calc_gk_rad_vars_nI_nu_advance(calc_gk_rad_vars, 
      &confLocal, &local, vnu_surf, vnu, vsqnu_surf, vsqnu, 
      m0_ion, nvnu_surf, nvnu, nvsqnu_surf, nvsqnu);


    gkyl_dg_updater_rad_gyrokinetic_advance(slvr, &local, f, cflrate, rhs);
    printf("After updater advance, n=%i\n",n);
  }

  // Take 2nd moment of rhs to find energy loss
  struct gkyl_dg_updater_moment *m2_calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, 0, GKYL_ELECTRON_MASS, gk_geom, "M2", false, use_gpu);
  struct gkyl_array *m2_final = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_dg_updater_moment_gyrokinetic_advance(m2_calc, &local, &confLocal, rhs, m2_final);

  double *m00 = gkyl_array_fetch(m0, 0+ghost[0]);
  double *m20 = gkyl_array_fetch(m2_final, 0+ghost[0]);
  //double *m21 = gkyl_array_fetch(rhs, 0+ghost[0]);
  double *m00_ion = gkyl_array_fetch(m0_ion, 0+ghost[0]);
  
  //  double cell_avg0 = m20[0]/pow(sqrt(2),cdim);
  double cell_avg0 = 1.0/2.0*GKYL_ELECTRON_MASS*m20[0]/(m00[0]*m00_ion[0]);

  double correct = 4.419192427285379e-32;
  //  for (int i=0; i<30; i++){
  printf("cell_avg=%e, correct energy=%e, density=%.10e, nI=%e, m2=%e\n",cell_avg0, correct, m00[0], m00_ion[0], m20[0]);
    //}
  
  TEST_CHECK( gkyl_compare( correct*1e30, cell_avg0*1e30, 1e-12));
  TEST_CHECK( cell_avg0>0);

  // Release memory
  gkyl_array_release(vnu);
  gkyl_array_release(vsqnu);
  if (use_gpu) {
    gkyl_array_release(nvnu_host);
    gkyl_array_release(nvsqnu_host);      
  }
  gkyl_dg_calc_gk_rad_vars_release(calc_gk_rad_vars);

  gkyl_array_release(m0);
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_array_release(prim_moms);
  if (use_gpu) {
    gkyl_array_release(m0_dev);
    gkyl_array_release(prim_moms_dev);      
  }
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_proj_maxwellian_on_basis_release(proj_max);  

  gkyl_array_release(f);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_dg_updater_rad_gyrokinetic_release(slvr);
  //gkyl_dg_updater_moment_gyrokinetic_release(m0_calc);
  gkyl_dg_updater_moment_gyrokinetic_release(m2_calc);
  //gkyl_array_release(m0_final);
  gkyl_array_release(m2_final);
}

struct eval_nu_ctx {
  double Abar;
  double alpha;
  double beta;
  double V0; 
  double gamma;
  double Crad;
  double B0;
};

void eval_nu_1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1], mu = xn[2];

  struct eval_nu_ctx *radFit = ctx;
  double Abar  = radFit->Abar ;
  double alpha = radFit->alpha;
  double beta  = radFit->beta ;
  double V0    = radFit->V0   ; 
  double gamma = radFit->gamma;
  double Crad  = radFit->Crad ;
  double B0    = radFit->B0   ;
  double me = GKYL_ELECTRON_MASS;
  double eV = GKYL_ELEMENTARY_CHARGE;

  double vSq = pow(vpar,2)+2.*mu*B0/me;
  double denomFac = (me/(2.*eV*V0))*vSq;

  fout[0] = vSq > 1e-18? (Abar/Crad)*(alpha+beta)*pow(vSq,gamma/2.) / 
    ( beta*pow(denomFac,-alpha/2.)+alpha*pow(denomFac,beta/2.) ) 
    : 0.;
}

void
projnu_1x2v(int poly_order, bool use_gpu)
{
  // Project nu onto an (x,vpar,mu) grid. Use eval_on_nodes now
  // though we'll have to write an updater for it if we keep the
  // dependence on B.
  int cdim = 1;
  int vdim = 2;
  int ndim = cdim + vdim;

  double me = GKYL_ELECTRON_MASS;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double qe = -eV;
  double Te = 100*2.8*eV;
  double B0 = 4.2;

  // H fitting parameters.
  double Abar = 0.166594984602572;
  double alpha = 8.000072245234049e+03;
  double beta = 0.875480342166318;
  double V0 = 3.071566470338830;
  double gamma = -3.956239523465181;
  double Crad = 8.*sqrt(M_PI)*pow(GKYL_ELEMENTARY_CHARGE,5./2.)/me;
  struct eval_nu_ctx radFit = {
    .Abar  = Abar ,
    .alpha = alpha,
    .beta  = beta ,
    .V0    = V0   ,
    .gamma = gamma,
    .Crad  = Crad ,
//    .B0    = 0.5*B0   ,
    .B0    = B0   ,
//    .B0    = 2.*B0   ,
  };

  double vtElc = sqrt(Te/me);
  double lower[] = {-M_PI, -4.0*vtElc, 0.};
  double upper[] = { M_PI,  4.0*vtElc, 0.5*me*pow(4.0*vtElc,2)/B0};
  int cells[] = {2, 32, 16};

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d]; 
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // projection nu
  struct gkyl_array *nu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_eval_on_nodes *evonnod = gkyl_eval_on_nodes_new(&grid, &basis, 1, eval_nu_1x2v, &radFit);

  gkyl_eval_on_nodes_advance(evonnod, 0.0, &local, nu);

  gkyl_grid_sub_array_write(&grid, &local, nu, "ctest_dg_rad_gyrokinetic_nu.gkyl");

  // Release memory
  gkyl_array_release(nu);
  gkyl_eval_on_nodes_release(evonnod);
}

void test_1x1v_p1() { test_1x(1, 1, false); }
void test_1x2v_p1() { test_1x(1, 2, false); }
void test_1x1v_p2() { test_1x(2, 1, false); }
void test_1x2v_p2() { test_1x(2, 2, false); }

void test_projnu_1x2v_p1() { projnu_1x2v(1, false); }

// #ifdef GKYL_HAVE_CUDA

// void test_1x1v_p1_gpu() { test_1x(1, 1, true); }
// void test_1x2v_p1_gpu() { test_1x(1, 2, true); }
// void test_1x1v_p2_gpu() { test_1x(2, 1, true); }
// void test_1x2v_p2_gpu() { test_1x(2, 2, true); }

// #endif

TEST_LIST = {
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x2v_p1", test_1x2v_p1 },
  { "test_1x1v_p2", test_1x1v_p2 },
  { "test_1x2v_p2", test_1x2v_p2 },
  { "test_projnu_1x2v_p1", test_projnu_1x2v_p1 },

// #ifdef GKYL_HAVE_CUDA
//   { "test_1x1v_p1_gpu", test_1x1v_p1_gpu },
//   { "test_1x2v_p1_gpu", test_1x2v_p1_gpu },
//   { "test_1x1v_p2_gpu", test_1x1v_p2_gpu },
//   { "test_1x2v_p2_gpu", test_1x2v_p2_gpu },

// #endif
  { NULL, NULL },
};
