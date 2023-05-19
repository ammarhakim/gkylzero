// Test the perpendicular FEM Helmholtz/Poisson solver.
//
#include <acutest.h>

#include <math.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <gkyl_fem_poisson_perp.h>

#define PERP_DIM 2

double poly_test_func_1x(double x, double a, double *c)
{
  // Function that can be used to produce homogeneous Dirichlet or Neumann
  // boundary values depending on the choice of a and c. It assumes x \in [0,1].
  return pow(x,2)/2.-a*pow(x,4)/12.+c[0]*x+c[1];
}

void evalFunc_dirichletx_dirichlety(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  double kz = 1.;
  double xp = x, yp = y, zp = z;
  fout[0] = -( (1 - b*pow(yp,2))*poly_test_func_1x(xp, a, c)
//              +(1 - a*pow(xp,2))*poly_test_func_1x(yp, b, d) )*sin(kz*z);
              +(1 - a*pow(xp,2))*poly_test_func_1x(yp, b, d) )*(1.+kz*z);
//              +(1 - a*pow(xp,2))*poly_test_func_1x(yp, b, d) )*(1.+kz*z+0.5*pow(z,2));
}
void evalFunc_dirichletx_dirichlety_sol(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  double a = 2.;
  double c[] = {a/12.-1./2., 0.};
  double b = 2.;
  double d[] = {b/12.-1./2., 0.};
  double kz = 1.;
  double xp = x, yp = y, zp = z;
//  fout[0] = poly_test_func_1x(xp, a, c)*poly_test_func_1x(yp, b, d)*sin(kz*z);
  fout[0] = poly_test_func_1x(xp, a, c)*poly_test_func_1x(yp, b, d)*(1.+kz*z);
//  fout[0] = poly_test_func_1x(xp, a, c)*poly_test_func_1x(yp, b, d)*(1.+kz*z+0.5*pow(z,2));
}

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

// Apply periodic BCs in one direction.
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, sgr.lower_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.upper_ghost[dir]);

  gkyl_array_copy_to_buffer(buff->data, fld, sgr.upper_skin[dir]);
  gkyl_array_copy_from_buffer(fld, buff->data, sgr.lower_ghost[dir]);
}

void
test_fem_poisson_perp_consteps(int poly_order, const int *cells, struct gkyl_poisson_bc bcs, bool use_gpu)
{
  double epsilon_0 = 1.0;
  double lower[] = {0.,0.,-M_PI}, upper[] = {1.,1.,M_PI};
//  double lower[] = {0.,0.,-M_PI}, upper[] = {1.,1.,M_PI};
  int dim = sizeof(lower)/sizeof(lower[0]);

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);
  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &localRange_ext, ghost);

  // Projection updater for DG field.
  gkyl_proj_on_basis *projob, *projob_sol;
  if ((bcs.lo_type[0]==GKYL_POISSON_DIRICHLET && bcs.up_type[0]==GKYL_POISSON_DIRICHLET) &&
      (bcs.lo_type[1]==GKYL_POISSON_DIRICHLET && bcs.up_type[1]==GKYL_POISSON_DIRICHLET)) {
    projob = gkyl_proj_on_basis_new(&grid, &basis,
      poly_order+1, 1, evalFunc_dirichletx_dirichlety, NULL);
    projob_sol = gkyl_proj_on_basis_new(&grid, &basis,
      2*(poly_order+1), 1, evalFunc_dirichletx_dirichlety_sol, NULL);
  }

  // Create DG field we wish to make continuous.
  struct gkyl_array *rho = mkarr(basis.num_basis, localRange_ext.volume);
  // Create array holding continuous field we'll compute.
  struct gkyl_array *phi = mkarr(basis.num_basis, localRange_ext.volume);
  // Create DG field for permittivity tensor.
  int epsnum = PERP_DIM+ceil((pow(3.,PERP_DIM-1)-PERP_DIM)/2);
  struct gkyl_array *eps = use_gpu? mkarr_cu(epsnum*basis.num_basis, localRange_ext.volume)
                                  : mkarr(   epsnum*basis.num_basis, localRange_ext.volume);
  // Analytic solution.
  struct gkyl_array *phisol = mkarr(basis.num_basis, localRange_ext.volume);
  // Device copies:
  struct gkyl_array *rho_cu, *phi_cu;
  if (use_gpu) {
    rho_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
    phi_cu = mkarr_cu(basis.num_basis, localRange_ext.volume);
  }

  // Project RHS charge density on basis.
  gkyl_proj_on_basis_advance(projob, 0.0, &localRange, rho);
  struct gkyl_array *perbuff = mkarr(basis.num_basis, skin_ghost.lower_skin[dim-1].volume);
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, rho, d, skin_ghost);
  if (use_gpu) gkyl_array_copy(rho_cu, rho);
//  gkyl_grid_sub_array_write(&grid, &localRange, rho, "ctest_fem_poisson_perp_2x_rho_1.gkyl");

  // Project the permittivity onto the basis.
  double dg0norm = pow(sqrt(2.),dim);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 0*basis.num_basis);
  gkyl_array_shiftc(eps,        0.*dg0norm, 1*basis.num_basis);
  gkyl_array_shiftc(eps, epsilon_0*dg0norm, 2*basis.num_basis);

  // Project the analytic solution.
  gkyl_proj_on_basis_advance(projob_sol, 0.0, &localRange, phisol);
  gkyl_grid_sub_array_write(&grid, &localRange, phisol, "ctest_fem_poisson_perp_phisol_8x8x8_p1.gkyl");

  // FEM poisson solver.
  gkyl_fem_poisson_perp *poisson = gkyl_fem_poisson_perp_new(&grid, basis, &bcs, eps, NULL, use_gpu);

  // Set the RHS source.
  if (use_gpu)
    gkyl_fem_poisson_perp_set_rhs(poisson, rho_cu);
  else
    gkyl_fem_poisson_perp_set_rhs(poisson, rho);

  // Solve the problem.
  if (use_gpu) {
    gkyl_fem_poisson_perp_solve(poisson, phi_cu);
    gkyl_array_copy(phi, phi_cu);
#ifdef GKYL_HAVE_CUDA
    cudaDeviceSynchronize();
#endif
  } else {
    gkyl_fem_poisson_perp_solve(poisson, phi);
  }
  for (int d=0; d<dim; d++)
    if (bcs.lo_type[d] == GKYL_POISSON_PERIODIC) apply_periodic_bc(perbuff, phi, d, skin_ghost);
  gkyl_grid_sub_array_write(&grid, &localRange, phi, "ctest_fem_poisson_perp_phi_8x8x8_p1.gkyl");

  if (poly_order == 1) {
    if ((bcs.lo_type[0] == GKYL_POISSON_DIRICHLET && bcs.up_type[0] == GKYL_POISSON_DIRICHLET) &&
        (bcs.lo_type[1] == GKYL_POISSON_DIRICHLET && bcs.up_type[1] == GKYL_POISSON_DIRICHLET)) {
      // Solution; checked convergence:
      const double sol[512] = {
        -1.4857902425065184e-04, -8.5782139646996372e-05, -8.5782139646865672e-05, 1.8914704406197892e-04, -4.9526341416772864e-05, 1.0920409680558287e-04,
        1.0920409680550554e-04, 6.3049014687256998e-05, -3.7869881831633913e-04, -2.1864186469671895e-04, -4.7077585402820190e-05, 4.8209875139245301e-04,
        -2.7180256605144123e-05, 2.7833984389239624e-04, 5.9931650281304105e-05, 3.4601554422903200e-05, -4.8238508035628331e-04, -2.7850515599677118e-04,
        -1.2785705897207671e-05, 6.1409551253425998e-04, -7.3818307415179943e-06, 3.5454820946987623e-04, 1.6276715296076932e-05, 9.3973659577562421e-06,
        -4.7821812488875947e-04, -2.7609936313588445e-04, 1.5191498758083683e-05, 6.0879081146090802e-04, 8.7708158973654205e-06, 3.5148553887719326e-04,
        -1.9339385888747259e-05, -1.1165599648880918e-05, -3.9208040614287438e-04, -2.2636772803058295e-04, 3.4540136347233381e-05, 4.9913404823207755e-04,
        1.9941757017920681e-05, 2.8817517710849280e-04, -4.3970975879881029e-05, -2.5386654760770837e-05, -2.5807022372215224e-04, -1.4899691313582957e-04,
        4.2830678547521361e-05, 3.2853372287024223e-04, 2.4728303788985857e-05, 1.8967903333702248e-04, -5.4525167891612326e-05, -3.1480120359827638e-05,
        -1.1891747993225567e-04, -6.8657039050248661e-05, 3.7509195538050490e-05, 1.5138671107812069e-04, 2.1655944140984915e-05, 8.7403158392697900e-05,
        -4.7750707052706831e-05, -2.7568883570882659e-05, -2.6974823754671981e-05, -1.5573921756102582e-05, 1.5573921756086665e-05, 3.4340030182754069e-05,
        8.9916079182148039e-06, 1.9826225669992933e-05, -1.9826225669985302e-05, -1.1446676727580283e-05, -3.7869881831637274e-04, -4.7077585402708903e-05,
        -2.1864186469670288e-04, 4.8209875139248608e-04, -2.7180256605210324e-05, 5.9931650281245863e-05, 2.7833984389239038e-04, 3.4601554422944427e-05,
        -9.6593915317034899e-04, -1.2040150071024966e-04, -1.2040150071024630e-04, 1.2296791992509439e-03, -1.5153325635635661e-05, 1.5327592891960440e-04,
        1.5327592891959123e-04, 1.9290789976223312e-05, -1.2313170233906195e-03, -1.5389090291881229e-04, -3.2814484098396489e-05, 1.5675158485681691e-03,
        -4.1817897444776702e-06, 1.9590927819016461e-04, 4.1774151505829872e-05, 5.3235857015899117e-06, -1.2211748055640191e-03, -1.5284687644820770e-04,
        3.8670096290763775e-05, 1.5546043993793627e-03, 4.7845587083223804e-06, 1.9458019071073231e-04, -4.9228580170615438e-05, -6.0909347156166817e-06,
        -1.0016098602570024e-03, -1.2554426638121464e-04, 8.8095783986847693e-05, 1.2750894369280114e-03, 1.0978610563435414e-05, 1.5982287543424525e-04,
        -1.1214945864320036e-04, -1.3976210615571825e-05, -6.5979684339352812e-04, -8.2940058872077652e-05, 1.0925005331179741e-04, 8.3994778696934870e-04,
        1.3618940110575512e-05, 1.0558601423795712e-04, -1.3907969009605144e-04, -1.7337455796105482e-05, -3.0491150826769746e-04, -3.8726663276881513e-05,
        9.5643090454551556e-05, 3.8816455270337567e-04, 1.1907675738098372e-05, 4.9300592207781632e-05, -1.2175748182276777e-04, -1.5158947764468026e-05,
        -6.9626408103776196e-05, -9.0509819492167675e-06, 4.0198825461396472e-05, 8.8637203992388735e-05, 5.2255868648144249e-06, 1.1522262245235528e-05,
        -5.1174713585175428e-05, -6.6523812089424644e-06, -4.8238508035620151e-04, -1.2785705897252288e-05, -2.7850515599680664e-04, 6.1409551253430368e-04,
        -7.3818307414815423e-06, 1.6276715296141225e-05, 3.5454820946987097e-04, 9.3973659577153728e-06, -1.2313170233906072e-03, -3.2814484098386040e-05,
        -1.5389090291881698e-04, 1.5675158485681897e-03, -4.1817897444823730e-06, 4.1774151505828591e-05, 1.9590927819015664e-04, 5.3235857015929212e-06,
        -1.5708705777981958e-03, -4.2150433122696279e-05, -4.2150433122698495e-05, 1.9997811123959371e-03, -1.2083229378437410e-06, 5.3659188242121851e-05,
        5.3659188242117866e-05, 1.5382434574345833e-06, -1.5587645637327890e-03, -4.2060661306192221e-05, 4.9139843602174565e-05, 1.9843696720032254e-03,
        1.2601527202679231e-06, 5.3544905127007844e-05, -6.2556987501449741e-05, -1.6042248447095148e-06, -1.2792161864866266e-03, -3.4731820796201924e-05,
        1.1225748725242052e-04, 1.6284934001327763e-03, 2.9711553210228856e-06, 4.4214997854706335e-05, -1.4290827386115453e-04, -3.7823996304679761e-06,
        -8.4357750896041524e-04, -2.3165757864813304e-05, 1.3925862048675021e-04, 1.0739079292105724e-03, 3.7065142258784857e-06, 2.9490936864656485e-05,
        -1.7728179706443526e-04, -4.7185409457017731e-06, -3.9114247992314450e-04, -1.1058811420873844e-05, 1.2195486531873429e-04, 4.9794003061819580e-04,
        3.2834345625943677e-06, 1.4078309516759583e-05, -1.5525342423245637e-04, -4.1799436025258761e-06, -8.9955228480428413e-05, -2.6858679675509360e-06,
        5.1935675378201352e-05, 1.1451660589930813e-04, 1.5506865940911167e-06, 3.4192174121908953e-06, -6.6116193242601623e-05, -1.9740860933638614e-06,
        -4.7821812488872613e-04, 1.5191498758100369e-05, -2.7609936313588413e-04, 6.0879081146098728e-04, 8.7708158973496454e-06, -1.9339385888790925e-05,
        3.5148553887716583e-04, -1.1165599648852878e-05, -1.2211748055639972e-03, 3.8670096290759126e-05, -1.5284687644821456e-04, 1.5546043993793848e-03,
        4.7845587083257541e-06, -4.9228580170613290e-05, 1.9458019071072680e-04, -6.0909347156182750e-06, -1.5587645637327831e-03, 4.9139843602175452e-05,
        -4.2060661306194789e-05, 1.9843696720032327e-03, 1.2601527202678208e-06, -6.2556987501449768e-05, 5.3544905127004523e-05, -1.6042248447091459e-06,
        -1.5474982642617466e-03, 4.8565262338569481e-05, 4.8565262338568776e-05, 1.9700272218948593e-03, -1.5918873674820159e-06, -6.1825522557911871e-05,
        -6.1825522557912508e-05, 2.0265363267635107e-06, -1.2706869384245629e-03, 3.9656184460819470e-05, 1.1125183114826664e-04, 1.6176353260059105e-03,
        -3.5517711434683191e-06, -5.0483909874733382e-05, -1.4162803339381479e-04, 4.5215468089140504e-06, -8.3879554668975795e-04, 2.5926625069031862e-05,
        1.3810078013050519e-04, 1.0678202998640151e-03, -4.3749936672351580e-06, -3.3005631306618052e-05, -1.7580781995370960e-04, 5.5695420273585997e-06,
        -3.8988499721200347e-04, 1.1784819402652368e-05, 1.2107784645253936e-04, 4.9633920480195821e-04, -3.7897816410764171e-06, -1.5002546732690018e-05,
        -1.5413694411715089e-04, 4.8245436976423641e-06, -9.0086007742574192e-05, 2.6103625253466520e-06, 5.2011180820413507e-05, 1.1468309313383139e-04,
        -1.5070935067046508e-06, -3.3230959625222529e-06, -6.6212314692308015e-05, 1.9185903485009741e-06, -3.9208040614284722e-04, 3.4540136347213093e-05,
        -2.2636772803059100e-04, 4.9913404823208253e-04, 1.9941757017931554e-05, -4.3970975879880263e-05, 2.8817517710850375e-04, -2.5386654760776736e-05,
        -1.0016098602569918e-03, 8.8095783986845972e-05, -1.2554426638121602e-04, 1.2750894369280275e-03, 1.0978610563435211e-05, -1.1214945864320599e-04,
        1.5982287543424081e-04, -1.3976210615569642e-05, -1.2792161864866210e-03, 1.1225748725241937e-04, -3.4731820796203604e-05, 1.6284934001327808e-03,
        2.9711553210233968e-06, -1.4290827386115594e-04, 4.4214997854704112e-05, -3.7823996304677249e-06, -1.2706869384245618e-03, 1.1125183114826646e-04,
        3.9656184460818718e-05, 1.6176353260059101e-03, -3.5517711434682683e-06, -1.4162803339381555e-04, -5.0483909874733890e-05, 4.5215468089141537e-06,
        -1.0441173998824859e-03, 9.1153799606617992e-05, 9.1153799606618073e-05, 1.3292032360397348e-03, -8.0518327772840401e-06, -1.1604243491016443e-04,
        -1.1604243491016380e-04, 1.0250305362999503e-05, -6.9005433895489503e-04, 5.9949151256281410e-05, 1.1326460393002962e-04, 8.7846678973585207e-04,
        -9.9641793477503907e-06, -7.6317668737877270e-05, -1.4419037369695810e-04, 1.2684799080065410e-05, -3.2161392308227690e-04, 2.7631503624010198e-05,
        9.9454569354364179e-05, 4.0942739519954688e-04, -8.6944232123168641e-06, -3.5176009937019723e-05, -1.2660964699911679e-04, 1.1068348703517426e-05,
        -7.4676777957819192e-05, 6.2861604388864368e-06, 4.3114657856159156e-05, 9.5066526934358737e-05, -3.6293164215592856e-06, -8.0025338133193964e-06,
        -5.4886684916474010e-05, 4.6202650513188232e-06, -2.5807022372216449e-04, 4.2830678547518969e-05, -1.4899691313582071e-04, 3.2853372287025562e-04,
        2.4728303788984804e-05, -5.4525167891608260e-05, 1.8967903333701682e-04, -3.1480120359831297e-05, -6.5979684339352498e-04, 1.0925005331179476e-04,
        -8.2940058872077706e-05, 8.3994778696934924e-04, 1.3618940110576470e-05, -1.3907969009605478e-04, 1.0558601423795538e-04, -1.7337455796106072e-05,
        -8.4357750896041372e-04, 1.3925862048674929e-04, -2.3165757864814100e-05, 1.0739079292105683e-03, 3.7065142258784857e-06, -1.7728179706443878e-04,
        2.9490936864655540e-05, -4.7185409457013005e-06, -8.3879554668975860e-04, 1.3810078013050441e-04, 2.5926625069031462e-05, 1.0678202998640096e-03,
        -4.3749936672350818e-06, -1.7580781995371169e-04, -3.3005631306617849e-05, 5.5695420273589242e-06, -6.9005433895489578e-04, 1.1326460393002904e-04,
        5.9949151256281458e-05, 8.7846678973584904e-04, -9.9641793477503653e-06, -1.4419037369695927e-04, -7.6317668737876267e-05, 1.2684799080065638e-05,
        -4.5700790975607651e-04, 7.4600267375337498e-05, 7.4600267375337959e-05, 5.8178935875592588e-04, -1.2358685770139086e-05, -9.4969125900869786e-05,
        -9.4969125900869000e-05, 1.5733101584855904e-05, -2.1405750025111815e-04, 3.4466226050632516e-05, 6.5667216952079626e-05, 2.7250376448509594e-04,
        -1.0812713792346902e-05, -4.3876885114389393e-05, -8.3596995207867956e-05, 1.3765017386719724e-05, -5.0159272049234077e-05, 7.8690281972934216e-06,
        2.8959469219980895e-05, 6.3854760712488500e-05, -4.5431855479681580e-06, -1.0017587816761052e-05, -3.6866563286394175e-05, 5.7836570226375672e-06,
        -1.1891747993226251e-04, 3.7509195538055992e-05, -6.8657039050244378e-05, 1.5138671107812595e-04, 2.1655944140983332e-05, -4.7750707052715593e-05,
        8.7403158392689037e-05, -2.7568883570880843e-05, -3.0491150826769691e-04, 9.5643090454552613e-05, -3.8726663276881581e-05, 3.8816455270336027e-04,
        1.1907675738097432e-05, -1.2175748182277376e-04, 4.9300592207778536e-05, -1.5158947764468221e-05, -3.9114247992314591e-04, 1.2195486531873341e-04,
        -1.1058811420874908e-05, 4.9794003061817552e-04, 3.2834345625941889e-06, -1.5525342423246220e-04, 1.4078309516759917e-05, -4.1799436025256102e-06,
        -3.8988499721200818e-04, 1.2107784645253757e-04, 1.1784819402651616e-05, 4.9633920480194390e-04, -3.7897816410767364e-06, -1.5413694411715395e-04,
        -1.5002546732686913e-05, 4.8245436976436923e-06, -3.2161392308228053e-04, 9.9454569354363149e-05, 2.7631503624011527e-05, 4.0942739519954103e-04,
        -8.6944232123161492e-06, -1.2660964699911720e-04, -3.5176009937017941e-05, 1.1068348703517651e-05, -2.1405750025111896e-04, 6.5667216952079626e-05,
        3.4466226050632868e-05, 2.7250376448509437e-04, -1.0812713792346967e-05, -8.3596995207868051e-05, -4.3876885114388695e-05, 1.3765017386719678e-05,
        -1.0163877852160976e-04, 3.0438753201852812e-05, 3.0438753201852917e-05, 1.2939023268193514e-04, -9.5264492369834572e-06, -3.8749750996843554e-05,
        -3.8749750996843378e-05, 1.2127551130928198e-05, -2.4458655728475410e-05, 6.9692295538358886e-06, 1.4121211468851588e-05, 3.1136847587379800e-05,
        -4.0236865589514017e-06, -8.8721081333432448e-06, -1.7976867336290102e-05, 5.1223140190652031e-06, -2.6974823754661329e-05, 1.5573921756091266e-05,
        -1.5573921756096677e-05, 3.4340030182743030e-05, 8.9916079182173179e-06, -1.9826225669985950e-05, 1.9826225669986620e-05, -1.1446676727580623e-05,
        -6.9626408103756491e-05, 4.0198825461406474e-05, -9.0509819492174502e-06, 8.8637203992336435e-05, 5.2255868648150289e-06, -5.1174713585190749e-05,
        1.1522262245218026e-05, -6.6523812089505900e-06, -8.9955228480439079e-05, 5.1935675378196955e-05, -2.6858679675677970e-06, 1.1451660589925285e-04,
        1.5506865940821917e-06, -6.6116193242616016e-05, 3.4192174122066818e-06, -1.9740860933552026e-06, -9.0086007742595225e-05, 5.2011180820405816e-05,
        2.6103625253575236e-06, 1.1468309313381435e-04, -1.5070935066976264e-06, -6.6212314692306511e-05, -3.3230959625159697e-06, 1.9185903485014870e-06,
        -7.4676777957820601e-05, 4.3114657856161548e-05, 6.2861604388869069e-06, 9.5066526934354780e-05, -3.6293164215604896e-06, -5.4886684916472533e-05,
        -8.0025338133181208e-06, 4.6202650513183065e-06, -5.0159272049234396e-05, 2.8959469219981146e-05, 7.8690281972935876e-06, 6.3854760712487497e-05,
        -4.5431855479681952e-06, -3.6866563286393782e-05, -1.0017587816760623e-05, 5.7836570226374503e-06, -2.4458655728475423e-05, 1.4121211468851677e-05,
        6.9692295538358971e-06, 3.1136847587379671e-05, -4.0236865589514525e-06, -1.7976867336290004e-05, -8.8721081333431685e-06, 5.1223140190651497e-06,
        -6.1937980258105384e-06, 3.5759909575078880e-06, 3.5759909575078880e-06, 7.8849527650921135e-06, -2.0645993419368465e-06, -4.5523796014734159e-06,
        -4.5523796014734159e-06, 2.6283175883640376e-06,
      };
      long i = 0;
      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &localRange);
      while (gkyl_range_iter_next(&iter)) {
        long loc = gkyl_range_idx(&localRange, iter.idx);
        const double *phi_p = gkyl_array_cfetch(phi, loc);
        if (iter.idx[2] == 3) {
          // Only check one cell in z:
          for (int m=0; m<basis.num_basis; m++) {
            TEST_CHECK( gkyl_compare(sol[i], phi_p[m], 1e-12) );
            TEST_MSG("Expected: %.13e in cell (%d,%d,%d)", sol[i], iter.idx[0], iter.idx[1], iter.idx[2]);
            TEST_MSG("Produced: %.13e", phi_p[m]);
            i += 1;
          }
        }
      }
    }
  } else if (poly_order == 2) {
  }

  gkyl_fem_poisson_perp_release(poisson);
  gkyl_proj_on_basis_release(projob);
  gkyl_proj_on_basis_release(projob_sol);
  gkyl_array_release(perbuff);
  gkyl_array_release(rho);
  gkyl_array_release(eps);
  gkyl_array_release(phi);
  gkyl_array_release(phisol);
  if (use_gpu) {
    gkyl_array_release(rho_cu);
    gkyl_array_release(phi_cu);
  }

}

void test_p1_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps(1, cells, bc_tv, false);
}

void test_p2_dirichletx_dirichlety_consteps() {
  int cells[] = {8,8,8};
  struct gkyl_poisson_bc bc_tv;
  bc_tv.lo_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[0] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.up_type[1] = GKYL_POISSON_DIRICHLET;
  bc_tv.lo_value[0].v[0] = 0.;
  bc_tv.up_value[0].v[0] = 0.;
  bc_tv.lo_value[1].v[0] = 0.;
  bc_tv.up_value[1].v[0] = 0.;
  test_fem_poisson_perp_consteps(2, cells, bc_tv, false);
}

TEST_LIST = {
  { "test_p1_dirichletx_dirichlety", test_p1_dirichletx_dirichlety_consteps },
//  { "test_p2_dirichletx_dirichlety", test_p2_dirichletx_dirichlety_consteps },
  { NULL, NULL },
};

