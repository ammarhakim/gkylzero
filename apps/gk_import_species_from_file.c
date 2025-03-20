#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_file_import_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_neut_species *gkns, bool is_gk, struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file. Intended options include importing:
  //   1) ICs with same grid.
  //   2) ICs one dimensionality lower (e.g. 2x2v for 3x2v sim).
  //   3) ICs with same grid extents but different resolution (NYI).

  struct gkyl_rect_grid grid = is_gk ? gks->grid : gkns->grid;
  const char *basis_id = is_gk ? gks->basis.id : gkns->basis.id;
  enum gkyl_basis_type b_type = is_gk ? gks->basis.b_type : gkns->basis.b_type;
  struct gkyl_range local_s = is_gk ? gks->local : gkns->local;
  struct gkyl_array *f_s = is_gk ? gks->f : gkns->f;
  struct gkyl_basis basis = is_gk ? gks->basis : gkns->basis;
  
  int pdim = grid.ndim;
  int cdim = app->cdim;
  int vdim = pdim - cdim;
  int poly_order = app->poly_order;

  struct gkyl_rect_grid grid_do; // Donor grid.
  struct gkyl_array_header_info hdr;
  int pdim_do, vdim_do, cdim_do;
  bool same_res = true;

  // Read the header of the input file, extract needed info an create a grid
  // and other things needed.
  FILE *fp;
  with_file(fp, inp.file_name, "r") {

    int status = gkyl_grid_sub_array_header_read_fp(&grid_do, &hdr, fp);

    pdim_do = grid_do.ndim;
    vdim_do = vdim; // Assume velocity space dimensionality is the same.
    cdim_do = pdim_do - vdim_do;

    // Perform some basic checks.
    if (pdim_do == pdim) {
      for (int d=0; d<pdim; d++) {
        assert(grid_do.lower[d] == grid.lower[d]);
        assert(grid_do.upper[d] == grid.upper[d]);
      }
      // Check if the grid resolution is the same.
      for (int d=0; d<pdim; d++)
        same_res = same_res && (grid_do.cells[d] == grid.cells[d]);
    }
    else {
      // Assume the loaded file has one lower conf-space dimension.
      // Primarily meant for loading:
      //   - 1x2v for a 2x2v sim / - 1x3v for a 2x3v sim.
      //   - 2x2v for a 3x2v sim / - 2x3v for a 3x3v sim.
      assert(pdim_do == pdim-1);
      for (int d=0; d<cdim_do-1; d++) {
        assert(grid_do.lower[d] == grid.lower[d]);
        assert(grid_do.upper[d] == grid.upper[d]);
        assert(grid_do.cells[d] == grid.cells[d]);
        assert(grid_do.dx[d] == grid.dx[d]);
      }
      assert(grid_do.lower[cdim_do-1] == grid.lower[cdim-1]);
      assert(grid_do.upper[cdim_do-1] == grid.upper[cdim-1]);
      assert(grid_do.cells[cdim_do-1] == grid.cells[cdim-1]);
      assert(grid_do.dx[cdim_do-1] == grid.dx[cdim-1]);
      for (int d=0; d<vdim; d++) {
        assert(grid_do.lower[cdim_do+d] == grid.lower[cdim+d]);
        assert(grid_do.upper[cdim_do+d] == grid.upper[cdim+d]);
        assert(grid_do.cells[cdim_do+d] == grid.cells[cdim+d]);
        assert(grid_do.dx[cdim_do+d] == grid.dx[cdim+d]);
      }
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }
      );
    assert(strcmp(basis_id, meta.basis_type_nm) == 0);
    assert(poly_order == meta.poly_order);
    gkyl_grid_sub_array_header_release(&hdr);
  }

  // Donor basis.
  struct gkyl_basis basis_do;
  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&basis_do, pdim_do, poly_order);
      }
      else if (poly_order == 1) {
        // p=2 in vparallel
        gkyl_cart_modal_gkhybrid(&basis_do, cdim_do, vdim_do); 
      }
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&basis_do, pdim_do, poly_order); // for canonical PB
      break;
    default:
      assert(false);
      break;
  }

  // Donor global range.
  int ghost_do[pdim_do];
  for (int d=0; d<cdim_do; d++) ghost_do[d] = 1;
  for (int d=0; d<vdim_do; d++) ghost_do[cdim_do+d] = 0;
  struct gkyl_range global_ext_do, global_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &global_ext_do, &global_do);

  // Create a donor communicator.
  int cuts_tar[GKYL_MAX_DIM] = {-1}, cuts_do[GKYL_MAX_CDIM] = {-1};
  gkyl_rect_decomp_get_cuts(app->decomp, cuts_tar);
  if (cdim_do == cdim-1) {
    for (int d=0; d<cdim_do-1; d++) {
      cuts_do[d] = cuts_tar[d];
    }
    cuts_do[cdim_do-1] = cuts_tar[cdim-1];
  }
  else {
    for (int d=0; d<cdim; d++) {
      cuts_do[d] = cuts_tar[d];
    }
  }
  // Set velocity space cuts to 1 as we do not use MPI in vel-space.
  for (int d=0; d<vdim; d++) {
    cuts_do[cdim_do+d] = 1;
  }

  struct gkyl_rect_decomp *decomp_do = gkyl_rect_decomp_new_from_cuts(pdim_do, cuts_do, &global_do);
  struct gkyl_comm* comm_do = is_gk ? gkyl_comm_split_comm(gks->comm, 0, decomp_do)
    : gkyl_comm_split_comm(gkns->comm, 0, decomp_do);

  // Donor local range.
  int my_rank = 0;
  gkyl_comm_get_rank(comm_do, &my_rank);

  struct gkyl_range local_ext_do, local_do;
  gkyl_create_ranges(&decomp_do->ranges[my_rank], ghost_do, &local_ext_do, &local_do);

  // Donor array.
  struct gkyl_array *fdo = mkarr(app->use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_host = app->use_gpu? mkarr(false, basis_do.num_basis, local_ext_do.volume)
                                            : gkyl_array_acquire(fdo);

  // Read donor field.
  struct gkyl_app_restart_status rstat;
  rstat.io_status = gkyl_comm_array_read(comm_do, &grid_do, &local_do, fdo_host, inp.file_name);
  if (app->use_gpu) {
    gkyl_array_copy(fdo, fdo_host);
  }

  if (pdim_do == pdim-1) {
    struct gkyl_translate_dim_gyrokinetic* transdim = gkyl_translate_dim_gyrokinetic_new(cdim_do,
      basis_do, cdim, basis, app->use_gpu);
    gkyl_translate_dim_gyrokinetic_advance(transdim, &local_do, &local_s, fdo, f_s);
    gkyl_translate_dim_gyrokinetic_release(transdim);
  }
  else {
    if (same_res) {
      gkyl_array_copy(f_s, fdo);
    }
    else {
      // Interpolate the donor distribution to the target grid.
      struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(app->cdim, &basis,
        &grid_do, &grid, &local_do, &local_s, ghost_do, app->use_gpu);
      gkyl_dg_interpolate_advance(interp, fdo, f_s);
      gkyl_dg_interpolate_release(interp);
    }
  }

  if (inp.type == GKYL_IC_IMPORT_AF || inp.type == GKYL_IC_IMPORT_AF_B) {
    // Scale f by a conf-space factor.
    gkyl_proj_on_basis *proj_conf_scale = gkyl_proj_on_basis_new(&app->grid, &app->basis,
      poly_order+1, 1, inp.conf_scale, inp.conf_scale_ctx);
    struct gkyl_array *xfac = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *xfac_ho = app->use_gpu? mkarr(false, app->basis.num_basis, app->local_ext.volume)
                                             : gkyl_array_acquire(xfac_ho);
    gkyl_proj_on_basis_advance(proj_conf_scale, 0.0, &app->local, xfac_ho);
    gkyl_array_copy(xfac, xfac_ho);
    gkyl_dg_mul_conf_phase_op_range(&app->basis, &basis, f_s, xfac, f_s, &app->local, &local_s);
    gkyl_proj_on_basis_release(proj_conf_scale);
    gkyl_array_release(xfac_ho);
    gkyl_array_release(xfac);
  }
  if (inp.type == GKYL_IC_IMPORT_F_B || inp.type == GKYL_IC_IMPORT_AF_B) {
    // Add a phase factor to f.
    struct gk_proj proj_phase_add;
    if (is_gk) {
      gk_species_projection_init(app, gks, inp.phase_add, &proj_phase_add);
      gk_species_projection_calc(app, gks, &proj_phase_add, gks->fnew, 0.0);
      gkyl_array_accumulate_range(f_s, 1.0, gks->fnew, &local_s);
      gk_species_projection_release(app, &proj_phase_add);
    }
    else {
      gk_neut_species_projection_init(app, gkns, inp.phase_add, &proj_phase_add);
      gk_neut_species_projection_calc(app, gkns, &proj_phase_add, gkns->fnew, 0.0);
      gkyl_array_accumulate_range(f_s, 1.0, gkns->fnew, &gkns->local);
      gk_neut_species_projection_release(app, &proj_phase_add);
    }
  }

  gkyl_rect_decomp_release(decomp_do);
  gkyl_comm_release(comm_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_host);
}
