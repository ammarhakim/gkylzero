#include <gkyl_array_rio.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkylt_mirrorgridgen.h>

#include <stc/cstr.h>

bool
gkylt_mirrorgridgen(const struct gkylt_mirrorgridgen_inp *inp)
{
  double clower[3], cupper[3];
  int cells[3];
  for (int d=0; d<3; ++d) {
    clower[d] = inp->lower[d];
    cupper[d] = inp->upper[d];
    cells[d] = inp->cells[d];
  }
  
  struct gkyl_rect_grid comp_grid;
  gkyl_rect_grid_init(&comp_grid, 3, clower, cupper, cells);

  bool status = true;  
  if (!gkyl_check_file_exists(inp->psiRZ_fname)) {
    fprintf(stderr, "Unable to find file %s!\n", inp->psiRZ_fname);
    status = false;
    goto cleanup;
  }

  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, inp->psiRZ_fname);

    // write DG projection of mapc2p to file
  cstr cubic_fileNm = cstr_from_fmt("%s-psi-cubic.gkyl", inp->out_prefix);
  
  // create mirror geometry
  struct gkyl_mirror_grid_gen *geom =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &comp_grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        .nrnodes = psi_grid.cells[0]-1, // cells and not nodes
        .nznodes = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = inp->fl_coord,
        .include_axis = inp->include_axis,
        .write_psi_cubic = inp->write_psi_cubic,
        .psi_cubic_fname = cubic_fileNm.str
      }
    );

  cstr_drop(&cubic_fileNm);

  // write out node coordinates
  struct gkyl_rect_grid nodal_grid;
  gkyl_rect_grid_init(&nodal_grid, 2,
    (double[]) { 0.0, 0.0 }, // arbitrary
    (double[]) { 1.0, 1.0 }, // arbitrary
    (int[]) { cells[0]+1, cells[2]+1 }
  );

  struct gkyl_range node_range;
  gkyl_range_init_from_shape(&node_range, 2, (int[2]) { cells[0]+1, cells[2]+1 });

  do {
    cstr fileNm = cstr_from_fmt("%s-nodal_coords.gkyl", inp->out_prefix);
    enum gkyl_array_rio_status io_status =
      gkyl_grid_sub_array_write(&nodal_grid, &node_range, 0, geom->nodesrz, fileNm.str);
    cstr_drop(&fileNm);
  } while (0);
  

  gkyl_mirror_grid_gen_release(geom);
  gkyl_array_release(psi);

  cleanup:
  
  return status;
}
