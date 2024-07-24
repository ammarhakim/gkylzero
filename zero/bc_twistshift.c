#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_mat.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_priv.h>
#include <gkyl_util.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_array_ops.h>

struct gkyl_bc_twistshift*
gkyl_bc_twistshift_new(const struct gkyl_bc_twistshift_inp *inp)
{

  // Allocate space for new updater.
  struct gkyl_bc_twistshift *up = gkyl_malloc(sizeof(struct gkyl_bc_twistshift));

  up->bc_dir = inp->bc_dir;
  up->shift_dir = inp->shift_dir;
  up->shear_dir = inp->shear_dir;
  up->edge = inp->edge;
  up->basis = inp->basis;
  up->grid = inp->grid;
  up->shift_func     = inp->shift_func; // Function defining the shift.
  up->shift_func_ctx = inp->shift_func_ctx; // Context for shift_func.
  up->use_gpu = inp->use_gpu;

  // Assume the poly order of the DG shift is the same as that of the field,
  // unless requested otherwise.
  up->shift_poly_order = inp->basis.poly_order;
  if (inp->shift_poly_order)
    up->shift_poly_order = inp->shift_poly_order;

  up->local_ext_r = inp->local_ext_r;
  const int ndim = inp->local_ext_r.ndim;

  // Create a range only extended in bc_dir.
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = inp->local_ext_r.lower[d] + inp->num_ghost[d];
    upper_bcdir_ext[d] = inp->local_ext_r.upper[d] - inp->num_ghost[d];
  }
  lower_bcdir_ext[up->bc_dir] = inp->local_ext_r.lower[up->bc_dir];
  upper_bcdir_ext[up->bc_dir] = inp->local_ext_r.upper[up->bc_dir];
  gkyl_sub_range_init(&up->local_bcdir_ext_r, &inp->local_ext_r, lower_bcdir_ext, upper_bcdir_ext);

  double lo1d[1], up1d[1];  int cells1d[1];

  // Create 1D grid and range in the diretion of the shear.
  gkyl_range_init(&up->shear_r, 1, (int[]) {up->local_bcdir_ext_r.lower[inp->shear_dir]},
                                   (int[]) {up->local_bcdir_ext_r.upper[inp->shear_dir]});
  lo1d[0] = inp->grid.lower[up->shear_dir];
  up1d[0] = inp->grid.upper[up->shear_dir];
  cells1d[0] = inp->grid.cells[up->shear_dir];
  gkyl_rect_grid_init(&up->shear_grid, 1, lo1d, up1d, cells1d);

  // Create 1D grid and range in the diretion of the shift.
  gkyl_range_init(&up->shift_r, 1, (int[]) {up->local_bcdir_ext_r.lower[inp->shift_dir]},
                                   (int[]) {up->local_bcdir_ext_r.upper[inp->shift_dir]});
  lo1d[0] = inp->grid.lower[up->shift_dir];
  up1d[0] = inp->grid.upper[up->shift_dir];
  cells1d[0] = inp->grid.cells[up->shift_dir];
  gkyl_rect_grid_init(&up->shift_grid, 1, lo1d, up1d, cells1d);

  // Create 2D grid (and range) the twist-shift takes place in.
  int dimlo, dimup;
  if (up->shift_dir < up->shear_dir) {
    dimlo = up->shift_dir;
    dimup = up->shear_dir;
    up->shift_dir_in_ts_grid = 0;
    up->shear_dir_in_ts_grid = 1;
  }
  else {
    dimlo = up->shear_dir;
    dimup = up->shift_dir;
    up->shift_dir_in_ts_grid = 1;
    up->shear_dir_in_ts_grid = 0;
  }
  gkyl_range_init(&up->ts_r, 2, (int[]) {up->local_bcdir_ext_r.lower[dimlo], up->local_bcdir_ext_r.lower[dimup]},
                                (int[]) {up->local_bcdir_ext_r.upper[dimlo], up->local_bcdir_ext_r.upper[dimup]});
  double lo2d[] = {inp->grid.lower[dimlo], inp->grid.lower[dimup]};
  double up2d[] = {inp->grid.upper[dimlo], inp->grid.upper[dimup]};
  int cells2d[] = {inp->grid.cells[dimlo], inp->grid.cells[dimup]};
  gkyl_rect_grid_init(&up->ts_grid, 2, lo2d, up2d, cells2d);

  // Project the shift onto the shift basis.
  gkyl_cart_modal_serendip(&up->shift_b, 1, up->shift_poly_order);
  up->shift = gkyl_array_new(GKYL_DOUBLE, up->shift_b.num_basis, up->shear_r.volume);
  gkyl_eval_on_nodes *evup = gkyl_eval_on_nodes_new(&up->shear_grid, &up->shift_b, 1,
    inp->shift_func, inp->shift_func_ctx);
  gkyl_eval_on_nodes_advance(evup, 0.0, &up->shear_r, up->shift);
  gkyl_eval_on_nodes_release(evup);

  // Find the donor cells for each target. Store the number of donors for each
  // shear_dir idx (num_do) & the shift_dir idx of each donor (shift_dir_idx_do).
  // i.e. allocates and assigns num_do and up->shift_dir_idx_do.
  ts_find_donors(up);

  // Array of cummulative number of donors at given shear_dir cell.
  const int num_do_cum_sz = up->grid.cells[up->shear_dir]+1;
  int num_do_cum_ho[num_do_cum_sz];
  num_do_cum_ho[0] = 0;
  for (int i=1; i<num_do_cum_sz; i++)
    num_do_cum_ho[i] = num_do_cum_ho[i-1] + up->num_do[i-1];

  if (!up->use_gpu) {
    up->num_do_cum = gkyl_malloc(num_do_cum_sz * sizeof(int));
    memcpy(up->num_do_cum, num_do_cum_ho, num_do_cum_sz * sizeof(int));
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->num_do_cum = gkyl_cu_malloc(num_do_cum_sz * sizeof(int));
    cudaMemcpy(up->num_do_cum, num_do_cum_ho, num_do_cum_sz * sizeof(int), GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Choose the kernels that do the subcell and full cell integrals
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_twistshift_kernels));
  gkyl_bc_twistshift_choose_kernels(inp->basis, inp->cdim, up->kernels);

  // The BC is applied as a set of matrix-matrix multiplications
  //   f_i = sum_{q}^{N_do(i)} A_q,i B_q,i
  // where i indicates the shear_dir cell index, A_q is a
  // num_basis x num_basis matrix containing the discretization
  // of subcell integrals, B_q is a num_basis x (Ny * Nvpar * Nmu) matrix with
  // the DG coefficients of f common to a given A_q matrix, and thus where f_i
  // is a num_basis x (Ny*Nvpar*Nmu) matrix.
  //
  // Naming scheme:
  //   A_q: scimat (subscell integral matrices).
  //   B_q: fmat (distribution function matrices).
  //   A_q . B_q: mm_contr (contributions from mat-mat multiplication).

  // Calculate the entries in the matrices used to apply the BC.
  up->scimat = ts_calc_mats(up);

  // Number of colums in fmat.
  int fmat_num_col = 1;
  for (int d=0; d<ndim; d++) {
    if (d != up->bc_dir && d != up->shear_dir)
      fmat_num_col *= up->grid.cells[d];
  }

  if (!up->use_gpu) {
    up->fmat = gkyl_nmat_new(up->scimat->num, up->scimat->nr, fmat_num_col);
    up->mm_contr = gkyl_nmat_new(up->scimat->num, up->scimat->nr, fmat_num_col);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->fmat = gkyl_nmat_cu_dev_new(up->scimat->num, up->scimat->nr, fmat_num_col);
    up->mm_contr = gkyl_nmat_cu_dev_new(up->scimat->num, up->scimat->nr, fmat_num_col);
  }
#endif

  // Index translation from num-numcol plane index to linear index into the
  // donor distribution function gkyl_array.
  up->num_numcol_fidx_do = ts_calc_num_numcol_fidx_do(up);

  // Index translation from num-numcol plane index to linear index into the
  // tar distribution function gkyl_array.
  up->num_numcol_fidx_tar = ts_calc_num_numcol_fidx_tar(up);

  // Permutted ghost range, for indexing into the target field.
  int lo4D[4] = { up->local_bcdir_ext_r.lower[1], up->local_bcdir_ext_r.lower[3],
                  up->local_bcdir_ext_r.lower[4], up->local_bcdir_ext_r.lower[0] };
  int up4D[4] = { up->local_bcdir_ext_r.upper[1], up->local_bcdir_ext_r.upper[3],
                  up->local_bcdir_ext_r.upper[4], up->local_bcdir_ext_r.upper[0] };
  gkyl_range_init(&up->permutted_ghost_r, ndim-1, lo4D, up4D);

  // Create a ghost range, to clear it before adding contributions from TS BC.
  struct gkyl_range tmp_skin;
  gkyl_skin_ghost_ranges(&tmp_skin, &up->ghost_r, up->bc_dir, up->edge, &inp->local_ext_r, inp->num_ghost);

  return up;
}

void gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_twistshift_advance_cu(up, fdo, ftar);
    return;
  }
#endif

  // Assign the distribution matrices.
  // This assumes that fdo->ncomp = fmat->nr.
  // Recall:
  //   fmat->num = sum_i^Nx num_do(i)
  //   fmat->nr = num_basis = ncomp
  //   fmat->nc = Ny*Nvpar*Nmu
  for (size_t i=0; i<up->fmat->num * up->fmat->nr * up->fmat->nc; i++) {

    long nc_idx = i / up->fmat->nr; // num-num_col index: current num-num_col plane.
    int row_idx = i % up->fmat->nr; // row index: current DG coeff.

    // This if-statement may only be needed in GPU kernel, not for CPUs.
    if ((nc_idx < up->fmat->num * up->fmat->nc) && (row_idx < fdo->ncomp)) {
      const double *fdo_c = (const double*) gkyl_array_cfetch(fdo, up->num_numcol_fidx_do[nc_idx]);
      struct gkyl_mat mcurr = gkyl_nmat_get(up->fmat, nc_idx % up->fmat->num);

      gkyl_mat_set(&mcurr, row_idx, nc_idx/up->fmat->num, fdo_c[row_idx]);
    }
  }

  // Perform the mat-mat multiplications.
  gkyl_nmat_mm(1.0, 0.0, GKYL_NO_TRANS, up->scimat, GKYL_NO_TRANS, up->fmat, up->mm_contr);

  gkyl_array_clear_range(ftar, 0.0, &up->ghost_r);

  // Perform reduction over num_do contributions from mat-mat mults (mm_contr).
  int num_cells_skin = up->grid.cells[up->shear_dir] * up->fmat->nc;
  for (size_t i=0; i<ftar->ncomp * num_cells_skin; i++) {

    long linidx_tar = i / ftar->ncomp;
    int row_idx = i % ftar->ncomp;

    int idx[GKYL_MAX_DIM];
    gkyl_sub_range_inv_idx(&up->permutted_ghost_r, linidx_tar, idx);
    int y_idx    = idx[0];
    int vpar_idx = idx[1];
    int mu_idx   = idx[2];
    int x_idx    = idx[3];

    // This if-statement may only be needed in GPU kernel, not for CPUs.
    if (linidx_tar<num_cells_skin && row_idx<ftar->ncomp) {
      double *ftar_c = (double*) gkyl_array_fetch(ftar, up->num_numcol_fidx_tar[linidx_tar]);
      int Nvpar = up->grid.cells[3], Nmu = up->grid.cells[4];
      int start = (mu_idx-1) * up->mm_contr->num
        + (vpar_idx-1) * Nmu * up->mm_contr->num
        + (y_idx-1) * Nvpar * Nmu * up->mm_contr->num;

      int do_start = up->num_do_cum[x_idx-1];
      int do_end   = up->num_do_cum[x_idx-1+1];
      for (int j=do_start; j<do_end; j++) { // Only loop over num_do[i] elements.
        int linidx_mm_contr = start + j;
        struct gkyl_mat mat = gkyl_nmat_get(up->mm_contr, linidx_mm_contr % up->mm_contr->num);
        ftar_c[row_idx] += gkyl_mat_get(&mat, row_idx, linidx_mm_contr / up->mm_contr->num);
      }
    }
  }
}

void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up) {
  // Release memory associated with this updater.
  if (!up->use_gpu) {
    gkyl_free(up->num_do_cum);
    gkyl_free(up->num_numcol_fidx_do);
    gkyl_free(up->num_numcol_fidx_tar);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->num_do_cum);
    gkyl_cu_free(up->num_numcol_fidx_do);
    gkyl_cu_free(up->num_numcol_fidx_tar);
  }
#endif

  gkyl_nmat_release(up->fmat);
  gkyl_nmat_release(up->mm_contr);

  gkyl_nmat_release(up->scimat);

  gkyl_free(up->kernels);

  gkyl_array_release(up->shift);

  gkyl_free(up->num_do);
  gkyl_free(up->shift_dir_idx_do);

  gkyl_free(up);
}
