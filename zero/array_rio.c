#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <gkyl_array_rio.h>
#include <gkyl_elem_type_priv.h>

void
gkyl_array_write(const struct gkyl_array *arr, FILE *fp)
{
  uint64_t esznc = arr->esznc, size = arr->size;
  fwrite(&esznc, sizeof(uint64_t), 1, fp);
  fwrite(&size, sizeof(uint64_t), 1, fp);
  fwrite(arr->data, arr->esznc*arr->size, 1, fp);
}

void
gkyl_sub_array_write(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    fwrite(_F(start), arr->esznc*skip.delta, 1, fp);
  }
#undef _F
}

int
gkyl_grid_sub_array_header_write_fp(const struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, FILE *fp)
{
  const char g0[5] = "gkyl0";

  // Version 1 header
  fwrite(g0, sizeof(char[5]), 1, fp);
  uint64_t version = 1;
  fwrite(&version, sizeof(uint64_t), 1, fp);
  fwrite(&hdr->file_type, sizeof(uint64_t), 1, fp);
  uint64_t meta_size = 0; // THIS WILL CHANGE ONCE METADATA IS EMBEDDED
  fwrite(&meta_size, sizeof(uint64_t), 1, fp);
  
  // Version 0 format is used for rest of the header
  uint64_t real_type = gkyl_array_data_type[hdr->etype];
  fwrite(&real_type, sizeof(uint64_t), 1, fp);
  gkyl_rect_grid_write(grid, fp);

  fwrite(&hdr->esznc, sizeof(uint64_t), 1, fp);
  fwrite(&hdr->tot_cells, sizeof(uint64_t), 1, fp);

  return errno;
}

int
gkyl_grid_sub_array_write_fp(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
  gkyl_grid_sub_array_header_write_fp(grid,
    &(struct gkyl_array_header_info) {
      .file_type = gkyl_file_type_int[GKYL_FIELD_DATA_FILE],
      .etype = arr->type,
      .esznc = arr->esznc,
      .tot_cells = range->volume
    },
    fp
  );

  gkyl_sub_array_write(range, arr, fp);
  return errno;
}

int
gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "w") {
    err = gkyl_grid_sub_array_write_fp(grid, range, arr, fp);
  }
  return err;
}

struct gkyl_array*
gkyl_array_new_from_file(enum gkyl_elem_type type, FILE *fp)
{
  struct gkyl_array* arr = 0;
  
  uint64_t esznc, size;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return 0;
  if (1 != fread(&size, sizeof(uint64_t), 1, fp))
    return 0;

  int ncomp = esznc/gkyl_elem_type_size[type];
  arr = gkyl_array_new(type, ncomp, size);
  if (1 != fread(arr->data, arr->esznc*arr->size, 1, fp)) {
    gkyl_array_release(arr);
    arr = 0;
  }
  return arr;
}

bool
gkyl_sub_array_read(const struct gkyl_range *range, struct gkyl_array *arr, FILE *fp)
{
#define _F(loc) gkyl_array_fetch(arr, loc)
  
  uint64_t esznc, size;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return false;
  if (1 != fread(&size, sizeof(uint64_t), 1, fp))
    return false;

  if ((size != range->volume) || (size > arr->size))
    return false;
  
  // construct skip iterator to allow reading (potentially) in chunks
  // rather than element by element or requiring a copy of data
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    if (1 != fread(_F(start), arr->esznc*skip.delta, 1, fp))
      return false;
  }

  return true;
#undef _F
}

int
gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char* fname)
{
  FILE *fp = 0;
  with_file (fp, fname, "r") {

    size_t frr;
    // Version 1 header
    
    char g0[6];
    frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
    g0[5] = '\0'; // add the NULL
    if (strcmp(g0, "gkyl0") != 0)
      return 1;

    uint64_t version;
    frr = fread(&version, sizeof(uint64_t), 1, fp);
    if (version != 1)
      return 1;

    uint64_t file_type;
    frr = fread(&file_type, sizeof(uint64_t), 1, fp);
    if (file_type != 1)
      return 1;

    uint64_t meta_size;
    frr = fread(&meta_size, sizeof(uint64_t), 1, fp);

    // read ahead by specified bytes: meta-data is not read in this
    // method
    fseek(fp, meta_size, SEEK_CUR);
    
    uint64_t real_type = 0;
    if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
      break;
    
    gkyl_rect_grid_read(grid, fp);
    gkyl_sub_array_read(range, arr, fp);
  }
  return errno;  
}

struct gkyl_array*
gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid, const char* fname)
{
  struct gkyl_array *arr = 0;
  FILE *fp = 0;
  with_file (fp, fname, "r") {

    size_t frr;
    // Version 1 header
    
    char g0[6];
    frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
    g0[5] = '\0'; // add the NULL
    if (strcmp(g0, "gkyl0") != 0)
      return 0;

    uint64_t version;
    frr = fread(&version, sizeof(uint64_t), 1, fp);
    if (version != 1)
      return 0;

    uint64_t file_type;
    frr = fread(&file_type, sizeof(uint64_t), 1, fp);
    if (file_type != 1)
      return 0;

    uint64_t meta_size;
    frr = fread(&meta_size, sizeof(uint64_t), 1, fp);

    // read ahead by specified bytes: READ META-DATA HERE
    fseek(fp, meta_size, SEEK_CUR);    
    
    uint64_t real_type = 0;
    if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
      break;
    
    gkyl_rect_grid_read(grid, fp);
    arr = gkyl_array_new_from_file(real_type, fp);
  }
  return arr;
}
