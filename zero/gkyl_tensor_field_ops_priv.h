#pragma once

// Private header, not for direct use in user code

#include <gkyl_tensor_field.h>

void 
gkyl_tensor_field_raise_or_lower_idx_in_place(struct gkyl_tensor_field *met, int raised_idx, struct gkyl_tensor_field *ten)
{

  // temporary tensor at a single point 
  enum gkyl_tensor_index_loc iloc = { GKYL_TENSOR_INDEX_LOWER };
  struct gkyl_tensor_field *tensor_out = gkyl_tensor_field_new(ten->rank,ten->ndim,1,&iloc);

  // iterate over the field of tensors
  for (long i=0; i<ten->size; ++i){
    const double *metric = gkyl_tensor_field_cfetch(met, i);

    // loop over the output tensor indices
    struct gkyl_range_iter iter_tf_out;
    gkyl_range_iter_init(&iter_tf_out, &ten->trange);
    while (gkyl_range_iter_next(&iter_tf_out)) {
      double ten_out_elem = 0; 
      int index_raised = iter_tf_out.idx[raised_idx];
      
      // summed over index, j
      for (int j=0; j<ten->ndim; ++j) {
        int idx_met[GKYL_MAX_DIM] = { index_raised, j };
        
        // Get the tensor element we are indexing from
        int idx_tf[GKYL_MAX_DIM];
        for (int k=0; k<GKYL_MAX_DIM; ++k) 
          idx_tf[k] = (raised_idx != k) ?  iter_tf_out.idx[k] : j;

        double met_elem = gkyl_tensor_field_elem_fetch(met, i, idx_met);
        double ten_elem = gkyl_tensor_field_elem_fetch(ten, i, idx_tf);
        ten_out_elem += met_elem*ten_elem;
      }

      // set the resulting sum in the temporary tensor
      gkyl_tensor_field_elem_set(tensor_out, i, iter_tf_out.idx, ten_out_elem);
    }

    // Set the old tensor to the new tensor values
    struct gkyl_range_iter iter_tf;
    gkyl_range_iter_init(&iter_tf, &ten->trange);
    while (gkyl_range_iter_next(&iter_tf)) {
      const double ten_out_elem = gkyl_tensor_field_elem_fetch(tensor_out, i, iter_tf.idx);
      gkyl_tensor_field_elem_set(ten, i, iter_tf.idx, ten_out_elem);
    }
  }

  // release
  gkyl_tensor_field_release(tensor_out);
}

void 
gkyl_tensor_field_raise_or_lower_idx_set(const struct gkyl_tensor_field *met, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *tensor_out)
{

  // iterate over the field of tensors
  for (long i=0; i<ten->size; ++i){
    const double *metric = gkyl_tensor_field_cfetch(met, i);

    // loop over the output tensor indices
    struct gkyl_range_iter iter_tf_out;
    gkyl_range_iter_init(&iter_tf_out, &ten->trange);
    while (gkyl_range_iter_next(&iter_tf_out)) {
      double ten_out_elem = 0; 
      int index_raised = iter_tf_out.idx[raised_idx];
      
      // summed over index, j
      for (int j=0; j<ten->ndim; ++j) {
        int idx_met[GKYL_MAX_DIM] = { index_raised, j };
        
        // Get the tensor element we are indexing from
        int idx_tf[GKYL_MAX_DIM];
        for (int k=0; k<GKYL_MAX_DIM; ++k) 
          idx_tf[k] = (raised_idx != k) ?  iter_tf_out.idx[k] : j;

        double met_elem = gkyl_tensor_field_elem_fetch(met, i, idx_met);
        double ten_elem = gkyl_tensor_field_elem_fetch(ten, i, idx_tf);
        ten_out_elem += met_elem*ten_elem;
      }

      // set the resulting sum in the temporary tensor
      gkyl_tensor_field_elem_set(tensor_out, i, iter_tf_out.idx, ten_out_elem);
    }
  }
}