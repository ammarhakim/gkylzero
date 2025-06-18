#include <gkyl_tensor_field_ops_priv.h>
#include <gkyl_tensor_field_ops.h>
#include <gkyl_tensor_field.h>

#include <assert.h>

void 
gkyl_tensor_field_lower_idx_in_place(struct gkyl_tensor_field *metric, int lowered_idx, struct gkyl_tensor_field *ten)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[lowered_idx] == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_LOWER );
  gkyl_tensor_field_raise_or_lower_idx_in_place(metric, lowered_idx, ten);
  ten->iloc[lowered_idx] = GKYL_TENSOR_INDEX_LOWER;
}

void 
gkyl_tensor_field_raise_idx_in_place(struct gkyl_tensor_field *metric, int raised_idx, struct gkyl_tensor_field *ten)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[raised_idx] == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_UPPER );
  gkyl_tensor_field_raise_or_lower_idx_in_place(metric, raised_idx, ten);
  ten->iloc[raised_idx] = GKYL_TENSOR_INDEX_UPPER;
}

void 
gkyl_tensor_field_lower_idx_set(const struct gkyl_tensor_field *metric, int lowered_idx,
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *ten_out)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[lowered_idx] == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_LOWER );
  gkyl_tensor_field_raise_or_lower_idx_set(metric, lowered_idx, ten, ten_out);
  ten_out->iloc[lowered_idx] = GKYL_TENSOR_INDEX_LOWER;
}

void 
gkyl_tensor_field_raise_idx_set(const struct gkyl_tensor_field *metric, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *ten_out)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[raised_idx] == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_UPPER );
  gkyl_tensor_field_raise_or_lower_idx_set(metric, raised_idx, ten, ten_out);
  ten_out->iloc[raised_idx] = GKYL_TENSOR_INDEX_UPPER;
}